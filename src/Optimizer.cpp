#ifndef _OPTIMIZER_CPP_
#define _OPTIMIZER_CPP_

#include "Optimizer.h"
#include "TaoInterfaceRegistration.h"


namespace reg
{




/********************************************************************
 * Name: Optimizer
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Optimizer"
Optimizer::Optimizer()
{
    this->Initialize();
}




/********************************************************************
 * Name: Optimizer
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~Optimizer"
Optimizer::~Optimizer(void)
{
    this->ClearMemory();
}




/********************************************************************
 * Name: Optimizer
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Optimizer"
Optimizer::Optimizer(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * Name: Initialize
 * Description: init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode Optimizer::Initialize(void)
{
    PetscFunctionBegin;

    this->m_Tao = NULL;
    this->m_OptimizationProblem = NULL;

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: ClearMemory
 * Description: clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode Optimizer::ClearMemory(void)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if(this->m_Tao != NULL){
        ierr=TaoDestroy(&this->m_Tao); CHKERRQ(ierr);
        this->m_Tao=NULL;
    }


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SetProblem
 * Description: set the optimization problem
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetProblem"
PetscErrorCode Optimizer::SetProblem(Optimizer::OptProbType* optprob)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr=Assert(optprob!=NULL,"opt prob is null"); CHKERRQ(ierr);

    this->m_OptimizationProblem = optprob;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: DoSetup
 * Description: set the optimization problem
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DoSetup"
PetscErrorCode Optimizer::DoSetup()
{
    PetscErrorCode ierr;
    IntType nlukwn,ngukwn;
    ScalarType gatol,grtol,gttol,reltol,abstol,divtol;
    IntType maxit;
    Vec v0;
    Mat HMatVec;
    KSP taoksp;
    PC taokktpc;
    PetscFunctionBegin;

    ierr=Assert(this->m_OptimizationProblem !=NULL, "optimizer: optimization problem not set"); CHKERRQ(ierr);

    nlukwn = 3*this->m_Opt->GetNLocal();
    ngukwn = 3*this->m_Opt->GetNGlobal();

    // allocate vector fields
    ierr=VecCreate(PETSC_COMM_WORLD,&v0); CHKERRQ(ierr);
    ierr=VecSetSizes(v0,nlukwn,ngukwn); CHKERRQ(ierr);
    ierr=VecSetFromOptions(v0); CHKERRQ(ierr);
    ierr=VecSet(v0,0.0); CHKERRQ(ierr);

    if(this->m_Tao != NULL){
        ierr=TaoDestroy(&this->m_Tao); CHKERRQ(ierr);
        this->m_Tao=NULL;
    }

    std::string method = "nls";
    ierr=TaoCreate(PETSC_COMM_WORLD,&this->m_Tao); CHKERRQ(ierr);
    ierr=TaoSetType(this->m_Tao,"nls"); CHKERRQ(ierr);
    ierr=TaoSetInitialVector(this->m_Tao,v0); CHKERRQ(ierr);

    // set the routine to evaluate the objective and compute the gradient
    ierr=TaoSetObjectiveRoutine(this->m_Tao,EvaluateObjective,(void*)this->m_OptimizationProblem); CHKERRQ(ierr);
    ierr=TaoSetGradientRoutine(this->m_Tao,EvaluateGradient,(void*)this->m_OptimizationProblem); CHKERRQ(ierr);
    ierr=TaoSetObjectiveAndGradientRoutine(this->m_Tao,EvaluateObjectiveGradient,(void*)this->m_OptimizationProblem); CHKERRQ(ierr);

    // set the monitor for the optimization process
    ierr=TaoCancelMonitors(this->m_Tao); CHKERRQ(ierr);
    ierr=TaoSetMonitor(this->m_Tao,OptimizationMonitor,this->m_OptimizationProblem,NULL); CHKERRQ(ierr);

    TaoLineSearch ls;
    ierr=TaoGetLineSearch(this->m_Tao,&ls); CHKERRQ(ierr);
    ierr=TaoLineSearchSetType(ls,"armijo"); CHKERRQ(ierr);


    // set tolearances for optimizer
    gatol = this->m_Opt->GetOptTol(0);
    grtol = this->m_Opt->GetOptTol(1);
    gttol = this->m_Opt->GetOptTol(2);
    ierr=TaoSetTolerances(this->m_Tao,gatol,grtol,gttol); CHKERRQ(ierr);

    ierr=TaoSetMaximumIterations(this->m_Tao,this->m_Opt->GetOptMaxit() - 1); CHKERRQ(ierr);

    ierr=MatCreateShell(PETSC_COMM_WORLD,nlukwn,nlukwn,ngukwn,ngukwn,static_cast<void*>(this->m_OptimizationProblem),&HMatVec); CHKERRQ(ierr);
    ierr=MatShellSetOperation(HMatVec,MATOP_MULT,(void(*)(void))HessianMatVec); CHKERRQ(ierr);
    ierr=MatSetOption(HMatVec,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);
    ierr=TaoSetHessianRoutine(this->m_Tao,HMatVec,HMatVec,EvaluateHessian,static_cast<void*>(&this->m_OptimizationProblem)); CHKERRQ(ierr);


    // get the ksp of the optimizer and set options
    ierr=TaoGetKSP(this->m_Tao,&taoksp); CHKERRQ(ierr);

    // ksp is only nonzero if we use a newton type method
    if (taoksp != NULL){

        if(this->m_Opt->GetVerbosity() >= 2){
            ierr=KSPMonitorSet(taoksp,KrylovMonitor,this->m_OptimizationProblem,NULL); CHKERRQ(ierr);
        }

        // set the preconditioner
        ierr=KSPGetPC(taoksp,&taokktpc); CHKERRQ(ierr);

        if(this->m_Opt->GetPrecondMeth() == NOPC){
            if (strcmp(method.c_str(),"nls") == 0){
                ierr=PetscOptionsSetValue(NULL,"-tao_nls_pc_type","none"); CHKERRQ(ierr);
                ierr=PetscOptionsSetValue(NULL,"-tao_nls_ksp_type","cg"); CHKERRQ(ierr);
                ierr=TaoSetFromOptions(this->m_Tao); CHKERRQ(ierr);
            }
            else if (strcmp(method.c_str(),"ntr") == 0){
                ierr=PetscOptionsSetValue(NULL,"-tao_ntr_pc_type","none"); CHKERRQ(ierr);
                ierr=PetscOptionsSetValue(NULL,"-tao_ntr_ksp_type","stcg"); CHKERRQ(ierr);
                ierr=TaoSetFromOptions(this->m_Tao); CHKERRQ(ierr);
            }
            ierr=PCSetType(taokktpc,PCNONE); CHKERRQ(ierr);
        }
        else if (  (this->m_Opt->GetPrecondMeth() == INVREG)
                || (this->m_Opt->GetPrecondMeth() == TWOLEVEL) ) {
            if (strcmp(method.c_str(),"nls") == 0){
                ierr=PetscOptionsSetValue(NULL,"-tao_nls_pc_type","petsc"); CHKERRQ(ierr);
                ierr=PetscOptionsSetValue(NULL,"-tao_nls_ksp_type","cg"); CHKERRQ(ierr);
                ierr=TaoSetFromOptions(this->m_Tao); CHKERRQ(ierr);
            }
            else if (strcmp(method.c_str(),"ntr") == 0){
                ierr=PetscOptionsSetValue(NULL,"-tao_ntr_pc_type","petsc"); CHKERRQ(ierr);
                ierr=PetscOptionsSetValue(NULL,"-tao_ntr_ksp_type","stcg"); CHKERRQ(ierr);
                ierr=TaoSetFromOptions(this->m_Tao); CHKERRQ(ierr);
            }
            ierr=PCSetType(taokktpc,PCSHELL); CHKERRQ(ierr);
            ierr=PCShellSetApply(taokktpc,PrecondMatVec); CHKERRQ(ierr);
            ierr=PCShellSetContext(taokktpc,this->m_OptimizationProblem); CHKERRQ(ierr);
            //ierr=PCShellSetName(taokktpc,"kktpc"); CHKERRQ(ierr);
            ierr=PCShellSetSetUp(taokktpc,PrecondSetup); CHKERRQ(ierr);
        }
        else{
            ierr=reg::ThrowError("preconditioner not defined"); CHKERRQ(ierr);
        }

        // set tolerances for krylov subspace method
        reltol = this->m_Opt->GetKKTSolverTol(0); // 1E-12;
        abstol = this->m_Opt->GetKKTSolverTol(1); // 1E-12;
        divtol = this->m_Opt->GetKKTSolverTol(2); // 1E+06;
        maxit  = this->m_Opt->GetKKTMaxit(); // 1000;
        ierr=KSPSetTolerances(taoksp,reltol,abstol,divtol,maxit); CHKERRQ(ierr);
        //ierr=KSPSetInitialGuessNonzero(taoksp,PETSC_FALSE); CHKERRQ(ierr);
        ierr=KSPSetInitialGuessNonzero(taoksp,PETSC_TRUE); CHKERRQ(ierr);

        //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
        //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
        //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
        ierr=KSPSetNormType(taoksp,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    }


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Run
 * Description: run the optimizer
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Run"
PetscErrorCode Optimizer::Run()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // check of optimization problem has been set
    ierr=Assert(this->m_OptimizationProblem !=NULL, "optimizer: optimization problem not set"); CHKERRQ(ierr);

    // do the setup
    if (this->m_Tao == NULL){
        ierr=this->DoSetup(); CHKERRQ(ierr);
    }
    ierr=Assert(this->m_Tao !=NULL, "optimizer: problem in tao setup"); CHKERRQ(ierr);

    // do the inversion
    ierr=this->m_Opt->StartTimer(T2SEXEC); CHKERRQ(ierr);

    ierr=Msg("starting optimization"); CHKERRQ(ierr);
    ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);
    ierr=Msg("optimization done"); CHKERRQ(ierr);

    ierr=this->m_Opt->StopTimer(T2SEXEC); CHKERRQ(ierr);

    // display info to user, once we're done
    ierr=TaoView(this->m_Tao,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Run
 * Description: run the optimizer
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunBetaCont"
PetscErrorCode Optimizer::RunBetaCont()
{
    PetscErrorCode ierr;
    std::string levelmsg;
    std::stringstream levelss;
    ScalarType beta;
    Vec x;
    unsigned int maxsteps = 10;
    bool stop;
    PetscFunctionBegin;

    // check of optimization problem has been set
    ierr=Assert(this->m_OptimizationProblem !=NULL, "optimizer: optimization problem not set"); CHKERRQ(ierr);

    // do the setup
    if (this->m_Tao == NULL){
        ierr=this->DoSetup(); CHKERRQ(ierr);
    }
    ierr=Assert(this->m_Tao !=NULL, "optimizer: problem in tao setup"); CHKERRQ(ierr);

    // do the inversion
    ierr=this->m_Opt->StartTimer(T2SEXEC); CHKERRQ(ierr);

    levelmsg = "level ";
    ierr=Msg("starting optimization (parameter continuation)"); CHKERRQ(ierr);

    beta = 1.0;

    for(unsigned int i = 0; i < maxsteps; ++i){

        this->m_Opt->SetRegularizationWeight(beta);

        levelss << std::setw(3) << i <<" (beta="<<beta<<")";
        ierr=Msg(levelmsg + levelss.str()); CHKERRQ(ierr);
        levelss.str("");

        // solve optimization probelm for current regularization parameter
        ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);

        ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);
        stop=false; // not yet we've not
        ierr=this->m_OptimizationProblem->CheckBounds(x,stop); CHKERRQ(ierr);

        if (stop) break; // if bound reached go home

        beta /= 10.0;

    }
    ierr=Msg("optimization done"); CHKERRQ(ierr);

    ierr=this->m_Opt->StopTimer(T2SEXEC); CHKERRQ(ierr);

    // display info to user, once we're done
    ierr=TaoView(this->m_Tao,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}







/********************************************************************
 * Name: GetSolution
 * Description: get the solution
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetSolution"
PetscErrorCode Optimizer::GetSolution(Vec &x)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(this->m_Tao !=NULL, "optimization has not been performed"); CHKERRQ(ierr);

    ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Finalize
 * Description: finalize
 * *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Finalize"
PetscErrorCode Optimizer::Finalize()
{
    PetscErrorCode ierr;
    Vec x;
    PetscFunctionBegin;

    ierr=Assert(this->m_Tao !=NULL, "optimization has not been performed"); CHKERRQ(ierr);
    ierr=Assert(this->m_OptimizationProblem !=NULL, "optimizer: optimization problem not set"); CHKERRQ(ierr);

    ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);

    // finalize the registration (write out all data)
    ierr=this->m_OptimizationProblem->Finalize(x); CHKERRQ(ierr);

    // display info to user, once we're done
    ierr=TaoView(this->m_Tao,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


} // end of name space



#endif
