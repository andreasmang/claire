/**
 * @file Optimizer.cpp
 *
 * @author Andreas Mang
 */




#ifndef _OPTIMIZER_CPP_
#define _OPTIMIZER_CPP_

#include "Optimizer.hpp"




namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Optimizer"
Optimizer::Optimizer()
{
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~Optimizer"
Optimizer::~Optimizer(void)
{
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 * @param opt base class for registration options and arguments
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Optimizer"
Optimizer::Optimizer(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode Optimizer::Initialize(void)
{
    PetscFunctionBegin;

    this->m_Tao = NULL;
    this->m_Solution = NULL;
    this->m_OptimizationProblem = NULL;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief clean up
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

    if (this->m_Solution != NULL){
        ierr=VecDestroy(&this->m_Solution); CHKERRQ(ierr);
        this->m_Solution = NULL;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set the initial guess (warm start)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetInitialGuess"
PetscErrorCode Optimizer::SetInitialGuess(VecField* x)
{
    PetscErrorCode ierr;
    IntType nlu,ngu;

    PetscFunctionBegin;

    // compute the number of unknowns
    nlu = 3*this->m_Opt->GetDomainPara().nlocal;
    ngu = 3*this->m_Opt->GetDomainPara().nglobal;

    if(this->m_Solution==NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_Solution); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_Solution,nlu,ngu); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_Solution); CHKERRQ(ierr);
        ierr=VecSet(this->m_Solution,0.0); CHKERRQ(ierr);
    }

    // the input better is not zero
    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=x->GetComponents(this->m_Solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief parse initial guess to tao
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetInitialGuess"
PetscErrorCode Optimizer::SetInitialGuess()
{
    PetscErrorCode ierr;
    IntType nlu,ngu;

    PetscFunctionBegin;

    nlu = 3*this->m_Opt->GetDomainPara().nlocal;
    ngu = 3*this->m_Opt->GetDomainPara().nglobal;

    // check if tao has been set up
    ierr=Assert(this->m_Tao!=NULL,"tao is null"); CHKERRQ(ierr);

    if(this->m_Solution==NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_Solution); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_Solution,nlu,ngu); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_Solution); CHKERRQ(ierr);
        ierr=VecSet(this->m_Solution,0.0); CHKERRQ(ierr);
    }

    // parse initial guess to tao
    ierr=TaoSetInitialVector(this->m_Tao,this->m_Solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them; currently this is only supported for
 * registration problems)
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
 * @brief set up optimization problem (this function sets up
 * and parses the default parameters to TAO)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DoSetup"
PetscErrorCode Optimizer::SetupTao()
{
    PetscErrorCode ierr;
    IntType nlu,ngu;
    ScalarType gatol,grtol,gttol,reltol,abstol,divtol;
    IntType maxit;
    Mat HMatVec;
    KSP taoksp;
    PC taokktpc;
    PetscFunctionBegin;

    ierr=Assert(this->m_OptimizationProblem !=NULL,"optimization problem not set"); CHKERRQ(ierr);

    // compute the number of unknowns
    nlu = 3*this->m_Opt->GetDomainPara().nlocal;
    ngu = 3*this->m_Opt->GetDomainPara().nglobal;

    // if tao exists, kill it
    if(this->m_Tao != NULL){
        ierr=TaoDestroy(&this->m_Tao); CHKERRQ(ierr);
        this->m_Tao=NULL;
    }

    std::string method = "nls";
    ierr=TaoCreate(PETSC_COMM_WORLD,&this->m_Tao); CHKERRQ(ierr);
    ierr=TaoSetType(this->m_Tao,"nls"); CHKERRQ(ierr);

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
    gatol = this->m_Opt->GetOptTol(0);  // ||g(x)||              <= gatol
    grtol = this->m_Opt->GetOptTol(1);  // ||g(x)|| / |J(x)|     <= grtol
    gttol = this->m_Opt->GetOptTol(2);  // ||g(x)|| / ||g(x0)||  <= gttol
    ierr=TaoSetTolerances(this->m_Tao,gatol,grtol,gttol); CHKERRQ(ierr);

    ierr=TaoSetMaximumIterations(this->m_Tao,this->m_Opt->GetOptMaxit() - 1); CHKERRQ(ierr);

    ierr=MatCreateShell(PETSC_COMM_WORLD,nlu,nlu,ngu,ngu,static_cast<void*>(this->m_OptimizationProblem),&HMatVec); CHKERRQ(ierr);
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
        ierr=KSPSetInitialGuessNonzero(taoksp,PETSC_FALSE); CHKERRQ(ierr);
        //ierr=KSPSetInitialGuessNonzero(taoksp,PETSC_TRUE); CHKERRQ(ierr);

        //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
        //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
        //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
        ierr=KSPSetNormType(taoksp,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
        //ierr=KSPSetNormType(taoksp,KSP_NORM_PRECONDITIONED); CHKERRQ(ierr);
    }


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief run the optimizer (main interface; calls specific functions
 * according to user settings (parameter continuation, grid
 * continuation, scale continuation, ...))
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Run"
PetscErrorCode Optimizer::Run()
{
    PetscErrorCode ierr;
    int rank;
    Vec x;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // check if optimization problem has been set
    ierr=Assert(this->m_OptimizationProblem !=NULL,"optimization problem not set"); CHKERRQ(ierr);

    // do setup
    if (this->m_Tao == NULL){
        ierr=this->SetupTao(); CHKERRQ(ierr);
    }
    ierr=Assert(this->m_Tao !=NULL,"problem in tao setup"); CHKERRQ(ierr);

    // set initial guess
    ierr=this->SetInitialGuess(); CHKERRQ(ierr);

    // solve optimization problem
    ierr=this->m_Opt->StartTimer(T2SEXEC); CHKERRQ(ierr);
    ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);
    ierr=this->m_Opt->StopTimer(T2SEXEC); CHKERRQ(ierr);

    // get solution
    ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);

    // copy solution into place holder
    ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}




/********************************************************************
 * @param x vector to hold solution
 * @brief get the solution
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetSolution"
PetscErrorCode Optimizer::GetSolution(Vec &x)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // check if we have solved the problem / set up tao
    ierr=Assert(this->m_Tao!=NULL, "optimization object not initialized"); CHKERRQ(ierr);

    // get solution
    ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief finalize optimization (displays information for user)
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Finalize"
PetscErrorCode Optimizer::Finalize()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(this->m_Tao !=NULL,"tao not set up"); CHKERRQ(ierr);

    // display info to user, once we're done
    ierr=TaoView(this->m_Tao,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of name space




#endif
