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
    ierr=TaoSetConvergenceTest(this->m_Tao,CheckConvergence,this->m_OptimizationProblem); CHKERRQ(ierr);

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
    int rank,indent,numindent,maxiter,iter,linelength;
    ScalarType gatol,grtol,gttol,gnorm,J,g0norm;
    bool stop[3],converged;
    std::string line,msg;
    TaoConvergedReason reason;
    std::stringstream ss;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    ierr=Assert(this->m_Tao !=NULL,"tao not set up"); CHKERRQ(ierr);

    ierr=TaoGetTolerances(this->m_Tao,&gatol,&grtol,&gttol); CHKERRQ(ierr);
    g0norm = this->m_OptimizationProblem->GetInitialGradNorm();

    ierr=TaoGetMaximumIterations(this->m_Tao,&maxiter); CHKERRQ(ierr);
    ierr=TaoGetSolutionStatus(this->m_Tao,&iter,&J,&gnorm,NULL,NULL,&reason); CHKERRQ(ierr);

    linelength = this->m_Opt->GetLineLength();
    line = std::string(linelength,'-');

    stop[0] = (gnorm < gttol*g0norm);// relative change of gradient
    stop[1] = (gnorm < gatol);       // absolute norm of gradient
    stop[2] = (iter  > maxiter);     // max number of iterations

    // check if we converged
    converged=false;
    for (int i = 0; i < 3; ++i){ if (stop[i]) converged=true; }

    indent    = 25;
    numindent = 5;
    if (rank == 0){

        if (converged){
            std::cout<< " convergence criteria" <<std::endl;

            // relative change of gradient
            ss << "[ " << stop[0] << "    ||g|| = " << std::setw(14) <<
                std::right << std::scientific << gnorm << " < " <<
                std::left << std::setw(14) << gttol*g0norm << " = " << "tol";
            std::cout << std::left << std::setw(100) << ss.str() << "]" << std::endl;
            ss.str(std::string()); ss.clear();

            // absolute norm of gradient
            ss << "[ " << stop[1] << "    ||g|| = " << std::setw(14) <<
                std::right << std::scientific << gnorm << " < "  <<
                std::left << std::setw(14) << gatol << " = " << "tol";
            std::cout << std::left << std::setw(100) << ss.str() << "]" << std::endl;
            ss.str(std::string()); ss.clear();

            // number of iterations
            ss << "[ " << stop[2] << "     iter = " << std::setw(14) <<
                std::right << iter  << " > " <<
                std::left << std::setw(14) << maxiter << " = " << "maxiter";
            std::cout << std::left << std::setw(100) << ss.str() << "]" << std::endl;
            ss.str(std::string()); ss.clear();
        }
    }

    if (!converged){
        switch(reason){
            case TAO_CONVERGED_STEPTOL:
            {
                msg="line search failed";
                break;
            }
            case TAO_CONVERGED_MINF:
            {
                msg="objective value to small";
                break;
            }
            case TAO_DIVERGED_NAN:
            {
                msg="numerical issues (NaN detected)";
                break;
            }
            case TAO_DIVERGED_MAXFCN:
            {
                msg="maximal number of function evaluations reached";
                break;
            }
            case TAO_DIVERGED_LS_FAILURE:
            {
                msg="line search failed";
                break;
            }
            case TAO_DIVERGED_TR_REDUCTION:
            {
                msg="trust region failed";
                break;
            }
            default:{
                msg="did not converge; reason not defined";
                break;
            }
        }
        ierr=WrngMsg(msg); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 1){

        if (converged && !rank) std::cout<< line <<std::endl;
        ss << std::left << std::setw(indent)
           << "outer iterations"
           << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(ITERATIONS) - 1;
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ss << std::left << std::setw(indent)
           << "objective evals"
           << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(OBJEVAL);
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ss << std::left << std::setw(indent)
           << "hessian matvecs"
           << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(HESSMATVEC);
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ss << std::left << std::setw(indent)
           << "pde solves"
           << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(PDESOLVE);
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // display info to user, once we're done
    //ierr=TaoView(this->m_Tao,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of name space




#endif
