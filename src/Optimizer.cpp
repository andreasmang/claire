/**
 * @file Optimizer.cpp
 *
 * @author Andreas Mang
 */


#ifndef _OPTIMIZER_CPP_
#define _OPTIMIZER_CPP_

#include "Optimizer.hpp"
#include "TaoInterfaceRegistration.hpp"




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
    this->m_InitialGuess = NULL;
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

    if (this->m_InitialGuess != NULL){
        ierr=VecDestroy(&this->m_InitialGuess); CHKERRQ(ierr);
        this->m_InitialGuess = NULL;
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
PetscErrorCode Optimizer::SetInitialGuess(Vec x0)
{
    PetscErrorCode ierr;
    IntType nlukwn,ngukwn;

    PetscFunctionBegin;

    ierr=Assert(x0!=NULL,"opt prob is null"); CHKERRQ(ierr);

    nlukwn = 3*this->m_Opt->GetNLocal();
    ngukwn = 3*this->m_Opt->GetNGlobal();

    // create all zero initial guess, if it has not been set already
    if (this->m_InitialGuess == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_InitialGuess); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_InitialGuess,nlukwn,ngukwn); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_InitialGuess); CHKERRQ(ierr);
    }

    ierr=VecCopy(x0,this->m_InitialGuess); CHKERRQ(ierr);

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
PetscErrorCode Optimizer::DoSetup()
{
    PetscErrorCode ierr;
    IntType nlukwn,ngukwn;
    ScalarType gatol,grtol,gttol,reltol,abstol,divtol;
    IntType maxit;
    Mat HMatVec;
    KSP taoksp;
    PC taokktpc;
    PetscFunctionBegin;

    ierr=Assert(this->m_OptimizationProblem !=NULL, "optimizer: optimization problem not set"); CHKERRQ(ierr);

    nlukwn = 3*this->m_Opt->GetNLocal();
    ngukwn = 3*this->m_Opt->GetNGlobal();

    // create all zero initial guess, if it has not been set already
    if (this->m_InitialGuess == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_InitialGuess); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_InitialGuess,nlukwn,ngukwn); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_InitialGuess); CHKERRQ(ierr);
        ierr=VecSet(this->m_InitialGuess,0.0); CHKERRQ(ierr);
    }
    if (this->m_Solution != NULL){
        ierr=VecDestroy(&this->m_Solution); CHKERRQ(ierr);
        this->m_Solution = NULL;
    }
    ierr=VecCreate(PETSC_COMM_WORLD,&this->m_Solution); CHKERRQ(ierr);
    ierr=VecSetSizes(this->m_Solution,nlukwn,ngukwn); CHKERRQ(ierr);
    ierr=VecSetFromOptions(this->m_Solution); CHKERRQ(ierr);

    // store the best we have so far in solution vector
    ierr=VecCopy(this->m_InitialGuess,this->m_Solution); CHKERRQ(ierr);

    if(this->m_Tao != NULL){
        ierr=TaoDestroy(&this->m_Tao); CHKERRQ(ierr);
        this->m_Tao=NULL;
    }

    std::string method = "nls";
    ierr=TaoCreate(PETSC_COMM_WORLD,&this->m_Tao); CHKERRQ(ierr);
    ierr=TaoSetType(this->m_Tao,"nls"); CHKERRQ(ierr);
    ierr=TaoSetInitialVector(this->m_Tao,this->m_InitialGuess); CHKERRQ(ierr);

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

    if( this->m_Opt->DoRegParaReductionSearch() ){

        // start the parameter continuation (we start with
        // a large regularization paramameter until we
        // hit bound on jacobian or have reached desired
        // regularization weight)
        ierr=this->RunRegParaReductionSearch(); CHKERRQ(ierr);

    }
    else if( this->m_Opt->DoRegParaBinarySearch() ){

        // start the parameter continuation (we first reduce
        // the regularization parameter by one order of magnitude
        // and from there start a binary search)
        ierr=this->RunRegParaBinarySearch(); CHKERRQ(ierr);

    }
    else if( this->m_Opt->DoRegParaContinuation() ){

        // start the parameter continuation (we reduce the
        // regularization parameter until we hit the user
        // defined target parameter)
        ierr=this->RunRegParaContinuation(); CHKERRQ(ierr);

    }
    else{

        // solve
        ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);

        // get solution
        ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);

        // copy solution into place holder
        ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);

    }

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ierr=Msg("optimization done"); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;

    ierr=this->m_Opt->StopTimer(T2SEXEC); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief run the optimizer; we search for an optimal
 * regularization weight using a binary search
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunRegParaBinarySearch"
PetscErrorCode Optimizer::RunRegParaBinarySearch()
{
    PetscErrorCode ierr;
    std::string levelmsg;
    std::stringstream ss;
    ScalarType beta,betamin,betascale,dbetascale,
                betastar,betahat,dbeta,dbetamin;
    Vec x;
    int maxsteps,level,rank;
    bool stop,boundreached;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // get parameters
    maxsteps  = this->m_Opt->GetMaxStepsParaCont();
    betamin   = this->m_Opt->GetBetaMinParaCont();
    betascale = this->m_Opt->GetBetaScaleParaCont();

    ierr=Assert(betascale < 1.0,"scale for beta > 1"); CHKERRQ(ierr);
    ierr=Assert(betascale > 0.0,"scale for beta <= 0.0"); CHKERRQ(ierr);
    ierr=Assert(betamin > 0.0,"lower bound for beta < 0"); CHKERRQ(ierr);
    ierr=Assert(betamin < 1.0,"lower bound for beta > 1"); CHKERRQ(ierr);

    // store initial guess
    ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);
    ierr=VecDuplicate(x,&this->m_Solution); CHKERRQ(ierr);

    // initialize parameters
    beta = 1.0;
    betastar = beta;

    levelmsg = "level ";
    ierr=Msg("starting optimization (parameter continuation)"); CHKERRQ(ierr);

    // reduce regularization parameter by one order of magnitude until
    // we hit the tolerance
    stop=false;
    level = 0;
    while(level < maxsteps){

        this->m_Opt->SetRegularizationWeight(beta);

        ss << std::scientific << std::setw(3)
            << level <<" ( betav="<<beta
            <<"; betav*="<<betastar<<" )";

        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ierr=Msg(levelmsg + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss.str( std::string() ); ss.clear();

        // solve optimization probelm for current regularization parameter
        ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);

        // get the solution vector
        ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);

        // check bounds on jacobian
        ierr=this->m_OptimizationProblem->CheckBounds(x,stop); CHKERRQ(ierr);

        if (stop) break; // if bound reached go home

        // remember regularization parameter
        betastar = beta;

        // if we got here, the solution is valid
        ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);

        // reduce beta
        beta *= betascale;

        // if regularization parameter is smaller than
        // lower bound, let's stop this
        if (beta < betamin){

            if (this->m_Opt->GetVerbosity() > 0){
                ss << std::scientific
                   <<"regularization parameter smaller than lower bound (betav="
                   <<beta<<" < " << betamin <<"=betavmin)";
                ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.str( std::string() ); ss.clear();
            }
            break;
        }

        ++level;

    } // until we hit the tolerance

    // now start binary search
    stop=false;

    // get scale for delta beta; this parameter determines how
    // accurate we solve (how many digits) with respect to the
    // order of magnitude of magnitude of the regularization
    // parameter)
    dbetascale = this->m_Opt->GetDeltaBetaScaleParaCont();
    ierr=Assert(dbetascale < 1.0,"scale for delta betav > 1"); CHKERRQ(ierr);
    ierr=Assert(dbetascale > 0.0,"scale for delta betav < 0"); CHKERRQ(ierr);

    //update beta
    dbetamin = dbetascale*betastar;
    betahat = betascale*betastar;
    dbeta = (betastar-betahat)/2.0;
    beta  = betastar-dbeta;
    ++level;

    while(!stop){

        // set regularization parameter
        this->m_Opt->SetRegularizationWeight(beta);

        // display regularization parameter to user
        ss << std::setw(3) << level <<" ( betav="<<beta<<"; betav*="<<betastar<<" )";
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ierr=Msg(levelmsg + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss.str( std::string() ); ss.clear();

        // solve optimization probelm for current regularization parameter
        ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);

        // get the solution vector
        ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);

        // check bounds on jacobian
        boundreached=false;
        ierr=this->m_OptimizationProblem->CheckBounds(x,boundreached); CHKERRQ(ierr);

        // if bound is reached, the lower bound is now beta
        // if not, beta is our new best estimate
        if (boundreached){ betahat = beta; }
        else{

            betastar = beta; // new best estimate

            // if we got here, the solution is valid
            ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);
        }

        // increase or reduce beta
        dbeta = (betastar - betahat)/2.0;
        beta  = betastar - dbeta;

        if (fabs(dbeta) < dbetamin){ stop = true; }

        ++level;
    }

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss <<std::scientific<<"estimated regularization parameter betav="<<betastar;
    ierr=Msg(ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss.str( std::string() ); ss.clear();

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief solves the optimization problem by simply reducing
 * the regularization parameter until the mapping becomes
 * non-diffeomorphic/breaches the user defined bound; stored
 * velocity field (solution) is last iterate that resulted in
 * diffeomorphic deformation map (as judged by the determinant
 * of the deformation gradient)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunRegParaReductionSearch"
PetscErrorCode Optimizer::RunRegParaReductionSearch()
{
    PetscErrorCode ierr;
    std::string levelmsg;
    std::stringstream ss;
    ScalarType beta,betamin,betastar,betascale;
    Vec x;
    int maxsteps, rank;
    bool stop;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // get parameters
    maxsteps  = this->m_Opt->GetMaxStepsParaCont();
    betamin   = this->m_Opt->GetBetaMinParaCont();
    betascale = this->m_Opt->GetBetaScaleParaCont();

    ierr=Assert(betascale < 1.0,"scale for beta > 1"); CHKERRQ(ierr);
    ierr=Assert(betascale > 0.0,"scale for beta <= 0.0"); CHKERRQ(ierr);
    ierr=Assert(betamin > 0.0,"lower bound for beta < 0"); CHKERRQ(ierr);
    ierr=Assert(betamin < 1.0,"lower bound for beta > 1"); CHKERRQ(ierr);

    levelmsg = "level ";
    ierr=Msg("starting optimization (parameter continuation)"); CHKERRQ(ierr);

    // copy solution
    ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);
    ierr=VecDuplicate(x,&this->m_Solution); CHKERRQ(ierr);

    // set initial regularization weight
    beta = 1.0;
    betastar = beta;

    // reduce regularization parameter by one order of magnitude until
    // we hit user defined tolerances (which either is a lower bound
    // on the regularization parameter or a lower bound on the
    // determinant of the deformation gradient)
    for(int i = 0; i < maxsteps; ++i){

        // set regularization weight
        this->m_Opt->SetRegularizationWeight(beta);

        // display message to user
        ss << std::scientific << std::setw(3) << i <<" (beta="<<beta<<")";
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ierr=Msg(levelmsg + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss.str( std::string() ); ss.clear();

        // solve optimization problem for current regularization parameter
        ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);

        // get solution
        ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);

        stop=false; // not yet we're not

        // check if we hit the bound
        ierr=this->m_OptimizationProblem->CheckBounds(x,stop); CHKERRQ(ierr);

        if (stop) break; // if bound reached go home

        // remember best estimate
        betastar = beta;

        // if we got here, the solution is valid
        ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);

        beta *= betascale; // reduce beta

        // if the regularization parameter is smaller than
        // the lower bound, we're done
        if (beta < betamin){

            if (this->m_Opt->GetVerbosity() > 0){
                ss <<"regularization parameter smaller than lower bound (betav="
                   <<beta<<" < " << betamin <<"=betavmin)";
                ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.str( std::string() ); ss.clear();
            }
            break;
        }

    } // parameter reduction


    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss <<std::scientific<<"estimated regularization parameter betav="<<betastar;
    ierr=Msg(ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss.str( std::string() ); ss.clear();

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief solves the optimization problem
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunRegParaContinuation"
PetscErrorCode Optimizer::RunRegParaContinuation()
{
    PetscErrorCode ierr;
    std::string levelmsg;
    std::stringstream ss;
    ScalarType beta,betastar;
    Vec x;
    int level,rank;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    levelmsg = "level ";
    ierr=Msg("starting optimization (parameter continuation)"); CHKERRQ(ierr);

    // get initial guess (best we have so far)
    ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);
    ierr=VecDuplicate(x,&this->m_Solution); CHKERRQ(ierr);

    // get target regularization weight
    betastar=this->m_Opt->GetRegularizationWeight(0);
    ierr=Assert(betastar > 0.0,"target beta < 0"); CHKERRQ(ierr);
    ierr=Assert(betastar < 1.0,"target beta > 1"); CHKERRQ(ierr);


    // reduce regularization parameter
    level = 0;
    beta=1.0; // set initial regularization weight
    while(beta > betastar){

        // set regularization weight
        this->m_Opt->SetRegularizationWeight(beta);

        // display message to user
        ss << std::scientific << std::setw(3) << level <<" (beta="<<beta<<"; beta*="<<betastar<<")";
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ierr=Msg(levelmsg + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss.str( std::string() ); ss.clear();

        // solve optimization problem for current regularization parameter
        ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);

        beta /= static_cast<ScalarType>(10); // reduce beta
        ++level;

    } // parameter reduction

    beta = betastar;

    // set regularization weight
    this->m_Opt->SetRegularizationWeight(beta);

    // display message to user
    ss << std::scientific << std::setw(3) << level <<" (beta="<<beta<<"; beta*="<<betastar<<")";
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ierr=Msg(levelmsg + ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss.str( std::string() ); ss.clear();


    // solve optimization problem for current regularization parameter
    ierr=TaoSolve(this->m_Tao); CHKERRQ(ierr);

    // get solution
    ierr=TaoGetSolutionVector(this->m_Tao,&x); CHKERRQ(ierr);
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

    ierr=Assert(x!=NULL, "input pointer is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_Tao!=NULL, "optimization object not initialized"); CHKERRQ(ierr);
    ierr=Assert(this->m_Solution!=NULL, "solution vector is null"); CHKERRQ(ierr);

    // copy the solution to input
    ierr=VecCopy(this->m_Solution,x); CHKERRQ(ierr);

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
    ierr=Assert(this->m_OptimizationProblem !=NULL,"optimizer: optimization problem not set"); CHKERRQ(ierr);

    // finalize the registration (write out all data)
    ierr=this->m_OptimizationProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // display info to user, once we're done
    ierr=TaoView(this->m_Tao,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    ierr=this->m_Opt->DisplayTimeToSolution(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of name space




#endif
