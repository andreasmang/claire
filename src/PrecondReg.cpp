#ifndef _PRECONDREG_CPP_
#define _PRECONDREG_CPP_

#include "PrecondReg.hpp"



namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PrecondReg"
PrecondReg::PrecondReg()
{
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~PrecondReg"
PrecondReg::~PrecondReg()
{
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistration"
PrecondReg::PrecondReg(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode PrecondReg::Initialize()
{
    PetscFunctionBegin;

    this->m_Opt = NULL;
    this->m_OptCoarse = NULL;
    this->m_MatVec = NULL;
    this->m_KrylovMethod = NULL;

    this->m_OptProbCoarse = NULL;
    this->m_VelocityFieldCoarse = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode PrecondReg::ClearMemory()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_KrylovMethod != NULL){
        ierr=KSPDestroy(&this->m_KrylovMethod); CHKERRQ(ierr);
        this->m_KrylovMethod = NULL;
    }

    if (this->m_MatVec != NULL){
        ierr=MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    if (this->m_OptProbCoarse != NULL){
        delete this->m_OptProbCoarse;
        this->m_OptProbCoarse=NULL;
    }
    if (this->m_OptCoarse != NULL){
        delete this->m_OptCoarse;
        this->m_OptCoarse = NULL;
    }

    if (this->m_VelocityFieldCoarse != NULL){
        delete this->m_VelocityFieldCoarse;
        this->m_VelocityFieldCoarse=NULL;
    }
    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetProblem"
PetscErrorCode PrecondReg::SetProblem(PrecondReg::OptProbType* optprob)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_OptProb = optprob;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MatVec"
PetscErrorCode PrecondReg::MatVec(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    // switch case for choice of preconditioner
    switch(this->m_Opt->GetKrylovSolverPara().pctype){
        case NOPC:
        {
            ierr=WrngMsg("no preconditioner used"); CHKERRQ(ierr);
            ierr=VecCopy(x,Px); CHKERRQ(ierr);
            break;
        }
        case INVREG:
        {
            ierr=this->ApplyInvRegPC(Px,x); CHKERRQ(ierr);
            break;
        }
        case TWOLEVEL:
        {
            ierr=this->Apply2LevelPC(Px,x); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner not defined"); CHKERRQ(ierr);
            break;
        }
    }

    // increment counter
    this->m_Opt->IncrementCounter(PCMATVEC);



    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief apply inverse of regularization operator as preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvRegPC"
PetscErrorCode PrecondReg::ApplyInvRegPC(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    // check if optimization problem is set up
    ierr=Assert(this->m_OptProb!=NULL,"null pointer"); CHKERRQ(ierr);

    // start timer
    ierr=this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);

    // apply inverse regularization operator
    ierr=this->m_OptProb->ApplyInvRegOp(Px,x); CHKERRQ(ierr);

    // stop timer
    ierr=this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Apply2LevelPC"
PetscErrorCode PrecondReg::Apply2LevelPC(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    // check if optimization problem is set up
    ierr=Assert(this->m_OptProb!=NULL,"null pointer"); CHKERRQ(ierr);

    // setup preconditioner
    ierr=this->SetUp2LevelPC(); CHKERRQ(ierr);

    // do setup
    if (this->m_KrylovMethod == NULL){
        ierr=this->SetupKrylovMethod(); CHKERRQ(ierr);
    }
    ierr=this->SetTolerancesKrylovMethod(); CHKERRQ(ierr);

    // start timer
    ierr=this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);


    ierr=KSPSolve(this->m_KrylovMethod,x,Px); CHKERRQ(ierr);


    // stop timer
    ierr=this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief setup 2 level preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetUp2LevelPC"
PetscErrorCode PrecondReg::SetUp2LevelPC()
{
    PetscErrorCode ierr;
    Vec m=NULL,lambda=NULL;
    PetscFunctionBegin;

    // check if optimization problem is set up
    ierr=Assert(this->m_OptProb!=NULL,"null pointer"); CHKERRQ(ierr);

    // start timer
    ierr=this->m_Opt->StartTimer(PMVSETUP); CHKERRQ(ierr);

    if (this->m_OptProbCoarse==NULL){

        try{ this->m_OptCoarse = new RegOpt(*this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=this->m_OptCoarse->DoSetup(false); CHKERRQ(ierr);

        // allocate class for registration
        if (this->m_Opt->GetRegModel() == COMPRESSIBLE){
            try{ this->m_OptProbCoarse = new OptimalControlRegistration(this->m_OptCoarse); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else if (this->m_Opt->GetRegModel() == STOKES){
            try{ this->m_OptProbCoarse = new OptimalControlRegistrationIC(this->m_OptCoarse); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else if (this->m_Opt->GetRegModel() == RELAXEDSTOKES){
            try{ this->m_OptProbCoarse = new OptimalControlRegistrationIC(this->m_OptCoarse); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else{
            try{ this->m_OptProbCoarse = new OptimalControlRegistration(this->m_OptCoarse); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

    }

    // allocate velocity field
    if(this->m_VelocityFieldCoarse==NULL){
        try{ this->m_VelocityFieldCoarse = new VecField(this->m_OptCoarse); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr=this->m_OptProb->GetControlVariable(this->m_VelocityFieldCoarse); CHKERRQ(ierr);
    ierr=this->m_OptProb->GetStateVariable(m); CHKERRQ(ierr);
    ierr=this->m_OptProb->GetAdjointVariable(lambda); CHKERRQ(ierr);

    ierr=Assert(m!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(lambda!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_VelocityFieldCoarse!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->m_OptProbCoarse->SetControlVariable(this->m_VelocityFieldCoarse); CHKERRQ(ierr);
    ierr=this->m_OptProbCoarse->SetStateVariable(m); CHKERRQ(ierr);
    ierr=this->m_OptProbCoarse->SetAdjointVariable(lambda); CHKERRQ(ierr);

    // stop timer
    ierr=this->m_Opt->StopTimer(PMVSETUP); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief do setup for krylov method
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupKrylovMethod"
PetscErrorCode PrecondReg::SetupKrylovMethod()
{
    PetscErrorCode ierr;
    PC pc=NULL;
    IntType nl,ng;

    PetscFunctionBegin;

    // get sizes
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;


    // create krylov method
    if (this->m_KrylovMethod != NULL){
        ierr=KSPDestroy(&this->m_KrylovMethod); CHKERRQ(ierr);
        this->m_KrylovMethod = NULL;
    }
    ierr=KSPCreate(PETSC_COMM_WORLD,&this->m_KrylovMethod); CHKERRQ(ierr);


    switch (this->m_Opt->GetKrylovSolverPara().pcsolver){
        case CHEB:
        {
            // chebyshev iteration
            ierr=KSPSetType(this->m_KrylovMethod,KSPCHEBYSHEV); CHKERRQ(ierr);
            break;
        }
        case PCG:
        {
            // preconditioned conjugate gradient
            ierr=KSPSetType(this->m_KrylovMethod,KSPCG); CHKERRQ(ierr);
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient
            ierr=KSPSetType(this->m_KrylovMethod,KSPFCG); CHKERRQ(ierr);
            break;
        }
        case GMRES:
        {
            // GMRES
            ierr=KSPSetType(this->m_KrylovMethod,KSPGMRES); CHKERRQ(ierr);
            break;
        }
        case FGMRES:
        {
            // flexible GMRES
            ierr=KSPSetType(this->m_KrylovMethod,KSPFGMRES); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }

    //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
    //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
    //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
    ierr=KSPSetNormType(this->m_KrylovMethod,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    ierr=KSPSetInitialGuessNonzero(this->m_KrylovMethod,PETSC_TRUE); CHKERRQ(ierr);

    // set up matvec for preconditioner
    if (this->m_MatVec != NULL){
        ierr=MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    ierr=MatCreateShell(PETSC_COMM_WORLD,3*nl,3*nl,3*ng,3*ng,this,&this->m_MatVec); CHKERRQ(ierr);
    ierr=MatShellSetOperation(this->m_MatVec,MATOP_MULT,(void(*)(void))InvertPrecondMatVec); CHKERRQ(ierr);
    ierr=KSPSetOperators(this->m_KrylovMethod,this->m_MatVec,this->m_MatVec);CHKERRQ(ierr);
    ierr=KSPMonitorSet(this->m_KrylovMethod,InvertPrecondKrylovMonitor,this,PETSC_NULL); CHKERRQ(ierr);

    // remove preconditioner
    ierr=KSPGetPC(this->m_KrylovMethod,&pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCNONE); CHKERRQ(ierr); ///< set no preconditioner

    // finish
    ierr=KSPSetUp(this->m_KrylovMethod); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set the tolerances for krylov method
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetTolerancesKrylovMethod"
PetscErrorCode PrecondReg::SetTolerancesKrylovMethod()
{
    PetscErrorCode ierr;
    ScalarType reltol,abstol,divtol,scale;
    IntType maxit;

    PetscFunctionBegin;

    divtol = 1E+06;
    abstol = 1E-16;
    reltol = 1E-16;
    maxit  = 1000;

    // check for null pointer
    ierr=Assert(this->m_KrylovMethod!=NULL,"null pointer"); CHKERRQ(ierr);

    scale = this->m_Opt->GetKrylovSolverPara().pcsolvertol;

    switch (this->m_Opt->GetKrylovSolverPara().pcsolver){
        case CHEB:
        {
            // chebyshev iteration
            maxit = this->m_Opt->GetKrylovSolverPara().pcsolvermaxit;
            break;
        }
        case PCG:
        {
            // preconditioned conjugate gradient
            reltol = scale*this->m_Opt->GetKrylovSolverPara().reltol;
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient
            maxit = this->m_Opt->GetKrylovSolverPara().pcsolvermaxit;
            break;
        }
        case GMRES:
        {
            // GMRES
            reltol = scale*this->m_Opt->GetKrylovSolverPara().reltol;
            break;
        }
        case FGMRES:
        {
            // flexible GMRES
            maxit = this->m_Opt->GetKrylovSolverPara().pcsolvermaxit;
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }
    reltol = std::max(reltol,1E-16); // make sure tolerance is non-zero
    reltol = std::min(reltol,0.5); // make sure tolerance smaller than 0.25

    // set tolerances
    ierr=KSPSetTolerances(this->m_KrylovMethod,reltol,abstol,divtol,maxit); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief do setup for two level preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode PrecondReg::HessianMatVec(Vec Hx, Vec x)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // check if optimization problem is set up
    ierr=Assert(this->m_OptProbCoarse!=NULL,"null pointer"); CHKERRQ(ierr);

    // apply hessian (hessian matvec)
    ierr=this->m_OptProbCoarse->HessianMatVec(Hx,x); CHKERRQ(ierr);
    //ierr=this->m_OptProb->HessianMatVec(Hx,x); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




} // end of namespace




#endif // _PRECONDREG_CPP_
