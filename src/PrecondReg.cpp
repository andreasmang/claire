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

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode PrecondReg::ClearMemory()
{
    //PetscErrorCode ierr;
    PetscFunctionBegin;


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
    this->m_OptimizationProblem = optprob;

    PetscFunctionReturn(0);
}



/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
/*#undef __FUNCT__
#define __FUNCT__ "Apply"
PetscErrorCode PrecondReg::Apply(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr=this->DoSetup(); CHKERRQ(ierr);

    ierr=KSPSolve(this->m_KrylovMethod,x,Px); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
*/




/********************************************************************
 * @brief
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DoSetup"
PetscErrorCode PrecondReg::DoSetup()
{
    PetscErrorCode ierr;
    PC pc=NULL;
    ScalarType reltol,abstol,divtol;
    IntType maxit;
    IntType nl,ng;

    PetscFunctionBegin;

    divtol = 1E+06;
    abstol = 1E-16;
    reltol = 1E-16;
    maxit  = 1000;

    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    if (this->m_KrylovMethod == NULL){
        ierr=KSPCreate(PETSC_COMM_WORLD,&this->m_KrylovMethod); CHKERRQ(ierr);
    }

    switch (this->m_Opt->GetKrylovSolverPara().pcsolver){
        case CHEB:
        {
            // chebyshev iteration
            ierr=KSPSetType(this->m_KrylovMethod,KSPCHEBYSHEV); CHKERRQ(ierr);
            maxit  = 10;
            break;
        }
        case PCG:
        {
            // preconditioned conjugate gradient
            ierr=KSPSetType(this->m_KrylovMethod,KSPCG); CHKERRQ(ierr);
            reltol = this->m_Opt->GetKrylovSolverPara().reltol;
            reltol *= 1E-1;
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient
            ierr=KSPSetType(this->m_KrylovMethod,KSPFCG); CHKERRQ(ierr);
            maxit  = 10;
            break;
        }
        case GMRES:
        {
            // GMRES
            ierr=KSPSetType(this->m_KrylovMethod,KSPGMRES); CHKERRQ(ierr);
            reltol = this->m_Opt->GetKrylovSolverPara().reltol;
            reltol *= 1E-1;
            break;
        }
        case FGMRES:
        {
            // flexible GMRES
            ierr=KSPSetType(this->m_KrylovMethod,KSPFGMRES); CHKERRQ(ierr);
            maxit  = 10;
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }
    reltol = std::max(reltol,1E-16); // make sure tolerance is non-zero
    reltol = std::min(reltol,0.25); // make sure tolerance smaller than 0.25


    ierr=KSPSetTolerances(this->m_KrylovMethod,reltol,abstol,divtol,maxit); CHKERRQ(ierr);

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
    ierr=MatShellSetOperation(this->m_MatVec,MATOP_MULT,(void(*)(void))TwoLevelPCMatVec); CHKERRQ(ierr);

    // set operator
    ierr=KSPSetOperators(this->m_KrylovMethod,this->m_MatVec,this->m_MatVec);CHKERRQ(ierr);
    ierr=KSPMonitorSet(this->m_KrylovMethod,PrecondMonitor,this,PETSC_NULL); CHKERRQ(ierr);
    // remove preconditioner
    ierr=KSPGetPC(this->m_KrylovMethod,&pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCNONE); CHKERRQ(ierr); ///< set no preconditioner

    // finish
    ierr=KSPSetUp(this->m_KrylovMethod); CHKERRQ(ierr);

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
            //ierr=this->Apply2LevelPrecond(Px,x); CHKERRQ(ierr);
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
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvRegPC"
PetscErrorCode PrecondReg::ApplyInvRegPC(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr=Assert(this->m_OptimizationProblem!=NULL,"null pointer"); CHKERRQ(ierr);

    // start timer
    ierr=this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);

    // apply inverse regularization operator
    ierr=this->m_OptimizationProblem->ApplyInvRegOp(Px,x); CHKERRQ(ierr);

    // stop timer
    ierr=this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


} // end of namespace


#endif // _PRECONDREG_CPP_
