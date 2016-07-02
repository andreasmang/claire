#ifndef _TWOLEVELPRECONDREG_CPP_
#define _TWOLEVELPRECONDREG_CPP_

#include "TwoLevelPrecondReg.hpp"

namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "TwoLevelPrecondReg"
TwoLevelPrecondReg::TwoLevelPrecondReg()
{
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~TwoLevelPrecondReg"
TwoLevelPrecondReg::~TwoLevelPrecondReg()
{
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistration"
TwoLevelPrecondReg::TwoLevelPrecondReg(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode TwoLevelPrecondReg::Initialize()
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
PetscErrorCode TwoLevelPrecondReg::ClearMemory()
{
    //PetscErrorCode ierr;
    PetscFunctionBegin;


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Apply"
PetscErrorCode TwoLevelPrecondReg::Apply(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr=this->DoSetup(); CHKERRQ(ierr);

    ierr=KSPSolve(this->m_PCKSP,x,Px); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief estimate hessian eigenvalues
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EstimateEigenValues"
PetscErrorCode TwoLevelPrecondReg::EstimateEigenValues()
{
    PetscErrorCode ierr;
    ScalarType emax,emin;
    PetscFunctionBegin;

    emax = 1.0;
    emin = 0.1;

    ierr=KSPChebyshevSetEigenvalues(this->m_PCKSP,emax,emin); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DoSetup"
PetscErrorCode TwoLevelPrecondReg::DoSetup()
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

    if (this->m_PCKSP == NULL){
        ierr=KSPCreate(PETSC_COMM_WORLD,&this->m_PCKSP); CHKERRQ(ierr);
    }

    switch (this->m_Opt->GetKrylovSolverPara().pcsolver){
        case CHEB:
        {
            // chebyshev iteration
            ierr=KSPSetType(this->m_PCKSP,KSPCHEBYSHEV); CHKERRQ(ierr);
            maxit  = 10;
            break;
        }
        case PCG:
        {
            // preconditioned conjugate gradient
            ierr=KSPSetType(this->m_PCKSP,KSPCG); CHKERRQ(ierr);
            reltol = this->m_Opt->GetKrylovSolverPara().reltol;
            reltol *= 1E-1;
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient
            ierr=KSPSetType(this->m_PCKSP,KSPFCG); CHKERRQ(ierr);
            maxit  = 10;
            break;
        }
        case GMRES:
        {
            // GMRES
            ierr=KSPSetType(this->m_PCKSP,KSPGMRES); CHKERRQ(ierr);
            reltol = this->m_Opt->GetKrylovSolverPara().reltol;
            reltol *= 1E-1;
            break;
        }
        case FGMRES:
        {
            // flexible GMRES
            ierr=KSPSetType(this->m_PCKSP,KSPFGMRES); CHKERRQ(ierr);
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


    ierr=KSPSetTolerances(this->m_PCKSP,reltol,abstol,divtol,maxit); CHKERRQ(ierr);

    //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
    //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
    //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
    ierr=KSPSetNormType(this->m_PCKSP,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    ierr=KSPSetInitialGuessNonzero(this->m_PCKSP,PETSC_TRUE); CHKERRQ(ierr);

    // set up matvec for preconditioner
    if (this->m_PCMatVec != NULL){
        ierr=MatDestroy(&this->m_PCMatVec); CHKERRQ(ierr);
        this->m_PCMatVec = NULL;
     }

    ierr=MatCreateShell(PETSC_COMM_WORLD,3*nl,3*nl,3*ng,3*ng,this,&this->m_PCMatVec); CHKERRQ(ierr);
    ierr=MatShellSetOperation(this->m_PCMatVec,MATOP_MULT,(void(*)(void))TwoLevelPCMatVec); CHKERRQ(ierr);

    // set operator
    ierr=KSPSetOperators(this->m_PCKSP,this->m_PCMatVec,this->m_PCMatVec);CHKERRQ(ierr);
    ierr=KSPMonitorSet(this->m_PCKSP,PrecondMonitor,this,PETSC_NULL); CHKERRQ(ierr);
    // remove preconditioner
    ierr=KSPGetPC(this->m_PCKSP,&pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCNONE); CHKERRQ(ierr); ///< set no preconditioner

    // finish
    ierr=KSPSetUp(this->m_PCKSP); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MatVec"
PetscErrorCode TwoLevelPrecondReg::MatVec(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr=VecCopy(x,Px); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of namespace


#endif
