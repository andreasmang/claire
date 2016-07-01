#ifndef _OPTIMIZATIONPROBLEM_CPP_
#define _OPTIMIZATIONPROBLEM_CPP_

#include "OptimizationProblem.hpp"


namespace reg
{




/****************************************************************************
 * Function: PrecondMatVec
 * Description: computes the matrix vector product Px
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "TwoLevelPCMatVec"
PetscErrorCode TwoLevelPCMatVec(Mat P, Vec x, Vec Px)
{
    PetscErrorCode ierr;
    void* ptr;
    OptimizationProblem *optprob = NULL;

    PetscFunctionBegin;

    ierr=MatShellGetContext(P,&ptr); CHKERRQ(ierr);
    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr=optprob->TwoLevelPrecondMatVec(Px,x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/****************************************************************************
 * Function: PrecondMonitor
 * Description: monitor evolution of krylov subspace method
 *****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PrecondMonitor"
PetscErrorCode PrecondMonitor(KSP ksp,IntType it,ScalarType rnorm,void* ptr)
{
    PetscErrorCode ierr;
    (void)ksp;
    OptimizationProblem* optprob=NULL;
    std::stringstream itss, rnss;
    std::string kspmeth, msg;

    PetscFunctionBegin;

    optprob = static_cast<OptimizationProblem*>(ptr);
    ierr=Assert(optprob!=NULL,"user is null pointer"); CHKERRQ(ierr);

    kspmeth=" >> PC  "; itss << std::setw(3) << it; rnss << std::scientific << rnorm;
    msg = kspmeth +  itss.str() + "  ||r||_2 = " + rnss.str();
    ierr=DbgMsg(msg); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: OptimizationProblem
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimizationProblem"
OptimizationProblem::OptimizationProblem()
{
    this->Initialize();
}




/********************************************************************
 * Name: OptimizationProblem
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimizationProblem"
OptimizationProblem::OptimizationProblem(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * Name: OptimizationProblem
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~OptimizationProblem"
OptimizationProblem::~OptimizationProblem(void)
{
}




/********************************************************************
 * Name: Initialize
 * Description: init class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode OptimizationProblem::Initialize(void)
{
    PetscFunctionBegin;

    this->m_Opt=NULL;

    this->m_InitGradNorm=0.0;
    this->m_InitObjectiveVal=0.0;
    this->m_InitDistanceVal=0.0;

    this->m_NumOuterIter=0;
    this->m_PCMatVec=0;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SetOptions
 * Description: set the registration options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetOptions"
PetscErrorCode OptimizationProblem::SetOptions(RegOpt* opt)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // check for null pointer
    ierr=Assert(opt != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Opt != NULL){
        delete this->m_Opt;
        this->m_Opt = NULL;
    }

    // overwrite
    this->m_Opt = opt;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: DerivativeCheck
 * Description: check gradient based on a taylor expansion
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DerivativeCheck"
PetscErrorCode OptimizationProblem::DerivativeCheck()
{
    PetscErrorCode ierr=0;
    Vec v=NULL,vtilde=NULL,w=NULL,dvJ=NULL;
    PetscRandom rctx;
    IntType nl,ng;
    ScalarType h=0.0,htilde=0.0,Jv=0.0,dvJw=0.0,
                Jvtilde=0.0,e[2],normv=0.0,normw=0.0;
    char buffer[256];

    PetscFunctionBegin;

    ierr=Assert(this->m_Opt!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=DbgMsg("performing derivative check (gradient)"); CHKERRQ(ierr);
    sprintf(buffer,"%-12s %-12s %-12s","h","e(h)","e(h^2)");
    ierr=DbgMsg(buffer); CHKERRQ(ierr);

    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    // create an extra array for initial guess (has to be flat for optimizer)
    ierr=VecCreate(PETSC_COMM_WORLD,&v); CHKERRQ(ierr);
    ierr=VecSetSizes(v,3*nl,3*ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(v); CHKERRQ(ierr);
    ierr=VecSet(v,0.0); CHKERRQ(ierr);

    ierr=VecDuplicate(v,&vtilde); CHKERRQ(ierr);
    ierr=VecDuplicate(v,&dvJ); CHKERRQ(ierr);
    ierr=VecDuplicate(v,&w); CHKERRQ(ierr);

    // create random vectors
    ierr=PetscRandomCreate(PETSC_COMM_WORLD,&rctx); CHKERRQ(ierr);
    ierr=PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);
    ierr=VecSetRandom(v,rctx); CHKERRQ(ierr);
    ierr=VecSetRandom(w,rctx); CHKERRQ(ierr);

    // compute norm of random vectors
    ierr=VecNorm(v,NORM_2,&normv); CHKERRQ(ierr);
    ierr=VecNorm(w,NORM_2,&normw); CHKERRQ(ierr);

    // normalize random number
    ierr=VecScale(v,1.0/normv); CHKERRQ(ierr);
    ierr=VecScale(w,1.0/normw); CHKERRQ(ierr);

    // compute value of objective functional
    ierr=this->EvaluateObjective(&Jv,v); CHKERRQ(ierr);
    ierr=this->EvaluateGradient(dvJ,v); CHKERRQ(ierr);

    // do the derivative check
    h = 1E-8;
    for (int i = 0; i < 10; ++i){

        // compute step size
        htilde = h*pow(10.0,i);

        // update velocity field
        ierr=VecCopy(v,vtilde); CHKERRQ(ierr);

        // perturb velocity field
        ierr=VecAXPY(vtilde,htilde,w); CHKERRQ(ierr);

        // evaluate objective
        ierr=this->EvaluateObjective(&Jvtilde,vtilde); CHKERRQ(ierr);

        // inner product between perturbation and gradient
        ierr=VecTDot(w,dvJ,&dvJw); CHKERRQ(ierr);
        dvJw*=htilde;

        e[0] = (Jvtilde - Jv);
        e[1] = (Jvtilde - Jv - dvJw);

        e[0] = std::abs(e[0]);
        e[1] = std::abs(e[1]);

        sprintf(buffer,"%e %e %e",htilde,e[0],e[1]);

        ierr=DbgMsg(buffer); CHKERRQ(ierr);
    }

    // clean up
    ierr=PetscRandomDestroy(&rctx); CHKERRQ(ierr);
    ierr=VecDestroy(&v); CHKERRQ(ierr);
    ierr=VecDestroy(&w); CHKERRQ(ierr);
    ierr=VecDestroy(&vtilde); CHKERRQ(ierr);
    ierr=VecDestroy(&dvJ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianSymmetryCheck
 * Description: check symmetry of hessian
 * the idea is to use the identity
 *   \langle A x, A x \rangle = \langle A^T*Ax, x \rangle
 * for the inner product
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianSymmetryCheck"
PetscErrorCode OptimizationProblem::HessianSymmetryCheck()
{
    PetscErrorCode ierr=0;
    IntType nl,ng;
    Vec v=NULL,Hv=NULL,HHv=NULL;
    ScalarType HvHv=0.0, HHvv=0.0, symerr=0.0, relsymerr=0.0, normHv=0.0;
    std::string msg;
    PetscRandom rctx;

    PetscFunctionBegin;

    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    // create an extra array for initial guess (has to be flat for optimizer)
    ierr=VecCreate(PETSC_COMM_WORLD,&v); CHKERRQ(ierr);
    ierr=VecSetSizes(v,3*nl,3*ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(v); CHKERRQ(ierr);

    // allocate hessian mat vec
    ierr=VecDuplicate(v,&Hv); CHKERRQ(ierr);
    ierr=VecDuplicate(v,&HHv); CHKERRQ(ierr);

    // create random vectors
    ierr=PetscRandomCreate(PETSC_COMM_WORLD,&rctx); CHKERRQ(ierr);
    ierr=PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);
    ierr=VecSetRandom(v,rctx); CHKERRQ(ierr);

    // apply hessian to vector field
    ierr=this->HessianMatVec(Hv,v); CHKERRQ(ierr);
    ierr=this->HessianMatVec(HHv,Hv); CHKERRQ(ierr);

    ierr=VecTDot(Hv,Hv,&HvHv); CHKERRQ(ierr);
    ierr=VecTDot(HHv,v,&HHvv); CHKERRQ(ierr);
    ierr=VecNorm(Hv,NORM_2,&normHv); CHKERRQ(ierr);

    symerr=std::abs(HvHv-HHvv);
    relsymerr=symerr/normHv;

    std::stringstream sserr, ssrelerr;
    sserr << symerr;
    ssrelerr << relsymerr;
    msg = "symmetry error of hessian: " + sserr.str()
        + " (relative " + ssrelerr.str() + ")";

    ierr=DbgMsg(msg); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of namespace




#endif // _OPTIMIZATIONPROBLEM_CPP_

