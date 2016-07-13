#ifndef _KRYLOVINTERFACEREG_CPP_
#define _KRYLOVINTERFACEREG_CPP_

#include "KrylovInterfaceReg.hpp"



namespace reg
{




/****************************************************************************
 * @brief monitor evolution of krylov subspace method
 *****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "KrylovMonitor"
PetscErrorCode KrylovMonitor(KSP krylovmethod,IntType it,ScalarType rnorm,void* ptr)
{
    PetscErrorCode ierr;
    (void)krylovmethod;
    KSPConvergedReason reason;
    OptimizationProblem* optprob;
    std::stringstream itss, rnss;
    std::string kspmeth, msg;

    PetscFunctionBegin;

    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    kspmeth = optprob->GetOptions()->GetKrylovSolverPara().name;
    itss << std::setw(5) << it;
    rnss << std::scientific << rnorm;
    msg = kspmeth + "   " + itss.str() + "  ||r||_2 = " + rnss.str();
    ierr=DbgMsg(msg); CHKERRQ(ierr);

    ierr=KSPGetConvergedReason(krylovmethod,&reason); CHKERRQ(ierr);
    ierr=DispKSPConvReason(reason); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/****************************************************************************
 * @brief preprocess right hand side and initial condition before entering
 * the krylov subspace method; in the context of numerical optimization this
 * means we preprocess the gradient and the incremental control variable
 * @para krylovmethod pointer to krylov method
 * @para b right hand side of equation
 * @para x solution vector
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PreKrylovSolve"
PetscErrorCode PreKrylovSolve(KSP krylovmethod,Vec b, Vec x,void* ptr)
{
    PetscErrorCode ierr;
    ScalarType gnorm=0.0,g0norm=0.0,reltol,abstol,divtol,
                uppergradbound,lowergradbound;
    IntType maxit;
    std::stringstream ss;
    std::string msg;
    OptimizationProblem* optprob=NULL;

    PetscFunctionBegin;

    uppergradbound=0.5;
    lowergradbound=1E-12;

    (void)krylovmethod;

    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    if (optprob->GetOptions()->GetHessianMatVecType() == PRECONDMATVEC){
        // get current gradient and compute norm
        // before we apply the preconditioner to
        // the right hand side
        ierr=VecNorm(b,NORM_2,&gnorm); CHKERRQ(ierr);
    }

    // apply pre processing to right hand side and initial condition
    ierr=optprob->PreKrylovSolve(b,x); CHKERRQ(ierr);

    // user forcing sequence to estimate adequate tolerance
    // for solution of KKT system (Eisenstat-Walker)
    if(optprob->GetOptions()->GetKrylovSolverPara().fseqtype != NOFS){

        if (gnorm == 0.0){
            // get current gradient and compute norm
            ierr=VecNorm(b,NORM_2,&gnorm); CHKERRQ(ierr);
        }

        if (!optprob->GetOptions()->GetKrylovSolverPara().g0normset){

            if(optprob->GetOptions()->GetVerbosity() > 1){
                ss << std::fixed << std::scientific << gnorm;
                msg = optprob->GetOptions()->GetKrylovSolverPara().name +
                    ": setting initial ||g|| in krylov method (" + ss.str() + ")";
                ierr=DbgMsg(msg); CHKERRQ(ierr);
                ss.str(std::string()); ss.clear();
            }

            optprob->GetOptions()->SetInitialGradNormKrylovMethod(gnorm);

        }

        // get initial value for gradient
        g0norm = optprob->GetOptions()->GetKrylovSolverPara().g0norm;
        ierr=Assert(g0norm > 0.0,"initial gradient is zero"); CHKERRQ(ierr);

        // normalize
        gnorm /= g0norm;

        // get current tolerances
        ierr=KSPGetTolerances(krylovmethod,&reltol,&abstol,&divtol,&maxit); CHKERRQ(ierr);

        if(optprob->GetOptions()->GetKrylovSolverPara().fseqtype == QDFS){
            // assuming quadratic convergence (we do not solver more
            // accurately than 12 digits)
            reltol=PetscMax(lowergradbound,PetscMin(uppergradbound,gnorm));
        }
        else{
            // assuming superlinear convergence (we do not solver
            // more accurately than 12 digitis)
            reltol=PetscMax(lowergradbound,PetscMin(uppergradbound,std::sqrt(gnorm)));
        }

        // overwrite tolerances with estimate
        ierr=KSPSetTolerances(krylovmethod,reltol,abstol,divtol,maxit); CHKERRQ(ierr);
    }

    // pass tolerance to optimization problem (for preconditioner)
    optprob->SetRelTolKrylovMethod(reltol);

    if(optprob->GetOptions()->GetVerbosity() > 0){
        ss << std::fixed << std::scientific << reltol;
        msg = optprob->GetOptions()->GetKrylovSolverPara().name +
              ": computing solution of hessian system (tol=" + ss.str() + ")";
        ierr=DbgMsg(msg); CHKERRQ(ierr);
    }

    // check symmetry of hessian
    #ifdef _REG_DEBUG_
    ierr=optprob->HessianSymmetryCheck(); CHKERRQ(ierr);
    #endif// _REG_DEBUG_

    PetscFunctionReturn(0);
}




/****************************************************************************
 * @brief postprocess right hand side and initial condition before entering
 * the krylov subspace method; in the context of numerical optimization this
 * means we postprocess the gradient and the incremental control variable
 * @para krylovmethod pointer to krylov method
 * @para b right hand side of equation
 * @para x solution vector
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PostKrylovSolve"
PetscErrorCode PostKrylovSolve(KSP krylovmethod,Vec b, Vec x,void* ptr)
{
    PetscErrorCode ierr;
    OptimizationProblem* optprob=NULL;
    KSPConvergedReason reason;
    std::string convmsg;

    PetscFunctionBegin;

    (void)krylovmethod;

    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr=optprob->PostKrylovSolve(b,x); CHKERRQ(ierr);

    if(optprob->GetOptions()->GetVerbosity() > 0){
        ierr=KSPGetConvergedReason(krylovmethod,&reason);
        ierr=DispKSPConvReason(reason); CHKERRQ(ierr);
    }


    PetscFunctionReturn(0);
}




/****************************************************************************
 * @briefdisplay the convergence reason of the KSP method
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DispKSPConvReason"
PetscErrorCode DispKSPConvReason(KSPConvergedReason flag)
{
    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    switch(flag){
        case KSP_CONVERGED_RTOL_NORMAL:
        {
            msg="krylov method converged: ||r||_2 < tol ||b||_2";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ATOL_NORMAL:
        {
            msg="krylov method converged: ||r||_2 < tol";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_RTOL:
        {
            msg="krylov method converged: ||r||_2 < tol ||b||_2";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ATOL:
        {
            msg="krylov method converged: ||r||_2 < tol";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ITS:
        {
            //used by the KSPPREONLY solver after the single iteration of
            //the preconditioner is applied
            msg="krylov method converged: k > maxit";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_CG_NEG_CURVE:
        {
            msg="krylov method: negative curvature detected";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_CG_CONSTRAINED:
        {
            msg="krylov method: convergence is reached";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_STEP_LENGTH:
        {
            msg="krylov method: converged step length";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_HAPPY_BREAKDOWN:
        {
            msg="krylov method: converged happy break down";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NULL:
        {
            msg="krylov method: divergence detected";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_ITS:
        {
            msg="krylov method: max number of iterations reached";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_DTOL:
        {
            msg="krylov method: divergence detected (||r||_2 increased by a factor of divtol)";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_BREAKDOWN:
        {
            //breakdown in Krylov method was detected
            //method could not continue to enlarge Krylov subspace;
            //could be due to a singlular matrix or preconditioner
            msg="krylov method: generic breakdown (potential singular operator)";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_BREAKDOWN_BICG:
        {
            msg="krylov method: initial ||r||_2 is orthogonal to preconditioned r";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NONSYMMETRIC:
        {
            msg="krylov method: operators (A or P) are not symmetric";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_INDEFINITE_PC:
        {
            msg="krylov method: preconditioner is indefinite";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NANORINF:
        {
            msg="krylov method: ||r||_2 is NAN or INF";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_INDEFINITE_MAT:
        {
            msg="krylov method: A is indefinite";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ITERATING:
        {
            // don't display anaything and continue iterating
            break;
        }
        default:
        {
            msg="krylov method: convergence reason not defined";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }

    }


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief monitor for evolution of krylov subspace method for the
 * inversion of the preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "InvertPrecondKrylovMonitor"
PetscErrorCode InvertPrecondKrylovMonitor(KSP krylovmethod,IntType it,ScalarType rnorm,void* ptr)
{
    PetscErrorCode ierr;
    (void)krylovmethod;
    KSPConvergedReason reason;
    PrecondReg *precond = NULL;
    std::stringstream itss, rnss;
    std::string kspmeth, msg;

    PetscFunctionBegin;

    precond = (PrecondReg*)ptr;
    ierr=Assert(precond!=NULL,"null pointer"); CHKERRQ(ierr);

    kspmeth = " >> PC " + precond->GetOptions()->GetKrylovSolverPara().pcname;
    itss << std::setw(5) << it;
    rnss << std::scientific << rnorm;
    msg = kspmeth +  itss.str() + "  ||r||_2 = " + rnss.str();
    ierr=DbgMsg(msg); CHKERRQ(ierr);

    ierr=KSPGetConvergedReason(krylovmethod,&reason); CHKERRQ(ierr);
    ierr=DispKSPConvReason(reason); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/****************************************************************************
 * @brief computes the hessian matrix vector product Hx = H*xtilde
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "InvertPrecondMatVec"
PetscErrorCode InvertPrecondMatVec(Mat P, Vec x, Vec Px)
{
    PetscErrorCode ierr;
    void* ptr;
    PrecondReg *preconditioner = NULL;

    PetscFunctionBegin;

    ierr=MatShellGetContext(P,(void**)&ptr); CHKERRQ(ierr);
    preconditioner = (PrecondReg*)ptr;
    ierr=Assert(preconditioner!=NULL,"null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr=preconditioner->HessianMatVec(Px,x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/****************************************************************************
 * @brief
 * @para krylovmethod pointer to krylov method
 * @para b right hand side of equation
 * @para x solution vector
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "InvertPrecondPreKrylovSolve"
PetscErrorCode InvertPrecondPreKrylovSolve(KSP krylovmethod,Vec b, Vec x,void* ptr)
{
    PetscErrorCode ierr;
    PrecondReg* precond=NULL;
    IntType maxits;
    ScalarType reltol,abstol,divtol,scale,eigmax,eigmin;
    std::stringstream itss,rnss;
    std::string msg;

    PetscFunctionBegin;

    precond = (PrecondReg*)ptr;
    ierr=Assert(precond!=NULL,"null pointer"); CHKERRQ(ierr);

    // setup preconditioner
    if (!precond->GetOptions()->GetKrylovSolverPara().pcsetupdone){
        ierr=precond->DoSetup(); CHKERRQ(ierr);
    }

    // set the default values
    maxits = 1E3;
    reltol = precond->GetOptions()->GetKrylovSolverPara().pctol[0];
    abstol = precond->GetOptions()->GetKrylovSolverPara().pctol[1];
    divtol = precond->GetOptions()->GetKrylovSolverPara().pctol[2];
    scale  = precond->GetOptions()->GetKrylovSolverPara().pctolscale;

    switch (precond->GetOptions()->GetKrylovSolverPara().pcsolver){
        case CHEB:
        {
            // chebyshev iteration
            maxits = precond->GetOptions()->GetKrylovSolverPara().pcmaxit;
            ierr=precond->EstimateEigenValues(); CHKERRQ(ierr);
            break;
        }
        case PCG:
        {
            // preconditioned conjugate gradient
            reltol = scale*precond->GetOptions()->GetKrylovSolverPara().reltol;
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient
            maxits = precond->GetOptions()->GetKrylovSolverPara().pcmaxit;
            break;
        }
        case GMRES:
        {
            // GMRES
            reltol = scale*precond->GetOptions()->GetKrylovSolverPara().reltol;
            break;
        }
        case FGMRES:
        {
            // flexible GMRES
            maxits = precond->GetOptions()->GetKrylovSolverPara().pcmaxit;
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }
    reltol = std::max( reltol, 1E-16 ); // make sure tolerance is non-zero
    reltol = std::min( reltol, 5E-1  ); // make sure tolerance smaller than 0.25

    // set tolerances
    ierr=KSPSetTolerances(krylovmethod,reltol,abstol,divtol,maxits); CHKERRQ(ierr);


    if (precond->GetOptions()->GetVerbosity() > 1){

        ierr=KSPGetTolerances(krylovmethod,&reltol,NULL,NULL,&maxits); CHKERRQ(ierr);

        rnss << std::fixed << std::scientific << reltol;
        itss << maxits;

        msg = precond->GetOptions()->GetKrylovSolverPara().pcname
            + ": inverting preconditioner (tol=" + rnss.str()
            + ", maxit=" + itss.str() + ")";
        ierr=DbgMsg(msg); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}




}


#endif // _KRYLOVINTERFACEREG_CPP_
