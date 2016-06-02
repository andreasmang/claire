#ifndef _TAOINTERFACEREGISTRATION_CPP
#define _TAOINTERFACEREGISTRATION_CPP_


#include "TaoInterfaceRegistration.hpp"
#include "OptimizationProblem.hpp"

namespace reg
{

/****************************************************************************
 * Function: EvaluateObjective
 * Description: general purpose function to evaluate objective
 * input:
 *  tao    pointer to tao solver
 *  x      iterate x objective is to be evaluated at
 *  Jx     objective value J(x)
 *  ptr    pointer to actual optimziation problem (has to be
 *         implemented by the user)
 * output:
 *  Jx     objective value J(x)
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateObjective"
PetscErrorCode EvaluateObjective(Tao tao,Vec x,ScalarType* Jx,void* ptr)
{
    PetscErrorCode ierr;
    OptimizationProblem *optprob = NULL;

    PetscFunctionBegin;

    (void)tao;
    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    // compute objective value
    ierr=optprob->EvaluateObjective(Jx,x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/****************************************************************************
 * Function: EvaluateGradient
 * Description: general purpose function to evaluate the gradient
 * input:
 *  tao    pointer to tao solver
 *  x      iterate x gradient is to be evaluated at
 *  gx     gradient of objective functional, i.e. g(x)
 *  ptr    pointer to actual optimziation problem
 * output:
 *  gx     gradient of objective functional
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradient"
PetscErrorCode EvaluateGradient(Tao tao, Vec x, Vec gx, void* ptr)
{
    PetscErrorCode ierr;
    OptimizationProblem* optprob = NULL;

    PetscFunctionBegin;

    (void)tao;
    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"user is null pointer"); CHKERRQ(ierr);

    // compute gradient
    ierr=optprob->EvaluateGradient(gx,x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



/****************************************************************************
 * Function: EvaluateObjectiveGradient
 * Description: general purpose function to evaluate the ogjective
 * and the gradient
 * input:
 *  tao     pointer to tao solver
 *  x       iterate x gradient is to be evaluated at
 *  Jx      objective value J(x)
 *  gx      gradient of objective functional, i.e. g(x)
 *  ptr     pointer to actual optimziation problem
 * output:
 *  gx      gradient of objective functional
 *  Jx      objective functional
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateObjectiveGradient"
PetscErrorCode EvaluateObjectiveGradient(Tao tao, Vec x, ScalarType* Jx, Vec gx, void* ptr)
{
    PetscErrorCode ierr;
    OptimizationProblem* optprob = NULL;
    (void)tao;

    PetscFunctionBegin;

    // cast pointer
    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"user is null pointer"); CHKERRQ(ierr);

    // evaluate objective and gradient
    ierr=optprob->EvaluateObjective(Jx,x); CHKERRQ(ierr);
    ierr=optprob->EvaluateGradient(gx,x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/****************************************************************************
 * Function: ConstructHessian
 * Description: general purpose function to compute the hessian matrix;
 * has to be se in TaoSetHessianRoutine; this is a general purpose template
 * for matrix based algorithms; if your code is matrix free, use the
 * evaluate hessian function below!
 *
 * input:
 *  tao     pointer to tao solver
 *  x       iterate x hessian is to be evaluated at
 *  ptr     pointer to actual optimziation problem (has to be
 *          implemented by the user
 * output:
 *  H       hessian matrix
 *  P       preconditioner
 *  flag    flag used to set the Hessian matrix and linear
 *          solver in the KSP routine
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ConstructHessian"
PetscErrorCode ConstructHessian(Tao tao,Vec x,Mat* H,Mat* P,MatStructure* flag,void* ptr)
{
    PetscErrorCode ierr;
    (void)tao; (void)x; (void)H; (void)P; (void)flag; (void)ptr; (void)ierr;

    PetscFunctionBegin;

    /* DO NOTHING CURRENTLY: TODO */

    PetscFunctionReturn(0);
}




/****************************************************************************
 * Function: evaluate hessian
 * Description: general purpose function to apply the hessian
 * input:
 *  tao      pointer to tao solver
 *  x        iterate x objective is to be evaluated at
 *  ptr      pointer to actual optimziation problem (has to be
 *              implemented by the user
 * output:
 *  H        hessian matrix (can be matrix free)
 *  P        preconditioner (can be matrix free)
 *  flag     flag used to set the Hessian matrix and linear
 *           solver in the KSP routine
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateHessian"
PetscErrorCode EvaluateHessian(Tao tao,Vec x,Mat H,Mat Hpre,void* ptr)
{
    PetscErrorCode ierr;
    (void)Hpre; (void)x;
    OptimizationProblem *optprob = NULL;
    std::string msg;
    Vec dJ;
    KSP ksp;
    ScalarType gnorm,gnorm0,reltol,abstol,divtol,uppergradbound,lowergradbound;
    IntType maxit;

    PetscFunctionBegin;

    uppergradbound=0.5;
    lowergradbound=1E-12;

    // get solver context
    ierr=MatShellGetContext(H,(void**)&ptr); CHKERRQ(ierr);
    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    // get krylov subspace object
    ierr=TaoGetKSP(tao,&ksp); CHKERRQ(ierr);

    // get current tolerances
    ierr=KSPGetTolerances(ksp,&reltol,&abstol,&divtol,&maxit); CHKERRQ(ierr);

    // user forcing sequence to estimate adequate tolerance
    // for solution of KKT system (Eisenstat-Walker)
    if(optprob->GetOptions()->GetFSeqType() != NOFS){

        // get initial value for gradient
        gnorm0 = optprob->GetInitialGradNorm();
        ierr=Assert(gnorm0 > 0.0, "initial gradient is zero"); CHKERRQ(ierr);

        // get current gradient and compute norm
        ierr=TaoGetGradientVector(tao,&dJ); CHKERRQ(ierr);
        ierr=VecNorm(dJ,NORM_2,&gnorm); CHKERRQ(ierr);

        gnorm /= gnorm0;

        if(optprob->GetOptions()->GetFSeqType() == QDFS){
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
        ierr=KSPSetTolerances(ksp,reltol,abstol,divtol,maxit); CHKERRQ(ierr);
    }

    // pass tolerance to optimization problem (for preconditioner)
    optprob->SetKSPTolerance(reltol);

    if(optprob->GetOptions()->GetVerbosity() >= 2){
        std::stringstream ss;
        ss << std::fixed <<std::scientific << reltol;
        msg = "computing solution of KKT system (tol=" + ss.str() + ")";
        ierr=DbgMsg(msg); CHKERRQ(ierr);
    }

    // check symmetry of hessian
//    ierr=optprob->HessianSymmetryCheck(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



/****************************************************************************
 * Function: HessianMatVec
 * Description: computes the hessian matrix vector product Hx = H*xtilde
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode HessianMatVec(Mat H, Vec x, Vec Hx)
{
    PetscErrorCode ierr;
    void* ptr;
    OptimizationProblem *optprob = NULL;

    PetscFunctionBegin;

    ierr=MatShellGetContext(H,(void**)&ptr); CHKERRQ(ierr);
    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr=optprob->HessianMatVec(Hx,x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



/****************************************************************************
 * Function: PrecondMatVec
 * Description: computes the matrix vector product Px
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PrecondMatVec"
PetscErrorCode PrecondMatVec(PC Hpre, Vec x, Vec Hprex)
{
    PetscErrorCode ierr;
    void* ptr;
    OptimizationProblem *optprob = NULL;

    PetscFunctionBegin;

    ierr=PCShellGetContext(Hpre,&ptr); CHKERRQ(ierr);
    optprob = (OptimizationProblem*)ptr;
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr=optprob->PrecondMatVec(Hprex,x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



/****************************************************************************
 * Function: PrecondSetup
 * Description: setup the preconditioner
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PrecondSetup"
PetscErrorCode PrecondSetup(PC Hpre)
{
    PetscErrorCode ierr;
    (void)Hpre; (void)ierr;

    PetscFunctionBegin;

    /* do nothing*/

    PetscFunctionReturn(0);
}



/****************************************************************************
 * Function: Monitor
 * Description: monitor the optimization process
 * input:
 *  tao    pointer to tao solver
 *  ptr    pointer to actual optimziation problem (has to be
 *         implemented by the user)
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimizationMonitor"
PetscErrorCode OptimizationMonitor(Tao tao, void* ptr)
{
    PetscErrorCode ierr;
    IntType iter;
    int iterdisp;
    char msg[256];
    ScalarType J=0.0,gnorm=0.0,cnorm=0.0,step=0.0,D=0.0,
                J0=0.0,D0=0.0,gnorm0=0.0,alpha,grtol,gatol,gttol;
    OptimizationProblem* optprob = NULL;
    TaoConvergedReason taoconvreason;
    TaoLineSearch ls=NULL;
    bool newtonkrylov;
    KSPConvergedReason kspconvreason;
    TaoLineSearchConvergedReason taolsconvreason;
    Vec x=NULL;
    KSP ksp=NULL;

    PetscFunctionBegin;


    optprob = static_cast<OptimizationProblem*>(ptr);
    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);

    // do we use a newton krylov method
    newtonkrylov =  (optprob->GetOptions()->GetOptMeth() == GAUSSNEWTON)
                 || (optprob->GetOptions()->GetOptMeth() == FULLNEWTON);

    if (newtonkrylov){

        if(optprob->GetOptions()->GetVerbosity() >= 2){

            std::string convmsg;
            ierr=TaoGetKSP(tao,&ksp); CHKERRQ(ierr);
            ierr=KSPGetConvergedReason(ksp,&kspconvreason);
            ierr=DispKSPConvReason(kspconvreason); CHKERRQ(ierr);

        }

    }

    ierr=TaoGetLineSearch(tao,&ls); CHKERRQ(ierr);
    ierr=TaoLineSearchGetSolution(ls,NULL,&J,NULL,&step,&taolsconvreason); CHKERRQ(ierr);

    if(optprob->GetOptions()->GetVerbosity() >= 2){
        ierr=DispLSConvReason(taolsconvreason); CHKERRQ(ierr);
    }


    // get current iteration, objective value, norm of gradient, norm of
    // contraint, step length / trust region radius and termination reason
    ierr=TaoGetSolutionStatus(tao,&iter,&J,&gnorm,&cnorm,&step,&taoconvreason); CHKERRQ(ierr);

    // remember current iterate
    optprob->SetNumOuterIter(iter);

    // tao: display convergence reason
    if(optprob->GetOptions()->GetVerbosity() >= 2){
        ierr=DispTaoConvReason(taoconvreason); CHKERRQ(ierr);
    }

    // compute l2 distance at current iteration
    ierr=optprob->EvaluateL2Distance(&D); CHKERRQ(ierr);

    // parse initial gradient and display header
    if ( optprob->InFirstIteration() == true ){

        ierr=PetscPrintf(MPI_COMM_WORLD,"%s\n",std::string(optprob->GetOptions()->GetLineLength(),'-').c_str());
        PetscPrintf(PETSC_COMM_WORLD," %s  %-20s %-20s %-20s %-20s %-20s\n",
                                    "iter","objective (rel)","mismatch (rel)",
                                    "||gradient||_2,rel","||gradient||_2","step");
        ierr=PetscPrintf(MPI_COMM_WORLD,"%s\n",std::string(optprob->GetOptions()->GetLineLength(),'-').c_str());

        optprob->SetInitialGradNorm(gnorm); // set the initial gradient norm
        optprob->SetInitialMismatchVal(D); // set the initial l2 distance
        optprob->SetInitialObjVal(J); // set the initial objective value

    }
    else{
        // handle multilevel/multiresolution problems
        if (iter == 0){

            // the tolerance of the solver has to be modified
            // so that the relative reduction is with respect
            // to the initial gradient

            gatol = optprob->GetOptions()->GetOptTol(0);  // ||g(x)||              <= gatol
            grtol = optprob->GetOptions()->GetOptTol(1);  // ||g(x)|| / |J(x)|     <= grtol
            gttol = optprob->GetOptions()->GetOptTol(2);  // ||g(x)|| / ||g(x0)||  <= gttol

            // get initial gradient
            gnorm0 = optprob->GetInitialGradNorm();
            gnorm0 = (gnorm0 > 0.0) ? gnorm0 : 1.0;
            alpha = PetscMax(1.0,gnorm0/gnorm);
            ierr=TaoSetTolerances(tao,gatol,grtol,alpha*gttol); CHKERRQ(ierr);
        }
    }

    // get initial gradient
    gnorm0 = optprob->GetInitialGradNorm();
    gnorm0 = (gnorm0 > 0.0) ? gnorm0 : 1.0;

    // get initial l2 distance
    D0 = optprob->GetInitialMismatchVal();
    D0 = (D0 > 0.0) ? D0 : 1.0;

    // get initial objective value
    J0 = optprob->GetInitialObjVal();
    J0 = (J0 > 0.0) ? J0 : 1.0;

    // get the solution vector and finalize the iteration
    ierr=TaoGetSolutionVector(tao,&x); CHKERRQ(ierr);
    ierr=optprob->FinalizeIteration(x); CHKERRQ(ierr);

    // display progress to user
    iterdisp = static_cast<int>(iter);
    sprintf(msg,"  %03d  %-20.12E %-20.12E %-20.12E %-20.12E %.6f",iterdisp,J/J0,D/D0,gnorm/gnorm0,gnorm,step);
    PetscPrintf(MPI_COMM_WORLD,"%-80s\n",msg);

    // make sure that everybody knows that we have performed the
    // first iteration already
    optprob->InFirstIteration(false);

    // go home
    PetscFunctionReturn(0);
}



/****************************************************************************
 * Function: KrylovMonitor
 * Description: monitor evolution of krylov subspace method
 *****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "KrylovMonitor"
PetscErrorCode KrylovMonitor(KSP ksp,IntType it,ScalarType rnorm,void* ptr)
{
    PetscErrorCode ierr;
    (void)ksp;
    OptimizationProblem* optprob=NULL;
    std::stringstream itss, rnss;
    std::string kspmeth, msg;

    PetscFunctionBegin;

    optprob = static_cast<OptimizationProblem*>(ptr);
    ierr=Assert(optprob!=NULL,"user is null pointer"); CHKERRQ(ierr);

    kspmeth="PCG  "; itss << std::setw(3) << it; rnss << std::scientific << rnorm;
    msg = kspmeth +  itss.str() + "  ||r||_2 = " + rnss.str();
    ierr=DbgMsg(msg); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/****************************************************************************
 * Function: DispKSPConvReason
 * Description:display the convergence reason of the KSP method
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
            msg="KSP convergence ||r||_2 < rtol ||b||_2";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ATOL_NORMAL:
        {
            msg="KSP convergence ||r||_2 < atol";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_RTOL:
        {
            msg="KSP convergence ||r||_2 < rtol ||b||_2";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ATOL:
        {
            msg="KSP convergence ||r||_2 < atol";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ITS:
        {
            //used by the KSPPREONLY solver after the single iteration of
            //the preconditioner is applied
            msg="KSP convergence k > maxit (KSPPREONLY)";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_CG_NEG_CURVE:
        {
            msg="KSP negative curvature detected";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_CG_CONSTRAINED:
        {
            msg="KSP convergence is reached along a constrained (full step goes beyond trust region)";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_STEP_LENGTH:
        {
            msg="KSP converged step length";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_HAPPY_BREAKDOWN:
        {
            msg="KSP converged happy break down";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NULL:
        {
            msg="KSP divergence detected";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_ITS:
        {
            msg="KSP max number of iterations reached";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_DTOL:
        {
            msg="KSP divergence detected (||r||_2 increased by a factor of divtol)";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_BREAKDOWN:
        {
            //breakdown in Krylov method was detected
            //method could not continue to enlarge Krylov subspace;
            //could be due to a singlular matrix or preconditioner
            msg="KSP generic breakdown; potential singular operator";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_BREAKDOWN_BICG:
        {
            msg="KSP initial ||r||_2 is orthogonal to preconditioned residual";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NONSYMMETRIC:
        {
            msg="KSP operators (A or P) are not symmetric";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_INDEFINITE_PC:
        {
            msg="KSP preconditioner is indefinite";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NANORINF:
        {
            msg="KSP: ||r||_2 is NAN or INF";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_INDEFINITE_MAT:
        {
            msg="KSP A is indefinite";
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
            msg="KSP convergence reason not defined";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }

    }


    PetscFunctionReturn(0);
}



/****************************************************************************
 * Function: DispLSConvReason
 * Description:display the convergence reason of the KSP method
 * Author: Andreas Mang
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DispLSConvReason"
PetscErrorCode DispLSConvReason(TaoLineSearchConvergedReason flag)
{
    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    switch(flag){
        case TAOLINESEARCH_FAILED_INFORNAN:
        {
            msg="LS: function evaluation gave INF or NAN";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_FAILED_BADPARAMETER:
        {
            msg="LS: bad parameter detected";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_FAILED_ASCENT:
        {
            msg="LS: search direction is not a descent direction";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_MAXFCN:
        {
            msg="LS: maximum number of function evaluations reached";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_UPPERBOUND:
        {
            msg="LS: step size reached upper bound";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_LOWERBOUND:
        {
            msg="LS: step size reached lower bound";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_RTOL:
        {
            msg="LS: range of uncertainty is smaller than given tolerance";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_OTHER:
        {
            msg="LS: line search stopped (other)";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_CONTINUE_ITERATING:
        {
//            msg="LS: continue iterating";
//            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_SUCCESS:
        {
            msg="LS: line serach was successfull";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        default:
        {
            msg="LS: reason not defined";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }

    }

    PetscFunctionReturn(0);
}



/****************************************************************************
 * Function: DispTaoConvReason
 * Description:display the convergence reason of the optimizer
 * Author: Andreas Mang
 ****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DispTaoConvReason"
PetscErrorCode DispTaoConvReason(TaoConvergedReason flag)
{
    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    switch(flag){
//        case TAO_CONVERGED_FATOL:
//        {
//            msg="TAO: convergence J(x) - J(x*) <= Jabstol";
//            ierr=DbgMsg(msg); CHKERRQ(ierr);
//            break;
//        }
//        case TAO_CONVERGED_FRTOL:
//        {
//            msg="TAO: convergence |J(x) - J(x*)|/|J(x)| <= Jreltol";
//            ierr=DbgMsg(msg); CHKERRQ(ierr);
//            break;
//        }
        case TAO_CONVERGED_GATOL:
        {
            msg="TAO: convergence ||g(x)|| <= gabstol";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_GRTOL:
        {
            msg="TAO: convergence ||g(x)||/J(x) <= greltol";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_GTTOL:
        {
            msg="TAO: convergence ||g(x)||/||g(x0)|| <= greltol";
            ierr=DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_STEPTOL:
        {
            msg="TAO: step size too small";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_MINF:
        {
            msg="TAO: objective value to small (J(x) < J_min)";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_USER:
        {
            msg="TAO: user defined convergence criteria met";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_MAXITS:
        {
            msg="TAO: maximum number of iterations reached";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_NAN:
        {
            msg="TAO: numerical problems (NAN detected)";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_MAXFCN:
        {
            msg="TAO: maximal number of function evaluations reached";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_LS_FAILURE:
        {
            msg="TAO: line search failed";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_TR_REDUCTION:
        {
            msg="TAO: trust region failed";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_USER:
        {
            msg="TAO: user defined divergence criterion met";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONTINUE_ITERATING:
        {
            // display nothing
            break;
        }
        default:
        {
            msg="TAO: convergence reason not defined";
            ierr=WrngMsg(msg); CHKERRQ(ierr);
            break;
        }

    }


    PetscFunctionReturn(0);
}

} // end of name space


#endif // _TAOINTERFACEREGISTRATION_H_
