/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CLAIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _KRYLOVINTERFACE_CPP_
#define _KRYLOVINTERFACE_CPP_

#include "KrylovInterface.hpp"


namespace reg {




/****************************************************************************
 * @brief monitor evolution of krylov subspace method
 ****************************************************************************/
PetscErrorCode KrylovMonitor(KSP krylovmethod, IntType it,
                             ScalarType rnorm, void* ptr) {
    PetscErrorCode ierr = 0;
    (void)krylovmethod;
    KSPConvergedReason reason;
    OptimizationProblem* optprob;
    std::stringstream itss, rnss;
    std::string kspmeth, msg;

    PetscFunctionBegin;

    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    if (optprob->GetOptions()->m_Verbosity > 0) {
        kspmeth = optprob->GetOptions()->m_KrylovMethod.name;
        itss << std::setw(5) << it;
        rnss << std::scientific << rnorm;
        msg = kspmeth + "   " + itss.str() + "  ||r||_2 = " + rnss.str();
        ierr = DbgMsg(msg); CHKERRQ(ierr);

        ierr = KSPGetConvergedReason(krylovmethod, &reason); CHKERRQ(ierr);
        ierr = DispKSPConvReason(reason); CHKERRQ(ierr);
    }

    if (optprob->GetOptions()->m_Log.enabled[LOGKSPRES]) {
        optprob->GetOptions()->LogKSPResidual(it, rnorm);
    }

    optprob->GetOptions()->m_KrylovMethod.iter = it;

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief preprocess right hand side and initial condition before entering
 * the krylov subspace method; in the context of numerical optimization this
 * means we preprocess the gradient and the incremental control variable
 * @para[in] krylovmethod pointer to krylov method
 * @para[in] b right hand side of equation
 * @para[in] x solution vector
 ****************************************************************************/
PetscErrorCode PreKrylovSolve(KSP krylovmethod, Vec b, Vec x, void* ptr) {
    PetscErrorCode ierr = 0;
    ScalarType gnorm = 0.0, g0norm = 0.0, reltol, abstol, divtol,
                uppergradbound, lowergradbound;
    IntType maxit;
    std::stringstream ss;
    std::string msg;
    OptimizationProblem* optprob = NULL;

    PetscFunctionBegin;

    uppergradbound = 0.5;
//#if defined(PETSC_USE_REAL_SINGLE)
//    lowergradbound = 1E-9;
//#else
    lowergradbound = 1E-16;
//#endif
    (void)krylovmethod;

    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    // set the iteration count to zero
    optprob->GetOptions()->m_KrylovMethod.iter = 0;

    if ((optprob->GetOptions()->m_KrylovMethod.matvectype == PRECONDMATVEC)
     || (optprob->GetOptions()->m_KrylovMethod.matvectype == PRECONDMATVECSYM))  {
        // get current gradient and compute norm
        // before we apply the preconditioner to the right hand side
        ierr = VecNorm(b, NORM_2, &gnorm); CHKERRQ(ierr);
    }

    // use default tolreance
    reltol = optprob->GetOptions()->m_KrylovMethod.tol[0];

    // apply preprocessing to right hand side and initial condition
    ierr = optprob->PreKrylovSolve(b, x); CHKERRQ(ierr);

    // user forcing sequence to estimate adequate tolerance
    // for solution of KKT system (Eisenstat-Walker)
    if (optprob->GetOptions()->m_KrylovMethod.fseqtype != NOFS) {
        if (gnorm == 0.0) {
            // get current gradient and compute norm
            ierr = VecNorm(b, NORM_2, &gnorm); CHKERRQ(ierr);
        }

        if (!optprob->GetOptions()->m_KrylovMethod.g0normset) {
            if (optprob->GetOptions()->m_Verbosity > 1) {
                ss << std::fixed << std::scientific << gnorm;
                msg = optprob->GetOptions()->m_KrylovMethod.name +
                    ": setting initial ||g|| in krylov method (" + ss.str() + ")";
                ierr = DbgMsg(msg); CHKERRQ(ierr);
                ss.str(std::string()); ss.clear();
            }
            optprob->GetOptions()->m_KrylovMethod.g0norm = gnorm;
            optprob->GetOptions()->m_KrylovMethod.g0normset = true;
        }

        // get initial value for gradient
        g0norm = optprob->GetOptions()->m_KrylovMethod.g0norm;
        ierr = Assert(g0norm > 0.0, "initial gradient is zero"); CHKERRQ(ierr);

        // normalize
        gnorm /= g0norm;

        // get current tolerances
        ierr = KSPGetTolerances(krylovmethod, &reltol, &abstol, &divtol, &maxit); CHKERRQ(ierr);

        if (optprob->GetOptions()->m_KrylovMethod.fseqtype == QDFS) {
            // assuming quadratic convergence (we do not solver more
            // accurately than 12 digits)
            reltol = PetscMax(lowergradbound, PetscMin(uppergradbound, gnorm));
        } else {
            // assuming superlinear convergence (we do not solver
            // more accurately than 12 digitis)
            reltol = PetscMax(lowergradbound, PetscMin(uppergradbound, std::sqrt(gnorm)));
        }

        // overwrite tolerances with estimate
        ierr = KSPSetTolerances(krylovmethod, reltol, abstol, divtol, maxit); CHKERRQ(ierr);
    }

    // pass tolerance to optimization problem (for preconditioner)
    optprob->GetOptions()->m_KrylovMethod.reltol = reltol;

    if (optprob->GetOptions()->m_Verbosity > 0) {
        ss << std::fixed << std::scientific << reltol;
        msg = optprob->GetOptions()->m_KrylovMethod.name +
              ": computing solution of hessian system (tol=" + ss.str() + ")";
        ierr = DbgMsg(msg); CHKERRQ(ierr);
    }

    // we might want to recompute eigenvalues at every Newton iteration;
    // if the user sets this flag, we pretend the eigenvalues have not
    // been estimated just yet
    if (optprob->GetOptions()->m_KrylovMethod.reesteigvals == 2) {
        optprob->GetOptions()->m_KrylovMethod.eigvalsestimated = false;
    }

    // check symmetry of hessian
    if (optprob->GetOptions()->m_KrylovMethod.checkhesssymmetry) {
        ierr = optprob->HessianSymmetryCheck(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief postprocess right hand side and initial condition before entering
 * the krylov subspace method; in the context of numerical optimization this
 * means we postprocess the gradient and the incremental control variable
 * @para[in] krylovmethod pointer to krylov method
 * @para[in] b right hand side of equation
 * @para[in] x solution vector
 ****************************************************************************/
PetscErrorCode PostKrylovSolve(KSP krylovmethod, Vec b, Vec x, void* ptr) {
    PetscErrorCode ierr = 0;
    OptimizationProblem* optprob = NULL;
    KSPConvergedReason reason;
    std::string convmsg;

    PetscFunctionBegin;

    (void)krylovmethod;

    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr = optprob->PostKrylovSolve(b, x); CHKERRQ(ierr);

    if (optprob->GetOptions()->m_Verbosity > 0) {
        ierr = KSPGetConvergedReason(krylovmethod, &reason);
        ierr = DispKSPConvReason(reason); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @briefdisplay the convergence reason of the KSP method
 ****************************************************************************/
PetscErrorCode DispKSPConvReason(KSPConvergedReason flag) {
    PetscErrorCode ierr = 0;
    std::string msg;
    PetscFunctionBegin;

    switch (flag) {
        case KSP_CONVERGED_ITERATING:
        {
            // don't display anaything and continue iterating
            break;
        }
        case KSP_CONVERGED_RTOL_NORMAL:
        {
            msg = "krylov method converged: ||r||_2 < tol ||b||_2";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ATOL_NORMAL:
        {
            msg = "krylov method converged: ||r||_2 < tol";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_RTOL:
        {
            msg = "krylov method converged: ||r||_2 < tol ||b||_2";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ATOL:
        {
            msg = "krylov method converged: ||r||_2 < tol";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_ITS:
        {
            //used by the KSPPREONLY solver after the single iteration of
            //the preconditioner is applied
            msg = "krylov method converged: k > maxit";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_CG_NEG_CURVE:
        {
            msg = "krylov method: negative curvature detected";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_CG_CONSTRAINED:
        {
            msg = "krylov method: convergence is reached";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_STEP_LENGTH:
        {
            msg = "krylov method: converged step length";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_CONVERGED_HAPPY_BREAKDOWN:
        {
            msg = "krylov method: converged happy break down";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NULL:
        {
            msg = "krylov method: divergence detected";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_ITS:
        {
            msg = "krylov method: max number of iterations reached";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_DTOL:
        {
            msg = "krylov method: divergence detected (||r||_2 increased by a factor of divtol)";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_BREAKDOWN:
        {
            //breakdown in Krylov method was detected
            //method could not continue to enlarge Krylov subspace;
            //could be due to a singlular matrix or preconditioner
            msg = "krylov method: generic breakdown (potential singular operator)";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_BREAKDOWN_BICG:
        {
            msg = "krylov method: initial ||r||_2 is orthogonal to preconditioned r";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NONSYMMETRIC:
        {
            msg = "krylov method: operators (A or P) are not symmetric";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_INDEFINITE_PC:
        {
            msg = "krylov method: preconditioner is indefinite";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_NANORINF:
        {
            msg = "krylov method: ||r||_2 is NAN or INF";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case KSP_DIVERGED_INDEFINITE_MAT:
        {
            msg = "krylov method: A is indefinite";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        default:
        {
            msg = "krylov method: convergence reason not defined";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }

    }
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief monitor for evolution of krylov subspace method for
 * inversion of preconditioner
 * @para[in] krylovmethod pointer to krylov method
 * @para[in] current iteration index
 * @para[in] current norm of residual
 * @para[in] pointer to preconditioner object
 *******************************************************************/
PetscErrorCode InvertPrecondKrylovMonitor(KSP krylovmethod, IntType it,
                                          ScalarType rnorm, void* ptr) {
    PetscErrorCode ierr = 0;
    (void)krylovmethod;
    KSPConvergedReason reason;
    Preconditioner* precond = NULL;
    std::stringstream itss, rnss;
    std::string kspmeth, msg;

    PetscFunctionBegin;

    precond = reinterpret_cast<Preconditioner*>(ptr);
    ierr = Assert(precond != NULL, "null pointer"); CHKERRQ(ierr);

    kspmeth = " >> PRECOND " + precond->GetOptions()->m_KrylovMethod.pcname;
    itss << std::setw(5) << it;
    rnss << std::scientific << rnorm;
    msg = kspmeth +  itss.str() + "  ||r||_2 = " + rnss.str();
    ierr = DbgMsg(msg); CHKERRQ(ierr);

    ierr = KSPGetConvergedReason(krylovmethod, &reason); CHKERRQ(ierr);
    ierr = DispKSPConvReason(reason); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief computes the hessian matrix vector product Hx = H*xtilde
 ****************************************************************************/
PetscErrorCode InvertPrecondMatVec(Mat P, Vec x, Vec Px) {
    PetscErrorCode ierr = 0;
    void* ptr;
    Preconditioner* precond = NULL;
    PetscFunctionBegin;

    ierr = MatShellGetContext(P, reinterpret_cast<void**>(&ptr)); CHKERRQ(ierr);
    precond = reinterpret_cast<Preconditioner*>(ptr);
    ierr = Assert(precond != NULL, "null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr = precond->HessianMatVec(Px, x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief initialization before applying the preconditioner
 * @para[in] krylovmethod pointer to krylov method
 * @para[in] b right hand side of equation
 * @para[in] x solution vector
 ****************************************************************************/
PetscErrorCode InvertPrecondPreKrylovSolve(KSP krylovmethod, Vec b,
                                           Vec x, void* ptr) {
    PetscErrorCode ierr = 0;
    Preconditioner* precond = NULL;
    IntType maxits;
    ScalarType reltol, abstol, divtol, scale, lowerbound, upperbound;
    std::stringstream itss, rnss;
    std::string msg;

    PetscFunctionBegin;

    precond = reinterpret_cast<Preconditioner*>(ptr);
    ierr = Assert(precond != NULL, "null pointer"); CHKERRQ(ierr);

    // setup preconditioner
    if (!precond->GetOptions()->m_KrylovMethod.pcsetupdone) {
        ierr = precond->DoSetup(); CHKERRQ(ierr);
    }

    // set the default values
    upperbound = 5E-1;
//#if defined(PETSC_USE_REAL_SINGLE)
//#else
    lowerbound = 1E-16;
//#endif
    maxits = 1E3;
    reltol = precond->GetOptions()->m_KrylovMethod.pctol[0];
    abstol = precond->GetOptions()->m_KrylovMethod.pctol[1];
    divtol = precond->GetOptions()->m_KrylovMethod.pctol[2];
    scale  = precond->GetOptions()->m_KrylovMethod.pctolscale;

    switch (precond->GetOptions()->m_KrylovMethod.pcsolver) {
        case CHEB:
        {
            // chebyshev iteration
            maxits = precond->GetOptions()->m_KrylovMethod.pcmaxit;
            // estimate the eigenvalues
            ierr = precond->EstimateEigenValues(); CHKERRQ(ierr);
            break;
        }
        case PCG:
        {
            // flexible conjugate gradient or gmres as outer
            if ( precond->GetOptions()->m_KrylovMethod.solver == FCG
              || precond->GetOptions()->m_KrylovMethod.solver == FGMRES) {
                // preconditioned conjugate gradient
                maxits = precond->GetOptions()->m_KrylovMethod.pcmaxit;
            } else {
                reltol = scale*precond->GetOptions()->m_KrylovMethod.reltol;
            }
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient or gmres as outer
            if ( precond->GetOptions()->m_KrylovMethod.solver == FCG
              || precond->GetOptions()->m_KrylovMethod.solver == FGMRES) {
                // preconditioned conjugate gradient
                maxits = precond->GetOptions()->m_KrylovMethod.pcmaxit;
            } else {
                reltol = scale*precond->GetOptions()->m_KrylovMethod.reltol;
            }
            break;
        }
        case GMRES:
        {
            // flexible conjugate gradient or gmres as outer
            if ( precond->GetOptions()->m_KrylovMethod.solver == FCG
              || precond->GetOptions()->m_KrylovMethod.solver == FGMRES) {
                // preconditioned conjugate gradient
                maxits = precond->GetOptions()->m_KrylovMethod.pcmaxit;
            } else {
                reltol = scale*precond->GetOptions()->m_KrylovMethod.reltol;
            }
            break;
        }
        case FGMRES:
        {
            // flexible conjugate gradient or gmres as outer
            if ( precond->GetOptions()->m_KrylovMethod.solver == FCG
              || precond->GetOptions()->m_KrylovMethod.solver == FGMRES) {
                // preconditioned conjugate gradient
                maxits = precond->GetOptions()->m_KrylovMethod.pcmaxit;
            } else {
                reltol = scale*precond->GetOptions()->m_KrylovMethod.reltol;
            }
            break;
        }
        default:
        {
            ierr = ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }

    reltol = std::max(reltol, lowerbound);  // make sure tolerance is non-zero
    reltol = std::min(reltol, upperbound);   // make sure tolerance smaller than 0.25

    // set tolerances
    ierr = KSPSetTolerances(krylovmethod, reltol, abstol, divtol, maxits); CHKERRQ(ierr);

    if (precond->GetOptions()->m_Verbosity > 1) {
        rnss << std::fixed << std::scientific << reltol;
        itss << maxits;
        msg = precond->GetOptions()->m_KrylovMethod.pcname
            + ": inverting preconditioner (tol=" + rnss.str()
            + ", maxit=" + itss.str() + ")";
        ierr = DbgMsg(msg); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _KRYLOVINTERFACE_CPP_
