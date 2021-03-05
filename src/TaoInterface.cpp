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

#ifndef _TAOINTERFACE_CPP_
#define _TAOINTERFACE_CPP_

#include "TaoInterface.hpp"




namespace reg {




/****************************************************************************
 * @brief general purpose function to evaluate objective
 * @param tao pointer to tao solver
 * @param x iterate x objective is to be evaluated at
 * @param Jx objective value J(x)
 * @param ptr pointer to optimziation problem (has to be implemented by user)
 ****************************************************************************/
PetscErrorCode EvaluateObjective(Tao tao, Vec x, ScalarType* Jx, void* ptr) {
    PetscErrorCode ierr;
    OptimizationProblem *optprob = NULL;

    PetscFunctionBegin;

    (void)tao;
    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);
    
    // compute objective value
    ierr = optprob->EvaluateObjective(Jx, x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief general purpose function to evaluate the gradient
 * @param tao pointer to tao solver
 * @param x iterate x gradient is to be evaluated at
 * @param gx gradient of objective functional, i.e. g(x)
 * @param ptr pointer to actual optimziation problem
 ****************************************************************************/
PetscErrorCode EvaluateGradient(Tao tao, Vec x, Vec gx, void* ptr) {
    PetscErrorCode ierr = 0;
    OptimizationProblem* optprob = NULL;

    PetscFunctionBegin;

    (void)tao;
    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    // compute gradient
    ierr = optprob->EvaluateGradient(gx, x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief general purpose function to evaluate objective and gradient
 * @param tao pointer to tao solver
 * @param x iterate x gradient is to be evaluated at
 * @param Jx objective value J(x)
 * @param gx gradient of objective functional, i.e. g(x)
 * @param  ptr pointer to actual optimziation problem
 ****************************************************************************/
PetscErrorCode EvaluateObjectiveGradient(Tao tao, Vec x, ScalarType* Jx, Vec gx, void* ptr) {
    PetscErrorCode ierr = 0;
    OptimizationProblem* optprob = NULL;
    (void)tao;

    PetscFunctionBegin;

    // cast pointer
    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    // evaluate objective and gradient
    ierr = optprob->EvaluateObjective(Jx, x); CHKERRQ(ierr);
    ierr = optprob->EvaluateGradient(gx, x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief this function is called before the hessian matvec; I used to have
 * the initialization of the krylov method here; it's now moved to the
 * preprocessing routine for the krylov method;
 * @param tao pointer to tao solver
 * @param x iterate x objective is to be evaluated at
 * @param H hessian matrix (can be matrix free)
 * @param Hpre preconditioner (can be matrix free)
 * @param ptr pointer to optimzation problem (has to be implemented by user)
 ****************************************************************************/
PetscErrorCode EvaluateHessian(Tao tao, Vec x, Mat H, Mat Hpre, void* ptr) {
    PetscErrorCode ierr = 0;
    (void)tao; (void)x; (void)H; (void)Hpre; (void)ptr;

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief computes the hessian matrix vector product Hx = H*xtilde
 ****************************************************************************/
PetscErrorCode HessianMatVec(Mat H, Vec x, Vec Hx) {
    PetscErrorCode ierr = 0;
    void* ptr;
    OptimizationProblem *optprob = NULL;

    PetscFunctionBegin;

    ierr = MatShellGetContext(H, reinterpret_cast<void**>(&ptr)); CHKERRQ(ierr);

    // cast pointer
    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr = optprob->HessianMatVec(Hx, x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief computes the matrix vector product Px
 ****************************************************************************/
PetscErrorCode PrecondMatVec(PC Hpre, Vec x, Vec Hprex) {
    PetscErrorCode ierr = 0;
    void* ptr;
    Preconditioner *preconditioner = NULL;

    PetscFunctionBegin;

    ierr = PCShellGetContext(Hpre, &ptr); CHKERRQ(ierr);

    preconditioner = reinterpret_cast<Preconditioner*>(ptr);
    ierr = Assert(preconditioner != NULL, "null pointer"); CHKERRQ(ierr);

    // apply hessian
    ierr = preconditioner->MatVec(Hprex, x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief convergence test for optimization
 * @param tao pointer to tao solver
 * @param optprob pointer to optimziation problem
 ****************************************************************************/
PetscErrorCode CheckConvergenceGradObj(Tao tao, void* ptr) {
    PetscErrorCode ierr = 0;
    IntType iter, maxiter, miniter, iterbound;
    OptimizationProblem* optprob = NULL;
    std::stringstream ss, sc;
    ScalarType jx, jxold, gnorm, step, gatol, grtol,
                gttol, g0norm, gtolbound, minstep, theta,
                normx, normdx, tolj, tolx, tolg;
    const int nstop = 7;
    bool stop[nstop];
    Vec x;

    PetscFunctionBegin;

    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    miniter   = optprob->GetOptions()->m_OptPara.miniter;     // lower bound for iterations
    iterbound = optprob->GetOptions()->m_OptPara.iterbound;   // upper bound for iterations
    gtolbound = optprob->GetOptions()->m_OptPara.gtolbound;   // lower bound for gradient

    jxold     = optprob->GetOptions()->m_Monitor.jvalold;     // value of objective from former iteration
    g0norm    = optprob->GetOptions()->m_Monitor.gradnorm0;   // initial norm of gradient

    // compute minimal step size
    minstep = std::pow(2.0, 10.0);
    minstep = 1.0 / minstep;

    // get initial gradient
    g0norm = (g0norm > 0.0) ? g0norm : 1.0;

    // get user defined tolerances
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
    ierr = TaoGetTolerances(tao, &gatol, &grtol, &gttol); CHKERRQ(ierr);
#else
    ierr = TaoGetTolerances(tao, NULL, NULL, &gatol, &grtol, &gttol); CHKERRQ(ierr);
#endif
    ierr = TaoGetMaximumIterations(tao, &maxiter); CHKERRQ(ierr);
    if (maxiter > iterbound) iterbound = maxiter;

    // compute tolerances for stopping conditions
    tolj = gttol;
    tolx = std::sqrt(gttol);
#if __cplusplus > 199711L
    tolg = std::cbrt(gttol);
#else
    tolg = std::pow(gttol, 1.0/3.0);
#endif

    // get solution status
    ierr = TaoGetSolutionStatus(tao, &iter, &jx, &gnorm, NULL, &step, NULL); CHKERRQ(ierr);
    ierr = TaoGetSolutionVector(tao, &x); CHKERRQ(ierr);

    // compute theta
    theta = 1.0 + std::abs(jx);

    ierr = optprob->ComputeUpdateNorm(x, normdx, normx); CHKERRQ(ierr);

    // check for NaN value
    if (PetscIsInfOrNanReal(jx)) {
        ierr = WrngMsg("objective is NaN"); CHKERRQ(ierr);
        ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_NAN); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    // check for NaN value
    if (PetscIsInfOrNanReal(gnorm)) {
        ierr = WrngMsg("||g|| is NaN"); CHKERRQ(ierr);
        ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_NAN); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    optprob->Converged(false);
    // initialize flags for stopping conditions
    for (int i = 0; i < nstop; ++i) stop[i] = false;

    // only check convergence criteria after a certain number of iterations
    if (iter >= miniter && iter > 1) {
        if (step < minstep) {
            ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_STEPTOL); CHKERRQ(ierr);
            PetscFunctionReturn(ierr);
        }

        // |j(x_{k-1}) - j(x_k)| < tolj*abs(1+J)
        if (std::abs(jxold-jx) < tolj*theta) {
            stop[0] = true;
        }
        ss << "[  " << stop[0] << "    |dJ|  = " << std::setw(14)
           << std::right << std::scientific << std::abs(jxold-jx) << "    <    "
           << std::left << std::setw(14) << tolj*theta << " = " << "tol*|1+J|";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        ss.str(std::string()); ss.clear();

        // ||dx|| < sqrt(tolj)*(1+||x||)
        if (normdx < tolx*(1+normx)) {
            stop[1] = true;
        }
        ss << "[  " << stop[1] << "    |dx|  = " << std::setw(14)
           << std::right << std::scientific << normdx << "    <    "
           << std::left << std::setw(14) << tolx*(1+normx) << " = " << "sqrt(tol)*(1+||x||)";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        ss.str(std::string()); ss.clear();

        // ||g_k||_2 < cbrt(tolj)*abs(1+Jc)
        if (gnorm < tolg*theta) {
            stop[2] = true;
        }
        ss  << "[  " << stop[2] << "    ||g|| = " << std::setw(14)
            << std::right << std::scientific << gnorm << "    <    "
            << std::left << std::setw(14) << tolg*theta << " = " << "cbrt(tol)*|1+J|";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        ss.str(std::string()); ss.clear();

        // ||g_k||_2 < tol
        if (gnorm < gatol) {
            stop[3] = true;
        }
        ss  << "[  " << stop[3] << "    ||g|| = " << std::setw(14)
            << std::right << std::scientific << gnorm << "    <    "
            << std::left << std::setw(14) << gatol << " = " << "tol";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        ss.str(std::string()); ss.clear();

        if (gnorm < gtolbound*g0norm) {
            stop[4] = true;
        }
        ss  << "[  " << stop[4] << "    ||g|| = " << std::setw(14)
            << std::right << gnorm  << "    >    "
            << std::left << std::setw(14) << gtolbound*g0norm << " = " << "kappa*||g0||";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        ss.str(std::string()); ss.clear();

        if (iter > maxiter) {
            stop[5] = true;
        }
        ss  << "[  " << stop[5] << "    iter  = " << std::setw(14)
            << std::right << iter  << "    >    "
            << std::left << std::setw(14) << maxiter << " = " << "maxiter";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        ss.str(std::string()); ss.clear();

        if (iter > iterbound) {
            stop[6] = true;
        }
        ss  << "[  " << stop[6] << "    iter  = " << std::setw(14)
            << std::right << iter  << "    >    "
            << std::left << std::setw(14) << iterbound << " = " << "iterbound";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        sc << std::endl;
        ss.str(std::string()); ss.clear();

        optprob->SetConvergenceMessage(sc.str());
        optprob->GetOptions()->m_Monitor.jvalold = jx;

        if (stop[0] && stop[1] && stop[2]) {
            ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_USER); CHKERRQ(ierr);
            optprob->Converged(true);
            PetscFunctionReturn(ierr);
        } else if (stop[3]) {
            ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL); CHKERRQ(ierr);
            optprob->Converged(true);
            PetscFunctionReturn(ierr);
        } else if (stop[4] && stop[5]) {
            ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_MAXITS); CHKERRQ(ierr);
            optprob->Converged(true);
            PetscFunctionReturn(ierr);
        } else if (stop[6]) {
            ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_MAXITS); CHKERRQ(ierr);
            optprob->Converged(true);
            PetscFunctionReturn(ierr);
        }
    } else {
        // if the gradient is zero, we should terminate immediately
        if (gnorm == 0) {
            ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL); CHKERRQ(ierr);
            PetscFunctionReturn(ierr);
        }
    }

    // perform derivative check
    if (optprob->GetOptions()->m_OptPara.derivativecheckenabled) {
//        ierr = optprob->DerivativeCheckGradient(); CHKERRQ(ierr);
        ierr = optprob->DerivativeCheckHessian(); CHKERRQ(ierr);
//        ierr = optprob->DerivativeCheckHessianFD(); CHKERRQ(ierr);
    }

    // if we're here, we're good to go
    ierr = TaoSetConvergedReason(tao, TAO_CONTINUE_ITERATING); CHKERRQ(ierr);

    // go home
    PetscFunctionReturn(ierr);
}






/****************************************************************************
 * @brief convergence test for optimization
 * @param tao pointer to tao solver
 * @param optprob pointer to optimziation problem
 ****************************************************************************/
PetscErrorCode CheckConvergenceGrad(Tao tao, void* ptr) {
    PetscErrorCode ierr = 0;
    IntType iter, maxiter, miniter;
    OptimizationProblem* optprob = NULL;
    ScalarType J, gnorm, step, gatol, grtol, gttol, g0norm, minstep;
    bool stop[3];
    int verbosity;
    std::stringstream ss, sc;

    PetscFunctionBegin;

    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    verbosity = optprob->GetOptions()->m_Verbosity;

    minstep = std::pow(2.0, 10.0);
    minstep = 1.0 / minstep;
    miniter = optprob->GetOptions()->m_OptPara.miniter;

    // get initial gradient
    g0norm = optprob->GetOptions()->m_Monitor.gradnorm0;
    g0norm = (g0norm > 0.0) ? g0norm : 1.0;

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
    ierr = TaoGetTolerances(tao, &gatol, &grtol, &gttol); CHKERRQ(ierr);
#else
    ierr = TaoGetTolerances(tao, NULL, NULL, &gatol, &grtol, &gttol); CHKERRQ(ierr);
#endif

    ierr = TaoGetMaximumIterations(tao, &maxiter); CHKERRQ(ierr);
    ierr = TaoGetSolutionStatus(tao, &iter, &J, &gnorm, NULL, &step, NULL); CHKERRQ(ierr);

    // check for NaN value
    if (PetscIsInfOrNanReal(J)) {
        ierr = WrngMsg("objective is NaN"); CHKERRQ(ierr);
        ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_NAN); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    // check for NaN value
    if (PetscIsInfOrNanReal(gnorm)) {
        ierr = WrngMsg("||g|| is NaN"); CHKERRQ(ierr);
        ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_NAN); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    // only check convergence criteria after a certain number of iterations
    stop[0] = false; stop[1] = false; stop[2] = false;
    optprob->Converged(false);
    if (iter >= miniter) {
        if (verbosity > 1) {
            ss << "step size in linesearch: " << std::scientific << step;
            ierr = DbgMsg1(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        if (step < minstep) {
            ss << "step  = " << std::scientific << step << " < " << minstep << " = " << "bound";
            ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
            ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_STEPTOL); CHKERRQ(ierr);
            PetscFunctionReturn(ierr);
        }

        // ||g_k||_2 < tol*||g_0||
        if (gnorm < gttol*g0norm) {
            ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_GTTOL); CHKERRQ(ierr);
            stop[0] = true;
        }
        ss << "[  " << stop[0] << "    ||g|| = " << std::setw(14)
           << std::right << std::scientific << gnorm << "    <    "
           << std::left << std::setw(14) << gttol*g0norm << " = " << "tol";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        ss.str(std::string()); ss.clear();

        // ||g_k||_2 < tol
        if (gnorm < gatol) {
            ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL); CHKERRQ(ierr);
            stop[1] = true;
        }
        ss  << "[  " << stop[1] << "    ||g|| = " << std::setw(14)
            << std::right << std::scientific << gnorm << "    <    "
            << std::left << std::setw(14) << gatol << " = " << "tol";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        ss.str(std::string()); ss.clear();

        // iteration number exceeds limit
        if (iter > maxiter) {
            ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_MAXITS); CHKERRQ(ierr);
            stop[2] = true;
        }
        ss  << "[  " << stop[2] << "     iter = " << std::setw(14)
            << std::right << iter  << "    >    "
            << std::left << std::setw(14) << maxiter << " = " << "maxiter";
        sc << std::left << std::setw(100) << ss.str() << "]" << std::endl;
        sc << std::endl;
        ss.str(std::string()); ss.clear();

        optprob->SetConvergenceMessage(sc.str());

        if (stop[0] || stop[1] || stop[2]) {
            optprob->Converged(true);
            PetscFunctionReturn(ierr);
        }

    } else {
        // if the gradient is zero, we should terminate immediately
        if (gnorm == 0) {
            ss << "||g|| = " << std::scientific << 0.0 << " < " << gatol  << " = " << "bound";
            ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
            ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL); CHKERRQ(ierr);
            PetscFunctionReturn(ierr);
        }
    }

    // perform derivative check
    if (optprob->GetOptions()->m_OptPara.derivativecheckenabled) {
//        ierr = optprob->DerivativeCheckGradient(); CHKERRQ(ierr);
        ierr = optprob->DerivativeCheckHessian(); CHKERRQ(ierr);
//        ierr = optprob->DerivativeCheckHessianFD(); CHKERRQ(ierr);
    }

    // if we're here, we're good to go
    ierr = TaoSetConvergedReason(tao, TAO_CONTINUE_ITERATING); CHKERRQ(ierr);

    // go home
    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief monitor the optimization process
 * @param tao pointer to tao solver
 * @param ptr pointer to optimziation problem (has to be implemented by user)
 ****************************************************************************/
PetscErrorCode OptimizationMonitor(Tao tao, void* ptr) {
    PetscErrorCode ierr = 0;
    IntType iter;
    int iterdisp, cnthess, cntgrad;
    char msg[256];
    std::string statusmsg;
    ScalarType J, gnorm, step, D, J0, D0, gnorm0;
    OptimizationProblem* optprob = NULL;
    Vec x = NULL;
    TaoConvergedReason convreason;

    PetscFunctionBegin;

    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    // reinit the initialization of the preconditioner
    optprob->GetOptions()->m_KrylovMethod.pcsetupdone = false;

    if (optprob->GetOptions()->m_Verbosity > 1) {
        ierr = GetLineSearchStatus(tao, optprob); CHKERRQ(ierr);
    }

    // get current iteration, objective value, norm of gradient, norm of
    // contraint, step length / trust region radius and termination reason
    ierr = TaoGetSolutionStatus(tao, &iter, &J, &gnorm, NULL, &step, &convreason); CHKERRQ(ierr);

    // save gradient norm
    optprob->GetOptions()->m_Monitor.gradnorm = gnorm;
    
    if (iter > 0) {
    // remember current iterate
      optprob->IncrementIterations();
    }

    // tao: display convergence reason
    if (optprob->GetOptions()->m_Verbosity > 0) {
        ierr = GetSolverStatus(convreason, statusmsg); CHKERRQ(ierr);
        optprob->GetOptions()->m_Monitor.solverstatus = statusmsg;
    }

    // compute l2 distance at current iteration
    D = optprob->GetOptions()->m_Monitor.dval;

    // get initial gradient
    gnorm0 = optprob->GetOptions()->m_Monitor.gradnorm0;
    gnorm0 = (gnorm0 > 0.0) ? gnorm0 : 1.0;

    // get initial l2 distance
    D0 = optprob->GetOptions()->m_Monitor.dval0;
    D0 = (D0 > 0.0) ? D0 : 1.0;

    // get initial objective value
    J0 = optprob->GetOptions()->m_Monitor.jval0;
    J0 = (J0 > 0.0) ? J0 : 1.0;

    // get the solution vector and finalize the iteration
    ierr = TaoGetSolutionVector(tao, &x); CHKERRQ(ierr);
    ierr = optprob->FinalizeIteration(x); CHKERRQ(ierr);
    
    cnthess = optprob->GetOptions()->GetCounter(HESSMATVEC);
    cntgrad = optprob->GetOptions()->GetCounter(OBJEVAL);

    // display progress to user
    iterdisp = static_cast<int>(iter);
    sprintf(msg, "  %03d %04d %04d  %-18.12E %-18.12E %-18.12E %-18.12E %.6f",
            iterdisp, cnthess, cntgrad, J/J0, D/D0, gnorm/gnorm0, gnorm, step);
    PetscPrintf(MPI_COMM_WORLD, "%-80s\n", msg);

    if (optprob->GetOptions()->m_Monitor.qmval != -1.0) {
        sprintf(msg, "  >> l2 inner product for q and m1: %-20.12E",optprob->GetOptions()->m_Monitor.qmval);
        PetscPrintf(MPI_COMM_WORLD, "%-80s\n", msg);
    }

    // go home
    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief display the convergence reason of the KSP method
 ****************************************************************************/
PetscErrorCode GetLineSearchStatus(Tao tao, void* ptr) {
    PetscErrorCode ierr = 0;
    std::string msg;
    ScalarType J, step;
    TaoLineSearchConvergedReason flag;
    OptimizationProblem* optprob = NULL;
    TaoLineSearch ls = NULL;

    PetscFunctionBegin;

    optprob = reinterpret_cast<OptimizationProblem*>(ptr);
    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = TaoGetLineSearch(tao, &ls); CHKERRQ(ierr);

#if PETSC_VERSION_GE(3,13,4)
    ierr = TaoLineSearchGetSolution(ls, 0, &J, 0, &step, &flag); CHKERRQ(ierr);
#else
		IntType nl, ng;
    Vec x = nullptr, g = nullptr;
    if (optprob->GetOptions()->m_Verbosity < 2 && !(optprob->GetOptions()->m_Log.enabled[LOGKSPRES]) ) {
      ierr = TaoLineSearchGetSolution(ls, 0, &J, 0, &step, &flag); CHKERRQ(ierr);
    } else {
      nl = optprob->GetOptions()->m_Domain.nl;
      ng = optprob->GetOptions()->m_Domain.ng;
      ierr = VecCreate(x, nl, ng); CHKERRQ(ierr);
      ierr = VecCreate(g, 3*nl, 3*ng); CHKERRQ(ierr);
      ierr = TaoLineSearchGetSolution(ls, x, &J, g, &step, &flag); CHKERRQ(ierr);
      if (x != NULL) {ierr = VecDestroy(&x); CHKERRQ(ierr); x = NULL;}
      if (g != NULL) {ierr = VecDestroy(&g); CHKERRQ(ierr); g = NULL;}
    }
#endif

    switch(flag) {
        case TAOLINESEARCH_FAILED_INFORNAN:
        {
            msg = "linesearch: function evaluation gave INF or NaN";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_FAILED_BADPARAMETER:
        {
            msg = "linesearch: bad parameter detected";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_FAILED_ASCENT:
        {
            msg = "linesearch: search direction is not a descent direction";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_MAXFCN:
        {
            msg = "linesearch: maximum number of function evaluations reached";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_UPPERBOUND:
        {
            msg = "linesearch: step size reached upper bound";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_LOWERBOUND:
        {
            msg = "linesearch: step size reached lower bound";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_RTOL:
        {
            msg = "linesearch: range of uncertainty is smaller than given tolerance";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_HALTED_OTHER:
        {
            msg = "linesearch: stopped (other)";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAOLINESEARCH_CONTINUE_ITERATING:
        {
            // do nothing, cause everything's fine
            break;
        }
        case TAOLINESEARCH_SUCCESS:
        {
            msg = "linesearch: successful";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        default:
        {
            msg = "linesearch: status not defined";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
    }

    PetscFunctionReturn(ierr);
}




/****************************************************************************
 * @brief display the convergence reason of the optimizer
 ****************************************************************************/
PetscErrorCode GetSolverStatus(TaoConvergedReason flag, std::string& msg) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    switch (flag) {
        case TAO_CONVERGED_GATOL:
        {
            msg = "solver converged: ||g(x)|| <= tol";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_GRTOL:
        {
            msg = "solver converged: ||g(x)||/J(x) <= tol";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_GTTOL:
        {
            msg = "solver converged: ||g(x)||/||g(x0)|| <= tol";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_STEPTOL:
        {
            msg = "step size too small";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_MINF:
        {
            msg = "objective value to small";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONVERGED_USER:
        {
            msg = "solver converged";
            ierr = DbgMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_MAXITS:
        {
            msg = "maximum number of iterations reached";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_NAN:
        {
            msg = "numerical problems (NAN detected)";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_MAXFCN:
        {
            msg = "maximal number of function evaluations reached";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_LS_FAILURE:
        {
            msg = "line search failed";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_TR_REDUCTION:
        {
            msg = "trust region failed";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_DIVERGED_USER:
        {
            msg = "user defined divergence criterion met";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
        case TAO_CONTINUE_ITERATING:
        {
            // display nothing
            break;
        }
        default:
        {
            msg = "convergence reason not defined";
            ierr = WrngMsg(msg); CHKERRQ(ierr);
            break;
        }
    }

    PetscFunctionReturn(ierr);
}




}   // namespace reg




#endif  // _TAOINTERFACE_CPP_
