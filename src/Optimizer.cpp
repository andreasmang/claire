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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _OPTIMIZER_CPP_
#define _OPTIMIZER_CPP_

#include <string>

#include "Optimizer.hpp"
#include "KrylovInterface.hpp"
#include "TaoInterface.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
Optimizer::Optimizer() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
Optimizer::~Optimizer(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 * @param opt base class for registration options and arguments
 *******************************************************************/
Optimizer::Optimizer(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Optimizer::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Tao = NULL;
    this->m_Precond = NULL;

    this->m_MatVec = NULL;
    this->m_PreProc = NULL;
    this->m_Solution = NULL;

    this->m_KrylovMethod = NULL;
    this->m_OptimizationProblem = NULL;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Optimizer::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // clean up tao
    if (this->m_Tao != NULL) {
        ierr = TaoDestroy(&this->m_Tao); CHKERRQ(ierr);
        this->m_Tao = NULL;
    }

    // delete solution vector
    if (this->m_Solution != NULL) {
        ierr = VecDestroy(&this->m_Solution); CHKERRQ(ierr);
        this->m_Solution = NULL;
    }

    if (this->m_MatVec != NULL) {
        ierr = MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the initial guess (warm start)
 *******************************************************************/
PetscErrorCode Optimizer::SetInitialGuess(VecField* x) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    ScalarType value;
    IntType nlu, ngu;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    if (this->m_Solution == NULL) {
        // compute number of unknowns
        nlu = 3*this->m_Opt->m_Domain.nl;
        ngu = 3*this->m_Opt->m_Domain.ng;
        ierr = VecCreate(this->m_Solution, nlu, ngu); CHKERRQ(ierr);
        ierr = VecSet(this->m_Solution, 0.0); CHKERRQ(ierr);
    }

    // the input better is not zero
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = x->GetComponents(this->m_Solution); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 1) {
        ierr = VecNorm(this->m_Solution, NORM_2, &value); CHKERRQ(ierr);
        ss << "norm of initial guess " << std::scientific << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief parse initial guess to tao
 *******************************************************************/
PetscErrorCode Optimizer::SetInitialGuess() {
    PetscErrorCode ierr = 0;
    IntType nlu, ngu;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    // check if tao has been set up
    ierr = Assert(this->m_Tao != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Solution == NULL) {
        nlu = 3*this->m_Opt->m_Domain.nl;
        ngu = 3*this->m_Opt->m_Domain.ng;
        ierr = VecCreate(this->m_Solution, nlu, ngu); CHKERRQ(ierr);
        ierr = VecSet(this->m_Solution, 0.0); CHKERRQ(ierr);
    }

    // parse initial guess to tao
    ierr = TaoSetInitialVector(this->m_Tao, this->m_Solution); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them; currently this is only supported for
 * registration problems)
 *******************************************************************/
PetscErrorCode Optimizer::SetProblem(Optimizer::OptProbType* optprob) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_OptimizationProblem = optprob;

    this->m_Opt->Exit(__func__);
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the preconditioner
 * @param[in] precond interface to preconditioner
 *******************************************************************/
PetscErrorCode Optimizer::SetPreconditioner(Preconditioner* precond) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    ierr = Assert(precond != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_Precond = precond;

    this->m_Opt->Exit(__func__);
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set up optimization problem (this function sets up
 * and parses the default parameters to TAO)
 *******************************************************************/
PetscErrorCode Optimizer::SetupTao() {
    PetscErrorCode ierr = 0;
    IntType nlu, ngu;
    ScalarType gatol, grtol, gttol, reltol, abstol, divtol;
    IntType maxit;
    void* optprob;
    PC preconditioner;
    TaoLineSearch linesearch;
    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_OptimizationProblem !=NULL, "optimization problem not set"); CHKERRQ(ierr);

    // compute the number of unknowns
    nlu = 3*this->m_Opt->m_Domain.nl;
    ngu = 3*this->m_Opt->m_Domain.ng;

    // if tao exists, kill it
    if (this->m_Tao != NULL) {
        ierr = TaoDestroy(&this->m_Tao); CHKERRQ(ierr);
        this->m_Tao = NULL;
    }

    std::string method = "nls";
    ierr = TaoCreate(PETSC_COMM_WORLD, &this->m_Tao); CHKERRQ(ierr);
    ierr = TaoSetType(this->m_Tao, "nls"); CHKERRQ(ierr);

    // get the ksp of the optimizer and set options
    ierr = TaoGetKSP(this->m_Tao, &this->m_KrylovMethod); CHKERRQ(ierr);

    // ksp is only nonzero if we use a newton type method
    if (this->m_KrylovMethod != NULL) {
        // switch of the standard preconditioner
        if (strcmp(method.c_str(), "nls") == 0) {
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7) && (PETSC_VERSION_MINOR < 10)
            ierr = PetscOptionsSetValue(NULL, "-tao_nls_pc_type", "petsc"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue(NULL, "-tao_nls_ksp_type", "petsc"); CHKERRQ(ierr);
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 10)
            ierr = PetscOptionsSetValue(NULL, "-tao_nls_pc_type", "none"); CHKERRQ(ierr);
            // TODO: CHECK IF WE CAN USE OTHER TYPES
            if (this->m_Opt->m_KrylovMethod.solver == PCG) {
                ierr = PetscOptionsSetValue(NULL, "-tao_nls_ksp_type", "cg"); CHKERRQ(ierr);
            } else if (this->m_Opt->m_KrylovMethod.solver == GMRES) {
                ierr = PetscOptionsSetValue(NULL, "-tao_nls_ksp_type", "gmres"); CHKERRQ(ierr);
            }
#else
            ierr = PetscOptionsSetValue("-tao_nls_pc_type", "petsc"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue("-tao_nls_ksp_type", "petsc"); CHKERRQ(ierr);
#endif
            ierr = TaoSetFromOptions(this->m_Tao); CHKERRQ(ierr);
        } else if (strcmp(method.c_str(), "ntr") == 0) {
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7) && (PETSC_VERSION_MINOR < 10)
            ierr = PetscOptionsSetValue(NULL, "-tao_ntr_pc_type", "petsc"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue(NULL, "-tao_ntr_ksp_type", "petsc"); CHKERRQ(ierr);
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 10)
            // TODO: CHECK IF WE CAN USE OTHER TYPES
            ierr = PetscOptionsSetValue(NULL, "-tao_ntr_pc_type", "none"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue(NULL, "-tao_ntr_ksp_type", "stcg"); CHKERRQ(ierr); // TODO: TEST ME
#else
            ierr = PetscOptionsSetValue("-tao_ntr_pc_type", "petsc"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue("-tao_ntr_ksp_type", "petsc"); CHKERRQ(ierr);
#endif
            ierr = TaoSetFromOptions(this->m_Tao); CHKERRQ(ierr);
        }

        // set tolerances for krylov subspace method
        reltol = this->m_Opt->m_KrylovMethod.tol[0];     // 1E-12;
        abstol = this->m_Opt->m_KrylovMethod.tol[1];     // 1E-12;
        divtol = this->m_Opt->m_KrylovMethod.tol[2];     // 1E+06;
        maxit  = this->m_Opt->m_KrylovMethod.maxiter;    // 1000;
        maxit  = std::max(static_cast<IntType>(0), maxit-1);
        ierr = KSPSetTolerances(this->m_KrylovMethod, reltol, abstol, divtol, maxit); CHKERRQ(ierr);
        ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethod, PETSC_FALSE); CHKERRQ(ierr);
//        ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethod, PETSC_TRUE); CHKERRQ(ierr);

        // KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
        // KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
        // KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
        ierr = KSPSetNormType(this->m_KrylovMethod, KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
//        ierr = KSPSetNormType(this->m_KrylovMethod,KSP_NORM_PRECONDITIONED); CHKERRQ(ierr);
//        ierr = KSPSetNormType(this->m_KrylovMethod,KSP_NORM_NATURAL); CHKERRQ(ierr);

        // set the kylov method
        if (this->m_Opt->m_KrylovMethod.solver == GMRES) {
            ierr = KSPSetType(this->m_KrylovMethod, KSPGMRES); CHKERRQ(ierr);
        } else if (this->m_Opt->m_KrylovMethod.solver == PCG) {
            ierr = KSPSetType(this->m_KrylovMethod, KSPCG); CHKERRQ(ierr);
        } else if (this->m_Opt->m_KrylovMethod.solver == FCG) {
            ierr = KSPSetType(this->m_KrylovMethod, KSPFCG); CHKERRQ(ierr);
        } else if (this->m_Opt->m_KrylovMethod.solver == FGMRES) {
            ierr = KSPSetType(this->m_KrylovMethod, KSPFGMRES); CHKERRQ(ierr);
        } else {
            ierr = ThrowError("interface for solver not provided"); CHKERRQ(ierr);
        }

        // apply projection operator to gradient and
        // solution if needed (two-level preconditioner)
        ierr = KSPSetPostSolve(this->m_KrylovMethod, PostKrylovSolve, this->m_OptimizationProblem); CHKERRQ(ierr);
        ierr = KSPSetPreSolve(this->m_KrylovMethod, PreKrylovSolve, this->m_OptimizationProblem); CHKERRQ(ierr);

        // set krylov monitor
        if (this->m_Opt->m_Verbosity > 0) {  /// || (this->m_Opt->GetLogger()->IsEnabled(LOGKSPRES))) {
            ierr = KSPMonitorSet(this->m_KrylovMethod, KrylovMonitor, this->m_OptimizationProblem, NULL); CHKERRQ(ierr);
        }

        // set the preconditioner
        ierr = KSPGetPC(this->m_KrylovMethod, &preconditioner); CHKERRQ(ierr);
        ierr = KSPSetFromOptions(this->m_KrylovMethod); CHKERRQ(ierr);

        // switch between different preconditioners
        if (this->m_Opt->m_KrylovMethod.pctype == NOPC) {
            ierr = PCSetType(preconditioner, PCNONE); CHKERRQ(ierr);
        } else {
            ierr = Assert(this->m_Precond != NULL, "null pointer"); CHKERRQ(ierr);

            // we have to create a shell object for the preconditioner,
            // since our solver is matrix free
            ierr = PCSetType(preconditioner, PCSHELL); CHKERRQ(ierr);
            ierr = PCShellSetApply(preconditioner, PrecondMatVec); CHKERRQ(ierr);
            ierr = PCShellSetContext(preconditioner, this->m_Precond); CHKERRQ(ierr);
//            ierr = PCShellSetName(taokktpc,"kktpc"); CHKERRQ(ierr);
//            ierr = PCShellSetSetUp(preconditioner, PrecondSetup); CHKERRQ(ierr);
        }
    }

    optprob = reinterpret_cast<void*>(this->m_OptimizationProblem);

    // set the routine to evaluate the objective and compute the gradient
    ierr = TaoSetObjectiveRoutine(this->m_Tao, EvaluateObjective, optprob); CHKERRQ(ierr);
    ierr = TaoSetGradientRoutine(this->m_Tao, EvaluateGradient, optprob); CHKERRQ(ierr);
    ierr = TaoSetObjectiveAndGradientRoutine(this->m_Tao, EvaluateObjectiveGradient, optprob); CHKERRQ(ierr);

    // set the monitor for the optimization process
    ierr = TaoCancelMonitors(this->m_Tao); CHKERRQ(ierr);
    ierr = TaoSetMonitor(this->m_Tao, OptimizationMonitor, this->m_OptimizationProblem, NULL); CHKERRQ(ierr);

    // set function to test stopping conditions
    if (this->m_Opt->m_OptPara.stopcond == GRAD) {
        ierr = TaoSetConvergenceTest(this->m_Tao, CheckConvergenceGrad, this->m_OptimizationProblem); CHKERRQ(ierr);
    } else if (this->m_Opt->m_OptPara.stopcond == GRADOBJ) {
        ierr = TaoSetConvergenceTest(this->m_Tao, CheckConvergenceGradObj, this->m_OptimizationProblem); CHKERRQ(ierr);
    } else {
        ierr = ThrowError("stop condition not defined"); CHKERRQ(ierr);
    }

    ierr = TaoGetLineSearch(this->m_Tao, &linesearch); CHKERRQ(ierr);

    switch(this->m_Opt->m_OptPara.glmethod) {
        case NOGM:
        {
            ierr = TaoLineSearchSetType(linesearch, "unit"); CHKERRQ(ierr);
            break;
        }
        case ARMIJOLS:
        {
            ierr = TaoLineSearchSetType(linesearch, "armijo"); CHKERRQ(ierr);
            break;
        }
        case OWARMIJOLS:
        {
            ierr = TaoLineSearchSetType(linesearch, "owarmijo"); CHKERRQ(ierr);
            break;
        }
        case MTLS:
        {
            ierr = TaoLineSearchSetType(linesearch, "more-thuente"); CHKERRQ(ierr);
            break;
        }
        case GPCGLS:
        {
            ierr = TaoLineSearchSetType(linesearch, "gpcg"); CHKERRQ(ierr);
            break;
        }
        case IPMLS:
        {
            ierr = TaoLineSearchSetType(linesearch, "ipm"); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("globalization method not defined"); CHKERRQ(ierr);
        }
    }

    //ierr = TaoLineSearchView(linesearch, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);

    // set tolerances for optimizer
    gatol = this->m_Opt->m_OptPara.tol[0];   // ||g(x)||             <= gatol
    grtol = this->m_Opt->m_OptPara.tol[1];   // ||g(x)|| / |J(x)|    <= grtol
    gttol = this->m_Opt->m_OptPara.tol[2];   // ||g(x)|| / ||g(x0)|| <= gttol

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
    ierr = TaoSetTolerances(this->m_Tao, gatol, grtol, gttol); CHKERRQ(ierr);
#else
    ierr = TaoSetTolerances(this->m_Tao, 1E-12, 1E-12, gatol, grtol, gttol); CHKERRQ(ierr);
#endif
    ierr = TaoSetMaximumIterations(this->m_Tao, this->m_Opt->m_OptPara.maxiter-1); CHKERRQ(ierr);
    ierr = TaoSetFunctionLowerBound(this->m_Tao, 1E-6); CHKERRQ(ierr);

    ierr = MatCreateShell(PETSC_COMM_WORLD, nlu, nlu, ngu, ngu, optprob, &this->m_MatVec); CHKERRQ(ierr);
    ierr = MatShellSetOperation(this->m_MatVec, MATOP_MULT, (void(*)(void))HessianMatVec); CHKERRQ(ierr);
    ierr = MatSetOption(this->m_MatVec, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(ierr);
    ierr = TaoSetHessianRoutine(this->m_Tao, this->m_MatVec, this->m_MatVec, EvaluateHessian, optprob); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);
    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief run the optimizer (main interface; calls specific functions
 * according to user settings (parameter continuation, grid
 * continuation, scale continuation, ...))
 *******************************************************************/
PetscErrorCode Optimizer::Run(bool presolve) {
    PetscErrorCode ierr;
    ScalarType gtol;
    IntType maxit;
    std::stringstream ss;
    Vec x = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // check if optimization problem has been set
    ierr = Assert(this->m_OptimizationProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // do setup
    if (this->m_Tao == NULL) {
        ierr = this->SetupTao(); CHKERRQ(ierr);
    }
    ierr = Assert(this->m_Tao != NULL, "null pointer"); CHKERRQ(ierr);

    // modify tolerance if requestged ||g(x)|| / ||g(x0)|| <= gttol
    if (presolve) {
        gtol = this->m_Opt->m_OptPara.presolvetol[2];
        maxit = this->m_Opt->m_OptPara.presolvemaxit;
        if (this->m_Opt->m_Verbosity > 1) {
            ss << "presolve: relative gradient tolerance: " << std::scientific << gtol;
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
    } else {
        gtol = this->m_Opt->m_OptPara.tol[2];
        maxit = this->m_Opt->m_OptPara.maxiter;
    }

    // set tolerance
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
    ierr = TaoSetTolerances(this->m_Tao, PETSC_DEFAULT, PETSC_DEFAULT, gtol); CHKERRQ(ierr);
#else
    ierr = TaoSetTolerances(this->m_Tao, PETSC_DEFAULT, PETSC_DEFAULT,
                                         PETSC_DEFAULT, PETSC_DEFAULT, gtol); CHKERRQ(ierr);
#endif
    ierr = TaoSetMaximumIterations(this->m_Tao, maxit-1); CHKERRQ(ierr);

    // set initial guess
    ierr = this->SetInitialGuess(); CHKERRQ(ierr);
    ierr = TaoSetUp(this->m_Tao); CHKERRQ(ierr);

    if (this->m_Opt->m_KrylovMethod.pctype != NOPC) {
        // in case we call the optimizer/solver several times
        // we have to make sure that the preconditioner is reset
        ierr = this->m_Precond->Reset(); CHKERRQ(ierr);
    }

    // solve optimization problem
    ierr = this->m_Opt->StartTimer(T2SEXEC); CHKERRQ(ierr);
    ierr = TaoSolve(this->m_Tao); CHKERRQ(ierr);
    ierr = this->m_Opt->StopTimer(T2SEXEC); CHKERRQ(ierr);

    // get solution
    ierr = TaoGetSolutionVector(this->m_Tao, &x); CHKERRQ(ierr);

    // copy solution into place holder
    ierr = VecCopy(x, this->m_Solution); CHKERRQ(ierr);

    // if we did not converge, we should display what's going on
    if (!this->m_OptimizationProblem->Converged()) {
        ierr = OptimizationMonitor(this->m_Tao, this->m_OptimizationProblem); CHKERRQ(ierr);
    }



    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief get the solution
 * @param x vector to hold solution
 *******************************************************************/
PetscErrorCode Optimizer::GetSolution(Vec &x) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // check if we have solved the problem / set up tao
    ierr = Assert(this->m_Tao != NULL, "null pointer"); CHKERRQ(ierr);

    // get solution
    ierr = TaoGetSolutionVector(this->m_Tao, &x); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief get the solution status (convergence reason)
 * @param get solution status
 *******************************************************************/
PetscErrorCode Optimizer::GetSolutionStatus(bool &converged) {
    PetscErrorCode ierr = 0;
    TaoConvergedReason reason;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // check if we have solved the problem / set up tao
    ierr = Assert(this->m_Tao != NULL, "null pointer"); CHKERRQ(ierr);

    // get solution
    ierr = TaoGetConvergedReason(this->m_Tao, &reason); CHKERRQ(ierr);

    converged = true;
    if (reason < 0) {
        if (reason != TAO_DIVERGED_MAXITS) {
            converged = false;
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief finalize optimization (displays information for user)
 ********************************************************************/
PetscErrorCode Optimizer::Finalize() {
    PetscErrorCode ierr = 0;
    int rank, indent, numindent, linelength;
    std::string line;
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << " convergence criteria" << std::endl;
        std::cout << std::endl;
        std::cout << this->m_OptimizationProblem->GetConvergenceMessage();
    }

    linelength = this->m_Opt->m_LineLength;
    line = std::string(linelength, '-');
    indent = 35; numindent = 5;
    if (this->m_Opt->m_Verbosity > 1) {
        if (this->m_OptimizationProblem->Converged() && !rank) {
            std::cout << line << std::endl;
        }
        ss << std::left << std::setw(indent)
           << "newton iterations" << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(ITERATIONS) - 1;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ss << std::left << std::setw(indent)
           << "objective evaluations" << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(OBJEVAL);
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ss << std::left << std::setw(indent)
           << "hessian matrix vector products" << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(HESSMATVEC);
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ss << std::left << std::setw(indent)
           << "pde solves" << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(PDESOLVE);
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // display info to user, once we're done
//    ierr = TaoView(this->m_Tao, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif
