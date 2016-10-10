/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _OPTIMIZER_CPP_
#define _OPTIMIZER_CPP_

#include <string>

#include "Optimizer.hpp"
#include "KrylovInterfaceReg.hpp"
#include "TaoInterfaceRegistration.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Optimizer"
Optimizer::Optimizer() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~Optimizer"
Optimizer::~Optimizer(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 * @param opt base class for registration options and arguments
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Optimizer"
Optimizer::Optimizer(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode Optimizer::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Tao = NULL;
    this->m_Precond = NULL;
    this->m_PreProc = NULL;
    this->m_Solution = NULL;
    this->m_KrylovMethod = NULL;
    this->m_OptimizationProblem = NULL;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
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

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the initial guess (warm start)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetInitialGuess"
PetscErrorCode Optimizer::SetInitialGuess(VecField* x) {
    PetscErrorCode ierr = 0;
    IntType nlu, ngu;

    PetscFunctionBegin;

    if (this->m_Solution == NULL) {
        // compute the number of unknowns
        nlu = 3*this->m_Opt->GetDomainPara().nlocal;
        ngu = 3*this->m_Opt->GetDomainPara().nglobal;

        ierr = VecCreate(this->m_Solution, nlu, ngu); CHKERRQ(ierr);
        ierr = VecSet(this->m_Solution, 0.0); CHKERRQ(ierr);
    }

    // the input better is not zero
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = x->GetComponents(this->m_Solution); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief parse initial guess to tao
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetInitialGuess"
PetscErrorCode Optimizer::SetInitialGuess() {
    PetscErrorCode ierr = 0;
    IntType nlu, ngu;

    PetscFunctionBegin;

    // check if tao has been set up
    ierr = Assert(this->m_Tao != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Solution == NULL) {
        nlu = 3*this->m_Opt->GetDomainPara().nlocal;
        ngu = 3*this->m_Opt->GetDomainPara().nglobal;
        ierr = VecCreate(this->m_Solution, nlu, ngu); CHKERRQ(ierr);
        ierr = VecSet(this->m_Solution, 0.0); CHKERRQ(ierr);
    }

    // parse initial guess to tao
    ierr = TaoSetInitialVector(this->m_Tao, this->m_Solution); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them; currently this is only supported for
 * registration problems)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetProblem"
PetscErrorCode Optimizer::SetProblem(Optimizer::OptProbType* optprob) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_OptimizationProblem = optprob;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the preconditioner
 * @param[in] precond interface to preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetPreconditioner"
PetscErrorCode Optimizer::SetPreconditioner(PrecondReg* precond) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(precond != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_Precond = precond;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set up optimization problem (this function sets up
 * and parses the default parameters to TAO)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupTao"
PetscErrorCode Optimizer::SetupTao() {
    PetscErrorCode ierr = 0;
    IntType nlu, ngu;
    ScalarType gatol, grtol, gttol, reltol, abstol, divtol;
    IntType maxit;
    Mat matvec;
    PC preconditioner;
    TaoLineSearch linesearch;
    PetscFunctionBegin;

    ierr = Assert(this->m_OptimizationProblem !=NULL, "optimization problem not set"); CHKERRQ(ierr);

    // compute the number of unknowns
    nlu = 3*this->m_Opt->GetDomainPara().nlocal;
    ngu = 3*this->m_Opt->GetDomainPara().nglobal;

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
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
            ierr = PetscOptionsSetValue(NULL, "-tao_nls_pc_type", "petsc"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue(NULL, "-tao_nls_ksp_type", "petsc"); CHKERRQ(ierr);
#else
            ierr = PetscOptionsSetValue("-tao_nls_pc_type", "petsc"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue("-tao_nls_ksp_type", "petsc"); CHKERRQ(ierr);
#endif
            ierr = TaoSetFromOptions(this->m_Tao); CHKERRQ(ierr);
        } else if (strcmp(method.c_str(), "ntr") == 0) {
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
            ierr = PetscOptionsSetValue(NULL, "-tao_ntr_pc_type", "petsc"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue(NULL, "-tao_ntr_ksp_type", "petsc"); CHKERRQ(ierr);
#else
            ierr = PetscOptionsSetValue("-tao_ntr_pc_type", "petsc"); CHKERRQ(ierr);
            ierr = PetscOptionsSetValue("-tao_ntr_ksp_type", "petsc"); CHKERRQ(ierr);
#endif
            ierr = TaoSetFromOptions(this->m_Tao); CHKERRQ(ierr);
        }

        // set tolerances for krylov subspace method
        reltol = this->m_Opt->GetKrylovSolverPara().tol[0];     // 1E-12;
        abstol = this->m_Opt->GetKrylovSolverPara().tol[1];     // 1E-12;
        divtol = this->m_Opt->GetKrylovSolverPara().tol[2];     // 1E+06;
        maxit  = this->m_Opt->GetKrylovSolverPara().maxit;      // 1000;
        ierr = KSPSetTolerances(this->m_KrylovMethod, reltol, abstol, divtol, maxit); CHKERRQ(ierr);
        ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethod, PETSC_FALSE); CHKERRQ(ierr);
//        ierr = KSPSetInitialGuessNonzero(krylovmethod,PETSC_TRUE); CHKERRQ(ierr);

        // KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
        // KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
        // KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
        ierr = KSPSetNormType(this->m_KrylovMethod, KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
//        ierr = KSPSetNormType(this->m_KrylovMethod,KSP_NORM_PRECONDITIONED); CHKERRQ(ierr);

        // set the kylov method
        if (this->m_Opt->GetKrylovSolverPara().solver == GMRES) {
            ierr = KSPSetType(this->m_KrylovMethod, KSPGMRES); CHKERRQ(ierr);
        } else if (this->m_Opt->GetKrylovSolverPara().solver == PCG) {
            ierr = KSPSetType(this->m_KrylovMethod, KSPCG); CHKERRQ(ierr);
        } else {
            ierr = ThrowError("interface for solver not provided"); CHKERRQ(ierr);
        }

        // apply projection operator to gradient and solution
        ierr = KSPSetPostSolve(this->m_KrylovMethod, PostKrylovSolve, this->m_OptimizationProblem); CHKERRQ(ierr);
        ierr = KSPSetPreSolve(this->m_KrylovMethod, PreKrylovSolve, this->m_OptimizationProblem); CHKERRQ(ierr);

        // set krylov monitor
        if (this->m_Opt->GetVerbosity() > 0) {  /// || (this->m_Opt->GetLogger()->IsEnabled(LOGKSPRES))) {
            ierr = KSPMonitorSet(this->m_KrylovMethod, KrylovMonitor, this->m_OptimizationProblem, NULL); CHKERRQ(ierr);
        }

        // set the preconditioner
        ierr = KSPGetPC(this->m_KrylovMethod, &preconditioner); CHKERRQ(ierr);
        ierr = KSPSetFromOptions(this->m_KrylovMethod); CHKERRQ(ierr);

        // switch between different preconditioners
        if (this->m_Opt->GetKrylovSolverPara().pctype == NOPC) {
            ierr = PCSetType(preconditioner, PCNONE); CHKERRQ(ierr);
        } else {
            ierr = Assert(this->m_Precond != NULL, "null pointer"); CHKERRQ(ierr);

            // we have to create a shell object for the preconditioner,
            // since our solver is matrix free
            ierr = PCSetType(preconditioner, PCSHELL); CHKERRQ(ierr);
            ierr = PCShellSetApply(preconditioner, PrecondMatVec); CHKERRQ(ierr);
            ierr = PCShellSetContext(preconditioner, this->m_Precond); CHKERRQ(ierr);
//            ierr = PCShellSetName(taokktpc,"kktpc"); CHKERRQ(ierr);
//            ierr = PCShellSetSetUp(preconditioner,PrecondSetup); CHKERRQ(ierr);
        }
    }

    // set the routine to evaluate the objective and compute the gradient
    ierr = TaoSetObjectiveRoutine(this->m_Tao, EvaluateObjective,
                                    reinterpret_cast<void*>(this->m_OptimizationProblem)); CHKERRQ(ierr);
    ierr = TaoSetGradientRoutine(this->m_Tao, EvaluateGradient,
                                    reinterpret_cast<void*>(this->m_OptimizationProblem)); CHKERRQ(ierr);
    ierr = TaoSetObjectiveAndGradientRoutine(this->m_Tao, EvaluateObjectiveGradient,
                                    reinterpret_cast<void*>(this->m_OptimizationProblem)); CHKERRQ(ierr);

    // set the monitor for the optimization process
    ierr = TaoCancelMonitors(this->m_Tao); CHKERRQ(ierr);
    ierr = TaoSetMonitor(this->m_Tao, OptimizationMonitor, this->m_OptimizationProblem, NULL); CHKERRQ(ierr);
    ierr = TaoSetConvergenceTest(this->m_Tao, CheckConvergence, this->m_OptimizationProblem); CHKERRQ(ierr);

    ierr = TaoGetLineSearch(this->m_Tao, &linesearch); CHKERRQ(ierr);
    ierr = TaoLineSearchSetType(linesearch, "armijo"); CHKERRQ(ierr);

    // set tolerances for optimizer
    gatol = this->m_Opt->GetOptPara().tol[0];   // ||g(x)||             <= gatol
    grtol = this->m_Opt->GetOptPara().tol[1];   // ||g(x)|| / |J(x)|    <= grtol
    gttol = this->m_Opt->GetOptPara().tol[2];   // ||g(x)|| / ||g(x0)|| <= gttol

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
    ierr = TaoSetTolerances(this->m_Tao, gatol, grtol, gttol); CHKERRQ(ierr);
#else
    ierr = TaoSetTolerances(this->m_Tao, 1E-12, 1E-12, gatol, grtol, gttol); CHKERRQ(ierr);
#endif
    ierr = TaoSetMaximumIterations(this->m_Tao, this->m_Opt->GetOptPara().maxit-1); CHKERRQ(ierr);
    ierr = TaoSetFunctionLowerBound(this->m_Tao, 1E-6); CHKERRQ(ierr);

    ierr = MatCreateShell(PETSC_COMM_WORLD, nlu, nlu, ngu, ngu,
                            reinterpret_cast<void*>(this->m_OptimizationProblem), &matvec); CHKERRQ(ierr);
    ierr = MatShellSetOperation(matvec, MATOP_MULT, (void(*)(void))HessianMatVec); CHKERRQ(ierr);
    ierr = MatSetOption(matvec, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(ierr);
    ierr = TaoSetHessianRoutine(this->m_Tao, matvec, matvec, EvaluateHessian,
                                reinterpret_cast<void*>(&this->m_OptimizationProblem)); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief run the optimizer (main interface; calls specific functions
 * according to user settings (parameter continuation, grid
 * continuation, scale continuation, ...))
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Run"
PetscErrorCode Optimizer::Run(bool presolve) {
    PetscErrorCode ierr;
    Vec x;
    ScalarType gtol;
    std::stringstream ss;
    PetscFunctionBegin;

    // check if optimization problem has been set
    ierr = Assert(this->m_OptimizationProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // do setup
    if (this->m_Tao == NULL) {
        ierr = this->SetupTao(); CHKERRQ(ierr);
    }
    ierr = Assert(this->m_Tao != NULL, "null pointer"); CHKERRQ(ierr);

    // modify tolerance if requestged ||g(x)|| / ||g(x0)|| <= gttol
    if (presolve) {
        gtol = this->m_Opt->GetOptPara().presolvetol[2];
        if (this->m_Opt->GetVerbosity() > 1) {
            ss << "presolve: relative gradient tolerance: " << std::scientific << gtol;
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
    } else {
        gtol = this->m_Opt->GetOptPara().tol[2];
    }

    // set tolerance
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
    ierr = TaoSetTolerances(this->m_Tao, PETSC_DEFAULT, PETSC_DEFAULT, gtol); CHKERRQ(ierr);
#else
    ierr = TaoSetTolerances(this->m_Tao, PETSC_DEFAULT,
                                         PETSC_DEFAULT,
                                         PETSC_DEFAULT,
                                         PETSC_DEFAULT, gtol); CHKERRQ(ierr);
#endif
    // set initial guess
    ierr = this->SetInitialGuess(); CHKERRQ(ierr);
    ierr = TaoSetUp(this->m_Tao); CHKERRQ(ierr);

    // in case we call the optimizer/solver several times
    // we have to make sure that the preconditioner is reset
    ierr = this->m_Precond->Reset(); CHKERRQ(ierr);

    // solve optimization problem
    ierr = this->m_Opt->StartTimer(T2SEXEC); CHKERRQ(ierr);
    ierr = TaoSolve(this->m_Tao); CHKERRQ(ierr);
    ierr = this->m_Opt->StopTimer(T2SEXEC); CHKERRQ(ierr);

    // get solution
    ierr = TaoGetSolutionVector(this->m_Tao, &x); CHKERRQ(ierr);

    // copy solution into place holder
    ierr = VecCopy(x, this->m_Solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief get the solution
 * @param x vector to hold solution
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetSolution"
PetscErrorCode Optimizer::GetSolution(Vec &x) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // check if we have solved the problem / set up tao
    ierr = Assert(this->m_Tao != NULL, "null pointer"); CHKERRQ(ierr);

    // get solution
    ierr = TaoGetSolutionVector(this->m_Tao, &x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief get the solution status (convergence reason)
 * @param get solution status
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetSolutionStatus"
PetscErrorCode Optimizer::GetSolutionStatus(bool &converged) {
    PetscErrorCode ierr = 0;
    TaoConvergedReason reason;
    PetscFunctionBegin;

    // check if we have solved the problem / set up tao
    ierr = Assert(this->m_Tao != NULL, "null pointer"); CHKERRQ(ierr);

    // get solution
    ierr = TaoGetConvergedReason(this->m_Tao, &reason); CHKERRQ(ierr);

    converged = true;
    if (reason < 0) {
        converged = false;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief finalize optimization (displays information for user)
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Finalize"
PetscErrorCode Optimizer::Finalize() {
    PetscErrorCode ierr = 0;
    int rank, indent, numindent, linelength;
    IntType maxiter, iter;
    ScalarType gatol, grtol, gttol, gnorm, J, g0norm;
    bool stop[3], converged;
    std::string line, msg;
    TaoConvergedReason reason;
    std::stringstream ss;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = Assert(this->m_Tao !=NULL, "null pointer"); CHKERRQ(ierr);

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
    ierr = TaoGetTolerances(this->m_Tao, &gatol, &grtol, &gttol); CHKERRQ(ierr);
#else
    ierr = TaoGetTolerances(this->m_Tao, NULL, NULL, &gatol, &grtol, &gttol); CHKERRQ(ierr);
#endif
    g0norm = this->m_OptimizationProblem->GetInitialGradNorm();

    ierr = TaoGetMaximumIterations(this->m_Tao, &maxiter); CHKERRQ(ierr);
    ierr = TaoGetSolutionStatus(this->m_Tao, &iter, &J, &gnorm, NULL, NULL, &reason); CHKERRQ(ierr);

    linelength = this->m_Opt->GetLineLength();
    line = std::string(linelength, '-');

    stop[0] = (gnorm < gttol*g0norm);   ///< relative change of gradient
    stop[1] = (gnorm < gatol);          ///< absolute norm of gradient
    stop[2] = (iter  > maxiter);        ///< max number of iterations

    // check if we converged
    converged = false;
    for (int i = 0; i < 3; ++i) {
        if (stop[i]) converged=true;
    }

    indent    = 25;
    numindent = 5;
    if (rank == 0) {
        if (converged) {
            std::cout << std::endl;
            std::cout << " convergence criteria" << std::endl;
            std::cout << std::endl;

            // relative change of gradient
            ss  << "[  " << stop[0] << "    ||g|| = " << std::setw(14)
                << std::right << std::scientific << gnorm << "    <    "
                << std::left << std::setw(14) << gttol*g0norm << " = " << "tol";
            std::cout << std::left << std::setw(100) << ss.str() << "]" << std::endl;
            ss.str(std::string()); ss.clear();

            // absolute norm of gradient
            ss  << "[  " << stop[1] << "    ||g|| = " << std::setw(14)
                << std::right << std::scientific << gnorm << "    <    "
                << std::left << std::setw(14) << gatol << " = " << "tol";
            std::cout << std::left << std::setw(100) << ss.str() << "]" << std::endl;
            ss.str(std::string()); ss.clear();

            // number of iterations
            ss  << "[  " << stop[2] << "     iter = " << std::setw(14)
                << std::right << iter  << "    >    "
                << std::left << std::setw(14) << maxiter << " = " << "maxiter";
            std::cout << std::left << std::setw(100) << ss.str() << "]" << std::endl;
            ss.str(std::string()); ss.clear();
            std::cout << std::endl;
        }
    }

    if (!converged) {
        switch (reason) {
            case TAO_CONVERGED_STEPTOL:
            {
                msg = "line search failed";
                break;
            }
            case TAO_CONVERGED_MINF:
            {
                msg = "objective value to small";
                break;
            }
            case TAO_DIVERGED_NAN:
            {
                msg = "numerical issues (NaN detected)";
                break;
            }
            case TAO_DIVERGED_MAXFCN:
            {
                msg = "maximal number of function evaluations reached";
                break;
            }
            case TAO_DIVERGED_LS_FAILURE:
            {
                msg = "line search failed";
                break;
            }
            case TAO_DIVERGED_TR_REDUCTION:
            {
                msg = "trust region failed";
                break;
            }
            default:
            {
                msg = "did not converge; reason not defined";
                break;
            }
        }
        ierr = WrngMsg(msg); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 1) {
        if (converged && !rank) std::cout<< line <<std::endl;
        ss << std::left << std::setw(indent)
           << "outer iterations" << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(ITERATIONS) - 1;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ss << std::left << std::setw(indent)
           << "objective evals" << std::right << std::setw(numindent)
           << this->m_Opt->GetCounter(OBJEVAL);
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ss << std::left << std::setw(indent)
           << "hessian matvecs" << std::right << std::setw(numindent)
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

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif
