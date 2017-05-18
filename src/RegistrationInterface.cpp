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
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGISTRATIONINTERFACE_CPP_
#define _REGISTRATIONINTERFACE_CPP_

#include "RegistrationInterface.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegistrationInterface::RegistrationInterface() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegistrationInterface::~RegistrationInterface(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 * @param[in] opt base class for registration options and arguments
 *******************************************************************/
RegistrationInterface::RegistrationInterface(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode RegistrationInterface::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;
    this->m_PreProc = NULL;
    this->m_ReadWrite = NULL;
    this->m_Optimizer = NULL;
    this->m_Precond = NULL;
    this->m_RegProblem = NULL;

    this->m_TemplateImage = NULL;
    this->m_ReferenceImage = NULL;

    this->m_TemplatePyramid = NULL;
    this->m_ReferencePyramid = NULL;

    this->m_Solution = NULL;

    this->m_IsTemplateSet = false;
    this->m_IsReferenceSet = false;

    this->m_DeleteSolution = true;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode RegistrationInterface::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // delete class for registration problem
    if (this->m_RegProblem != NULL) {
        delete this->m_RegProblem;
        this->m_RegProblem = NULL;
    }

    // delete class for optimizer
    if (this->m_Optimizer != NULL) {
        delete this->m_Optimizer;
        this->m_Optimizer = NULL;
    }

    // delete class for pre-processing
    if (this->m_PreProc != NULL) {
        delete this->m_PreProc;
        this->m_PreProc = NULL;
    }

    // delete class for preconditioner
    if (this->m_Precond != NULL) {
        delete this->m_Precond;
        this->m_Precond = NULL;
    }

    if (this->m_DeleteSolution) {
        if (this->m_Solution != NULL) {
            delete this->m_Solution;
            this->m_Solution = NULL;
        }
    }

    if (this->m_ReferencePyramid != NULL) {
        delete this->m_ReferencePyramid;
        this->m_ReferencePyramid = NULL;
    }

    if (this->m_TemplatePyramid != NULL) {
        delete this->m_TemplatePyramid;
        this->m_TemplatePyramid = NULL;
    }

    // if we did not read/set the images, we can
    // destroy the containers
    if (!this->m_IsReferenceSet) {
        if (this->m_ReferenceImage != NULL) {
            ierr = VecDestroy(&this->m_ReferenceImage); CHKERRQ(ierr);
            this->m_ReferenceImage = NULL;
        }
    }

    if (!this->m_IsTemplateSet) {
        if (this->m_TemplateImage != NULL) {
            ierr = VecDestroy(&this->m_TemplateImage); CHKERRQ(ierr);
            this->m_TemplateImage = NULL;
        }
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set initial guess
 * @param[in] x    vector field that contains initial guess
 *******************************************************************/
PetscErrorCode RegistrationInterface::SetInitialGuess(VecField* x, bool copy) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // the input better is not zero
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);

    if (copy) {
        // if we have not setup initial guess, do so
        if (this->m_Solution == NULL) {
            try {this->m_Solution = new VecField(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        ierr = this->m_Solution->Copy(x); CHKERRQ(ierr);
        this->m_DeleteSolution = true;
    } else {
        this->m_Solution = x;
        this->m_DeleteSolution = false;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get solution of optimization problem
 * @param[out] x    vector field for solution; needs to be allocated
 * elswhere
 *******************************************************************/
PetscErrorCode RegistrationInterface::GetSolution(VecField* x, bool copy) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // the input better is not zero
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_Solution != NULL, "null pointer"); CHKERRQ(ierr);

    if (copy) {
        ierr = x->Copy(this->m_Solution); CHKERRQ(ierr);
    } else {
        x = this->m_Solution;
    }
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
PetscErrorCode RegistrationInterface::SetReadWrite(ReadWriteReg* readwrite) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(readwrite != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite = readwrite;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode RegistrationInterface::SetReferenceImage(Vec mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);

    this->m_ReferenceImage = mR;
    this->m_IsReferenceSet = true;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set template image
 *******************************************************************/
PetscErrorCode RegistrationInterface::SetTemplateImage(Vec mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);

    this->m_TemplateImage = mT;
    this->m_IsTemplateSet = true;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get final state at t=1 for current iterate v
 * (stored in state variable)
 *******************************************************************/
PetscErrorCode RegistrationInterface::GetFinalState(Vec m1) {
    PetscErrorCode ierr = 0;
    Vec m = NULL;
    IntType nl, nc, nt;
    ScalarType *p_m = NULL, *p_m1 = NULL;
    PetscFunctionBegin;

    ierr = Assert(m1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = this->m_RegProblem->GetStateVariable(m); CHKERRQ(ierr);

    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    nt = this->m_Opt->GetDomainPara().nt;

    ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(m1, &p_m1); CHKERRQ(ierr);

    try {std::copy(p_m+nt*nc*nl, p_m+(nt+1)*nc*nl, p_m1);}
    catch (std::exception&) {
        ierr = ThrowError("copy failed"); CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);
    ierr = VecRestoreArray(m1, &p_m1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get the residual at t=1 for current iterate v
 * (stored in adjoint variable)
 *******************************************************************/
PetscErrorCode RegistrationInterface::GetResidual(Vec residual) {
    PetscErrorCode ierr = 0;
    Vec lambda = NULL, m = NULL;
    IntType nl, nc, nt, l0, l1, l;
    ScalarType *p_l = NULL, *p_res = NULL, *p_m = NULL, *p_mr = NULL;
    PetscFunctionBegin;

    ierr = Assert(residual != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);


    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    nt = this->m_Opt->GetDomainPara().nt;

    // get array for residual
    ierr = VecGetArray(residual, &p_res); CHKERRQ(ierr);

    if (this->m_Opt->GetOptPara().method == FULLNEWTON) {
        // get lambda
        ierr = this->m_RegProblem->GetAdjointVariable(lambda); CHKERRQ(ierr);
        ierr = VecGetArray(lambda, &p_l); CHKERRQ(ierr);
        l0 =  nt*nc*nl; l1 = (nt+1)*nc*nl;
        try {std::copy(p_l+l0, p_l+l1, p_res);}
        catch (std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(lambda, &p_l); CHKERRQ(ierr);
    } else {
        l = nt*nc*nl;  // index for final condition
        // compute terminal condition \lambda_1 = -(m_1 - m_R) = m_R - m_1
        ierr = this->m_RegProblem->GetStateVariable(m); CHKERRQ(ierr);
        ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nc*nl; ++i) {
            p_l[i] = p_mr[i] - p_m[l+i];  // compute initial condition
        }
}  // omp
        ierr = VecRestoreArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
        ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);
    }


    ierr = VecRestoreArray(residual, &p_res); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read/write object
 *******************************************************************/
PetscErrorCode RegistrationInterface::DispLevelMsg(std::string msg, int rank) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (rank == 0) {
        std::cout << std::string(this->m_Opt->GetLineLength(), '-') << std::endl;
    }

    ierr = Msg(msg); CHKERRQ(ierr);

    if (rank == 0) {
        std::cout << std::string(this->m_Opt->GetLineLength(), '-') << std::endl;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evalute the regularization functional
 * @param[in] v velocity field
 * @param[out] value sobolev norm of velocity v
 *******************************************************************/
PetscErrorCode RegistrationInterface
::EvaluateRegularizationFunctional(ScalarType* value, VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // reset registration problem
    if (this->m_RegProblem == NULL) {
        ierr = this->SetupRegProblem(); CHKERRQ(ierr);
    }

    ierr = this->m_RegProblem->EvaluateRegularizationFunctional(value, v); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set up the registration problem and optimizer
 ********************************************************************/
PetscErrorCode RegistrationInterface::SetupSolver() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    ScalarType vn1, vn2, vn3;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("setting up solver"); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg(" >> allocation of optimizer"); CHKERRQ(ierr);
    }
    // reset optimizer
    if (this->m_Optimizer != NULL) {
        delete this->m_Optimizer;
        this->m_Optimizer = NULL;
    }
    try {this->m_Optimizer = new OptimizerType(this->m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // set up optimization/registration problem
    ierr = this->SetupRegProblem(); CHKERRQ(ierr);

    // reset/setup preprocessing
    if (this->m_PreProc != NULL) {
        delete this->m_PreProc; this->m_PreProc = NULL;
    }
    try {this->m_PreProc = new Preprocessing(this->m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetKrylovSolverPara().pctype != NOPC) {
        if (this->m_Opt->GetVerbosity() > 2) {
            ierr = DbgMsg("allocating preconditioner"); CHKERRQ(ierr);
        }
        // reset/setup preconditioner
        if (this->m_Precond != NULL) {
            delete this->m_Precond; this->m_Precond = NULL;
        }
        try {this->m_Precond = new PrecondReg(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_Precond->SetPreProc(this->m_PreProc); CHKERRQ(ierr);
        ierr = this->m_Precond->SetProblem(this->m_RegProblem); CHKERRQ(ierr);
        ierr = this->m_Optimizer->SetPreconditioner(this->m_Precond); CHKERRQ(ierr);
    }

    // set up initial condition
    if (this->m_Solution == NULL) {
        if (this->m_Opt->GetVerbosity() > 2) {
            ierr = DbgMsg("allocating solution vector"); CHKERRQ(ierr);
        }
        try {this->m_Solution = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_Solution->SetValue(0.0); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 1) {
        ierr = VecNorm(this->m_Solution->m_X1, NORM_2, &vn1); CHKERRQ(ierr);
        ierr = VecNorm(this->m_Solution->m_X2, NORM_2, &vn2); CHKERRQ(ierr);
        ierr = VecNorm(this->m_Solution->m_X3, NORM_2, &vn3); CHKERRQ(ierr);
        ss  << "norm of initial guess: "
            << std::scientific << "(" << vn1 << " " << vn2 << " " << vn3 << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set up the data (e.g., apply preprocessing)
 ********************************************************************/
PetscErrorCode RegistrationInterface::SetupData(Vec& mR, Vec& mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // presmoothing, if necessary
    if (this->m_IsTemplateSet && this->m_IsReferenceSet) {
        ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

        // allocate
        ierr = VecDuplicate(this->m_TemplateImage, &mT); CHKERRQ(ierr);
        ierr = VecDuplicate(this->m_ReferenceImage, &mR); CHKERRQ(ierr);

        if (this->m_Opt->GetRegFlags().applysmoothing) {
            if (this->m_PreProc == NULL) {
                try{this->m_PreProc = new Preprocessing(this->m_Opt);}
                catch (std::bad_alloc&) {
                    ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
                }
            }
            ierr = this->m_PreProc->Smooth(mR, this->m_ReferenceImage); CHKERRQ(ierr);
            ierr = this->m_PreProc->Smooth(mT, this->m_TemplateImage); CHKERRQ(ierr);
        } else {
            ierr = VecCopy(this->m_ReferenceImage, mR); CHKERRQ(ierr);
            ierr = VecCopy(this->m_TemplateImage, mT); CHKERRQ(ierr);
        }
        ierr = this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr = this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);
    } else {
        // set up synthetic test problem
        ierr = this->m_RegProblem->SetupSyntheticProb(this->m_ReferenceImage, this->m_TemplateImage); CHKERRQ(ierr);
        ierr = this->m_RegProblem->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
        ierr = this->m_RegProblem->SetTemplateImage(this->m_TemplateImage); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set up the registration problem, which essentially is
 * equivalent to allocating the class
 ********************************************************************/
PetscErrorCode RegistrationInterface::SetupRegProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // reset registration problem
    if (this->m_RegProblem != NULL) {
        delete this->m_RegProblem; this->m_RegProblem = NULL;
    }

    // allocate class for registration
    if (this->m_Opt->GetRegModel() == COMPRESSIBLE) {
        try {this->m_RegProblem = new OptimalControlRegistration(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    } else if (this->m_Opt->GetRegModel() == STOKES) {
        try {this->m_RegProblem = new OptimalControlRegistrationIC(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    } else if (this->m_Opt->GetRegModel() == RELAXEDSTOKES) {
        try {this->m_RegProblem = new OptimalControlRegistrationRelaxedIC(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    } else {
        ierr = ThrowError("registration model not available"); CHKERRQ(ierr);
    }

    ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = this->m_RegProblem->SetReadWrite(this->m_ReadWrite); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief main function to call in order to solve the optimization
 * problem
 ********************************************************************/
PetscErrorCode RegistrationInterface::Run() {
    PetscErrorCode ierr = 0;
    IntType nxmax, nx;
    std::stringstream ss;
    int rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = Msg("starting optimization"); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(), '-') << std::endl;
    ierr = PetscPrintf(PETSC_COMM_WORLD," %s  %-20s %-20s %-20s %-20s %-20s\n",
                       "iter", "objective (rel)", "mismatch (rel)",
                       "||gradient||_2,rel", "||gradient||_2", "step"); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(), '-') << std::endl;

    // switch between solvers we have to solve optimization problem
    if (this->m_Opt->GetParaCont().enabled) {
        ierr = this->RunSolverRegParaCont(); CHKERRQ(ierr);
    } else if (this->m_Opt->GetScaleContPara().enabled) {
        ierr = this->RunSolverScaleCont(); CHKERRQ(ierr);
    } else if (this->m_Opt->GetGridContPara().enabled) {
        nxmax = PETSC_MIN_INT;
        for (int i = 0; i < 3; ++i) {
            nx = this->m_Opt->GetDomainPara().nx[i];
            nxmax = nx > nxmax ? nx : nxmax;
        }

        // run grid continuation
        if (nxmax >= 32) {
            ierr = this->RunSolverGridCont(); CHKERRQ(ierr);
        } else {
            ss << "max(nx) = " << nxmax << " too small for grid continuation; switching to default solver";
            ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
            ierr = this->RunSolver(); CHKERRQ(ierr);
        }
    } else {
        ierr = this->RunSolver(); CHKERRQ(ierr);
    }

    ierr = this->DispLevelMsg("optimization done", rank); CHKERRQ(ierr);
    ierr = this->Finalize(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief run single level solver (no grid, scale, or parameter
 * continuation is performed)
 ********************************************************************/
PetscErrorCode RegistrationInterface::RunSolver() {
    PetscErrorCode ierr = 0;
    Vec mT = NULL, mR = NULL, x = NULL;
    ScalarType beta, betastar;
    bool boundreached, monitor;
    int level = 0;
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // do the setup
    ierr = this->SetupSolver(); CHKERRQ(ierr);

    ierr = Assert(this->m_Optimizer != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = SetupData(mR, mT); CHKERRQ(ierr);

    // initialize registration problem
    ierr = this->m_RegProblem->SetControlVariable(this->m_Solution); CHKERRQ(ierr);
    ierr = this->m_RegProblem->InitializeSolver(); CHKERRQ(ierr);
    ierr = this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

    monitor = this->m_Opt->GetRegMonitor().detdgradenabled;

    // init solver
    ierr = this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // reset all the clocks we have used so far
    ierr = this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr = this->m_Opt->ResetCounters(); CHKERRQ(ierr);
    this->m_Opt->SetSolutionStatus(0); // assume everything is good

    if (monitor) {
        boundreached = true; // enter search
        betastar = this->m_Opt->GetBeta(3);
        beta     = this->m_Opt->GetBeta(0);

        // if not set by user, assume something
        if (betastar == 0.0 || betastar == beta) betastar = beta/1E-1;
        ierr = Assert(betastar > beta, "parameter error"); CHKERRQ(ierr);

        // we will always enter here
        while (boundreached && level < 10) {
            // set initial guess
            ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

            // run the optimization
            ierr = this->m_Optimizer->Run(); CHKERRQ(ierr);

            // get the solution
            ierr = this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);

            // check bounds on determinant of deformation gradient
            ierr = this->m_RegProblem->CheckBounds(x, boundreached); CHKERRQ(ierr);

            if (boundreached) {
                beta += (betastar - beta) / 2.0;

                ss << "increasing beta from " << this->m_Opt->GetBeta(0) << " to " << beta;
                ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
                ss.clear(); ss.str(std::string());

                this->m_Opt->SetBeta(0, beta);
                ++ level;
            }
        }

        if (level == 0) {
            // nothing happened (solution and beta are valid)
            if (boundreached == false) {
                this->m_Opt->SetSolutionStatus(0);
                ierr = this->m_Solution->SetComponents(x); CHKERRQ(ierr);
            } else {
                ierr = ThrowError("should never happen"); CHKERRQ(ierr);
            }
        } else {
            if (boundreached == false) {
                // search for beta worked, so accept the solution
                this->m_Opt->SetSolutionStatus(1);
                ierr = this->m_Solution->SetComponents(x); CHKERRQ(ierr);
            } else {
                this->m_Opt->SetSolutionStatus(2);
                this->m_Opt->SetBeta(0, betastar);
            }
        }

    } else {
        // set initial guess
        ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run optimization
        ierr = this->m_Optimizer->Run(); CHKERRQ(ierr);

        // get solution
        ierr = this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
        ierr = this->m_Solution->SetComponents(x); CHKERRQ(ierr);
    }

    // finalize registration
    ierr = this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // destroy vectors if they were allocated
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr);}
    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr);}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief main function to run the parameter continuation
 ********************************************************************/
PetscErrorCode RegistrationInterface::RunSolverRegParaCont() {
    PetscErrorCode ierr = 0;
    Vec mT = NULL, mR = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // do the setup
    ierr = this->SetupSolver(); CHKERRQ(ierr);

    // check if setup was complete
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_Optimizer != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = SetupData(mR, mT); CHKERRQ(ierr);

    // reset all the clocks we have used so far
    ierr = this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr = this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // switch between the different strategies for
    // doing the parameter continuation (default one
    // is binary search)
    switch (this->m_Opt->GetParaCont().strategy) {
        case PCONTBINSEARCH:
        {
            ierr = this->RunSolverRegParaContBinarySearch(); CHKERRQ(ierr);
            break;
        }
        case PCONTREDUCESEARCH:
        {
            ierr = this->RunSolverRegParaContReductSearch(); CHKERRQ(ierr);
            break;
        }
        case PCONTINUATION:
        {
            ierr = this->RunSolverRegParaContReduction(); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("parameter continuation strategy not valid"); CHKERRQ(ierr);
            break;
        }
    }

    // destroy vector
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr);}
    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr);}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief run the optimizer; we search for an optimal
 * regularization weight using a binary search; we reduce/lift the
 * regularization parameter until we found a deformation map that
 * is diffeomorphic and results in a map that is close to the bound
 * on jacobian set by user
 *******************************************************************/
PetscErrorCode RegistrationInterface::RunSolverRegParaContBinarySearch() {
    PetscErrorCode ierr = 0;
    int maxsteps, level, rank;
    bool stop, boundreached, converged;
    std::ofstream logwriter;
    std::stringstream ss;
    std::string filename, msg;
    ScalarType beta, betamin, betascale, dbetascale,
                betastar, betahat, dbeta, dbetamin;
    Vec x = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_Optimizer != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // get parameters
    betamin = this->m_Opt->GetBetaMinParaCont();
    betascale = this->m_Opt->GetParaCont().betascale;
    maxsteps = this->m_Opt->GetParaCont().maxsteps;

    ierr = Assert(betascale < 1.0 && betascale > 0.0, "scale for beta not in (0,1)"); CHKERRQ(ierr);
    ierr = Assert(betamin > 0.0 && betamin < 1.0, "lower bound for beta in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr = this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess for current level
    ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    // initialize parameters (can be user defined)
    beta = this->m_Opt->GetParaCont().beta0;
    betastar = beta;

    // initialize registration problem (evaluate objective and gradient
    // for zero velocity field)
    this->m_Opt->SetBeta(0, beta);
//    this->m_Opt->SetBeta(1, beta);
    ierr = this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 0) {
        ierr = DbgMsg("starting coarse search for regularization weight"); CHKERRQ(ierr);
    }

    // reduce regularization parameter by one order of magnitude until
    // we hit tolerance
    stop = false; level = 0;
    while (level < maxsteps) {
        this->m_Opt->SetBeta(0, beta);
//        this->m_Opt->SetBeta(1, beta);
        //this->m_Opt->InitialGradNormSet(false);

        ss << std::scientific << std::setw(3)
            << "level " << level << " ( betav=" << beta
            << "; betav*=" << betastar << " )";
        ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        if (this->m_Opt->GetOptPara().fastsolve) {
//            ierr = this->m_RegProblem->InitializeOptimization(this->m_Solution); CHKERRQ(ierr);
            ierr = this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);
        }

        // set initial guess for current level
        ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run the optimization
        ierr = this->m_Optimizer->Run(); CHKERRQ(ierr);

        // if we did not converge, we do not want to decrease
        // the regularization parameter
        ierr = this->m_Optimizer->GetSolutionStatus(converged); CHKERRQ(ierr);
        if (!converged && level > 0) break;

        // get the solution
        ierr = this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);

        // check bounds on jacobian
        ierr = this->m_RegProblem->CheckBounds(x, stop); CHKERRQ(ierr);

        // we have to make sure that the initial parameter was
        // not too small
        if (stop) {
            if (level > 0) {
                break;  ///< if bound reached go home
            } else {
                // we reached bound in first step -> increase beta
                beta /= betascale;
                level = -1;  ///< reset level to 0
            }
        } else {
            // remember regularization parameter
            betastar = beta;
            // if we got here, the solution is valid
            ierr = this->m_Solution->SetComponents(x); CHKERRQ(ierr);
            // reduce beta
            beta *= betascale;
        }

        // if regularization parameter is smaller than
        // lower bound, let's stop this
        if (beta < betamin) {
            if (this->m_Opt->GetVerbosity() > 0) {
                ss << std::scientific
                   << "regularization parameter smaller than lower bound (betav="
                   << beta << " < " << betamin << "=betavmin)";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.str(std::string()); ss.clear();
            }
            break;
        }
        ++level;
    }  ///< until we hit the tolerance

    if (this->m_Opt->GetVerbosity() > 0) {
        ierr = DbgMsg("starting fine search for regularization weight"); CHKERRQ(ierr);
    }

    stop = false;  ///< reset search

    // compute new trial beta
    betahat = betascale*betastar;
    dbeta = (betastar-betahat)/2.0;
    beta = betastar-dbeta;
    ++level;

    // get scale for delta beta; this parameter determines how
    // accurate we solve (how many digits) with respect to the
    // order of magnitude of magnitude of the regularization
    // parameter)
    dbetascale = this->m_Opt->GetParaCont().dbetascale;
    msg = "scale for delta betav not in (0,1)";
    ierr = Assert(dbetascale < 1.0 && dbetascale > 0.0, msg); CHKERRQ(ierr);
    dbetamin = dbetascale*betastar;

    while (!stop) {
        // set regularization parameter
        this->m_Opt->SetBeta(0, beta);
//        this->m_Opt->SetBeta(1, beta);

        // display regularization parameter to user
        ss << std::setw(3) << "level " << level << " ( betav="
           << beta << "; betav*=" << betastar << " )";
        ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        if (this->m_Opt->GetOptPara().fastsolve) {
//            ierr = this->m_RegProblem->InitializeOptimization(this->m_Solution); CHKERRQ(ierr);
            ierr = this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);
        }

        // set initial guess for current level
        ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run the optimization
        ierr = this->m_Optimizer->Run(); CHKERRQ(ierr);

        // if we did not converge, beta is too small, also
        ierr = this->m_Optimizer->GetSolutionStatus(converged); CHKERRQ(ierr);
        if (converged) {
            // get the solution
            ierr = this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);

            // check bounds on jacobian
            boundreached = false;
            ierr = this->m_RegProblem->CheckBounds(x, boundreached); CHKERRQ(ierr);

            // if bound is reached, the lower bound is now beta
            // if not, beta is our new best estimate
            if (boundreached) {
                betahat = beta;
            } else {
                betastar = beta;  // new best estimate

                // if we got here, the solution is valid
                ierr = this->m_Solution->SetComponents(x); CHKERRQ(ierr);
            }
        } else {
            ierr = WrngMsg("solver did not converge"); CHKERRQ(ierr);
            betahat = beta;
        }

        // increase or reduce beta
        dbeta = (betastar - betahat)/2.0;
        beta = betastar - dbeta;
        if (fabs(dbeta) < dbetamin) {
            stop = true;
            if (this->m_Opt->GetVerbosity() > 0) {
                ss  << std::setw(3) << "update for beta too small ( dbeta="
                    << fabs(dbeta) << " < " << dbetamin << "=dbetamin )";
                ierr = DbgMsg(ss.str());
                ss.str(std::string()); ss.clear();
            }
        }
        ++level;
    }

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(), '-') << std::endl;
    ss << std::scientific << "estimated regularization parameter betav=" << betastar;
    ierr = Msg(ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(), '-') << std::endl;
    ss.str(std::string()); ss.clear();

    // if output folder is set
    if (!this->m_Opt->GetReadWriteFlags().xfolder.empty()) {
        if (rank == 0) {
            filename = this->m_Opt->GetReadWriteFlags().xfolder;
            filename += "parameter-continuation-estimated-beta.log";
            // create output file or append to output file
            logwriter.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
            ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);
            ss << std::scientific << "betav " << std::setw(3) << std::right << betastar;
            logwriter << ss.str() << std::endl;
            ss.str(std::string()); ss.clear();
        }
    }

    // wrap up
    ierr = this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief we solves the optimization problem by simply reducing
 * the regularization parameter until the mapping becomes
 * non-diffeomorphic/breaches the user defined bound; stored
 * velocity field (solution) is last iterate that resulted in
 * diffeomorphic deformation map (as judged by the determinant
 * of the deformation gradient)
 *******************************************************************/
PetscErrorCode RegistrationInterface::RunSolverRegParaContReductSearch() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    ScalarType beta, betamin, betastar, betascale;
    Vec x = NULL;
    int maxsteps, rank, level;
    bool stop;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // get parameters
    betamin = this->m_Opt->GetBetaMinParaCont();
    maxsteps = this->m_Opt->GetParaCont().maxsteps;
    betascale = this->m_Opt->GetParaCont().betascale;

    ierr = Assert(betascale < 1.0 && betascale > 0.0, "scale for beta not in (0,1)"); CHKERRQ(ierr);
    ierr = Assert(betamin > 0.0 && betamin < 1.0, "lower bound for beta in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr = this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess for current level
    ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    // set initial regularization weight
    beta = this->m_Opt->GetParaCont().beta0;
    betastar = beta;

    // initialize registration problem (evaluate objective and gradient
    // for zero velocity field)
    this->m_Opt->SetBeta(0, beta);
//    this->m_Opt->SetBeta(1, beta);
    ierr = this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

    // reduce regularization parameter by one order of magnitude until
    // we hit user defined tolerances (which either is a lower bound
    // on the regularization parameter or a lower bound on the
    // determinant of the deformation gradient)
    level = 0;
    while (level < maxsteps) {
        // set regularization weight
        this->m_Opt->SetBeta(0, beta);
//        this->m_Opt->SetBeta(1, beta);

        // display message to user
        ss << std::scientific << std::setw(3) << "level " << level << " (beta=" << beta << ")";
        ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        // set initial guess for current level
        ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run the optimization
        ierr = this->m_Optimizer->Run(); CHKERRQ(ierr);

        // get the solution
        ierr = this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);

        // check bounds on jacobian
        stop = false;
        ierr = this->m_RegProblem->CheckBounds(x, stop); CHKERRQ(ierr);

        if (stop) break; // if bound reached go home

        // remember best estimate
        betastar = beta;

        // if we got here, the solution is valid
        ierr = this->m_Solution->SetComponents(x); CHKERRQ(ierr);

        beta *= betascale; // reduce beta

        // if the regularization parameter is smaller than
        // the lower bound, we're done
        if (beta < betamin) {
            if (this->m_Opt->GetVerbosity() > 0) {
                ss << "regularization parameter smaller than lower bound (betav="
                   << beta << " < " << betamin << "=betavmin)";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.str(std::string()); ss.clear();
            }
            break;
        }
        ++level;
    }  ///< parameter reduction

    ss << std::scientific << "estimated regularization parameter betav=" << betastar;
    ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();

    // wrap up
    ierr = this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief we solves the optimization problem by simply reducing
 * the regularization parameter until we have reached the
 * target regularization weight set by the user
 *******************************************************************/
PetscErrorCode RegistrationInterface::RunSolverRegParaContReduction() {
    PetscErrorCode ierr = 0;
    int level, rank;
    std::stringstream ss;
    ScalarType beta, betastar;
    Vec x;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // get target regularization weight
    betastar = this->m_Opt->GetParaCont().targetbeta;
    ierr = Assert(betastar > 0.0 && betastar < 1.0, "target beta not in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr = this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess for current level
    ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    // reduce regularization parameter
    level = 0; beta = 1.0;

    // initialize registration problem (evaluate objective and gradient
    // for zero velocity field)
    this->m_Opt->SetBeta(0, beta);
//    this->m_Opt->SetBeta(1, beta);
    //ierr = this->m_RegProblem->InitializeOptimization(this->m_Solution); CHKERRQ(ierr);
    ierr = this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

    while (beta > betastar) {
        // set regularization weight
        this->m_Opt->SetBeta(0, beta);
//        this->m_Opt->SetBeta(1, beta);

        // display message to user
        ss << std::scientific << std::setw(3) << "level "
           << level << " (beta=" << beta << "; beta*=" << betastar << ")";
        ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        // run the optimization (TODO: make the fast solve modular)
        ierr = this->m_Optimizer->Run(true); CHKERRQ(ierr);

        // reduce by one order of magnitude
        beta /= static_cast<ScalarType>(10); // reduce beta
        ++level;
    } // parameter reduction

    beta = betastar;

    // set regularization weight
    this->m_Opt->SetBeta(0, beta);
//    this->m_Opt->SetBeta(1, beta);

    // display message to user
    ss << std::scientific << std::setw(3)
        << "level " << level << " (beta="
        << beta << "; beta*=" << betastar << ")";
    ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();

    // solve optimization problem for user defined regularization parameter
    ierr = this->m_Optimizer->Run(); CHKERRQ(ierr);

    // get the solution
    ierr = this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
    ierr = this->m_Solution->SetComponents(x); CHKERRQ(ierr);

    // wrap up
    ierr = this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief run solver using a scale continuation scheme; that is, we
 * will successively reduce the smoothing of the images to be
 * registered to get to finer and finer scales; this is supposed to
 * reduce the non-linearity in the problem
 *******************************************************************/
PetscErrorCode RegistrationInterface::RunSolverScaleCont() {
    PetscErrorCode ierr = 0;
    Vec mT = NULL, mR = NULL, x = NULL;
    std::stringstream ss;
    int level, maxlevel = 6, rank;
    bool solve;
    ScalarType sigma[3], nxhalf;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // set up preprocessing
    if (this->m_PreProc == NULL) {
        try {this->m_PreProc = new Preprocessing(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // do the setup
    ierr = this->SetupSolver(); CHKERRQ(ierr);

    // check if everything has been set up correctly
    ierr = Assert(this->m_Optimizer != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // set up synthetic problem if we did not read images
    if (!this->m_IsTemplateSet && !this->m_IsReferenceSet) {
        ierr = this->m_RegProblem->SetupSyntheticProb(this->m_ReferenceImage, this->m_TemplateImage); CHKERRQ(ierr);
    }

    // check if images have been set
    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate local images
    ierr = VecDuplicate(this->m_TemplateImage, &mT); CHKERRQ(ierr);
    ierr = VecDuplicate(this->m_ReferenceImage, &mR); CHKERRQ(ierr);

    // set images
    ierr = this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
    ierr = this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

    // reset all the clocks we have used so far
    ierr = this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr = this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // set optimization problem
    ierr = this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess for current level
    ierr = this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    level=0;
    while (level < maxlevel) {
        solve = true;
        for (int i = 0; i < 3; ++i) {
            // get and set sigma for current level
            sigma[i] = this->m_Opt->GetScaleContPara().sigma[i][level];
            this->m_Opt->SetSigma(i,sigma[i]);

            // if sigma is bigger than half of the grid size, don't compute
            nxhalf = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[i])/2.0;
            if (nxhalf <= sigma[i]) solve = false;
        }

        // solve problem
        if (solve) {
            ierr = this->m_PreProc->Smooth(mT, this->m_TemplateImage); CHKERRQ(ierr);
            ierr = this->m_PreProc->Smooth(mR, this->m_ReferenceImage); CHKERRQ(ierr);

            // rescale images (TODO: has to be removed, cause gaussian should be
            // scale invariant)
            //ierr = Rescale(mR, 0.0, 1.0); CHKERRQ(ierr);
            //ierr = Rescale(mT, 0.0, 1.0); CHKERRQ(ierr);

            // display message to user
            ss << std::scientific << std::setw(3)
               << "level " << level << " sigma=("
               << sigma[0] << "," << sigma[1] << "," << sigma[2] << ")";
            ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();

            // compute gradient, distance measure, and initial objective
            // value for zero velocity field, but updated images
            //ierr = this->m_RegProblem->InitializeOptimization(this->m_Solution); CHKERRQ(ierr);
            ierr = this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

            // run the optimization
            ierr = this->m_Optimizer->Run(); CHKERRQ(ierr);
        } else {
            ss << std::scientific << std::setw(3)
                << "skipping level " << level << " sigma=("
                << sigma[0] << "," << sigma[1] << "," << sigma[2] << ")";
            ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        ++level;  // increment counter
    }

    // get the solution
    ierr = this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
    ierr = this->m_Solution->SetComponents(x); CHKERRQ(ierr);

    // wrap up
    ierr = this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // destroy vector
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr);}
    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr);}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief run solver using a grid continuation scheme; that is, we
 * will successively increase the number of grid points of the
 * template and reference image
 ********************************************************************/
PetscErrorCode RegistrationInterface::RunSolverGridCont() {
    PetscErrorCode ierr = 0;
    int rank, level, nlevels, computelevel;
    std::stringstream ss;
    std::string ext;
    IntType nl, ng, isize[3], nx[3];
    Vec mT = NULL, mR = NULL, xstar = NULL;
    VecField *v = NULL;
    ScalarType greltol, tolscale = 10;
    bool solve;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // set up preprocessing
    if (this->m_PreProc == NULL) {
        try{this->m_PreProc = new Preprocessing(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    this->m_PreProc->ResetGridChangeOps(true);

    if (!this->m_IsTemplateSet && !this->m_IsReferenceSet) {
        // do the setup
        ierr = this->SetupSolver(); CHKERRQ(ierr);

        // check if everything has been set up correctly
        ierr = Assert(this->m_Optimizer != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

        // set up synthetic test problem
        ierr = this->m_RegProblem->SetupSyntheticProb(this->m_ReferenceImage,
                                                      this->m_TemplateImage); CHKERRQ(ierr);
    }

    // make sure images have not been set
    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate multilevel pyramid for reference image
    if (this->m_Opt->GetVerbosity() > 1) {
        ierr = DbgMsg("setup: reference image multilevel pyramid"); CHKERRQ(ierr);
    }
    if (this->m_ReferencePyramid == NULL) {
        try {this->m_ReferencePyramid = new MultiLevelPyramid(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    // setup multilevel pyramid for reference image
    ierr = this->m_ReferencePyramid->SetPreProc(this->m_PreProc); CHKERRQ(ierr);
    ierr = this->m_ReferencePyramid->DoSetup(this->m_ReferenceImage); CHKERRQ(ierr);

    // allocate multilevel pyramid for template image
    if (this->m_Opt->GetVerbosity() > 1) {
        ierr = DbgMsg("setup: template image multilevel pyramid"); CHKERRQ(ierr);
    }
    if (this->m_TemplatePyramid == NULL) {
        try {this->m_TemplatePyramid = new MultiLevelPyramid(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    // setup multilevel pyramid for template image
    ierr = this->m_TemplatePyramid->SetPreProc(this->m_PreProc); CHKERRQ(ierr);
    ierr = this->m_TemplatePyramid->DoSetup(this->m_TemplateImage); CHKERRQ(ierr);

    // get grid size
    for (int i = 0; i < 3; ++i) {
        nx[i] = this->m_Opt->GetGridContPara().nx[0][i];
    }
    ierr = this->m_Opt->GetSizes(nx, nl, ng); CHKERRQ(ierr);

    // TODO: allow for warm start
    try {v = new VecField(nl, ng);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    ierr = v->SetValue(0.0); CHKERRQ(ierr);

    // reset tolerance for gradient (optimization); we do not want
    // to solve as accurately when we solve on the coarse grid
    greltol = this->m_Opt->GetOptPara().tol[2];

    if (greltol < 0.01) {
        if (this->m_Opt->GetVerbosity() > 1) {
            ss  << std::scientific << "increasing tolerance for gradient: "
                << greltol << " >> " << tolscale*greltol;
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        this->m_Opt->SetOptTol(2, tolscale*greltol);
    }

    // reset all the clocks we have used so far
    ierr = this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr = this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // get number of levels
    nlevels = this->m_Opt->GetGridContPara().nlevels;

    // run multi-level solver
    computelevel = 0;
    level = 0; // this->m_Opt->GetGridContPara().minlevel;
    while (level < nlevels) {
        // get number of grid points for current level
        for (int i = 0; i < 3; ++i) {
            nx[i] = this->m_Opt->GetGridContPara().nx[level][i];
            isize[i] = this->m_Opt->GetGridContPara().isize[level][i];
        }
        nl = this->m_Opt->GetGridContPara().nl[level];
        ng = this->m_Opt->GetGridContPara().ng[level];

        // display user message
        ss << std::scientific << "level " << std::setw(3) << level
           << "    nx=(" << nx[0] << "," << nx[1] << "," << nx[2]
           << "); (nl,ng)=(" << nl << "," << ng
           << "); isize=(" << isize[0] << ","
           << isize[1] << "," << isize[2] << ")";
         ierr = this->DispLevelMsg(ss.str(), rank); CHKERRQ(ierr);
         ss.str(std::string()); ss.clear();

        solve = true;
        if (this->m_Opt->GetPDESolverPara().type == SL) {
            if (isize[0] < 6 || isize[1] < 6) solve = false;
        }

        if (solve) {
            // get the individual images from the pyramid
            ierr = this->m_ReferencePyramid->GetLevel(&mR, level); CHKERRQ(ierr);
            ierr = Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);
            ierr = this->m_TemplatePyramid->GetLevel(&mT, level); CHKERRQ(ierr);
            ierr = Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);

            // initialize
            for (int i = 0; i < 3; ++i) {
                this->m_Opt->SetNumGridPoints(i, nx[i]);
            }
            ierr = this->m_Opt->DoSetup(false); CHKERRQ(ierr);

            // store intermediate results
            if (this->m_Opt->GetReadWriteFlags().iterates) {
                ext = this->m_Opt->GetReadWriteFlags().extension;
                ss << "reference-image-level=" << level << ext;
                ierr = this->m_ReadWrite->Write(mR, ss.str()); CHKERRQ(ierr);
                ss.str(std::string()); ss.clear();

                ss << "template-image-level=" << level << ext;
                ierr = this->m_ReadWrite->Write(mT, ss.str()); CHKERRQ(ierr);
                ss.str(std::string()); ss.clear();
            }

            if (!this->m_Opt->GetGridContPara().maxit.empty()) {
                int l = this->m_Opt->GetGridContPara().maxit.size() - 1;
                l = computelevel > l ? l : computelevel;
                int maxit = this->m_Opt->GetGridContPara().maxit[l];
                if (this->m_Opt->GetVerbosity() > 1) {
                    ss  <<"setting max number of iterations to " << maxit;
                    ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                    ss.str(std::string()); ss.clear();
                }
                this->m_Opt->SetOptMaxIter(maxit);
            }

            // do the setup
            ierr = this->SetupSolver(); CHKERRQ(ierr);

            // set images
            ierr = this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);
            ierr = this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);

            // compute initial gradient, objective and
            // distance mesure for zero velocity field
            ierr = this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

            // set initial guess and registraiton problem
            ierr = this->m_Optimizer->SetInitialGuess(v); CHKERRQ(ierr);
            ierr = this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

            // reset tolerances
            if ((level == (nlevels-1)) && (greltol < 0.01)) {
                if (this->m_Opt->GetVerbosity() > 1) {
                    ss  << std::scientific << "reseting tolerance for gradient: "
                        << tolscale*greltol << " >> " << greltol;
                    ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                    ss.str(std::string()); ss.clear();
                }
                this->m_Opt->SetOptTol(2, greltol);
            }

            // run the optimizer
            ierr = this->m_Optimizer->Run(); CHKERRQ(ierr);

            // get and parse solution
            ierr = this->m_Optimizer->GetSolution(xstar); CHKERRQ(ierr);
            ierr = v->SetComponents(xstar); CHKERRQ(ierr);

            ++computelevel;
        } else {
            ss << "skipping level " << level;
            ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        ++level;  // increment level

        if (level < nlevels) {
            ierr = this->ProlongVelocityField(v, level); CHKERRQ(ierr);
            if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr); mR = NULL;}
            if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr); mT = NULL;}
        }
    }

    // wrap up
    ierr = this->m_RegProblem->Finalize(v); CHKERRQ(ierr);

    if (v != NULL) {delete v; v = NULL;};
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr); mR = NULL;}
    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr); mT = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief prolong velocity field
 ********************************************************************/
PetscErrorCode RegistrationInterface::ProlongVelocityField(VecField*& v, int level) {
    PetscErrorCode ierr = 0;
    IntType nx_f[3], nx_c[3], nl, ng;
    VecField *v_f = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("prolonging velocity field"); CHKERRQ(ierr);
    }

    // set up preprocessing
    if (this->m_PreProc == NULL) {
        try{this->m_PreProc = new Preprocessing(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    this->m_PreProc->ResetGridChangeOps(true);

    // get number of grid points for current level
    for (int i = 0; i < 3; ++i) {
        nx_f[i] = this->m_Opt->GetGridContPara().nx[level  ][i];
        nx_c[i] = this->m_Opt->GetGridContPara().nx[level-1][i];
    }

    // get number of points to allocate
    ierr = this->m_Opt->GetSizes(nx_f, nl, ng); CHKERRQ(ierr);

    // allocate container for velocity field
    try {v_f = new reg::VecField(nl, ng);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // apply prolongation operator
    ierr = this->m_PreProc->Prolong(v_f, v, nx_f, nx_c); CHKERRQ(ierr);

    // allocate container for velocity field
    if (v != NULL) {delete v; v = NULL;}
    try {v = new reg::VecField(nl, ng);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    ierr = v->Copy(v_f); CHKERRQ(ierr);

    this->m_PreProc->ResetGridChangeOps(false);

    if (v_f != NULL) {delete v_f; v_f = NULL;}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief finalize optimization (displays information for user)
 ********************************************************************/
PetscErrorCode RegistrationInterface::Finalize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // finalize optimizer (show tao output)
    ierr = this->m_Optimizer->Finalize(); CHKERRQ(ierr);

    // display time to solution
    ierr = this->m_Opt->DisplayTimeToSolution(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief run postprocessing of input data
 ********************************************************************/
PetscErrorCode RegistrationInterface::RunPostProcessing() {
    PetscErrorCode ierr = 0;
    Vec mR = NULL, mT = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = this->SetupRegProblem(); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // user needs to set template and reference image and the solution
    ierr = Assert(this->m_Solution != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate image containers
    ierr = VecDuplicate(this->m_TemplateImage, &mT); CHKERRQ(ierr);
    ierr = VecDuplicate(this->m_ReferenceImage, &mR); CHKERRQ(ierr);

    if (this->m_Opt->GetRegFlags().applysmoothing) {
        // allocate preprocessing class
        if (this->m_PreProc == NULL) {
            try{this->m_PreProc = new Preprocessing(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        // apply smoothing
        ierr = this->m_PreProc->Smooth(mR, this->m_ReferenceImage); CHKERRQ(ierr);
        ierr = this->m_PreProc->Smooth(mT, this->m_TemplateImage); CHKERRQ(ierr);
    } else {
        // copy input images
        ierr = VecCopy(this->m_ReferenceImage, mR); CHKERRQ(ierr);
        ierr = VecCopy(this->m_TemplateImage, mT); CHKERRQ(ierr);
    }

    // set reference and template images
    ierr = this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
    ierr = this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

    // compute stuff
    ierr = this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // destroy vectors
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr);}
    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr);}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem, given some trial velocity
 * field v and a template image m0
 * @param[in] m0 initial condition/template image
 * @param[out] m1 transported template image
 ********************************************************************/
PetscErrorCode RegistrationInterface::SolveForwardProblem(Vec m1, Vec m0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m0 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(m1 != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_RegProblem == NULL) {
        ierr = this->SetupRegProblem(); CHKERRQ(ierr);
    }
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_Solution != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Opt->GetRegFlags().applysmoothing) {
        // allocate preprocessing class
        if (this->m_PreProc == NULL) {
            try {this->m_PreProc = new Preprocessing(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        // apply smoothing
        ierr = this->m_PreProc->Smooth(m1, m0); CHKERRQ(ierr);
        ierr = VecCopy(m1, m0); CHKERRQ(ierr);
        ierr = VecSet(m1, 0.0); CHKERRQ(ierr);
    }

    // set the control variable
    ierr = this->m_RegProblem->SetControlVariable(this->m_Solution); CHKERRQ(ierr);
    ierr = this->m_RegProblem->SolveForwardProblem(m1, m0); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem, given some trial velocity
 * field v and a deformed template image m1
 * @param[in] m1 deformed template image
 * @param[out] l0 adjoint variable at t=0
 ********************************************************************/
PetscErrorCode RegistrationInterface::SolveAdjointProblem(Vec l0, Vec m1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = this->SetupRegProblem(); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // user needs to set template and reference image and the solution
    ierr = Assert(this->m_Solution != NULL, "null pointer"); CHKERRQ(ierr);

    // set the control variable
    ierr = this->m_RegProblem->SetControlVariable(this->m_Solution); CHKERRQ(ierr);
    ierr = this->m_RegProblem->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
    ierr = this->m_RegProblem->SolveAdjointProblem(l0, m1); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map or deformation gradient
 ********************************************************************/
PetscErrorCode RegistrationInterface::ComputeDefFields() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = this->SetupRegProblem(); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // user needs to set template and reference image and the solution
    ierr = Assert(this->m_Solution != NULL, "null pointer"); CHKERRQ(ierr);

    // compute stuff
    ierr = this->m_RegProblem->SetControlVariable(this->m_Solution); CHKERRQ(ierr);

    if (this->m_Opt->GetReadWriteFlags().detdefgrad) {
        ierr = Msg("computing determinant of deformation gradient"); CHKERRQ(ierr);
        ierr = this->m_RegProblem->ComputeDetDefGrad(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->GetReadWriteFlags().defgrad) {
        ierr = Msg("computing deformation gradient"); CHKERRQ(ierr);
        ierr = this->m_RegProblem->ComputeDefGrad(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->GetReadWriteFlags().defmap) {
        ierr = Msg("computing deformation map"); CHKERRQ(ierr);
        ierr = this->m_RegProblem->ComputeDeformationMap(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->GetReadWriteFlags().deffield) {
        ierr = Msg("computing displacement field"); CHKERRQ(ierr);
        ierr = this->m_RegProblem->ComputeDisplacementField(true); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 ********************************************************************/
PetscErrorCode RegistrationInterface::ComputeDetDefGrad(Vec detj) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = this->SetupRegProblem(); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // user needs to set template and reference image and the solution
    ierr = Assert(this->m_Solution != NULL, "null pointer"); CHKERRQ(ierr);

    // compute stuff
    ierr = this->m_RegProblem->SetControlVariable(this->m_Solution); CHKERRQ(ierr);

    ierr = Msg("computing determinant of deformation gradient"); CHKERRQ(ierr);
    ierr = this->m_RegProblem->ComputeDetDefGrad(false, detj); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map
 ********************************************************************/
PetscErrorCode RegistrationInterface::ComputeDeformationMap(VecField* y) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = this->SetupRegProblem(); CHKERRQ(ierr);
    ierr = Assert(this->m_RegProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // user needs to set template and reference image and the solution
    ierr = Assert(this->m_Solution != NULL, "null pointer"); CHKERRQ(ierr);

    // compute stuff
    ierr = this->m_RegProblem->SetControlVariable(this->m_Solution); CHKERRQ(ierr);

    ierr = Msg("computing deformation map"); CHKERRQ(ierr);
    ierr = this->m_RegProblem->ComputeDeformationMap(false, y); CHKERRQ(ierr);

    if (this->m_Opt->GetRegFlags().checkdefmapsolve) {
        ierr = this->m_RegProblem->CheckDefMapConsistency(); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  //  namespace reg




#endif  // _REGISTRATIONINTERFACE_CPP_
