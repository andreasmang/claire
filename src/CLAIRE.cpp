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

#ifndef _CLAIRE_CPP_
#define _CLAIRE_CPP_

// global includes
#include <string>
#include <algorithm>

// local includes
#include "CLAIRE.hpp"

#ifdef REG_HAS_CUDA
#include "cuda_helper.hpp"
#endif

#include "PreconditionerKernel.hpp"
#include "TwoLevel.hpp"

#ifdef ZEITGEIST_DEV
std::string ZeitGeist::prefix = "";
#endif

namespace reg {

/********************************************************************
 * @brief default constructor
 *******************************************************************/
CLAIRE::CLAIRE() : SuperClass() {
    this->Initialize();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
CLAIRE::CLAIRE(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}

/********************************************************************
 * @brief default destructor
 *******************************************************************/
CLAIRE::~CLAIRE() {
    this->ClearMemory();
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode CLAIRE::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_StateVariable = NULL;       ///< state variable
    this->m_AdjointVariable = NULL;     ///< adjoint variable
    this->m_IncStateVariable = NULL;    ///< incremental state variable
    this->m_IncAdjointVariable = NULL;  ///< incremental adjoint variable
    
    this->m_CoarseReg = nullptr;
    this->m_CoarseRegOpt = nullptr;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode CLAIRE::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
#ifdef ZEITGEIST
    if (this->m_Opt->m_Domain.level == 0) {
      Msg("-----------------------------------------------------------------------------------------------------");
      Msg("ZeitGeist:");
      for (auto zg : ZeitGeist::zgMap()) {
        char txt[120];
        double global_runtime;
        double local_runtime = zg.second.Total_s();
        MPI_Reduce(&local_runtime, &global_runtime, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
        global_runtime /= this->m_Opt->rank_cnt;
        sprintf(txt, "  %-60s: %5lix, %0.10lf",zg.first.c_str(), zg.second.Count(), global_runtime);
        Msg(txt);
      }
      Msg("-----------------------------------------------------------------------------------------------------");
    }
#endif
#ifdef REG_HAS_CUDA
  if (this->m_Opt->m_Verbosity > 2) {
    cudaPrintDeviceMemory();
  }
#endif
    
    ierr = Free(this->m_CoarseReg); CHKERRQ(ierr);
    if (this->m_CoarseRegOpt) {
      this->m_CoarseRegOpt->m_FFT.fft = nullptr;
    }
    ierr = Free(this->m_CoarseRegOpt); CHKERRQ(ierr);

    // delete all variables
    ierr = this->ClearVariables(); CHKERRQ(ierr);
  
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode CLAIRE::ClearVariables(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // delete all variables
    if (this->m_Opt->m_Verbosity > 2) {
      std::stringstream ss;
      size_t total = 0;
      if (this->m_StateVariable) total += this->m_StateVariable->GetSize();
      if (this->m_AdjointVariable) total += this->m_AdjointVariable->GetSize();
      if (this->m_IncStateVariable) total += this->m_IncStateVariable->GetSize();
      if (this->m_IncAdjointVariable) total += this->m_IncAdjointVariable->GetSize();
      ss << "memory allocated: "<< std::scientific << total;
      ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
      ss.clear(); ss.str(std::string());
    }
    
    Free(this->m_StateVariable);
    Free(this->m_AdjointVariable);
    Free(this->m_IncStateVariable);
    Free(this->m_IncAdjointVariable);

    PetscFunctionReturn(ierr);
}

PetscErrorCode CLAIRE::CreateCoarseReg() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = Assert(this->m_CoarseRegOpt != nullptr, "coarse grid RegOpt not initialized"); CHKERRQ(ierr);
  
  ierr = AllocateOnce<CLAIRE>(this->m_CoarseReg, this->m_CoarseRegOpt); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode CLAIRE::InitializeCoarseReg() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = Assert(this->m_Opt != nullptr, "fine grid RegOpt not initialized"); CHKERRQ(ierr);
  
  //printf("%p\n", this->m_Opt->m_FFT_coarse.fft);
  //ierr = Assert(this->m_Opt->m_FFT.fft != nullptr, "no FFT"); CHKERRQ(ierr);
  if (this->m_Opt->m_FFT_coarse.fft) {
    
    ierr = AllocateOnce(this->m_CoarseRegOpt, *this->m_Opt); CHKERRQ(ierr);
  
    this->m_CoarseRegOpt->m_FFT = this->m_Opt->m_FFT_coarse;
    this->m_CoarseRegOpt->m_FFT_coarse.fft = nullptr;
    this->m_CoarseRegOpt->m_Domain.level = this->m_Opt->m_Domain.level;
    this->m_CoarseRegOpt->m_Domain.nl = 1;
    this->m_CoarseRegOpt->m_Domain.ng = 1;
    this->m_CoarseRegOpt->m_Domain.mpicomm = this->m_Opt->m_Domain.mpicomm;
    this->m_CoarseRegOpt->m_Domain.rowcomm = this->m_Opt->m_Domain.rowcomm;
    this->m_CoarseRegOpt->m_Domain.colcomm = this->m_Opt->m_Domain.colcomm;
    for (int i=0; i<3; ++i) {
      this->m_CoarseRegOpt->m_Domain.nx[i] = this->m_Opt->m_FFT_coarse.nx[i];
      this->m_CoarseRegOpt->m_Domain.istart[i] = this->m_Opt->m_FFT_coarse.istart[i];
      this->m_CoarseRegOpt->m_Domain.isize[i] = this->m_Opt->m_FFT_coarse.isize[i];
      this->m_CoarseRegOpt->m_Domain.ng *= this->m_CoarseRegOpt->m_Domain.nx[i];
      this->m_CoarseRegOpt->m_Domain.nl *= this->m_CoarseRegOpt->m_Domain.isize[i];
      this->m_CoarseRegOpt->m_Domain.hx[i] = 2.*PETSC_PI/static_cast<ScalarType>(this->m_CoarseRegOpt->m_Domain.nx[i]);
    }
    this->m_CoarseRegOpt->m_SetupDone = true;
    
    ierr = this->CreateCoarseReg(); CHKERRQ(ierr);
    
    ierr = this->m_CoarseReg->SetControlVariable(nullptr); CHKERRQ(ierr);
    
    ierr = this->m_CoarseReg->InitializeSolver();
  }
  
  PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief setup solver (to not have to allocate things on the
 * fly; this allows us to essentially do a warm start)
 *******************************************************************/
PetscErrorCode CLAIRE::InitializeSolver(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = AllocateOnce(this->m_StateVariable, this->m_Opt, true, true); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_AdjointVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_IncAdjointVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_IncStateVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField3, this->m_Opt); CHKERRQ(ierr);
    
    
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);

    /*if (this->m_Opt->m_PDESolver.type == SL) {
        ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
        ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }*/

    //ierr = AllocateOnce(this->m_DifferentiationFD, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce<DifferentiationSM>(this->m_Differentiation, this->m_Opt); CHKERRQ(ierr);
    if (this->m_DeformationFields != NULL) {
      this->m_DeformationFields->SetDifferentiation(this->m_Differentiation);
    }

    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }
    
    if (this->m_TransportProblem == nullptr) {
        ierr = this->SetupTransportProblem(); CHKERRQ(ierr);
    }
    
    if (this->m_Opt->m_PDESolver.type == SL) {
      ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
      ierr = this->m_TransportProblem->InitializeControlVariable(this->m_VelocityField); CHKERRQ(ierr);
    }
    
    //ierr = this->InitializeCoarseReg(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief initialize the optimization (we essentially evaluate
 * the objective functional and the gradient for a given initial
 * guess)
 *******************************************************************/
PetscErrorCode CLAIRE::InitializeOptimization() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    std::stringstream ss;
    ScalarType value, alpha, jvt, jv, lsred, descent;
    Vec g = NULL, dv = NULL, v = NULL, vtilde = NULL;
    bool lssuccess, restoreinitialguess = false;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // if velocity field is null pointer, we did not set
    // any initial guess
    if (this->m_VelocityField == NULL) {
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = DbgMsg2("allocating velocity field"); CHKERRQ(ierr);
        }
        ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
    } else { // user might have provided initial guess
        ierr = this->IsVelocityZero(); CHKERRQ(ierr);
        if (!this->m_VelocityIsZero) {
            restoreinitialguess = true;
            ierr = AllocateOnce(this->m_WorkVecField3, this->m_Opt); CHKERRQ(ierr);
            ierr = this->m_WorkVecField3->Copy(this->m_VelocityField); CHKERRQ(ierr);
        }
    }
    ierr = VecCreate(g, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecCreate(v, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecSet(v, 0.0); CHKERRQ(ierr);

    // if we use a non-zero initial guess, we compute
    // the first velocity using a steepest descent approach
    if (!this->m_Opt->m_OptPara.usezeroinitialguess) {
        ierr = VecCreate(dv, 3*nl, 3*ng); CHKERRQ(ierr);
        ierr = VecCreate(vtilde, 3*nl, 3*ng); CHKERRQ(ierr);

        lsred = 1E-4;  // reduction rate for line search
        for (int l = 0; l < 1; ++l) {
            // evaluate objective function
            ierr = this->EvaluateObjective(&jv, v); CHKERRQ(ierr);

            // compute gradient
            ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);

            // compute search direction (gradient in sobolev space)
            ierr = this->EvaluateGradient(dv, v); CHKERRQ(ierr);

            // inner product between gradient and search direction
            ierr = VecTDot(g, dv, &descent); CHKERRQ(ierr);

            alpha = 1.0; lssuccess = false;
            for (int i = 0; i < 20; ++i) {
                // compute trial velocity
                ierr = VecWAXPY(vtilde, -alpha, dv, v); CHKERRQ(ierr);

                // evaluate objective function
                ierr = this->EvaluateObjective(&jvt, vtilde); CHKERRQ(ierr);

                // armijo rule
                if (jvt < jv + alpha*lsred*descent) {
                    lssuccess = true;
                    break;
                }
                alpha /= 2.0;
            }
            if (lssuccess) {
                if (this->m_Opt->m_Verbosity > 1) {
                    ss << "line search successful (initialization; alpha=" << std::scientific << alpha << ")";
                    ierr = DbgMsg1(ss.str()); CHKERRQ(ierr);
                    ss.clear(); ss.str(std::string());
                }
                ierr = VecCopy(vtilde, v); CHKERRQ(ierr);
            } else {
                if (this->m_Opt->m_Verbosity > 1) {
                    ss << "line search failed (initialization; alpha=" << std::scientific << alpha << ")";
                    ierr = DbgMsg1(ss.str()); CHKERRQ(ierr);
                    ss.clear(); ss.str(std::string());
                }
                break;
            }
        }
    }
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // evaluate distance measure
    ierr = this->EvaluateDistanceMeasure(&value); CHKERRQ(ierr);
    this->m_Opt->m_Monitor.dval0 = value;

    // evaluate objective functional
    ierr = this->EvaluateObjective(&value, v); CHKERRQ(ierr);
    this->m_Opt->m_Monitor.jval0 = value;

    // compute gradient
    ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);

    // compute gradient norm
    ierr = VecNorm(g, NORM_2, &value); CHKERRQ(ierr);
    this->m_Opt->m_Monitor.gradnorm0 = value;

    if (this->m_Opt->m_Verbosity > 0) {
        ss << "initial gradient norm: "<< std::scientific << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    // if we had a non-zero initial velocity, we'll restore it
    if (restoreinitialguess) {
        ierr = this->m_VelocityField->Copy(this->m_WorkVecField3); CHKERRQ(ierr);
    }

    // clean up
    if (v != NULL) {ierr = VecDestroy(&v); CHKERRQ(ierr); v = NULL;}
    if (g != NULL) {ierr = VecDestroy(&g); CHKERRQ(ierr); g = NULL;}
    if (dv != NULL) {ierr = VecDestroy(&dv); CHKERRQ(ierr); dv = NULL;}
    if (vtilde != NULL) {ierr = VecDestroy(&vtilde); CHKERRQ(ierr); vtilde = NULL;}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set m at t=0
 * @param[in] m0 density/image at t=0 (initial condition of forward
 * problem)
 *******************************************************************/
PetscErrorCode CLAIRE::SetInitialState(Vec m0) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m0 != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate state variable
    ierr = AllocateOnce(this->m_StateVariable, this->m_Opt, true, this->m_Opt->m_RegFlags.runinversion); CHKERRQ(ierr);

    // copy m_0 to m(t=0)
    ierr = this->m_StateVariable->SetFrame(m0, 0); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief get m at t=1
 * @param[out] m1 density/image at t=1 (solution of forward problem)
 *******************************************************************/
PetscErrorCode CLAIRE::GetFinalState(Vec m1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    IntType nt = this->m_Opt->m_Domain.nt;

    if (!this->m_Opt->m_RegFlags.runinversion) {
        nt = 0; // we did not store the time history
    }

    // copy m(t=1) to m_1
    ierr = this->m_StateVariable->GetFrame(m1, nt);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set lambda at t=1
 * @param[out] l1 adjoint variable at t=1 (final condition of adjoint
 * problem)
 *******************************************************************/
PetscErrorCode CLAIRE::SetFinalAdjoint(Vec l1) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(l1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    IntType nt = this->m_Opt->m_Domain.nt;

    // we do not store the time history for a gauss-newton approximation
    if (this->m_Opt->m_OptPara.method == GAUSSNEWTON) {
        nt = 0;
    }

    // allocate pointer if not done so already
    //ierr = AllocateOnce(this->m_AdjointVariable, this->m_Opt, true, true); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_AdjointVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);

    // copy l1 to lambda(t=1)
    ierr = this->m_AdjointVariable->GetFrame(l1, nt); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the forward problem (we assume the user has
 * set the velocity field)
 * @param[in] m0 density/image at t=0
 * @param[out] m1 density/image at t=1 (solution of transport
 * equation)
 *******************************************************************/
PetscErrorCode CLAIRE::SolveForwardProblem(Vec m1, Vec m0) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    DebugGPUStartEvent(__FUNCTION__);

    ierr = Assert(m0 != NULL, "null pointer"); CHKERRQ(ierr);

    // set initial condition
    ierr = AllocateOnce(this->m_TemplateImage, this->m_Opt, m0, true); CHKERRQ(ierr);
    ierr = this->m_TemplateImage->SetVector(m0); CHKERRQ(ierr);
    //this->m_TemplateImage = m0;

    // compute solution of state equation
    ierr = this->SolveStateEquation(); CHKERRQ(ierr);

    // only copy if output is necessary
    if (m1 != NULL) {
        ierr = this->GetFinalState(m1); CHKERRQ(ierr);
    }

    DebugGPUStopEvent();
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem (we assume the user has
 * set the velocity field)
 * @param[in] m1 density/image at t=1
 * @param[out] l0 adjoint variable at t=0 (solution of transport
 * equation)
 *******************************************************************/
PetscErrorCode CLAIRE::SolveAdjointProblem(Vec l0, Vec m1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    DebugGPUStartEvent(__FUNCTION__);

    ierr = Assert(m1 != NULL, "null pointer"); CHKERRQ(ierr);

    IntType nt = this->m_Opt->m_Domain.nt;

    // allocate state variable
    ierr = AllocateOnce(this->m_StateVariable, this->m_Opt, 0.0, true, true); CHKERRQ(ierr);
    
    // copy memory for m_1
    ierr = this->m_StateVariable->SetFrame(m1, nt); CHKERRQ(ierr);

    // compute solution of state equation
    ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);

    // copy memory for lambda0
    ierr = this->m_AdjointVariable->GetFrame(l0, 0); CHKERRQ(ierr);
    
    DebugGPUStopEvent();

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set state variable from externally
 *******************************************************************/
PetscErrorCode CLAIRE::SetStateVariable(Vec m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != NULL, "null pointer"); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("setting state variable"); CHKERRQ(ierr);
    }

    // we have to allocate the variable, because we delete it
    // at the end once we're done; since it comes from external
    // we need to make sure that we don't delete the external
    // pointer
    ierr = AllocateOnce(this->m_StateVariable, this->m_Opt, true, true); CHKERRQ(ierr);

    ierr = this->m_StateVariable->Copy(m); CHKERRQ(ierr);

    // if semi lagrangian pde solver is used,
    // we have to initialize it here
    /*if (this->m_Opt->m_PDESolver.type == SL) {
        ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
        // compute trajectory
        ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    }*/

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set state variable from externally
 *******************************************************************/
PetscErrorCode CLAIRE::GetStateVariable(Vec& m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m == NULL, "null pointer expected"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    m = *this->m_StateVariable;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set adjoint variable from externally
 *******************************************************************/
PetscErrorCode CLAIRE::SetAdjointVariable(Vec lambda) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda != NULL, "null pointer"); CHKERRQ(ierr);

    // we have to allocate the variable, because we delete it
    // at the end once we're done; since it comes from external
    // we need to make sure that we don't delete the external pointer
    ierr = AllocateOnce(this->m_AdjointVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);

    ierr = this->m_AdjointVariable->Copy(lambda); CHKERRQ(ierr);

    /*if (this->m_Opt->m_PDESolver.type == SL) {
        ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
        // compute trajectory
        ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }*/

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief get adjoint variable
 *******************************************************************/
PetscErrorCode CLAIRE::GetAdjointVariable(Vec& lambda) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda == NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    lambda = *this->m_AdjointVariable;
//    ierr = VecCopy(this->m_AdjointVariable, lambda); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief evaluate the l2 distance between m_R and m_1
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateDistanceMeasure(ScalarType* D) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    //ZeitGeist_define(EVAL_DIST);
    //ZeitGeist_tick(EVAL_DIST);
    
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
    
    // compute solution of state equation
    ierr = this->SolveStateEquation(); CHKERRQ(ierr);
    
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate distance measure
    if (this->m_Opt->m_Distance.reset) {
        ierr = this->SetupDistanceMeasure(); CHKERRQ(ierr);
        this->m_Opt->m_Distance.reset = false;
    }
    if (this->m_DistanceMeasure == NULL) {
        ierr = this->SetupDistanceMeasure(); CHKERRQ(ierr);
    }
    
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);

    // set state variable
    ierr = this->m_DistanceMeasure->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetWorkScaField(this->m_WorkScaField1); CHKERRQ(ierr);

    // evaluate distance measure
    ierr = this->m_DistanceMeasure->EvaluateFunctional(D); CHKERRQ(ierr);
    
    //ZeitGeist_tock(EVAL_DIST);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief evaluates the objective value
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateObjective(ScalarType* J, Vec v) {
    PetscErrorCode ierr = 0;
    ScalarType D = 0.0, R = 0.0;
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("evaluating objective"); CHKERRQ(ierr);
    }

    // allocate
    //ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
    //if (v) {
    //  ierr = Free(this->m_VelocityField); CHKERRQ(ierr);
    //  ierr = AllocateOnce(this->m_VelocityField, this->m_Opt, v); CHKERRQ(ierr);
    //} else {
      ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
      ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);
    //}
    
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // start timer
    ZeitGeist_define(EVAL_OBJ);
    ZeitGeist_tick(EVAL_OBJ);
    ierr = this->m_Opt->StartTimer(OBJEXEC); CHKERRQ(ierr);

    // set components of velocity field
    //ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // evaluate the regularization model
    ierr = this->EvaluateDistanceMeasure(&D); CHKERRQ(ierr);
    
    ZeitGeist_define(OBJ_REG);
    ZeitGeist_tick(OBJ_REG);
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (!this->m_VelocityIsZero) {
        // evaluate the regularization model
        ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
        ierr = this->m_Regularization->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
        ierr = this->m_Regularization->SetWorkScaField(this->m_WorkScaField1); CHKERRQ(ierr);
        //ierr = this->m_Regularization->SetDifferentiation(this->m_DifferentiationFD); CHKERRQ(ierr);
        ierr = this->m_Regularization->EvaluateFunctional(&R, this->m_VelocityField); CHKERRQ(ierr);
        //ierr = this->m_Regularization->SetDifferentiation(this->m_Differentiation); CHKERRQ(ierr);
    }
    ZeitGeist_tock(OBJ_REG);

    // add up the contributions
    *J = D + R;
    
    // store for access (e.g., used in coupling)
    this->m_Opt->m_Monitor.jval = *J;
    this->m_Opt->m_Monitor.dval = D;
    this->m_Opt->m_Monitor.rval = R;

    if (this->m_Opt->m_Verbosity > 1) {
        ss << "J(v) = D(v) + R(v) = " << std::scientific
           << this->m_Opt->m_Monitor.dval << " + "
           << this->m_Opt->m_Monitor.rval;
        ierr = DbgMsg1(ss.str()); CHKERRQ(ierr);
    }

    // stop timer
    ierr = this->m_Opt->StopTimer(OBJEXEC); CHKERRQ(ierr);
    ZeitGeist_tock(EVAL_OBJ);

    // increment counter for objective evaluations
    this->m_Opt->IncrementCounter(OBJEVAL);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief evaluates the reduced gradient of the lagrangian
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateGradient(Vec g, Vec v) {
    PetscErrorCode ierr = 0;
    ScalarType value, nvx1, nvx2, nvx3;
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("evaluating gradient"); CHKERRQ(ierr);
    }
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate
    //ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);

    // start timer
    ZeitGeist_define(EVAL_GRAD);
    ZeitGeist_tick(EVAL_GRAD);
    ierr = this->m_Opt->StartTimer(GRADEXEC); CHKERRQ(ierr);

    // parse input arguments
    //if (v != NULL) {
    //    ierr = Free(this->m_VelocityField); CHKERRQ(ierr);
    //    ierr = AllocateOnce(this->m_VelocityField, this->m_Opt, v); CHKERRQ(ierr);
    //} else {
      ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
    if (v != NULL) {
      ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);
    }

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = this->m_VelocityField->Norm(nvx1, nvx2, nvx3); CHKERRQ(ierr);
        ss  << "||v||_2 = (" << std::scientific
            << nvx1 << "," << nvx2 << "," << nvx3 << ")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }
    
    /*if (this->m_Opt->m_Verbosity > 1) {
      CFLStatKernel kernel;
      kernel.dt = 1.0/static_cast<ScalarType>(this->m_Opt->m_Domain.nt);
      //kernel.h = this->m_Opt->m_Domain.isize[0]*this->m_Opt->m_Domain.hx[0];
      kernel.h = 2.*M_PI/static_cast<ScalarType>(this->m_Opt->m_Domain.nx[0]);
      kernel.h *= 4; //this->m_Opt->m_Domain.isize[0];
      kernel.ng = static_cast<ScalarType>(this->m_Opt->m_Domain.ng);
      kernel.nl = this->m_Opt->m_Domain.nl;
      
      ierr = this->m_VelocityField->GetArraysRead(kernel.pV); CHKERRQ(ierr);
      
      ScalarType res = 0;
      
      kernel.CFLx(res);
      MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
      
      res /= static_cast<ScalarType>(this->m_Opt->m_Domain.ng);
      
      ss  << "Vx*dt > slab_size: " << res*100 << "%%";
      ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
      ss.clear(); ss.str(std::string());
      
      ierr = this->m_VelocityField->RestoreArrays(); CHKERRQ(ierr);
    }*/

    // compute solution of adjoint equation (i.e., \lambda(x,t))
    // and compute body force \int_0^1 grad(m)\lambda dt
    // which is assigned to work vecfield 2
    ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);
    
    ierr = this->m_WorkVecField2->DebugInfo("adjoint grad", __LINE__, __FILE__); CHKERRQ(ierr);

    // evaluate gradient of regularization model
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    
    if (this->m_VelocityIsZero) {
        // \vect{g}_v = \D{K}[\vect{b}]
        if (g != NULL) {
            ierr = this->m_WorkVecField2->GetComponents(g); CHKERRQ(ierr);
        }
    } else {
        switch (this->m_Opt->m_OptPara.gradtype) {
            case L2GRAD:
            {
                // evaluate l2 gradient
                ierr = this->EvaluateL2Gradient(g); CHKERRQ(ierr);
                break;
            }
/*
            case SGRAD:
            {
                // evaluate sobolev gradient
                ierr = this->EvaluateSobolevGradient(g, false); CHKERRQ(ierr);
                break;
            }
            case SYMSGRAD:
            {
                // evaluate sobolev gradient
                ierr = this->EvaluateSobolevGradient(g, true); CHKERRQ(ierr);
                break;
            }
*/
            default:
            {
                ierr = ThrowError("operator not implemented"); CHKERRQ(ierr);
                break;
            }
        }
    }

    // parse to output
    if (g != NULL) {
        // get and scale by lebesque measure
//        hd = this->m_Opt->GetLebesgueMeasure();
//        ierr = VecScale(g, hd); CHKERRQ(ierr);

        if (this->m_Opt->m_Verbosity > 2) {
            ierr = VecNorm(g, NORM_2, &value); CHKERRQ(ierr);
            ss << "||g||_2 = " << std::scientific << value;
            ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
    }
    
    if (this->m_Opt->m_Verbosity > 3) {
      if(this->m_StateVariable) this->m_StateVariable->DebugInfo("state       ",__LINE__,__FILE__);
      else DbgMsg3("state            : nullptr");
      if(this->m_AdjointVariable) this->m_AdjointVariable->DebugInfo("adjoint     ",__LINE__,__FILE__);
      else DbgMsg3("adjoint          : nullptr");
      if(this->m_IncStateVariable) this->m_IncStateVariable->DebugInfo("inc state   ",__LINE__,__FILE__);
      else DbgMsg3("inc state        : nullptr");
      if(this->m_IncAdjointVariable) this->m_IncAdjointVariable->DebugInfo("inc adjoint ",__LINE__,__FILE__);
      else DbgMsg3("inc adjoint      : nullptr");
      if(this->m_ReferenceImage) this->m_ReferenceImage->DebugInfo("reference   ",__LINE__,__FILE__);
      else DbgMsg3("reference        : nullptr");
      if(this->m_TemplateImage) this->m_TemplateImage->DebugInfo("template    ",__LINE__,__FILE__);
      else DbgMsg3("template         : nullptr");
      if(this->m_VelocityField) this->m_VelocityField->DebugInfo("velocity    ",__LINE__,__FILE__);
      else DbgMsg3("velocity         : nullptr");
      if(this->m_IncVelocityField) this->m_IncVelocityField->DebugInfo("inc velocity",__LINE__,__FILE__);
      else DbgMsg3("inc velocity     : nullptr");
    }

    // stop timer
    ierr = this->m_Opt->StopTimer(GRADEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(GRADEVAL);
    ZeitGeist_tock(EVAL_GRAD);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief evaluates the reduced gradient of the lagrangian (l2)
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateL2Gradient(Vec g) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // evaluate / apply gradient operator for regularization
    //ierr = this->m_Regularization->EvaluateGradient(this->m_WorkVecField1, this->m_VelocityField); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Diff.diffReg == FINITE && 
        (this->m_Opt->m_RegNorm.type == H1 || this->m_Opt->m_RegNorm.type == H1SN)) {
      ierr = AllocateOnce(this->m_DifferentiationFD, this->m_Opt); CHKERRQ(ierr);
      ierr = this->m_Regularization->SetDifferentiation(this->m_DifferentiationFD); CHKERRQ(ierr);
    }
    ierr = this->m_Regularization->EvaluateGradient(this->m_WorkVecField1, this->m_VelocityField); CHKERRQ(ierr);
    if (this->m_Opt->m_Diff.diffReg == FINITE && 
        (this->m_Opt->m_RegNorm.type == H1 || this->m_Opt->m_RegNorm.type == H1SN)) {
      ierr = AllocateOnce<DifferentiationSM>(this->m_Differentiation, this->m_Opt); CHKERRQ(ierr);
      ierr = this->m_Regularization->SetDifferentiation(this->m_Differentiation); CHKERRQ(ierr);
    }
    
    ierr = this->m_WorkVecField1->DebugInfo("reg grad", __LINE__, __FILE__); CHKERRQ(ierr);
    if (this->m_Opt->m_ReadWriteFlags.iterates) {
        std::stringstream ss;
        std::string ext;
        IntType iter = 0;
        std::string ver = "cpu-";
#ifdef REG_HAS_CUDA
        ver = "gpu-";
#endif
        
        // parse extension
        ext = this->m_Opt->m_FileNames.extension;
        // log objective values
        iter = this->m_Opt->GetCounter(ITERATIONS);
        
        ss  << ver << "reg-grad-field-i=" << std::setw(3) << std::setfill('0') << iter << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1, ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
        ss  << ver << "adjoint-grad-field-i=" << std::setw(3) << std::setfill('0') << iter << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField2, ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }
    
    // \vect{g}_v = \beta_v \D{A}[\vect{v}] + \D{K}[\vect{b}]
    ierr = this->m_WorkVecField1->AXPY(1.0, this->m_WorkVecField2); CHKERRQ(ierr);
    
    ierr = this->m_WorkVecField1->DebugInfo("gradient    ", __LINE__, __FILE__); CHKERRQ(ierr);
    if (this->m_Opt->m_ReadWriteFlags.iterates) {
        std::stringstream ss;
        std::string ext;
        IntType iter = 0;
        std::string ver = "cpu-";
#ifdef REG_HAS_CUDA
        ver = "gpu-";
#endif
        
        // parse extension
        ext = this->m_Opt->m_FileNames.extension;
        // log objective values
        iter = this->m_Opt->GetCounter(ITERATIONS);
        
        ss  << ver << "grad-field-i=" << std::setw(3) << std::setfill('0') << iter << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1, ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }
    
    // copy
    if (g != NULL) {
      ierr = this->m_WorkVecField1->GetComponents(g); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief evaluates the reduced gradient of the lagrangian (in
 * sobolev space incuded by regularization operator)
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateSobolevGradient(Vec g, bool flag) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // evaluate / apply gradient operator for regularization
    ierr = this->m_Regularization->ApplyInverse(this->m_WorkVecField1, this->m_WorkVecField2, flag); CHKERRQ(ierr);

    // \vect{g}_v = \vect{v} + (\beta_v \D{A})^{-1}\D{K}[\vect{b}]
    ierr = this->m_WorkVecField1->AXPY(1.0, this->m_VelocityField); CHKERRQ(ierr);

    // copy to output
    if (g != NULL) {
      ierr = this->m_WorkVecField1->GetComponents(g); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief applies the hessian to a vector
 * @param[in] vtilde incremental velocity field
 * @param[in] scale flag to switch on scaling by lebesgue measure
 * @param[out] Hvtilde hessian applied to vector
 *******************************************************************/
PetscErrorCode CLAIRE::HessianMatVec(Vec Hvtilde, Vec vtilde, bool scale) {
    PetscErrorCode ierr = 0;
    ScalarType hd; //, gamma;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("computing hessian matvec"); CHKERRQ(ierr);
    }
    
    if (this->m_Opt->m_Verbosity > 3) {
      if(this->m_StateVariable) this->m_StateVariable->DebugInfo("state       ",__LINE__,__FILE__);
      else DbgMsg3("state            : nullptr");
      if(this->m_AdjointVariable) this->m_AdjointVariable->DebugInfo("adjoint     ",__LINE__,__FILE__);
      else DbgMsg3("adjoint          : nullptr");
      if(this->m_IncStateVariable) this->m_IncStateVariable->DebugInfo("inc state   ",__LINE__,__FILE__);
      else DbgMsg3("inc state        : nullptr");
      if(this->m_IncAdjointVariable) this->m_IncAdjointVariable->DebugInfo("inc adjoint ",__LINE__,__FILE__);
      else DbgMsg3("inc adjoint      : nullptr");
      if(this->m_ReferenceImage) this->m_ReferenceImage->DebugInfo("reference   ",__LINE__,__FILE__);
      else DbgMsg3("reference        : nullptr");
      if(this->m_TemplateImage) this->m_TemplateImage->DebugInfo("template    ",__LINE__,__FILE__);
      else DbgMsg3("template         : nullptr");
      if(this->m_VelocityField) this->m_VelocityField->DebugInfo("velocity    ",__LINE__,__FILE__);
      else DbgMsg3("velocity         : nullptr");
      if(this->m_IncVelocityField) this->m_IncVelocityField->DebugInfo("inc velocity",__LINE__,__FILE__);
      else DbgMsg3("inc velocity     : nullptr");
    }

    ZeitGeist_define(EVAL_HESS);
    ZeitGeist_tick(EVAL_HESS);
    ierr = this->m_Opt->StartTimer(HMVEXEC); CHKERRQ(ierr);
    
    /*{
    VecField *g_vec = new VecField(m_Opt, vtilde);
    
    TwoLevelFFT op(this->m_Opt);
    op.Restrict(this->m_WorkVecField1, g_vec);
    op.Prolong(g_vec, this->m_WorkVecField1);
        
    delete g_vec;
    }*/
    
    /*VecCopy(vtilde, Hvtilde);
    hd = this->m_Opt->GetLebesgueMeasure();
    VecScale(Hvtilde, hd);*/

    // switch between hessian operators
    switch (this->m_Opt->m_KrylovMethod.matvectype) {
        case DEFAULTMATVEC:
        {
            // apply hessian H to \tilde{v}
            ierr = this->HessMatVec(Hvtilde, vtilde); CHKERRQ(ierr);
            //ierr = VecCopy(vtilde, Hvtilde); CHKERRQ(ierr);
            break;
        }
        case H0MATVEC:
        {
            // apply hessian H to \tilde{v}
            ierr = this->H0HessMatVec(Hvtilde, vtilde); CHKERRQ(ierr);
            //ierr = VecCopy(vtilde, Hvtilde); CHKERRQ(ierr);
            break;
        }
        case PRECONDMATVEC:
        {
            // apply analytically preconditioned hessian H to \tilde{v}
            ierr = this->PrecondHessMatVec(Hvtilde, vtilde); CHKERRQ(ierr);
            break;
        }
        case PRECONDMATVECSYM:
        {
            // apply analytically preconditioned hessian H to \tilde{v}
            // hessian operator is symmetrized
            ierr = this->PrecondHessMatVecSym(Hvtilde, vtilde); CHKERRQ(ierr);
            //ierr = this->SymTwoLevelHessMatVec(Hvtilde, vtilde); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("operator not implemented"); CHKERRQ(ierr);
            break;
        }
    }
    
    /*{
    VecField *g_vec = new VecField(m_Opt, Hvtilde);
    
    TwoLevelFFT op(this->m_Opt);
    op.Restrict(this->m_WorkVecField1, g_vec);
    op.Prolong(g_vec, this->m_WorkVecField1);
        
    delete g_vec;
    }*/

    if (Hvtilde != NULL) {
        // TODO @ Andreas: fix for two-level precond
        // scale by lebesgue measure
        if (scale == false) {
            hd = this->m_Opt->GetLebesgueMeasure();
            ierr = VecScale(Hvtilde, 1/hd); CHKERRQ(ierr);
        }

//        gamma = this->m_Opt->m_KrylovMethod.hessshift;
//        if (gamma > 0.0) {
//            ierr = VecAXPY(Hvtilde, gamma, vtilde); CHKERRQ(ierr);
//        }
    }

    // stop hessian matvec timer
    ierr = this->m_Opt->StopTimer(HMVEXEC); CHKERRQ(ierr);
    ZeitGeist_tock(EVAL_HESS);

    // increment matvecs
    this->m_Opt->IncrementCounter(HESSMATVEC);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief applies inverse of H(v=0) (used as preconditioner for our problem)
 * This is a PCG solver solving the H0x = r system and prolonging and restricting it.
 *******************************************************************/
PetscErrorCode CLAIRE::ApplyInvHessian(Vec precx, Vec x, VecField* gradM, bool first, bool twolevel, Preprocessing* preproc) {
    PetscErrorCode ierr = 0;
    H0PrecondKernel kernel;
    const ScalarType *ptr = nullptr;
    DifferentiationSM* diff;
    
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    ZeitGeist_define(PC_H0);
    ZeitGeist_tick(PC_H0);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(gradM != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField3, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkScaField4, this->m_Opt); CHKERRQ(ierr);
    
    //ScalarType beta = sqrt(this->m_Opt->m_RegNorm.beta[0]);
    
    VecField *vecR, *vecP, *vecM, *vecX, *vecV;
    vecR = this->m_WorkVecField3;
    vecP = this->m_WorkVecField2;
    vecM = this->m_WorkVecField1;
    vecV = this->m_WorkVecField4;
    vecX = nullptr;
    //vecX = this->m_WorkVecField1;
    
    ierr = VecCopy(x, precx); CHKERRQ(ierr);
    ierr = AllocateOnce(vecX, this->m_Opt, precx); CHKERRQ(ierr);
    
    //FILE *handle = fopen("inreg_res.txt", "a");
    //ScalarType norm;
    
    //int nit =  this->m_Opt->GetCounter(ITERATIONS) - 1;
    //int kit = this->m_Opt->GetCounter(HESSMATVEC);
    
    //vecX->Norm2(norm);
    
    //fprintf(handle, "%i,%i,%e", nit, kit, sqrt(norm));
    
    ScalarType beta, betav, betaw;
    
    beta = this->m_Opt->m_RegNorm.beta[0];
    betav = this->m_Opt->m_RegNorm.beta[0];
    betaw = this->m_Opt->m_RegNorm.beta[2];
    if (this->m_Opt->m_RegNorm.beta[3] > 0 && this->m_Opt->m_RegNorm.beta[3] > beta) {
      beta = this->m_Opt->m_RegNorm.beta[3];
    }
    
#if 0
     ierr = this->m_Differentiation->InvRegLapOp(vecX, vecX, false, betav); CHKERRQ(ierr);
     
     /*TwoLevelFFT op(this->m_Opt);
     op.Restrict(vecR, vecX);
     op.Prolong(vecX, vecR);*/
     
     
     {
     ScalarType *ptr[3];
     
     ScalarType n2, n3, n4;
     
     kernel.nl = this->m_Opt->m_Domain.nl;
     vecX->GetArraysRead(kernel.pGmt);
     kernel.Norm(n3);
     MPI_Allreduce(MPI_IN_PLACE, &n3, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
     
     vecX->RestoreArrays();
     
     vecX->GetArraysReadWrite(ptr);
     
     ComplexType *cptr;
     
     cptr = this->m_Opt->m_FFT.fft->m_WorkSpace;
     
     vecM->SetValue(-100000);
     
     this->m_Opt->m_FFT.fft->m_WorkVecField = vecM;
     
     this->m_Opt->m_FFT.fft->FFT_R2C(ptr[0], cptr);
     this->m_Opt->m_FFT.fft->Scale(cptr, 1./static_cast<ScalarType>(this->m_Opt->m_Domain.ng));
     this->m_Opt->m_FFT.fft->Norm(norm, cptr);
     //this->m_Opt->m_FFT.fft->LowPassFilter(cptr, 0.5);
     this->m_Opt->m_FFT.fft->FFT_C2R(cptr, ptr[0]);
     
     this->m_Opt->m_FFT.fft->FFT_R2C(ptr[1], cptr);
     //this->m_Opt->m_FFT.fft->LowPassFilter(cptr, 0.5);
     this->m_Opt->m_FFT.fft->Scale(cptr, 1./static_cast<ScalarType>(this->m_Opt->m_Domain.ng));
     this->m_Opt->m_FFT.fft->FFT_C2R(cptr, ptr[1]);
     
     this->m_Opt->m_FFT.fft->FFT_R2C(ptr[2], cptr);
     //this->m_Opt->m_FFT.fft->LowPassFilter(cptr, 0.5);
     this->m_Opt->m_FFT.fft->Scale(cptr, 1./static_cast<ScalarType>(this->m_Opt->m_Domain.ng));
     this->m_Opt->m_FFT.fft->FFT_C2R(cptr, ptr[2]);
     
     vecX->RestoreArrays();
    
     kernel.nl = this->m_Opt->m_Domain.nl;
     vecX->GetArraysRead(kernel.pGmt);
     kernel.Norm(n2);
     MPI_Allreduce(MPI_IN_PLACE, &n2, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
     
     vecX->RestoreArrays();
     //vecX->Norm(n2, n3, n4);
     
     printf("%e %e %e\n", n3, n2, norm*this->m_Opt->m_Domain.ng);
    }

#else

    /*ierr = this->m_Differentiation->InvRegLapOp(vecR, vecX, false, betav); CHKERRQ(ierr);
    
    ScalarType norm2;
    vecR->Norm2(norm2);
    fprintf(handle, ",%e", sqrt(norm2));
    
    if (sqrt(norm2/norm) < this->m_Opt->m_KrylovMethod.reltol) {
      vecX->Copy(vecR);
    } else {*/
    
    bool solve_h0 = true;
    bool sym_h0 = true;
    
    kernel.beta = 1; //betav;
    
    diff = (DifferentiationSM*)this->m_Differentiation;
    
    kernel.nl = this->m_Opt->m_Domain.nl;
    
    IntType *nx = this->m_Opt->m_FFT.nx;
    
    kernel.diag = beta*((nx[0]-2)*(nx[0]-1) + (nx[1]-2)*(nx[1]-1) + (nx[2]-2)*(nx[2]-1))/12.;
    
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }
    
    ScalarType normref, cg_a, cg_b, cg_r, cg_p;
    
    if (first) {
      if (this->m_GradientState) {
        ierr = gradM->Copy(this->m_GradientState[this->m_Opt->m_Domain.nt]); CHKERRQ(ierr);
      } else {
        ierr = this->m_StateVariable->GetArrayRead(ptr, 0, this->m_Opt->m_Domain.nt); CHKERRQ(ierr);
        if (this->m_Opt->m_Diff.diffPDE == FINITE) {
          ierr = this->m_DifferentiationFD->Gradient(gradM, ptr); CHKERRQ(ierr);
        } else {
          ierr = this->m_Differentiation->Gradient(gradM, ptr); CHKERRQ(ierr);
        }
        ierr = this->m_StateVariable->RestoreArray();
      }
      
      if (twolevel) {
        TwoLevelFFT twolevel_op(this->m_Opt);
        ierr = twolevel_op.Restrict(gradM, gradM); CHKERRQ(ierr);
        /*TwoLevelFinite twolevel_op(this->m_Opt);
        ierr = vecM->Copy(gradM);
        ierr = twolevel_op.Restrict(gradM, vecM); CHKERRQ(ierr);*/
      }
    }
    
    if (!twolevel) {
      ScalarType norm, norm2;
      ierr = vecX->Norm2(norm); CHKERRQ(ierr);
      ierr = this->m_Differentiation->InvRegLapOp(vecR, vecX, sym_h0, betav); CHKERRQ(ierr);
      ierr = vecR->Norm2(norm2); CHKERRQ(ierr);
      if (sqrt(norm2/norm) < this->m_Opt->m_KrylovMethod.reltol) {
        ierr = vecX->Copy(vecR); CHKERRQ(ierr);
        solve_h0 = false;
        ZeitGeist_define(PC_H0_INVREG);
        ZeitGeist_add(PC_H0_INVREG,1);
      }
    } else {
      if (this->m_Opt->m_Verbosity > 2) {
          ScalarType norm;
          ierr = vecX->Norm2(norm);
          std::stringstream ss;
          ss << "KSP res: " << sqrt(norm);
          ierr = DbgMsgCall(ss.str()); CHKERRQ(ierr);
      }
      
      TwoLevelRegFFT twolevel_op(this->m_Opt, betav, this->m_Opt->m_KrylovMethod.reltol, sym_h0);      
      //TwoLevelFinite twolevel_op(this->m_Opt);      
      this->m_Opt->m_FFT.fft->m_WorkVecField = vecM;
      
      ierr = twolevel_op.Restrict(vecR, vecX); CHKERRQ(ierr);
      
      solve_h0 = twolevel_op.restricted;
      //solve_h0 = true; //twolevel_op.restricted;
      
      if (!solve_h0) {
        ierr = twolevel_op.Prolong(vecX, vecR); CHKERRQ(ierr);
        ZeitGeist_define(PC_H0_INVREG);
        ZeitGeist_add(PC_H0_INVREG,1);
      }
      
      //ierr = this->m_Differentiation->InvRegLapOp(vecR, vecR, false, betav); CHKERRQ(ierr);
    }
    
    if (solve_h0) {
    
    if (!twolevel) {
      ierr = gradM->GetArraysRead(kernel.pGmt); CHKERRQ(ierr);
    } else {
      ierr = gradM->GetArraysRead(kernel.pGmt); CHKERRQ(ierr);
      ierr = diff->SetFFT(&this->m_Opt->m_FFT_coarse); CHKERRQ(ierr);
      kernel.nl = this->m_Opt->m_FFT_coarse.isize[0]*this->m_Opt->m_FFT_coarse.isize[1]*this->m_Opt->m_FFT_coarse.isize[2];
    }
    
    ierr = vecX->Copy(vecR); CHKERRQ(ierr);
    //ierr = vecX->SetValue(0);
    
    IntType innerloop = this->m_Opt->m_KrylovMethod.pcmaxit;
    
    ScalarType cg_eps;
    
    cg_eps = this->m_Opt->m_KrylovMethod.pctolint[0]*this->m_Opt->m_KrylovMethod.reltol;
    
    ierr = this->m_WorkScaField1->GetArrayWrite(kernel.pWS); CHKERRQ(ierr);
    
    normref = 1.;
    int rval;

    {
      if (sym_h0) {
        ierr = diff->InvRegLapOp(vecV, vecX, true, beta); CHKERRQ(ierr);
      } else {
        vecV = vecX;
      }
      
      ierr = vecV->GetArraysReadWrite(kernel.pVhat); CHKERRQ(ierr);
      ierr = vecM->GetArraysReadWrite(kernel.pM); CHKERRQ(ierr);
      ierr = kernel.gMgMT(); CHKERRQ(ierr);
      ierr = vecM->RestoreArrays(); CHKERRQ(ierr);
      ierr = vecV->RestoreArrays(); CHKERRQ(ierr);
      
      if (!sym_h0) {
        ierr = diff->InvRegLerayOp(vecM, vecM, betav, betaw, beta); CHKERRQ(ierr);
      } else {
        ierr = diff->LerayOperator(vecM, vecM, betav, betaw); CHKERRQ(ierr);
        ierr = diff->InvRegLapOp(vecM, vecM, true, beta); CHKERRQ(ierr);
      }
      
      ierr = vecX->GetArraysReadWrite(kernel.pVhat); CHKERRQ(ierr);
      ierr = vecM->GetArraysReadWrite(kernel.pM); CHKERRQ(ierr);
      ierr = vecR->GetArraysReadWrite(kernel.pRes); CHKERRQ(ierr);
      ierr = vecP->GetArraysReadWrite(kernel.pP); CHKERRQ(ierr);
      ierr = kernel.res(cg_r);
      ierr = vecX->RestoreArrays(); CHKERRQ(ierr);
      ierr = vecM->RestoreArrays(); CHKERRQ(ierr);
      ierr = vecR->RestoreArrays(); CHKERRQ(ierr);
      ierr = vecP->RestoreArrays(); CHKERRQ(ierr);
      
      rval = MPI_Allreduce(MPI_IN_PLACE, &cg_r, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
      ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
  
      normref = sqrt(cg_r);
      
      if (this->m_Opt->m_Verbosity > 2) {
        std::stringstream ss;
        ss << "PC res: " << sqrt(cg_r);
        ss << ", " << sqrt(cg_r)/normref;
        ierr = DbgMsgCall(ss.str()); CHKERRQ(ierr);
      }
      int i;
      for (i = 0; i<innerloop; ++i) {
        if (sym_h0) {
          ierr = diff->InvRegLapOp(vecV, vecP, true, beta); CHKERRQ(ierr);
        } else {
          vecV = vecP;
        }
        
        ierr = vecV->GetArraysReadWrite(kernel.pVhat); CHKERRQ(ierr);
        ierr = vecM->GetArraysReadWrite(kernel.pM); CHKERRQ(ierr);
        ierr = kernel.gMgMT(); CHKERRQ(ierr);
        ierr = vecM->RestoreArrays(); CHKERRQ(ierr);
        ierr = vecV->RestoreArrays(); CHKERRQ(ierr);
        
        if (!sym_h0) {
          ierr = diff->InvRegLerayOp(vecM, vecM, betav, betaw, beta); CHKERRQ(ierr);
        } else {
          ierr = diff->LerayOperator(vecM, vecM, betav, betaw); CHKERRQ(ierr);
          ierr = diff->InvRegLapOp(vecM, vecM, true, beta); CHKERRQ(ierr);
        }
        
        ierr = vecM->GetArraysReadWrite(kernel.pM); CHKERRQ(ierr);
        ierr = vecP->GetArraysReadWrite(kernel.pP); CHKERRQ(ierr);
        ierr = kernel.pTAp(cg_p);
        ierr = vecM->RestoreArrays(); CHKERRQ(ierr);
        ierr = vecP->RestoreArrays(); CHKERRQ(ierr);
        
        rval = MPI_Allreduce(MPI_IN_PLACE, &cg_p, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
        
        cg_a = cg_r/cg_p;
        ierr = vecX->GetArraysReadWrite(kernel.pVhat); CHKERRQ(ierr);
        ierr = vecM->GetArraysReadWrite(kernel.pM); CHKERRQ(ierr);
        ierr = vecR->GetArraysReadWrite(kernel.pRes); CHKERRQ(ierr);
        ierr = vecP->GetArraysReadWrite(kernel.pP); CHKERRQ(ierr);
        ierr = kernel.CGres(cg_a);
        ierr = vecX->RestoreArrays(); CHKERRQ(ierr);
        ierr = vecM->RestoreArrays(); CHKERRQ(ierr);
        ierr = vecR->RestoreArrays(); CHKERRQ(ierr);
        ierr = vecP->RestoreArrays(); CHKERRQ(ierr);
        
        rval = MPI_Allreduce(MPI_IN_PLACE, &cg_a, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
        
        if (sqrt(cg_a) < cg_eps*normref) {
          cg_r = cg_a;
          ++i;
          break;
        }
        cg_b = cg_a/cg_r;
        cg_r = cg_a;
        ierr = vecR->GetArraysReadWrite(kernel.pRes); CHKERRQ(ierr);
        ierr = vecP->GetArraysReadWrite(kernel.pP); CHKERRQ(ierr);
        ierr = kernel.CGp(cg_b);
        ierr = vecR->RestoreArrays(); CHKERRQ(ierr);
        ierr = vecP->RestoreArrays(); CHKERRQ(ierr);
        
        rval = MPI_Allreduce(MPI_IN_PLACE, &cg_b, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
      }
      
      if (this->m_Opt->m_Verbosity > 2)  {
        std::stringstream ss;
        ss << "PC " << (twolevel?"coarse: ": "fine: ") << i;
        ss << " ,kernel size: " << kernel.nl;
        ierr = DbgMsgCall(ss.str()); CHKERRQ(ierr);
      }
      
      if (twolevel) {
        ZeitGeist_define(PC_H0_coarse);
        ZeitGeist_add(PC_H0_coarse, i);
      } else {
        ZeitGeist_define(PC_H0_fine);
        ZeitGeist_add(PC_H0_fine,i);
      }
      
      if (twolevel) {
        ierr = diff->SetFFT(&this->m_Opt->m_FFT); CHKERRQ(ierr);
        
        TwoLevelRegFFT twolevel_op(this->m_Opt, betav);
        twolevel_op.restricted = true;
        ierr = twolevel_op.Prolong(vecX, vecX);
        
        /*TwoLevelFinite twolevel_op(this->m_Opt);
        ierr = vecR->Copy(vecX); CHKERRQ(ierr);
        ierr = twolevel_op.Prolong(vecX, vecR);*/
              
        ierr = gradM->RestoreArrays(); CHKERRQ(ierr);
        
        kernel.nl = this->m_Opt->m_Domain.nl;
        
        cg_eps = this->m_Opt->m_KrylovMethod.pctolint[2];
      } else {
        ierr = gradM->RestoreArrays(); CHKERRQ(ierr);
      }
    }
    
    if (this->m_Opt->m_Verbosity > 2) {
      std::stringstream ss;
      ss << "PC res: " << sqrt(cg_r) << ", " << sqrt(cg_r)/normref;
      ierr = DbgMsgCall(ss.str()); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();
    }
    if (this->m_Opt->m_Verbosity > 2) {
        ScalarType norm;
        ierr = vecX->Norm2(norm);
        std::stringstream ss;
        ss << "KSP res: " << sqrt(norm);
        ierr = DbgMsgCall(ss.str()); CHKERRQ(ierr);
    }
    
    ierr = this->m_WorkScaField1->RestoreArray(); CHKERRQ(ierr);
    
    }
  //}
#endif
    
    if (sym_h0) {
      ierr = diff->InvRegLapOp(vecX, vecX, true, betav); CHKERRQ(ierr);
    }

    //vecX->Norm2(norm);
    
    //fprintf(handle, ",%e\n", sqrt(norm));
    //fclose(handle);
        
    //ierr = this->m_WorkVecField1->GetComponents(precx); CHKERRQ(ierr);
    ierr = Free(vecX); CHKERRQ(ierr);
    
    ZeitGeist_tock(PC_H0);
  
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief applies the hessian to a vector (default way of doing this)
 *******************************************************************/
PetscErrorCode CLAIRE::H0HessMatVec(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    H0PrecondKernel kernel;
    ScalarType hd;
    const ScalarType *ptr = nullptr;
    
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
  
    ierr = AllocateOnce(this->m_IncVelocityField, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    
    hd  = this->m_Opt->GetLebesgueMeasure();
        
    kernel.nl = this->m_Opt->m_Domain.nl;
    
    if (this->m_GradientState) {
      ierr = this->m_GradientState[0]->GetArraysRead(kernel.pGmt); CHKERRQ(ierr);
    } else {
      ierr = this->m_StateVariable->GetArrayRead(ptr, 0, 0); CHKERRQ(ierr);
      ierr = this->m_Differentiation->Gradient(this->m_WorkVecField1, ptr); CHKERRQ(ierr);
      ierr = this->m_StateVariable->RestoreArray();
      ierr = this->m_WorkVecField1->GetArraysRead(kernel.pGmt); CHKERRQ(ierr);
    }
    
    // parse input
    if (vtilde != NULL) {
        ierr = this->m_IncVelocityField->SetVector(vtilde); CHKERRQ(ierr);
    }
    
    ierr = this->m_IncVelocityField->GetArraysReadWrite(kernel.pVhat); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArraysReadWrite(kernel.pM); CHKERRQ(ierr);
    ierr = kernel.gMgMT(); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->RestoreArrays(); CHKERRQ(ierr);
    
    if (this->m_GradientState) {
      ierr = this->m_GradientState[0]->RestoreArrays(); CHKERRQ(ierr);
    } else {
      ierr = this->m_WorkVecField1->RestoreArrays(); CHKERRQ(ierr);
    }
    
    // apply K[\tilde{b}]
    ierr = this->ApplyProjection(); CHKERRQ(ierr);
    // scale result by hd
    ierr = this->m_WorkVecField2->Scale(hd); CHKERRQ(ierr);
        
    // apply 2nd variation of regularization model to
    // incremental control variable: \beta*\D{A}[\vect{\tilde{v}}]
    ierr = this->m_Regularization->HessianMatVec(this->m_WorkVecField1, this->m_IncVelocityField); CHKERRQ(ierr);
    
    // \D{H}\vect{\tilde{v}} = \beta*\D{A}[\vect{\tilde{v}}] + \D{K}[\vect{\tilde{b}}]
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    ierr = this->m_WorkVecField1->AXPY(1.0, this->m_WorkVecField2); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 3) {
      if(this->m_WorkVecField1) this->m_WorkVecField1->DebugInfo("hessian     ",__LINE__,__FILE__);
      else DbgMsg3("hessian          : nullptr");
    }

    // pass to output
    if (Hvtilde != NULL) {
        ierr = this->m_WorkVecField1->GetComponents(Hvtilde); CHKERRQ(ierr);
    }
    
    ierr = this->m_VelocityField->DebugInfo("velocity", __LINE__, __FILE__); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->DebugInfo("inc velocity", __LINE__, __FILE__); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief applies the hessian to a vector (default way of doing this)
 *******************************************************************/
PetscErrorCode CLAIRE::HessMatVec(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("regular hessian matvec"); CHKERRQ(ierr);
    }

    this->m_Opt->Enter(__func__);
    
    VecField* Hv = nullptr;
    
    if (Hvtilde != NULL) {
      ierr = AllocateOnce(Hv, this->m_Opt, Hvtilde); CHKERRQ(ierr);
    } else {
      ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
      Hv = this->m_WorkVecField1;
    }
    

    // allocate container for incremental velocity field
    //ierr = AllocateOnce(this->m_IncVelocityField, this->m_Opt); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }
    
    // parse input
    //if (vtilde != NULL) {
    //    ierr = Free(this->m_IncVelocityField); CHKERRQ(ierr);
    //    ierr = AllocateOnce(this->m_IncVelocityField, this->m_Opt, vtilde); CHKERRQ(ierr);
    //} else {
      ierr = AllocateOnce(this->m_IncVelocityField, this->m_Opt); CHKERRQ(ierr);
    if (vtilde != NULL) {
      ierr = this->m_IncVelocityField->SetVector(vtilde); CHKERRQ(ierr);
    }
    
    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t)
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // compute incremental body force
//    ierr = this->ComputeIncBodyForce(); CHKERRQ(ierr);

    // apply 2nd variation of regularization model to
    // incremental control variable: \beta*\D{A}[\vect{\tilde{v}}]
    //ierr = this->m_Regularization->HessianMatVec(Hv, this->m_IncVelocityField); CHKERRQ(ierr);
    
    /*if (this->m_Opt->m_Diff.diffReg == FINITE && 
        (this->m_Opt->m_RegNorm.type == H1 || this->m_Opt->m_RegNorm.type == H1SN)) {
      ierr = AllocateOnce(this->m_DifferentiationFD, this->m_Opt); CHKERRQ(ierr);
      ierr = this->m_Regularization->SetDifferentiation(this->m_DifferentiationFD); CHKERRQ(ierr);
    }*/
    ierr = this->m_Regularization->HessianMatVec(Hv, this->m_IncVelocityField); CHKERRQ(ierr);
    /*if (this->m_Opt->m_Diff.diffReg == FINITE && 
        (this->m_Opt->m_RegNorm.type == H1 || this->m_Opt->m_RegNorm.type == H1SN)) {
      ierr = AllocateOnce<DifferentiationSM>(this->m_Differentiation, this->m_Opt); CHKERRQ(ierr);
      ierr = this->m_Regularization->SetDifferentiation(this->m_Differentiation); CHKERRQ(ierr);
    }*/
    

    // \D{H}\vect{\tilde{v}} = \beta*\D{A}[\vect{\tilde{v}}] + \D{K}[\vect{\tilde{b}}]
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    ierr = Hv->AXPY(1.0, this->m_WorkVecField2); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 3) {
      if(Hv) Hv->DebugInfo("hessian     ",__LINE__,__FILE__);
      else DbgMsg3("hessian          : nullptr");
    }

    // pass to output
    if (Hvtilde != NULL) {
        Free(Hv);
        //ierr = this->m_WorkVecField1->GetComponents(Hvtilde); CHKERRQ(ierr);
    }
    
    ierr = this->m_VelocityField->DebugInfo("velocity", __LINE__, __FILE__); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->DebugInfo("inc velocity", __LINE__, __FILE__); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief applies the analytically (spectrally) preconditioned
 * hessian, i.e.
 * P(H[\tilde{v}]) = (\beta A)^{-1}(\beta A[\tilde{v}] + b[\tilde{v}])
 *                 = \tilde{v} + (\beta A)^{-1}(b[\tilde{v}])
 * it is important to note, that this matrix is no longer symmetric;
 * we therefore can't use pcg
 *******************************************************************/
PetscErrorCode CLAIRE::PrecondHessMatVec(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType hd;
    PetscFunctionBegin;
    
    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("preconditioned hessian matvec"); CHKERRQ(ierr);
    }

    this->m_Opt->Enter(__func__);

    // allocate container for incremental velocity field
    ierr = AllocateOnce(this->m_IncVelocityField, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // parse input
    if (vtilde != NULL) {
        ierr = this->m_IncVelocityField->SetVector(vtilde); CHKERRQ(ierr);
    }

    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t) and incremental body force
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // apply inverse of 2nd variation of regularization model to
    // incremental body force: (\beta \D{A})^{-1}\D{K}[\vect{\tilde{b}}]
    ierr = this->m_Regularization->ApplyInverse(this->m_WorkVecField1, this->m_WorkVecField2, false); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \vect{\tilde{v}} + (\beta \D{A})^{-1} \D{K}[\vect{\tilde{b}}]
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    hd = this->m_Opt->GetLebesgueMeasure();
    ierr = this->m_WorkVecField2->WAXPY(hd, this->m_IncVelocityField, this->m_WorkVecField1); CHKERRQ(ierr);

    // pass to output
    if (Hvtilde != NULL) {
        ierr = this->m_WorkVecField2->GetComponents(Hvtilde); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief applies the analytically (spectrally) preconditioned
 * hessian, i.e.
 * P(H[\tilde{v}]) = \tilde{v} + P^{1/2}(b[\tilde{v}])P^{1/2}
 * it is important to note, that this matrix is no longer symmetric;
 * we therefore can't use pcg
 *******************************************************************/
PetscErrorCode CLAIRE::PrecondHessMatVecSym(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType hd;
    PetscFunctionBegin;
    
    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("symmetric preconditioned hessian matvec"); CHKERRQ(ierr);
    }

    this->m_Opt->Enter(__func__);

    // allocate container for incremental velocity field
    ierr = AllocateOnce(this->m_IncVelocityField, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // allocate work vec field 5 (1,2,3, and 4 are overwritten
    // during the computation of the incremental forward and adjoint
    // solve and the computation of the incremental body force)
    ierr = Assert(Hvtilde != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vtilde != NULL, "null pointer"); CHKERRQ(ierr);
    VecField *hv = nullptr;
    VecField *v = nullptr;
    ierr = AllocateOnce(hv, this->m_Opt, Hvtilde); CHKERRQ(ierr);
    ierr = AllocateOnce(v, this->m_Opt, vtilde); CHKERRQ(ierr);
    
    /*TwoLevelFinite op(this->m_Opt);
    op.Restrict(this->m_WorkVecField1, hv);
    op.Prolong(hv, this->m_WorkVecField1);*/
    
    // apply inverse of 2nd variation of regularization model to
    // incremental body force: (\beta\D{A})^{-1/2}\D{K}[\vect{\tilde{b}}](\beta\D{A})^{-1/2}
  
    // apply (\beta\D{A})^{-1/2} to incremental velocity field
    ierr = this->m_Regularization->ApplyInverse(this->m_IncVelocityField, v, true); CHKERRQ(ierr);
    
    
    /*{
      VecField *g_vec = this->m_IncVelocityField;
    
      TwoLevelFFT op(this->m_Opt);
      op.Restrict(this->m_WorkVecField1, g_vec);
      op.Prolong(g_vec, this->m_WorkVecField1);
    }*/

    // now solve the PDEs given the preconditioned incremental velocity field

    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t) and compute incremental body force
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // apply (\beta\D{A})^{-1/2} to incremental body force
    ierr = this->m_Regularization->ApplyInverse(hv, this->m_WorkVecField2, true); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \vect{\tilde{v}} + (\beta \D{A})^{-1/2}\D{K}[\vect{\tilde{b}}](\beta \D{A})^{-1/2}
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    hd = this->m_Opt->GetLebesgueMeasure();
    //ierr = hv->Scale(hd); CHKERRQ(ierr);
    ierr = hv->AXPY(hd, v); CHKERRQ(ierr);
    
    /*{
      VecField *g_vec = hv;
    
      TwoLevelFFT op(this->m_Opt);
      op.Restrict(this->m_WorkVecField1, g_vec);
      op.Prolong(g_vec, this->m_WorkVecField1);
    }*/

    //TwoLevelFFT op(this->m_Opt);
    /*op.Restrict(this->m_WorkVecField1, hv);
    op.Prolong(hv, this->m_WorkVecField1);*/
    
    ierr = Free(hv); CHKERRQ(ierr);
    ierr = Free(v); CHKERRQ(ierr);
    
    // pass to output
    /*if (Hvtilde != NULL) {
        ierr = this->m_WorkVecField5->GetComponents(Hvtilde); CHKERRQ(ierr);
    }*/
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

PetscErrorCode CLAIRE::SymTwoLevelHessMatVec(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType hd;
    PetscFunctionBegin;
    
    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("symmetric preconditioned two-level hessian matvec"); CHKERRQ(ierr);
    }

    this->m_Opt->Enter(__func__);

    ierr = Assert(Hvtilde != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vtilde != NULL, "null pointer"); CHKERRQ(ierr);
    VecField *hv = nullptr;
    VecField *v = nullptr;
    ierr = AllocateOnce(hv, this->m_Opt, Hvtilde); CHKERRQ(ierr);
    ierr = AllocateOnce(v, this->m_Opt, vtilde); CHKERRQ(ierr);

    // parse input (store incremental velocity field \tilde{v})
    ierr = AllocateOnce(this->m_CoarseReg->m_IncVelocityField, this->m_CoarseRegOpt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_CoarseReg->m_WorkVecField5, this->m_CoarseRegOpt); CHKERRQ(ierr);

    ScalarType beta;
    beta = this->m_Opt->m_RegNorm.beta[0];
    /*if (this->m_Opt->m_RegNorm.beta[3] > 0 && this->m_Opt->m_RegNorm.beta[3] > beta) {
          beta = this->m_Opt->m_RegNorm.beta[3];
    }*/
    
    TwoLevelRegFFT op(this->m_Opt, beta, 0, true);
    //TwoLevelFinite op(this->m_Opt);
    
    op.Restrict(this->m_CoarseReg->m_WorkVecField5, v);

    //this->m_CoarseRegOpt->m_RegNorm.beta[0] = beta;
    //this->m_CoarseReg->m_Regularization->ApplyInverse(this->m_CoarseReg->m_WorkVecField5,this->m_CoarseReg->m_WorkVecField5, true);
    
    if (beta <= 0.1) {
    
    VecField p(this->m_CoarseRegOpt);
    VecField r(this->m_CoarseRegOpt);
    VecField x(this->m_CoarseRegOpt);
    VecField m(this->m_CoarseRegOpt);
    VecField *b = this->m_CoarseReg->m_WorkVecField5;
    
    x.SetValue(0);
    
    ScalarType rTr, pTAp, rel;
    
    r.Copy(b);
    p.Copy(&r);
    
    VecDot(r.m_X, r.m_X, &rTr);
        
    rel=sqrt(rTr)*this->m_Opt->m_KrylovMethod.reltol*this->m_Opt->m_KrylovMethod.pctolint[0];
    
    this->m_CoarseRegOpt->m_KrylovMethod.matvectype = PRECONDMATVECSYM;
    
    for (int k=0;k<this->m_Opt->m_KrylovMethod.pcmaxit;++k) {
      this->m_CoarseReg->HessianMatVec(m.m_X, p.m_X, false);
      VecDot(p.m_X, m.m_X, &pTAp);
      ScalarType cgalpha = rTr/pTAp;
      x.AXPY(cgalpha, &p);
      r.AXPY(-cgalpha, &m);
      ScalarType rTr2;
      VecDot(r.m_X, r.m_X, &rTr2);
      if (sqrt(rTr2) < rel) {
        rTr = rTr2;
        break;
      }
      ScalarType cgbeta = rTr2/rTr;
      p.Scale(cgbeta);
      p.AXPY(1., &r);
      rTr = rTr2;
    }
    
    //this->m_CoarseReg->m_Regularization->ApplyInverse(&x,&x, true);
    
    op.Prolong(hv, &x);
    } else {
      op.Prolong(hv, this->m_CoarseReg->m_WorkVecField5);
    }
    this->m_Regularization->ApplyInverse(hv, hv, true);
    
   /* //this->m_CoarseReg->m_Regularization->ApplyInverse(this->m_CoarseReg->m_IncVelocityField, this->m_CoarseReg->m_WorkVecField5, true);
    
    this->m_CoarseReg->SolveIncStateEquation();
    this->m_CoarseReg->SolveIncAdjointEquation();
    
    this->m_CoarseReg->m_Regularization->ApplyInverse(this->m_CoarseReg->m_WorkVecField1, this->m_CoarseReg->m_WorkVecField2, true);
    
    //hd = this->m_CoarseRegOpt->GetLebesgueMeasure();
    //this->m_CoarseReg->m_WorkVecField1->Scale(this->m_Opt->GetLebesgueMeasure()/hd);
    
    hd = this->m_CoarseRegOpt->GetLebesgueMeasure();
    ierr = this->m_CoarseReg->m_WorkVecField1->AXPY(hd, this->m_CoarseReg->m_WorkVecField5); CHKERRQ(ierr);
    
    op.Prolong(hv, this->m_CoarseReg->m_WorkVecField1);  
    */
    
    ierr = Free(hv); CHKERRQ(ierr);
    ierr = Free(v); CHKERRQ(ierr);
    
    // pass to output
    /*if (Hvtilde != NULL) {
        ierr = this->m_WorkVecField5->GetComponents(Hvtilde); CHKERRQ(ierr);
    }*/
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute initial condition given some initial guess for
 * the state variable $m$ and the adjoint variable $\lambda$
 * @param[in] m state variable
 * @param[in] lambda adjoint variable
 *******************************************************************/
PetscErrorCode CLAIRE::ComputeInitialCondition(Vec m, Vec lambda) {
    PetscErrorCode ierr = 0;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

    ext = this->m_Opt->m_FileNames.extension;

    // allocate container for incremental velocity field
    ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);

    // allocate regularization model
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // allocate state and adjoint variables
    ierr = AllocateOnce(this->m_StateVariable, this->m_Opt, true, true); CHKERRQ(ierr);

    // allocate state and adjoint variables
    ierr = AllocateOnce(this->m_AdjointVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("piccard iteration"); CHKERRQ(ierr);
    }

    // copy input state and adjoint variable to class variables
    ierr = this->m_StateVariable->Copy(m); CHKERRQ(ierr);
    ierr = this->m_AdjointVariable->Copy(lambda); CHKERRQ(ierr);
    //ierr = VecCopy(m, *this->m_StateVariable); CHKERRQ(ierr);
    //ierr = VecCopy(lambda, *this->m_AdjointVariable); CHKERRQ(ierr);

    // compute body force (assigned to work vec field 2)
    // TODO: this will crash for GAUSS NEWTON
    //ierr = this->ComputeBodyForce(); CHKERRQ(ierr);

    // piccard step: solve A[v] = - ht \sum_j \lambda^j grad(m^j)
    ierr = this->m_WorkVecField2->Scale(-1.0); CHKERRQ(ierr);

    // apply inverse regularization operator / spectral preconditioning
    ierr = this->m_Regularization->ApplyInverse(this->m_VelocityField, this->m_WorkVecField2); CHKERRQ(ierr);

    // reset the adjoint variables
    //ierr = VecSet(*this->m_StateVariable, 0.0); CHKERRQ(ierr);
    //ierr = VecSet(*this->m_AdjointVariable, 0.0); CHKERRQ(ierr);
    ierr = this->m_StateVariable->Set(0.0); CHKERRQ(ierr);
    ierr = this->m_AdjointVariable->Set(0.0); CHKERRQ(ierr);

    ierr = this->m_ReadWrite->Write(this->m_VelocityField, "initial-condition"+ext); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::StoreStateVariable() {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt;
    ScalarType *p_m = NULL, *p_mj = NULL;
    std::stringstream ss;
    std::string ext;

    PetscFunctionBegin;
    
    ierr = DebugGPUNotImplemented(); CHKERRQ(ierr);

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    ierr = Assert(nt > 0, "nt <= 0"); CHKERRQ(ierr);
    ext = this->m_Opt->m_FileNames.extension;

    /*if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }*/
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
    ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

    // store time history
    ierr = VecGetArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
    // store individual time points
    for (IntType j = 0; j <= nt; ++j) {
        for (IntType k = 0; k < nc; ++k) {
            ierr = VecGetArray(*this->m_WorkScaField1, &p_mj); CHKERRQ(ierr);
            try {std::copy(p_m+j*nl*nc + k*nl, p_m+j*nl*nc + (k+1)*nl, p_mj);}
            catch (std::exception& err) {
                ierr = ThrowError(err); CHKERRQ(ierr);
            }
            ierr = VecRestoreArray(*this->m_WorkScaField1, &p_mj); CHKERRQ(ierr);
            // write out
            ss.str(std::string()); ss.clear();

            if (nc > 1) {
                ss << "state-variable-k=" << std::setw(3) << std::setfill('0')  << k
                   << "-j=" << std::setw(3) << std::setfill('0') << j << ext;
            } else {
                ss << "state-variable-j=" << std::setw(3) << std::setfill('0') << j << ext;
            }
            ierr = this->m_ReadWrite->Write(*this->m_WorkScaField1, ss.str()); CHKERRQ(ierr);
        }  // for number of vector components
    }  // for number of time points
    ierr = VecRestoreArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveStateEquation() {
    PetscErrorCode ierr = 0;
    IntType nc, nt;
    std::stringstream ss;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    ZeitGeist_define(PDE_STATE);
    ZeitGeist_tick(PDE_STATE);
    
    DebugGPUStartEvent(__FUNCTION__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    
    if (this->m_TransportProblem == nullptr) {
      ierr = this->SetupTransportProblem(); CHKERRQ(ierr);
      //ierr = this->m_TransportProblem->SetDifferentiation(Differentiation::Type::Spectral); CHKERRQ(ierr);
    }

    // check cfl condition / update time step
    if (this->m_Opt->m_PDESolver.monitorcflnumber ||
        this->m_Opt->m_PDESolver.adapttimestep) {
        ierr = this->ComputeCFLCondition(); CHKERRQ(ierr);
    }

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    ierr = Assert(nt > 0, "nt <= 0"); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "solving state equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->m_Domain.nx[0]
           << "," << this->m_Opt->m_Domain.nx[1]
           << "," << this->m_Opt->m_Domain.nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // set initial condition m_0 = m_T
    ierr = this->SetInitialState(this->m_TemplateImage->m_X); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 3) {
      if(this->m_StateVariable) this->m_StateVariable->DebugInfo("state       ",__LINE__,__FILE__);
      else DbgMsg3("state            : nullptr");
      if(this->m_TemplateImage) this->m_TemplateImage->DebugInfo("template    ",__LINE__,__FILE__);
      else DbgMsg3("template         : nullptr");
    }
    
    if (this->m_Opt->m_Mem.savestategrad) {
      if (!this->m_GradientState) {
        ierr = AllocateArrayOnce(this->m_GradientState, nc*(this->m_Opt->m_RegFlags.runinversion?nt+1:1)); CHKERRQ(ierr);
        for (int i=0;i<nc*(this->m_Opt->m_RegFlags.runinversion?nt+1:1);++i) {
          this->m_GradientState[i] = nullptr;
          ierr = AllocateOnce(this->m_GradientState[i], this->m_Opt); CHKERRQ(ierr);
        }
      }
    }
    if (this->m_Opt->m_Mem.savestategradx) {
      if (!this->m_GradientXState) {
        ierr = AllocateArrayOnce(this->m_GradientXState, nc*(this->m_Opt->m_RegFlags.runinversion?nt:1)); CHKERRQ(ierr);
        for (int i=0;i<nc*(this->m_Opt->m_RegFlags.runinversion?nt:1);++i) {
          this->m_GradientXState[i] = nullptr;
          ierr = AllocateOnce(this->m_GradientXState[i], this->m_Opt); CHKERRQ(ierr);
        }
      }
    }

    ierr = this->m_TransportProblem->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetControlVariable(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetGradientState(this->m_GradientState); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetGradientXState(this->m_GradientXState); CHKERRQ(ierr);

    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);
    
    ierr = this->m_TransportProblem->SolveForwardProblem(); CHKERRQ(ierr);
    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ScalarType maxval, minval, nvx1, nvx2, nvx3;
        ierr = VecMax(*this->m_StateVariable, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(*this->m_StateVariable, NULL, &minval); CHKERRQ(ierr);
        ss  << "state variable: [" << std::scientific
            << minval << "," << maxval << "]";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = this->m_VelocityField->Norm(nvx1, nvx2, nvx3); CHKERRQ(ierr);
        ss  << "velocity norm: (" << std::scientific
            << nvx1 << "," << nvx2 << "," << nvx3 <<")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
        
        ierr = VecMax(this->m_VelocityField->m_X1, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(this->m_VelocityField->m_X1, NULL, &minval); CHKERRQ(ierr);
        ss  << "velocity min/max: [" << std::scientific
            << minval << "," << maxval << "]";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
        ierr = VecMax(this->m_VelocityField->m_X2, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(this->m_VelocityField->m_X2, NULL, &minval); CHKERRQ(ierr);
        ss  << "                  [" << std::scientific
            << minval << "," << maxval << "] ";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
        ierr = VecMax(this->m_VelocityField->m_X3, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(this->m_VelocityField->m_X3, NULL, &minval); CHKERRQ(ierr);
        ss  << "                  [" << std::scientific
            << minval << "," << maxval << "]";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // store time series
    if (this->m_Opt->m_ReadWriteFlags.timeseries) {
        ierr = this->StoreStateVariable(); CHKERRQ(ierr);
    }
    
    if (this->m_Opt->m_Verbosity > 3) {
      if(this->m_StateVariable) this->m_StateVariable->DebugInfo("state       ",__LINE__,__FILE__);
      else DbgMsg3("state            : nullptr");
    }

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);
    
    DebugGPUStopEvent();
    
    ZeitGeist_tock(PDE_STATE);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveAdjointEquation() {
    PetscErrorCode ierr = 0;
    IntType nc, nt;
    ScalarType hd;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    ZeitGeist_define(PDE_ADJOINT);
    ZeitGeist_tick(PDE_ADJOINT);
    
    DebugGPUStartEvent(__FUNCTION__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);
    nc = this->m_Opt->m_Domain.nc;
    hd  = this->m_Opt->GetLebesgueMeasure();

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "solving adjoint equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->m_Domain.nx[0]
           << "," << this->m_Opt->m_Domain.nx[1]
           << "," << this->m_Opt->m_Domain.nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ierr = AllocateOnce(this->m_AdjointVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);
    
    if (this->m_DistanceMeasure == NULL) {
        ierr = this->SetupDistanceMeasure(); CHKERRQ(ierr);
    }
    ierr = this->m_DistanceMeasure->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetAdjointVariable(this->m_AdjointVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetFinalConditionAE(); CHKERRQ(ierr);
    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 3) {
      if(this->m_StateVariable) this->m_StateVariable->DebugInfo("state       ",__LINE__,__FILE__);
      else DbgMsg3("state            : nullptr");
      if(this->m_AdjointVariable) this->m_AdjointVariable->DebugInfo("adjoint     ",__LINE__,__FILE__);
      else DbgMsg3("adjoint          : nullptr");
      if(this->m_ReferenceImage) this->m_ReferenceImage->DebugInfo("reference   ",__LINE__,__FILE__);
      else DbgMsg3("reference        : nullptr");
      if(this->m_TemplateImage) this->m_TemplateImage->DebugInfo("template    ",__LINE__,__FILE__);
      else DbgMsg3("template         : nullptr");
    }
    
    if (this->m_TransportProblem == nullptr) {
      ierr = this->SetupTransportProblem(); CHKERRQ(ierr);
    }
    ierr = this->m_TransportProblem->SetTemplateImage(this->m_TemplateImage); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetControlVariable(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetAdjointVariable(this->m_AdjointVariable); CHKERRQ(ierr);

    ierr = this->m_TransportProblem->SolveAdjointProblem(); CHKERRQ(ierr);
    
    ierr = this->m_WorkVecField2->DebugInfo("post adj grad", __LINE__, __FILE__); CHKERRQ(ierr);
        
    // apply projection
    ierr = this->ApplyProjection(); CHKERRQ(ierr);
    
    // scale result by hd
    ierr = this->m_WorkVecField2->Scale(hd); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 3) {
      if(this->m_AdjointVariable) this->m_AdjointVariable->DebugInfo("adjoint     ",__LINE__,__FILE__);
      else DbgMsg3("adjoint          : nullptr");
    }

    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);
        
    DebugGPUStopEvent();
    ZeitGeist_tock(PDE_ADJOINT);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \idiv m\vect{v} = 0  with initial condition m_0 = m_T
 * (solved forward in time)
 *******************************************************************/
PetscErrorCode CLAIRE::SolveContinuityEquationSL() {
    PetscErrorCode ierr = 0;
    ScalarType mx, rhs0, rhs1, ht, hthalf;
    ScalarType *p_m = nullptr;
    ScalarType *p_divv = nullptr;
    ScalarType *p_divvx = nullptr;
    ScalarType *p_mx = nullptr;
    IntType nl, nc, nt, l, lnext;
    bool store;

    PetscFunctionBegin;

    ierr = DebugGPUNotImplemented(); CHKERRQ(ierr);

    this->m_Opt->Enter(__func__);

    // flag to identify if we store the time history
    store = this->m_Opt->m_RegFlags.runinversion;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField3, &p_mx); CHKERRQ(ierr);

    // compute divergence of velocity field
    this->m_Differentiation->Divergence(p_divv, this->m_VelocityField);

    // interpolate velocity field v(X)
//    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField1, this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    // compute divergence of velocity field at X
//    ierr = this->m_WorkVecField1->GetArrays(p_vec1, p_vec2, p_vec3); CHKERRQ(ierr);
//    this->m_Opt->StartTimer(FFTSELFEXEC);
//    accfft_divergence_t(p_divvx, p_vec1, p_vec2, p_vec3, this->m_Opt->m_FFT.plan, timer);
//    this->m_Opt->StopTimer(FFTSELFEXEC);
//    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    // evaluate div(v) along characteristic X
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_divvx, p_divv, "state"); CHKERRQ(ierr);

    // perform numerical time integration for state variable
    for (IntType j = 0; j < nt; ++j) {
        if (store) {
            l = j*nl*nc; lnext = (j+1)*nl*nc;
        } else {
            l = 0; lnext = 0;
        }

        // scaling for trapezoidal rule (for body force)
        for (IntType k = 0; k < nc; ++k) {
            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_mx, p_m + l + k*nl, "state"); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                mx = p_mx[i];

                rhs0 = -mx*p_divvx[i];
                rhs1 = -(mx + ht*rhs0)*p_divv[i];
                //if (std::abs(p_divv[i]) > 0.1) { std::cout << p_divv[i] << " ";}
                // compute \lambda(x,t^{j+1})
                p_m[lnext + k*nl + i] = mx + hthalf*(rhs0 + rhs1);
            }
        }
}  // omp
    }

    ierr = RestoreRawPointer(this->m_WorkScaField3, &p_mx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the incremental state equation
 * \p_t \tilde{m} + \igrad m \cdot \vect{\tilde{v}}
 *                + \igrad \tilde{m} \cdot \vect{v} = 0
 * subject to \tilde{m}_0 = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncStateEquation(void) {
    PetscErrorCode ierr = 0;
    IntType nc, nt;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    ZeitGeist_define(PDE_INCSTATE);
    ZeitGeist_tick(PDE_INCSTATE);
    
    DebugGPUStartEvent(__FUNCTION__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "solving incremental state equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->m_Domain.nx[0]
           << "," << this->m_Opt->m_Domain.nx[1]
           << "," << this->m_Opt->m_Domain.nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // allocate variables
    ierr = AllocateOnce(this->m_IncStateVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);
    
    
    if (this->m_TransportProblem == nullptr) {
      ierr = this->SetupTransportProblem(); CHKERRQ(ierr);
      ierr = this->m_TransportProblem->SetControlVariable(this->m_VelocityField); CHKERRQ(ierr);
    }
    ierr = this->m_TransportProblem->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetIncStateVariable(this->m_IncStateVariable); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetIncControlVariable(this->m_IncVelocityField); CHKERRQ(ierr);

    // set initial value
    ierr = this->m_IncStateVariable->Set(0.0); CHKERRQ(ierr);

    // start timer
    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    ierr = this->m_TransportProblem->SolveIncForwardProblem(); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 2) {
        ScalarType maxval, minval;
        ierr = VecMax(*this->m_IncStateVariable, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(*this->m_IncStateVariable, NULL, &minval); CHKERRQ(ierr);
        ss << "incremental state variable: [" << std::scientific << minval << "," << maxval << "]";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ScalarType mtilde_norm;
    ierr = VecNorm(this->m_IncStateVariable->m_X, NORM_2, &mtilde_norm);
    //PetscPrintf(PETSC_COMM_WORLD, "mtilde norm = %0.4e\n", mtilde_norm);


    // stop timer
    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    DebugGPUStopEvent();
    
    ZeitGeist_tock(PDE_INCSTATE);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncAdjointEquation(void) {
    PetscErrorCode ierr = 0;
    IntType nc, nt;
    ScalarType hd;
    std::stringstream ss;

    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ZeitGeist_define(PDE_INCADJOINT);
    ZeitGeist_tick(PDE_INCADJOINT);
    
    DebugGPUStartEvent(__FUNCTION__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    hd  = this->m_Opt->GetLebesgueMeasure();
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "solving incremental adjoint equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->m_Domain.nx[0]
           << "," << this->m_Opt->m_Domain.nx[1]
           << "," << this->m_Opt->m_Domain.nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // allocate state and adjoint variables
    ierr = AllocateOnce(this->m_IncAdjointVariable, this->m_Opt, true, this->m_Opt->m_OptPara.method == FULLNEWTON); CHKERRQ(ierr);
    
    if (this->m_DistanceMeasure == NULL) {
        ierr = this->SetupDistanceMeasure(); CHKERRQ(ierr);
    }
    //ierr = this->m_DistanceMeasure->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
    //ierr = this->m_DistanceMeasure->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetIncStateVariable(this->m_IncStateVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetIncAdjointVariable(this->m_IncAdjointVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetFinalConditionIAE(); CHKERRQ(ierr);
    
    if (this->m_TransportProblem == nullptr) {
      ierr = this->SetupTransportProblem(); CHKERRQ(ierr);
    }
    ierr = this->m_TransportProblem->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetControlVariable(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetIncAdjointVariable(this->m_IncAdjointVariable); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SolveIncAdjointProblem(); CHKERRQ(ierr);
    
    // apply K[\tilde{b}]
    ierr = this->ApplyProjection(); CHKERRQ(ierr);
    // scale result by hd
    ierr = this->m_WorkVecField2->Scale(hd); CHKERRQ(ierr);
    if (this->m_Opt->m_Verbosity > 2) {
        ScalarType maxval, minval;
        ierr = VecMax(*this->m_IncAdjointVariable, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(*this->m_IncAdjointVariable, NULL, &minval); CHKERRQ(ierr);
        ss << "incremental adjoint variable: [" << std::scientific << minval << "," << maxval << "]";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    DebugGPUStopEvent();
    
    ZeitGeist_tock(PDE_INCADJOINT);
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief apply projection to map \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
PetscErrorCode CLAIRE::ApplyProjection() {
    PetscErrorCode ierr = 0;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief finalize the current iteration
 *******************************************************************/
PetscErrorCode CLAIRE::FinalizeIteration(Vec v) {
    PetscErrorCode ierr = 0;
    int rank;
    IntType nl, nc, nt, iter;
    std::string filename, ext;
    std::stringstream ss;
    std::ofstream logwriter;
    ScalarType *p_m1 = NULL, *p_m = NULL, rval, dval, jval, gval;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get number of time points and grid points
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    // parse extension
    ext = this->m_Opt->m_FileNames.extension;

    // if not yet allocted, do so
    ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
    // set velocity field
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // log objective values
    if (this->m_Opt->m_Log.enabled[LOGCONV]) {
        iter = this->m_Opt->GetCounter(ITERATIONS);
        jval = this->m_Opt->m_Monitor.jval;
        dval = this->m_Opt->m_Monitor.dval;
        rval = this->m_Opt->m_Monitor.rval;
        this->m_Opt->LogConvergence(iter, jval, dval, rval);
    }

    // log norm of gradient
    if (this->m_Opt->m_Log.enabled[LOGGRAD]) {
        gval = this->m_Opt->m_Monitor.gradnorm;
        this->m_Opt->m_Log.gradnorm.push_back(gval);
    }

    // store iterates
    if (this->m_Opt->m_ReadWriteFlags.iterates) {
        // allocate
        ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
        ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
        if (this->m_WorkScaFieldMC == NULL) {
            ierr = AllocateOnce(this->m_WorkScaFieldMC, this->m_Opt); CHKERRQ(ierr);
            //ierr = VecCreate(this->m_WorkScaFieldMC, nl*nc, ng*nc); CHKERRQ(ierr);
        }

        iter = this->m_Opt->GetCounter(ITERATIONS);
        ierr = Assert(iter >= 0, "problem in counter"); CHKERRQ(ierr);
        
        std::string ver = "cpu-";
#ifdef REG_HAS_CUDA
        ver = "gpu-";
#endif

        // copy memory for m_1
        ierr = VecGetArray(*this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecGetArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception& err) {
            ierr = ThrowError(err); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(*this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecRestoreArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);

        ss  << ver << "deformed-template-image-i=" << std::setw(3) << std::setfill('0') << iter << ext;
        ierr = this->m_ReadWrite->Write(*this->m_WorkScaFieldMC, ss.str(), nc > 1); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
        
        // copy memory for m_1
        ierr = VecGetArray(*this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecGetArray(*this->m_AdjointVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m, p_m+nl*nc, p_m1);}
        catch (std::exception& err) {
            ierr = ThrowError(err); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(*this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecRestoreArray(*this->m_AdjointVariable, &p_m); CHKERRQ(ierr);

        ss  << ver << "adjoint-image-i=" << std::setw(3) << std::setfill('0') << iter << ext;
        ierr = this->m_ReadWrite->Write(*this->m_WorkScaFieldMC, ss.str(), nc > 1); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        // construct file names for velocity field components
        ss  << ver << "velocity-field-i=" << std::setw(3) << std::setfill('0') << iter << ext;
        filename = ss.str();
        ss.str(std::string()); ss.clear();

        // velocity field out
        ierr = this->m_ReadWrite->Write(this->m_VelocityField, filename); CHKERRQ(ierr);
    }  // store iterates

    if (this->m_Opt->m_StoreCheckPoints) {
        ierr = this->m_ReadWrite->Write(this->m_VelocityField, "velocity-field-checkpoint"+ext); CHKERRQ(ierr);
    }

    // compute determinant of deformation gradient and write it to file
    if (this->m_Opt->m_Monitor.detdgradenabled) {
        ierr = this->ComputeDetDefGrad(); CHKERRQ(ierr);

        // if user enabled the logger
        if (this->m_Opt->m_Log.enabled[LOGJAC]) {
            if (rank == 0) {
                filename  = this->m_Opt->m_FileNames.xfolder;
                filename += "registration-performance-detdefgrad.log";

                // create output file or append to output file
                logwriter.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
                ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);
                ss  << std::scientific
                    <<  "iter = "     << this->m_Opt->GetCounter(ITERATIONS)
                    <<  "   betav = " << this->m_Opt->m_RegNorm.beta[0] << "    "
                    << std::left << std::setw(20) << this->m_Opt->m_Monitor.detdgradmin << " "
                                 << std::setw(20) << this->m_Opt->m_Monitor.detdgradmean <<" "
                                 << std::setw(20) << this->m_Opt->m_Monitor.detdgradmax;
                logwriter << ss.str() << std::endl;
                ss.str(std::string()); ss.clear();
            }  // if on master rank
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief finalize the registration
 *******************************************************************/
PetscErrorCode CLAIRE::Finalize(VecField* v) {
    PetscErrorCode ierr = 0;
    std::string filename, fn, ext;
    IntType nl, nc, nt;
    int rank, nproc;
    std::ofstream logwriter;
    std::stringstream ss, ssnum;
    ScalarType value;
    ScalarType *p_m1 = NULL, *p_mt = NULL, *p_mr = NULL, *p_m = NULL, *p_dr = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

    // get sizes
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    if (this->m_Opt->m_Verbosity > 1) {
        ierr = DbgMsg1("finalizing registration"); CHKERRQ(ierr);
    }

    // if not yet allocted, do so and copy input
    ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
    ierr = this->m_VelocityField->Copy(v); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);

    // process timer
    ierr = this->m_Opt->ProcessTimers(); CHKERRQ(ierr);

    // parse extension
    ext = this->m_Opt->m_FileNames.extension;

    // compute residuals
    if (this->m_Opt->m_Log.enabled[LOGDIST]) {
        ierr = AllocateOnce(this->m_WorkScaFieldMC, this->m_Opt, true); CHKERRQ(ierr);
        ierr = VecWAXPY(*this->m_WorkScaFieldMC, -1.0, *this->m_TemplateImage, *this->m_ReferenceImage); CHKERRQ(ierr);

        ierr = VecNorm(*this->m_WorkScaFieldMC, NORM_2, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(0, value);

        ierr = VecNorm(*this->m_WorkScaFieldMC, NORM_INFINITY, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(1, value);

        // deformed template out (compute solution of state equation)
        ierr = this->SolveStateEquation(); CHKERRQ(ierr);

        // copy memory for m_1
        ierr = VecGetArray(*this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecGetArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception& err) {
            ierr = ThrowError(err); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
        ierr = VecRestoreArray(*this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);

        // ||m_R - m_1||
        ierr = VecAXPY(*this->m_WorkScaFieldMC, -1.0, *this->m_ReferenceImage); CHKERRQ(ierr);

        ierr = VecNorm(*this->m_WorkScaFieldMC, NORM_2, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(2, value);

        ierr = VecNorm(*this->m_WorkScaFieldMC, NORM_INFINITY, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(3, value);
    }

    // write deformed template image to file
    if (this->m_Opt->m_ReadWriteFlags.deftemplate) {
        ierr = AllocateOnce(this->m_WorkScaFieldMC, this->m_Opt, true); CHKERRQ(ierr);
        // copy memory for m_1
        ierr = VecGetArray(*this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecGetArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception& err) {
            ierr = ThrowError(err); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
        ierr = VecRestoreArray(*this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);

        if (this->m_Opt->m_RegFlags.registerprobmaps) {
            ierr = EnsurePartitionOfUnity(*this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = ShowValues(*this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = ComputeBackGround(*this->m_WorkScaField1, *this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->WriteT(*this->m_WorkScaField1, "background-image" + ext, false); CHKERRQ(ierr);
        }

        ierr = this->m_ReadWrite->WriteT(*this->m_WorkScaFieldMC, "deformed-template-image" + ext, nc > 1); CHKERRQ(ierr);
    }

    // write residual images to file
    if (this->m_Opt->m_ReadWriteFlags.residual || this->m_Opt->m_ReadWriteFlags.invresidual) {
        ierr = AllocateOnce(this->m_WorkScaFieldMC, this->m_Opt, true); CHKERRQ(ierr);
        ierr = VecGetArray(*this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);

        // write residual at t = 0 to file
        ierr = VecGetArray(*this->m_TemplateImage, &p_mt); CHKERRQ(ierr);
        if (this->m_Opt->m_ReadWriteFlags.residual) {
            ierr = VecGetArray(*this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            for (IntType i = 0; i < nl*nc; ++i) {
                p_dr[i] = PetscAbs(p_mt[i] - p_mr[i]);
            }
            ierr = VecRestoreArray(*this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->Write(*this->m_WorkScaFieldMC, "residual-t=0" + ext, nc > 1); CHKERRQ(ierr);
        }

        if (this->m_Opt->m_ReadWriteFlags.invresidual) {
            ierr = VecGetArray(*this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            for (IntType i = 0; i < nl*nc; ++i) {
                p_dr[i] = 1.0 - PetscAbs(p_mt[i] - p_mr[i]);
            }
            ierr = VecRestoreArray(*this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->Write(*this->m_WorkScaFieldMC, "inv-residual-t=0" + ext, nc > 1); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(*this->m_TemplateImage, &p_mt); CHKERRQ(ierr);

        // write residual at t = 1 to file
        ierr = VecGetArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
        if (this->m_Opt->m_ReadWriteFlags.residual) {
            ierr = VecGetArray(*this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            for (IntType i = 0; i < nl*nc; ++i) {
                p_dr[i] = PetscAbs(p_m[nt*nl*nc + i] - p_mr[i]);
            }
            ierr = VecRestoreArray(*this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
//            ierr = ShowValues(this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->Write(*this->m_WorkScaFieldMC, "residual-t=1" + ext, nc > 1); CHKERRQ(ierr);
        }
        if (this->m_Opt->m_ReadWriteFlags.invresidual) {
            ierr = VecGetArray(*this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            for (IntType i = 0; i < nl*nc; ++i) {
                p_dr[i] = 1.0 - PetscAbs(p_m[nt*nl*nc + i] - p_mr[i]);
            }
            ierr = VecRestoreArray(*this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);

//            ierr = ShowValues(this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->Write(*this->m_WorkScaFieldMC, "inv-residual-t=1" + ext, nc > 1); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);

        // restore reference image
        ierr = VecRestoreArray(*this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    }

    // write velocity field to file
    if (this->m_Opt->m_ReadWriteFlags.velocity) {
        ierr = this->m_ReadWrite->Write(this->m_VelocityField, "velocity-field" + ext); CHKERRQ(ierr);
    }

    // write norm of velocity field to file
    if (this->m_Opt->m_ReadWriteFlags.velnorm) {
        ierr = this->m_VelocityField->Norm(*this->m_WorkScaField1); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(*this->m_WorkScaField1, "velocity-field-norm" + ext); CHKERRQ(ierr);
    }

    // write determinant of deformation gradient to file
    if (this->m_Opt->m_ReadWriteFlags.detdefgrad) {
        ierr = this->ComputeDetDefGrad(true); CHKERRQ(ierr);
    }

    // write determinant of deformation gradient to file
    if (this->m_Opt->m_ReadWriteFlags.defgrad) {
        if (this->m_DeformationFields == NULL) {
            ierr = this->SetupDeformationField(); CHKERRQ(ierr);
        }
        ierr = this->m_DeformationFields->ComputeDefGrad(true); CHKERRQ(ierr);
    }

    // write deformation map to file
    if (this->m_Opt->m_ReadWriteFlags.defmap) {
        if (this->m_DeformationFields == NULL) {
            ierr = this->SetupDeformationField(); CHKERRQ(ierr);
        }
        ierr = this->m_DeformationFields->ComputeDeformationMap(true); CHKERRQ(ierr);
    }

    // write deformation field to file
    if (this->m_Opt->m_ReadWriteFlags.deffield) {
        if (this->m_DeformationFields == NULL) {
            ierr = this->SetupDeformationField(); CHKERRQ(ierr);
        }
        ierr = this->m_DeformationFields->ComputeDisplacementField(true); CHKERRQ(ierr);
    }

    // write template and reference image
    if (this->m_Opt->m_ReadWriteFlags.templateim) {
//        ierr = ShowValues(this->m_TemplateImage, nc); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->WriteT(*this->m_TemplateImage, "template-image"+ext, nc > 1); CHKERRQ(ierr);
    }
    if (this->m_Opt->m_ReadWriteFlags.referenceim) {
//        ierr = ShowValues(this->m_ReferenceImage, nc); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->WriteR(*this->m_ReferenceImage, "reference-image"+ext, nc > 1); CHKERRQ(ierr);
    }

    // write log file
    ierr = this->m_Opt->WriteLogFile(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

}  // namespace reg

#endif  // _CLAIRE_CPP_
