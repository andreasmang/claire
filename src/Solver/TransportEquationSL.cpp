/*************************************************************************
 *  Copyright (c) 2018.
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

#ifndef _TRANSPORTEQUATIONSL_CPP_
#define _TRANSPORTEQUATIONSL_CPP_

#include "TransportEquationSL.hpp"

namespace reg {

/********************************************************************
 * @brief default constructor
 *******************************************************************/
TransportEquationSL::TransportEquationSL() : SuperClass() {
    this->Initialize();
}

/********************************************************************
 * @brief default destructor
 *******************************************************************/
TransportEquationSL::~TransportEquationSL() {
    this->ClearMemory();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
TransportEquationSL::TransportEquationSL(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode TransportEquationSL::Initialize() {
    PetscFunctionBegin;
    
    this->m_SemiLagrangianMethod = nullptr;
    this->m_GradientState = nullptr;

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode TransportEquationSL::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
        
    Free(this->m_SemiLagrangianMethod);

    PetscFunctionReturn(ierr);
}

PetscErrorCode TransportEquationSL::InitializeControlVariable(VecField *field) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_VelocityField = field;
    
    ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
  
    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the forward problem (i.e., the continuity equation)
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveForwardProblem() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
    if (VelocityIsZero) {
        ierr = SuperClass::SolveForwardProblem(); CHKERRQ(ierr);
    } else {
        ierr = this->SolveStateEquation(this->m_VelocityField); CHKERRQ(ierr);
    }
    
    ierr = this->ComputeGradientState(VelocityIsZero); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

PetscErrorCode TransportEquationSL::ComputeGradientState (bool V0) {
  PetscErrorCode ierr = 0;
  IntType nt, nc;
  const ScalarType *pM = nullptr;
  
  PetscFunctionBegin;
  
  this->m_Opt->Enter(__func__);
  
  nt = this->m_Opt->m_Domain.nt;
  nc = this->m_Opt->m_Domain.nc;
    
  ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
  ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    
  if (V0) {
    if (this->m_GradientState) {
      for (IntType k = 0; k < nc; ++k) {
        ierr = this->m_StateVariable->GetArrayRead(pM, k); CHKERRQ(ierr);
        ierr = this->m_Differentiation->Gradient(this->m_GradientState[k], pM); CHKERRQ(ierr);
        if (this->m_GradientXState) {
          ierr = this->m_GradientXState[k]->Copy(this->m_GradientState[k]); CHKERRQ(ierr);
        }
      }
      ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
      for (IntType j = 1; j < (this->m_Opt->m_RegFlags.runinversion?nt+1:1); ++j) {
        for (IntType k = 0; k < nc; ++k) {
          this->m_GradientState[j*nc+k]->Copy(this->m_GradientState[k]);
          if (this->m_GradientXState && j < nt) {
            ierr = this->m_GradientXState[j*nc + k]->Copy(this->m_GradientState[k]); CHKERRQ(ierr);
          }
        }
      }
    } else if (this->m_GradientXState) {
      for (IntType k = 0; k < nc; ++k) {
        ierr = this->m_StateVariable->GetArrayRead(pM, k); CHKERRQ(ierr);
        ierr = this->m_Differentiation->Gradient(this->m_WorkVecField[0], pM);
        for (IntType j = 0; j < (this->m_Opt->m_RegFlags.runinversion?nt:1); ++j) {
          ierr = this->m_GradientXState[k]->Copy(this->m_WorkVecField[0]); CHKERRQ(ierr);
        }
      }
      ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    }
  } else {
    ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
        
    if (this->m_GradientState) {
      for (IntType j = 0; j < (this->m_Opt->m_RegFlags.runinversion?nt+1:1); ++j) {
        for (IntType k = 0; k < nc; ++k) {
          ierr = this->m_StateVariable->GetArrayRead(pM, k, j); CHKERRQ(ierr);
          ierr = this->m_Differentiation->Gradient(this->m_GradientState[j*nc + k], pM);
        }
      }
      ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    }
    
    if (this->m_GradientXState) {
      // compute trajectory
      //ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
      //ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
      
      for (IntType j = 0; j < (this->m_Opt->m_RegFlags.runinversion?nt:1); ++j) {
        for (IntType k = 0; k < nc; ++k) {
          if (this->m_GradientState) {
            ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_GradientXState[j*nc + k], this->m_GradientState[j*nc + k], "state"); CHKERRQ(ierr);
          } else {
            ierr = this->m_StateVariable->GetArrayRead(pM, k, j); CHKERRQ(ierr);
            ierr = this->m_Differentiation->Gradient(this->m_WorkVecField[0], pM);
            ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_GradientXState[j*nc + k], this->m_WorkVecField[0], "state"); CHKERRQ(ierr);
          }
        }
      }
      if (!this->m_GradientState) {
        ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
      }
    }
  }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the inverse forward problem (i.e., the continuity equation)
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveInverseProblem() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
    if (VelocityIsZero) {
        ierr = SuperClass::SolveInverseProblem(); CHKERRQ(ierr);
    } else {
        ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
  
        // compute trajectory
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
        
        ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField[0], this->m_VelocityField, "state"); CHKERRQ(ierr);
        ierr = this->m_WorkVecField[1]->Copy(this->m_WorkVecField[0]); CHKERRQ(ierr);
        ierr = this->m_WorkVecField[1]->Scale(-1.); CHKERRQ(ierr);
        
        ierr = this->SolveStateEquation(this->m_WorkVecField[1]); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the forward problem (i.e., the continuity equation)
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveStateEquation(VecField *v) {
    PetscErrorCode ierr = 0;
    IntType nc, nt, l, lnext;
    ScalarType *pM = nullptr, *pMnext = nullptr;
    bool store = true;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(v != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);

    // flag to identify if we store the time history
    store = this->m_Opt->m_RegFlags.runinversion;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;

    ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(v, "state"); CHKERRQ(ierr);
    
    // get state variable m
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        if (store) {
            l = j; lnext = (j+1);
        } else {
            l = 0; lnext = 0;
        }
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            ierr = this->m_StateVariable->GetArrayReadWrite(pM, k, l); CHKERRQ(ierr);
            ierr = this->m_StateVariable->GetArrayReadWrite(pMnext, k, lnext); CHKERRQ(ierr);
            // compute m(X,t^{j+1}) (interpolate state variable)
            ierr = this->m_SemiLagrangianMethod->Interpolate(pMnext, pM, "state"); CHKERRQ(ierr);
        }
    }
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveAdjointProblem() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
    if (VelocityIsZero) {
        ierr = SuperClass::SolveAdjointProblem(); CHKERRQ(ierr);
    } else {
        ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveAdjointEquation() {
    PetscErrorCode ierr = 0;
    const ScalarType *pM = nullptr;
    IntType nc, nt, ll, llnext;
    TransportKernelAdjointSL kernel;
    bool fullnewton = false;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    kernel.scale = kernel.ht/static_cast<ScalarType>(nc);

    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[2] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);
  
    ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);

    // compute trajectory for adjoint equations
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    // for full newton we store the adjoint variable
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        fullnewton = true;
    }
    
    ierr = this->m_WorkScaField[0]->GetArrayWrite(kernel.pDivV); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->GetArrayWrite(kernel.pDivVx); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[2]->GetArrayWrite(kernel.pLx); CHKERRQ(ierr);

    // compute divergence of velocity field
    ierr = this->m_Differentiation->Divergence(kernel.pDivV,  this->m_VelocityField); CHKERRQ(ierr);
        
    // compute divergence of velocity field at X
    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pGm); CHKERRQ(ierr);
    }

    // evaluate div(v) along characteristic X
    ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pDivVx, kernel.pDivV, "adjoint"); CHKERRQ(ierr);
    
    // init body force for numerical integration
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArraysReadWrite(kernel.pB); CHKERRQ(ierr);

    // perform numerical time integration for adjoint variable and
    // add up body force
    for (IntType j = 0; j < nt; ++j) {
        if (fullnewton) {
            ll = (nt-j); llnext = (nt-(j+1));
        } else {
            ll = 0; llnext = 0;
        }

        // scaling for trapezoidal rule (for body force)
        if (j == 0) kernel.scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {
            ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL, k, ll); CHKERRQ(ierr);
            ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pLnext, k, llnext); CHKERRQ(ierr);
            ierr = this->m_StateVariable->GetArrayRead(pM, k, nt-j); CHKERRQ(ierr);

            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pLx, kernel.pL, "adjoint"); CHKERRQ(ierr);

            // compute gradient of m (for body force)
            if (this->m_GradientState) {
              ierr = this->m_GradientState[(nt-j)*nc + k]->GetArraysReadWrite(kernel.pGm); CHKERRQ(ierr);
            } else {
              ierr = this->m_Differentiation->Gradient(kernel.pGm, pM); CHKERRQ(ierr);
            }
            
            // compute \lambda(x,t^{j+1}) and bodyforce
            ierr = kernel.ComputeBodyForcePart1(); CHKERRQ(ierr);
            
            if (this->m_GradientState) {
              ierr = this->m_GradientState[(nt-j)*nc + k]->RestoreArrays(); CHKERRQ(ierr);
            }
        }
        // trapezoidal rule (revert scaling; for body force)
        if (j == 0) kernel.scale *= 2.0;
    }
    
    kernel.scale *= 0.5;
    
    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
      ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL, k); CHKERRQ(ierr);
      ierr = this->m_StateVariable->GetArrayRead(pM, k); CHKERRQ(ierr);

      // compute gradient of m (for incremental body force)
      if (this->m_GradientState) {
        ierr = this->m_GradientState[k]->GetArraysReadWrite(kernel.pGm); CHKERRQ(ierr);
      } else {
        ierr = this->m_Differentiation->Gradient(kernel.pGm, pM); CHKERRQ(ierr);
      }
      
      ierr = kernel.ComputeBodyForcePart2(); CHKERRQ(ierr);
      
      if (this->m_GradientState) {
        ierr = this->m_GradientState[k]->RestoreArrays(); CHKERRQ(ierr);
      }
    }
        
    ierr = this->m_WorkVecField[1]->RestoreArrays(); CHKERRQ(ierr);
    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[0]->RestoreArrays(); CHKERRQ(ierr);
    }
    
    ierr = this->m_WorkScaField[0]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[2]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_AdjointVariable->RestoreArray(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveIncForwardProblem() {
    PetscErrorCode ierr = 0;
    IntType nt, nc, lmt, lmtnext;
    ScalarType *pMtilde = nullptr;
    const ScalarType *pM = nullptr;
    TransportKernelIncStateSL kernel;
    bool fullnewton = false;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.hthalf = 0.5*this->m_Opt->GetTimeStepSize();
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[2] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);

    if (this->m_SemiLagrangianMethod == nullptr) {
        ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    } else if (this->m_Opt->m_KrylovMethod.pctype == TWOLEVEL) {
      ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
      ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    }
    //ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);

    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {   // gauss newton
        fullnewton = true;
    }
    
    ierr = this->m_VelocityField->DebugInfo("velocity", __LINE__, __FILE__); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->DebugInfo("inc velocity", __LINE__, __FILE__); CHKERRQ(ierr);
    ierr = this->m_IncStateVariable->DebugInfo("inc state pre ", __LINE__, __FILE__); CHKERRQ(ierr);
    ierr = this->m_StateVariable->DebugInfo("state pre ", __LINE__, __FILE__); CHKERRQ(ierr);

    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[2]->GetArraysWrite(kernel.pGm); CHKERRQ(ierr);
    }
    if (!this->m_GradientXState) {
      ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pGmx); CHKERRQ(ierr);
    }

    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField[1], this->m_IncVelocityField, "state"); CHKERRQ(ierr);
    
    ierr = this->m_WorkVecField[1]->DebugInfo("work vec", __LINE__, __FILE__); CHKERRQ(ierr);

    ierr = this->m_WorkVecField[1]->GetArraysRead(kernel.pVtildex); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArraysRead(kernel.pVtilde); CHKERRQ(ierr);

    for (IntType k = 0; k < nc; ++k) {  // for all image components
      if (this->m_GradientState) {
        ierr = this->m_GradientState[k]->GetArraysReadWrite(kernel.pGm); CHKERRQ(ierr);
      } else {
        ierr = this->m_StateVariable->GetArrayRead(pM, k, 0, 0); CHKERRQ(ierr);
        ierr = this->m_Differentiation->Gradient(kernel.pGm, pM); CHKERRQ(ierr);
      }
      for (IntType j = 0; j < nt; ++j) {  // for all time points
        if (fullnewton) {   // full newton
            lmt = j; lmtnext = (j+1);
        } else {
            lmt = 0; lmtnext = 0;
        }
          ierr = this->m_IncStateVariable->GetArrayReadWrite(kernel.pMtilde, k, lmtnext, 0); CHKERRQ(ierr);
          ierr = this->m_IncStateVariable->GetArrayReadWrite(pMtilde, k, lmt, 0); CHKERRQ(ierr);
          //ierr = this->m_StateVariable->GetArrayRead(pM, k, j); CHKERRQ(ierr);

          // interpolate incremental adjoint variable \tilde{m}^j(X)
          ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pMtilde, pMtilde, "state"); CHKERRQ(ierr);
          // compute gradient for state variable
          //ierr = this->m_Differentiation->Gradient(kernel.pGm, pM); CHKERRQ(ierr);

          if (this->m_GradientXState) {
            ierr = this->m_GradientXState[j*nc + k]->GetArraysReadWrite(kernel.pGmx); CHKERRQ(ierr);
          } else {
            ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pGmx[0], kernel.pGmx[1], kernel.pGmx[2], kernel.pGm[0], kernel.pGm[1], kernel.pGm[2], "state"); CHKERRQ(ierr);
          }
          
          // first part of time integration
          //ierr = kernel.TimeIntegrationPart1(); CHKERRQ(ierr);
          
          if (this->m_GradientState) {
            ierr = this->m_GradientState[j*nc + k]->RestoreArrays(); CHKERRQ(ierr);
            ierr = this->m_GradientState[(j+1)*nc + k]->GetArraysReadWrite(kernel.pGm); CHKERRQ(ierr);
          } else {
            ierr = this->m_StateVariable->GetArrayRead(pM, k, j+1, 0); CHKERRQ(ierr);
            // compute gradient for state variable at next time time point
            ierr = this->m_Differentiation->Gradient(kernel.pGm, pM); CHKERRQ(ierr);
          }

          // second part of time integration
          //ierr = kernel.TimeIntegrationPart2(); CHKERRQ(ierr);
          ierr = kernel.TimeIntegrationAll(); CHKERRQ(ierr);
          
          if (this->m_GradientXState) {
            ierr = this->m_GradientXState[j*nc + k]->RestoreArrays(); CHKERRQ(ierr);
          }
      } // for all time points
      if (this->m_GradientState) {
        ierr = this->m_GradientState[nt*nc + k]->RestoreArrays(); CHKERRQ(ierr);
      }
    }

    ierr = this->m_IncVelocityField->RestoreArrays(); CHKERRQ(ierr);
    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[2]->RestoreArrays(); CHKERRQ(ierr);
      ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    }
    if (!this->m_GradientXState) {
      ierr = this->m_WorkVecField[0]->RestoreArrays(); CHKERRQ(ierr);
    }
    ierr = this->m_WorkVecField[1]->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_IncStateVariable->RestoreArray(); CHKERRQ(ierr);
    
    ierr = this->m_IncStateVariable->DebugInfo("inc state post", __LINE__, __FILE__); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveIncAdjointProblem() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    
    if (this->m_Opt->m_OptPara.method == GAUSSNEWTON) {   // gauss newton
        ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
        if (VelocityIsZero) {
            ierr = SuperClass::SolveIncAdjointProblem(); CHKERRQ(ierr);
        } else {
            ierr = this->SolveIncAdjointEquationGN(); CHKERRQ(ierr);
        }
    } else if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
    } else {
        ierr = ThrowError("update method not defined"); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem with Gauss-Newton
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveIncAdjointEquationGN() {
    PetscErrorCode ierr = 0;
    const ScalarType *pM = nullptr;
    TransportKernelAdjointSL kernel;
    IntType nt, nc;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[2] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    kernel.scale = kernel.ht/static_cast<ScalarType>(nc);

    if (this->m_SemiLagrangianMethod == nullptr) {
        ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    } else if (this->m_Opt->m_Domain.level > 0) {
      ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
      ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }

    // compute divergence of velocity field
    ierr = this->m_WorkScaField[0]->GetArrayWrite(kernel.pDivV); CHKERRQ(ierr);
    ierr = this->m_Differentiation->Divergence(kernel.pDivV,this->m_VelocityField); CHKERRQ(ierr);
    
    ierr = this->m_WorkScaField[1]->GetArrayWrite(kernel.pDivVx); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pDivVx, kernel.pDivV, "adjoint"); CHKERRQ(ierr);
    
    ierr = this->m_WorkScaField[2]->GetArrayWrite(kernel.pLx); CHKERRQ(ierr);
    
    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pGm); CHKERRQ(ierr);
    }

    // initialize work vec field
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArraysReadWrite(kernel.pB); CHKERRQ(ierr);
    
    for (IntType j = 0; j < nt; ++j) {
        if (j == 0) kernel.scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {
            ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pL, k); CHKERRQ(ierr);
            if (!this->m_GradientState) {
              ierr = this->m_StateVariable->GetArrayRead(pM, k, nt-j); CHKERRQ(ierr);
            }
            kernel.pLnext = kernel.pL;
            
            ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pLx, kernel.pL, "adjoint"); CHKERRQ(ierr);

            // compute gradient of m^j
            if (this->m_GradientState) {
              ierr = this->m_GradientState[(nt-j)*nc + k]->GetArraysReadWrite(kernel.pGm); CHKERRQ(ierr);
            } else {
              ierr = this->m_Differentiation->Gradient(kernel.pGm, pM); CHKERRQ(ierr);
            }
            
            // compute body force
            ierr = kernel.ComputeBodyForcePart1();
            
            if (this->m_GradientState) {
              ierr = this->m_GradientState[(nt-j)*nc + k]->RestoreArrays();
            }
        }  // for all image components
        if (j == 0) kernel.scale *= 2.0;
    }  // for all time points
    
    kernel.scale *= 0.5;

    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pL, k); CHKERRQ(ierr);
        
        // compute gradient of m (for incremental body force)
        if (this->m_GradientState) {
          ierr = this->m_GradientState[k]->GetArraysReadWrite(kernel.pGm); CHKERRQ(ierr);
        } else {
          ierr = this->m_StateVariable->GetArrayRead(pM, k); CHKERRQ(ierr);
          ierr = this->m_Differentiation->Gradient(kernel.pGm,pM); CHKERRQ(ierr);
        }
        
        ierr = kernel.ComputeBodyForcePart2(); CHKERRQ(ierr);
        
        if (this->m_GradientState) {
          ierr = this->m_GradientState[k]->RestoreArrays();
        }
    }

    ierr = this->m_IncAdjointVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[2]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[0]->RestoreArray(); CHKERRQ(ierr);
    
    if (!this->m_GradientState) {
      ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
      ierr = this->m_WorkVecField[0]->RestoreArrays(); CHKERRQ(ierr);
    }

    ierr = this->m_WorkVecField[1]->RestoreArrays(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

}  // namespace reg

#endif
