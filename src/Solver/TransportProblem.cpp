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

#ifndef _CONTINUITYEQUATION_CPP_
#define _CONTINUITYEQUATION_CPP_

#include "TransportProblem.hpp"
#include "DifferentiationSM.hpp"
#include "DifferentiationFD.hpp"

namespace reg {

/********************************************************************
 * @brief default constructor
 *******************************************************************/
TransportProblem::TransportProblem() {
    this->Initialize();
}

/********************************************************************
 * @brief default destructor
 *******************************************************************/
TransportProblem::~TransportProblem() {
    this->ClearMemory();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
TransportProblem::TransportProblem(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode TransportProblem::Initialize() {
    PetscFunctionBegin;

    this->m_Opt = nullptr;
    this->m_TemplateImage = nullptr;
    this->m_ReferenceImage = nullptr;
    this->m_StateVariable = nullptr;
    this->m_AdjointVariable = nullptr;
    this->m_IncStateVariable = nullptr;
    this->m_IncAdjointVariable = nullptr;
    
    this->m_VelocityField = nullptr;
    this->m_IncVelocityField = nullptr;
    
    this->m_GradientState = nullptr;
    this->m_GradientXState = nullptr;
    
    for (int i = 0;i < 5;++i) {
      this->m_WorkScaField[i] = nullptr;
      this->m_WorkVecField[i] = nullptr;
    }
    
    this->m_Differentiation = nullptr;
    this->m_DiffAllocated = false;

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode TransportProblem::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (this->m_DiffAllocated) {
      Free(this->m_Differentiation);
    }
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode TransportProblem::SetReferenceImage(ScaField* mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mR != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_ReferenceImage = mR;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set template image (i.e., the image to be deformed)
 *******************************************************************/
PetscErrorCode TransportProblem::SetTemplateImage(ScaField* mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mT != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_TemplateImage = mT;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set state variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetStateVariable(ScaField* m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_StateVariable = m;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set gradient state variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetGradientState(VecField** g) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    this->m_GradientState = g;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set interpolated gradient state variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetGradientXState(VecField** g) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    this->m_GradientXState = g;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set adjoint variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetAdjointVariable(ScaField* lambda) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_AdjointVariable = lambda;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set velocity field
 *******************************************************************/
PetscErrorCode TransportProblem::SetControlVariable(VecField *field) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_VelocityField = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

PetscErrorCode TransportProblem::InitializeControlVariable(VecField *field) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_VelocityField = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set incremental adjoint variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetIncAdjointVariable(ScaField* m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_IncAdjointVariable = m;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set incremental state variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetIncStateVariable(ScaField* m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_IncStateVariable = m;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set incremental velocity field
 *******************************************************************/
PetscErrorCode TransportProblem::SetIncControlVariable(VecField *field) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_IncVelocityField = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set working scalar field
 ********************************************************************/
PetscErrorCode TransportProblem::SetWorkScaField(ScaField* field, IntType idx) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(idx >= 1 && idx <= 5, "index out of range"); CHKERRQ(ierr);
    
    this->m_WorkScaField[idx-1] = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set working scalar field
 ********************************************************************/
PetscErrorCode TransportProblem::SetWorkVecField(VecField* field, IntType idx) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(idx >= 1 && idx <= 5, "index out of range"); CHKERRQ(ierr);
    
    this->m_WorkVecField[idx-1] = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set the differentiation method
 *******************************************************************/
PetscErrorCode TransportProblem::SetDifferentiation(Differentiation::Type type) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);

    if (this->m_Differentiation != nullptr && this->m_Differentiation->m_Type != type && this->m_DiffAllocated) {
      Free(this->m_Differentiation);
    }
    
    if (this->m_Differentiation == nullptr) {
      switch (type) {
      case Differentiation::Type::Spectral:
        ierr = AllocateOnce<DifferentiationSM>(this->m_Differentiation, this->m_Opt); CHKERRQ(ierr);
        this->m_DiffAllocated = true;
        break;
      case Differentiation::Type::Finite:
        ierr = AllocateOnce<DifferentiationFD>(this->m_Differentiation, this->m_Opt); CHKERRQ(ierr);
        this->m_DiffAllocated = true;
        break;
      default:
        ierr = ThrowError("no valid differentiation method"); CHKERRQ(ierr);
        break;
      };
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set the differentiation method
 *******************************************************************/
PetscErrorCode TransportProblem::SetDifferentiation(Differentiation *diff) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);

    if (this->m_DiffAllocated) {
      Free(this->m_Differentiation);
      this->m_DiffAllocated = false;
    }
    
    this->m_Differentiation = diff;
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the forward problem for zero velocity fields
 *******************************************************************/
PetscErrorCode TransportProblem::SolveForwardProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    
    // m and \lambda are constant in time
    
    if (this->m_Opt->m_RegFlags.runinversion) {
        ierr = this->m_StateVariable->CopyFrame(0); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the inverse forward problem for zero velocity fields
 *******************************************************************/
PetscErrorCode TransportProblem::SolveInverseProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    IntType nt = this->m_Opt->m_Domain.nt;
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    
    // m and \lambda are constant in time
    
    if (this->m_Opt->m_RegFlags.runinversion) {
        ierr = this->m_StateVariable->CopyFrame(nt); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem for zero velocity fields
 *******************************************************************/
PetscErrorCode TransportProblem::SolveAdjointProblem() {
    PetscErrorCode ierr = 0;
    const ScalarType *pM = nullptr;
    TransportKernelAdjoint kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    IntType nt = this->m_Opt->m_Domain.nt;
    IntType nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.scale = 1.0/static_cast<ScalarType>(nc);
    
    ierr = Assert(this->m_TemplateImage != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);
    
    //ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL); CHKERRQ(ierr);
    // adjoint variable is constant in time
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
      this->m_AdjointVariable->CopyFrame(nt); CHKERRQ(ierr);
        /*for (IntType j = 1; j <= nt; ++j) {
            TransportKernelCopy(kernel.pL + nt*nc*nl, kernel.pL + (nt-j)*nl*nc, nc*nl);
        }*/
    }

    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    
    // m and \lambda are constant in time
    ierr = this->m_WorkVecField[1]->GetArraysReadWrite(kernel.pB); CHKERRQ(ierr);
    
    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pGm); CHKERRQ(ierr);
    }

    
    for (IntType k = 0; k < nc; ++k) {  // for all components
        ierr = this->m_TemplateImage->GetArrayRead(pM, k); CHKERRQ(ierr);
        ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL, k); CHKERRQ(ierr);
        
        // compute gradient of m
        if (this->m_GradientState) {
          ierr = this->m_GradientState[k]->GetArraysReadWrite(kernel.pGm); CHKERRQ(ierr);
        } else {
          ierr = this->m_Differentiation->Gradient(kernel.pGm, pM); CHKERRQ(ierr);
        }

        // b = \sum_k\int_0^1 \lambda_k \grad m_k dt
        ierr = kernel.ComputeBodyForce(); CHKERRQ(ierr);
        
        if (this->m_GradientState) {
          ierr = this->m_GradientState[k]->RestoreArrays(); CHKERRQ(ierr);
        }
    }
   
    ierr = this->m_WorkVecField[1]->RestoreArrays(); CHKERRQ(ierr);
    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[0]->RestoreArrays(); CHKERRQ(ierr);
    }
    ierr = this->m_TemplateImage->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_AdjointVariable->RestoreArray(); CHKERRQ(ierr);
      
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the incremental adjoint problem for zero velocity fields with Gauss-Newton-scheme
 *******************************************************************/
PetscErrorCode TransportProblem::SolveIncAdjointProblem() {
    PetscErrorCode ierr = 0;
    const ScalarType *pM = nullptr;
    TransportKernelAdjoint kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    IntType nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.scale = 1.0/static_cast<ScalarType>(nc);
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);

    // m and \lambda are constant in time
    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pGm); CHKERRQ(ierr);
    }

    // init body force for numerical integration
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArraysReadWrite(kernel.pB); CHKERRQ(ierr);

    // $m$ and $\tilde{\lambda}$ are constant
    for (IntType k = 0; k < nc; ++k) {  // for all components
        ierr = this->m_StateVariable->GetArrayRead(pM, k);
        ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pL, k); CHKERRQ(ierr);
        
        // compute gradient of m
        if (this->m_GradientState) {
          ierr = this->m_GradientState[k]->GetArraysReadWrite(kernel.pGm); CHKERRQ(ierr);
        } else {
          ierr = this->m_Differentiation->Gradient(kernel.pGm, pM); CHKERRQ(ierr);
        }

        ierr = kernel.ComputeBodyForce(); CHKERRQ(ierr);
        
        if (this->m_GradientState) {
          ierr = this->m_GradientState[k]->RestoreArrays(); CHKERRQ(ierr);
        }
    }
    ierr = this->m_WorkVecField[1]->RestoreArrays(); CHKERRQ(ierr);
    if (!this->m_GradientState) {
      ierr = this->m_WorkVecField[0]->RestoreArrays(); CHKERRQ(ierr);
    }
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_IncAdjointVariable->RestoreArray(); CHKERRQ(ierr);
  
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

}  // namespace reg

#endif  // _CONTINUITYEQUATION_CPP_
