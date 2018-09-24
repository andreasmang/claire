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
    
    for (int i = 0;i < 5;++i) {
      this->m_WorkScaField[i] = nullptr;
      this->m_WorkVecField[i] = nullptr;
    }
    
    this->m_Differentiation = nullptr;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode TransportProblem::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (this->m_Differentiation != nullptr) {
      delete this->m_Differentiation;
      this->m_Differentiation = nullptr;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode TransportProblem::SetReferenceImage(Vec mR) {
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
PetscErrorCode TransportProblem::SetTemplateImage(Vec mT) {
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
PetscErrorCode TransportProblem::SetStateVariable(Vec m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_StateVariable = m;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set adjoint variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetAdjointVariable(Vec lambda) {
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

/********************************************************************
 * @brief set incremental adjoint variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetIncAdjointVariable(Vec m) {
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
PetscErrorCode TransportProblem::SetIncStateVariable(Vec m) {
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
PetscErrorCode TransportProblem::SetWorkScaField(Vec field, IntType idx) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(idx >= 0 && idx < 5, "index out of range"); CHKERRQ(ierr);
    
    this->m_WorkScaField[idx] = field;

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
    ierr = Assert(idx >= 0 && idx < 5, "index out of range"); CHKERRQ(ierr);
    
    this->m_WorkVecField[idx] = field;

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

    if (this->m_Differentiation != nullptr && this->m_Differentiation->m_Type != type) {
      delete this->m_Differentiation;
      this->m_Differentiation = nullptr;
    }
    
    if (this->m_Differentiation == nullptr) {
      switch (type) {
      case Differentiation::Type::Spectral:
        try {
          this->m_Differentiation = new DifferentiationSM(this->m_Opt);
        } catch (std::bad_alloc& err) {
          ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        break;
      case Differentiation::Type::Finite:
        try {
          this->m_Differentiation = new DifferentiationFD(this->m_Opt);
        } catch (std::bad_alloc& err) {
          ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
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
 * @brief solve the forward problem for zero velocity fields
 *******************************************************************/
PetscErrorCode TransportProblem::SolveForwardProblem() {
    PetscErrorCode ierr = 0;
    ScalarType *pM = nullptr, *pL = nullptr;
    IntType nc, nl, nt;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    nt = this->m_Opt->m_Domain.nt;
    
    ierr = Assert(this->m_TemplateImage != nullptr, "null pointer"); CHKERRQ(ierr);
    
    // m and \lambda are constant in time
    
    if (this->m_Opt->m_RegFlags.runinversion) {
        ierr = GetRawPointer(this->m_StateVariable, &pM); CHKERRQ(ierr);
        for (IntType j = 1; j <= nt; ++j) {
            TransportKernelCopy(pM, pM + j*nc*nl, nc*nl);
        }
        ierr = RestoreRawPointer(this->m_StateVariable, &pM); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem for zero velocity fields
 *******************************************************************/
PetscErrorCode TransportProblem::SolveAdjointProblem() {
    PetscErrorCode ierr = 0;
    ScalarType *pM = nullptr, *pL = nullptr;
    IntType nc, nl, nt;
    TransportKernelAdjoint kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    nt = this->m_Opt->m_Domain.nt;
    kernel.nl = nl;
    kernel.scale = 1.0/static_cast<ScalarType>(nc);
    
    ierr = Assert(this->m_TemplateImage != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = GetRawPointer(this->m_AdjointVariable, &pL); CHKERRQ(ierr);
    // adjoint variable is constant in time
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        for (IntType j = 1; j <= nt; ++j) {
            TransportKernelCopy(pL + nt*nc*nl, pL + (nt-j)*nl*nc, nc*nl);
        }
    }

    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);

    // m and \lambda are constant in time
    ierr = GetRawPointer(this->m_TemplateImage, &pM); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArrays(kernel.pGm); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArrays(kernel.pB); CHKERRQ(ierr);
    
    for (IntType k = 0; k < nc; ++k) {  // for all components
        kernel.pL = pL + k*nl;
        
        // compute gradient of m
        ierr = this->m_Differentiation->Gradient(kernel.pGm,pM+k*nl); CHKERRQ(ierr);

        // b = \sum_k\int_{\Omega} \lambda_k \grad m_k dt
        ierr = kernel.ComputeBodyForce(); CHKERRQ(ierr);
    }  // pragma
   
    ierr = this->m_WorkVecField[1]->RestoreArrays(kernel.pB); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->RestoreArrays(kernel.pGm); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_TemplateImage, &pM); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &pL); CHKERRQ(ierr);
  
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the incremental adjoint problem for zero velocity fields with Gauss-Newton-scheme
 *******************************************************************/
PetscErrorCode TransportProblem::SolveIncAdjointProblem() {
    PetscErrorCode ierr = 0;
    ScalarType *pM = nullptr, *pLtilde = nullptr;
    IntType nc, nl;
    TransportKernelAdjoint kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    kernel.nl = nl;
    kernel.scale = 1.0/static_cast<ScalarType>(nc);
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);

    // copy terminal condition \tilde{\lambda}_1 = -\tilde{m}_1 to all time points
    ierr = GetRawPointer(this->m_IncAdjointVariable, &pLtilde); CHKERRQ(ierr);

    // m and \lambda are constant in time
    ierr = GetRawPointer(this->m_StateVariable, &pM); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArrays(kernel.pGm); CHKERRQ(ierr);

    // init body force for numerical integration
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArrays(kernel.pB); CHKERRQ(ierr);

    // $m$ and $\tilde{\lambda}$ are constant
    for (IntType k = 0; k < nc; ++k) {  // for all components
        kernel.pL = pLtilde + k*nl;
        // compute gradient of m
        ierr = this->m_Differentiation->Gradient(kernel.pGm,pM + k*nl); CHKERRQ(ierr);

        ierr = kernel.ComputeBodyForce(); CHKERRQ(ierr);
    }
    ierr = this->m_WorkVecField[1]->RestoreArrays(kernel.pB); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->RestoreArrays(kernel.pGm); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &pM); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &pLtilde); CHKERRQ(ierr);
  
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


}  // namespace reg




#endif  // _CONTINUITYEQUATION_CPP_
