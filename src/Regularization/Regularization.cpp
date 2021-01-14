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

#ifndef _REGULARIZATION_CPP_
#define _REGULARIZATION_CPP_

#include "Regularization.hpp"
#include "DifferentiationSM.hpp"
#include "DifferentiationFD.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
Regularization::Regularization() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
Regularization::~Regularization(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
Regularization::Regularization(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Regularization::Initialize(void) {
    PetscFunctionBegin;

    this->m_Opt = nullptr;
    this->m_WorkVecField = nullptr;
    this->m_Differentiation = nullptr;

    this->m_v1hat = nullptr;
    this->m_v2hat = nullptr;
    this->m_v3hat = nullptr;
    
    this->m_DiffAllocated = false;

    PetscFunctionReturn(0);
}


/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Regularization::SetSpectralData(ComplexType* xhat1,
                                               ComplexType* xhat2,
                                               ComplexType* xhat3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    //ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);
    //ierr = this->m_Differentiation->SetupData(xhat1, xhat2, xhat3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Regularization::SetWorkVecField(VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(v != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_WorkVecField = v;

    PetscFunctionReturn(ierr);
}

PetscErrorCode Regularization::SetWorkScaField(ScaField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(v != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_WorkScaField = v;

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief set Differentiation interface
 *******************************************************************/
PetscErrorCode Regularization::SetDifferentiation(Differentiation *diff) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Differentiation != nullptr && this->m_DiffAllocated) {
      delete this->m_Differentiation;
      this->m_Differentiation = nullptr;
      this->m_DiffAllocated = false;
    }
    
    this->m_Differentiation = diff;

    PetscFunctionReturn(ierr);
}
/********************************************************************
 * @brief set Differentiation interface
 *******************************************************************/
PetscErrorCode Regularization::SetDifferentiation(Differentiation::Type type) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Differentiation != nullptr && this->m_Differentiation->m_Type != type) {
      delete this->m_Differentiation;
      this->m_Differentiation = nullptr;
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

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Regularization::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (this->m_DiffAllocated) {
      ierr = Free(this->m_Differentiation); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _REGULARIZATIONREGISTRATION_CPP_
