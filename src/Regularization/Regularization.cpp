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

    ierr = Assert(xhat1 != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xhat2 != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xhat3 != nullptr, "null pointer"); CHKERRQ(ierr);

    this->m_v1hat = xhat1;
    this->m_v2hat = xhat2;
    this->m_v3hat = xhat3;

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

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Regularization::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (this->m_Differentiation != nullptr) {
      delete this->m_Differentiation;
      this->m_Differentiation = nullptr;
    }

    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _REGULARIZATIONREGISTRATION_CPP_
