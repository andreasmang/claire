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

    this->m_Opt = NULL;
    this->m_WorkVecField = NULL;

    this->m_v1hat = NULL;
    this->m_v2hat = NULL;
    this->m_v3hat = NULL;

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

    ierr = Assert(xhat1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xhat2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xhat3 != NULL, "null pointer"); CHKERRQ(ierr);

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

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_WorkVecField = v;

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Regularization::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _REGULARIZATIONREGISTRATION_CPP_
