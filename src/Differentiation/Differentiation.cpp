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

#ifndef _DIFFERENTIATION_CPP_
#define _DIFFERENTIATION_CPP_

#include "Differentiation.hpp"




namespace reg {





/********************************************************************
 * @brief default constructor
 *******************************************************************/
Differentiation::Differentiation() : m_Type(Type::None) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
Differentiation::~Differentiation() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
Differentiation::Differentiation(RegOpt* opt, Type type) : m_Type(type) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Differentiation::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Differentiation::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




}  // end of namespace




#endif  // _DIFFERENTIATION_CPP_
