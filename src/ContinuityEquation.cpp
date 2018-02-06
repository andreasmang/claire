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

#include "ContinuityEquation.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
ContinuityEquation::ContinuityEquation() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
ContinuityEquation::~ContinuityEquation() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
ContinuityEquation::ContinuityEquation(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode ContinuityEquation::Initialize() {
    PetscFunctionBegin;

    this->m_TemplateImage = NULL;
    this->m_ReferenceImage = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode ContinuityEquation::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem (i.e., the continuity equation)
 *******************************************************************/
PetscErrorCode ContinuityEquation::SolveForwardProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode ContinuityEquation::SolveAdjointProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode ContinuityEquation::SolveIncForwardProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode ContinuityEquation::SolveIncAdjointProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _CONTINUITYEQUATION_CPP_
