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

#ifndef _DIFFERENTIATIONSMGPU_CPP_
#define _DIFFERENTIATIONSMGPU_CPP_




#ifdef REG_HAS_CUDA




#include "DifferentiationSMGPU.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DifferentiationSM::DifferentiationSM() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DifferentiationSM::DifferentiationSM(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DifferentiationSM::~DifferentiationSM() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DifferentiationSM::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DifferentiationSM::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Gradient(ScalarType *g1,
                                           ScalarType *g2,
                                           ScalarType *g3,
                                           ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply laplacian operator to scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Laplacian(ScalarType *l,
                                            ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply laplacian operator to vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Laplacian(ScalarType *l,
                                            ScalarType *v1,
                                            ScalarType *v2,
                                            ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply biharmonic operator to vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Biharmonic(ScalarType *b,
                                            ScalarType *v1,
                                            ScalarType *v2,
                                            ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




}  // end of name space




#endif  // REG_HAS_CUDA




#endif  // _DIFFERENTIATIONSMGPU_CPP_
