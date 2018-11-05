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

#ifndef _DIFFERENTIATIONKERNEL_HPP_
#define _DIFFERENTIATIONKERNEL_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"

namespace reg {
namespace DifferentiationKernel {

struct VectorField {
  ComplexType *pXHat[3];
  
  IntType nx[3];
  IntType nl[3];
  IntType nstart[3];
  
  ScalarType scale;
  
  PetscErrorCode Laplacian(ScalarType);
  PetscErrorCode Laplacian(ScalarType, ScalarType);
  PetscErrorCode Bilaplacian(ScalarType);
  PetscErrorCode Bilaplacian(ScalarType, ScalarType);
  PetscErrorCode InverseBilaplacian(ScalarType);
  PetscErrorCode InverseBilaplacianSqrt(ScalarType);
  PetscErrorCode InverseBilaplacian(ScalarType, ScalarType);
  PetscErrorCode InverseBilaplacianSqrt(ScalarType, ScalarType);
};

} // namespace DifferentiationKernel
} // namespace reg

#endif