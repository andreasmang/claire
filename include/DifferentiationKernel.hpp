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

struct DifferentiationKernel {
  ComplexType *pXHat[3];
  
  IntType nx[3];
  IntType nl[3];
  IntType nstart[3];
  
  ScalarType scale;
  
  ScalarType tol;
  
  PetscErrorCode ScalarLaplacian(ScalarType);
  PetscErrorCode LaplacianMod(ScalarType);
  PetscErrorCode Laplacian(ScalarType, ScalarType=0.0);
  PetscErrorCode LaplacianTol(ScalarType, ScalarType=0.0);
  PetscErrorCode Bilaplacian(ScalarType, ScalarType=0.0);
  PetscErrorCode Trilaplacian(ScalarType, ScalarType=0.0);
  PetscErrorCode InverseLaplacian(bool, ScalarType, ScalarType=0.0);
  PetscErrorCode InverseBilaplacian(bool, ScalarType, ScalarType=0.0);
  PetscErrorCode InverseTrilaplacian(bool, ScalarType, ScalarType=0.0);
  PetscErrorCode TrilaplacianFunctional(ScalarType, ScalarType=0.0);
  
  PetscErrorCode GaussianFilter(const ScalarType*);
  
  PetscErrorCode Gradient();
  PetscErrorCode Divergence();
  
  PetscErrorCode Leray(ScalarType, ScalarType);
  PetscErrorCode InvRegLeray(ScalarType, ScalarType, ScalarType);
};

} // namespace reg

#endif
