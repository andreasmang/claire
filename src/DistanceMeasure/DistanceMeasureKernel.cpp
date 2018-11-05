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

#include "DistanceMeasureKernel.hpp"

namespace reg {
namespace DistanceMeasureKernel {
  
PetscErrorCode EvaluateFunctionalSL2::ComputeFunctionalMask() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  value = 0.0;
  for (IntType k = 0; k < nc; ++k) {  // for all image components
      for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
          ScalarType dr = (pMr[k*nl+i] - pM[k*nl+i]);
          value += pW[i]*dr*dr;
      }
  }
  
  PetscFunctionReturn(ierr);
}
  
PetscErrorCode EvaluateFunctionalSL2::ComputeFunctional() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  value = 0.0;
  for (IntType i = 0; i < nc*nl; ++i) {
      ScalarType dr = (pMr[i] - pM[i]);
      value += dr*dr;
  }
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionSL2::ComputeFinalConditionAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nc*nl; ++i) {
    pL[i] = pMr[i] - pM[i];
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionSL2::ComputeFinalConditionMaskAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType k = 0; k < nc; ++k) {  // for all image components
    for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
      pL[k*nl+i] = pW[i]*(pMr[k*nl+i] - pM[k*nl+i]);
    }
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionSL2::ComputeFinalConditionIAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl*nc; ++i) {
    pL[i] = -pM[i]; // / static_cast<ScalarType>(nc);
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionSL2::ComputeFinalConditionMaskIAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType k = 0; k < nc; ++k) {  // for all image components
    for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
      pL[k*nl+i] = -pW[i]*pM[k*nl+i];
    }
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}

} // namespace DistanceMeasureKernel
} // namespace reg
