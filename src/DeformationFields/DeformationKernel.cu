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

#include "DeformationKernel.hpp"
#include "cuda_helper.hpp"
#include "DeformationKernel.txx"

namespace reg {
  
PetscErrorCode DetDefGradKernel::IntegrateSL() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = KernelUtils::KernelCallGPU<DetDefGradSLKernel>(nl, pJ, pJx, pDivV, pDivVx, alpha, ht); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DetDefGradKernel::InitSL(ScalarType val) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = KernelUtils::KernelCallGPU<DetDefGradSLKernel>(nl, pJ, val); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

} // namespace reg
