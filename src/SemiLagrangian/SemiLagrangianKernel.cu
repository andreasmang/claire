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

#include "SemiLagrangianKernel.hpp"
#include "cuda_helper.hpp"
#include "SemiLagrangianKernel.txx"

namespace reg {
  
using KernelUtils::SpacialKernelCallGPU;
  
PetscErrorCode TrajectoryKernel::RK2_Step1() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = SpacialKernelCallGPU<RK2Kernel>(istart, isize,
                                         pX[0], pX[1], pX[2],
                                         pV[0], pV[1], pV[2],
                                         ix[0], ix[1], ix[2],
                                         hx[0], hx[1], hx[2]); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode TrajectoryKernel::RK2_Step2() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  const ScalarType half = 0.5;

  ierr = SpacialKernelCallGPU<RK2Kernel>(istart, isize,
                                         pX[0], pX[1], pX[2],
                                         pV[0], pV[1], pV[2],
                                         pVx[0], pVx[1], pVx[2],
                                         ix[0], ix[1], ix[2],
                                         hx[0]*half, hx[1]*half, hx[2]*half); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

} // namespace reg
