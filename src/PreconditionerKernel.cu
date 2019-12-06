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

#include "PreconditionerKernel.hpp"
#include "cuda_helper.hpp"

#include "PreconditionerKernel.txx"

namespace reg {
  
using KernelUtils::KernelCallGPU;

PetscErrorCode H0PrecondKernel::Iteration () {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = KernelCallGPU<H0Kernel>(nl, pVhat, pGmt, pRHS, pReg, omg); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode H0PrecondKernel::IterationPart1 () {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = KernelCallGPU<H0Kernel>(nl, pM[0], pM[1], pM[2], pVhat[0], pVhat[1], pVhat[2], pGmt[0], pGmt[1], pGmt[2]); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode H0PrecondKernel::IterationPart2 () {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = KernelCallGPU<H0Kernel>(nl, pVhat, pM, omg, beta); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode H0PrecondKernel::Diagonal () {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = KernelCallGPU<H0Kernel>(nl, pVhat, pGmt, beta); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

} // namespace reg
