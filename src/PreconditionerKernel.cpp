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
#include "PreconditionerKernel.txx"

namespace reg {
  
using KernelUtils::KernelCall;
using KernelUtils::ReductionKernelCall;

PetscErrorCode H0PrecondKernel::gMgMT2 () {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = KernelCall<H0Kernel2>(nl, 
                                 pM[0], pM[1], pM[2], 
                                 pVhat[0], pVhat[1], pVhat[2], 
                                 pGmt[0], pGmt[1], pGmt[2]); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}
PetscErrorCode H0PrecondKernel::res2 (ScalarType &res) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ReductionKernelCall<H0Kernel2>(res, nl, 
                                          pM[0], pM[1], pM[2],
                                          pP[0], pP[1], pP[2],
                                          pRes[0], pRes[1], pRes[2],
                                          pGmt[0], pGmt[1], pGmt[2],
                                          diag); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}
PetscErrorCode H0PrecondKernel::pTAp2 (ScalarType &res) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ReductionKernelCall<H0Kernel2>(res, nl, 
                                          pM[0], pM[1], pM[2],
                                          pP[0], pP[1], pP[2],
                                          pGmt[0], pGmt[1], pGmt[2],
                                          diag); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode H0PrecondKernel::gMgMT () {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = KernelCall<H0Kernel>(nl, 
                                 pM[0], pM[1], pM[2], 
                                 pVhat[0], pVhat[1], pVhat[2], 
                                 pGmt[0], pGmt[1], pGmt[2]); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode H0PrecondKernel::res (ScalarType &res) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ReductionKernelCall<H0Kernel>(res, nl, 
                                          pM[0], pM[1], pM[2],
                                          pP[0], pP[1], pP[2],
                                          pRes[0], pRes[1], pRes[2],
                                          pVhat[0], pVhat[1], pVhat[2],
                                          beta); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode H0PrecondKernel::pTAp (ScalarType &res) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ReductionKernelCall<H0Kernel>(res, nl, 
                                          pM[0], pM[1], pM[2],
                                          pP[0], pP[1], pP[2],
                                          beta); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode H0PrecondKernel::CGres (ScalarType &res) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ScalarType alpha = res;
  
  ierr = ReductionKernelCall<H0KernelCG>(res, nl, 
                                            pM[0], pM[1], pM[2],
                                            pP[0], pP[1], pP[2],
                                            pRes[0], pRes[1], pRes[2],
                                            pVhat[0], pVhat[1], pVhat[2],
                                            alpha); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode H0PrecondKernel::CGp (ScalarType alpha) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
    
  ierr = KernelCall<H0KernelCG>(nl, 
                                   pP[0], pP[1], pP[2],
                                   pRes[0], pRes[1], pRes[2],
                                   alpha); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

} // namespace reg
