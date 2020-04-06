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

#ifndef _SPECTRALKERNEL_CPP_
#define _SPECTRALKERNEL_CPP_

#include "SpectralKernel.hpp"
#include "cuda_helper.hpp"

#include "SpectralKernel.txx"

using KernelUtils::SpectralKernelCall;

namespace reg {

PetscErrorCode SpectralKernel::LowPassFilter(ComplexType *pXHat, ScalarType pct) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ScalarType l1, l2, l3;
  l1 = static_cast<ScalarType>(nx[0])*0.5*pct;
  l2 = static_cast<ScalarType>(nx[1])*0.5*pct;
  l3 = static_cast<ScalarType>(nx[2])*0.5*pct;
  ierr = SpectralKernelCall<LowPassFilterKernel>(nstart, nx, nl, pXHat, l1, l2, l3, scale); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode SpectralKernel::HighPassFilter(ComplexType *pXHat, ScalarType pct) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ScalarType l1, l2, l3;
  l1 = static_cast<ScalarType>(nx[0])*0.5*pct;
  l2 = static_cast<ScalarType>(nx[1])*0.5*pct;
  l3 = static_cast<ScalarType>(nx[2])*0.5*pct;
  ierr = SpectralKernelCall<HighPassFilterKernel>(nstart, nx, nl, pXHat, l1, l2, l3, scale); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode SpectralKernel::Scale(ComplexType *pX, ScalarType val) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = SpectralKernelCall<ScaleKernel>(nstart, nx, nl, pX, scale*val); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode SpectralKernel::Restrict(ComplexType *pXc, const ComplexType *pXf, const IntType nxc[3], const IntType ostart_c[3]) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode SpectralKernel::Prolong(ComplexType *pXf, const ComplexType *pXc, const IntType nxc[3],  const IntType ostart_c[3]) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  PetscFunctionReturn(ierr);
}

} // namepsace reg

#endif
