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

#ifndef _DIFFERENTIATIONKERNEL_CPP_
#define _DIFFERENTIATIONKERNEL_CPP_

#include "DifferentiationKernel.hpp"
#include "CLAIREUtils.hpp"

#include "DifferentiationKernel.txx"

namespace reg {
  
using KernelUtils::SpectralKernelCall;
  
PetscErrorCode DifferentiationKernel::LaplacianTol(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ScalarType lognx = 0.;
  lognx += log2(static_cast<ScalarType>(nx[0]));
  lognx += log2(static_cast<ScalarType>(nx[1]));
  lognx += log2(static_cast<ScalarType>(nx[2]));
  
  KernelUtils::array3_t<ComplexType*> v;
  v.x = pXHat[0];
  v.y = pXHat[1];
  v.z = pXHat[2];
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCall<NLaplacianFilterKernel<1> >(nstart, nx, nl, v, 
      b0*scale, tol*lognx); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCall<RelaxedNLaplacianKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::ScalarLaplacian(ScalarType b0) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = SpectralKernelCall<NLaplacianKernel<1> >(nstart, nx, nl, 
    pXHat[0], b0*scale); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationKernel::LaplacianMod(ScalarType b0) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = SpectralKernelCall<NLaplacianModKernel<1> >(nstart, nx, nl, 
    pXHat[0], pXHat[1], pXHat[2], scale, b0); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationKernel::Laplacian(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCall<NLaplacianKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCall<RelaxedNLaplacianKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::Bilaplacian(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCall<NLaplacianKernel<2> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCall<RelaxedNLaplacianKernel<2> >(nstart, nx, nl,
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::Trilaplacian(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCall<NLaplacianKernel<3> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCall<RelaxedNLaplacianKernel<3> >(nstart, nx, nl,
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::TrilaplacianFunctional(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ThrowError("relaxed trilaplacian operator not implemented"); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::InverseLaplacian(bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (usesqrt) {
    if (b1 == 0.0) {
      ierr = SpectralKernelCall<InverseNLaplacianSqrtKernel<1> >(nstart, nx, nl,
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCall<RelaxedInverseNLaplacianSqrtKernel<1> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  } else {
    if (b1 == 0.0) {
      ierr = SpectralKernelCall<InverseNLaplacianKernel<1> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2], 
        scale, b0); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCall<RelaxedInverseNLaplacianKernel<1> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::InverseBilaplacian(bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (usesqrt) {
    if (b1 == 0.0) {
      /// scale/sqrt(b0*|lapik|^2) = scale/(sqrt(b0)*|lapik|)
      ierr = SpectralKernelCall<InverseNLaplacianKernel<1> >(nstart, nx, nl,
        pXHat[0], pXHat[1], pXHat[2],
        scale, sqrt(b0)); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCall<RelaxedInverseNLaplacianSqrtKernel<2> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  } else {
    if (b1 == 0.0) {
      ierr = SpectralKernelCall<InverseNLaplacianKernel<2> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2], 
        scale, b0); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCall<RelaxedInverseNLaplacianKernel<2> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::InverseTrilaplacian(bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (usesqrt) {
    if (b1 == 0.0) {
      ierr = SpectralKernelCall<InverseNLaplacianSqrtKernel<3> >(nstart, nx, nl,
        pXHat[0], pXHat[1], pXHat[2],
        scale, sqrt(b0)); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCall<RelaxedInverseNLaplacianSqrtKernel<3> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  } else {
    if (b1 == 0.0) {
      ierr = SpectralKernelCall<InverseNLaplacianKernel<3> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2], 
        scale, b0); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCall<RelaxedInverseNLaplacianKernel<3> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::Leray(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = SpectralKernelCall<LerayKernel>(nstart, nx, nl, 
    pXHat[0], pXHat[1], pXHat[2], 
    scale, b0, b1); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::GaussianFilter(const ScalarType c[3]) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = SpectralKernelCall<GaussianFilterKernel>(nstart, nx, nl, 
    pXHat[0], c[0], c[1], c[2], scale); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::Gradient() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = SpectralKernelCall<GradientKernel>(nstart, nx, nl, 
    pXHat[0], pXHat[1], pXHat[2], scale); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationKernel::Divergence() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = SpectralKernelCall<DivergenceKernel>(nstart, nx, nl, 
    pXHat[0], pXHat[1], pXHat[2], scale); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

} // namespace reg

#endif
