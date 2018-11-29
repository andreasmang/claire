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

template<typename KernelFn, typename ... Args>
PetscErrorCode SpectralKernelCall(IntType nstart[3], IntType nx[3], IntType nl[3], Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  int3 nx3, nl3;
  nx3.x = nx[0];
  nx3.y = nx[1];
  nx3.z = nx[2];
  nl3.x = nl[0];
  nl3.y = nl[1];
  nl3.z = nl[2];
  
#pragma omp parallel
{
#pragma omp for
    for (IntType i1 = 0; i1 < nl[0]; ++i1) {
        for (IntType i2 = 0; i2 < nl[1]; ++i2) {
            for (IntType i3 = 0; i3 < nl[2]; ++i3) {
                int3 w;
                w.x = i1 + nstart[0];
                w.y = i2 + nstart[1];
                w.z = i3 + nstart[2];
                
                ComputeWaveNumber(w, nx3);
                IntType i = GetLinearIndex(i1, i2, i3, nl3);
                
                KernelFn::call(i, w, args...);
            }
        }
    }
}  // pragma omp parallel
  
  PetscFunctionReturn(ierr);
}

namespace reg {
namespace DifferentiationKernel {
  
PetscErrorCode VectorField::LaplacianTol(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ScalarType lognx = 0.;
  lognx += log2(static_cast<ScalarType>(nx[0]));
  lognx += log2(static_cast<ScalarType>(nx[1]));
  lognx += log2(static_cast<ScalarType>(nx[2]));
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCall<NLaplacianFilterKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, tol*lognx); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCall<RelaxedNLaplacianKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}
  
PetscErrorCode VectorField::Laplacian(ScalarType b0, ScalarType b1) {
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

PetscErrorCode VectorField::Bilaplacian(ScalarType b0, ScalarType b1) {
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

PetscErrorCode VectorField::Trilaplacian(ScalarType b0, ScalarType b1) {
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

PetscErrorCode VectorField::TrilaplacianFunctional(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ThrowError("relaxed trilaplacian operator not implemented"); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseLaplacian(bool usesqrt, ScalarType b0, ScalarType b1) {
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

PetscErrorCode VectorField::InverseBilaplacian(bool usesqrt, ScalarType b0, ScalarType b1) {
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

PetscErrorCode VectorField::InverseTrilaplacian(bool usesqrt, ScalarType b0, ScalarType b1) {
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

PetscErrorCode VectorField::Leray(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = SpectralKernelCall<LerayKernel>(nstart, nx, nl, 
    pXHat[0], pXHat[1], pXHat[2], 
    scale, b0, b1); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

} // namespace DifferentiationKernel
} // namespace reg

#endif
