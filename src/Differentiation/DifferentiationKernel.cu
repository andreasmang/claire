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
#include "cuda_helper.hpp"

#include "DifferentiationKernel.txx"

template<typename KernelFn, typename ... Args>
__global__ void SpectralKernelGPU(int3 wave, int3 nx, int3 nl, Args ... args) {
  int i1 = threadIdx.x + blockIdx.x*blockDim.x;
  int i2 = blockIdx.y;
  int i3 = blockIdx.z;
  
  if (i1 < nl.x) {
    wave.x += i1;
    wave.y += i2;
    wave.z += i3;

    ComputeWaveNumber(wave, nx);
    int i = GetLinearIndex(i1, i2, i3, nl);

    KernelFn::call(i, wave, args...);
  }
}
template<typename KernelFn, typename ... Args>
PetscErrorCode SpectralKernelCallGPU(IntType nstart[3], IntType nx[3], IntType nl[3], 
    Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  dim3 block(256,1,1);
  dim3 grid((nl[0] + 255)/256,nl[1],nl[2]);
  int3 wave, nx3, nl3;
  wave.x = nstart[0]; wave.y = nstart[1]; wave.z = nstart[2];
  nx3.x = nx[0]; nx3.y = nx[1]; nx3.z = nx[2];
  nl3.x = nl[0]; nl3.y = nl[1]; nl3.z = nl[2];
  
  if (nl[0]*nl[1]*nl[2] > 0) {
    SpectralKernelGPU<KernelFn><<<grid, block>>>(wave, nx3, nl3, args...);
    ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
    ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);
  }
  
  PetscFunctionReturn(ierr);
}

namespace reg {
namespace DifferentiationKernel {
  
PetscErrorCode VectorField::Laplacian(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCallGPU<NLaplacianKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCallGPU<RelaxedNLaplacianKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::LaplacianTol(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ScalarType lognx = 0.;
  lognx += log2(static_cast<ScalarType>(nx[0]));
  lognx += log2(static_cast<ScalarType>(nx[1]));
  lognx += log2(static_cast<ScalarType>(nx[2]));
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCallGPU<NLaplacianFilterKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, tol*lognx); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCallGPU<RelaxedNLaplacianKernel<1> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Bilaplacian(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCallGPU<NLaplacianKernel<2> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCallGPU<RelaxedNLaplacianKernel<2> >(nstart, nx, nl,
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Trilaplacian(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (b1 == 0.0) {
    ierr = SpectralKernelCallGPU<NLaplacianKernel<3> >(nstart, nx, nl, 
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale); CHKERRQ(ierr);
  } else {
    ierr = SpectralKernelCallGPU<RelaxedNLaplacianKernel<3> >(nstart, nx, nl,
      pXHat[0], pXHat[1], pXHat[2], 
      b0*scale, b1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::TrilaplacianFunctional(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ThrowError("trilaplacian operator not implemented"); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseLaplacian(bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (usesqrt) {
    if (b1 == 0.0) {
      ierr = SpectralKernelCallGPU<InverseNLaplacianSqrtKernel<1> >(nstart, nx, nl,
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCallGPU<RelaxedInverseNLaplacianSqrtKernel<1> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  } else {
    if (b1 == 0.0) {
      ierr = SpectralKernelCallGPU<InverseNLaplacianKernel<1> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2], 
        scale, b0); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCallGPU<RelaxedInverseNLaplacianKernel<1> >(nstart, nx, nl, 
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
      ierr = SpectralKernelCallGPU<InverseNLaplacianKernel<1> >(nstart, nx, nl,
        pXHat[0], pXHat[1], pXHat[2],
        scale, sqrt(b0)); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCallGPU<RelaxedInverseNLaplacianSqrtKernel<2> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  } else {
    if (b1 == 0.0) {
      ierr = SpectralKernelCallGPU<InverseNLaplacianKernel<2> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2], 
        scale, b0); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCallGPU<RelaxedInverseNLaplacianKernel<2> >(nstart, nx, nl, 
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
      ierr = SpectralKernelCallGPU<InverseNLaplacianSqrtKernel<3> >(nstart, nx, nl,
        pXHat[0], pXHat[1], pXHat[2],
        scale, sqrt(b0)); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCallGPU<RelaxedInverseNLaplacianSqrtKernel<3> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  } else {
    if (b1 == 0.0) {
      ierr = SpectralKernelCallGPU<InverseNLaplacianKernel<3> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2], 
        scale, b0); CHKERRQ(ierr);
    } else {
      ierr = SpectralKernelCallGPU<RelaxedInverseNLaplacianKernel<3> >(nstart, nx, nl, 
        pXHat[0], pXHat[1], pXHat[2],
        scale, b0, b1); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Leray(ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = SpectralKernelCallGPU<LerayKernel>(nstart, nx, nl, 
    pXHat[0], pXHat[1], pXHat[2], 
    scale, b0, b1); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

} // namespace DifferentiationKernel
} // namespace reg

#endif
