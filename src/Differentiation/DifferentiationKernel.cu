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

/********************************************************************
 * @brief computes wave number on GPU
 *******************************************************************/
__device__ inline void ComputeWaveNumber(dim3 &w, dim3 &n) {
  if (w.x >  n.x/2) w.x -= n.x;
  if (w.y >  n.y/2) w.y -= n.y;
  if (w.z >  n.z/2) w.z -= n.z;
}

/********************************************************************
 * @brief computes linear array index on GPU
 *******************************************************************/
__device__ inline int GetLinearIndex(int i1, int i2, int i3, dim3 size) {
  return i1*size.y*size.z + i2*size.z + i3;
}

template<int N> __device__ inline ScalarType pow(ScalarType x) {
  return x*pow<N-1>(x);
}
template<> __device__ inline ScalarType pow<0>(ScalarType x) {
  return 1;
}

/********************************************************************
 * @brief computes linear array index on GPU
 *******************************************************************/
template<int N> 
__device__ inline ScalarType ComputeNLaplaceNumber(dim3 w) {
  ScalarType norm = static_cast<ScalarType>(w.x*w.x + w.y*w.y + w.z*w.z);
  return pow<N>(norm);
}

template<int N> struct NLaplacianKernelGPU {
  static __device__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, ScalarType b0) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w);
    ScalarType regop = b0*lapik;

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NLaplacianRegularizationKernelGPU {
  static __device__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, ScalarType b0, ScalarType b1) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w) + b1;
    ScalarType regop = b0*lapik;

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NInvLaplacianKernelGPU {
  static __device__ inline void call (int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, ScalarType scale, ScalarType b0) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w);
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NInvLaplacianSqrtKernelGPU {
  static __device__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, ScalarType scale, ScalarType b0) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w);
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/sqrt(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NInvLaplacianRegularizationKernelGPU {
  static __device__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0, ScalarType b1) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w) + b1;
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> 
struct NInvLaplacianRegularizationSqrtKernelGPU {
  static __device__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0, ScalarType b1) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w) + b1;
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/sqrt(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};

template<typename KernelFn, typename ... Args>
__global__ void SpectralKernelGPU(dim3 wave, dim3 nx, dim3 nl, Args ... args) {
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
PetscErrorCode SpectralKernelCallGPU(IntType nstart[3], IntType nx[3], IntType nl[3], Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  dim3 block(256,1,1);
  dim3 grid((nl[0] + 255)/256,nl[1],nl[2]);
  dim3 wave(nstart[0],nstart[1],nstart[2]);
  dim3 nx3(nx[0],nx[1],nx[2]);
  dim3 nl3(nl[0],nl[1],nl[2]);
  
  SpectralKernelGPU<KernelFn><<<grid,block>>>(wave, nx3, nl3, args...);
  ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
  ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);
  
  PetscFunctionReturn(ierr);
}

namespace reg {
namespace DifferentiationKernel {
  
PetscErrorCode VectorField::Laplacian(ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCallGPU<NLaplacianKernelGPU<1> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], b0*scale); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Laplacian(ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCallGPU<NLaplacianRegularizationKernelGPU<1> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], b0*scale, b1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Bilaplacian(ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCallGPU<NLaplacianKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], b0*scale); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Bilaplacian(ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCallGPU<NLaplacianRegularizationKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], b0*scale, b1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseBilaplacian(ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCallGPU<NInvLaplacianKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], scale, b0); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseBilaplacianSqrt(ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCallGPU<NInvLaplacianSqrtKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], scale, b0); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseBilaplacian(ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCallGPU<NInvLaplacianRegularizationKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], scale, b0, b1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseBilaplacianSqrt(ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCallGPU<NInvLaplacianRegularizationSqrtKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], scale, b0, b1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

} // namespace DifferentiationKernel
} // namespace reg

#endif
