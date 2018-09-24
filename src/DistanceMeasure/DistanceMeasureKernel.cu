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
#include "cuda_helper.hpp"

__global__ void VecSubMulGPU(ScalarType *pL, const ScalarType *pW, const ScalarType *pMr, 
    const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  if (i < nl) {
    pL[k] = pW[i] * (pMr[k] - pM[k]);
  }
}

__global__ void VecSubGPU(ScalarType *pL, const ScalarType *pMr, const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  if (i < nl) {
    pL[k] = pMr[k] - pM[k];
  }
}

__global__ void VecMulGPU(ScalarType *pL, const ScalarType *pW, const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  if (i < nl) {
    pL[k] = - pW[i] * pM[k];
  }
}

__global__ void VecNegGPU(ScalarType *pL, const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  if (i < nl) {
    pL[k] = - pM[k];
  }
}

namespace reg {

PetscErrorCode DistanceMeasureKernelSL2::ComputeFinalConditionAE() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  PetscFunctionBegin;
  
  VecSubGPU<<<grid, block>>>(pL, pMr, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode DistanceMeasureKernelSL2::ComputeFinalConditionMaskAE() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  PetscFunctionBegin;
  
  VecSubMulGPU<<<grid, block>>>(pL, pW, pMr, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode DistanceMeasureKernelSL2::ComputeFinalConditionIAE() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  PetscFunctionBegin;
  
  VecNegGPU<<<grid, block>>>(pL, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode DistanceMeasureKernelSL2::ComputeFinalConditionMaskIAE() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  PetscFunctionBegin;
  
  VecSubGPU<<<grid, block>>>(pL, pW, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

} // namespace reg
