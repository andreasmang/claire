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

__global__ void VecSubMulGPU(ScalarType *pL, const ScalarType *pW, const ScalarType *pWts,
    const ScalarType *pMr, const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  if (i < nl) {
    pL[k] = pW[i] * pWts[blockIdx.y] * (pMr[k] - pM[k]);
  }
}

__global__ void VecSubGPU(ScalarType *pL, const ScalarType *pWts, const ScalarType *pMr, 
    const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  if (i < nl) {
    pL[k] = pWts[blockIdx.y] * (pMr[k] - pM[k]);
  }
}

__global__ void VecMulGPU(ScalarType *pL, const ScalarType *pW, const ScalarType *pWts,
    const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  if (i < nl) {
    pL[k] = - pW[i] * pWts[blockIdx.y] * pM[k];
  }
}

__global__ void VecNegGPU(ScalarType *pL, const ScalarType *pWts, const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  if (i < nl) {
    pL[k] = - pWts[blockIdx.y] * pM[k];
  }
}

// TODO I would rather use cuBlas functions for reduction for performace reasons
template<int N> inline __device__ void LocalReductionSum(ScalarType *shared) {
  if (threadIdx.x < N) {
    shared[threadIdx.x] += shared[threadIdx.x + N];
  }
  __syncthreads();
  LocalReductionSum<N/2>(shared);
}

template<> inline __device__ void LocalReductionSum<1>(ScalarType *shared) {
  if (threadIdx.x == 0) {
    shared[0] += shared[1];
  }
  __syncthreads();
}

template<int N> __global__ void ReductionSum(ScalarType *res, int n) {
  __shared__ ScalarType value[N];
  
  value[threadIdx.x] = 0.0;
  for (int i=threadIdx.x; i<n; i+=N) {
    value[threadIdx.x] += res[i];
  }
  
  __syncthreads();
  
  LocalReductionSum<N/2>(value);
  
  if (threadIdx.x == 0) {
    res[0] = value[0];
  }
}

template<int N>
__global__ void DistanceMeasureFunctionalGPU(ScalarType *res, 
    const ScalarType *pW, const ScalarType *pWts, 
    const ScalarType *pMr,const ScalarType *pM, 
    int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  __shared__ ScalarType value[N];
  
  ScalarType tmp;
  
  value[threadIdx.x] = 0.0;
  
  if (i < nl) {
    tmp = pMr[k] - pM[k];
    value[threadIdx.x] = tmp*tmp*pW[i]*pWts[blockIdx.y];
  }
  
  __syncthreads();
    
  LocalReductionSum<N/2>(value);
  
  if (threadIdx.x == 0) {
    res[blockIdx.x + blockIdx.y*gridDim.x] = value[0];
  }
}

template<int N>
__global__ void DistanceMeasureFunctionalGPU(ScalarType *res, 
    const ScalarType *pWts, const ScalarType *pMr, const ScalarType *pM, int nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int k = blockIdx.y*nl + i;
  
  __shared__ ScalarType value[N];
  ScalarType tmp;
  
  value[threadIdx.x] = 0.0;  
  if (i < nl) {
    tmp = (pMr[k] - pM[k]);
    value[threadIdx.x] = tmp*tmp*pWts[blockIdx.y];
  }
  
  __syncthreads();
    
  LocalReductionSum<N/2>(value);
  
  if (threadIdx.x == 0) {
    res[blockIdx.x + blockIdx.y*gridDim.x] = value[0];
  }
}

namespace reg {
namespace DistanceMeasureKernel {

/* Compute Masked Registration Functional */
PetscErrorCode EvaluateFunctionalSL2::ComputeFunctionalMask() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  ScalarType *res = nullptr;
  PetscFunctionBegin;
  
  ierr = AllocateMemoryOnce(res, grid.x*grid.y*sizeof(ScalarType)); CHKERRQ(ierr);
  
  DistanceMeasureFunctionalGPU<256><<<grid, block>>>(res, pW, pWts, pMr, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  ReductionSum<1024><<<1, 1024>>>(res, grid.x*grid.y);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  ierr = cudaMemcpy(reinterpret_cast<void*>(&value), reinterpret_cast<void*>(res), sizeof(ScalarType), cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
    
  FreeMemory(res);
  
  PetscFunctionReturn(ierr);
}

/* Compute the Registration Functional */
PetscErrorCode EvaluateFunctionalSL2::ComputeFunctional() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  ScalarType *res = nullptr;
  PetscFunctionBegin;
  
  ierr = AllocateMemoryOnce(res, grid.x*grid.y*sizeof(ScalarType)); CHKERRQ(ierr);
  
  DistanceMeasureFunctionalGPU<256><<<grid, block>>>(res, pWts, pMr, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  ReductionSum<1024><<<1, 1024>>>(res, grid.x*grid.y);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  ierr = cudaMemcpy(reinterpret_cast<void*>(&value), reinterpret_cast<void*>(res), sizeof(ScalarType), cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
    
  FreeMemory(res);
  
  PetscFunctionReturn(ierr);
}

/* Final Condition for Adjoint Equation */
PetscErrorCode FinalConditionSL2::ComputeFinalConditionAE() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  PetscFunctionBegin;
  
  VecSubGPU<<<grid, block>>>(pL, pWts, pMr, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

/* Final Condition for Masked Adjoint Equation */
PetscErrorCode FinalConditionSL2::ComputeFinalConditionMaskAE() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  PetscFunctionBegin;
  
  VecSubMulGPU<<<grid, block>>>(pL, pW, pWts, pMr, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

/* Final Condition for Incremental Adjoint Equation */
PetscErrorCode FinalConditionSL2::ComputeFinalConditionIAE() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  PetscFunctionBegin;
  
  VecNegGPU<<<grid, block>>>(pL, pWts, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

/* Final Condition for Masked Incremental Adjoint Equation */
PetscErrorCode FinalConditionSL2::ComputeFinalConditionMaskIAE() {
  PetscErrorCode ierr = 0;
  dim3 block(256, 1, 1);
  dim3 grid((nl + 255)/256, nc, 1);
  PetscFunctionBegin;
 
  VecMulGPU<<<grid, block>>>(pL, pW, pWts, pM, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

} // namespace DistanceMeasureKernel
} // namespace reg
