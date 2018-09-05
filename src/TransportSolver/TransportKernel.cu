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

#include "TransportKernel.hpp"
#include "cuda_helper.hpp"

__global__ void TransformKernelAdjointSLGPU(ScalarType *pL, ScalarType* pLnext,
    ScalarType *pLx, ScalarType *pDivV, ScalarType *pDivVx, 
    ScalarType *pVec1, ScalarType *pVec2, ScalarType *pVec3,
    ScalarType *pB1, ScalarType *pB2, ScalarType *pB3, 
    ScalarType scale, ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    ScalarType lambda  = pL[i];
    ScalarType lambdax = pLx[i];

    ScalarType rhs0 = lambdax*pDivVx[i];
    ScalarType rhs1 = (lambdax + ht*rhs0)*pDivV[i];

    // compute \lambda(x,t^{j+1})
    pLnext[i] = lambdax + 0.5*ht*(rhs0 + rhs1);

    // compute bodyforce
    pB1[i] += scale*pVec1[i]*lambda;
    pB2[i] += scale*pVec2[i]*lambda;
    pB3[i] += scale*pVec3[i]*lambda;
  }
}
__global__ void TransformKernelAdjointSLGPU(ScalarType *pL,
    ScalarType *pVec1, ScalarType *pVec2, ScalarType *pVec3,
    ScalarType *pB1, ScalarType *pB2, ScalarType *pB3, 
    ScalarType scale, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    ScalarType lambda  = pL[i];

    // compute bodyforce
    pB1[i] += scale*pVec1[i]*lambda;
    pB2[i] += scale*pVec2[i]*lambda;
    pB3[i] += scale*pVec3[i]*lambda;
  }
}
__global__ void TransformKernelIncStateSLGPU(ScalarType *pM,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    ScalarType *pG1, ScalarType *pG2, ScalarType *pG3, 
    ScalarType hthalf, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    pM[i] -= hthalf*(pG1[i]*pV1[i] + pG2[i]*pV2[i] + pG3[i]*pV3[i]);
  }
}
#define TransformKernelIncAdjointGNGPU TransformKernelAdjointSLGPU

namespace reg {
  
PetscErrorCode TransportKernelAdjointSL::ComputeBodyForcePart1() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;

  TransformKernelAdjointSLGPU<<<grid, block>>>(pL, pLnext, pLx, pDivV, pDivVx, 
    pGm[0], pGm[1], pGm [2], 
    pB[0], pB[1], pB[2], 
    scale, ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointSL::ComputeBodyForcePart2() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;

  TransformKernelAdjointSLGPU<<<grid, block>>>(pL, 
    pGm[0], pGm[1], pGm[2], 
    pB[0], pB[1], pB[2], 
    scale, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateSL::TimeIntegrationPart1() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelIncStateSLGPU<<<grid, block>>>(pMtilde, 
    pVtildex[0], pVtildex[1], pVtildex[2], 
    pGm[0], pGm[1], pGm[2], 
    hthalf, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateSL::TimeIntegrationPart2() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelIncStateSLGPU<<<grid, block>>>(pMtilde, 
    pVtilde[0], pVtilde[1], pVtilde[2], 
    pGm[0], pGm[1], pGm[2], 
    hthalf, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointGN::ComputeBodyForce() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelIncAdjointGNGPU<<<grid, block>>>(pLtilde, 
      pGm[0], pGm[1], pGm[2], 
      pBtilde[0], pBtilde[1], pBtilde[2], 
      scale, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

} // namespace reg

