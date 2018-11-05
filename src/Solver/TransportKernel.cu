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
__global__ void TransformKernelAdjointSLGPU(const ScalarType *pL,
    const ScalarType *pVec1, const ScalarType *pVec2, const ScalarType *pVec3,
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
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    pM[i] -= ht*(pG1[i]*pV1[i] + pG2[i]*pV2[i] + pG3[i]*pV3[i]);
  }
}
#define TransformKernelIncAdjointGPU TransformKernelAdjointSLGPU

__global__ void TransformKernelEulerGPU(const ScalarType *pM,
    const ScalarType *pG1, const ScalarType *pG2, const ScalarType *pG3,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    ScalarType *pMbar, ScalarType *pRHS, 
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  ScalarType rhs;
  
  if (i < nl) {
    rhs = - pG1[i]*pV1[i] - pG2[i]*pV2[i] - pG3[i]*pV2[i];
    pRHS[i] = rhs;
    pMbar[i] = pM[i] +  ht*rhs;
  }
}
__global__ void TransformKernelEulerGPU(const ScalarType *pM,
    const ScalarType *pG1, const ScalarType *pG2, const ScalarType *pG3,
    const ScalarType *pGt1, const ScalarType *pGt2, const ScalarType *pGt3,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    const ScalarType *pVt1, const ScalarType *pVt2, const ScalarType *pVt3,
    ScalarType *pMbar, ScalarType *pRHS, 
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  ScalarType rhs;
  
  if (i < nl) {
    rhs = - pGt1[i]*pV1[i] - pGt2[i]*pV2[i] - pGt3[i]*pV2[i]
          - pG1[i]*pVt1[i] - pG2[i]*pVt2[i] - pG3[i]*pVt2[i];
    pRHS[i] = rhs;
    pMbar[i] = pM[i] +  ht*rhs;
  }
}
__global__ void TransformKernelRK2GPU(const ScalarType *pM, const ScalarType *pRHS,
    const ScalarType *pG1, const ScalarType *pG2, const ScalarType *pG3,
    const ScalarType *pGt1, const ScalarType *pGt2, const ScalarType *pGt3,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    const ScalarType *pVt1, const ScalarType *pVt2, const ScalarType *pVt3,
    ScalarType *pMbar,
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  ScalarType rhs;
  
  if (i < nl) {
    rhs = - pGt1[i]*pV1[i] - pGt2[i]*pV2[i] - pGt3[i]*pV2[i]
          - pG1[i]*pVt1[i] - pG2[i]*pVt2[i] - pG3[i]*pVt2[i];
    pMbar[i] = pM[i] +  ht*(rhs + pRHS[i]);
  }
}
__global__ void TransformKernelEulerGPU(const ScalarType *pM,
    const ScalarType *pG1, const ScalarType *pG2, const ScalarType *pG3,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    ScalarType *pMbar,
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  ScalarType rhs;
  
  if (i < nl) {
    rhs = - pG1[i]*pV1[i] - pG2[i]*pV2[i] - pG3[i]*pV2[i];
    pMbar[i] = pM[i] +  ht*rhs;
  }
}
__global__ void TransformKernelEulerGPU(const ScalarType *pM, const ScalarType *pRHS,
    ScalarType *pMbar, 
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (i < nl) {
    pMbar[i] = pM[i] +  ht*pRHS[i];
  }
}
__global__ void TransformKernelRK2GPU(const ScalarType *pM, const ScalarType *pRHS, 
    const ScalarType *pG1, const ScalarType *pG2, const ScalarType *pG3,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    ScalarType *pMnext, 
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  ScalarType rhs;
  
  if (i < nl) {
    rhs = - pG1[i]*pV1[i] - pG2[i]*pV2[i] - pG3[i]*pV2[i];
    pMnext[i] = pM[i] +  ht*(rhs + pRHS[i]);
  }
}
__global__ void TransformKernelRK2GPU(const ScalarType *pM, const ScalarType *pRHS1, const ScalarType *pRHS2,
    ScalarType *pMbar, 
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (i < nl) {
    pMbar[i] = pM[i] +  ht*(pRHS1[i] + pRHS2[i]); 
  }
}
__global__ void TransformKernelScaleGPU(const ScalarType *pL,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    ScalarType *pB1, ScalarType *pB2, ScalarType *pB3, 
    IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    ScalarType lambda  = pL[i];

    pB1[i] = pV1[i]*lambda;
    pB2[i] = pV2[i]*lambda;
    pB3[i] = pV3[i]*lambda;
  }
}
__global__ void TransformKernelScaleGPU(const ScalarType *pL, const ScalarType *pLt,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    const ScalarType *pVt1, const ScalarType *pVt2, const ScalarType *pVt3,
    ScalarType *pB1, ScalarType *pB2, ScalarType *pB3, 
    IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    ScalarType lambda  = pL[i];
    ScalarType lambdat  = pLt[i];

    pB1[i] = pV1[i]*lambdat + pVt1[i]*lambda;
    pB2[i] = pV2[i]*lambdat + pVt2[i]*lambda;
    pB3[i] = pV3[i]*lambdat + pVt3[i]*lambda;
  }
}
__global__ void TransformKernelScaleGPU(const ScalarType *pL, const ScalarType *pLt, const ScalarType* pRHS,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    const ScalarType *pVt1, const ScalarType *pVt2, const ScalarType *pVt3,
    ScalarType *pB1, ScalarType *pB2, ScalarType *pB3, 
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    ScalarType lambdat = pLt[i] + ht * pRHS[i];
    ScalarType lambda  = pL[i];

    pB1[i] = pV1[i]*lambdat + pVt1[i]*lambda;
    pB2[i] = pV2[i]*lambdat + pVt2[i]*lambda;
    pB3[i] = pV3[i]*lambdat + pVt3[i]*lambda;
  }
}
__global__ void TransformKernelScaleEulerGPU(const ScalarType *pL, const ScalarType *pRHS,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    ScalarType *pB1, ScalarType *pB2, ScalarType *pB3, 
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    ScalarType lambda  = pL[i] + ht*pRHS[i];

    pB1[i] = pV1[i]*lambda;
    pB2[i] = pV2[i]*lambda;
    pB3[i] = pV3[i]*lambda;
  }
}
__global__ void TransformKernelScaleRK2GPU(const ScalarType *pL, 
    const ScalarType *pRHS1, const ScalarType * pRHS2,
    const ScalarType *pV1, const ScalarType *pV2, const ScalarType *pV3,
    ScalarType *pB1, ScalarType *pB2, ScalarType *pB3, 
    ScalarType *pLnext,
    ScalarType ht, ScalarType scale, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    ScalarType lambda = pL[i];
    pLnext[i] = lambda + ht*(pRHS1[i]+ pRHS2[i]);

    pB1[i] += scale*pV1[i]*lambda;
    pB2[i] += scale*pV2[i]*lambda;
    pB3[i] += scale*pV3[i]*lambda;
  }
}
#define TransformKernelScaleAddGPU TransformKernelAdjointSLGPU
__global__ void TransformKernelContinuityGPU(const ScalarType *pMx,
    const ScalarType *pDivV, const ScalarType *pDivVx, ScalarType *pMnext,
    ScalarType ht, IntType nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    ScalarType mx = pMx[i];
    ScalarType rhs0 = - mx * pDivVx[i];
    ScalarType rhs1 = - (mx + ht*rhs0)*pDivV[i];
    
    pMnext[i] = mx + 0.5*ht*(rhs0 + rhs1);
  }
}

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

PetscErrorCode TransportKernelAdjoint::ComputeBodyForce() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelIncAdjointGPU<<<grid, block>>>(pL, 
      pGm[0], pGm[1], pGm[2], 
      pB[0], pB[1], pB[2], 
      scale, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelStateRK2::TimeIntegrationPart1() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelEulerGPU<<<grid, block>>>(pM,
    pGmx[0], pGmx[1], pGmx[2], 
    pVx[0], pVx[1], pVx[2],
    pMbar, pRHS,
    ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelStateRK2::TimeIntegrationPart2() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelRK2GPU<<<grid, block>>>(pM, pRHS,
    pGmx[0], pGmx[1], pGmx[2], 
    pVx[0], pVx[1], pVx[2],
    pMnext, 
    0.5*ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointRK2::TimeIntegrationPart1() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelScaleGPU<<<grid, block>>>(pL,
    pV[0], pV[1], pV[2],
    pVec[0], pVec[1], pVec[2],
    nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointRK2::TimeIntegrationPart2() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelScaleEulerGPU<<<grid, block>>>(pL, pRHS[0],
    pV[0], pV[1], pV[2],
    pVec[0], pVec[1], pVec[2],
    ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointRK2::TimeIntegrationPart3() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelScaleRK2GPU<<<grid, block>>>(pL, pRHS[0], pRHS[1],
    pVec[0], pVec[1], pVec[2],
    pB[0], pB[1], pB[2],
    pLnext,
    0.5*ht, scale, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointRK2::TimeIntegrationPart4() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelScaleAddGPU<<<grid, block>>>(pL,
    pVec[0], pVec[1], pVec[2],
    pB[0], pB[1], pB[2],
    0.5*scale, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart1a() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelScaleGPU<<<grid, block>>>(pL,
    pVtx[0], pVtx[1], pVtx[2],
    pLtjVx[0], pLtjVx[1], pLtjVx[2],
    nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart1b() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelScaleGPU<<<grid, block>>>(pL, pLt,
    pVx[0], pVx[1], pVx[2],
    pVtx[0], pVtx[1], pVtx[2],
    pLtjVx[0], pLtjVx[1], pLtjVx[2],
    nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart2a() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelEulerGPU<<<grid, block>>>(pLt, pRHS[0], pLtnext, ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart2b() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelScaleGPU<<<grid, block>>>(pL, pLt, pRHS[0],
    pVx[0], pVx[1], pVx[2],
    pVtx[0], pVtx[1], pVtx[2],
    pLtjVx[0], pLtjVx[1], pLtjVx[2],
    ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart3b() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelRK2GPU<<<grid, block>>>(pLt, pRHS[0], pRHS[1], pLtnext, 0.5*ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateRK2::TimeIntegrationEuler() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelEulerGPU<<<grid, block>>>(pMt,
    pGmx[0], pGmx[1], pGmx[2],
    pVtx[0], pVtx[1], pVtx[2],
    pMtnext,
    ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateRK2::TimeIntegrationPart1() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelEulerGPU<<<grid, block>>>(pMt,
    pGmx[0], pGmx[1], pGmx[2],
    pGmtx[0], pGmtx[1], pGmtx[2],
    pVx[0], pVx[1], pVx[2],
    pVtx[0], pVtx[1], pVtx[2],
    pMtbar, pRHS,
    ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateRK2::TimeIntegrationPart2() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelRK2GPU<<<grid, block>>>(pMt, pRHS,
    pGmx[0], pGmx[1], pGmx[2],
    pGmtx[0], pGmtx[1], pGmtx[2],
    pVx[0], pVx[1], pVx[2],
    pVtx[0], pVtx[1], pVtx[2],
    pMtnext,
    0.5*ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelContinuity::TimeIntegration() {
  PetscErrorCode ierr = 0;
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  PetscFunctionBegin;
  
  TransformKernelContinuityGPU<<<grid, block>>>(pMx, pDivV, pDivVx, pMnext, ht, nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();

  PetscFunctionReturn(ierr);
}

template<typename T>
PetscErrorCode TransportKernelCopy(T* org, T* dest, IntType ne) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  cudaMemcpy((void*)dest,(void*)org,sizeof(T)*ne,cudaMemcpyDeviceToDevice);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}
template PetscErrorCode TransportKernelCopy<ScalarType>(ScalarType*, ScalarType*, IntType);


} // namespace reg

