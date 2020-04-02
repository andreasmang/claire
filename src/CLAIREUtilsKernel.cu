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

#ifndef _CLAIREUTILSKERNEL_CU_
#define _CLAIREUTILSKERNEL_CU_

#include "CLAIREUtils.hpp"
#include "cuda_helper.hpp"

// CUDA kernel to evaluate point-wise norm of a vector field
__global__ void VecFieldPointWiseNormKernel(ScalarType *p_m, const ScalarType *p_X1, const ScalarType *p_X2, const ScalarType *p_X3, IntType nl) {
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    if (i < nl) {
        p_m[i] = sqrtf(p_X1[i]*p_X1[i] + p_X2[i]*p_X2[i] + p_X3[i]*p_X3[i]);
    }
}

__global__ void CopyStridedToFlatVecKernel(ScalarType *pX, const ScalarType *p_x1, const ScalarType *p_x2, const ScalarType *p_x3, IntType nl) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    if (i < nl) {
        pX[3*i + 0] = p_x1[i];
        pX[3*i + 1] = p_x2[i];
        pX[3*i + 2] = p_x3[i];
    }

}

__global__ void CopyStridedFromFlatVecKernel(ScalarType *p_x1, ScalarType *p_x2, ScalarType *p_x3, const ScalarType* pX, IntType nl) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    if (i < nl) {
        p_x1[i] = pX[3*i + 0];
        p_x2[i] = pX[3*i + 1];
        p_x3[i] = pX[3*i + 2];
    }
}

namespace reg {
  

/********************************************************************
 * @brief compute pointwise norm of vector field
 *******************************************************************/
PetscErrorCode VecFieldPointWiseNormGPU(ScalarType* p_m, const ScalarType* p_X1, const ScalarType* p_X2, const ScalarType* p_X3, IntType nl) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    dim3 block(256, 1, 1);
    dim3 grid((nl + 255)/256, 1, 1);
    
    VecFieldPointWiseNormKernel<<<grid, block>>>(p_m, p_X1, p_X2, p_X3, nl);
    cudaDeviceSynchronize();
    cudaCheckKernelError();

    PetscFunctionReturn(ierr);

}


/********************************************************************
 * @brief Copy vector field to a flat array in strided fashion
 *******************************************************************/
PetscErrorCode CopyStridedToFlatVec(ScalarType* pX, const ScalarType* p_x1, const ScalarType* p_x2, const ScalarType* p_x3, IntType nl) {
    PetscFunctionBegin;
    PetscErrorCode ierr = 0;
    
    int threads = 256;
    int blocks = (nl + 255)/threads;

    CopyStridedToFlatVecKernel<<<blocks,threads>>>(pX, p_x1, p_x2, p_x3, nl);
    cudaDeviceSynchronize();
    cudaCheckKernelError();

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Copy vector field to a flat array in strided fashion
 *******************************************************************/
PetscErrorCode CopyStridedFromFlatVec(ScalarType* p_x1, ScalarType* p_x2, ScalarType* p_x3, const ScalarType* pX, IntType nl) {
    PetscFunctionBegin;
    PetscErrorCode ierr = 0;
    
    int threads = 256;
    int blocks = (nl + 255)/threads;

    CopyStridedFromFlatVecKernel<<<blocks,threads>>>(p_x1, p_x2, p_x3, pX, nl);
    cudaDeviceSynchronize();
    cudaCheckKernelError();

    PetscFunctionReturn(ierr);
}

}

#endif
