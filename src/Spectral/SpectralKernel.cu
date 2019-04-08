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
#include <cuda.h>
#include <cuda_runtime.h>

#include "SpectralKernel.txx"

using KernelUtils::SpectralKernelCallGPU;

namespace reg {
  
PetscErrorCode SpectralKernel::LowPassFilter(ComplexType *pXHat, ScalarType pct) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ScalarType l1, l2, l3;
  l1 = static_cast<ScalarType>(nx[0])*0.5*pct;
  l2 = static_cast<ScalarType>(nx[1])*0.5*pct;
  l3 = static_cast<ScalarType>(nx[2])*0.5*pct;
  ierr = SpectralKernelCallGPU<LowPassFilterKernel>(nstart, nx, nl, pXHat, l1, l2, l3, scale); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode SpectralKernel::HighPassFilter(ComplexType *pXHat, ScalarType pct) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ScalarType l1, l2, l3;
  l1 = static_cast<ScalarType>(nx[0])*0.5*pct;
  l2 = static_cast<ScalarType>(nx[1])*0.5*pct;
  l3 = static_cast<ScalarType>(nx[2])*0.5*pct;
  ierr = SpectralKernelCallGPU<HighPassFilterKernel>(nstart, nx, nl, pXHat, l1, l2, l3, scale); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode SpectralKernel::Scale(ComplexType *pX, ScalarType val) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = SpectralKernelCallGPU<ScaleKernel>(nstart, nx, nl, pX, val); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

__global__ void FilterKernel(int3 wave, int3 nl, ComplexType *x, int3 nxc) {
  int i3 = threadIdx.x + blockIdx.x*blockDim.x;
  int i2 = blockIdx.y;
  int i1 = blockIdx.z;
  
  if (i3 < nl.z) {
    wave.x += i1;
    wave.y += i2;
    wave.z += i3;

    int i = i1*nl.y + i2*nl.z + i3;
    
    if (wave.x > nxc.x || wave.y > nxc.y || wave.z > nxc.z) {
      x[i][0] = 0.;
      x[i][1] = 0.;
    }
  }
}

PetscErrorCode SpectralKernel::Restrict(ComplexType *pXc, const ComplexType *pXf, const IntType nxc[3]) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  cudaMemset(pXc, 0, sizeof(ComplexType)*nxc[0]*nxc[1]*nxc[2]);
  /*for (int x = 0; x < nxc[0]; ++x) {
    for (int y = 0; y < nxc[1]; ++y) {
      int idx_c, idx_f;
      idx_c = y*nxc[2] + x*nxc[2]*nxc[1];
      idx_f = y*nl[2] + x*nl[2]*nl[1];
      if (y > nxc[1]/2)
        idx_f += (nl[1] - nxc[1])*nl[2];
      if (x > nxc[0]/2)
        idx_f += (nl[0] - nxc[0])*nl[2]*nl[1];
      if (x != nxc[0]/2 && y != nxc[1]/2)
        cudaMemcpyAsync(&pXc[idx_c], &pXf[idx_f], sizeof(ComplexType)*(nxc[2]-1), cudaMemcpyDeviceToDevice);
    }
  }*/
  
  cudaMemcpy3DParms params = {0};
  params.kind = cudaMemcpyDeviceToDevice;
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(reinterpret_cast<const void*>(pXf)), nl[2]*sizeof(ComplexType), nl[2], nl[1]);
  params.dstPtr = make_cudaPitchedPtr(reinterpret_cast<void*>(pXc), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.extent =  make_cudaExtent(nxc[2]-1, nxc[1]/2, nxc[0]/2);
  cudaMemcpy3D(&params);
  size_t offset_f, offset_c;
  offset_f = (nl[1]-(nxc[1]/2) + 1)*nl[2];
  offset_c = (nxc[1]/2+1)*nxc[2];
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(reinterpret_cast<const void*>(&pXf[offset_f])), nl[2]*sizeof(ComplexType), nl[2], nl[1]);
  params.dstPtr = make_cudaPitchedPtr(reinterpret_cast<void*>(&pXc[offset_c]), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.extent =  make_cudaExtent(nxc[2]-1, nxc[1]/2-1, nxc[0]/2-1);
  cudaMemcpy3D(&params);
  offset_f = (nl[0]-(nxc[0]/2) + 1)*nl[2]*nl[1];
  offset_c = (nxc[0]/2+1)*nxc[2]*nxc[1];
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(reinterpret_cast<const void*>(&pXf[offset_f])), nl[2]*sizeof(ComplexType), nl[2], nl[1]);
  params.dstPtr = make_cudaPitchedPtr(reinterpret_cast<void*>(&pXc[offset_c]), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.extent =  make_cudaExtent(nxc[2]-1, nxc[1]/2, nxc[0]/2-1);
  cudaMemcpy3D(&params);
  offset_f = (nl[1]-(nxc[1]/2) + 1)*nl[2] + (nl[0]-(nxc[0]/2) + 1)*nl[2]*nl[1];
  offset_c = (nxc[1]/2+1)*nxc[2] + (nxc[0]/2+1)*nxc[2]*nxc[1];
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(reinterpret_cast<const void*>(&pXf[offset_f])), nl[2]*sizeof(ComplexType), nl[2], nl[1]);
  params.dstPtr = make_cudaPitchedPtr(reinterpret_cast<void*>(&pXc[offset_c]), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.extent =  make_cudaExtent(nxc[2]-1, nxc[1]/2-1, nxc[0]/2-1);
  cudaMemcpy3D(&params);
  
  cudaDeviceSynchronize();
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode SpectralKernel::Prolong(ComplexType *pXf, const ComplexType *pXc, const IntType nxc[3]) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  cudaMemset(pXf, 0, sizeof(ComplexType)*nl[0]*nl[1]*nl[2]);
  /*cudaMemcpy3DParms params = {0};
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(static_cast<const void*>(pXc)), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.dstPtr = make_cudaPitchedPtr(static_cast<void*>(pXf), nx[2]*sizeof(ComplexType), nx[2], nx[1]);
  params.extent =  make_cudaExtent(nxc[2], nxc[1], nxc[0]);
  params.kind = cudaMemcpyDeviceToDevice;
  cudaMemcpy3D(&params);*/
  /*
  for (int x = 0; x < nxc[0]; ++x) {
    for (int y = 0; y < nxc[1]; ++y) {
      int idx_c, idx_f;
      idx_c = y * nxc[2] + x * nxc[2]*nxc[1];
      idx_f = y * nl[2] + x * nl[2]*nl[1];
      if (y > nxc[1]/2)
        idx_f += (nl[1] - nxc[1])*nl[2];
      if (x > nxc[0]/2)
        idx_f += (nl[0] - nxc[0])*nl[2]*nl[1];
      if (x != nxc[0]/2 && y != nxc[1]/2)
      cudaMemcpyAsync(&pXf[idx_f], &pXc[idx_c], sizeof(ComplexType)*(nxc[2]-1), cudaMemcpyDeviceToDevice);
    }
  }*/
  
  cudaMemcpy3DParms params = {0};
  params.kind = cudaMemcpyDeviceToDevice;
  params.dstPtr = make_cudaPitchedPtr(reinterpret_cast<void*>(pXf), nl[2]*sizeof(ComplexType), nl[2], nl[1]);
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(reinterpret_cast<const void*>(pXc)), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.extent =  make_cudaExtent(nxc[2]-1, nxc[1]/2, nxc[0]/2);
  cudaMemcpy3D(&params);
  size_t offset_f, offset_c;
  offset_f = (nl[1]-(nxc[1]/2) + 1)*nl[2];
  offset_c = (nxc[1]/2+1)*nxc[2];
  params.dstPtr = make_cudaPitchedPtr(reinterpret_cast<void*>(&pXf[offset_f]), nl[2]*sizeof(ComplexType), nl[2], nl[1]);
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(reinterpret_cast<const void*>(&pXc[offset_c])), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.extent =  make_cudaExtent(nxc[2]-1, nxc[1]/2-1, nxc[0]/2-1);
  cudaMemcpy3D(&params);
  offset_f = (nl[0]-(nxc[0]/2) + 1)*nl[2]*nl[1];
  offset_c = (nxc[0]/2+1)*nxc[2]*nxc[1];
  params.dstPtr = make_cudaPitchedPtr(reinterpret_cast<void*>(&pXf[offset_f]), nl[2]*sizeof(ComplexType), nl[2], nl[1]);
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(reinterpret_cast<const void*>(&pXc[offset_c])), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.extent =  make_cudaExtent(nxc[2]-1, nxc[1]/2, nxc[0]/2-1);
  cudaMemcpy3D(&params);
  offset_f = (nl[1]-(nxc[1]/2) + 1)*nl[2] + (nl[0]-(nxc[0]/2) + 1)*nl[2]*nl[1];
  offset_c = (nxc[1]/2+1)*nxc[2] + (nxc[0]/2+1)*nxc[2]*nxc[1];
  params.dstPtr = make_cudaPitchedPtr(reinterpret_cast<void*>(&pXf[offset_f]), nl[2]*sizeof(ComplexType), nl[2], nl[1]);
  params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(reinterpret_cast<const void*>(&pXc[offset_c])), nxc[2]*sizeof(ComplexType), nxc[2], nxc[1]);
  params.extent =  make_cudaExtent(nxc[2]-1, nxc[1]/2-1, nxc[0]/2-1);
  cudaMemcpy3D(&params);
  
  cudaDeviceSynchronize();
  
  /*dim3 block, grid;
  if (nl[2] <= 1024 && nl[2] >= 32) {
    block.x = (nl[2] + 31)/32;
    block.x *= 32;
    grid.x = 1;
  } else {
    block.x = 128; // 128 threads per block
    grid.x = (nl[2] + 127)/128;  // $\lceil nl_2 / 128 \rceil$
  }
  grid.y = nl[1];
  grid.z = nl[0];
  int3 nl3, wave, nxc3;
  wave.x = nstart[0]; wave.y = nstart[1]; wave.z = nstart[2];
  nl3.x = nl[0]; nl3.y = nl[1]*nl[2]; nl3.z = nl[2];
  nxc3.x = nxc[0]; nxc3.y = nxc[1]; nxc3.z = nxc[2];
  
  if (nl[0]*nl[1]*nl[2] > 0) {
    FilterKernel<<<grid, block>>>(wave, nl3, pXf, nxc3);
    ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
    ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);
  }*/
  
  PetscFunctionReturn(ierr);
}

} // namepsace reg

#endif
