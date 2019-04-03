/*************************************************************************
 *  Copyright (c) 2018.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you c1n redistribute it and/or modify
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


#ifndef _TEXTUREDIFFERENTIATIONKERNEL_CPP_
#define _TEXTUREDIFFERENTIATIONKERNEL_CPP_

#include "TextureDifferentiationKernel.hpp"
#include "cuda_helper.hpp"
#include "cuda_profiler_api.h"
#define HALO 4 
//#define HALO 2 
#define spencil 4
#define lpencil 32
#define sharedrows 64
#define perthreadcomp 8

const float h_c[HALO] = {4.f / 5.f , -1.f / 5.f , 4.f / 105.f, -1.f / 280.f};
//const float h_c[HALO] = {2.f / 3.f , -1.f / 12.f};


// device constants
__constant__ int d_nx, d_ny, d_nz;
__constant__ float d_invnx, d_invny, d_invnz;
__constant__ float d_invhx, d_invhy, d_invhz;
__constant__ float d_cx[HALO], d_cy[HALO], d_cz[HALO];

const int sx = spencil;
const int sy = sharedrows;
const int sxx = lpencil;
const int syy = sharedrows;


__device__ inline int getLinearIdx(int i, int j, int k) {
    return i*d_ny*d_nz + j*d_nz + k;
}


__global__ void gradient_z(ScalarType* dfz, const ScalarType* f) {
  __shared__ float s_f[sx][sy+2*HALO]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i   = blockIdx.z;
  int j   = blockIdx.y*blockDim.y + threadIdx.y;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x + HALO;       // local k for shared memory ac3ess + halo offset
  int sj = threadIdx.y; // local j for shared memory ac3ess

  int globalIdx = getLinearIdx(i,j,k);

  s_f[sj][sk] = f[globalIdx];

  __syncthreads();
    
  int lid,rid;
  // fill in periodic images in shared memory array 
  if (threadIdx.x < HALO) {
    lid = k%d_nz-HALO;
    if (lid<0) lid+=d_nz;
    s_f[sj][sk-HALO] = f[i*d_ny*d_nz + j*d_nz + lid];
    rid = (k+sy)%d_nz;
    s_f[sj][sy+sk] = f[i*d_ny*d_nz + j*d_nz + rid];
  }

  __syncthreads();
/* 
  dfz[globalIdx] = 
    (  c1 * ( s_f[sj][sk+1] - s_f[sj][sk-1] )
    +  c2 * ( s_f[sj][sk+2] - s_f[sj][sk-2] )
    +  c3 * ( s_f[sj][sk+3] - s_f[sj][sk-3] )
    +  c4 * ( s_f[sj][sk+4] - s_f[sj][sk-4] ) );
*/
    dfz[globalIdx] = 0;
    for(int l=0; l<HALO; l++) {
        dfz[globalIdx] += d_cz[l] * (s_f[sj][sk+1+l] - s_f[sj][sk-1-l]);
    }
}

__global__ void gradient_y(ScalarType* dfy, const ScalarType* f) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
    
    
  for(int j = threadIdx.y; j < syy; j += blockDim.y) {
    int globalIdx = getLinearIdx(i, blockIdx.y*syy + j ,k);
    int sj = j + HALO;
    s_f[sj][sk] = f[globalIdx];
  }

  __syncthreads();

  
  int lid,rid, sj = threadIdx.y + HALO;
  int y = syy*blockIdx.y + threadIdx.y;
  // fill in periodic images in shared memory array 
  if (threadIdx.y < HALO) {
    lid = y%d_ny-HALO;
    if (lid<0) lid+=d_ny;
    s_f[sj-HALO][sk] = f[i*d_ny*d_nz + lid*d_nz + k];
    rid = (y+syy)%d_ny;
    s_f[sj+syy][sk] = f[i*d_ny*d_nz + rid*d_nz + k];
  }

  __syncthreads();
  
  for(int j = threadIdx.y; j < syy; j += blockDim.y) {
    int globalIdx = getLinearIdx(i, blockIdx.y*syy + j ,k);
    int sj = j + HALO;
    dfy[globalIdx] = 0;
    for( int l=0; l<HALO; l++) {
        dfy[globalIdx] += d_cy[l] * ( s_f[sj+1+l][sk] - s_f[sj-1-l][sk]);
    }
    /*
  dfy[globalIdx] = 
    (  c1 * ( s_f[sj+1][sk] - s_f[sj-1][sk] )
    +  c2 * ( s_f[sj+2][sk] - s_f[sj-2][sk] )
    +  c3 * ( s_f[sj+3][sk] - s_f[sj-3][sk] )
    +  c4 * ( s_f[sj+4][sk] - s_f[sj-4][sk] ) );
    */
  }

}

__global__ void gradient_x(ScalarType* dfx, const ScalarType* f) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int j  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
    
    
  for(int i = threadIdx.y; i < syy; i += blockDim.y) {
    int globalIdx = getLinearIdx(blockIdx.y*syy + i, j ,k);
    int si = i + HALO;
    s_f[si][sk] = f[globalIdx];
  }

  __syncthreads();

  
  int lid,rid, si = threadIdx.y + HALO;
  int x = syy*blockIdx.y + threadIdx.y;
  // fill in periodic images in shared memory array 
  if (threadIdx.y < HALO) {
    lid = x%d_nx-HALO;
    if (lid<0) lid+=d_nx;
    s_f[si-HALO][sk] = f[lid*d_ny*d_nz + j*d_nz + k];
    rid = (x+syy)%d_nx;
    s_f[si+syy][sk] = f[rid*d_ny*d_nz + j*d_nz + k];
  }

  __syncthreads();
  
  for(int i = threadIdx.y; i < syy; i += blockDim.y) {
    int globalIdx = getLinearIdx(blockIdx.y*syy + i , j, k);
    int si = i + HALO;
    dfx[globalIdx] = 0;
    for( int l=0; l<HALO; l++) {
        dfx[globalIdx] += d_cx[l] * ( s_f[si+1+l][sk] - s_f[si-1-l][sk]);
    }
  /*
  dfx[globalIdx] = 
    (  c1 * ( s_f[si+1][sk] - s_f[si-1][sk] )
    +  c2 * ( s_f[si+2][sk] - s_f[si-2][sk] )
    +  c3 * ( s_f[si+3][sk] - s_f[si-3][sk] )
    +  c4 * ( s_f[si+4][sk] - s_f[si-4][sk] ) );
    */
  }

}

__global__ void TextureDivXComputeKernel(cudaTextureObject_t tex, ScalarType* div) {
    const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
    const int tidy = blockDim.y * blockIdx.y + threadIdx.y;
    const int tidz = blockDim.z * blockIdx.z + threadIdx.z;

    if (tidx < d_nx && tidy < d_ny && tidz < d_nz) {    
      // global index
      const int gid = tidz + tidy*d_nz + tidx*d_ny*d_nz;
      float3 id = make_float3( tidz*d_invnz, tidy*d_invny, tidx*d_invnx);
      
      float dfx=0;
      for(int l=1; l<HALO+1; l++) {
          dfx += (tex3D<float>(tex, id.x, id.y, id.z + l*d_invnx) - tex3D<float>(tex, id.x, id.y, id.z - l*d_invnx))*d_cx[l-1];
      }
      div[gid] = dfx;
    }
}

__global__ void TextureDivYComputeKernel(cudaTextureObject_t tex, ScalarType* div) {
    const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
    const int tidy = blockDim.y * blockIdx.y + threadIdx.y;
    const int tidz = blockDim.z * blockIdx.z + threadIdx.z;

    if (tidx < d_nx && tidy < d_ny && tidz < d_nz) {    
      // global index
      const int gid = tidz + tidy*d_nz + tidx*d_ny*d_nz;
      float3 id = make_float3( tidz*d_invnz, tidy*d_invny, tidx*d_invnx);
      
      float dfy=0;
      for(int l=1; l<HALO+1; l++) {
          dfy += (tex3D<float>(tex, id.x, id.y + l*d_invny, id.z) - tex3D<float>(tex, id.x, id.y - l*d_invny, id.z))*d_cy[l-1];
      }
      div[gid] += dfy;
    }
}

__global__ void TextureDivZComputeKernel(cudaTextureObject_t tex, ScalarType* div) {
    const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
    const int tidy = blockDim.y * blockIdx.y + threadIdx.y;
    const int tidz = blockDim.z * blockIdx.z + threadIdx.z;

    if (tidx < d_nx && tidy < d_ny && tidz < d_nz) {    
      // global index
      const int gid = tidz + tidy*d_nz + tidx*d_ny*d_nz;
      float3 id = make_float3( tidz*d_invnz, tidy*d_invny, tidx*d_invnx);
      
      float dfz=0;
      for(int l=1; l<HALO+1; l++) {
          dfz += (tex3D<float>(tex, id.x + l*d_invnz, id.y, id.z) - tex3D<float>(tex, id.x - l*d_invnz, id.y, id.z))*d_cz[l-1];
      }
      div[gid] += dfz;
    }
}


__global__ void TextureGradientComputeKernel(cudaTextureObject_t tex, ScalarType* dmx, ScalarType* dmy, ScalarType* dmz) {
    const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
    const int tidy = blockDim.y * blockIdx.y + threadIdx.y;
    const int tidz = blockDim.z * blockIdx.z + threadIdx.z;
    
    if (tidx < d_nx && tidy < d_ny && tidz < d_nz) {
      // global index
      const int gid = tidz + tidy*d_nz + tidx*d_ny*d_nz;
      float3 id = make_float3( tidz*d_invnz, tidy*d_invny, tidx*d_invnx);
    
      float dfx=0,dfy=0,dfz=0;
      for(int l=1; l<HALO+1; l++) {
          dfz += (tex3D<float>(tex, id.x + l*d_invnz, id.y, id.z) - tex3D<float>(tex, id.x - l*d_invnz, id.y, id.z))*d_cz[l-1];
          dfy += (tex3D<float>(tex, id.x, id.y + l*d_invny, id.z) - tex3D<float>(tex, id.x, id.y - l*d_invny, id.z))*d_cy[l-1];
          dfx += (tex3D<float>(tex, id.x, id.y, id.z + l*d_invnx) - tex3D<float>(tex, id.x, id.y, id.z - l*d_invnx))*d_cx[l-1];
      }
      dmz[gid] = dfz;
      dmy[gid] = dfy;
      dmx[gid] = dfx;
    }
}


 /*   
__global__ void TextureGradientComputeKernel_old(cudaTextureObject_t tex, ScalarType* dmx, ScalarType* dmy, ScalarType* dmz) {
    const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
    const int tidy = blockDim.y * blockIdx.y + threadIdx.y;
    const int tidz = blockDim.z * blockIdx.z + threadIdx.z;
    
    
    // global index
    const int gid = tidz + tidy*d_nz + tidx*d_ny*d_nz;
    //float3 id = make_float3( tidx*inv_nx.x, tidy*inv_nx.y, tidz*inv_nx.z);
    //float3 id = make_float3( tidx/nx.x, tidy/nx.y, tidz/nx.z);
    float3 id = make_float3( tidz*d_invnz, tidy*d_invny, tidx*d_invnx);
    //float3 id = make_float3( tidz*inv_nx.z, tidy*inv_nx.y, tidx*inv_nx.x);
    dmz[gid] = (tex3D<float>(tex, id.x + d_invnz, id.y, id.z) - tex3D<float>(tex, id.x - d_invnz, id.y, id.z))*d_cz[0];
    dmy[gid] = (tex3D<float>(tex, id.x, id.y + d_invny, id.z) - tex3D<float>(tex, id.x, id.y - d_invny, id.z))*d_cy[0];
    dmx[gid] = (tex3D<float>(tex, id.x, id.y, id.z + d_invnx) - tex3D<float>(tex, id.x, id.y, id.z - d_invnx))*d_cz[0];

    // print value
    //float3 hx = make_float3( 2*M_PI/(nx.x-1), 2*M_PI/(nx.y-1), 2*M_PI/(nx.z-1)); 
    float3 hx = make_float3( 2*M_PI/(nx.x), 2*M_PI/(nx.y), 2*M_PI/(nx.z)); 
    float x1 = hx.x*tidx - M_PI;
    float x2 = hx.y*tidy - M_PI;
    float x3 = hx.z*tidz - M_PI;
    float m,mtrue,gxtrue,gytrue,gztrue;
        
    mtrue = sinf(M_PI*sinf(x1)) * sinf(M_PI*sinf(x2)) * sinf(M_PI*sinf(x3));
    
    gxtrue = M_PI*cosf(x1)*cosf(M_PI*sinf(x1)) *
                              sinf(M_PI*sinf(x2)) *
                              sinf(M_PI*sinf(x3));
    gytrue = sinf(M_PI*sinf(x1)) *
                              M_PI*cosf(x2)*cosf(M_PI*sinf(x2)) *
                              sinf(M_PI*sinf(x3));
    gztrue = sinf(M_PI*sinf(x1)) *
                              sinf(M_PI*sinf(x2)) *
                              M_PI*cosf(x3)*cosf(M_PI*sinf(x3));
    //mtrue = sinf(M_PI*x1) + sinf(M_PI*x2) + sinf(M_PI*x3);
    //mtrue = gid;
    //gxtrue = 0;
    m = tex3D<float>(tex, id.x, id.y, id.z);
    
    if (gid >=0) {
        //printf("gid=%d,\t tidx=%d,\t tidy=%d,\t tidz=%d,\t x=%f,\t y=%f,\t z=%f,\t mtrue=%f,\t m=%f,\t ref=%f,\t gxt=%f,\t gxc=%f\n\n", gid, tidx, tidy, tidz, x1, x2,
        //x3, mtrue, m, ref[gid], gxtrue, gx[gid]);
    }

}

*/   

void printFloat3(float3 a){
    printf("x = %f\t y = %f\t z = %f\n",a.x,a.y,a.z);
}


namespace reg {

cudaTextureObject_t gpuInitEmptyGradientTexture(IntType *nx) {
   cudaError_t err = cudaSuccess;
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
   cudaExtent extent = make_cudaExtent(nx[2], nx[1], nx[0]);
   cudaArray* cuArray;
   err = cudaMalloc3DArray(&cuArray, &channelDesc, extent, 0);
   if (err != cudaSuccess){
        fprintf(stderr, "Failed to allocate 3D cudaArray (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
   }
    
    /* create cuda resource description */
    struct cudaResourceDesc resDesc;
    memset( &resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = cuArray;

    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.addressMode[0] = cudaAddressModeWrap;
    texDesc.addressMode[1] = cudaAddressModeWrap;
    texDesc.addressMode[2] = cudaAddressModeWrap;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.normalizedCoords = 1;

    cudaTextureObject_t texObj = 0;
    err = cudaCreateTextureObject( &texObj, &resDesc, &texDesc, NULL);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to create texture (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    return texObj;
}


/********************************************************************************
 * @brief update texture object by copying volume data to 3D cudaArray container
 *******************************************************************************/
void updateTextureFromVolume(cudaPitchedPtr volume, cudaExtent extent, cudaTextureObject_t texObj) {
    cudaError_t err = cudaSuccess;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

    /* create cuda resource description */
    struct cudaResourceDesc resDesc;
    memset( &resDesc, 0, sizeof(resDesc));
    cudaGetTextureObjectResourceDesc( &resDesc, texObj);

    cudaMemcpy3DParms p = {0};
    p.srcPtr = volume;
    p.dstArray = resDesc.res.array.array;
    p.extent = extent;
    p.kind = cudaMemcpyDeviceToDevice;
    err = cudaMemcpy3D(&p);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to copy 3D memory to cudaArray (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

}

PetscErrorCode initConstants(IntType* nx) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  const float3 inv_nx = make_float3(  1.0f/static_cast<float>(nx[0]),
                                        1.0f/static_cast<float>(nx[1]), 
                                        1.0f/static_cast<float>(nx[2]));
  const float3 hx = make_float3(2*M_PI/(nx[0]), 2*M_PI/(nx[1]), 2*M_PI/(nx[2]));
  const float3 inv_hx = make_float3(0.5f/hx.x, 0.5f/hx.y, 0.5f/hx.z);
  
  cudaMemcpyToSymbol(d_nx, &nx[0], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_ny, &nx[1], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_nz, &nx[2], sizeof(int), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(d_invnx, &inv_nx.x, sizeof(float), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_invny, &inv_nx.y, sizeof(float), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_invnz, &inv_nx.z, sizeof(float), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(d_invhx, &inv_hx.x, sizeof(float), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_invhy, &inv_hx.y, sizeof(float), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_invhz, &inv_hx.z, sizeof(float), 0, cudaMemcpyHostToDevice);
  
  float h_ct[HALO];
  for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx.x;
  cudaMemcpyToSymbol(d_cx, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);
  for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx.y;
  cudaMemcpyToSymbol(d_cy, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);
  for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx.z;
  cudaMemcpyToSymbol(d_cz, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode computeTextureDivergence(ScalarType* l, const ScalarType* g1, const ScalarType* g2, const ScalarType* g3, cudaTextureObject_t mtex, IntType* nx) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
   // create a cudaExtent for input resolution
  cudaExtent extent = make_cudaExtent(nx[2], nx[1], nx[0]);
  
  // Texture gradient
  dim3 threadsPerBlock(1,8,32);
  dim3 numBlocks(nx[0]/threadsPerBlock.x, (nx[1]+7)/threadsPerBlock.y, (nx[2]+31)/threadsPerBlock.z);
  
  cudaPitchedPtr m_cudaPitchedPtr;
  
  // make input image a cudaPitchedPtr for m
  m_cudaPitchedPtr = make_cudaPitchedPtr((void*)(g1), nx[2]*sizeof(ScalarType), nx[2], nx[1]);
  // update texture object
  updateTextureFromVolume(m_cudaPitchedPtr, extent, mtex);

  TextureDivXComputeKernel<<<numBlocks, threadsPerBlock>>>(mtex, l);
  cudaCheckKernelError();
  cudaDeviceSynchronize();
  
  // make input image a cudaPitchedPtr for m
  m_cudaPitchedPtr = make_cudaPitchedPtr((void*)(g2), nx[2]*sizeof(ScalarType), nx[2], nx[1]);
  // update texture object
  updateTextureFromVolume(m_cudaPitchedPtr, extent, mtex);

  TextureDivYComputeKernel<<<numBlocks, threadsPerBlock>>>(mtex, l);
  cudaCheckKernelError();
  cudaDeviceSynchronize();
  
  // make input image a cudaPitchedPtr for m
  m_cudaPitchedPtr = make_cudaPitchedPtr((void*)(g3), nx[2]*sizeof(ScalarType), nx[2], nx[1]);
  // update texture object
  updateTextureFromVolume(m_cudaPitchedPtr, extent, mtex);

  TextureDivZComputeKernel<<<numBlocks, threadsPerBlock>>>(mtex, l);
  cudaCheckKernelError();
  cudaDeviceSynchronize();
  
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode computeTextureGradient(ScalarType* gx, ScalarType* gy, ScalarType* gz, const ScalarType* m, cudaTextureObject_t mtex, IntType* nx) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    //float time=0, dummy_time=0;
    //int repcount = 1;
    //cudaEvent_t startEvent, stopEvent;
    //cudaEventCreate(&startEvent);
    //cudaEventCreate(&stopEvent);

    /*const float3 inv_nx = make_float3(  1.0f/static_cast<float>(nx[0]),
                                        1.0f/static_cast<float>(nx[1]), 
                                        1.0f/static_cast<float>(nx[2]));
    const float3 hx = make_float3(2*M_PI/(nx[0]), 2*M_PI/(nx[1]), 2*M_PI/(nx[2]));
    const float3 inv_hx = make_float3(0.5f/hx.x, 0.5f/hx.y, 0.5f/hx.z);
    long int nq = nx[0]*nx[1]*nx[2]; 

    
    cudaMemcpyToSymbol(d_nx, &nx[0], sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_ny, &nx[1], sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_nz, &nx[2], sizeof(int), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(d_invnx, &inv_nx.x, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_invny, &inv_nx.y, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_invnz, &inv_nx.z, sizeof(float), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(d_invhx, &inv_hx.x, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_invhy, &inv_hx.y, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_invhz, &inv_hx.z, sizeof(float), 0, cudaMemcpyHostToDevice);

    float h_c[HALO] = {4.f / 5.f , -1.f / 5.f , 4.f / 105.f, -1.f / 280.f};
    float h_ct[HALO];
    for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx.x;
    cudaMemcpyToSymbol(d_cx, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);
    for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx.y;
    cudaMemcpyToSymbol(d_cy, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);
    for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx.z;
    cudaMemcpyToSymbol(d_cz, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);*/


    // create a common cudaResourceDesc objects
    //struct cudaResourceDesc resDesc;
    //memset(&resDesc, 0, sizeof(resDesc));
   
#if 0
    // Texture Kernel
    // make input image a cudaPitchedPtr for m
    cudaPitchedPtr m_cudaPitchedPtr = make_cudaPitchedPtr((void*)(m), nx[2]*sizeof(ScalarType), nx[2], nx[1]);
    
    // create a cudaExtent for input resolution
    cudaExtent extent = make_cudaExtent(nx[2], nx[1], nx[0]);
    
    // update texture object
    updateTextureFromVolume(m_cudaPitchedPtr, extent, mtex);
   
 
    // Texture gradient
    dim3 threadsPerBlock(1,8,32);
    dim3 numBlocks(nx[0]/threadsPerBlock.x, (nx[1]+7)/threadsPerBlock.y, (nx[2]+31)/threadsPerBlock.z);
    TextureGradientComputeKernel<<<numBlocks, threadsPerBlock>>>(mtex, gx, gy, gz);
    cudaCheckKernelError();
    cudaDeviceSynchronize();
#else
    // Shared Texture Kernel
    // Z-Gradient
    dim3 threadsPerBlock_z(sy, sx, 1);
    dim3 numBlocks_z(nx[2]/sy, nx[1]/sx, nx[0]);
    gradient_z<<<numBlocks_z, threadsPerBlock_z>>>(gz,m);
    cudaCheckKernelError();
    
    // Y-Gradient 
    dim3 threadsPerBlock_y(sxx, syy/perthreadcomp, 1);
    dim3 numBlocks_y(nx[2]/sxx, nx[1]/syy, nx[0]);
    gradient_y<<<numBlocks_y, threadsPerBlock_y>>>(gy, m);
    cudaCheckKernelError();
    
    // X-Gradient
    dim3 threadsPerBlock_x(sxx, syy/perthreadcomp, 1);
    dim3 numBlocks_x(nx[2]/sxx, nx[0]/syy, nx[1]);
    gradient_x<<<numBlocks_x, threadsPerBlock_x>>>(gx, m);
    cudaCheckKernelError();
    cudaDeviceSynchronize();
#endif
    // start recording the interpolation kernel
    /*time = 0; dummy_time = 0; 
    cudaEventRecord(startEvent,0); 
    

    for (int rep=0; rep<repcount; rep++) { 
        TextureGradientComputeKernel<<<numBlocks, threadsPerBlock>>>(mtex, gx, gy, gz);  */
/*        gradient_z<<<numBlocks_z, threadsPerBlock_z>>>(gz, m);
        gradient_y<<<numBlocks_y, threadsPerBlock_y>>>(gy, m);
        gradient_x<<<numBlocks_x, threadsPerBlock_x>>>(gx, m);
*/
/*        if ( cudaSuccess != cudaGetLastError())
                    printf("Error in running gradx kernel\n");
        cudaCheckKernelError();
    }

    cudaEventRecord(stopEvent,0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dummy_time, startEvent, stopEvent);
    time+=dummy_time;
    cudaDeviceSynchronize();
    cudaEventDestroy(startEvent);
    cudaEventDestroy(stopEvent);
    
    // print interpolation time and number of interpolations in Mvoxels/sec
    printf("> gradient avg eval time = %fmsec\n", time/repcount);*/
    //*interp_time += time;
    
    PetscFunctionReturn(ierr);

}



}


#endif
