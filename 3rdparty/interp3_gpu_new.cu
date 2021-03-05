/*--------------------------------------------------------------------------*\
Copyright (c) 2008-2010, Danny Ruijters. All rights reserved.
http://www.dannyruijters.nl/cubicinterpolation/
This file is part of CUDA Cubic B-Spline Interpolation (CI).

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
*  Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
*  Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
*  Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are
those of the authors and should not be interpreted as representing official
policies, either expressed or implied.

When using this code in a scientific project, please cite one or all of the
following papers:
*  Daniel Ruijters and Philippe Thévenaz,
   GPU Prefilter for Accurate Cubic B-Spline Interpolation, 
   The Computer Journal, vol. 55, no. 1, pp. 15-20, January 2012.
   http://dannyruijters.nl/docs/cudaPrefilter3.pdf
*  Daniel Ruijters, Bart M. ter Haar Romeny, and Paul Suetens,
   Efficient GPU-Based Texture Interpolation using Uniform B-Splines,
   Journal of Graphics Tools, vol. 13, no. 4, pp. 61-69, 2008.
\*--------------------------------------------------------------------------*/

#include <stdio.h>
#include "petsc.h"
#include "petscconf.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <chrono>
#include <ctime> 
#include <algorithm>
#include <memcpy.cu>
#include <bspline_kernel.cu>
#include <lagrange_kernel.cu>
#include "interp3_gpu_new.hpp"
#include "cuda_helper.hpp"
#include "cuda_profiler_api.h"
#include "zeitgeist.hpp"
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>


#define PI ((double)3.14159265358979323846264338327950288419716939937510)
#define TWO_PI 2.0f*CUDART_PI
#define KERNEL_DIM 4
#define MAX_BLOCKS 1024
#define HALO 7
#define spencil 4
#define lpencil 32
#define sharedrows 64
#define perthreadcomp 8

const int sx = spencil;
const int sy = sharedrows;
const int sxx = lpencil;
const int syy = sharedrows;

// device constants
/*
__constant__ float d_c[HALO+1];
__constant__ int d_nx, d_ny, d_nz;
__constant__ float d_h[3];
__constant__ float d_iX0[3];
__constant__ float d_iX1[3];
__constant__ int d_istart[3];
__constant__ int d_isize[3];
*/

enum CUBIC_INTERP_TYPE {
    FAST_SPLINE = 0,
    FAST_LAGRANGE = 1,
    SLOW_LAGRANGE = 2,
    SLOW_SPLINE = 3
};


template <typename T>
__host__ __device__
inline T rec3_fmaf(T a, T b, T c, T d, T e, T f) {
    return fmaf(a, b, fmaf(c, d, e*f));
}


template <typename T>
__host__ __device__
inline T rec4_fmaf(T a, T b, T c, T d, T e, T f, T g, T h) {
    return fmaf(a, b, fmaf(c, d, fmaf( e, f, g*h)));
}

__global__ void interp0gpu(float* m, float* q1, float* q2, float *q3, float *q, dim3 nx) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int x1 = q1[i];
  int x2 = q2[i];
  int x3 = q3[i];
  q[i] = m[x3 + x2*nx.z + x1*nx.z*nx.y];
}

void interp0(float* m, float* q1, float* q2, float* q3, float* q, int nx[3]) {
  dim3 n(nx[0], nx[1], nx[2]);
  int nl = nx[0]*nx[1]*nx[2];
  interp0gpu<<<nl/256,256>>>(m,q1,q2,q3,q,n);
}

__global__ void printVectorKernel(ScalarType* m, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<n) printf("m[%d] = %f\n", i , m[i]);
}

__global__ void print3DVectorKernel(ScalarType* arr1, ScalarType* arr2, ScalarType* arr3, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<n) printf("m[%d] = %f\nm[%d] = %f\nm[%d] = %f\n", 3*i , arr1[i], 3*i+1, arr2[i], 3*i+2, arr3[i]);
}

/********************************************************************
 * @brief device function for computing the linear index from given 3D indices
 *******************************************************************/
__device__ inline int getLinearIdx(int i, int j, int k, int3 nl) {
    return i*nl.y*nl.z + j*nl.z + k;
}

/********************************************************************
 * @brief function to add 1 to a an array
 *******************************************************************/
__global__ void load_store(float* f1, float* f2) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    f2[i] = f1[i] + 1;
}

/********************************************************************
 * @brief get coordinates from coordinate array (strided or block)
 *******************************************************************/
__device__ void getCoordinatesBlock(ScalarType* x, ScalarType* y, ScalarType* z, const ScalarType** xq, const int tid) {
  // block access - more memory efficient
  *x = xq[0][tid];
  *y = xq[1][tid];
  *z = xq[2][tid];
}

/********************************************************************
 * @brief get coordinates from coordinate array (strided or block)
 *******************************************************************/
__device__ void getCoordinates(ScalarType* x, ScalarType* y, ScalarType* z, const ScalarType** xq, const int tid) {
#if defined(BLOCK_COORDINATES)
  getCoordinatesBlock(x, y, z, xq, tid);
#else
  // strided
  *x = xq[0][3*tid+0];
  *y = xq[0][3*tid+1];
  *z = xq[0][3*tid+2];
#endif
}

template<typename T, int I> struct __prefilter_value {
  static const T value = __prefilter_value<T, I-1>::value*static_cast<T>(1.732050807568877293527446341505872366942805253810380628055806 - 2.);//(sqrt(3.) - 2.);
  static const T sum = __prefilter_value<T, I-1>::sum + 2*value;
};
template<typename T> struct __prefilter_value<T, 0> {
  static const T value = static_cast<T>(1.732050807568877293527446341505872366942805253810380628055806);//sqrt(3.);;
  static const T sum = value;
};
template<typename T, int I, int L> struct prefilter_value {
  static const T value = __prefilter_value<T,I>::value / __prefilter_value<T,L>::sum;
};

/*template<typename T, int I, int L> inline __device__ T prefilter_value() {
  T value = static_cast<T>(1.732050807568877293527446341505872366942805253810380628055806);//sqrt(3.);
  T sum = value;
  T rval = value;
  for (int i=1; i<=L; ++i) {
    value *= static_cast<T>(1.732050807568877293527446341505872366942805253810380628055806 - 2.);//(sqrt(3.) - 2.);
    sum += value*2;
    if (i == I) rval = value;
  }
  return rval/sum;
};*/

/********************************************************************
 * @brief prefilte for z-direction
 *******************************************************************/
__global__ void prefilter_z(float* dfz, float* f, int3 nl) {
  __shared__ float s_f[sx][sy+2*HALO]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i = blockIdx.z;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x + HALO;       // local k for shared memory ac3ess + halo offset
  int sj = threadIdx.y; // local j for shared memory ac3ess
  int zblock_width, id;
  
  const float d_c[HALO + 1] = { 
    prefilter_value<float, 0, HALO>::value,
    prefilter_value<float, 1, HALO>::value,
    prefilter_value<float, 2, HALO>::value,
    prefilter_value<float, 3, HALO>::value,
    prefilter_value<float, 4, HALO>::value,
    prefilter_value<float, 5, HALO>::value,
    prefilter_value<float, 6, HALO>::value,
    prefilter_value<float, 7, HALO>::value
  };
      
  
  if (blockIdx.x < gridDim.x - 1) {
    zblock_width = blockDim.x;
  }
  else {
    zblock_width = nl.z - blockIdx.x*blockDim.x;
  }
  
  bool internal = (j < nl.y) && (threadIdx.x < zblock_width);
  
  if (internal) {
    id = getLinearIdx(i,j,k, nl);
    s_f[sj][sk] = f[id];
  }

  __syncthreads();
    
  int lid,rid;
  // fill in periodic images in shared memory array 
  if (threadIdx.x < HALO) {
    lid = k%nl.z-HALO;
    if (lid<0) lid+=nl.z;
    s_f[sj][sk-HALO] = f[i*nl.y*nl.z + j*nl.z + lid];
    rid = (k+zblock_width)%nl.z;
    s_f[sj][zblock_width+sk] = f[i*nl.y*nl.z + j*nl.z + rid];
  }

  __syncthreads();
  
  ScalarType result = 0;
  if (internal) {
    result = d_c[0]*s_f[sj][sk];
    for(int l=0; l<HALO; l++) {
        result += d_c[l+1] * (s_f[sj][sk+1+l] + s_f[sj][sk-1-l]);
    }
    dfz[id] = result;
  }

}


/********************************************************************
 * @brief prefilter for y-direction
 *******************************************************************/
__global__ void prefilter_y(float* dfy, float* f, int3 nl) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
    
  int yblock_width, globalIdx, sj;
  bool internal;
  
  if ( blockIdx.y < gridDim.y - 1) {
    yblock_width = syy;
  }
  else {
    yblock_width = nl.y - syy*blockIdx.y;
  }
  
  const float d_c[HALO + 1] = { 
    prefilter_value<float, 0, HALO>::value,
    prefilter_value<float, 1, HALO>::value,
    prefilter_value<float, 2, HALO>::value,
    prefilter_value<float, 3, HALO>::value,
    prefilter_value<float, 4, HALO>::value,
    prefilter_value<float, 5, HALO>::value,
    prefilter_value<float, 6, HALO>::value,
    prefilter_value<float, 7, HALO>::value
  };
  
    
  for(int j = threadIdx.y; j < yblock_width; j += blockDim.y) {
    internal = ((blockIdx.y*syy+j) < nl.y) && (k < nl.z);
    if (internal) {
        globalIdx = getLinearIdx(i, blockIdx.y*syy + j ,k, nl);
        sj = j + HALO;
        s_f[sj][sk] = f[globalIdx];
    }
  }

  __syncthreads();

  
  int lid,rid;
  sj = threadIdx.y + HALO;
  int y = syy*blockIdx.y + threadIdx.y;
  // fill in periodic images in shared memory array 
  if (threadIdx.y < HALO) {
    lid = y%nl.y-HALO;
    if (lid<0) lid+=nl.y;
    s_f[sj-HALO][sk] = f[i*nl.y*nl.z + lid*nl.z + k];
    rid = (y+yblock_width)%nl.y;
    s_f[sj+yblock_width][sk] = f[i*nl.y*nl.z + rid*nl.z + k];
  }

  __syncthreads();
    
  
  ScalarType result;
  for(int j = threadIdx.y; j < yblock_width; j += blockDim.y) {
    result = 0;
    internal = ((blockIdx.y*syy+j) < nl.y) && (k < nl.z);
    if (internal) {
      globalIdx = getLinearIdx(i, blockIdx.y*syy + j ,k, nl);
      sj = j + HALO;
      result = d_c[0]*s_f[sj][sk];
      for( int l=0; l<HALO; l++) {
          result += d_c[l+1] * ( s_f[sj+1+l][sk] + s_f[sj-1-l][sk]);
      }
      dfy[globalIdx] =  result;
    }
  }

}


/********************************************************************
 * @brief prefilter for x-direction
 *******************************************************************/
__global__ void prefilter_x(float* dfx, float* f, int3 nl) {
  
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int j  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
    
  int xblock_width, globalIdx, si;
  bool internal;
  
  if ( blockIdx.y < gridDim.y - 1) {
    xblock_width = syy;
  }
  else {
    xblock_width = nl.x - syy*blockIdx.y;
  }
  
  const float d_c[HALO + 1] = { 
    prefilter_value<float, 0, HALO>::value,
    prefilter_value<float, 1, HALO>::value,
    prefilter_value<float, 2, HALO>::value,
    prefilter_value<float, 3, HALO>::value,
    prefilter_value<float, 4, HALO>::value,
    prefilter_value<float, 5, HALO>::value,
    prefilter_value<float, 6, HALO>::value,
    prefilter_value<float, 7, HALO>::value
  };
    
  for(int i = threadIdx.y; i < xblock_width; i += blockDim.y) {
    internal = ((blockIdx.y*syy + i) < nl.x) && (k < nl.z);
    if (internal) {
        globalIdx = getLinearIdx(blockIdx.y*syy + i, j ,k, nl);
        si = i + HALO;
        s_f[si][sk] = f[globalIdx];
    }
  }

  __syncthreads();

  
  int lid,rid;
  si = threadIdx.y + HALO;
  int x = syy*blockIdx.y + threadIdx.y;
  // fill in periodic images in shared memory array 
  if (threadIdx.y < HALO) {
    lid = x%nl.x-HALO;
    if (lid<0) lid+=nl.x;
    s_f[si-HALO][sk] = f[lid*nl.y*nl.z + j*nl.z + k];
    rid = (x+xblock_width)%nl.x;
    s_f[si+xblock_width][sk] = f[rid*nl.y*nl.z + j*nl.z + k];
  }

  __syncthreads();
    
  
  for(int i = threadIdx.y; i < syy; i += blockDim.y) {
    ScalarType result = 0;
    internal = ((blockIdx.y*syy + i) < nl.x) && (k < nl.z);
    if (internal) {
      int globalIdx = getLinearIdx(blockIdx.y*syy + i , j, k, nl);
      int si = i + HALO;
      result = d_c[0]*s_f[si][sk];
      for( int l=0; l<HALO; l++) {
        result += d_c[l+1] * ( s_f[si+1+l][sk] + s_f[si-1-l][sk]);
      }
      dfx[globalIdx] = result;
    } 
  }

}

/********************************************************************
 * @brief device function to do the interpolation of a single point using the Fast Lagrange Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__global__ void mgpu_cubicTex3DFastLagrange(cudaTextureObject_t tex,
                                        const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
                                        PetscScalar* yo,
                                        const float3 inv_ext, int nq)
{
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const PetscScalar *xq[3] = {xq1, xq2, xq3};
    
    if (tid < nq) {
        //float3 coord_grid = make_float3(zq[tid], yq[tid], xq[tid]);
        float3 coord_grid;
        getCoordinates(&coord_grid.z, &coord_grid.y, &coord_grid.x, xq, tid);
        const float3 index = floor(coord_grid);
        float3 w0, w1, w2, w3;
        lagrange_weights(coord_grid - index, w0, w1, w2, w3);
        
        // compute the locations for the trilinear, bilinear and linear interps
        const float3 g0 = w1 + w2;
        const float3 h0 = (w2/g0 + index + 0.5f)*inv_ext;
        const float idx[2] = { (index.x-0.5f)*inv_ext.x, (index.x+2.5f)*inv_ext.x};
        const float idy[2] = { (index.y-0.5f)*inv_ext.y, (index.y+2.5f)*inv_ext.y};
        const float idz[2] = { (index.z-0.5f)*inv_ext.z, (index.z+2.5f)*inv_ext.z};

        
        /////////////////////// slice 1 ////////////////////////////////////////
        // x weighting
        float Z0,Z2;
        //float row0,row1,row2;
        float point0 = tex3D<float>( tex, idx[0], idy[0], idz[0]); // tex000
        float point1 = tex3D<float>( tex, idx[1], idy[0], idz[0]); // tex100
        float lin0 = tex3D<float>( tex, h0.x, idy[0], idz[0]); // y0z0
        // (w0.x * tex000) + (g0.x * y0z0) + (w3.x * tex100)
        Z0 = w0.y*rec3_fmaf( w0.x,  point0,  g0.x,  lin0,  w3.x,  point1);
        
        point0 = tex3D<float>( tex, idx[0], idy[1], idz[0]); // tex010
        point1 = tex3D<float>( tex, idx[1], idy[1], idz[0]); // tex110
        lin0 = tex3D<float>( tex, h0.x, idy[1], idz[0]); // y1z0
        // (w0.x * tex010) + (g0.x * y1z0) + (w3.x * tex110)
        Z0 = fmaf(w3.y, rec3_fmaf( w0.x,  point0,  g0.x,  lin0,  w3.x,  point1), Z0);

        lin0 = tex3D<float>( tex, idx[0], h0.y, idz[0]); // x0z0
        float lin1 = tex3D<float>( tex, idx[1], h0.y, idz[0]); // x1z0
        float bi0 = tex3D<float>( tex, h0.x, h0.y, idz[0]); // z0
        // (w0.x * x0z0) + (g0.x * z0) + (w3.x * x1z0)
        Z0 = fmaf(g0.y, rec3_fmaf( w0.x,  lin0,    g0.x,  bi0,    w3.x,  lin1), Z0); 
        
        ////////////////////////////////////////// slice 3 //////////////////////////////////////////////////////////
        // x weighting
        point0 = tex3D<float>( tex, idx[0], idy[0], idz[1]); // tex001
        point1 = tex3D<float>( tex, idx[1], idy[0], idz[1]); // tex101
        lin0 = tex3D<float>( tex, h0.x, idy[0], idz[1]); // y0z1
        // (w0.x * tex001) + (g0.x * y0z1) + (w3.x * tex101)
        Z2 = w0.y * rec3_fmaf( w0.x, point0, g0.x, lin0, w3.x, point1);

        point0 = tex3D<float>( tex, idx[0], idy[1], idz[1]); // tex011
        point1 = tex3D<float>( tex, idx[1], idy[1], idz[1]); // tex111
        lin0 = tex3D<float>( tex, h0.x, idy[1], idz[1]); // y1z1
        // (w0.x * tex011) + (g0.x * y1z1) + (w3.x * tex111)
        Z2 = fmaf(w3.y, rec3_fmaf( w0.x, point0, g0.x, lin0,  w3.x, point1), Z2);

        lin0 = tex3D<float>( tex, idx[0], h0.y, idz[1]); // x0z1
        lin1 = tex3D<float>( tex, idx[1], h0.y, idz[1]); // x1z1
        bi0 = tex3D<float>( tex, h0.x, h0.y, idz[1]); // z1
        // (w0.x * x0z1) + (g0.x * z1) + (w3.x * x1z1);
        Z2 = fmaf(g0.y, rec3_fmaf( w0.x, lin0,   g0.x, bi0,   w3.x, lin1), Z2);
        Z2 = fmaf(w0.z, Z0, w3.z*Z2);

        //////////////////////////////////////////////////// slice 2 ////////////////////////////////////////////////
        // single trilinear lookup
        lin0 = tex3D<float>( tex, idx[0], idy[0], h0.z); // x0y0
        lin1 = tex3D<float>( tex, idx[1], idy[0], h0.z); // x1y0
        bi0 = tex3D<float>( tex, h0.x, idy[0], h0.z); // y0
        // (w0.x * x0y0) + (g0.x * y0) + (w3.x * x1y0)
        Z0 = w0.y* rec3_fmaf( w0.x, lin0, g0.x, bi0, w3.x, lin1);

        bi0 = tex3D<float>( tex, idx[0], h0.y, h0.z); // x0
        float bi1 = tex3D<float>( tex, idx[1], h0.y, h0.z); // x1
        float core = tex3D<float>( tex, h0.x, h0.y, h0.z); // core
        // (w0.x * x0) + (g0.x * core) + (w3.x * x1)
        Z0 = fmaf(g0.y, rec3_fmaf( w0.x, bi0, g0.x, core, w3.x, bi1), Z0);


        lin0 = tex3D<float>( tex, idx[0], idy[1], h0.z); // x0y1
        lin1 = tex3D<float>( tex, idx[1], idy[1], h0.z); // x1y1
        bi0 = tex3D<float>( tex, h0.x, idy[1], h0.z); // y1
        // (w0.x * x0y1) + (g0.x * y1) + (w3.x * x1y1)
        Z0 = fmaf(w3.y, rec3_fmaf( w0.x, lin0, g0.x, bi0, w3.x, lin1), Z0);
        yo[tid] = fmaf(g0.z, Z0, Z2);
    }
}
/********************************************************************
 * @brief device function to do the interpolation of a single point using the Fast Lagrange Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__global__ void cubicTex3DFastLagrange(cudaTextureObject_t tex,
                                        const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
                                        PetscScalar* yo,
                                        const float3 inv_ext, int nq)
{
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const PetscScalar *xq[3] = {xq1, xq2, xq3};
    
    if (tid < nq) {
        //float3 coord_grid = make_float3(zq[tid], yq[tid], xq[tid]);
        float3 coord_grid;
        getCoordinatesBlock(&coord_grid.z, &coord_grid.y, &coord_grid.x, xq, tid);
        const float3 index = floor(coord_grid);
        float3 w0, w1, w2, w3;
        lagrange_weights(coord_grid - index, w0, w1, w2, w3);
        
        // compute the locations for the trilinear, bilinear and linear interps
        const float3 g0 = w1 + w2;
        const float3 h0 = (w2/g0 + index + 0.5f)*inv_ext;
        const float idx[2] = { (index.x-0.5f)*inv_ext.x, (index.x+2.5f)*inv_ext.x};
        const float idy[2] = { (index.y-0.5f)*inv_ext.y, (index.y+2.5f)*inv_ext.y};
        const float idz[2] = { (index.z-0.5f)*inv_ext.z, (index.z+2.5f)*inv_ext.z};

        
        /////////////////////// slice 1 ////////////////////////////////////////
        // x weighting
        float Z0,Z2;
        //float row0,row1,row2;
        float point0 = tex3D<float>( tex, idx[0], idy[0], idz[0]); // tex000
        float point1 = tex3D<float>( tex, idx[1], idy[0], idz[0]); // tex100
        float lin0 = tex3D<float>( tex, h0.x, idy[0], idz[0]); // y0z0
        // (w0.x * tex000) + (g0.x * y0z0) + (w3.x * tex100)
        Z0 = w0.y*rec3_fmaf( w0.x,  point0,  g0.x,  lin0,  w3.x,  point1);
        
        point0 = tex3D<float>( tex, idx[0], idy[1], idz[0]); // tex010
        point1 = tex3D<float>( tex, idx[1], idy[1], idz[0]); // tex110
        lin0 = tex3D<float>( tex, h0.x, idy[1], idz[0]); // y1z0
        // (w0.x * tex010) + (g0.x * y1z0) + (w3.x * tex110)
        Z0 = fmaf(w3.y, rec3_fmaf( w0.x,  point0,  g0.x,  lin0,  w3.x,  point1), Z0);

        lin0 = tex3D<float>( tex, idx[0], h0.y, idz[0]); // x0z0
        float lin1 = tex3D<float>( tex, idx[1], h0.y, idz[0]); // x1z0
        float bi0 = tex3D<float>( tex, h0.x, h0.y, idz[0]); // z0
        // (w0.x * x0z0) + (g0.x * z0) + (w3.x * x1z0)
        Z0 = fmaf(g0.y, rec3_fmaf( w0.x,  lin0,    g0.x,  bi0,    w3.x,  lin1), Z0); 
        
        ////////////////////////////////////////// slice 3 //////////////////////////////////////////////////////////
        // x weighting
        point0 = tex3D<float>( tex, idx[0], idy[0], idz[1]); // tex001
        point1 = tex3D<float>( tex, idx[1], idy[0], idz[1]); // tex101
        lin0 = tex3D<float>( tex, h0.x, idy[0], idz[1]); // y0z1
        // (w0.x * tex001) + (g0.x * y0z1) + (w3.x * tex101)
        Z2 = w0.y * rec3_fmaf( w0.x, point0, g0.x, lin0, w3.x, point1);

        point0 = tex3D<float>( tex, idx[0], idy[1], idz[1]); // tex011
        point1 = tex3D<float>( tex, idx[1], idy[1], idz[1]); // tex111
        lin0 = tex3D<float>( tex, h0.x, idy[1], idz[1]); // y1z1
        // (w0.x * tex011) + (g0.x * y1z1) + (w3.x * tex111)
        Z2 = fmaf(w3.y, rec3_fmaf( w0.x, point0, g0.x, lin0,  w3.x, point1), Z2);

        lin0 = tex3D<float>( tex, idx[0], h0.y, idz[1]); // x0z1
        lin1 = tex3D<float>( tex, idx[1], h0.y, idz[1]); // x1z1
        bi0 = tex3D<float>( tex, h0.x, h0.y, idz[1]); // z1
        // (w0.x * x0z1) + (g0.x * z1) + (w3.x * x1z1);
        Z2 = fmaf(g0.y, rec3_fmaf( w0.x, lin0,   g0.x, bi0,   w3.x, lin1), Z2);
        Z2 = fmaf(w0.z, Z0, w3.z*Z2);

        //////////////////////////////////////////////////// slice 2 ////////////////////////////////////////////////
        // single trilinear lookup
        lin0 = tex3D<float>( tex, idx[0], idy[0], h0.z); // x0y0
        lin1 = tex3D<float>( tex, idx[1], idy[0], h0.z); // x1y0
        bi0 = tex3D<float>( tex, h0.x, idy[0], h0.z); // y0
        // (w0.x * x0y0) + (g0.x * y0) + (w3.x * x1y0)
        Z0 = w0.y* rec3_fmaf( w0.x, lin0, g0.x, bi0, w3.x, lin1);

        bi0 = tex3D<float>( tex, idx[0], h0.y, h0.z); // x0
        float bi1 = tex3D<float>( tex, idx[1], h0.y, h0.z); // x1
        float core = tex3D<float>( tex, h0.x, h0.y, h0.z); // core
        // (w0.x * x0) + (g0.x * core) + (w3.x * x1)
        Z0 = fmaf(g0.y, rec3_fmaf( w0.x, bi0, g0.x, core, w3.x, bi1), Z0);


        lin0 = tex3D<float>( tex, idx[0], idy[1], h0.z); // x0y1
        lin1 = tex3D<float>( tex, idx[1], idy[1], h0.z); // x1y1
        bi0 = tex3D<float>( tex, h0.x, idy[1], h0.z); // y1
        // (w0.x * x0y1) + (g0.x * y1) + (w3.x * x1y1)
        Z0 = fmaf(w3.y, rec3_fmaf( w0.x, lin0, g0.x, bi0, w3.x, lin1), Z0);
        yo[tid] = fmaf(g0.z, Z0, Z2);
    }
}

/********************************************************************
 * Fixed departure point
 * @brief device function to do the interpolation of a single point using the Vanilla Lagrange Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__global__ void fixedpointLagrange(PetscScalar* f, 
                                   const PetscScalar** xq, 
                                   PetscScalar* fq, 
                                   const float3 inv_ext, int3 nl) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  float3 qcoord;
  getCoordinates(&qcoord.z, &qcoord.y, &qcoord.x, xq, tid);
	__shared__ float fs[64];
	const float3 index = floor(qcoord);
	const float3 fraction = qcoord - index;
	float3 w0, w1, w2, w3;
	lagrange_weights(fraction, w0, w1, w2, w3);
    
    float wx[KERNEL_DIM] = {w0.x, w1.x, w2.x, w3.x};
    float wy[KERNEL_DIM] = {w0.y, w1.y, w2.y, w3.y};
    float wz[KERNEL_DIM] = {w0.z, w1.z, w2.z, w3.z};
    
    if (threadIdx.x == 0) {
        PetscScalar *fp = &f[9 + nl.z*9 + nl.y*nl.z*9];
    // indices for the source points to be loaded in Shared Memory
        for (int k=0; k<KERNEL_DIM; k++) 
            for (int j=0; j<KERNEL_DIM; j++) 
                for (int i=0; i<KERNEL_DIM; i++) 
                    fs[k + j*4 + i*16] = fp[k + j*4 + i*16];
    }

    __syncthreads();
    
    float sk,sj,temp=0;
    // computation begins
    for (int k=0; k<KERNEL_DIM; k++) {
        sk = 0;
        for (int j=0; j<KERNEL_DIM; j++)  {
           sj = wx[0]*fs[k + j*4 + 0] +
                wx[1]*fs[k + j*4 + 1*16] + 
                wx[2]*fs[k + j*4 + 2*16] + 
                wx[3]*fs[k + j*4 + 3*16];
           sk = fmaf(wy[j], sj, sk);
        }
        temp = fmaf(wz[k], sk, temp);
    }
   fq[tid] = temp; 
}

/********************************************************************
 * @brief device function to do the interpolation of a single point using the Vanilla Lagrange Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__global__ void mgpu_cubicTex3DSlowLagrange(cudaTextureObject_t tex,
                                        const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
                                        PetscScalar* yo,
                                        const float3 inv_ext, int nq)
{
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const PetscScalar *xq[3] = {xq1, xq2, xq3};
    
    if (tid < nq) {
        //float3 coord_grid = make_float3(zq[tid], yq[tid], xq[tid]);
        float3 coord_grid;
        getCoordinates(&coord_grid.z, &coord_grid.y, &coord_grid.x, xq, tid);
        const float3 index = floor(coord_grid);
        const float3 fraction = coord_grid - index;
        float3 w0, w1, w2, w3;
        lagrange_weights(fraction, w0, w1, w2, w3);

        float idx[KERNEL_DIM] = { (index.x-0.5f)*inv_ext.x, 
                        (index.x+0.5f)*inv_ext.x,
                        (index.x+1.5f)*inv_ext.x,
                        (index.x+2.5f)*inv_ext.x};
        
        float wx[KERNEL_DIM] = {w0.x, w1.x, w2.x, w3.x};

        float idy[KERNEL_DIM] = { (index.y-0.5f)*inv_ext.y, 
                        (index.y+0.5f)*inv_ext.y,
                        (index.y+1.5f)*inv_ext.y,
                        (index.y+2.5f)*inv_ext.y};
        float wy[KERNEL_DIM] = {w0.y, w1.y, w2.y, w3.y};

        float idz[KERNEL_DIM] = { (index.z-0.5f)*inv_ext.z, 
                        (index.z+0.5f)*inv_ext.z,
                        (index.z+1.5f)*inv_ext.z,
                        (index.z+2.5f)*inv_ext.z};
        float wz[KERNEL_DIM] = {w0.z, w1.z, w2.z, w3.z};
        
        float res = 0;
        int j,k;
        float sj,sk;
        
        for(k=0; k<KERNEL_DIM; k++){
            sk = 0;
#pragma unroll
            for(j=0; j<KERNEL_DIM; j++){
                sj = rec4_fmaf( wx[0],  tex3D<float>(tex, idx[0], idy[j], idz[k]), 
                                wx[1],  tex3D<float>(tex, idx[1], idy[j], idz[k]),
                                wx[2],  tex3D<float>(tex, idx[2], idy[j], idz[k]),
                                wx[3],  tex3D<float>(tex, idx[3], idy[j], idz[k]));
        
                sk = fmaf(wy[j], sj, sk);
            }
            res = fmaf(wz[k], sk, res);
        }
        yo[tid] = res;
    }
}

/********************************************************************
 * @brief device function to do the interpolation of a single point using the Vanilla Lagrange Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__global__ void cubicTex3DSlowLagrange(cudaTextureObject_t tex,
                                        const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
                                        PetscScalar* yo,
                                        const float3 inv_ext, int nq)
{
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const PetscScalar *xq[3] = {xq1, xq2, xq3};
    
    if (tid < nq) {
        //float3 coord_grid = make_float3(zq[tid], yq[tid], xq[tid]);
        float3 coord_grid;
        getCoordinatesBlock(&coord_grid.z, &coord_grid.y, &coord_grid.x, xq, tid);
        const float3 index = floor(coord_grid);
        const float3 fraction = coord_grid - index;
        float3 w0, w1, w2, w3;
        lagrange_weights(fraction, w0, w1, w2, w3);

        float idx[KERNEL_DIM] = { (index.x-0.5f)*inv_ext.x, 
                        (index.x+0.5f)*inv_ext.x,
                        (index.x+1.5f)*inv_ext.x,
                        (index.x+2.5f)*inv_ext.x};
        
        float wx[KERNEL_DIM] = {w0.x, w1.x, w2.x, w3.x};

        float idy[KERNEL_DIM] = { (index.y-0.5f)*inv_ext.y, 
                        (index.y+0.5f)*inv_ext.y,
                        (index.y+1.5f)*inv_ext.y,
                        (index.y+2.5f)*inv_ext.y};
        float wy[KERNEL_DIM] = {w0.y, w1.y, w2.y, w3.y};

        float idz[KERNEL_DIM] = { (index.z-0.5f)*inv_ext.z, 
                        (index.z+0.5f)*inv_ext.z,
                        (index.z+1.5f)*inv_ext.z,
                        (index.z+2.5f)*inv_ext.z};
        float wz[KERNEL_DIM] = {w0.z, w1.z, w2.z, w3.z};
        
        float res = 0;
        int j,k;
        float sj,sk;
        
        for(k=0; k<KERNEL_DIM; k++){
            sk = 0;
#pragma unroll
            for(j=0; j<KERNEL_DIM; j++){
                sj = rec4_fmaf( wx[0],  tex3D<float>(tex, idx[0], idy[j], idz[k]), 
                                wx[1],  tex3D<float>(tex, idx[1], idy[j], idz[k]),
                                wx[2],  tex3D<float>(tex, idx[2], idy[j], idz[k]),
                                wx[3],  tex3D<float>(tex, idx[3], idy[j], idz[k]));
        
                sk = fmaf(wy[j], sj, sk);
            }
            res = fmaf(wz[k], sk, res);
        }
        yo[tid] = res;
    }
}


/********************************************************************
 * @brief device function to do the interpolation of a single point using the Vanilla Spline Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__global__ void cubicTex3DSlowSpline(cudaTextureObject_t tex,
                                        const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
                                        PetscScalar* yo,
                                        const float3 inv_ext, int nq)
{
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const PetscScalar *xq[3] = {xq1, xq2, xq3};
    
    if (tid < nq) {
        //float3 coord_grid = make_float3(zq[tid], yq[tid], xq[tid]);
        float3 coord_grid;
        getCoordinatesBlock(&coord_grid.z, &coord_grid.y, &coord_grid.x, xq, tid);
        // transform the coordinate from [0,extent] to [-0.5, extent-0.5]
        float3 index = floor(coord_grid);
        const float3 fraction = coord_grid - index;
        index = index + 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]
    
        float res = 0;
        for (float z=-1; z < 2.5f; z++)  //range [-1, 2]
        {
            float bsplineZ = bspline(z-fraction.z);
            float w = (index.z + z)*inv_ext.z;
            for (float y=-1; y < 2.5f; y++)
            {
                float bsplineYZ = bspline(y-fraction.y) * bsplineZ;
                float v = (index.y + y)*inv_ext.y;
                for (float x=-1; x < 2.5f; x++)
                {
                    float bsplineXYZ = bspline(x-fraction.x) * bsplineYZ;
                    float u = (index.x + x)*inv_ext.z;
                    res = fmaf(bsplineXYZ , tex3D<float>(tex, u, v, w), res);
                }
            }
        }
        yo[tid] = res;
    } 
}



/********************************************************************
 * @brief device function to do the interpolation of a single point using the Fast Spline Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__global__ void cubicTex3DFastSpline(cudaTextureObject_t tex,
                                        const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
                                        PetscScalar* yo,
                                        const float3 inv_ext, int nq)
{
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const PetscScalar *xq[3] = {xq1, xq2, xq3};
    
    if (tid < nq) {
        //float3 coord_grid = make_float3(zq[tid], yq[tid], xq[tid]);
        float3 coord_grid;
        getCoordinatesBlock(&coord_grid.z, &coord_grid.y, &coord_grid.x, xq, tid);
        // shift the coordinate from [0,extent] to [-0.5, extent-0.5]
        const float3 index = floor(coord_grid);
        const float3 fraction = coord_grid - index;
        float3 w0, w1, w2, w3;
        bspline_weights(fraction, w0, w1, w2, w3);

        const float3 g0 = w0 + w1;
        const float3 g1 = 1.0f - g0;
        const float3 h0 = ((w1 / g0) - 0.5f + index)*inv_ext;
        const float3 h1 = ((w3 / g1) + 1.5f + index)*inv_ext;
        

        // fetch the eight linear interpolations
        // weighting and fetching is interleaved for performance and stability reasons
        float tex000 = tex3D<float>(tex, h0.x, h0.y, h0.z);
        float tex100 = tex3D<float>(tex, h1.x, h0.y, h0.z);
        //tex000 = g0.x * tex000 + g1.x * tex100;  //weigh along the x-direction
        tex000 = lerp(tex100, tex000, g0.x);
        
        float tex010 = tex3D<float>(tex, h0.x, h1.y, h0.z);
        float tex110 = tex3D<float>(tex, h1.x, h1.y, h0.z);
        //tex010 = g0.x * tex010 + g1.x * tex110;  //weigh along the x-direction
        tex010 = lerp( tex110, tex010, g0.x);
        //tex000 = g0.y * tex000 + g1.y * tex010;  //weigh along the y-direction
        tex000 = lerp( tex010, tex000, g0.y);
        
        float tex001 = tex3D<float>(tex, h0.x, h0.y, h1.z);
        float tex101 = tex3D<float>(tex, h1.x, h0.y, h1.z);
        //tex001 = g0.x * tex001 + g1.x * tex101;  //weigh along the x-direction
        tex001 = lerp( tex101, tex001, g0.x);
        
        float tex011 = tex3D<float>(tex, h0.x, h1.y, h1.z);
        float tex111 = tex3D<float>(tex, h1.x, h1.y, h1.z);
        //tex011 = g0.x * tex011 + g1.x * tex111;  //weigh along the x-direction
        tex011 = lerp( tex111, tex011, g0.x);
        //tex001 = g0.y * tex001 + g1.y * tex011;  //weigh along the y-direction
        tex001 = lerp( tex011, tex001, g0.y);

        //return (g0.z * tex000 + g1.z * tex001);  //weigh along the z-direction
        yo[tid] = lerp( tex001, tex000, g0.z);
    }
}


__device__ float linTex3D(cudaTextureObject_t tex, const float3 coord_grid, const float3 inv_reg_extent)
{
  const float3 coord = (coord_grid+0.5f)*inv_reg_extent;
  return tex3D<float>(tex, coord.x, coord.y, coord.z);
}

void getThreadBlockDimensionsX(dim3& threads, dim3& blocks, int* nx) {
  threads.x = sxx;
  threads.y = syy/perthreadcomp;
  threads.z = 1;
  blocks.x = (nx[2]+sxx-1)/sxx;
  blocks.y = (nx[0]+syy-1)/syy;
  blocks.z = nx[1];
}

void getThreadBlockDimensionsY(dim3& threads, dim3& blocks, int* nx) {
  threads.x = sxx;
  threads.y = syy/perthreadcomp;
  threads.z = 1;
  blocks.x = (nx[2]+sxx-1)/sxx;
  blocks.y = (nx[1]+syy-1)/syy;
  blocks.z = nx[0];
}

void getThreadBlockDimensionsZ(dim3& threads, dim3& blocks, int* nx) {
  threads.x = sy;
  threads.y = sx;
  threads.z = 1;
  blocks.x = (nx[2]+sy-1)/sy;
  blocks.y = (nx[1]+sx-1)/sx;
  blocks.z = nx[0];
}

// Fast prefilter for B-Splines
void CubicBSplinePrefilter3D_fast(float *m, int* nx, float *mtemp1, float *mtemp2) {
    
    cudaEvent_t startEvent, stopEvent;
    cudaEventCreate(&startEvent);
    cudaEventCreate(&stopEvent);

    /*float h_c[HALO+1];
    h_c[0] = sqrt(3);
    float sum = h_c[0];
    for(int l=1; l<HALO+1; l++) {
        h_c[l] = h_c[l-1]*(sqrt(3) - 2);
        sum += h_c[l]*2;
    }
    for(int l=0; l<HALO; l++) h_c[l] /= sum;
    cudaMemcpyToSymbol(d_c, h_c, sizeof(float)*(HALO+1), 0, cudaMemcpyHostToDevice);*/
    
    dim3 threadsPerBlock_x, numBlocks_x;
    dim3 threadsPerBlock_y, numBlocks_y;
    dim3 threadsPerBlock_z, numBlocks_z;
    getThreadBlockDimensionsX(threadsPerBlock_x, numBlocks_x, nx);
    getThreadBlockDimensionsY(threadsPerBlock_y, numBlocks_y, nx);
    getThreadBlockDimensionsZ(threadsPerBlock_z, numBlocks_z, nx);
    
    int3 nl;
    nl.x = nx[0]; nl.y = nx[1]; nl.z = nx[2];

    // X
    prefilter_x<<<numBlocks_x, threadsPerBlock_x>>>(mtemp1, m, nl);
    if ( cudaSuccess != cudaGetLastError())
        printf("Error in running gradx kernel\n");
    cudaCheckKernelError();
    // Y 
    prefilter_y<<<numBlocks_y, threadsPerBlock_y>>>(mtemp2, mtemp1, nl);
    if ( cudaSuccess != cudaGetLastError())
        printf("Error in running grady kernel\n");
    cudaCheckKernelError();
    // Z
    prefilter_z<<<numBlocks_z, threadsPerBlock_z>>>(mtemp1, mtemp2, nl);
    if ( cudaSuccess != cudaGetLastError())
        printf("Error in running gradz kernel\n");
    cudaCheckKernelError();
}




/********************************************************************
 * @brief function to create a 3D texture from the given cuda Pitched Pointer denoting volume (3D) data
 *******************************************************************/
extern "C" cudaTextureObject_t initTextureFromVolume(cudaPitchedPtr volume, cudaExtent extent) {
   cudaError_t err = cudaSuccess;
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
   cudaArray* cuArray;
   err = cudaMalloc3DArray(&cuArray, &channelDesc, extent, 0);
   if (err != cudaSuccess){
        fprintf(stderr, "Failed to allocate 3D cudaArray (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
   }
   
   cudaMemcpy3DParms p = {0};
   p.srcPtr = volume;
   p.dstArray = cuArray;
   p.extent = extent;
   p.kind = cudaMemcpyDeviceToDevice;
   err = cudaMemcpy3D(&p);
   if (err != cudaSuccess){
        fprintf(stderr, "Failed to copy 3D memory to cudaArray (error code %s)!\n", cudaGetErrorString(err));
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
    texDesc.filterMode = cudaFilterModeLinear;
    texDesc.normalizedCoords = 1;

    cudaTextureObject_t texObj = 0;
    err = cudaCreateTextureObject( &texObj, &resDesc, &texDesc, NULL);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to create texture (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    return texObj;
}

/********************************************************************
 * @brief create texture object with empty data (cudaArray)
 *******************************************************************/
extern "C" cudaTextureObject_t gpuInitEmptyTexture(IntType* nx) {

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
    texDesc.filterMode = cudaFilterModeLinear;
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
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy 3D memory to cudaArray (error name %s = %s)!\n", cudaGetErrorName(err), cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

/********************************************************************
 * @brief linear interpolation kernel for scalar field
 * @parm[in] yi_tex 3D texture used for interpolation
 * @parm[in] xq,yq,zq query coordinates
 * @parm[in] nx array denoting number of query coordinates in each dimension 
 * @parm[out] yo memory for storing interpolated values
 *******************************************************************/
__global__ void mgpu_interp3D_kernel_linear(
        cudaTextureObject_t  yi_tex,
        const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
        PetscScalar* yo,
        const float3 inv_nx, int nq) {
    // Get thread index 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const PetscScalar *xq[3] = {xq1, xq2, xq3};
    if (tid < nq) {
      float3 coord_grid;
      getCoordinates(&coord_grid.z, &coord_grid.y, &coord_grid.x, xq, tid);
      yo[tid] = linTex3D(yi_tex, coord_grid, inv_nx);
    }
}

/********************************************************************
 * @brief linear interpolation kernel for scalar field
 * @parm[in] yi_tex 3D texture used for interpolation
 * @parm[in] xq,yq,zq query coordinates
 * @parm[in] nx array denoting number of query coordinates in each dimension 
 * @parm[out] yo memory for storing interpolated values
 *******************************************************************/
__global__ void interp3D_kernel_linear(
        cudaTextureObject_t  yi_tex,
        const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
        PetscScalar* yo,
        const float3 inv_nx, int nq) {
    // Get thread index 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const PetscScalar *xq[3] = {xq1, xq2, xq3};
    if (tid < nq) {
      float3 coord_grid;
      getCoordinatesBlock(&coord_grid.z, &coord_grid.y, &coord_grid.x, xq, tid);
      yo[tid] = linTex3D(yi_tex, coord_grid, inv_nx);
    }
}


/********************************************************************
 * @brief host function to do interpolation of a scalar field
 * @parm[in] yi input data values 
 * @parm[in] xq1,yq1,zq1 query coordinates
 * @parm[in] yi_tex texture object
 * @parm[in] nx array denoting number of query coordinates in each dimension 
 * @parm[out] yo interpolated values
 * @parm[out] interp_time time for computing the interpolation
 *******************************************************************/
void gpuInterp3Dkernel(
           PetscScalar* yi,
           const PetscScalar** xq,
           PetscScalar* yo,
           float *tmp1, float* tmp2,
           int*  nx,
           cudaTextureObject_t yi_tex,
           int iporder,
           cudaExtent yi_extent,
           const float3 inv_nx,
           long int nq,
           cudaStream_t* stream)
{
    
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // SET cubic interpolation type here
    enum CUBIC_INTERP_TYPE interp_type = FAST_LAGRANGE;
    if (nprocs == 1)
      interp_type = FAST_SPLINE;

    cudaPitchedPtr yi_cudaPitchedPtr;
    if (iporder == 3) {
      //cudaMemcpyToSymbol(d_nx, &nx[0], sizeof(int), 0, cudaMemcpyHostToDevice);
      //cudaMemcpyToSymbol(d_ny, &nx[1], sizeof(int), 0, cudaMemcpyHostToDevice);
      //cudaMemcpyToSymbol(d_nz, &nx[2], sizeof(int), 0, cudaMemcpyHostToDevice);
      switch (interp_type) {
        case FAST_SPLINE:
          if (nprocs == 1) {
            CubicBSplinePrefilter3D_fast(yi, nx, tmp1, tmp2);
            yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(tmp1), nx[2]*sizeof(float), nx[2], nx[1]);
          }
            break;
        case SLOW_SPLINE:
          if (nprocs == 1 ) {
            CubicBSplinePrefilter3D_fast(yi, nx, tmp1, tmp2);
            yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(tmp1), nx[2]*sizeof(float), nx[2], nx[1]);
          }
            break;
        case FAST_LAGRANGE:
            yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(yi), nx[2]*sizeof(float), nx[2], nx[1]);
            break;
        case SLOW_LAGRANGE:
            yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(yi), nx[2]*sizeof(float), nx[2], nx[1]);
            break;
      };
    } else {
      // make input image a cudaPitchedPtr for fi
      yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(yi), nx[2]*sizeof(float), nx[2], nx[1]);
    }
    updateTextureFromVolume(yi_cudaPitchedPtr, yi_extent, yi_tex);
  
    int threads = 256;
    int blocks = (nq+255)/threads;
     
    // launch the interpolation kernel
    if (nprocs == 1 ) {
      switch (iporder) {
      case 1:
        interp3D_kernel_linear<<<blocks,threads, 0, *stream>>>(yi_tex, xq[0], xq[1], xq[2], yo, inv_nx, nq);
        break;
      case 3:
        switch (interp_type) {
          case FAST_SPLINE:
              cubicTex3DFastSpline<<<blocks,threads, 0, *stream>>>(yi_tex, xq[0], xq[1], xq[2], yo, inv_nx, nq);
              break;
          case SLOW_SPLINE:
              cubicTex3DSlowSpline<<<blocks,threads, 0, *stream>>>(yi_tex, xq[0], xq[1], xq[2], yo, inv_nx, nq);
              break;
          case FAST_LAGRANGE:
              cubicTex3DFastLagrange<<<blocks,threads, 0, *stream>>>(yi_tex, xq[0], xq[1], xq[2], yo, inv_nx, nq);
              break;
          case SLOW_LAGRANGE:
              cubicTex3DSlowLagrange<<<blocks,threads, 0, *stream>>>(yi_tex, xq[0], xq[1], xq[2], yo, inv_nx, nq);
              break;
        };
        break;
      };
    } else {
      switch (iporder) {
      case 1:
        mgpu_interp3D_kernel_linear<<<blocks,threads, 0, *stream>>>(yi_tex, xq[0], xq[1], xq[2], yo, inv_nx, nq);
        break;
      case 3:
        switch (interp_type) {
          case FAST_SPLINE:
              break;
          case SLOW_SPLINE:
              break;
          case FAST_LAGRANGE:
              mgpu_cubicTex3DFastLagrange<<<blocks,threads, 0, *stream>>>(yi_tex, xq[0], xq[1], xq[2], yo, inv_nx, nq);
              break;
          case SLOW_LAGRANGE:
              mgpu_cubicTex3DSlowLagrange<<<blocks,threads, 0, *stream>>>(yi_tex, xq[0], xq[1], xq[2], yo, inv_nx, nq);
              break;
        };
        break;
      };
    }
}

/********************************************************************
 * @brief host function to do interpolation of a scalar field
 * @parm[in] yi input data values 
 * @parm[in] xq1,yq1,zq1 query coordinates
 * @parm[in] yi_tex texture object
 * @parm[in] nx array denoting number of query coordinates in each dimension 
 * @parm[out] yo interpolated values
 * @parm[out] interp_time time for computing the interpolation
 *******************************************************************/
void gpuInterp3D(
           PetscScalar* yi,
           const PetscScalar** xq,
           PetscScalar* yo,
           float *tmp1, float* tmp2,
           IntType*  inx,
           long int nq,
           cudaTextureObject_t yi_tex,
           int iporder,
           float* interp_time)
{
    int nx[3];
    nx[0] = inx[0]; nx[1] = inx[1]; nx[2] = inx[2];
    // define inv of nx for normalizing in texture interpolation
    const float3 inv_nx = make_float3(  1.0f/static_cast<float>(nx[2]),
                                        1.0f/static_cast<float>(nx[1]), 
                                        1.0f/static_cast<float>(nx[0]));
    
    // create a cudaExtent for input resolution
    cudaExtent yi_extent = make_cudaExtent(nx[2], nx[1], nx[0]);
    
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    gpuInterp3Dkernel(yi,xq,yo,tmp1,tmp2,nx,yi_tex,iporder,yi_extent,inv_nx,nq, &stream);

    cudaStreamSynchronize(stream);
    cudaStreamDestroy(stream);


    cudaDeviceSynchronize();
}

/********************************************************************
 * @brief host function to do interpolation of a scalar field
 * @parm[in] yi input data values 
 * @parm[in] xq1,yq1,zq1 query coordinates
 * @parm[in] yi_tex texture object
 * @parm[in] nx array denoting number of query coordinates in each dimension  (this will be isize when using multi-GPU)
 * @parm[out] yo interpolated values
 * @parm[out] interp_time time for computing the interpolation
 *******************************************************************/
void gpuInterpVec3D(
           PetscScalar* yi1, PetscScalar* yi2, PetscScalar* yi3,
           const PetscScalar** xq,
           PetscScalar* yo1, PetscScalar* yo2, PetscScalar* yo3,
           float *tmp1, float* tmp2,
           IntType*  inx, long int nq, cudaTextureObject_t yi_tex, int iporder, float* interp_time)
{
    int nx[3];
    nx[0] = inx[0]; nx[1] = inx[1]; nx[2] = inx[2];
    // define inv of nx for normalizing in texture interpolation
    const float3 inv_nx = make_float3(  1.0f/static_cast<float>(nx[2]),
                                        1.0f/static_cast<float>(nx[1]), 
                                        1.0f/static_cast<float>(nx[0]));

    // create a cudaExtent for input resolution
    cudaExtent yi_extent = make_cudaExtent(nx[2], nx[1], nx[0]);
    
    // launch the 3 interpolations in 3 different streams to have some speed up
    // in case there is not enough work for the GPU
    cudaStream_t stream[3];
    for (int i=0; i<3; i++) {
      cudaStreamCreate(&stream[i]);
    }

    gpuInterp3Dkernel(yi1,xq,yo1,tmp1,tmp2,nx,yi_tex,iporder,yi_extent,inv_nx,nq, &stream[0]);
    gpuInterp3Dkernel(yi2,xq,yo2,tmp1,tmp2,nx,yi_tex,iporder,yi_extent,inv_nx,nq, &stream[1]);
    gpuInterp3Dkernel(yi3,xq,yo3,tmp1,tmp2,nx,yi_tex,iporder,yi_extent,inv_nx,nq, &stream[2]);

    for (int i=0; i<3; i++) 
      cudaStreamSynchronize(stream[i]);
    
    for (int i=0; i<3; i++)
      cudaStreamDestroy(stream[i]);

    
    cudaDeviceSynchronize();
}

void getMax(ScalarType* x, int nl, ScalarType* max) {
    thrust::device_ptr<ScalarType> x_thrust;
    x_thrust = thrust::device_pointer_cast (x);
    // find the max itr
    thrust::device_vector<ScalarType>::iterator it = thrust::max_element(x_thrust, x_thrust + nl);
    *max = *it;
}
    
__global__ void normalizeQueryPointsKernel(ScalarType* xq1, ScalarType* xq2, ScalarType* xq3, ScalarType* all_query_points, int nq, const float3 ng, const float3 offset) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x; 

    if (tid<nq) {
#ifdef BLOCK_COORDINATES
        xq1[tid] = fmaf(all_query_points[tid*3+0], ng.x, offset.x);
        xq2[tid] = fmaf(all_query_points[tid*3+1], ng.y, offset.y);
        xq3[tid] = fmaf(all_query_points[tid*3+2], ng.z, offset.z);
#else
        all_query_points[tid*3+0] = fmaf(all_query_points[tid*3+0], ng.x, offset.x);
        all_query_points[tid*3+1] = fmaf(all_query_points[tid*3+1], ng.y, offset.y);
        all_query_points[tid*3+2] = fmaf(all_query_points[tid*3+2], ng.z, offset.z);
#endif
    }
}


void normalizeQueryPoints(ScalarType* xq1, ScalarType* xq2, ScalarType* xq3, ScalarType* all_query_points, int nq, IntType* isize, IntType* nx, int* procid, int nghost) {
    
    const float3 offset = make_float3( static_cast<float>(nghost-procid[0]*isize[0]),
                                       static_cast<float>(0*nghost-0*procid[1]*isize[1]),
                                       static_cast<float>(0*nghost));
    const float3 ng = make_float3( nx[0], nx[1], nx[2] );

    int threads = 256;
    int blocks = (nq+threads-1)/threads;
    normalizeQueryPointsKernel<<<blocks,threads>>>(xq1, xq2, xq3, all_query_points, nq, ng, offset);
    cudaDeviceSynchronize();
}

void printGPUVector(ScalarType* arr, int nq) {
    int threads = 256;
    int blocks = (nq+threads-1)/threads;
    printVectorKernel<<<blocks, threads>>>(arr, nq);
}


void printGPU3DVector(ScalarType* arr1, ScalarType* arr2, ScalarType* arr3, int nq) {
    int threads = 256;
    int blocks = (nq+255)/threads;
    print3DVectorKernel<<<blocks, threads>>>(arr1, arr2, arr3, nq);
}


__global__ void copyQueryValuesKernel(ScalarType* dst, ScalarType* src, int* index, int len) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    int ind;
    if (tid < len) {
        ind = index[tid];
        dst[ind] = src[tid];
    }
}
    

void copyQueryValues(ScalarType* dst, ScalarType* src, int* index, int len) {
    int threads = 256;
    int blocks = (len+threads-1)/threads;

    copyQueryValuesKernel<<<blocks, threads>>>(dst, src, index, len);
    cudaDeviceSynchronize();
}

__global__ void enforcePeriodicityKernel(ScalarType* xq, ScalarType* yq, ScalarType* zq, int len, float3 dh) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    ScalarType x, y, z;
    if (tid < len) {
        x = xq[tid];
        y = yq[tid];
        z = zq[tid];
        while (x <= -dh.x) { x += 1; }
        while (y <= -dh.y) { y += 1; }
        while (z <= -dh.z) { z += 1; }
        
        while (x >= 1) { x -= 1; }
        while (y >= 1) { y -= 1; }
        while (z >= 1) { z -= 1; }

        xq[tid] = x;
        yq[tid] = y;
        zq[tid] = z;
    }
}

void enforcePeriodicity(ScalarType* xq, ScalarType* yq, ScalarType* zq, ScalarType* h, int len) {
    int threads = 256;
    int blocks = (len+threads-1)/threads;
    
    // copy constant h to device
    //cudaMemcpyToSymbol(d_h, h, sizeof(float)*(3), 0, cudaMemcpyHostToDevice);
    
    float3 dh;
    dh.x = h[0]; dh.y = h[1]; dh.z = h[2];

    enforcePeriodicityKernel<<<blocks, threads>>>(xq, yq, zq, len, dh);
    cudaDeviceSynchronize();
}

__global__ void checkDomainKernel(short* which_proc, ScalarType* xq, ScalarType* yq, ScalarType* zq, int len, int procid, const int2 isize, int c_dim1, float3 iX0, float3 iX1, float3 dh) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    ScalarType x,y,z;
    int dproc0, dproc1;
    int proc;
    if (tid < len) {
        x = xq[tid];
        y = yq[tid];
        z = zq[tid];
        if ( iX0.x-dh.x <= x && x <= iX1.x+dh.x &&
             iX0.y-dh.y <= y && y <= iX1.y+dh.y &&
             iX0.z-dh.z <= z && z <= iX1.z+dh.z ) {
            which_proc[tid] = static_cast<short>(procid);
        } else {
            dproc0=(int)(x/dh.x)/isize.x;
            dproc1=(int)(y/dh.y)/isize.y;
            proc=dproc0*c_dim1+dproc1; 
            which_proc[tid] = static_cast<short>(proc);
        }
    }
}

void checkDomain(short* which_proc, ScalarType* xq, ScalarType* yq, ScalarType* zq, ScalarType* iX0, ScalarType* iX1, ScalarType* h, int len, int procid, int isize0, int isize1, int c_dim1) {
    
    int threads = 256;
    int blocks = (len+threads-1)/threads;
    
    // copy constant h to device
    //cudaMemcpyToSymbol(d_h, h, 3*sizeof(float), 0, cudaMemcpyHostToDevice);
    //cudaMemcpyToSymbol(d_iX0, iX0, 3*sizeof(float), 0, cudaMemcpyHostToDevice);
    //cudaMemcpyToSymbol(d_iX1, iX1, 3*sizeof(float), 0, cudaMemcpyHostToDevice);
    const int2 isize = make_int2(isize0, isize1);
    
    float3 diX0, diX1, dh;
    diX0.x = iX0[0]; diX0.y = iX0[1]; diX0.z = iX0[2];
    diX1.x = iX1[0]; diX1.y = iX1[1]; diX1.z = iX1[2];
    dh.x = h[0]; dh.y = h[1]; dh.z = h[2];


    checkDomainKernel<<<blocks, threads>>>(which_proc, xq, yq, zq, len, procid, isize, c_dim1, diX0, diX1, dh);
    cudaDeviceSynchronize();
}


__host__ __device__ 
void TestFunction(ScalarType *val, const ScalarType x, const ScalarType y, const ScalarType z, int caseid) {
      *val = (caseid+1)*( sinf(8*x)*sinf(8*x) + sinf(2*y)*sinf(2*y) + sinf(4*z)*sinf(4*z) )/3.0;
      //*val = x;
}

__global__ void setup_kernel(curandState *state, const int3 size, const int3 start, const int3 n) {

    int ix = threadIdx.x + blockDim.x * blockIdx.x;
    int iy = threadIdx.y + blockDim.y * blockIdx.y;
    int iz = threadIdx.z + blockDim.z * blockIdx.z;
    

    int idx = ix*size.y*size.z + iy*size.z + iz; // global index local to GPU
    if (idx < size.x*size.y*size.z) {
      int i = ( (ix+start.x) * n.y * n.z ) + ( (iy+start.y) * n.z ) + iz+start.z; // global index to the grid
      curand_init(1234, i*3, 0, &state[idx]);
    }
}

__global__ void initializeGridKernel(ScalarType* xq, ScalarType* yq, ScalarType* zq, ScalarType* f, ScalarType* ref, const float3 h, const int3 size, const int3 start, const int3 n, int caseid) {
    int ix = threadIdx.x + blockDim.x * blockIdx.x;
    int iy = threadIdx.y + blockDim.y * blockIdx.y;
    int iz = threadIdx.z + blockDim.z * blockIdx.z;

    int i = (ix*size.y*size.z) + (iy*size.z) + iz;
    
    ScalarType x,y,z;

    if (i < size.x*size.y*size.z) {
        x = h.x*(float)(ix+start.x);
        y = h.y*(float)(iy+start.y);
        z = h.z*(float)(iz+start.z);
        
        TestFunction(&f[i], x, y, z, caseid);

        //x = 0;
        //y = 0;
        //z = 0;
        
        //float perturb=sinf(x)*sinf(y)*sinf(z); 
        //x += h.x*perturb;
        //y += h.y*perturb;
        //z += h.z*perturb;
        
        x -= 0.25*0.5*sinf(z)*cosf(y)*sinf(y);
        y -= 0.25*0.5*sinf(x)*cosf(z)*sinf(z);
        z -= 0.25*0.5*sinf(y)*cosf(x)*sinf(x);
        
        //x -= 0.25*0.5*sinf(z)*cosf(y)*sinf(y);
        //y -= 0.25*0.5*sinf(x)*cosf(z)*sinf(z);
        //z -= 0.25*0.5*sinf(y)*cosf(x)*sinf(x);
        
        //x += h.x*(dhx[i*3]*2-1);
        //y += h.y*(dhx[i*3+1]*2-1);
        //z += h.z*(dhx[i*3+2]*2-1);

        //x += h.x*(curand_uniform(&state[i])*2.0-1.0);
        //y += h.y*(curand_uniform(&state[i])*2.0-1.0);
        //z += h.z*(curand_uniform(&state[i])*2.0-1.0);
        
        //x = curand_uniform(&state)*2*PI;
        //y = curand_uniform(&state)*2*PI;
        //z = curand_uniform(&state)*2*PI;


        if (x > 2.*PI) x -= 2.*PI;
        if (y > 2.*PI) y -= 2.*PI;
        if (z > 2.*PI) z -= 2.*PI;
        
        if (x < 0) x += 2.*PI;
        if (y < 0) y += 2.*PI;
        if (z < 0) z += 2.*PI;
        
        TestFunction(&ref[i], x, y, z, caseid);
        
        xq[i] = x/(2.*PI);
        yq[i] = y/(2.*PI);
        zq[i] = z/(2.*PI);
    }
}

void initializeGrid(ScalarType* xq, ScalarType* yq, ScalarType* zq, ScalarType* f, ScalarType* ref, ScalarType* h, IntType* isize, IntType* istart, IntType* nx, int caseid) {
    dim3 threads(1,32,32);
    dim3 blocks((isize[0]+threads.x-1)/threads.x, (isize[1]+threads.y-1)/threads.y, (isize[2]+threads.z-1)/threads.z);


    const float3 hx = make_float3(h[0], h[1], h[2]);
    const int3 size = make_int3(isize[0], isize[1], isize[2]);
    const int3 start = make_int3(istart[0], istart[1], istart[2]);
    const int3 n = make_int3(nx[0], nx[1], nx[2]);

    initializeGridKernel<<<blocks, threads>>>(xq, yq, zq, f, ref, hx, size, start, n, caseid);
    cudaDeviceSynchronize();
}

__global__ void testKernel(ScalarType* f) {
  int tid = threadIdx.x + blockIdx.x*blockDim.x;
  
  ScalarType x = 0;
  for (int i=0; i<500; i++) {
    x += sinf((float)i);
  }
  f[tid] = x;
}

void test(ScalarType* f, int nq) {

  int threads = 256;
  int blocks = (nq+threads-1)/threads;

  testKernel<<<blocks, threads>>>(f);
  cudaDeviceSynchronize();
}


