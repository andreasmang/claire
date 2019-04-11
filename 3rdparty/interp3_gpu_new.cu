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
#include "petsccuda.h"
#include <cuda.h>
#include <cuda_runtime.h>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/extrema.h>
#include <thrust/pair.h>
#include <algorithm>
#include <thrust/device_ptr.h>

#include <memcpy.cu>
#include <cubicPrefilter3D.cu>
#include <bspline_kernel.cu>
#include <lagrange_kernel.cu>
#include "interp3_gpu_new.hpp"

#include "cuda_helper.hpp"
#include "cuda_profiler_api.h"

#include "zeitgeist.hpp"


#define PI ((double)3.14159265358979323846264338327950288419716939937510)
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
__constant__ float d_c[HALO+1];
__constant__ int d_nx, d_ny, d_nz;
//__constant__ float d_invnx, d_invny, d_invnz;
//__constant__ float d_invhx, d_invhy, d_invhz;

template <typename T>
__host__ __device__
inline T rec3_fmaf(T a, T b, T c, T d, T e, T f) {
    return fmaf(a, b, fmaf(c, d, e*f));
    //return a*b + (c*d + e*f);
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


__global__ void printVector(float *m, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<n) printf("m[%d] = %f\n", i , m[i]);
}

/********************************************************************
 * @brief device function for computing the linear index from given 3D indices
 *******************************************************************/
__device__ inline int getLinearIdx(int i, int j, int k) {
    return i*d_ny*d_nz + j*d_nz + k;
}



/********************************************************************
 * @brief function to add 1 to a an array
 *******************************************************************/
__global__ void load_store(float* f1, float* f2) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    f2[i] = f1[i] + 1;
}



/********************************************************************
 * @brief prefilte for z-direction
 *******************************************************************/
__global__ void prefilter_z(float* dfz, float* f) {
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
    
    dfz[globalIdx] = d_c[0]*s_f[sj][sk];
    for(int l=0; l<HALO; l++) {
        dfz[globalIdx] += d_c[l+1] * (s_f[sj][sk+1+l] + s_f[sj][sk-1-l]);
    }
}


/********************************************************************
 * @brief prefilter for y-direction
 *******************************************************************/
__global__ void prefilter_y(float* dfy, float* f) {
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
    dfy[globalIdx] = d_c[0]*s_f[sj][sk];
    for( int l=0; l<HALO; l++) {
        dfy[globalIdx] += d_c[l+1] * ( s_f[sj+1+l][sk] + s_f[sj-1-l][sk]);
    }
  }

}


/********************************************************************
 * @brief prefilter for x-direction
 *******************************************************************/
__global__ void prefilter_x(float* dfx, float* f) {
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
    dfx[globalIdx] = d_c[0]*s_f[si][sk];
    for( int l=0; l<HALO; l++) {
        dfx[globalIdx] += d_c[l+1] * ( s_f[si+1+l][sk] + s_f[si-1-l][sk]);
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
__device__ float cubicTex3D_lagrangeFast( cudaTextureObject_t tex, const float3 coord_grid, const float3 inv_ext)
{
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
    return fmaf(g0.z, Z0, Z2);
    
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
                                   const PetscScalar* xq, 
                                   const PetscScalar* yq, 
                                   const PetscScalar* zq,
                                   PetscScalar* fq, 
                                   const float3 inv_ext) {
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const float3 qcoord = make_float3(zq[tid], yq[tid], xq[tid]);
	__shared__ float fs[64];
	const float3 index = floor(qcoord);
	const float3 fraction = qcoord - index;
	float3 w0, w1, w2, w3;
	lagrange_weights(fraction, w0, w1, w2, w3);
    
    float wx[KERNEL_DIM] = {w0.x, w1.x, w2.x, w3.x};
    float wy[KERNEL_DIM] = {w0.y, w1.y, w2.y, w3.y};
    float wz[KERNEL_DIM] = {w0.z, w1.z, w2.z, w3.z};
    
    if (threadIdx.x == 0) {
        PetscScalar *fp = &f[9 + d_nz*9 + d_ny*d_nz*9];
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
__device__ float cubicTex3D_lagrangeSimple(cudaTextureObject_t tex, float3 coord, const float3 inv_ext)
{
	const float3 coord_grid = coord;
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
    
    float yq = 0;
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
        yq = fmaf(wz[k], sk, yq);
    }
    return yq;
}


/********************************************************************
 * @brief device function to do the interpolation of a single point using the Vanilla Spline Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__device__ float cubicTex3D_splineSimple(cudaTextureObject_t tex, float3 coord, const float3 inv_extent)
{
	// transform the coordinate from [0,extent] to [-0.5, extent-0.5]
	const float3 coord_grid = coord;
	float3 index = floor(coord_grid);
	const float3 fraction = coord_grid - index;
	index = index + 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]

	float result = 0.0f;
	for (float z=-1; z < 2.5f; z++)  //range [-1, 2]
	{
		float bsplineZ = bspline(z-fraction.z);
		float w = (index.z + z)*inv_extent.z;
		for (float y=-1; y < 2.5f; y++)
		{
			float bsplineYZ = bspline(y-fraction.y) * bsplineZ;
			float v = (index.y + y)*inv_extent.y;
			for (float x=-1; x < 2.5f; x++)
			{
				float bsplineXYZ = bspline(x-fraction.x) * bsplineYZ;
				float u = (index.x + x)*inv_extent.z;
				result = fmaf(bsplineXYZ , tex3D<float>(tex, u, v, w), result);
			}
		}
	}
	return result;
}



/********************************************************************
 * @brief device function to do the interpolation of a single point using the Fast Spline Method
 * @parm[in] tex input data texture used for interpolation
 * @parm[in] coord_grid query coordinate
 * @parm[in] inv_reg_extent inverse of the dimension of the 3D grid (1/nx, 1/ny, 1/nz)
 * @parm[out] interpolated value
 *******************************************************************/
__device__ float cubicTex3D_splineFast(cudaTextureObject_t tex, const float3 coord_grid, const float3 inv_reg_extent)
{
	// shift the coordinate from [0,extent] to [-0.5, extent-0.5]
	const float3 index = floor(coord_grid);
	const float3 fraction = coord_grid - index;
	float3 w0, w1, w2, w3;
	bspline_weights(fraction, w0, w1, w2, w3);

	const float3 g0 = w0 + w1;
	const float3 g1 = 1.0f - g0;
	const float3 h0 = ((w1 / g0) - 0.5f + index)*inv_reg_extent;
	const float3 h1 = ((w3 / g1) + 1.5f + index)*inv_reg_extent;
    

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
    return lerp( tex001, tex000, g0.z);
}

__device__ float linTex3D(cudaTextureObject_t tex, const float3 coord_grid, const float3 inv_reg_extent)
{
  const float3 coord = (coord_grid+0.5f)*inv_reg_extent;
  return tex3D<float>(tex, coord.x, coord.y, coord.z);
}

// Fast prefilter for B-Splines
void CubicBSplinePrefilter3D_fast(float *m, int* nx, float *mtemp1, float *mtemp2) {
    
    float time=0, dummy_time=0;
    int repcount = 1;
    cudaEvent_t startEvent, stopEvent;
    cudaEventCreate(&startEvent);
    cudaEventCreate(&stopEvent);

    float h_c[HALO+1];
    h_c[0] = sqrt(3);
    float sum = h_c[0];
    for(int l=1; l<HALO+1; l++) {
        h_c[l] = h_c[l-1]*(sqrt(3) - 2);
        sum += h_c[l]*2;
    }
    for(int l=0; l<HALO; l++) h_c[l] /= sum;
    cudaMemcpyToSymbol(d_c, h_c, sizeof(float)*(HALO+1), 0, cudaMemcpyHostToDevice);
    
    // temporary storage for intermediate results
    //float* mtemp;
    //cudaMalloc((void**) &mtemp, sizeof(float)*nx[0]*nx[1]*nx[2]);

    // Z-Prefilter - WARM UP
    dim3 threadsPerBlock_z(sy, sx, 1);
    dim3 numBlocks_z(nx[2]/sy, nx[1]/sx, nx[0]);
    /*prefilter_z<<<numBlocks_z, threadsPerBlock_z>>>(mtemp,m);
    if ( cudaSuccess != cudaGetLastError())
                printf("Error in running warmup gradz kernel\n");
    cudaCheckKernelError();
    
    // check 1-load, 1-add, 1-store time
    int threads = 256;
    int blocks = nx[2]*nx[1]*nx[0]/threads;
    cudaEventRecord(startEvent,0); 
    load_store<<<blocks, threads>>>(m, mtemp);
    cudaEventRecord(stopEvent,0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dummy_time, startEvent, stopEvent);
    time+=dummy_time;
    cudaDeviceSynchronize();
    printf("> laod-store avg time = %fmsec\n", time);
    time = 0;*/


    
    // Y-Gradient
    dim3 threadsPerBlock_y(sxx, syy/perthreadcomp, 1);
    dim3 numBlocks_y(nx[2]/sxx, nx[1]/syy, nx[0]);
    
    // X-Gradient
    dim3 threadsPerBlock_x(sxx, syy/perthreadcomp, 1);
    dim3 numBlocks_x(nx[2]/sxx, nx[0]/syy, nx[1]);

    // start recording the interpolation kernel
    //cudaEventRecord(startEvent,0); 
    

    //for (int rep=0; rep<repcount; rep++) { 
        // X
        prefilter_x<<<numBlocks_x, threadsPerBlock_x>>>(mtemp1, m);
        if ( cudaSuccess != cudaGetLastError())
            printf("Error in running gradx kernel\n");
        cudaCheckKernelError();
        // Y 
        prefilter_y<<<numBlocks_y, threadsPerBlock_y>>>(mtemp2, mtemp1);
        if ( cudaSuccess != cudaGetLastError())
            printf("Error in running gradx kernel\n");
        cudaCheckKernelError();
        // Z
        prefilter_z<<<numBlocks_z, threadsPerBlock_z>>>(mtemp1, mtemp2);
        if ( cudaSuccess != cudaGetLastError())
            printf("Error in running gradx kernel\n");
        cudaCheckKernelError();
    //}

    /*cudaEventRecord(stopEvent,0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dummy_time, startEvent, stopEvent);
    time+=dummy_time;
    cudaDeviceSynchronize();
    cudaEventDestroy(startEvent);
    cudaEventDestroy(stopEvent);*/
    
    //cudaMemcpy((void*)m, (void*)mtemp, sizeof(float)*nx[0]*nx[1]*nx[2], cudaMemcpyDeviceToDevice);
    //if ( cudaSuccess != cudaGetLastError())
    //            printf("Error in copying data\n");
    
    //if ( mtemp != NULL) cudaFree(mtemp);
    
    // print interpolation time and number of interpolations in Mvoxels/sec
    //printf("> prefilter avg eval time = %fmsec\n", time/repcount);
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
extern "C" cudaTextureObject_t gpuInitEmptyTexture(int* nx) {
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


/********************************************************************
 * @brief interpolation kernel for scalar field
 * @parm[in] yi_tex 3D texture used for interpolation
 * @parm[in] xq,yq,zq query coordinates
 * @parm[in] nx array denoting number of query coordinates in each dimension 
 * @parm[out] yo memory for storing interpolated values
 *******************************************************************/
__global__ void interp3D_kernel(
        cudaTextureObject_t  yi_tex,
        const PetscScalar* xq,
        const PetscScalar* yq,
        const PetscScalar* zq, 
        PetscScalar* yo,
        const float3 inv_nx, int nq)
{
    // Get thread index 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    
    if (tid < nq) {
      float3 qcoord = make_float3(zq[tid], yq[tid], xq[tid]);

      yo[tid] = cubicTex3D_splineFast(yi_tex, qcoord, inv_nx);
        //yo[tid] = cubicTex3D_splineSimple(yi_tex, qcoord, inv_nx);
        //yo[tid] = cubicTex3D_lagrangeSimple(yi_tex, qcoord, inv_nx);
      //yo[tid] = cubicTex3D_lagrangeFast(yi_tex, qcoord, inv_nx);

/*    const float h = 2*PI*inv_nx.x;
      const float3 q = qcoord*h;
      float votrue = computeVx(q.z, q.y, q.x);
      if (tid>=60 && tid<70) {
        printf("tidz = %d  x = %f  y = %f  z = %f  vi = %f  vo = %f  votrue  = %f\n",tid, qcoord.x, qcoord.y, qcoord.z, *((float*)(yi.ptr)+tid), yo[tid], votrue);
      }
*/
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
        const PetscScalar* xq,
        const PetscScalar* yq,
        const PetscScalar* zq, 
        PetscScalar* yo,
        const float3 inv_nx, int nq) {
    // Get thread index 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid < nq) {
      float3 qcoord = make_float3(zq[tid], yq[tid], xq[tid]);
      yo[tid] = linTex3D(yi_tex, qcoord, inv_nx);
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
           const PetscScalar* xq1,
           const PetscScalar* xq2,
           const PetscScalar* xq3,
           PetscScalar* yo,
           float *tmp1, float* tmp2,
           int*  nx,
           cudaTextureObject_t yi_tex,
           int iporder,
           cudaExtent yi_extent,
           const float3 inv_nx,
           long int nq)
{
    if (iporder == 3) {
      cudaMemcpyToSymbol(d_nx, &nx[0], sizeof(int), 0, cudaMemcpyHostToDevice);
      cudaMemcpyToSymbol(d_ny, &nx[1], sizeof(int), 0, cudaMemcpyHostToDevice);
      cudaMemcpyToSymbol(d_nz, &nx[2], sizeof(int), 0, cudaMemcpyHostToDevice);
      CubicBSplinePrefilter3D_fast(yi, nx, tmp1, tmp2);
      cudaPitchedPtr yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(tmp1), nx[2]*sizeof(float), nx[2], nx[1]);
      updateTextureFromVolume(yi_cudaPitchedPtr, yi_extent, yi_tex);
    } else {
      // make input image a cudaPitchedPtr for fi
      cudaPitchedPtr yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(yi), nx[2]*sizeof(float), nx[2], nx[1]);
      //CubicBSplinePrefilter3D_Periodic((float*)yi_cudaPitchedPtr.ptr, (uint)yi_cudaPitchedPtr.pitch, nx[2], nx[1], nx[0]);
      // update texture object
      updateTextureFromVolume(yi_cudaPitchedPtr, yi_extent, yi_tex);
    }
  
    int threads = 256;
    int blocks = (nq+255)/threads;
    
    // launch the interpolation kernel
    switch (iporder) {
    case 1:
      interp3D_kernel_linear<<<blocks,threads>>>(yi_tex, xq1, xq2, xq3, yo, inv_nx, nq);
      break;
    case 3:
      interp3D_kernel<<<blocks,threads>>>(yi_tex, xq1, xq2, xq3, yo, inv_nx, nq);
      break;
    //default:
      // Not implemented
    };
    cudaCheckKernelError();
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
           const PetscScalar* xq1,
           const PetscScalar* xq2,
           const PetscScalar* xq3,
           PetscScalar* yo,
           float *tmp1, float* tmp2,
           int*  nx,
           cudaTextureObject_t yi_tex,
           int iporder,
           float* interp_time)
{
   
    // timing variables
    //float time=0, dummy_time=0;
    //cudaEvent_t startEvent, stopEvent;
    //cudaEventCreate(&startEvent);
    //cudaEventCreate(&stopEvent);

    // define inv of nx for normalizing in texture interpolation
    const float3 inv_nx = make_float3(  1.0f/static_cast<float>(nx[2]),
                                        1.0f/static_cast<float>(nx[1]), 
                                        1.0f/static_cast<float>(nx[0]));
    long int nq = nx[0]*nx[1]*nx[2]; 

    /*cudaMemcpyToSymbol(d_invnx, &inv_nx.x, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_invny, &inv_nx.y, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_invnz, &inv_nx.z, sizeof(float), 0, cudaMemcpyHostToDevice);*/

    // define nxq, the dimensions of the grid
    //const float3 nxq = make_float3( nx[0], nx[1], nx[2]);

    // create a common cudaResourceDesc objects
    //struct cudaResourceDesc resDesc;
    //memset(&resDesc, 0, sizeof(resDesc));
   

    // initiate by computing the bspline coefficients for mt (in-place computation, updates mt)
    //if (iporder == 3) {
    //  cudaMemcpyToSymbol(d_nx, &nx[0], sizeof(int), 0, cudaMemcpyHostToDevice);
    //  cudaMemcpyToSymbol(d_ny, &nx[1], sizeof(int), 0, cudaMemcpyHostToDevice);
    //  cudaMemcpyToSymbol(d_nz, &nx[2], sizeof(int), 0, cudaMemcpyHostToDevice);
    //  CubicBSplinePrefilter3D_fast(yi, nx);
    //}

    // make input image a cudaPitchedPtr for fi
//    cudaPitchedPtr yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(yi), nx[2]*sizeof(float), nx[2], nx[1]);
    //CubicBSplinePrefilter3D_Periodic((float*)yi_cudaPitchedPtr.ptr, (uint)yi_cudaPitchedPtr.pitch, nx[2], nx[1], nx[0]);
    
    // create a cudaExtent for input resolution
    cudaExtent yi_extent = make_cudaExtent(nx[2], nx[1], nx[0]);
  
  
    gpuInterp3Dkernel(yi,xq1,xq2,xq3,yo,tmp1,tmp2,nx,yi_tex,iporder,yi_extent,inv_nx,nq);
    cudaDeviceSynchronize();
    // update texture object
/*    updateTextureFromVolume(yi_cudaPitchedPtr, yi_extent, yi_tex);

    int threads = 256;
    int blocks = (nq+255)/threads;
    
    // start recording the interpolation kernel
    
    //time = 0; dummy_time = 0; 
    //cudaEventRecord(startEvent,0); 
    
    // launch the interpolation kernel
    switch (iporder) {
    case 1:
      interp3D_kernel_linear<<<blocks,threads>>>(yi_tex, xq1, xq2, xq3, yo, inv_nx, nq);
      break;
    case 3:
      interp3D_kernel<<<blocks,threads>>>(yi_tex, xq1, xq2, xq3, yo, inv_nx, nq);
      break;
    //default:
      // Not implemented
    };
    cudaCheckKernelError();

    //cudaEventRecord(stopEvent,0);
    //cudaEventSynchronize(stopEvent);
    //cudaEventElapsedTime(&dummy_time, startEvent, stopEvent);
    //time+=dummy_time;
    //cudaDeviceSynchronize();
    
    //cudaEventDestroy(startEvent);
    //cudaEventDestroy(stopEvent);
    
    // print interpolation time and number of interpolations in Mvoxels/sec
    //printf("> interp time = %fmsec ==> %f MVoxels/sec\n", time, (nq/1E6)/(time/1000));
    //*interp_time += time;
*/
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
void gpuInterpVec3D(
           PetscScalar* yi1, PetscScalar* yi2, PetscScalar* yi3,
           const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3,
           PetscScalar* yo1, PetscScalar* yo2, PetscScalar* yo3,
           float *tmp1, float* tmp2,
           int*  nx, cudaTextureObject_t yi_tex, int iporder, float* interp_time)
{
    // define inv of nx for normalizing in texture interpolation
    const float3 inv_nx = make_float3(  1.0f/static_cast<float>(nx[2]),
                                        1.0f/static_cast<float>(nx[1]), 
                                        1.0f/static_cast<float>(nx[0]));
    long int nq = nx[0]*nx[1]*nx[2]; 
    
    // create a cudaExtent for input resolution
    cudaExtent yi_extent = make_cudaExtent(nx[2], nx[1], nx[0]);
  
    gpuInterp3Dkernel(yi1,xq1,xq2,xq3,yo1,tmp1,tmp2,nx,yi_tex,iporder,yi_extent,inv_nx,nq);
    gpuInterp3Dkernel(yi2,xq1,xq2,xq3,yo2,tmp1,tmp2,nx,yi_tex,iporder,yi_extent,inv_nx,nq);
    gpuInterp3Dkernel(yi3,xq1,xq2,xq3,yo3,tmp1,tmp2,nx,yi_tex,iporder,yi_extent,inv_nx,nq);
    
    cudaDeviceSynchronize();
}
