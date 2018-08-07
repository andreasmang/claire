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


#define PI ((double)3.14159265358979323846264338327950288419716939937510)
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
            _a > _b ? _a : _b; })

#define KERNEL_DIM 4
#define MAX_BLOCKS 1024

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
    float idx[2] = { (index.x-0.5f)*inv_ext.x, (index.x+2.5f)*inv_ext.x};
    float idy[2] = { (index.y-0.5f)*inv_ext.y, (index.y+2.5f)*inv_ext.y};
    float idz[2] = { (index.z-0.5f)*inv_ext.z, (index.z+2.5f)*inv_ext.z};

    // single trilinear lookup
    float core = tex3D<float>( tex, h0.x, h0.y, h0.z);

    // 6 bilinear lookups
    float z0 = tex3D<float>( tex, h0.x, h0.y, idz[0]);
    float z1 = tex3D<float>( tex, h0.x, h0.y, idz[1]);
    float y0 = tex3D<float>( tex, h0.x, idy[0], h0.z);
    float y1 = tex3D<float>( tex, h0.x, idy[1], h0.z);
    float x0 = tex3D<float>( tex, idx[0], h0.y, h0.z);
    float x1 = tex3D<float>( tex, idx[1], h0.y, h0.z);

    // 12 linear lookups
    // along z-axis
    float x0y0 = tex3D<float>( tex, idx[0], idy[0], h0.z);
    float x1y1 = tex3D<float>( tex, idx[1], idy[1], h0.z);
    float x0y1 = tex3D<float>( tex, idx[0], idy[1], h0.z);
    float x1y0 = tex3D<float>( tex, idx[1], idy[0], h0.z);
    // aling y-axis
    float x0z0 = tex3D<float>( tex, idx[0], h0.y, idz[0]);
    float x1z1 = tex3D<float>( tex, idx[1], h0.y, idz[1]);
    float x0z1 = tex3D<float>( tex, idx[0], h0.y, idz[1]);
    float x1z0 = tex3D<float>( tex, idx[1], h0.y, idz[0]);
    // along x-axis
    float y0z0 = tex3D<float>( tex, h0.x, idy[0], idz[0]);
    float y1z1 = tex3D<float>( tex, h0.x, idy[1], idz[1]);
    float y0z1 = tex3D<float>( tex, h0.x, idy[0], idz[1]);
    float y1z0 = tex3D<float>( tex, h0.x, idy[1], idz[0]);

    // 8 single point look ups
    float tex000 = tex3D<float>( tex, idx[0], idy[0], idz[0]);
    float tex100 = tex3D<float>( tex, idx[1], idy[0], idz[0]);
    float tex010 = tex3D<float>( tex, idx[0], idy[1], idz[0]);
    float tex110 = tex3D<float>( tex, idx[1], idy[1], idz[0]);
    float tex001 = tex3D<float>( tex, idx[0], idy[0], idz[1]);
    float tex101 = tex3D<float>( tex, idx[1], idy[0], idz[1]);
    float tex011 = tex3D<float>( tex, idx[0], idy[1], idz[1]);
    float tex111 = tex3D<float>( tex, idx[1], idy[1], idz[1]);

    // weighting in x direction
    // slice 1 (z=0)
    float row0 = rec3_fmaf( w0.x,  tex000,  g0.x,  y0z0,  w3.x,  tex100);
    float row1 = rec3_fmaf( w0.x,  x0z0,    g0.x,  z0,    w3.x,  x1z0);
    float row2 = rec3_fmaf( w0.x,  tex010,  g0.x,  y1z0,  w3.x,  tex110);
    // weighting along y direction
    float Z0 = rec3_fmaf( w0.y, row0, g0.y, row1, w3.y, row2);
    // slice 3 (z=1), weighing along x direction
    row0 = rec3_fmaf( w0.x, tex001, g0.x, y0z1, w3.x, tex101);
    row1 = rec3_fmaf( w0.x, x0z1,   g0.x, z1,   w3.x, x1z1);
    row2 = rec3_fmaf( w0.x, tex011, g0.x, y1z1, w3.x, tex111);
    // weighting along y direction
    float Z2 = rec3_fmaf( w0.y, row0, g0.y, row1, w3.y, row2);

    // slice 2 (z in middle, 4 bilinear, 4 linear and 1 trilinear lookup), weighing along x-direction
    row0 = rec3_fmaf( w0.x, x0y0, g0.x, y0, w3.x, x1y0);
    row1 = rec3_fmaf( w0.x, x0, g0.x, core, w3.x, x1);
    row2 = rec3_fmaf( w0.x, x1y0, g0.x, y1, w3.x, x1y1);
    // weighting along y direction
    float Z1 = rec3_fmaf( w0.y, row0, g0.y, row1, w3.y, row2);
    
    // weighting along z-direction
    return rec3_fmaf( w0.z, Z0, g0.z, Z1, w3.z, Z2);
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
 * @brief device function for computing the linear index from given 3D indices
 *******************************************************************/
__device__ int getLinearIdxfrom3DCoord(int x, int y, int z, int width, int height) {
    
    // width will be the pitch in case of pitched memory
    return  x +  width*y + width*height*z;
    
}

/********************************************************************
 * @brief device functions for computing the true velocity and function 
 * values at given coordinates for debugging purposes
 *******************************************************************/
__device__ float computeVx(float x, float y, float z) {
    return cosf(y)*cosf(z);
}

__device__ float computeVy(float x, float y, float z) {
    return sinf(x)*sinf(z);
}

__device__ float computeVz(float x, float y, float z) {
    return cosf(x)*cosf(y);
}

__device__ float computeM(float x, float y, float z) {
    return (sinf(x)*sinf(x) + sinf(y)*sinf(y) + sinf(z)*sinf(z))/3.0f;
}


/********************************************************************
 * @brief kernel function for computing the linear index from given 3D indices
 *******************************************************************/
__global__ void compareMemoryWithTexture( float* m, cudaPitchedPtr mptr, cudaTextureObject_t mtex, const float3 nx, const float3 inv_extent) {
    const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
    const int tidy = blockIdx.y * blockDim.y + threadIdx.y;
    const int tidz = blockIdx.z * blockDim.z + threadIdx.z;
    const int tid = tidx*nx.y*nx.z + tidy*nx.z + tidz;
    
    float3 q = make_float3(tidx+0.5, tidy+0.5, tidz+0.5);
    float3 qcoord = make_float3(tidx, tidy, tidz);
    qcoord *= 2*PI*inv_extent.x;
    q = q*inv_extent;
    float mval = m[tid];
    char* mptr1 = (char*)mptr.ptr;
    size_t pitch = mptr.pitch;
    size_t slicePitch = pitch*nx.y;
    char* slice = mptr1 + tidx*slicePitch;
    float* depth = (float*)(slice + tidy*pitch);
    float mptrval = depth[tidz]; 
    float mtexval = tex3D<float>(mtex, q.z, q.y, q.x);
    float mtrue = computeVx(qcoord.x, qcoord.y, qcoord.z);
    if (tid>=1 && tid<20) {
        printf("[%d,%d,%d] \t tid = %d \t  mval = %f \t mptrval = %f \t mtexval = %f\t mtrue = %f\n", tidx, tidy, tidz, tid, mval, mptrval, mtexval, mtrue);
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
        const float3 inv_nx)
{
    // Get thread index 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    float3 qcoord = make_float3(zq[tid], yq[tid], xq[tid]);
    // do single point interpolation - 4 methods

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



/********************************************************************
 * @brief host function to do interpolation of a scalar field
 * @parm[in] yi input data values 
 * @parm[in] xq1,yq1,zq1 query coordinates
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
           int*  nx,
           float* interp_time)
{
   
    // timing variables
    float time=0, dummy_time=0;
    cudaEvent_t startEvent, stopEvent;
    cudaEventCreate(&startEvent);
    cudaEventCreate(&stopEvent);

    // define inv of nx for normalizing in texture interpolation
    const float3 inv_nx = make_float3(  1.0f/static_cast<float>(nx[2]),
                                        1.0f/static_cast<float>(nx[1]), 
                                        1.0f/static_cast<float>(nx[0]));
    // define nxq, the dimensions of the grid
    const float3 nxq = make_float3( nx[0], nx[1], nx[2]);
    long int nq = nx[0]*nx[1]*nx[2]; 

    // create a common cudaResourceDesc objects
    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
   
    // create cudaPitchedPointers for query points by copying them to another device location
    //cudaPitchedPtr xq, yq, zq;
    //xq = CopyVolumeDeviceToDevice(xq1, nq, 1, 1);
    //yq = CopyVolumeDeviceToDevice(xq2, nq, 1, 1);
    //zq = CopyVolumeDeviceToDevice(xq3, nq, 1, 1);
    //////////////////////////////////////////////////////////////////////////////////////////////////    
    
   

    // make input image a cudaPitchedPtr for fi
    cudaPitchedPtr yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(yi), nx[2]*sizeof(float), nx[2], nx[1]);
    // initiate by computing the bspline coefficients for mt (in-place computation, updates mt)
    CubicBSplinePrefilter3D_Periodic((float*)yi_cudaPitchedPtr.ptr, (uint)yi_cudaPitchedPtr.pitch, nx[2], nx[1], nx[0]);
    // create a cudaExtent for input resolution
    cudaExtent yi_extent = make_cudaExtent(nx[2], nx[1], nx[0]);
    // create a texture from the spline coefficients
    cudaTextureObject_t yi_tex = initTextureFromVolume(yi_cudaPitchedPtr,  yi_extent);

    int threads = 256;
    int blocks = nq/threads;
    
    // check the correctness of the input data by checking the consistency across
    // the linear memory, pitched pointer and the texture
/*	dim3 dimBlock(1, 16, 16);
	dim3 dimGrid(nx[0] / dimBlock.x, nx[1] / dimBlock.y, nx[2] / dimBlock.z);
    compareMemoryWithTexture<<<dimBlock, dimGrid>>>(yi, yi_cudaPitchedPtr, yi_tex, nxq, inv_nx);
    if ( cudaSuccess != cudaGetLastError())
        printf("Error in running the interp3D kernel\n");
    cudaDeviceSynchronize(); 
    printf("\n----------------------------------------------------------------------------------------------------------\n");
*/

    // start recording the interpolation kernel
    time = 0; dummy_time = 0; 
    cudaEventRecord(startEvent,0); 
    
    // launch the interpolation kernel
    interp3D_kernel<<<blocks,threads>>>(yi_tex, xq1, xq2, xq3, yo, inv_nx);
    //interp3D_kernel<<<blocks,threads>>>(yi_tex, (PetscScalar*)xq.ptr, (PetscScalar*)yq.ptr, (PetscScalar*)zq.ptr, yo, inv_nx);
    if ( cudaSuccess != cudaGetLastError())
        printf("Error in running the interp3D kernel\n");

    cudaEventRecord(stopEvent,0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dummy_time, startEvent, stopEvent);
    time+=dummy_time;
    cudaDeviceSynchronize();
    
    // free texture and cudaArray from device memory
    cudaGetTextureObjectResourceDesc( &resDesc, yi_tex);
    cudaDestroyTextureObject(yi_tex);
    cudaFreeArray( resDesc.res.array.array);
    cudaEventDestroy(startEvent);
    cudaEventDestroy(stopEvent);
    
    printf("interp time = %fmsec\n", time);
    *interp_time += time;
    
}


/********************************************************************
 * @brief kernel function to do get the initial condition i.e grid indices for SemiLagrangian 
 * @parm[out] x,y,z memory space for storing coordinates
 * @parm[in] nxq dimensions for the grid
 *******************************************************************/
__global__ void getSMLInitialCondition_kernel(PetscScalar* x, PetscScalar* y, PetscScalar* z, const int3 nx) {
     int tidx = blockIdx.x * blockDim.x + threadIdx.x;
     int tidy = blockIdx.y * blockDim.y + threadIdx.y;
     int tidz = blockIdx.z * blockDim.z + threadIdx.z;
     int tid = tidx*nx.y*nx.z + tidy*nx.z + tidz;

    x[tid] = static_cast<PetscScalar>(tidx);
    y[tid] = static_cast<PetscScalar>(tidy);
    z[tid] = static_cast<PetscScalar>(tidz);

    if (tid>=4094 && tid<4098) {
        printf("tid = %d \t x = %0.1f \t y = %0.1f \t z = %0.1f\n", tid, x[tid], y[tid], z[tid]);
    }
    
}


/********************************************************************
 * @brief host function to do get the initial condition i.e grid indices for SemiLagrangian 
 * @parm[out] x,y,z memory space for storing coordinates
 * @parm[out] compute_time time required to do the initialization
 *******************************************************************/
void getSemiLagrangianInitialCondition(PetscScalar* x, PetscScalar* y, PetscScalar* z, int* nx, PetscScalar* compute_time) {
    float time=0, dummy_time=0;
    cudaEvent_t startEvent, stopEvent;
    cudaEventCreate(&startEvent);
    cudaEventCreate(&stopEvent);
    
    const int3 nxq = make_int3(nx[0], nx[1], nx[2]);

    dim3 dimBlock(1,16,16);
	dim3 dimGrid( nx[0] / dimBlock.x, nx[1] / dimBlock.y, nx[2] / dimBlock.z );
    
/*	uint dimX = min(min(PowTwoDivider(nx[0]), PowTwoDivider(nx[1])), 64);
	uint dimY = min(min(PowTwoDivider(nx[2]), PowTwoDivider(nx[1])), 512/dimX);
	dim3 dimBlock(dimX, dimY);
    
	// Replace the voxel values by the b-spline coefficients
	dim3 dimGrid(nx[1] / dimBlock.x, nx[2] / dimBlock.y);
*/    
    printf("blockx = %d, blocky = %d, gridx - %d, gridy = %d\n", dimBlock.x, dimBlock.y, dimGrid.x, dimGrid.y);
  getSMLInitialCondition_kernel<<<dimBlock, dimGrid>>>(x, y, z, nxq);
    
    if ( cudaSuccess != cudaGetLastError())
        printf("Error in running the SML initial condition kernel\n");
    
    cudaEventRecord(stopEvent,0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dummy_time, startEvent, stopEvent);
    time+=dummy_time/1000;
    cudaDeviceSynchronize();

    *compute_time += time;
}
