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
}


template <typename T>
__host__ __device__
inline T rec4_fmaf(T a, T b, T c, T d, T e, T f, T g, T h) {
    return fmaf(a, b, fmaf(c, d, fmaf( e, f, g*h)));
}


__device__ float cubicTex3D_lagrangeFast( cudaTextureObject_t tex, const float3 coord, const float3 inv_ext)
{
	const float3 coord_grid = coord;
	const float3 index = floor(coord_grid);
	const float3 fraction = coord_grid - index;
	float3 w0, w1, w2, w3;
	lagrange_weights(fraction, w0, w1, w2, w3);

    // compute the locations for the trilinear, bilinear and linear interps
    const float3 g0 = w1 + w2;
    const float3 h0 = ((w2/g0) + index + 0.5f)*inv_ext;
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


extern "C" cudaTextureObject_t initTextureFrom1DVector(cudaPitchedPtr arr, uint length) {
    cudaError_t err = cudaSuccess;
    /* create cuda resource description */
    struct cudaResourceDesc resDesc;
    memset( &resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = (float*)arr.ptr;
    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 32;
    resDesc.res.linear.desc.w = 32;
    resDesc.res.linear.sizeInBytes = length*sizeof(float);

    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.readMode = cudaReadModeElementType;

    cudaTextureObject_t texObj = 0;
    err = cudaCreateTextureObject( &texObj, &resDesc, &texDesc, NULL);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to create texture (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    return texObj;

}


__device__ int getLinearIdxfrom3DCoord(int x, int y, int z, int width, int height) {
    
    return  x +  width*y + width*height*z;
    
}


/*
 * Computes the Euler departure point using simple euler integration
 * @param[in] [vx,vy,vz] velocity vector field as 3D pitched pointers
 * @param[in] extent length of the regular grid in each dimension
*/
__global__ void getEulerPoint(
    cudaPitchedPtr vx,
    cudaPitchedPtr vy,
    cudaPitchedPtr vz,
    int3 extent,          
    float3 minlim,        
    float3 h, float dt, 
    cudaPitchedPtr xstar,
    cudaPitchedPtr x_sml) 
{
    // 3D grid of 3D block
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;
 
    // get the linear index of the 3D array
    int tid_vx = getLinearIdxfrom3DCoord( x, y, z, (int)vx.pitch, extent.y);
    int tid_vy = getLinearIdxfrom3DCoord( x, y, z, (int)vy.pitch, extent.y);
    int tid_vz = getLinearIdxfrom3DCoord( x, y, z, (int)vz.pitch, extent.y);
    int tid_xstar = getLinearIdxfrom3DCoord( x, y, z, (int)xstar.pitch, extent.y);
    int tid_xsml = getLinearIdxfrom3DCoord( x, y, z, (int)x_sml.pitch, extent.y); 

    float* vx1 = (float*) vx.ptr;
    float* vy1 = (float*) vy.ptr;
    float* vz1 = (float*) vz.ptr;
    float3 v = make_float3( *(vx1 + tid_vx), *(vy1 + tid_vy), *(vz1 + tid_vz));

    // get the regular grid coordinate
    float3 id = make_float3( (float)x, (float)y, (float)z);
    float3 x0 = minlim + h*id;
    
    // compute the departure point using euler
    *((float3*)xstar.ptr + tid_xstar) = x0 - dt*v;
    // also partially compute the semi-lagrangian point for RK2 scheme
    *((float3*)x_sml.ptr + tid_xsml) = x0 - dt*v/2.0f;

}   


/****************************************************************************/
__global__ void getSmlCoord(
    cudaTextureObject_t vxtex,
    cudaTextureObject_t vytex,
    cudaTextureObject_t vztex,
    cudaPitchedPtr xstar,
    float3 extent,                  
    float3 h, float dt, 
    cudaPitchedPtr x_sml) 
{
    // 3D grid of 3D block
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;
 
    // get the linear index of the 3D array
    int tid_xstar = getLinearIdxfrom3DCoord( x, y, z, (int)xstar.pitch, (int)extent.y);
    int tid_xsml = getLinearIdxfrom3DCoord( x, y, z, (int)x_sml.pitch, (int)extent.y); 

    // get the points where you want to interpolate
    float3 xq = *((float3*)xstar.ptr + tid_xstar);
    
    float3 v_interp;
    // interpolate ***************************** correct extent
    v_interp.x = cubicTex3D_splineFast(vxtex, xq, extent);
    v_interp.y = cubicTex3D_splineFast(vytex, xq, extent);
    v_interp.z = cubicTex3D_splineFast(vztex, xq, extent);
    // compute the semi-lagrangian point for RK2 scheme
    *((float3*)x_sml.ptr + tid_xsml) -=  dt*v_interp/2.0f;

}   


__global__ void interp3D_kernel(
        cudaTextureObject_t  yi_tex,
        const PetscScalar* xq,
        const PetscScalar* yq,
        const PetscScalar* zq, 
        PetscScalar* yo,
        const float3 inv_nx)
{
    // 3D grid of 3D block
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    float3 qcoord = make_float3(xq[tid], yq[tid], zq[tid]);
    
    // do single point interpolation
    yo[tid] = cubicTex3D_splineFast(yi_tex, qcoord, inv_nx);
    //mt_interp[mtid] = cubicTex3D_splineSimple( mt_tex, qcoord_normalised, inv_reg_extent);
    //mt_interp[mtid] = cubicTex3D_lagrangeSimple( mt_tex, qcoord_normalised, inv_reg_extent);
    // mt_interp[tid] = cubicTex3D_lagrangeFast( mt_tex, qcoord_normalised, inv_reg_extent);
}


__global__ void printSliceFromVolume(
        float* volume,
        int pitch,
        int height,
        int slice)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    
    int tid = getLinearIdxfrom3DCoord( x, y, slice, pitch, height);

    printf("x = %d, y = %d, id = %d, val = %f\n", x, y, tid, volume[tid]);
}

/*
static void printFloat3(float3 var) {
    printf("x = %f \t y = %f \t z = %f\n", var.x, var.y, var.z);
}
*/

/**********************************************************************************/
void gpuInterp3D(
           PetscScalar* yi,
           const PetscScalar* xq1,
           const PetscScalar* xq2,
           const PetscScalar* xq3,
           PetscScalar* yo,
           int*  nx)
{
    //cudaError_t err = cudaSuccess;
    
    // timing variables
    float time=0, dummy_time=0;
    cudaEvent_t startEvent, stopEvent;
    cudaEventCreate(&startEvent);
    cudaEventCreate(&stopEvent);

    // define inv of nx for normalizing in texture interpolation
    const float3 inv_nx = make_float3(  1.0f/static_cast<float>(nx[0]),
                                        1.0f/static_cast<float>(nx[1]), 
                                        1.0f/static_cast<float>(nx[2]));
    
    // create a common cudaResourceDesc objects
    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    
    // make input image a cudaPitchedPtr for fi
    cudaPitchedPtr yi_cudaPitchedPtr = make_cudaPitchedPtr(static_cast<void*>(yi), nx[0]*sizeof(float), nx[0], nx[1]);
    // initiate by computing the bspline coefficients for mt (in-place computation, updates mt)
    // CubicBSplinePrefilter3D_Periodic((float*)yi_cudaPitchedPtr.ptr, (uint)yi_cudaPitchedPtr.pitch, nx[0], ny[1], nx[2]);
    // create a cudaExtent for input resolution
    cudaExtent yi_extent = make_cudaExtent(nx[0], nx[1], nx[2]);
    // create a texture from the spline coefficients
    cudaTextureObject_t yi_tex = initTextureFromVolume(yi_cudaPitchedPtr,  yi_extent);

    // start recording the interpolation kernel
    time = 0; dummy_time = 0; 
    cudaEventRecord(startEvent,0); 
    int threads = 256;
    long int nq = nx[0]*nx[1]*nx[2];
    int blocks = nq/threads;
    // launch the kernel
    interp3D_kernel<<<blocks,threads>>>(yi_tex, xq1, xq2, xq3, yo, inv_nx); 
    
    if ( cudaSuccess != cudaGetLastError())
        printf("Error in running the interp3D kernel\n");

    cudaEventRecord(stopEvent,0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dummy_time, startEvent, stopEvent);
    time+=dummy_time/1000;
    cudaDeviceSynchronize();
    // print time
    printf("\n 3D interpolation of Q=%d query points on a grid N=%dx%dx%d took %0.2E sec\n\n", nx[0]*nx[1]*nx[2], nx[0], nx[1], nx[2], time);
    
    // free texture and cudaArray from device memory
    cudaGetTextureObjectResourceDesc( &resDesc, yi_tex);
    cudaDestroyTextureObject(yi_tex);
    cudaFreeArray( resDesc.res.array.array);
    cudaEventDestroy(startEvent);
    cudaEventDestroy(stopEvent);
    
}

