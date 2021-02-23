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
#define spencil 4
#define lpencil 32
#define sharedrows 32
#define perthreadcomp 8


const float h_c[HALO] = {4.f / 5.f , -1.f / 5.f , 4.f / 105.f, -1.f / 280.f};
const float h2_c[HALO+1] = {-205.f/72.f, 8.f / 5.f , -1.f / 5.f , 8.f / 315.f, -1.f / 560.f};


// device constants
__constant__ int d_nx, d_ny, d_nz;
__constant__ int d_isize0, d_isize1, d_isize2;
__constant__ int d_isize_g0, d_isize_g1, d_isize_g2;
__constant__ int d_halox, d_haloy, d_haloz;

__constant__ float d_invnx, d_invny, d_invnz;
__constant__ float d_invhx, d_invhy, d_invhz;
__constant__ float d_cx[HALO], d_cy[HALO], d_cz[HALO];
__constant__ float d_cxx[HALO+1], d_cyy[HALO+1], d_czz[HALO+1];

const int sx = spencil;
const int sy = sharedrows;
const int sxx = lpencil;
const int syy = sharedrows;



__device__ inline int getLinearIdx(int i, int j, int k, const dim3 nl) {
    return i*nl.y*nl.z + j*nl.z + k;
}

__device__ inline int getLinearIdx(int i, int j, int k) {
    return i*d_ny*d_nz + j*d_nz + k;
}

__device__ inline int getLinearIdx_ghost(int i, int j, int k) {
    return i*d_isize_g1*d_isize_g2 + j*d_isize_g2 + k;
}

__device__ inline int getLinearIdx_local(int i, int j, int k) {
    return i*d_isize1*d_isize2 + j*d_isize2 + k;
}

/**********************************************************************************
 * @brief compute z-gradient using 8th order finite differencing
 * @param[in] f = input scalar field with ghost padding
 * @param[out] dfz = z-component of gradient
**********************************************************************************/
__global__ void mgpu_gradient_z(ScalarType* dfz, const ScalarType* f) {
  __shared__ float s_f[sx][sy+2*HALO]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i = blockIdx.z;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x + HALO;       // local k for shared memory ac3ess + halo offset
  int sj = threadIdx.y; // local j for shared memory ac3ess
  int zblock_width;
  int id,lid,rid;
  
  if (blockIdx.x < gridDim.x - 1) {
    zblock_width = blockDim.x;
  }
  else {
    zblock_width = d_isize2 - blockIdx.x*blockDim.x;
  }
  
  bool internal = (j < d_isize1) && (threadIdx.x < zblock_width);
  
  if (internal) {
    id = getLinearIdx_ghost(i+d_halox,j+d_haloy,k+d_haloz);
    s_f[sj][sk] = f[id];
  }
  
  __syncthreads();
    
  // fill in periodic images in shared memory array 
  if (threadIdx.x < HALO) {
    if (d_haloz == 0) {
      lid = k%d_isize2-HALO;
      if (lid<0) lid+=d_isize2;
      id = getLinearIdx_ghost(i+d_halox, j+d_haloy, lid);
      s_f[sj][sk-HALO] = f[id];
      rid = (k+zblock_width)%d_isize2;
      id = getLinearIdx_ghost(i+d_halox, j+d_haloy, rid);
      s_f[sj][zblock_width+sk] = f[id];
    } else {
      id = getLinearIdx_ghost(i+d_halox, j+d_haloy, k);
      s_f[sj][sk-HALO] = f[id];
      id = getLinearIdx_ghost(i+d_halox, j+d_haloy, k+zblock_width+d_haloz);
      s_f[sj][zblock_width+sk] = f[id];
    }
  }
  
  __syncthreads();
  
  ScalarType result = 0;
  if (internal) {
    id = getLinearIdx_local(i,j,k);
    for(int l=0; l<HALO; l++) {
      result += d_cz[l] * (s_f[sj][sk+1+l] - s_f[sj][sk-1-l]);
    }
    dfz[id] += result;
  }
}



/**********************************************************************************
 * @brief compute y-component of gradient using 8th order finite differencing
 * @param[in]    f = scalar field with ghost layer padding
 * @param[out] dfy = y-component of gradient
**********************************************************************************/
__global__ void mgpu_gradient_y(ScalarType* dfy, const ScalarType* f) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
  int yblock_width, id, lid, rid, sj;
  bool internal;
  
  if ( blockIdx.y < gridDim.y - 1) {
    yblock_width = syy;
  }
  else {
    yblock_width = d_isize1 - syy*blockIdx.y;
  }
    
  for(int j = threadIdx.y; j < yblock_width; j += blockDim.y) {
    internal = ((blockIdx.y*syy+j) < d_isize1) && (k < d_isize2);
    if (internal) {
      id = getLinearIdx_ghost(i+d_halox, blockIdx.y*syy + j + d_haloy, k+d_haloz);
      sj = j + HALO;
      s_f[sj][sk] = f[id];
    }
  }
  
  __syncthreads();

  
  sj = threadIdx.y + HALO;
  int y = syy*blockIdx.y + threadIdx.y;
  // fill in periodic images in shared memory array 
  if (threadIdx.y < HALO) {
    if (d_haloy == 0) {
      lid = y%d_isize1-HALO;
      if (lid<0) lid+=d_isize1;
      id = getLinearIdx_ghost(i+d_halox, lid, k+d_haloz);
      s_f[sj-HALO][sk] = f[id];
      rid = (y+yblock_width)%d_isize1;
      id = getLinearIdx_ghost(i+d_halox, rid, k+d_haloz);
      s_f[sj+yblock_width][sk] = f[id];
    } else {
      id = getLinearIdx_ghost(i+d_halox, y, k+d_haloz);
      s_f[sj-HALO][sk] = f[id];
      id = getLinearIdx_ghost(i+d_halox, y+yblock_width+d_haloy, k+d_haloz);
      s_f[sj+yblock_width][sk] = f[id];
    }
  }

  __syncthreads();
  
  ScalarType result;
  for(int j = threadIdx.y; j < yblock_width; j += blockDim.y) {
    internal = ((blockIdx.y*syy+j) < d_isize1) && (k < d_isize2);
    result = 0;
    if (internal) {
      int id = getLinearIdx_local(i, blockIdx.y*syy + j ,k);
      int sj = j + HALO;
      for( int l=0; l<HALO; l++) {
        result += d_cy[l] * ( s_f[sj+1+l][sk] - s_f[sj-1-l][sk]);
      }
      dfy[id] += result;
    }
  }
}



/**********************************************************************************
 * @brief compute x-component of gradient using 8th order finite differencing
 * @param[in]    f = scalar field with ghost layer padding
 * @param[out] dfx = x-component of gradient
**********************************************************************************/
__global__ void mgpu_gradient_x(ScalarType* dfx, const ScalarType* f) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int j  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
    
  int id, lid, rid, si;
  int xblock_width;
  bool internal;
    
  if ( blockIdx.y < gridDim.y - 1) {
    xblock_width = syy;
  }
  else {
    xblock_width = d_isize0 - syy*blockIdx.y;
  }
  
  for(int i = threadIdx.y; i < xblock_width; i += blockDim.y) {
    internal = ((blockIdx.y*syy + i) < d_nx) && (k < d_nz);
    if (internal) {
      id = getLinearIdx_ghost(blockIdx.y*syy + i + d_halox, j + d_haloy, k + d_haloz);
      si = i + HALO;
      s_f[si][sk] = f[id];
    }
  }

  __syncthreads();

  
  si = threadIdx.y + HALO;
  int x = syy*blockIdx.y + threadIdx.y;
  // fill in periodic images in shared memory array 
  if (threadIdx.y < HALO) {
    if (d_halox == 0) {
      lid = x%d_isize0-HALO;
      if (lid<0) lid+=d_isize0;
      id = getLinearIdx_ghost(lid, j+d_haloy, k+d_haloz);
      s_f[si-HALO][sk] = f[id];
      rid = (x+xblock_width)%d_isize0;
      id = getLinearIdx_ghost(rid, j+d_haloy, k+d_haloz);
      s_f[si+xblock_width][sk] = f[id];
    } else {
      id = getLinearIdx_ghost(x, j+d_haloy, k+d_haloz);
      s_f[si-HALO][sk] = f[id];
      id = getLinearIdx_ghost(x+xblock_width+d_halox, j+d_haloy, k+d_haloz);
      s_f[si+xblock_width][sk] = f[id];
    }
  }

  __syncthreads();

  ScalarType result;
  for(int i = threadIdx.y; i < xblock_width; i += blockDim.y) {
    internal = ((blockIdx.y*syy + i) < d_isize0) && (k < d_isize2);
    result = 0;
    if (internal) {
      id = getLinearIdx_local(blockIdx.y*syy + i , j, k);
      si = i + HALO;
      for( int l=0; l<HALO; l++) {
          result +=  d_cx[l] * ( s_f[si+1+l][sk] - s_f[si-1-l][sk]);
      }
      dfx[id] += result;
    }
  }
}

/**********************************************************************************
 * @brief compute laplacian using 8th order finite differencing
 * @param[out]   ddf  = laplacian of scalar field f 
 * @param[in]    f    = scalar field f with ghost padding
 * @param[in]    beta = some constant which needs to be defined TODO
**********************************************************************************/
__global__ void mgpu_d_zz(ScalarType* ddf, const ScalarType* f, const ScalarType beta) {
  __shared__ float s_f[sx][sy+2*HALO]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i = blockIdx.z;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x + HALO;       // local k for shared memory ac3ess + halo offset
  int sj = threadIdx.y; // local j for shared memory ac3ess
  int zblock_width;
  int id,lid,rid;
  
  if (blockIdx.x < gridDim.x - 1) {
    zblock_width = blockDim.x;
  }
  else {
    zblock_width = d_isize2 - blockIdx.x*blockDim.x;
  }
  
  bool internal = (j < d_isize1) && (threadIdx.x < zblock_width);
  
  if (internal) {
    id = getLinearIdx_ghost(i+d_halox,j+d_haloy,k+d_haloz);
    s_f[sj][sk] = f[id];
  }
  
  __syncthreads();
    
  // fill in periodic images in shared memory array 
  if (threadIdx.x < HALO) {
    if (d_haloz == 0) {
      lid = k%d_isize2-HALO;
      if (lid<0) lid+=d_isize2;
      id = getLinearIdx_ghost(i+d_halox, j+d_haloy, lid);
      s_f[sj][sk-HALO] = f[id];
      rid = (k+zblock_width)%d_isize2;
      id = getLinearIdx_ghost(i+d_halox, j+d_haloy, rid);
      s_f[sj][zblock_width+sk] = f[id];
    } else {
      id = getLinearIdx_ghost(i+d_halox, j+d_haloy, k);
      s_f[sj][sk-HALO] = f[id];
      id = getLinearIdx_ghost(i+d_halox, j+d_haloy, k+zblock_width+d_haloz);
      s_f[sj][zblock_width+sk] = f[id];
    }
  }
  
  __syncthreads();
  
  ScalarType lval = d_czz[0]*s_f[sj][sk];
  if (internal) {
    id = getLinearIdx_local(i,j,k);
    for(int l=0; l<HALO; l++) {
      lval += d_czz[l] * (s_f[sj][sk+l] + s_f[sj][sk-l]);
    }
    ddf[id] = lval;
  }
}

/**********************************************************************************
 * @brief compute laplacian using 8th order finite differencing
 * @param[inout] ddf  = partial laplacian of f
 * @param[in]    f    = scalar field f with ghost padding
 * @param[in]    beta = some constant which needs to be defined TODO
**********************************************************************************/
__global__ void mgpu_d_yy(ScalarType* ddf, const ScalarType* f, const ScalarType beta) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
  int yblock_width, id, lid, rid, sj;
  bool internal;
  
  if ( blockIdx.y < gridDim.y - 1) {
    yblock_width = syy;
  }
  else {
    yblock_width = d_isize1 - syy*blockIdx.y;
  }
    
  for(int j = threadIdx.y; j < yblock_width; j += blockDim.y) {
    internal = ((blockIdx.y*syy+j) < d_isize1) && (k < d_isize2);
    if (internal) {
      id = getLinearIdx_ghost(i+d_halox, blockIdx.y*syy + j + d_haloy, k+d_haloz);
      sj = j + HALO;
      s_f[sj][sk] = f[id];
    }
  }
  
  __syncthreads();

  
  sj = threadIdx.y + HALO;
  int y = syy*blockIdx.y + threadIdx.y;
  // fill in periodic images in shared memory array 
  if (threadIdx.y < HALO) {
    if (d_haloy == 0) {
      lid = y%d_isize1-HALO;
      if (lid<0) lid+=d_isize1;
      id = getLinearIdx_ghost(i+d_halox, lid, k+d_haloz);
      s_f[sj-HALO][sk] = f[id];
      rid = (y+yblock_width)%d_isize1;
      id = getLinearIdx_ghost(i+d_halox, rid, k+d_haloz);
      s_f[sj+yblock_width][sk] = f[id];
    } else {
      id = getLinearIdx_ghost(i+d_halox, y, k+d_haloz);
      s_f[sj-HALO][sk] = f[id];
      id = getLinearIdx_ghost(i+d_halox, y+yblock_width+d_haloy, k+d_haloz);
      s_f[sj+yblock_width][sk] = f[id];
    }
  }

  __syncthreads();
  
  for(int j = threadIdx.y; j < yblock_width; j += blockDim.y) {
    internal = ((blockIdx.y*syy+j) < d_isize1) && (k < d_isize2);
    if (internal) {
      ScalarType lval = d_cyy[0]*s_f[sj][sk];
      int id = getLinearIdx_local(i, blockIdx.y*syy + j ,k);
      int sj = j + HALO;
      for( int l=0; l<HALO; l++) {
        lval += d_cyy[l] * ( s_f[sj+l][sk] + s_f[sj-l][sk]);
      }
      ddf[id] += lval;
    }
  }
}

/**********************************************************************************
 * @brief compute laplacian using 8th order finite differencing
 * @param[inout] ddf  = parital laplacian of scalar field f
 * @param[in]    f    = scalar field f with ghost padding
 * @param[in]    beta = some constant which needs to be defined TODO
**********************************************************************************/
__global__ void mgpu_d_xx(ScalarType* ddf, const ScalarType* f, const ScalarType beta) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int j  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
    
  int id, lid, rid, si;
  int xblock_width;
  bool internal;
    
  if ( blockIdx.y < gridDim.y - 1) {
    xblock_width = syy;
  }
  else {
    xblock_width = d_isize0 - syy*blockIdx.y;
  }
  
  for(int i = threadIdx.y; i < xblock_width; i += blockDim.y) {
    internal = ((blockIdx.y*syy + i) < d_nx) && (k < d_nz);
    if (internal) {
      id = getLinearIdx_ghost(blockIdx.y*syy + i + d_halox, j + d_haloy, k + d_haloz);
      si = i + HALO;
      s_f[si][sk] = f[id];
    }
  }

  __syncthreads();

  
  si = threadIdx.y + HALO;
  int x = syy*blockIdx.y + threadIdx.y;
  // fill in periodic images in shared memory array 
  if (threadIdx.y < HALO) {
    if (d_halox == 0) {
      lid = x%d_isize0-HALO;
      if (lid<0) lid+=d_isize0;
      id = getLinearIdx_ghost(lid, j+d_haloy, k+d_haloz);
      s_f[si-HALO][sk] = f[id];
      rid = (x+xblock_width)%d_isize0;
      id = getLinearIdx_ghost(rid, j+d_haloy, k+d_haloz);
      s_f[si+xblock_width][sk] = f[id];
    } else {
      id = getLinearIdx_ghost(x, j+d_haloy, k+d_haloz);
      s_f[si-HALO][sk] = f[id];
      id = getLinearIdx_ghost(x+xblock_width+d_halox, j+d_haloy, k+d_haloz);
      s_f[si+xblock_width][sk] = f[id];
    }
  }

  __syncthreads();

  for(int i = threadIdx.y; i < xblock_width; i += blockDim.y) {
    internal = ((blockIdx.y*syy + i) < d_isize0) && (k < d_isize2);
    if (internal) {
      id = getLinearIdx_local(blockIdx.y*syy + i , j, k);
      si = i + HALO;
      ScalarType lval = d_cxx[0]*s_f[si][sk];
      for( int l=0; l<HALO; l++) {
          lval += d_cxx[l] * ( s_f[si+l][sk] + s_f[si-l][sk]);
      }
      ddf[id] += lval;
      ddf[id] *= beta;
    }
  }
}

inline __device__ void add_op (ScalarType& a, ScalarType b) { a += b; }
inline __device__ void replace_op (ScalarType& a, ScalarType b) { a = b; }

/**************************************************************************************************
 * @brief compute z-component of gradient using 8th order finite differencing (single GPU version)
 * @param[in]    f = scalar field with no ghost layer padding
 * @param[out] dfz = z-component of gradient
**************************************************************************************************/
template<void(*Op)(ScalarType&,ScalarType)=add_op>
__global__ void gradient_z(ScalarType* dfz, const ScalarType* f, const dim3 nl, const float ih) {
  __shared__ float s_f[sx][sy+2*HALO]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i = blockIdx.z;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x + HALO;       // local k for shared memory ac3ess + halo offset
  int sj = threadIdx.y; // local j for shared memory ac3ess
  int zblock_width, id;
  
  const float cz[HALO] = {4.f / 5.f , -1.f / 5.f , 4.f / 105.f, -1.f / 280.f};
  
  if (blockIdx.x < gridDim.x - 1) {
    zblock_width = blockDim.x;
  }
  else {
    zblock_width = nl.z - blockIdx.x*blockDim.x;
  }
  
  bool internal = (j < nl.y) && (threadIdx.x < zblock_width);
  
  if (internal) {
    id = getLinearIdx(i,j,k);
    s_f[sj][sk] = f[id];
  }

  __syncthreads();
    
  int lid,rid;
  // fill in periodic images in shared memory array 
  if (threadIdx.x < HALO) {
    lid = k%nl.z-HALO;
    if (lid<0) lid+=nl.z;
    s_f[sj][sk-HALO] = f[i*nl.y*nl.z + j*nl.z + lid];
    rid = (k+zblock_width)%d_nz;
    s_f[sj][zblock_width+sk] = f[i*nl.y*nl.z + j*nl.z + rid];
  }

  __syncthreads();
  
  ScalarType result = 0;
  if (internal) {
    for(int l=0; l<HALO; l++) {
        result += cz[l] * (s_f[sj][sk+1+l] - s_f[sj][sk-1-l]);
    }
    Op(dfz[id], result*ih);
  }
}


/**************************************************************************************************
 * @brief compute y-component of gradient using 8th order finite differencing (single GPU version)
 * @param[in]    f = scalar field with no ghost layer padding
 * @param[out] dfy = y-component of gradient
**************************************************************************************************/
template<void(*Op)(ScalarType&,ScalarType)=add_op>
__global__ void gradient_y(ScalarType* dfy, const ScalarType* f, const dim3 nl, const float ih) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int i  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
    
  int yblock_width, globalIdx, sj;
  bool internal;
  
  const float cy[HALO] = {4.f / 5.f , -1.f / 5.f , 4.f / 105.f, -1.f / 280.f};
  
  if ( blockIdx.y < gridDim.y - 1) {
    yblock_width = syy;
  }
  else {
    yblock_width = nl.y - syy*blockIdx.y;
  }
  
    
  for(int j = threadIdx.y; j < yblock_width; j += blockDim.y) {
    internal = ((blockIdx.y*syy+j) < nl.y) && (k < nl.z);
    if (internal) {
        globalIdx = getLinearIdx(i, blockIdx.y*syy + j ,k);
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
      globalIdx = getLinearIdx(i, blockIdx.y*syy + j ,k);
      sj = j + HALO;
      for( int l=0; l<HALO; l++) {
          result +=  cy[l] * ( s_f[sj+1+l][sk] - s_f[sj-1-l][sk]);
      }
      Op(dfy[globalIdx], result*ih);
    }
  }
}

/**************************************************************************************************
 * @brief compute x-component of gradient using 8th order finite differencing (single GPU version)
 * @param[in]    f = scalar field with no ghost layer padding
 * @param[out] dfz = x-component of gradient
**************************************************************************************************/
template<void(*Op)(ScalarType&,ScalarType)=add_op>
__global__ void gradient_x(ScalarType* dfx, const ScalarType* f, const dim3 nl, const float ih) {
  __shared__ float s_f[syy+2*HALO][sxx]; // HALO-wide halo for central diferencing scheme
    
  // note i and k have been exchanged to ac3ount for k being the fastest changing index
  int j  = blockIdx.z;
  int k  = blockIdx.x*blockDim.x + threadIdx.x;
  int sk = threadIdx.x;       // local k for shared memory ac3ess, fixed
    
  int xblock_width, globalIdx, si;
  bool internal;
  
  const float cx[HALO] = {4.f / 5.f , -1.f / 5.f , 4.f / 105.f, -1.f / 280.f};
  
  if ( blockIdx.y < gridDim.y - 1) {
    xblock_width = syy;
  }
  else {
    xblock_width = nl.x - syy*blockIdx.y;
  }
    
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
      for( int l=0; l<HALO; l++) {
        result += cx[l] * ( s_f[si+1+l][sk] - s_f[si-1-l][sk]);
      }
      Op(dfx[globalIdx],result*ih);
    } 
  }
}

/**********************************************************************************
 * @brief compute laplacian using 8th order finite differencing (single GPU code)
 * @param[inout] dfz  = parital laplacian of scalar field f
 * @param[in]    f    = scalar field f
 * @param[in]    beta = some constant which needs to be defined TODO
**********************************************************************************/
__global__ void d_zz(ScalarType* ddf, const ScalarType* f, const ScalarType beta) {
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

    ScalarType lval = d_czz[0]*s_f[sj][sk];
    for(int l=1; l<=HALO; l++) {
      lval += d_czz[l] * (s_f[sj][sk+l] + s_f[sj][sk-l]);
    }
    ddf[globalIdx] = lval;
}

/**********************************************************************************
 * @brief compute laplacian using 8th order finite differencing (single GPU code)
 * @param[inout] dfy  = parital laplacian of scheme field f 
 * @param[in]    f    = scalar field f
 * @param[in]    beta = some constant which needs to be defined TODO
**********************************************************************************/
__global__ void d_yy(ScalarType* ddf, const ScalarType* f, const ScalarType beta) {
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
    ScalarType lval =d_cyy[0]*s_f[sj][sk];
    for( int l=1; l<=HALO; l++) {
        lval += d_cyy[l] * ( s_f[sj+l][sk] + s_f[sj-l][sk]);
    }
    ddf[globalIdx] += lval;
  }

}

/**********************************************************************************
 * @brief compute laplacian using 8th order finite differencing (single GPU code)
 * @param[inout] ddf  = partial laplacian of scalar field f
 * @param[in]    f    = scalar field f
 * @param[in]    beta = some constant which needs to be defined TODO
**********************************************************************************/
__global__ void d_xx(ScalarType* ddf, const ScalarType* f, const ScalarType beta) {
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
    ScalarType lval = d_cxx[0]*s_f[si][sk];
    for( int l=1; l<=HALO; l++) {
        lval += d_cxx[l] * ( s_f[si+l][sk] + s_f[si-l][sk]);
    }
    ddf[globalIdx] += lval;
    ddf[globalIdx] *= beta;
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

void printFloat3(float3 a){
    printf("x = %f\t y = %f\t z = %f\n",a.x,a.y,a.z);
}

namespace reg {

cudaTextureObject_t gpuInitEmptyGradientTexture(IntType *nx) {
   cudaTextureObject_t texObj = 0;
#if defined(USE_TEXTURE_GRADIENT)
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

    err = cudaCreateTextureObject( &texObj, &resDesc, &texDesc, NULL);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to create texture (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
#endif
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
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to copy 3D memory to cudaArray (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

PetscErrorCode initConstants(const IntType* iisize, const IntType* iisize_g, const ScalarType* hx, const IntType* ihalo) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  int halo[3], isize[3], isize_g[3];
  halo[0] = ihalo[0]; halo[1] = ihalo[1]; halo[2] = ihalo[2];
  isize[0] = iisize[0]; isize[1] = iisize[1]; isize[2] = iisize[2];
  isize_g[0] = iisize_g[0]; isize_g[1] = iisize_g[1]; isize_g[2] = iisize_g[2];
  
  float3 inv_nx = make_float3(  1.0f/static_cast<float>(isize[0]),
                                1.0f/static_cast<float>(isize[1]), 
                                1.0f/static_cast<float>(isize[2]));
  float3 inv_hx = make_float3(0.5f/hx[0], 0.5f/hx[1], 0.5f/hx[2]);

  cudaMemcpyToSymbol(d_halox, &halo[0], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_haloy, &halo[1], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_haloz, &halo[2], sizeof(int), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(d_isize0, &isize[0], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_isize1, &isize[1], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_isize2, &isize[2], sizeof(int), 0, cudaMemcpyHostToDevice);
  
  cudaMemcpyToSymbol(d_nx, &isize[0], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_ny, &isize[1], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_nz, &isize[2], sizeof(int), 0, cudaMemcpyHostToDevice);
  
  cudaMemcpyToSymbol(d_isize_g0, &isize_g[0], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_isize_g1, &isize_g[1], sizeof(int), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_isize_g2, &isize_g[2], sizeof(int), 0, cudaMemcpyHostToDevice);
  
  cudaMemcpyToSymbol(d_invnx, &inv_nx.x, sizeof(float), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_invny, &inv_nx.y, sizeof(float), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_invnz, &inv_nx.z, sizeof(float), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(d_invhx, &inv_hx.x, sizeof(float), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_invhy, &inv_hx.y, sizeof(float), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_invhz, &inv_hx.z, sizeof(float), 0, cudaMemcpyHostToDevice);
  
  float h_ct[HALO+1];
  for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx[0];
  cudaMemcpyToSymbol(d_cx, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);
  for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx[1];
  cudaMemcpyToSymbol(d_cy, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);
  for(int l=0; l<HALO; l++) h_ct[l] = h_c[l]/hx[2];
  cudaMemcpyToSymbol(d_cz, h_ct, sizeof(float)*HALO, 0, cudaMemcpyHostToDevice);
  
  for(int l=0; l<=HALO; l++) h_ct[l] = h2_c[l]/(hx[0]*hx[0]);
  cudaMemcpyToSymbol(d_cxx, h_ct, sizeof(float)*(HALO+1), 0, cudaMemcpyHostToDevice);
  for(int l=0; l<=HALO; l++) h_ct[l] = h2_c[l]/(hx[1]*hx[1]);
  cudaMemcpyToSymbol(d_cyy, h_ct, sizeof(float)*(HALO+1), 0, cudaMemcpyHostToDevice);
  for(int l=0; l<=HALO; l++) h_ct[l] = h2_c[l]/(hx[2]*hx[2]);
  cudaMemcpyToSymbol(d_czz, h_ct, sizeof(float)*(HALO+1), 0, cudaMemcpyHostToDevice);
  
  PetscFunctionReturn(ierr);
}

void getThreadBlockDimensionsX(dim3& threads, dim3& blocks, IntType* nx) {
  threads.x = sxx;
  threads.y = syy/perthreadcomp;
  threads.z = 1;
  blocks.x = (nx[2]+sxx-1)/sxx;
  blocks.y = (nx[0]+syy-1)/syy;
  blocks.z = nx[1];
}

void getThreadBlockDimensionsY(dim3& threads, dim3& blocks, IntType* nx) {
  threads.x = sxx;
  threads.y = syy/perthreadcomp;
  threads.z = 1;
  blocks.x = (nx[2]+sxx-1)/sxx;
  blocks.y = (nx[1]+syy-1)/syy;
  blocks.z = nx[0];
}

void getThreadBlockDimensionsZ(dim3& threads, dim3& blocks, IntType* nx) {
  threads.x = sy;
  threads.y = sx;
  threads.z = 1;
  blocks.x = (nx[2]+sy-1)/sy;
  blocks.y = (nx[1]+sx-1)/sx;
  blocks.z = nx[0];
}

PetscErrorCode computeDivergence(ScalarType* l, const ScalarType* g1, const ScalarType* g2, const ScalarType* g3, cudaTextureObject_t mtex, IntType* nx, ScalarType* hx, bool mgpu) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  size_t count = sizeof(ScalarType)*nx[0]*nx[1]*nx[2];
  if (mgpu) {
    ierr = cudaMemset((void*)l, 0, count); CHKERRCUDA(ierr);
  }

#if defined(USE_TEXTURE_GRADIENT)
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
#else
  const dim3 nl (nx[0], nx[1], nx[2]);
  // Shared memory implementation
  // Z-Gradient
  dim3 threadsPerBlock_z, numBlocks_z;
  getThreadBlockDimensionsZ(threadsPerBlock_z, numBlocks_z, nx);
  if (mgpu)
    mgpu_gradient_z<<<numBlocks_z, threadsPerBlock_z>>>(l,g3);
  else
    gradient_z<replace_op><<<numBlocks_z, threadsPerBlock_z>>>(l,g3,nl,1./hx[2]);
  cudaCheckKernelError();
    
  // Y-Gradient 
  dim3 threadsPerBlock_y, numBlocks_y;
  getThreadBlockDimensionsY(threadsPerBlock_y, numBlocks_y, nx);
  if (mgpu)
    mgpu_gradient_y<<<numBlocks_y, threadsPerBlock_y>>>(l, g2);
  else
    gradient_y<<<numBlocks_y, threadsPerBlock_y>>>(l, g2,nl,1./hx[1]);
  cudaCheckKernelError();
    
  // X-Gradient
  dim3 threadsPerBlock_x, numBlocks_x;
  getThreadBlockDimensionsX(threadsPerBlock_x, numBlocks_x, nx);
  if (mgpu)
    mgpu_gradient_x<<<numBlocks_x, threadsPerBlock_x>>>(l, g1);
  else
    gradient_x<<<numBlocks_x, threadsPerBlock_x>>>(l, g1,nl,1./hx[0]);
  cudaCheckKernelError();
  cudaDeviceSynchronize();
#endif
  
  PetscFunctionReturn(ierr);
}


PetscErrorCode computeDivergenceZ(ScalarType* l, ScalarType* gz, IntType* nx, ScalarType* hx, bool mgpu) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  const dim3 nl(nx[0], nx[1], nx[2]);
  
  dim3 threadsPerBlock_z, numBlocks_z;
  getThreadBlockDimensionsZ(threadsPerBlock_z, numBlocks_z, nx);
  if (mgpu)
    mgpu_gradient_z<<<numBlocks_z, threadsPerBlock_z>>>(l,gz);
  else
    gradient_z<<<numBlocks_z, threadsPerBlock_z>>>(l,gz,nl,1./hx[2]);
  cudaCheckKernelError();
  cudaDeviceSynchronize();
      
  PetscFunctionReturn(ierr);
}

PetscErrorCode computeDivergenceY(ScalarType* l, ScalarType* gy, IntType* nx, ScalarType* hx, bool mgpu) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  const dim3 nl (nx[0], nx[1], nx[2]);
  
  dim3 threadsPerBlock_y, numBlocks_y;
  getThreadBlockDimensionsY(threadsPerBlock_y, numBlocks_y, nx);
  if (mgpu)
    mgpu_gradient_y<<<numBlocks_y, threadsPerBlock_y>>>(l, gy);
  else
    gradient_y<<<numBlocks_y, threadsPerBlock_y>>>(l, gy, nl, 1./hx[1]);
  cudaCheckKernelError();
  cudaDeviceSynchronize();
  
  PetscFunctionReturn(ierr);
}


PetscErrorCode computeDivergenceX(ScalarType* l, ScalarType* gx, IntType* nx, ScalarType* hx, bool mgpu) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  const dim3 nl(nx[0], nx[1], nx[2]);
  
  dim3 threadsPerBlock_x, numBlocks_x;
  getThreadBlockDimensionsX(threadsPerBlock_x, numBlocks_x, nx);
  if (mgpu)
    mgpu_gradient_x<<<numBlocks_x, threadsPerBlock_x>>>(l, gx);
  else
    gradient_x<<<numBlocks_x, threadsPerBlock_x>>>(l, gx, nl, 1./hx[0]);
  cudaCheckKernelError();
  cudaDeviceSynchronize();

  PetscFunctionReturn(ierr);
}

PetscErrorCode computeGradient(ScalarType* gx, ScalarType* gy, ScalarType* gz, const ScalarType* m, cudaTextureObject_t mtex, IntType* nx, ScalarType* hx, bool mgpu) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // set all values to zero first
    size_t count = sizeof(ScalarType)*nx[0]*nx[1]*nx[2];
    if (mgpu) {
      ierr = cudaMemset((void*)gz, 0, count); CHKERRCUDA(ierr);
      ierr = cudaMemset((void*)gy, 0, count); CHKERRCUDA(ierr);
      ierr = cudaMemset((void*)gx, 0, count); CHKERRCUDA(ierr);
    }

#if defined(USE_TEXTURE_GRADIENT)
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
    dim3 nl (nx[0], nx[1], nx[2]);
    
    // Shared Memory implementation of Gradient Kernel
    // Z-Gradient
    dim3 threadsPerBlock_z, numBlocks_z;
    getThreadBlockDimensionsZ(threadsPerBlock_z, numBlocks_z, nx);
    if (mgpu)
      mgpu_gradient_z<<<numBlocks_z, threadsPerBlock_z>>>(gz,m);
    else
      gradient_z<replace_op><<<numBlocks_z, threadsPerBlock_z>>>(gz,m, nl, 1./hx[2]);
    cudaCheckKernelError();
    
    // Y-Gradient 
    dim3 threadsPerBlock_y, numBlocks_y;
    getThreadBlockDimensionsY(threadsPerBlock_y, numBlocks_y, nx);
    if (mgpu)
      mgpu_gradient_y<<<numBlocks_y, threadsPerBlock_y>>>(gy, m);
    else
      gradient_y<replace_op><<<numBlocks_y, threadsPerBlock_y>>>(gy, m, nl, 1./hx[1]);
    cudaCheckKernelError();
    
    // X-Gradient
    dim3 threadsPerBlock_x, numBlocks_x;
    getThreadBlockDimensionsX(threadsPerBlock_x, numBlocks_x, nx);
    if (mgpu)
      mgpu_gradient_x<<<numBlocks_x, threadsPerBlock_x>>>(gx, m);
    else
      gradient_x<replace_op><<<numBlocks_x, threadsPerBlock_x>>>(gx, m, nl, 1./hx[0]);
    cudaCheckKernelError();
    cudaDeviceSynchronize();
#endif
    PetscFunctionReturn(ierr);

}


PetscErrorCode computeLaplacian(ScalarType* ddm, const ScalarType* m, cudaTextureObject_t mtex, IntType* nx, ScalarType beta, bool mgpu) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    // Z-component
    dim3 threadsPerBlock_z, numBlocks_z;
    getThreadBlockDimensionsZ(threadsPerBlock_z, numBlocks_z, nx);
    if (mgpu)
      mgpu_d_zz<<<numBlocks_z, threadsPerBlock_z>>>(ddm, m, beta);
    else
      d_zz<<<numBlocks_z, threadsPerBlock_z>>>(ddm, m, beta);
    cudaCheckKernelError();
    
    // Y-component
    dim3 threadsPerBlock_y, numBlocks_y;
    getThreadBlockDimensionsY(threadsPerBlock_y, numBlocks_y, nx);
    if (mgpu)
      mgpu_d_yy<<<numBlocks_y, threadsPerBlock_y>>>(ddm, m, beta);
    else
      d_yy<<<numBlocks_y, threadsPerBlock_y>>>(ddm, m, beta);
    cudaCheckKernelError();
    
    // X-component
    dim3 threadsPerBlock_x, numBlocks_x;
    getThreadBlockDimensionsX(threadsPerBlock_x, numBlocks_x, nx);
    if (mgpu)
      mgpu_d_xx<<<numBlocks_x, threadsPerBlock_x>>>(ddm, m, beta);
    else
      d_xx<<<numBlocks_x, threadsPerBlock_x>>>(ddm, m, beta);
    cudaCheckKernelError();
    cudaDeviceSynchronize();
    
    PetscFunctionReturn(ierr);

}

}


#endif
