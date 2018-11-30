/*************************************************************************
 *  Copyright (c) 2018, Malte Brunn
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

#ifndef __KERNELUTILS_HPP__
#define __KERNELUTILS_HPP__

namespace KernelUtils {

#if defined(__CUDACC__) && defined(REG_HAS_CUDA) // for CUDA compiler
  #define __hostdevice__ __host__ __device__
#else // for regular compiler
  #define __hostdevice__
  struct int3 { int x, y, z; };
#endif

struct real3 { ScalarType x, y, z; };

/********************************************************************
 * @brief compute kernels are static methods
 * the first two arguments for spectral kernels are of type 
 *  int and real3, for local index and wave number, respectively
 * the first argument for general kernels is of type int
 *  for the local index
 * the first argument for reduction kernels is of type int
 *  for the local index. The kernel must return a ScalarType
 * 
 * A typical kernel looks like
 * struct KernelName {
 *    KernelOperator (int i, TypeName1 argument1, ...) {
 *    }
 * };
 * A reduction kernel looks like
 * struct KernelName {
 *    ReductionFunctional (int i, TypeName1 argument1, ...) {
 *      return localvalue;
 *    }
 * };
 * Regular, spectral, and reduction kernels can be combined in a
 * single struct/class.
*******************************************************************/
#define KernelOperator static __hostdevice__ inline void call
#define ReductionFunctional static __hostdevice__ inline ScalarType func

/********************************************************************
 * @brief computes wave number
 * TODO Note: The Nyquist is set to zero here! This is wrong for even (2nd, 4th, ...)
 * derivatives. Filtering the data with Nyquist = 0 changes the function space,
 * but makes this valid.
 *******************************************************************/
__hostdevice__ inline void ComputeWaveNumber(real3 &w, real3 n) {
  if (w.x > n.x/2) w.x -= n.x;
  else if (w.x == n.x/2) w.x = 0;
  if (w.y > n.y/2) w.y -= n.y;
  else if (w.y == n.y/2) w.y = 0;
  if (w.z > n.z/2) w.z -= n.z;
  else if (w.z == n.z/2) w.z = 0;
}

/********************************************************************
 * @brief computes linear array index in Z-Y-X layout
 * nl contains precomputed local sizes: 
 *  nl.x = nl0      local size in x
 *  nl.y = nl1*nl2  size of local pencil slice in Z-Y-X order
 *  nl.z = nl2      size of local pencil width in Z-Y-X order
 *******************************************************************/
__hostdevice__ inline int GetLinearIndex(int i1, int i2, int i3, int3 nl) {
  return i1*nl.y + i2*nl.z + i3;
}

/********************************************************************
 * @brief pre-processor evaluated power
 *******************************************************************/
template<int N> struct _pow {
  template<typename T> static __hostdevice__ inline T call(T x) {
    return x*_pow<N-1>::call(x);
  }
};
template<> struct _pow<0> {
  template<typename T> static __hostdevice__ inline T call(T x) {
    return 1;
  }
};
template<int N, typename T> __hostdevice__ inline T pow(T x) {
  return _pow<N>::call(x);
}

/********************************************************************
 * @brief computes absolute laplacian operator
 *******************************************************************/
__hostdevice__ inline ScalarType LaplaceNumber(real3 w) {
  return -(w.x*w.x + w.y*w.y + w.z*w.z);
}

#if defined(__CUDACC__) && defined(REG_HAS_CUDA) // compiled by CUDA compiler
/********************************************************************
 * @brief tree reduction for shared memory
 *******************************************************************/
template<int N> inline __device__ void LocalReductionSum(ScalarType *shared) {
  if (threadIdx.x < N) {
    shared[threadIdx.x] += shared[threadIdx.x + N];
  }
  __syncthreads();
  LocalReductionSum<N/2>(shared);
}
template<> inline __device__ void LocalReductionSum<1>(ScalarType *shared) {
  if (threadIdx.x == 0) {
    shared[0] += shared[1];
  }
  __syncthreads();
}

/********************************************************************
 * @brief tree reduction using N threads all in one block
 *******************************************************************/
template<int N> __global__ void ReductionSum(ScalarType *res, int n) {
  __shared__ ScalarType value[N];
  
  value[threadIdx.x] = 0.0;
  for (int i=threadIdx.x; i<n; i+=N) {
    value[threadIdx.x] += res[i];
  }
  
  __syncthreads();
  
  LocalReductionSum<N/2>(value);
  
  if (threadIdx.x == 0) {
    res[0] = value[0];
  }
}

/********************************************************************
 * @brief GPU kernel function wrapper for reduction kernels
 * this wrapper performs a reduction over all N threads in the block
 *******************************************************************/
template<int N, typename KernelFn, typename ... Args>
__global__ void ReductionKernelGPU(ScalarType *res, int nl, Args ... args) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  __shared__ ScalarType value[N];
  ScalarType tmp;
  
  value[threadIdx.x] = 0.0;
  if (i < nl) { // execute kernel
    value[threadIdx.x] = KernelFn::func(i, args...);
  }
  
  // start the tree reduction
  __syncthreads();
  LocalReductionSum<N/2>(value);
  if (threadIdx.x == 0) { // final value write back
    res[blockIdx.x] = value[0];
  }
}

/********************************************************************
 * @brief GPU kernel function wrapper for spectral operators
 * wave: must be initialized with the offset of local subdomain
 * nx: absolute dimensions of the domain
 * nl: contains precomputed local domain sizes: 
 *  nl.x = nl0      local size in x
 *  nl.y = nl1*nl2  size of local pencil slice in Z-Y-X order
 *  nl.z = nl2      size of local pencil width in Z-Y-X order
 *******************************************************************/
template<typename KernelFn, typename ... Args>
__global__ void SpectralKernelGPU(real3 wave, real3 nx, int3 nl, Args ... args) {
  int i1 = threadIdx.x + blockIdx.x*blockDim.x;
  int i2 = blockIdx.y;
  int i3 = blockIdx.z;
  
  if (i1 < nl.x) {
    wave.x += i1;
    wave.y += i2;
    wave.z += i3;

    ComputeWaveNumber(wave, nx);
    int i = GetLinearIndex(i1, i2, i3, nl);

    KernelFn::call(i, wave, args...);
  }
}

/********************************************************************
 * @brief GPU kernel function wrapper
 *******************************************************************/
template<typename KernelFn, typename ... Args>
__global__ void KernelGPU(int nl, Args ... args) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (i < nl) {
    KernelFn::call(i, args...);
  }
}

/********************************************************************
 * @brief Starts a GPU kernel for spectral operators
 * Note: Local index must not exceed $2^{31} - 1$
 *******************************************************************/
template<typename KernelFn, typename ... Args>
PetscErrorCode SpectralKernelCallGPU(IntType nstart[3], IntType nx[3], IntType nl[3], 
    Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  dim3 block(256,1,1); // 256 threads per block
  dim3 grid((nl[0] + 255)/256,nl[1],nl[2]); // $\lceil nl_0 / 256 \rceil, nl_1, nl_2 $
  real3 wave, nx3;
  int3 nl3;
  wave.x = nstart[0]; wave.y = nstart[1]; wave.z = nstart[2];
  nx3.x = nx[0]; nx3.y = nx[1]; nx3.z = nx[2];
  nl3.x = nl[0]; nl3.y = nl[1]*nl[2]; nl3.z = nl[2];
  
  if (nl[0]*nl[1]*nl[2] > 0) {
    SpectralKernelGPU<KernelFn><<<grid, block>>>(wave, nx3, nl3, args...);
    ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
    ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);
  }
  
  PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Starts a GPU kernel
 * Note: Local index must not exceed $2^{31} - 1$
 *******************************************************************/
template<typename KernelFn, typename ... Args>
PetscErrorCode KernelCallGPU(IntType nl, Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  dim3 block(256,1,1); // 256 threads per block
  dim3 grid((nl + 255)/256,1,1); // $\lceil nl_0 / 256 \rceil, nl_1, nl_2 $
  
  if (nl > 0) {
    KernelGPU<KernelFn><<<grid, block>>>(nl, args...);
    ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
    ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);
  }
  
  PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Starts a GPU kernel with value reduction
 * Note: Local index must not exceed $2^{31} - 1$
 * TODO: Work array is allocated each call
 *******************************************************************/
template<typename KernelFn, typename ... Args>
PetscErrorCode ReductionKernelCallGPU(ScalarType &value, IntType nl, Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  dim3 block(256, 1, 1); // 256 threads per block
  dim3 grid((nl + 255)/256, 1, 1);  // $\lceil nl_0 / 256 \rceil, nl_1, nl_2 $
  ScalarType *res = nullptr;

  ierr = reg::AllocateMemoryOnce(res, grid.x*sizeof(ScalarType)); CHKERRQ(ierr);
  
  value = 0.;
  
  if (nl > 0) {
    // execute the kernel and reduce over threads
    ReductionKernelGPU<256, KernelFn><<<grid, block>>>(res, nl, args...);
    ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
    ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);
    
    // reduce over work array
    ReductionSum<1024><<<1, 1024>>>(res, grid.x);
    ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
    ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);

    // copy result to cpu
    ierr = cudaMemcpy(reinterpret_cast<void*>(&value), reinterpret_cast<void*>(res), 
                      sizeof(ScalarType), cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
  }

  ierr = reg::FreeMemory(res); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

#endif

/********************************************************************
 * @brief Starts a CPU kernel for spectral operators
 * Note: Local index must not exceed $2^{31} - 1$
 *******************************************************************/
template<typename KernelFn, typename ... Args>
PetscErrorCode SpectralKernelCall(IntType nstart[3], IntType nx[3], IntType nl[3], 
    Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  real3 nx3;
  nx3.x = nx[0]; nx3.y = nx[1]; nx3.z = nx[2];
  int3 nl3;
  nl3.x = nl[0]; nl3.y = nl[1]*nl[2]; nl3.z = nl[2];
  
#pragma omp parallel for collapse(3)
  for (int i1 = 0; i1 < nl[0]; ++i1) {
      for (int i2 = 0; i2 < nl[1]; ++i2) {
          for (int i3 = 0; i3 < nl[2]; ++i3) {
              real3 w;
              w.x = static_cast<IntType>(i1) + nstart[0];
              w.y = static_cast<IntType>(i2) + nstart[1];
              w.z = static_cast<IntType>(i3) + nstart[2];
              
              ComputeWaveNumber(w, nx3);
              int i = GetLinearIndex(i1, i2, i3, nl3);
              
              KernelFn::call(i, w, args...);
          }
      }
  }
  
  PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Starts a CPU kernel
 * Note: Local index must not exceed $2^{31} - 1$
 *******************************************************************/
template<typename KernelFn, typename ... Args>
PetscErrorCode KernelCall(IntType nl, Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel for
  for (int i = 0; i < nl; ++i) {
      KernelFn::call(i, args...);
  }
  
  PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Starts a CPU kernel with value reduction
 * Note: Local index must not exceed $2^{31} - 1$
 *******************************************************************/
template<typename KernelFn, typename ... Args>
PetscErrorCode ReductionKernelCall(ScalarType &value, IntType nl, Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  value = 0.;
  
#pragma omp parallel for reduction(+:value)
  for (int i = 0; i < nl; ++i) {
    value += KernelFn::func(i, args...);
  }

  PetscFunctionReturn(ierr);
}

} // namespace KernelUtils

#endif // __KERNELUTILS_HPP__
