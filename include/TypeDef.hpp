#ifndef __TYPEDEF_HPP__
#define __TYPEDEF_HPP__

// local includes
#include "petsc.h"
#include "petscsys.h"

#if defined(REG_HAS_CUDA) || defined(REG_FFT_CUDA)
  #include "petsccuda.h"
  #include "cuda.h"
  #include <petsc/private/vecimpl.h>
#endif

#ifdef REG_FFT_CUDA
  #include "accfft_gpu.h"
  #include "accfft_gpuf.h"
  #include "accfft_operators_gpu.h"
#else
  #include "accfft.h"
  #include "accfftf.h"
  #include "accfft_operators.h"
#endif

#if defined(REG_DBG_CUDA) && defined(REG_HAS_CUDA)
  #include "nvToolsExt.h"
  #define DebugGPUStartEvent(str) nvtxRangePushA(str)
  #define DebugGPUStopEvent() nvtxRangePop()
  #define DebugGPUNotImplemented() WrngMsg("Not implemented for GPU")
#else
  #define DebugGPUStartEvent(str)
  #define DebugGPUStopEvent()
  #define DebugGPUNotImplemented() 0
#endif

#define DebugNotImplemented() ThrowError("Not implemented")

#define WrngMsg(str) WrngMsgCall(str, __LINE__, __FILE__);
#define ThrowError(str) ThrowErrorMsg(str, __LINE__, __FILE__)
#define DbgMsg(str) DbgMsgCall(str, __LINE__, __FILE__)


using IntType = PetscInt;
using ScalarType = PetscReal;

#if defined(PETSC_USE_REAL_SINGLE)
  using ComplexType = Complexf;
  using FFTWPlanType = fftwf_plan;
#else
  using ComplexType = Complex;
  using FFTWPlanType = fftw_plan;
#endif

#ifdef REG_FFT_CUDA // GPU FFT
  #ifdef PETSC_USE_REAL_SINGLE
    using FFTPlanType = accfft_plan_gpuf;
    const auto accfft_plan_dft_3d_r2c = accfft_plan_dft_3d_r2c_gpuf;
    const auto accfft_cleanup = accfft_cleanup_gpuf;
  #else
    using FFTPlanType = accfft_plan_gpu;
    const auto accfft_plan_dft_3d_r2c = accfft_plan_dft_3d_r2c_gpu;
    const auto accfft_cleanup = accfft_cleanup_gpu;
  #endif
  // dummy function replacing fftw init (not needed by GPU fft)
  inline int accfft_init(int nthreads) { return 0; }
  // AccFFT wrapper for GPU kernels
  #define accfft_execute_r2c_t accfft_execute_r2c_gpu_t
  #define accfft_execute_c2r_t accfft_execute_c2r_gpu_t
  #define accfft_grad_t accfft_grad_gpu_t
  #define accfft_laplace_t accfft_laplace_gpu_t
  #define accfft_divergence_t accfft_divergence_gpu_t
  #define accfft_biharmonic_t accfft_biharmonic_gpu_t
#else // CPU FFT:
  using FFTPlanType = accfft_plan_t<ScalarType, ComplexType, FFTWPlanType>;
#endif

/*
#include <cuda_runtime.h>
#undef PetscFunctionBegin
#define PetscFunctionBegin { \
  size_t free_mem; \
  size_t total_mem; \
  cudaMemGetInfo(&free_mem, &total_mem); \
  printf("\t\x1b[34m%10lf MiB \x1b[33m%s:\x1b[32m%i\x1b[0m\n", static_cast<double>(total_mem - free_mem)/1048576.0,  __FILE__, __LINE__); \
}*/

#endif // __TYPEDEF_HPP__
