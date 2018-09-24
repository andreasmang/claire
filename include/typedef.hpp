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

#if defined(REG_FFT_CUDA)
#include "accfft_gpuf.h"
#include "accfft_operators_gpu.h"
#else
#include "accfft.h"
#include "accfftf.h"
#include "accfft_operators.h"
#endif

#if defined(REG_DBG_CUDA)
#include "nvToolsExt.h"
#define DebugGPUStartEvent(str) nvtxRangePushA(str)
#define DebugGPUStopEvent() nvtxRangePop()
#else
#define DebugGPUStartEvent(str)
#define DebugGPUStopEvent()
#endif

using IntType = PetscInt;
using ScalarType = PetscReal;

#if defined(PETSC_USE_REAL_SINGLE)
using ComplexType = Complexf;
using FFTWPlanType = fftwf_plan;
#else
using ComplexType = Complex;
using FFTWPlanType = fftw_plan;
#endif

#if defined(REG_FFT_CUDA) // GPU FFT
#if defined(PETSC_USE_REAL_SINGLE)
using FFTPlanType = accfft_plan_gpuf;
//#define accfft_plan_dft_3d_r2c accfft_plan_dft_3d_r2c_gpuf
const auto accfft_plan_dft_3d_r2c = accfft_plan_dft_3d_r2c_gpuf;
const auto accfft_cleanup = accfft_cleanup_gpuf;
//inline void accfft_cleanup() { accfft_cleanup_gpuf(); }
#else
using FFTPlanType = accfft_plan_gpu;
//#define accfft_plan_dft_3d_r2c accfft_plan_dft_3d_r2c_gpu
const auto accfft_plan_dft_3d_r2c = accfft_plan_dft_3d_r2c_gpu;
const auto accfft_cleanup = accfft_cleanup_gpu;
//inline void accfft_cleanup() { accfft_cleanup_gpu(); }
#endif
// dummy function replacing fftw init (not needed by GPU fft)
inline int accfft_init(int nthreads) { return 0; }
// AccFFT wrapper for GPU kernels
#define accfft_execute_r2c_t accfft_execute_r2c_gpu_t
//const auto accfft_execute_r2c_t = accfft_execute_r2c_gpu_t<ScalarType>;
#define accfft_execute_c2r_t accfft_execute_c2r_gpu_t
//const auto accfft_execute_c2r_t = accfft_execute_c2r_gpu_t<ScalarType>;
#define accfft_grad_t accfft_grad_gpu_t
#define accfft_laplace_t accfft_laplace_gpu_t
#define accfft_divergence_t accfft_divergence_gpu_t
#define accfft_biharmonic_t accfft_biharmonic_gpu_t
#else // CPU FFT:
using FFTPlanType = accfft_plan_t<ScalarType, ComplexType, FFTWPlanType>;
#endif

#endif
