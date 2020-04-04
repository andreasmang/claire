#ifndef _GHOST_PLAN_HPP_
#define _GHOST_PLAN_HPP_

#include <stdlib.h>
#include <math.h> 
#include <mpi.h>
#include <vector>
#include <CLAIREUtils.hpp>
//#include <interp3.hpp>
#include <interp3_gpu_new.hpp>

#if defined(PETSC_USE_REAL_SINGLE)
  typedef float ScalarType;
  #define MPI_T MPI_FLOAT
  #define TC Complexf
  #define PL fftwf_plan
#else
  typedef double ScalarType;
  #define MPI_T MPI_DOUBLE
  #define TC Complex
  #define PL fftw_plan
#endif

namespace reg {

struct GhostPlan {
  
  public:
  GhostPlan(FFTPlanType*, int);
  ~GhostPlan();
  
  FFTPlanType* plan;
  int g_size;

  ScalarType* padded_data;
  ScalarType* ghost_data_xy;
  // top bottom
  ScalarType *GB, *GT;
  // left right
  ScalarType *RS, *GL;
  ScalarType *LS, *GR;

  size_t g_alloc_max;
  int isize_g[3], isize[3];
  int istart_g[3], istart[3];

  const int num_streams = 1;
  cudaStream_t* stream;
  
  void allocate(); 

  void share_ghost_xy(ScalarType* data, ScalarType* ghost_data, double* timers);
  
  void share_left_right(ScalarType* data, double* timers);

  void share_top_bottom(ScalarType* ghost_data, double* timers);

  size_t get_ghost_local_size(int * isize_g, int* istart_g);
};

}
#endif
