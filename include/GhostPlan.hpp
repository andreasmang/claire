#ifndef _GHOST_PLAN_HPP_
#define _GHOST_PLAN_HPP_

#include <stdlib.h>
#include <math.h> 
#include <mpi.h>
#include <vector>
#include <RegOpt.hpp>
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
    GhostPlan(RegOpt*, IntType);
    ~GhostPlan();
    
    RegOpt* m_Opt;

    ScalarType* padded_data;
    ScalarType* ghost_data_xy;
    // top bottom
    ScalarType *GB, *GT;
    // left right
    ScalarType *RS, *GL;
    ScalarType *LS, *GR;

    size_t g_alloc_max;
    int g_size;
    int isize_g[3], istart_g[3];
    int isize[3], istart[3];

    const int num_streams = 1;
    cudaStream_t* stream;
    
    void share_ghost_xy(const ScalarType* data, ScalarType* ghost_data);
    
    void share_ghost_x(const ScalarType* data, ScalarType* ghost_data);
    
    void share_left_right(const ScalarType* data, double* timers);

    void share_top_bottom(ScalarType* ghost_data, double* timers);

    size_t get_ghost_local_size_x(IntType * isize_g, IntType* istart_g);
    size_t get_ghost_local_size_xy(IntType * isize_g, IntType* istart_g);
};

}
#endif
