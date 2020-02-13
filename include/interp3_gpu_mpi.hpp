
#ifndef _INTERP3_GPU_MPI_HPP_
#define _INTERP3_GPU_MPI_HPP_


//#include <interp3_common.hpp>
//#include <interp3_gpu.hpp>
#include <interp3_gpu_new.hpp>
#include <mpi.h>
#include <vector>

//#if defined(PETSC_USE_REAL_SINGLE)
  typedef float Real;
  #define MPI_T MPI_FLOAT
  #define TC Complexf
  #define PL fftwf_plan
//#else
//  typedef double Real;
//  #define MPI_T MPI_DOUBLE
//  #define TC Complex
//  #define PL fftw_plan
//#endif

#define INTERP_PINNED // if defined will use pinned memory for GPU
struct Interp3_Plan_GPU{
  public:
  Interp3_Plan_GPU(size_t g_alloc_max);
  Real * query_points;
  void allocate (int N_pts, int data_dof);
  //void slow_run( Real* ghost_reg_grid_vals, int data_dof,
  //            int* N_reg, int * isize, int* istart, const int N_pts, const int g_size, Real* query_points_in,
  //            Real* query_values,int* c_dims, MPI_Comm c_comm);
  void scatter(int data_dof,
              int* N_reg, int * isize, int* istart, const int N_pts, const int g_size, Real* query_points_in,
              int* c_dims, MPI_Comm c_comm, double * timings);
  //void interpolate( Real* ghost_reg_grid_vals, int data_dof,
  //            int* N_reg, int * isize, int* istart, const int N_pts, const int g_size,
  //            Real* query_values,int* c_dims, MPI_Comm c_comm,double * timings);
  
  // current interpolation function
  void interpolate( Real* ghost_reg_grid_vals, 
                    int data_dof,
                    int* N_reg, 
                    int * isize, 
                    int* istart, 
                    int* isize_g,
                    const int nlghost,
                    const int N_pts, 
                    const int g_size,
                    Real* query_values,
                    int* c_dims, 
                    MPI_Comm c_comm,
                    double * timings, 
                    float *tmp1, 
                    float* tmp2,
                    cudaTextureObject_t yi_tex, 
                    int iporder, 
                    ScalarType* interp_time);

  size_t g_alloc_max; // size in bytes of the ghost input
  int N_reg_g[3];
  int isize_g[3];
  int total_query_points;
  int data_dof;
  size_t all_query_points_allocation;
  MPI_Datatype *stype,*rtype;
    
  // CPU memory pointers
  Real * all_query_points;
  Real* all_f_cubic;
  Real * f_cubic_unordered;
  int* f_index_procs_others_offset; // offset in the all_query_points array
  int* f_index_procs_self_offset  ; // offset in the query_outside array
  int* f_index_procs_self_sizes   ; // sizes of the number of interpolations that need to be sent to procs
  int* f_index_procs_others_sizes ; // sizes of the number of interpolations that need to be received from procs

  MPI_Request * s_request;
  MPI_Request * request;

  std::vector <int> *f_index;
  std::vector <Real> *query_outside;

  bool allocate_baked;
  bool scatter_baked;

  // GPU memory pointers
  Real* all_f_cubic_d;
  Real* all_query_points_d;
  Real* ghost_reg_grid_vals_d;
  Real *xq1, *xq2, *xq3; 

  ~Interp3_Plan_GPU();

};

void gpu_par_interp3_ghost_xyz_p( Real* reg_grid_vals, int data_dof,
              int* N_reg, int * isize, int* istart, const int N_pts,int g_size, Real* query_points,
              Real* query_values,int* c_dims,MPI_Comm c_comm);

#endif
