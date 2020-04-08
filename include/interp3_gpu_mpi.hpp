
#ifndef _INTERP3_GPU_MPI_HPP_
#define _INTERP3_GPU_MPI_HPP_


#include <mpi.h>
#include <set>
#include "interp3_common.hpp"
#include "interp3_gpu_new.hpp"
#include "CLAIREUtils.hpp"

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

//#define INTERP_PINNED // if defined will use pinned memory for GPU
struct Interp3_Plan_GPU{
  public:
  Interp3_Plan_GPU(size_t g_alloc_max);
  void allocate (int N_pts, int* data_dof, int nplans);
  void scatter(int* N_reg, int * isize, int* istart, const int N_pts, const int g_size, Real* query_points_in, int* c_dims, MPI_Comm c_comm, double * timings);
  
  void interpolate( Real* ghost_reg_grid_vals, 
                    int* N_reg, 
                    int * isize, 
                    int* istart, 
                    int* isize_g,
                    const IntType nlghost,
                    const IntType N_pts, 
                    const int g_size,
                    Real* query_values,
                    int* c_dims, 
                    MPI_Comm c_comm,
                    double * timings, 
                    float *tmp1, 
                    float* tmp2,
                    cudaTextureObject_t yi_tex, 
                    int iporder, 
                    ScalarType* interp_time,
                    int version);

  size_t g_alloc_max; // size in bytes of the ghost input
  int N_reg_g[3];
  int isize_g[3];
  int total_query_points;
	int data_dof_max;
  int nplans_;
  pvfmm::Iterator<int> data_dofs;
  size_t all_query_points_allocation;
  pvfmm::Iterator<MPI_Datatype> stypes, rtypes;
  
  pvfmm::Iterator<Real> all_query_points;
  pvfmm::Iterator<Real> all_f_cubic;
  pvfmm::Iterator<Real> f_cubic_unordered;
  pvfmm::Iterator<Real> query_points;

  pvfmm::Iterator<int> f_index_procs_others_offset; // offset in the all_query_points array
	pvfmm::Iterator<int> f_index_procs_self_offset; // offset in the query_outside array
	pvfmm::Iterator<int> f_index_procs_self_sizes; // sizes of the number of interpolations that need to be sent to procs
	pvfmm::Iterator<int> f_index_procs_others_sizes; // sizes of the number of interpolations that need to be received from procs

  pvfmm::Iterator<MPI_Request> s_request;
	pvfmm::Iterator<MPI_Request> request;

  std::vector <int> *f_index;
  std::vector <Real> *query_outside;

  bool allocate_baked;
  bool scatter_baked;
  
  std::vector<int> procs_i_send_to_; // procs who i have to send my q
  std::vector<int> procs_i_recv_from_; // procs whose q I have to recv
  int procs_i_send_to_size_, procs_i_recv_from_size_;

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
