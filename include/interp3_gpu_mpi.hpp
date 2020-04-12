
#ifndef _INTERP3_GPU_MPI_HPP_
#define _INTERP3_GPU_MPI_HPP_


#include <mpi.h>
#include <set>
#include "interp3_common.hpp"
#include "interp3_gpu_new.hpp"
#include <vector>
#include <iostream>
#include "CLAIREUtils.hpp"

#if defined(PETSC_USE_REAL_SINGLE)
  typedef float Real;
  #define MPI_T MPI_FLOAT
  #define TC Complexf
  #define PL fftwf_plan
#else
  typedef double Real;
  #define MPI_T MPI_DOUBLE
  #define TC Complex
  #define PL fftw_plan
#endif

//#define INTERP_PINNED // if defined will use pinned memory for GPU
struct Interp3_Plan_GPU{
  public:
  Interp3_Plan_GPU(size_t g_alloc_max, bool cuda_aware);
  void allocate (int N_pts, int* data_dofs, int nplans);
  void scatter(int* N_reg, int * isize, int* istart, const int N_pts, const int g_size, Real* query_points_in_x, Real* query_points_in_y, Real* query_points_in_z, int* c_dims, MPI_Comm c_comm, double * timings);
  void test_kernel(Real* f, int nq);
  void interpolate( Real* ghost_reg_grid_vals, 
                    int* isize_g,
                    const IntType nlghost,
                    const IntType N_pts, 
                    Real* query_values,
                    MPI_Comm c_comm,
                    float *tmp1, 
                    float* tmp2,
                    cudaTextureObject_t yi_tex, 
                    int iporder, 
                    ScalarType* interp_time,
                    int version);

  size_t g_alloc_max; // size in bytes of the ghost input
  int N_reg_g[3];
  int isize_g[3];
  size_t total_query_points;
	int data_dof_max;
  int nplans_;
  int* data_dofs;
  size_t all_query_points_allocation;
  MPI_Datatype* stypes, *rtypes;
  bool cuda_aware;

  // CPU memory pointers
  Real* all_query_points;
  Real* all_f_cubic;
  Real* f_cubic_unordered;
  Real* query_points;
    
  int* f_index_procs_others_offset; // offset in the all_query_points array
  int* f_index_procs_self_offset  ; // offset in the query_outside array
  int* f_index_procs_self_sizes   ; // sizes of the number of interpolations that need to be sent to procs
  int* f_index_procs_others_sizes ; // sizes of the number of interpolations that need to be received from procs

  MPI_Request* s_request;
	MPI_Request* request;

  thrust::device_ptr<Real> query_outside;
  thrust::device_ptr<int> f_index;
  
  // vector to store number of query points in each proc
  int* num_query_per_proc;
  int* query_outside_offset;
  int* f_index_offset;
  int neighbour_query_width;
  size_t max_query_points_capacity;


  int* which_proc;
  Real *query_points_x, *query_points_y, *query_points_z;
  bool allocate_baked;
  bool scatter_baked;
  
  std::vector<int> procs_i_send_to_; // procs who i have to send my q
  std::vector<int> procs_i_recv_from_; // procs whose q I have to recv
  int procs_i_send_to_size_, procs_i_recv_from_size_;

  // GPU memory pointers
  Real* all_f_cubic_d;
  Real* all_query_points_d;
  Real * f_cubic_unordered_d;
  //Real* ghost_reg_grid_vals_d;
  Real *xq1, *xq2, *xq3; 

  ~Interp3_Plan_GPU();

};

void get_count(int* arr, int size, int val, int* count);
void print_norm(ScalarType* arr, int N);
void print_max(ScalarType* arr, int N);
void test_count();

#endif
