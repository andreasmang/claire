
#ifndef _INTERP3_GPU_MPI_HPP_
#define _INTERP3_GPU_MPI_HPP_


#include <mpi.h>
#include <set>
#include "interp3_common.hpp"
#include "interp3_gpu_new.hpp"
#include <vector>
#include <iostream>
#include "CLAIREUtils.hpp"
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/extrema.h>

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

//#define INTERP_PINNED // if defined will use pinned memory for GPU
struct Interp3_Plan_GPU{
  public:
  Interp3_Plan_GPU(size_t g_alloc_max, bool cuda_aware);
  void allocate (int N_pts, int* data_dofs, int nplans, int gsize, IntType* isize_g);
  void scatter(IntType* N_reg, IntType * isize, IntType* istart, const int N_pts, const int g_size, ScalarType* query_points_in_x, ScalarType* query_points_in_y, ScalarType* query_points_in_z, int* c_dims, MPI_Comm c_comm, double * timings, std::string flag);
  void test_kernel(ScalarType* f, int nq);
  void interpolate( ScalarType* ghost_reg_grid_vals, 
                    IntType* isize_g,
                    const IntType nlghost,
                    const IntType N_pts, 
                    ScalarType** query_values,
                    MPI_Comm c_comm,
                    float *tmp1, 
                    float* tmp2,
                    cudaTextureObject_t yi_tex, 
                    int iporder, 
                    ScalarType* interp_time,
                    int version,
                    std::string flag);
  
  // number of equations = state and adjoint 
  static const int numeq = 2;
  // size in bytes of the ghost input
  size_t g_alloc_max; 
  //int N_reg_g[3];
  //int isize_g[3];
  //size_t total_query_points;
	int data_dof_max;
	// for scalar and vector field
  int nplans_;
  int* data_dofs;
  //size_t all_query_points_allocation;
  //MPI_Datatype* stypes, *rtypes;
  bool cuda_aware;
  bool output_baked;
  size_t all_f_capacity;
  
  /*
  int* f_index_procs_others_offset; // offset in the all_query_points array
  int* f_index_procs_self_offset  ; // offset in the query_send array
  int* f_index_procs_self_sizes   ; // sizes of the number of interpolations that need to be sent to procs
  int* f_index_procs_others_sizes ; // sizes of the number of interpolations that need to be received from procs
  */

  MPI_Request* s_request;
	MPI_Request* request;
  
  std::vector<int> proc_neighbours; // for slab decomp

  thrust::device_ptr<ScalarType> query_send;
  size_t query_points_send_capacity;
  int neighbour_query_send_width;
  bool query_baked;
  int* query_send_offset;
  //int* num_query_per_proc;
  //thrust::device_ptr<int> f_index;
  //int* f_index_offset;
  
  //int neighbour_query_recv_width;
  //size_t query_points_recv_capacity;

  short* which_proc;
  //int* which_proc;
  //ScalarType *query_points_x, *query_points_y, *query_points_z;
  //bool allocate_baked;
  //bool scatter_baked;
  
  //std::vector<int> procs_i_send_to_; // procs who i have to send my q
  //std::vector<int> procs_i_recv_from_; // procs whose q I have to recv
  //int procs_i_send_to_size_, procs_i_recv_from_size_;

  // GPU memory pointers for state and adjoint fields
  ScalarType* all_f;
  //ScalarType* all_query_points;
  ScalarType* f_unordered;
  //ScalarType *xq1, *xq2, *xq3; 

  ~Interp3_Plan_GPU();

  struct Equation {
    bool allocate_baked;
    bool scatter_baked;
    size_t total_query_points;
    size_t all_query_points_allocation;
    MPI_Datatype* stypes, *rtypes;

    int* f_index_procs_others_offset; // offset in the all_query_points array
    int* f_index_procs_self_offset  ; // offset in the query_send array
    int* f_index_procs_self_sizes   ; // sizes of the number of interpolations that need to be sent to procs
    int* f_index_procs_others_sizes ; // sizes of the number of interpolations that need to be received from procs

    thrust::device_ptr<int> f_index;
    int* f_index_offset;
    int* num_query_per_proc;
    int neighbour_query_recv_width;
    size_t query_points_recv_capacity;
  
    //ScalarType* all_f;
    ScalarType* all_query_points;
    ScalarType *xq1, *xq2, *xq3; 
  };

  Equation Eq[numeq];

};

void get_count(short* arr, int size, int val, int* count);
void print_norm(ScalarType* arr, int N);
void print_max(ScalarType* arr, int N);
void test_count();
//void print_vector(thrust::device_ptr<ScalarType> arr, int N, int stride);

#endif
