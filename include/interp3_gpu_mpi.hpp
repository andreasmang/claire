
#ifndef _INTERP3_GPU_MPI_HPP_
#define _INTERP3_GPU_MPI_HPP_

#include <interp3_gpu.hpp>
#include <interp3_common.hpp>
#include <mpi.h>
#include <vector>

#define INTERP_PINNED // if defined will use pinned memory for GPU
struct Interp3_Plan_GPU{
  public:
  Interp3_Plan_GPU(size_t g_alloc_max);
  Real * query_points;
  void allocate (int N_pts, int data_dof);
  void slow_run( Real* ghost_reg_grid_vals, int data_dof,
              int* N_reg, int * isize, int* istart, const int N_pts, const int g_size, Real* query_points_in,
              Real* query_values,int* c_dims, MPI_Comm c_comm);
  void scatter(int data_dof,
              int* N_reg, int * isize, int* istart, const int N_pts, const int g_size, Real* query_points_in,
              int* c_dims, MPI_Comm c_comm, double * timings);
  void interpolate( Real* ghost_reg_grid_vals, int data_dof,
              int* N_reg, int * isize, int* istart, const int N_pts, const int g_size,
              Real* query_values,int* c_dims, MPI_Comm c_comm,double * timings);

  size_t g_alloc_max; // size in bytes of the ghost input
  int N_reg_g[3];
  int isize_g[3];
  int total_query_points;
  int data_dof;
  size_t all_query_points_allocation;
  MPI_Datatype *stype,*rtype;

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


  Real* all_f_cubic_d;
  Real* all_query_points_d;
  Real* ghost_reg_grid_vals_d;

  ~Interp3_Plan_GPU();

};

void gpu_par_interp3_ghost_xyz_p( Real* reg_grid_vals, int data_dof,
              int* N_reg, int * isize, int* istart, const int N_pts,int g_size, Real* query_points,
              Real* query_values,int* c_dims,MPI_Comm c_comm);

#endif
