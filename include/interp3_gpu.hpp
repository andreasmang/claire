
#ifndef _INTERP3_GPU_HPP_
#define _INTERP3_GPU_HPP_
typedef float Real;
//#define MPI_T MPI_DOUBLE

#define COORD_DIM 3

#include <cuda.h>
#include <cuda_runtime_api.h>

void gpu_interp3_p( Real* reg_grid_vals, int data_dof,
              int* N_reg, const int N_pts, Real* query_points,
              Real* query_values);

void gpu_interp3_ghost_xyz_p( Real* reg_grid_vals_d, int data_dof,
    int* N_reg, int* istart, int * isize, const int N_pts, const int g_size, Real* query_points_d,
    Real* query_values_d,bool query_values_already_scaled=false);

__global__ void gpu_interp3_ghost_xyz_p_kernel( Real* reg_grid_vals, int data_dof,
    int* N_reg, int* isize,int* istart, const int N_pts, const int g_size, Real* query_points,
    Real* query_values,bool query_values_already_scaled=false);
#endif
