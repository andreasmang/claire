#ifndef _INTERP3_HPP_
#define _INTERP3_HPP_

#include "petsc.h"

//#include <accfft.h>
//#include <accfftf.h>
#include "TypeDef.hpp"
// #include <glog/logging.h>
// #define ParLOG if(procid==0) LOG(INFO)
#undef PCOUT
#define PCOUT if(procid==0) std::cerr
#define FAST_INTERP

#ifdef PETSC_USE_REAL_SINGLE
#ifndef POWER9
// enable ONLY for single precision
#define FAST_INTERPV
#endif
#endif

#define FAST_INTERP_BINNING
//#define HASWELL
//#define KNL

#if defined(KNL)
#define INTERP_USE_MORE_MEM_L1
#endif




#if defined(PETSC_USE_REAL_SINGLE)
  typedef float Real;
  #define MPI_T MPI_FLOAT
#else
  typedef double Real;
  #define MPI_T MPI_DOUBLE
#endif

#ifndef __INTEL_COMPILER
#define  _mm256_loadu2_m128(HIADDR, LOADDR) _mm256_set_m128(_mm_loadu_ps(HIADDR),_mm_loadu_ps(LOADDR))
//#undef FAST_INTERPV
#endif

#define COORD_DIM 3
#include <mpi.h>
#include <vector>
#include <set>
#include <interp3_common.hpp>

void rescale_xyz(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		int* isize, int* isize_g, const int N_pts, Real* Q_);
void rescale_xyzgrid(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		int* isize, int* isize_g, const int N_pts, pvfmm::Iterator<Real> Q_);
void interp3_p_col(Real* reg_grid_vals, int data_dof, int N_reg,
		const int N_pts, Real* query_points, Real* query_values);
void interp3_p(Real* reg_grid_vals, int data_dof, int N_reg, const int N_pts,
		Real* query_points, Real* query_values);

void interp3_p(Real* reg_grid_vals, int data_dof, int* N_reg, const int N_pts,
		Real* query_points, Real* query_values);

void vectorized_interp3_ghost_xyz_p(Real* __restrict reg_grid_vals, int data_dof, const int* __restrict N_reg,
		const int* __restrict N_reg_g, const int * __restrict isize_g, const int* __restrict istart, const int N_pts,
		const int g_size, Real* __restrict query_points, Real* __restrict query_values,
		bool query_values_already_scaled = false);

void optimized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * N_reg_g, int* isize_g, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values,
		bool query_values_already_scaled = false); // cubic interpolation

void interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * N_reg_g, int* isize_g, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values,
		bool query_values_already_scaled = false); // cubic interpolation

void interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * N_reg_g, int* isize_g, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values, int interp_order,
		bool query_values_already_scaled = false); // higher order interpolation

void interp3_ghost_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * N_reg_g, int* isize_g, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values);

void par_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * isize, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values, int* c_dims, MPI_Comm c_comm);

struct Interp3_Plan {
public:
	Interp3_Plan();
  pvfmm::Iterator<Real> query_points;
	void allocate(int N_pts, int* data_dofs = NULL,int nplans = 1);
	void slow_run(Real* ghost_reg_grid_vals, int data_dof, int* N_reg,
			int * isize, int* istart, const int N_pts, const int g_size,
			Real* query_points_in, Real* query_values, int* c_dims,
			MPI_Comm c_comm);
	void fast_scatter(int* N_reg, int * isize, int* istart,
			const int N_pts, const int g_size, Real* query_points_in,
			int* c_dims, MPI_Comm c_comm, double * timings);
	void scatter(int* N_reg, int * isize, int* istart,
			const int N_pts, const int g_size, Real* query_points_in,
			int* c_dims, MPI_Comm c_comm, double * timings);
	// void interpolate(Real* ghost_reg_grid_vals, int data_dof, int* N_reg,
	// 		int * isize, int* istart, const int N_pts, const int g_size,
	// 		Real* query_values, int* c_dims, MPI_Comm c_comm, double * timings);
/*
  void interpolate(Real* __restrict ghost_reg_grid_vals,
		int*__restrict N_reg, int *__restrict isize, int*__restrict istart, const int N_pts, const int g_size,
		Real*__restrict query_values, int*__restrict c_dims, MPI_Comm c_comm, double *__restrict timings, int version =0);
  */
void interpolate(Real* __restrict ghost_reg_grid_vals,
		int*__restrict N_reg, int *__restrict isize, int*__restrict istart, const int N_pts, const int g_size,
		Real*__restrict query_values, int*__restrict c_dims, MPI_Comm c_comm, double *__restrict timings, int version =0);

	void high_order_interpolate(Real* ghost_reg_grid_vals, int data_dof, int* N_reg,
			int * isize, int* istart, const int N_pts, const int g_size,
			Real* query_values, int* c_dims, MPI_Comm c_comm, double * timings, int interp_order);
	int N_reg_g[3];
	int isize_g[3];
	size_t total_query_points;
	int data_dof_max;
  int nplans_;
  pvfmm::Iterator<int> data_dofs_;
	size_t all_query_points_allocation;
  pvfmm::Iterator<MPI_Datatype> stypes, rtypes;

  pvfmm::Iterator<Real> all_query_points;
  pvfmm::Iterator<Real> all_f_cubic;
  pvfmm::Iterator<Real> f_cubic_unordered;

  pvfmm::Iterator<int> f_index_procs_others_offset; // offset in the all_query_points array
	pvfmm::Iterator<int> f_index_procs_self_offset; // offset in the query_outside array
	pvfmm::Iterator<int> f_index_procs_self_sizes; // sizes of the number of interpolations that need to be sent to procs
	pvfmm::Iterator<int> f_index_procs_others_sizes; // sizes of the number of interpolations that need to be received from procs

  pvfmm::Iterator<MPI_Request> s_request;
	pvfmm::Iterator<MPI_Request> request;

	std::vector<int> *f_index;
	std::vector<Real> *query_outside;

	bool allocate_baked;
	bool scatter_baked;

  std::vector<int> procs_i_send_to_; // procs who i have to send my q
  std::vector<int> procs_i_recv_from_; // procs whose q I have to recv
  int procs_i_send_to_size_, procs_i_recv_from_size_;

	~Interp3_Plan();

};

void par_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * isize, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values, int* c_dims, MPI_Comm c_comm,
		Interp3_Plan* interp_plan);

void gpu_interp3_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		const int N_pts, Real* query_points, Real* query_values);

void gpu_interp3_ghost_xyz_p(Real* reg_grid_vals_d, int data_dof, int* N_reg,
		int* istart, int * isize, const int N_pts, const int g_size,
		Real* query_points_d, Real* query_values_d);

void gpu_par_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * isize, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values, int* c_dims, MPI_Comm c_comm);

//class Interp_Mem_Mgr{
//	public:
//		Mem_Mgr();
//
//    bool PINNED;
//
//    Real * buffer;
//    Real * buffer_2;
//    Real * buffer_d;
//
//		~Mem_Mgr();
//};


// GHOST FUNCTIONS

namespace reg {
class RegOpt;
}

size_t ghost_local_size(reg::RegOpt* m_Opt, int g_size,
		int * isize_g, int* istart_g);
//size_t accfft_ghost_local_size_dft_r2c(accfft_plan* plan, int g_size,
//		int * isize_g, int* istart_g);
void get_ghost(reg::RegOpt* m_Optn, int g_size, int* isize_g, Real* data,
		Real* ghost_data);
//void accfft_get_ghost(accfft_plan* plan, int g_size, int* isize_g, Real* data,
//		Real* ghost_data);

size_t ghost_xyz_local_size(reg::RegOpt* m_Opt, int g_size,
		int * isize_g, int* istart_g);
//size_t accfft_ghost_xyz_local_size_dft_r2c(accfft_plan* plan, int g_size,
//		int * isize_g, int* istart_g);

void get_ghost_xyz(reg::RegOpt* m_Opt, int g_size, int* isize_g,
		Real* data, Real* ghost_data);

void share_ghost_layer(reg::RegOpt* m_Opt, int g_size, int* isize_g,
		Real* data, Real* ghost_data, pvfmm::Iterator<Real> padded_data, pvfmm::Iterator<Real> ghost_data_xy);

//void accfft_get_ghost_xyz(accfft_plan* plan, int g_size, int* isize_g,
//		Real* data, Real* ghost_data);

#endif
