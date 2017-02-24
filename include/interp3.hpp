#ifndef _INTERP3_HPP_
#define _INTERP3_HPP_
typedef double Real;
#define MPI_T MPI_DOUBLE

//typedef double Real;
//#define MPI_T MPI_DOUBLE
#define COORD_DIM 3
#include <mpi.h>
#include <vector>
#include <interp3_common.hpp>
void rescale_xyz(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		const int N_pts, Real* query_points);
void interp3_p_col(Real* reg_grid_vals, int data_dof, int N_reg,
		const int N_pts, Real* query_points, Real* query_values);
void interp3_p(Real* reg_grid_vals, int data_dof, int N_reg, const int N_pts,
		Real* query_points, Real* query_values);

void interp3_p(Real* reg_grid_vals, int data_dof, int* N_reg, const int N_pts,
		Real* query_points, Real* query_values);

void interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * N_reg_g, int* isize_g, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values,
		bool query_values_already_scaled = false);

void interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * N_reg_g, int* isize_g, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values, int interp_order,
		bool query_values_already_scaled = false);

void interp3_ghost_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * N_reg_g, int* isize_g, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values);

void par_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int * isize, int* istart, const int N_pts, int g_size,
		Real* query_points, Real* query_values, int* c_dims, MPI_Comm c_comm);

struct Interp3_Plan {
public:
	Interp3_Plan();
	Real * query_points;
	void allocate(int N_pts, int data_dof);
	void slow_run(Real* ghost_reg_grid_vals, int data_dof, int* N_reg,
			int * isize, int* istart, const int N_pts, const int g_size,
			Real* query_points_in, Real* query_values, int* c_dims,
			MPI_Comm c_comm);
	void scatter(int data_dof, int* N_reg, int * isize, int* istart,
			const int N_pts, const int g_size, Real* query_points_in,
			int* c_dims, MPI_Comm c_comm, double * timings);
	void interpolate(Real* ghost_reg_grid_vals, int data_dof, int* N_reg,
			int * isize, int* istart, const int N_pts, const int g_size,
			Real* query_values, int* c_dims, MPI_Comm c_comm, double * timings);
	void high_order_interpolate(Real* ghost_reg_grid_vals, int data_dof, int* N_reg,
			int * isize, int* istart, const int N_pts, const int g_size,
			Real* query_values, int* c_dims, MPI_Comm c_comm, double * timings, int interp_order);

	int N_reg_g[3];
	int isize_g[3];
	int total_query_points;
	int data_dof;
	size_t all_query_points_allocation;
	MPI_Datatype *stype, *rtype;

	Real * all_query_points;
	Real* all_f_cubic;
	Real * f_cubic_unordered;
	int* f_index_procs_others_offset; // offset in the all_query_points array
	int* f_index_procs_self_offset; // offset in the query_outside array
	int* f_index_procs_self_sizes; // sizes of the number of interpolations that need to be sent to procs
	int* f_index_procs_others_sizes; // sizes of the number of interpolations that need to be received from procs

	MPI_Request * s_request;
	MPI_Request * request;

	std::vector<int> *f_index;
	std::vector<Real> *query_outside;

	bool allocate_baked;
	bool scatter_baked;

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

#endif
