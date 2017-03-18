
#include <interp3.hpp>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <limits.h>
#define FAST_INTERP
#define FAST_INTERP_BINNING
Interp3_Plan::Interp3_Plan() {
	this->allocate_baked = false;
	this->scatter_baked = false;
}

void Interp3_Plan::allocate(int N_pts, int data_dof) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));

	f_index_procs_others_offset = (int*) malloc(nprocs * sizeof(int)); // offset in the all_query_points array
	f_index_procs_self_offset = (int*) malloc(nprocs * sizeof(int)); // offset in the query_outside array
	f_index_procs_self_sizes = (int*) malloc(nprocs * sizeof(int)); // sizes of the number of interpolations that need to be sent to procs
	f_index_procs_others_sizes = (int*) malloc(nprocs * sizeof(int)); // sizes of the number of interpolations that need to be received from procs

	s_request = new MPI_Request[nprocs];
	request = new MPI_Request[nprocs];

	f_index = new std::vector<int>[nprocs];
	query_outside = new std::vector<Real>[nprocs];

	f_cubic_unordered = (Real*) malloc(N_pts * sizeof(Real) * data_dof); // The reshuffled semi-final interpolated values are stored here

	stype = new MPI_Datatype[nprocs];
	rtype = new MPI_Datatype[nprocs];
	this->data_dof = data_dof;
	this->allocate_baked = true;
}

class Trip {
public:
	Trip() {
	}
	;
	Real x;
	Real y;
	Real z;
	int ind;
	int* N;
	Real* h;

};
static bool ValueCmp(Trip const & a, Trip const & b) {
	return a.z + a.y / a.h[1] * a.N[2] + a.x / a.h[0] * a.N[1] * a.N[2]
			< b.z + b.y / b.h[1] * b.N[2] + b.x / b.h[0] * b.N[1] * b.N[2];
}
static void sort_queries(std::vector<Real>* query_outside,
		std::vector<int>* f_index, int* N_reg, Real* h, MPI_Comm c_comm) {

	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);
	for (int proc = 0; proc < nprocs; ++proc) {
		int qsize = query_outside[proc].size() / COORD_DIM;
		Trip* trip = new Trip[qsize];

		for (int i = 0; i < qsize; ++i) {
			trip[i].x = query_outside[proc][i * COORD_DIM + 0];
			trip[i].y = query_outside[proc][i * COORD_DIM + 1];
			trip[i].z = query_outside[proc][i * COORD_DIM + 2];
			trip[i].ind = f_index[proc][i];
			trip[i].N = N_reg;
			trip[i].h = h;
		}

		std::sort(trip, trip + qsize, ValueCmp);

		query_outside[proc].clear();
		f_index[proc].clear();

		for (int i = 0; i < qsize; ++i) {
			query_outside[proc].push_back(trip[i].x);
			query_outside[proc].push_back(trip[i].y);
			query_outside[proc].push_back(trip[i].z);
			f_index[proc].push_back(trip[i].ind);
		}
		delete (trip);
	}
	return;
}


#include "libmorton/libmorton/include/morton.h"
class zTrip {
public:
	zTrip() {
	}
	;
  int mortid_;
  int i_;

};
inline uint64_t mortonEncode_for(unsigned int x, unsigned int y, unsigned int z){
    uint64_t answer = 0;
    for (uint64_t i = 0; i < (sizeof(uint64_t)* CHAR_BIT)/3; ++i) {
        answer |= ((x & ((uint64_t)1 << i)) << 2*i) | ((y & ((uint64_t)1 << i)) << (2*i + 1)) | ((z & ((uint64_t)1 << i)) << (2*i + 2));
    }
    return answer;
}
static bool zValueCmp(zTrip const & a, zTrip const & b) {

  return (a.mortid_ < b.mortid_);
}

static void zsort_queries(std::vector<Real>* query_outside,
		std::vector<int>* f_index, int* N_reg, Real* h, MPI_Comm c_comm) {

	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);
  const Real h0 = h[0];
  const Real h1 = h[1];
  const Real h2 = h[2];
  int N_reg_g[3];
  N_reg_g[0] = N_reg[0] + 2 * 3;
  N_reg_g[1] = N_reg[1] + 2 * 3;
  N_reg_g[2] = N_reg[2] + 2 * 3;

#ifdef FAST_INTERP_BINNING
	for (int proc = 0; proc < nprocs; ++proc) {
		int qsize = query_outside[proc].size() / COORD_DIM;
    std::vector<Real> tmp_query(query_outside[proc]); // tol hold xyz coordinates
    std::vector<int> tmp_f_index(f_index[proc]); // tol hold xyz coordinates
    int bsize_xyz[COORD_DIM];
    const int bsize = 16;
    bsize_xyz[0] = std::ceil(N_reg[0] / (Real)bsize);
    bsize_xyz[1] = std::ceil(N_reg[1] / (Real)bsize);
    bsize_xyz[2] = std::ceil(N_reg[2] / (Real)bsize);
    const int total_bsize = bsize_xyz[0] * bsize_xyz[1] * bsize_xyz[2];
    std::vector<int> bins[total_bsize];
		zTrip* trip = new zTrip[total_bsize];
    //std::cout << "------------------ total bin size = " << total_bsize << std::endl;

	  Real* x_ptr = &query_outside[proc][0];
		for (int i = 0; i < qsize; ++i) {
			int x = (int) std::abs(std::floor(x_ptr[0] / h0) / bsize);
			int y = (int) std::abs(std::floor(x_ptr[1] / h1) / bsize);
			int z = (int) std::abs(std::floor(x_ptr[2] / h2) / bsize);
      int indx = z + y * bsize_xyz[2] + x * bsize_xyz[2] * bsize_xyz[1];
      //std::cout << "x = " << x << " y = " << y << " z = " << z << " indx = " << indx << std::endl;
      bins[indx].push_back(i);
      x_ptr += 3;
		}

    for(int i = 0; i < bsize_xyz[0]; ++i) {
    for(int j = 0; j < bsize_xyz[1]; ++j) {
    for(int k = 0; k < bsize_xyz[2]; ++k) {
      int indx = k + j * bsize_xyz[2] + i * bsize_xyz[2] * bsize_xyz[1];
      trip[indx].mortid_ =  morton3D_32_encode(i, j, k);
      trip[indx].i_ = indx;
    }
    }
    }

		std::sort(trip, trip + total_bsize, zValueCmp);

    double time = -MPI_Wtime();
		query_outside[proc].clear();
		f_index[proc].clear();
		// query_outside[proc].reserve(qsize * COORD_DIM);
		// f_index[proc].reserve(qsize);
		for (int i = 0; i < total_bsize; ++i) {
      // push back contents of bins based on the order given by trip[i].i_
      int bindx = trip[i].i_;
      if(!bins[bindx].empty()){
        // std::cout << "bindx = " << bindx << " mid = " << trip[i].mortid_<< std::endl;
        std::vector<Real>* query_ptr= &query_outside[proc];
        std::vector<int>* f_index_ptr= &f_index[proc];
        for(int j = 0; j < bins[bindx].size(); ++j){
          int ind_ = bins[bindx][j];
			    query_ptr->push_back(tmp_query[ind_ * COORD_DIM + 0]);
			    query_ptr->push_back(tmp_query[ind_ * COORD_DIM + 1]);
			    query_ptr->push_back(tmp_query[ind_ * COORD_DIM + 2]);
			    f_index_ptr->push_back(tmp_f_index[ind_]);
        }
      }
		}
    time+=MPI_Wtime();
    // std::cout << "*** time = " << time << std::endl;
		delete (trip);
	}
#else
	for (int proc = 0; proc < nprocs; ++proc) {
		int qsize = query_outside[proc].size() / COORD_DIM;
		zTrip* trip = new zTrip[qsize];
    std::vector<Real> tmp_query(query_outside[proc]); // tol hold xyz coordinates
    std::vector<int> tmp_f_index(f_index[proc]); // tol hold xyz coordinates

	  //double* x_ptr = &query_outside[proc][i * COORD_DIM + 0];
	  Real* x_ptr = &query_outside[proc][0];
		for (int i = 0; i < qsize; ++i) {
			int x = (int) std::abs(std::floor(x_ptr[0] / h0));
			int y = (int) std::abs(std::floor(x_ptr[1] / h1));
			int z = (int) std::abs(std::floor(x_ptr[2] / h2));
      if(0){
        trip[i].mortid_ =  mortonEncode_for(x, y, z);
      }
      else{
        //trip[i].mortid_ =  mortonEncode_LUT(x, y, z);
        trip[i].mortid_ =  morton3D_32_encode(x, y, z);
        //trip[i].mortid_ =  mortonEncode_LUT(x, y, z);
      }
			trip[i].i_ = i;
      x_ptr += 3;
		}

		std::sort(trip, trip + qsize, zValueCmp);

		query_outside[proc].clear();
		f_index[proc].clear();

		for (int i = 0; i < qsize; ++i) {
      // std::cout << "bindx = " << trip[i].i_ << " mid = " << trip[i].mortid_<< std::endl;
			query_outside[proc].push_back(tmp_query[trip[i].i_ * COORD_DIM + 0]);
			query_outside[proc].push_back(tmp_query[trip[i].i_ * COORD_DIM + 1]);
			query_outside[proc].push_back(tmp_query[trip[i].i_ * COORD_DIM + 2]);
			f_index[proc].push_back(tmp_f_index[trip[i].i_]);
		}
		delete (trip);
	}
#endif
	return;
}
/*
 * Phase 1 of the parallel interpolation: This function computes which query_points needs to be sent to
 * other processors and which ones can be interpolated locally. Then a sparse alltoall is performed and
 * all the necessary information is sent/received including the coordinates for the query_points.
 * At the end, each process has the coordinates for interpolation of its own data and those of the others.
 *
 * IMPORTANT: This function must be called just once for a specific query_points. The reason is because of the
 * optimizations performed which assumes that the query_points do not change. For repeated interpolation you should
 * just call this function once, and instead repeatedly call Interp3_Plan::interpolate function.
 */
void Interp3_Plan::fast_scatter(int data_dof, int* N_reg, int * isize, int* istart,
		const int N_pts, const int g_size, Real* query_points_in, int* c_dims,
		MPI_Comm c_comm, double * timings) {
	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);

	if (this->allocate_baked == false) {
		std::cout
				<< "ERROR Interp3_Plan Scatter called before calling allocate.\n";
		return;
	}
	if (this->scatter_baked == true) {
		for (int proc = 0; proc < nprocs; ++proc) {
			std::vector<int>().swap(f_index[proc]);
			std::vector<Real>().swap(query_outside[proc]);
		}
	}
	all_query_points_allocation = 0;

	{

		//int N_reg_g[3], isize_g[3];
		N_reg_g[0] = N_reg[0] + 2 * g_size;
		N_reg_g[1] = N_reg[1] + 2 * g_size;
		N_reg_g[2] = N_reg[2] + 2 * g_size;

		isize_g[0] = isize[0] + 2 * g_size;
		isize_g[1] = isize[1] + 2 * g_size;
		isize_g[2] = isize[2] + 2 * g_size;

		Real h[3]; // original grid size along each axis
		h[0] = 1. / N_reg[0];
		h[1] = 1. / N_reg[1];
		h[2] = 1. / N_reg[2];

		// We copy query_points_in to query_points to aviod overwriting the input coordinates
		Real* query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		// Enforce periodicity
		for (int i = 0; i < N_pts; i++) {
			while (query_points[i * COORD_DIM + 0] <= -h[0]) {
				query_points[i * COORD_DIM + 0] =
						query_points[i * COORD_DIM + 0] + 1;
			}
			while (query_points[i * COORD_DIM + 1] <= -h[1]) {
				query_points[i * COORD_DIM + 1] =
						query_points[i * COORD_DIM + 1] + 1;
			}
			while (query_points[i * COORD_DIM + 2] <= -h[2]) {
				query_points[i * COORD_DIM + 2] =
						query_points[i * COORD_DIM + 2] + 1;
			}

			while (query_points[i * COORD_DIM + 0] >= 1) {
				query_points[i * COORD_DIM + 0] =
						query_points[i * COORD_DIM + 0] - 1;
			}
			while (query_points[i * COORD_DIM + 1] >= 1) {
				query_points[i * COORD_DIM + 1] =
						query_points[i * COORD_DIM + 1] - 1;
			}
			while (query_points[i * COORD_DIM + 2] >= 1) {
				query_points[i * COORD_DIM + 2] =
						query_points[i * COORD_DIM + 2] - 1;
			}
		}

		// Compute the start and end coordinates that this processor owns
		Real iX0[3], iX1[3];
		for (int j = 0; j < 3; j++) {
			iX0[j] = istart[j] * h[j];
			iX1[j] = iX0[j] + (isize[j] - 1) * h[j];
		}

		// Now march through the query points and split them into nprocs parts.
		// These are stored in query_outside which is an array of vectors of size nprocs.
		// That is query_outside[i] is a vector that contains the query points that need to
		// be sent to process i. Obviously for the case of query_outside[procid], we do not
		// need to send it to any other processor, as we own the necessary information locally,
		// and interpolation can be done locally.
		int Q_local = 0, Q_outside = 0;

		// This is needed for one-to-one correspondence with output f. This is becaues we are reshuffling
		// the data according to which processor it land onto, and we need to somehow keep the original
		// index to write the interpolation data back to the right location in the output.

		// This is necessary because when we want to compute dproc0 and dproc1 we have to divide by
		// the max isize. If the proc grid is unbalanced, the last proc's isize will be different
		// than others. With this approach we always use the right isize0 for all procs.
		int isize0 = std::ceil(N_reg[0] * 1. / c_dims[0]);
		int isize1 = std::ceil(N_reg[1] * 1. / c_dims[1]);
		for (int i = 0; i < N_pts; i++) {
			// The if condition checks whether the query points fall into the locally owned domain or not
			if (iX0[0] - h[0] <= query_points[i * COORD_DIM + 0]
					&& query_points[i * COORD_DIM + 0] <= iX1[0] + h[0]
					&& iX0[1] - h[1] <= query_points[i * COORD_DIM + 1]
					&& query_points[i * COORD_DIM + 1] <= iX1[1] + h[1]
					&& iX0[2] - h[2] <= query_points[i * COORD_DIM + 2]
					&& query_points[i * COORD_DIM + 2] <= iX1[2] + h[2]) {
				query_outside[procid].push_back(
						query_points[i * COORD_DIM + 0]);
				query_outside[procid].push_back(
						query_points[i * COORD_DIM + 1]);
				query_outside[procid].push_back(
						query_points[i * COORD_DIM + 2]);
				f_index[procid].push_back(i);
				++Q_local;
				//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
				continue;
			} else {
				// If the point does not reside in the processor's domain then we have to
				// first compute which processor owns the point. After computing that
				// we add the query point to the corresponding vector.
				int dproc0 = (int) (query_points[i * COORD_DIM + 0] / h[0])
						/ isize0;
				int dproc1 = (int) (query_points[i * COORD_DIM + 1] / h[1])
						/ isize1;
				int proc = dproc0 * c_dims[1] + dproc1; // Compute which proc has to do the interpolation
				//PCOUT<<"proc="<<proc<<std::endl;
				query_outside[proc].push_back(query_points[i * COORD_DIM + 0]);
				query_outside[proc].push_back(query_points[i * COORD_DIM + 1]);
				query_outside[proc].push_back(query_points[i * COORD_DIM + 2]);
				f_index[proc].push_back(i);
				++Q_outside;
				//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
				continue;
			}

		}

		// Now sort the query points in zyx order
#ifdef SORT_QUERIES
		timings[3]+=-MPI_Wtime();
		zsort_queries(query_outside,f_index,N_reg,h,c_comm);
		timings[3]+=+MPI_Wtime();
		//if(procid==0) std::cout<<"Sorting time="<<s_time<<std::endl;;
		//if(procid==0) std::cout<<"Sorting Queries\n";
#endif

		// Now we need to send the query_points that land onto other processor's domain.
		// This done using a sparse alltoallv.
		// Right now each process knows how much data to send to others, but does not know
		// how much data it should receive. This is a necessary information both for the MPI
		// command as well as memory allocation for received data.
		// So we first do an alltoall to get the f_index[proc].size from all processes.

		for (int proc = 0; proc < nprocs; proc++) {
			if (!f_index[proc].empty())
				f_index_procs_self_sizes[proc] = f_index[proc].size();
			else
				f_index_procs_self_sizes[proc] = 0;
		}
		timings[0] += -MPI_Wtime();
		MPI_Alltoall(f_index_procs_self_sizes, 1, MPI_INT,
				f_index_procs_others_sizes, 1, MPI_INT, c_comm);
		timings[0] += +MPI_Wtime();

		// Now we need to allocate memory for the receiving buffer of all query
		// points including ours. This is simply done by looping through
		// f_index_procs_others_sizes and adding up all the sizes.
		// Note that we would also need to know the offsets.
		f_index_procs_others_offset[0] = 0;
		f_index_procs_self_offset[0] = 0;
		for (int proc = 0; proc < nprocs; ++proc) {
			// The reason we multiply by COORD_DIM is that we have three coordinates per interpolation request
			all_query_points_allocation += f_index_procs_others_sizes[proc]
					* COORD_DIM;
			if (proc > 0) {
				f_index_procs_others_offset[proc] =
						f_index_procs_others_offset[proc - 1]
								+ f_index_procs_others_sizes[proc - 1];
				f_index_procs_self_offset[proc] = f_index_procs_self_offset[proc
						- 1] + f_index_procs_self_sizes[proc - 1];
			}
		}
		total_query_points = all_query_points_allocation / COORD_DIM;

		// This if condition is to allow multiple calls to scatter fucntion with different query points
		// without having to create a new plan
		if (this->scatter_baked == true) {
			free(this->all_query_points);
			free(this->all_f_cubic);
			all_query_points = (Real*) malloc(
					all_query_points_allocation * sizeof(Real));
			all_f_cubic = (Real*) malloc(
					total_query_points * sizeof(Real) * data_dof);
		} else {
			all_query_points = (Real*) malloc(
					all_query_points_allocation * sizeof(Real));
			all_f_cubic = (Real*) malloc(
					total_query_points * sizeof(Real) * data_dof);
		}

		// Now perform the allotall to send/recv query_points
		timings[0] += -MPI_Wtime();
		{
			int dst_r, dst_s;
			for (int i = 0; i < nprocs; ++i) {
				dst_r = i;    //(procid+i)%nprocs;
				dst_s = i;    //(procid-i+nprocs)%nprocs;
				s_request[dst_s] = MPI_REQUEST_NULL;
				request[dst_r] = MPI_REQUEST_NULL;
				int roffset = f_index_procs_others_offset[dst_r] * COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
				int soffset = f_index_procs_self_offset[dst_s] * COORD_DIM;
				if (f_index_procs_others_sizes[dst_r] != 0)
					MPI_Irecv(&all_query_points[roffset],
							f_index_procs_others_sizes[dst_r] * COORD_DIM,
							MPI_T, dst_r, 0, c_comm, &request[dst_r]);
				if (!query_outside[dst_s].empty())
					MPI_Isend(&query_outside[dst_s][0],
							f_index_procs_self_sizes[dst_s] * COORD_DIM, MPI_T,
							dst_s, 0, c_comm, &s_request[dst_s]);
			}
			// Wait for all the communication to finish
			MPI_Status ierr;
			for (int proc = 0; proc < nprocs; ++proc) {
				if (request[proc] != MPI_REQUEST_NULL)
					MPI_Wait(&request[proc], &ierr);
				if (s_request[proc] != MPI_REQUEST_NULL)
					MPI_Wait(&s_request[proc], &ierr);
			}
		}
		timings[0] += +MPI_Wtime();

		// Now perform the interpolation on all query points including those that need to
		// be sent to other processors and store them into all_f_cubic
		free(query_points);
	}

	for (int i = 0; i < nprocs; ++i) {
		MPI_Type_vector(data_dof, f_index_procs_self_sizes[i], N_pts, MPI_T,
				&rtype[i]);
		MPI_Type_vector(data_dof, f_index_procs_others_sizes[i],
				total_query_points, MPI_T, &stype[i]);
		MPI_Type_commit(&stype[i]);
		MPI_Type_commit(&rtype[i]);
	}

	rescale_xyz(g_size, N_reg, N_reg_g, istart, total_query_points,
			all_query_points);
	this->scatter_baked = true;
	return;
}

void Interp3_Plan::scatter(int data_dof, int* N_reg, int * isize, int* istart,
		const int N_pts, const int g_size, Real* query_points_in, int* c_dims,
		MPI_Comm c_comm, double * timings) {
#ifdef FAST_INTERP
  return fast_scatter(data_dof, N_reg, isize, istart, N_pts,
      g_size, query_points_in, c_dims, c_comm, timings);
#else
	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);

	if (this->allocate_baked == false) {
		std::cout
				<< "ERROR Interp3_Plan Scatter called before calling allocate.\n";
		return;
	}
	if (this->scatter_baked == true) {
		for (int proc = 0; proc < nprocs; ++proc) {
			std::vector<int>().swap(f_index[proc]);
			std::vector<Real>().swap(query_outside[proc]);
		}
	}
	all_query_points_allocation = 0;

	{

		//int N_reg_g[3], isize_g[3];
		N_reg_g[0] = N_reg[0] + 2 * g_size;
		N_reg_g[1] = N_reg[1] + 2 * g_size;
		N_reg_g[2] = N_reg[2] + 2 * g_size;

		isize_g[0] = isize[0] + 2 * g_size;
		isize_g[1] = isize[1] + 2 * g_size;
		isize_g[2] = isize[2] + 2 * g_size;

		Real h[3]; // original grid size along each axis
		h[0] = 1. / N_reg[0];
		h[1] = 1. / N_reg[1];
		h[2] = 1. / N_reg[2];

		// We copy query_points_in to query_points to aviod overwriting the input coordinates
		Real* query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		// Enforce periodicity
		for (int i = 0; i < N_pts; i++) {
			while (query_points[i * COORD_DIM + 0] <= -h[0]) {
				query_points[i * COORD_DIM + 0] =
						query_points[i * COORD_DIM + 0] + 1;
			}
			while (query_points[i * COORD_DIM + 1] <= -h[1]) {
				query_points[i * COORD_DIM + 1] =
						query_points[i * COORD_DIM + 1] + 1;
			}
			while (query_points[i * COORD_DIM + 2] <= -h[2]) {
				query_points[i * COORD_DIM + 2] =
						query_points[i * COORD_DIM + 2] + 1;
			}

			while (query_points[i * COORD_DIM + 0] >= 1) {
				query_points[i * COORD_DIM + 0] =
						query_points[i * COORD_DIM + 0] - 1;
			}
			while (query_points[i * COORD_DIM + 1] >= 1) {
				query_points[i * COORD_DIM + 1] =
						query_points[i * COORD_DIM + 1] - 1;
			}
			while (query_points[i * COORD_DIM + 2] >= 1) {
				query_points[i * COORD_DIM + 2] =
						query_points[i * COORD_DIM + 2] - 1;
			}
		}

		// Compute the start and end coordinates that this processor owns
		Real iX0[3], iX1[3];
		for (int j = 0; j < 3; j++) {
			iX0[j] = istart[j] * h[j];
			iX1[j] = iX0[j] + (isize[j] - 1) * h[j];
		}

		// Now march through the query points and split them into nprocs parts.
		// These are stored in query_outside which is an array of vectors of size nprocs.
		// That is query_outside[i] is a vector that contains the query points that need to
		// be sent to process i. Obviously for the case of query_outside[procid], we do not
		// need to send it to any other processor, as we own the necessary information locally,
		// and interpolation can be done locally.
		int Q_local = 0, Q_outside = 0;

		// This is needed for one-to-one correspondence with output f. This is becaues we are reshuffling
		// the data according to which processor it land onto, and we need to somehow keep the original
		// index to write the interpolation data back to the right location in the output.

		// This is necessary because when we want to compute dproc0 and dproc1 we have to divide by
		// the max isize. If the proc grid is unbalanced, the last proc's isize will be different
		// than others. With this approach we always use the right isize0 for all procs.
		int isize0 = std::ceil(N_reg[0] * 1. / c_dims[0]);
		int isize1 = std::ceil(N_reg[1] * 1. / c_dims[1]);
		for (int i = 0; i < N_pts; i++) {
			// The if condition checks whether the query points fall into the locally owned domain or not
			if (iX0[0] - h[0] <= query_points[i * COORD_DIM + 0]
					&& query_points[i * COORD_DIM + 0] <= iX1[0] + h[0]
					&& iX0[1] - h[1] <= query_points[i * COORD_DIM + 1]
					&& query_points[i * COORD_DIM + 1] <= iX1[1] + h[1]
					&& iX0[2] - h[2] <= query_points[i * COORD_DIM + 2]
					&& query_points[i * COORD_DIM + 2] <= iX1[2] + h[2]) {
				query_outside[procid].push_back(
						query_points[i * COORD_DIM + 0]);
				query_outside[procid].push_back(
						query_points[i * COORD_DIM + 1]);
				query_outside[procid].push_back(
						query_points[i * COORD_DIM + 2]);
				f_index[procid].push_back(i);
				++Q_local;
				//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
				continue;
			} else {
				// If the point does not reside in the processor's domain then we have to
				// first compute which processor owns the point. After computing that
				// we add the query point to the corresponding vector.
				int dproc0 = (int) (query_points[i * COORD_DIM + 0] / h[0])
						/ isize0;
				int dproc1 = (int) (query_points[i * COORD_DIM + 1] / h[1])
						/ isize1;
				int proc = dproc0 * c_dims[1] + dproc1; // Compute which proc has to do the interpolation
				//PCOUT<<"proc="<<proc<<std::endl;
				query_outside[proc].push_back(query_points[i * COORD_DIM + 0]);
				query_outside[proc].push_back(query_points[i * COORD_DIM + 1]);
				query_outside[proc].push_back(query_points[i * COORD_DIM + 2]);
				f_index[proc].push_back(i);
				++Q_outside;
				//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
				continue;
			}

		}

		// Now sort the query points in zyx order
#ifdef SORT_QUERIES
		timings[3]+=-MPI_Wtime();
		sort_queries(query_outside,f_index,N_reg,h,c_comm);
		timings[3]+=+MPI_Wtime();
		//if(procid==0) std::cout<<"Sorting time="<<s_time<<std::endl;;
		//if(procid==0) std::cout<<"Sorting Queries\n";
#endif

		// Now we need to send the query_points that land onto other processor's domain.
		// This done using a sparse alltoallv.
		// Right now each process knows how much data to send to others, but does not know
		// how much data it should receive. This is a necessary information both for the MPI
		// command as well as memory allocation for received data.
		// So we first do an alltoall to get the f_index[proc].size from all processes.

		for (int proc = 0; proc < nprocs; proc++) {
			if (!f_index[proc].empty())
				f_index_procs_self_sizes[proc] = f_index[proc].size();
			else
				f_index_procs_self_sizes[proc] = 0;
		}
		timings[0] += -MPI_Wtime();
		MPI_Alltoall(f_index_procs_self_sizes, 1, MPI_INT,
				f_index_procs_others_sizes, 1, MPI_INT, c_comm);
		timings[0] += +MPI_Wtime();

		// Now we need to allocate memory for the receiving buffer of all query
		// points including ours. This is simply done by looping through
		// f_index_procs_others_sizes and adding up all the sizes.
		// Note that we would also need to know the offsets.
		f_index_procs_others_offset[0] = 0;
		f_index_procs_self_offset[0] = 0;
		for (int proc = 0; proc < nprocs; ++proc) {
			// The reason we multiply by COORD_DIM is that we have three coordinates per interpolation request
			all_query_points_allocation += f_index_procs_others_sizes[proc]
					* COORD_DIM;
			if (proc > 0) {
				f_index_procs_others_offset[proc] =
						f_index_procs_others_offset[proc - 1]
								+ f_index_procs_others_sizes[proc - 1];
				f_index_procs_self_offset[proc] = f_index_procs_self_offset[proc
						- 1] + f_index_procs_self_sizes[proc - 1];
			}
		}
		total_query_points = all_query_points_allocation / COORD_DIM;

		// This if condition is to allow multiple calls to scatter fucntion with different query points
		// without having to create a new plan
		if (this->scatter_baked == true) {
			free(this->all_query_points);
			free(this->all_f_cubic);
			all_query_points = (Real*) malloc(
					all_query_points_allocation * sizeof(Real));
			all_f_cubic = (Real*) malloc(
					total_query_points * sizeof(Real) * data_dof);
		} else {
			all_query_points = (Real*) malloc(
					all_query_points_allocation * sizeof(Real));
			all_f_cubic = (Real*) malloc(
					total_query_points * sizeof(Real) * data_dof);
		}

		// Now perform the allotall to send/recv query_points
		timings[0] += -MPI_Wtime();
		{
			int dst_r, dst_s;
			for (int i = 0; i < nprocs; ++i) {
				dst_r = i;    //(procid+i)%nprocs;
				dst_s = i;    //(procid-i+nprocs)%nprocs;
				s_request[dst_s] = MPI_REQUEST_NULL;
				request[dst_r] = MPI_REQUEST_NULL;
				int roffset = f_index_procs_others_offset[dst_r] * COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
				int soffset = f_index_procs_self_offset[dst_s] * COORD_DIM;
				if (f_index_procs_others_sizes[dst_r] != 0)
					MPI_Irecv(&all_query_points[roffset],
							f_index_procs_others_sizes[dst_r] * COORD_DIM,
							MPI_T, dst_r, 0, c_comm, &request[dst_r]);
				if (!query_outside[dst_s].empty())
					MPI_Isend(&query_outside[dst_s][0],
							f_index_procs_self_sizes[dst_s] * COORD_DIM, MPI_T,
							dst_s, 0, c_comm, &s_request[dst_s]);
			}
			// Wait for all the communication to finish
			MPI_Status ierr;
			for (int proc = 0; proc < nprocs; ++proc) {
				if (request[proc] != MPI_REQUEST_NULL)
					MPI_Wait(&request[proc], &ierr);
				if (s_request[proc] != MPI_REQUEST_NULL)
					MPI_Wait(&s_request[proc], &ierr);
			}
		}
		timings[0] += +MPI_Wtime();

		// Now perform the interpolation on all query points including those that need to
		// be sent to other processors and store them into all_f_cubic
		free(query_points);
	}

	for (int i = 0; i < nprocs; ++i) {
		MPI_Type_vector(data_dof, f_index_procs_self_sizes[i], N_pts, MPI_T,
				&rtype[i]);
		MPI_Type_vector(data_dof, f_index_procs_others_sizes[i],
				total_query_points, MPI_T, &stype[i]);
		MPI_Type_commit(&stype[i]);
		MPI_Type_commit(&rtype[i]);
	}

	rescale_xyz(g_size, N_reg, N_reg_g, istart, total_query_points,
			all_query_points);
	this->scatter_baked = true;
	return;
#endif
}

/*
 * Phase 2 of the parallel interpolation: This function must be called after the scatter function is called.
 * It performs local interpolation for all the points that the processor has for itself, as well as the interpolations
 * that it has to send to other processors. After the local interpolation is performed, a sparse
 * alltoall is performed so that all the interpolated results are sent/received.
 *
 */
void Interp3_Plan::interpolate(Real* ghost_reg_grid_vals, int data_dof,
		int* N_reg, int * isize, int* istart, const int N_pts, const int g_size,
		Real* query_values, int* c_dims, MPI_Comm c_comm, double * timings) {
	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);
	if (this->allocate_baked == false) {
		std::cout
				<< "ERROR Interp3_Plan interpolate called before calling allocate.\n";
		return;
	}
	if (this->scatter_baked == false) {
		std::cout
				<< "ERROR Interp3_Plan interpolate called before calling scatter.\n";
		return;
	}

	timings[1] += -MPI_Wtime();
#ifdef FAST_INTERP
	optimized_interp3_ghost_xyz_p(ghost_reg_grid_vals, data_dof, N_reg, N_reg_g, isize_g,
			istart, total_query_points, g_size, all_query_points, all_f_cubic,
			true);
#else
	 interp3_ghost_xyz_p(ghost_reg_grid_vals, data_dof, N_reg, N_reg_g, isize_g,
			istart, total_query_points, g_size, all_query_points, all_f_cubic,
			true);
#endif
	timings[1] += +MPI_Wtime();

	// Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to
	// f_cubic_unordered.
	timings[0] += -MPI_Wtime();
	{
		int dst_r, dst_s;
		for (int i = 0; i < nprocs; ++i) {
			dst_r = i;  //(procid+i)%nprocs;
			dst_s = i;  //(procid-i+nprocs)%nprocs;
			s_request[dst_s] = MPI_REQUEST_NULL;
			request[dst_r] = MPI_REQUEST_NULL;
			// Notice that this is the adjoint of the first comm part
			// because now you are sending others f and receiving your part of f
			int soffset = f_index_procs_others_offset[dst_r];
			int roffset = f_index_procs_self_offset[dst_s];

			if (f_index_procs_self_sizes[dst_r] != 0)
				MPI_Irecv(&f_cubic_unordered[roffset], 1, rtype[i], dst_r, 0,
						c_comm, &request[dst_r]);
			if (f_index_procs_others_sizes[dst_s] != 0)
				MPI_Isend(&all_f_cubic[soffset], 1, stype[i], dst_s, 0, c_comm,
						&s_request[dst_s]);
		}
		MPI_Status ierr;
		for (int proc = 0; proc < nprocs; ++proc) {
			if (request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&request[proc], &ierr);
			if (s_request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&s_request[proc], &ierr);
		}
	}
	timings[0] += +MPI_Wtime();

	// Now copy back f_cubic_unordered to f_cubic in the correct f_index
	for (int dof = 0; dof < data_dof; ++dof) {
		for (int proc = 0; proc < nprocs; ++proc) {
			if (!f_index[proc].empty())
				for (int i = 0; i < f_index[proc].size(); ++i) {
					int ind = f_index[proc][i];
					//f_cubic[ind]=all_f_cubic[f_index_procs_others_offset[proc]+i];
					query_values[ind + dof * N_pts] =
							f_cubic_unordered[f_index_procs_self_offset[proc]
									+ i + dof * N_pts];
				}
		}
	}

	return;
}

void Interp3_Plan::high_order_interpolate(Real* ghost_reg_grid_vals, int data_dof,
		int* N_reg, int * isize, int* istart, const int N_pts, const int g_size,
		Real* query_values, int* c_dims, MPI_Comm c_comm, double * timings, int interp_order) {
	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);
	if (this->allocate_baked == false) {
		std::cout
				<< "ERROR Interp3_Plan interpolate called before calling allocate.\n";
		return;
	}
	if (this->scatter_baked == false) {
		std::cout
				<< "ERROR Interp3_Plan interpolate called before calling scatter.\n";
		return;
	}

	timings[1] += -MPI_Wtime();
	interp3_ghost_xyz_p(ghost_reg_grid_vals, data_dof, N_reg, N_reg_g, isize_g,
			istart, total_query_points, g_size, all_query_points, all_f_cubic, interp_order,
			true);
	timings[1] += +MPI_Wtime();

	// Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to
	// f_cubic_unordered.
	timings[0] += -MPI_Wtime();
	{
		int dst_r, dst_s;
		for (int i = 0; i < nprocs; ++i) {
			dst_r = i;  //(procid+i)%nprocs;
			dst_s = i;  //(procid-i+nprocs)%nprocs;
			s_request[dst_s] = MPI_REQUEST_NULL;
			request[dst_r] = MPI_REQUEST_NULL;
			// Notice that this is the adjoint of the first comm part
			// because now you are sending others f and receiving your part of f
			int soffset = f_index_procs_others_offset[dst_r];
			int roffset = f_index_procs_self_offset[dst_s];

			if (f_index_procs_self_sizes[dst_r] != 0)
				MPI_Irecv(&f_cubic_unordered[roffset], 1, rtype[i], dst_r, 0,
						c_comm, &request[dst_r]);
			if (f_index_procs_others_sizes[dst_s] != 0)
				MPI_Isend(&all_f_cubic[soffset], 1, stype[i], dst_s, 0, c_comm,
						&s_request[dst_s]);
		}
		MPI_Status ierr;
		for (int proc = 0; proc < nprocs; ++proc) {
			if (request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&request[proc], &ierr);
			if (s_request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&s_request[proc], &ierr);
		}
	}
	timings[0] += +MPI_Wtime();

	// Now copy back f_cubic_unordered to f_cubic in the correct f_index
	for (int dof = 0; dof < data_dof; ++dof) {
		for (int proc = 0; proc < nprocs; ++proc) {
			if (!f_index[proc].empty())
				for (int i = 0; i < f_index[proc].size(); ++i) {
					int ind = f_index[proc][i];
					//f_cubic[ind]=all_f_cubic[f_index_procs_others_offset[proc]+i];
					query_values[ind + dof * N_pts] =
							f_cubic_unordered[f_index_procs_self_offset[proc]
									+ i + dof * N_pts];
				}
		}
	}

	return;
}


Interp3_Plan::~Interp3_Plan() {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (this->allocate_baked) {
		free(query_points);

		free(f_index_procs_others_offset);
		free(f_index_procs_self_offset);
		free(f_index_procs_self_sizes);
		free(f_index_procs_others_sizes);

		delete (s_request);
		delete (request);
		//vectors
		for (int proc = 0; proc < nprocs; ++proc) {
			std::vector<int>().swap(f_index[proc]);
			std::vector<Real>().swap(query_outside[proc]);
		}
		free(f_cubic_unordered);

	}

	if (this->scatter_baked) {
		for (int i = 0; i < nprocs; ++i) {
			MPI_Type_free(&stype[i]);
			MPI_Type_free(&rtype[i]);
		}
		free(all_query_points);
		free(all_f_cubic);
	}

	if (this->allocate_baked) {
		delete (stype);
		delete (rtype);
	}
	return;
}

/*
 * A dummy function that performs the whole process of scatter and interpolation without any planning or
 * prior allocation. Use only for debugging.
 */
void Interp3_Plan::slow_run(Real* ghost_reg_grid_vals, int data_dof, int* N_reg,
		int * isize, int* istart, const int N_pts, const int g_size,
		Real* query_points_in, Real* query_values, int* c_dims,
		MPI_Comm c_comm) {
	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);

	int N_reg_g[3], isize_g[3];
	N_reg_g[0] = N_reg[0] + 2 * g_size;
	N_reg_g[1] = N_reg[1] + 2 * g_size;
	N_reg_g[2] = N_reg[2] + 2 * g_size;

	isize_g[0] = isize[0] + 2 * g_size;
	isize_g[1] = isize[1] + 2 * g_size;
	isize_g[2] = isize[2] + 2 * g_size;

	Real h[3]; // original grid size along each axis
	h[0] = 1. / N_reg[0];
	h[1] = 1. / N_reg[1];
	h[2] = 1. / N_reg[2];

	// We copy query_points_in to query_points to aviod overwriting the input coordinates
	Real* query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
	memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
	// Enforce periodicity
	for (int i = 0; i < N_pts; i++) {
		while (query_points[i * COORD_DIM + 0] <= -h[0]) {
			query_points[i * COORD_DIM + 0] = query_points[i * COORD_DIM + 0]
					+ 1;
		}
		while (query_points[i * COORD_DIM + 1] <= -h[1]) {
			query_points[i * COORD_DIM + 1] = query_points[i * COORD_DIM + 1]
					+ 1;
		}
		while (query_points[i * COORD_DIM + 2] <= -h[2]) {
			query_points[i * COORD_DIM + 2] = query_points[i * COORD_DIM + 2]
					+ 1;
		}

		while (query_points[i * COORD_DIM + 0] >= 1) {
			query_points[i * COORD_DIM + 0] = query_points[i * COORD_DIM + 0]
					- 1;
		}
		while (query_points[i * COORD_DIM + 1] >= 1) {
			query_points[i * COORD_DIM + 1] = query_points[i * COORD_DIM + 1]
					- 1;
		}
		while (query_points[i * COORD_DIM + 2] >= 1) {
			query_points[i * COORD_DIM + 2] = query_points[i * COORD_DIM + 2]
					- 1;
		}
	}

	// Compute the start and end coordinates that this processor owns
	Real iX0[3], iX1[3];
	for (int j = 0; j < 3; j++) {
		iX0[j] = istart[j] * h[j];
		iX1[j] = iX0[j] + (isize[j] - 1) * h[j];
	}

	// Now march through the query points and split them into nprocs parts.
	// These are stored in query_outside which is an array of vectors of size nprocs.
	// That is query_outside[i] is a vector that contains the query points that need to
	// be sent to process i. Obviously for the case of query_outside[procid], we do not
	// need to send it to any other processor, as we own the necessary information locally,
	// and interpolation can be done locally.
	int Q_local = 0, Q_outside = 0;

	// This is needed for one-to-one correspondence with output f. This is becaues we are reshuffling
	// the data according to which processor it land onto, and we need to somehow keep the original
	// index to write the interpolation data back to the right location in the output.
	std::vector<int> f_index[nprocs];
	std::vector<Real> query_outside[nprocs];
	for (int i = 0; i < N_pts; i++) {
		// The if condition check whether the query points fall into the locally owned domain or not
		if (iX0[0] - h[0] <= query_points[i * COORD_DIM + 0]
				&& query_points[i * COORD_DIM + 0] <= iX1[0] + h[0]
				&& iX0[1] - h[1] <= query_points[i * COORD_DIM + 1]
				&& query_points[i * COORD_DIM + 1] <= iX1[1] + h[1]
				&& iX0[2] - h[2] <= query_points[i * COORD_DIM + 2]
				&& query_points[i * COORD_DIM + 2] <= iX1[2] + h[2]) {
			query_outside[procid].push_back(query_points[i * COORD_DIM + 0]);
			query_outside[procid].push_back(query_points[i * COORD_DIM + 1]);
			query_outside[procid].push_back(query_points[i * COORD_DIM + 2]);
			f_index[procid].push_back(i);
			Q_local++;
			//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
			continue;
		} else {
			// If the point does not reside in the processor's domain then we have to
			// first compute which processor owns the point. After computing that
			// we add the query point to the corresponding vector.
			int dproc0 = (int) (query_points[i * COORD_DIM + 0] / h[0])
					/ isize[0];
			int dproc1 = (int) (query_points[i * COORD_DIM + 1] / h[1])
					/ isize[1];
			int proc = dproc0 * c_dims[1] + dproc1; // Compute which proc has to do the interpolation
			//PCOUT<<"proc="<<proc<<std::endl;
			query_outside[proc].push_back(query_points[i * COORD_DIM + 0]);
			query_outside[proc].push_back(query_points[i * COORD_DIM + 1]);
			query_outside[proc].push_back(query_points[i * COORD_DIM + 2]);
			f_index[proc].push_back(i);
			Q_outside++;
			//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
			continue;
		}

	}

	// Now we need to send the query_points that land onto other processor's domain.
	// This done using a sparse alltoallv.
	// Right now each process knows how much data to send to others, but does not know
	// how much data it should receive. This is a necessary information both for the MPI
	// command as well as memory allocation for received data.
	// So we first do an alltoall to get the f_index[proc].size from all processes.
	int f_index_procs_self_sizes[nprocs]; // sizes of the number of interpolations that need to be sent to procs
	int f_index_procs_others_sizes[nprocs]; // sizes of the number of interpolations that need to be received from procs

	for (int proc = 0; proc < nprocs; proc++) {
		if (!f_index[proc].empty())
			f_index_procs_self_sizes[proc] = f_index[proc].size();
		else
			f_index_procs_self_sizes[proc] = 0;
	}
	MPI_Alltoall(f_index_procs_self_sizes, 1, MPI_INT,
			f_index_procs_others_sizes, 1, MPI_INT, c_comm);

#ifdef VERBOSE2
	sleep(1);
	if(procid==0) {
		std::cout<<"procid="<<procid<<std::endl;
		std::cout<<"f_index_procs_self[0]="<<f_index_procs_self_sizes[0]<<" [1]= "<<f_index_procs_self_sizes[1]<<std::endl;
		std::cout<<"f_index_procs_others[0]="<<f_index_procs_others_sizes[0]<<" [1]= "<<f_index_procs_others_sizes[1]<<std::endl;
	}
	sleep(1);
	if(procid==1) {
		std::cout<<"procid="<<procid<<std::endl;
		std::cout<<"f_index_procs_self[0]="<<f_index_procs_self_sizes[0]<<" [1]= "<<f_index_procs_self_sizes[1]<<std::endl;
		std::cout<<"f_index_procs_others[0]="<<f_index_procs_others_sizes[0]<<" [1]= "<<f_index_procs_others_sizes[1]<<std::endl;
	}
#endif

	// Now we need to allocate memory for the receiving buffer of all query
	// points including ours. This is simply done by looping through
	// f_index_procs_others_sizes and adding up all the sizes.
	// Note that we would also need to know the offsets.
	size_t all_query_points_allocation = 0;
	int f_index_procs_others_offset[nprocs]; // offset in the all_query_points array
	int f_index_procs_self_offset[nprocs]; // offset in the query_outside array
	f_index_procs_others_offset[0] = 0;
	f_index_procs_self_offset[0] = 0;
	for (int proc = 0; proc < nprocs; ++proc) {
		// The reason we multiply by COORD_DIM is that we have three coordinates per interpolation request
		all_query_points_allocation += f_index_procs_others_sizes[proc]
				* COORD_DIM;
		if (proc > 0) {
			f_index_procs_others_offset[proc] = f_index_procs_others_offset[proc
					- 1] + f_index_procs_others_sizes[proc - 1];
			f_index_procs_self_offset[proc] =
					f_index_procs_self_offset[proc - 1]
							+ f_index_procs_self_sizes[proc - 1];
		}
	}
	int total_query_points = all_query_points_allocation / COORD_DIM;
	Real * all_query_points = (Real*) malloc(
			all_query_points_allocation * sizeof(Real));
#ifdef VERBOSE2
	if(procid==0) {
		std::cout<<"procid="<<procid<<std::endl;
		for (int proc=0;proc<nprocs;++proc)
		std::cout<<"proc= "<<proc<<" others_offset= "<<f_index_procs_others_offset[proc]<<" others_sizes= "<<f_index_procs_others_sizes[proc]<<std::endl;
		for (int proc=0;proc<nprocs;++proc)
		std::cout<<"proc= "<<proc<<" self_offset= "<<f_index_procs_self_offset[proc]<<" self_sizes= "<<f_index_procs_self_sizes[proc]<<std::endl;
	}
#endif

	MPI_Request * s_request = new MPI_Request[nprocs];
	MPI_Request * request = new MPI_Request[nprocs];

	// Now perform the allotall to send/recv query_points
	{
		int dst_r, dst_s;
		for (int i = 0; i < nprocs; ++i) {
			dst_r = i; //(procid+i)%nprocs;
			dst_s = i; //(procid-i+nprocs)%nprocs;
			s_request[dst_s] = MPI_REQUEST_NULL;
			request[dst_r] = MPI_REQUEST_NULL;
			int roffset = f_index_procs_others_offset[dst_r] * COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
			int soffset = f_index_procs_self_offset[dst_s] * COORD_DIM;
			if (f_index_procs_others_sizes[dst_r] != 0)
				MPI_Irecv(&all_query_points[roffset],
						f_index_procs_others_sizes[dst_r] * COORD_DIM, MPI_T,
						dst_r, 0, c_comm, &request[dst_r]);
			if (!query_outside[dst_s].empty())
				MPI_Isend(&query_outside[dst_s][0],
						f_index_procs_self_sizes[dst_s] * COORD_DIM, MPI_T,
						dst_s, 0, c_comm, &s_request[dst_s]);
			//if(procid==1){
			//std::cout<<"soffset="<<soffset<<" roffset="<<roffset<<std::endl;
			//std::cout<<"f_index_procs_self_sizes[0]="<<f_index_procs_self_sizes[0]<<std::endl;
			//std::cout<<"f_index_procs_others_sizes[0]="<<f_index_procs_others_sizes[0]<<std::endl;
			//std::cout<<"q_outside["<<dst_s<<"]="<<query_outside[dst_s][0]<<std::endl;
			//}
		}
		// Wait for all the communication to finish
		MPI_Status ierr;
		for (int proc = 0; proc < nprocs; ++proc) {
			if (request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&request[proc], &ierr);
			if (s_request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&s_request[proc], &ierr);
		}
	}

	//if(procid==1){
	//  std::cout<<"total_query_points="<<total_query_points<<std::endl;
	//  std::cout<<"----- procid="<<procid<<" Q="<<all_query_points[0]<<" "<<all_query_points[1]<<" "<<all_query_points[2]<<std::endl;
	//  std::cout<<"----- procid="<<procid<<" Q="<<all_query_points[3]<<" "<<all_query_points[4]<<" "<<all_query_points[5]<<std::endl;
	//}
	//PCOUT<<"**** Q_local="<<Q_local<<" f_index_procs_self_sizes[procid]="<<f_index_procs_self_sizes[procid]<<std::endl;
	//int dum=0;
	//for (int i=0;i<Q_local;i++){
	//  dum+=query_local[i*COORD_DIM+0]-all_query_points[f_index_procs_others_offset[procid]*3+i*COORD_DIM+0];
	//  dum+=query_local[i*COORD_DIM+1]-all_query_points[f_index_procs_others_offset[procid]*3+i*COORD_DIM+1];
	//  dum+=query_local[i*COORD_DIM+2]-all_query_points[f_index_procs_others_offset[procid]*3+i*COORD_DIM+2];
	//}

	// Now perform the interpolation on all query points including those that need to
	// be sent to other processors and store them into all_f_cubic
	Real* all_f_cubic = (Real*) malloc(
			total_query_points * sizeof(Real) * data_dof);
	interp3_ghost_xyz_p(ghost_reg_grid_vals, data_dof, N_reg, N_reg_g, isize_g,
			istart, total_query_points, g_size, all_query_points, all_f_cubic);

	//if(procid==0){
	//  std::cout<<"total_query_points="<<total_query_points<<std::endl;
	//  std::cout<<"procid="<<procid<<" Q="<<all_query_points[0]<<" "<<all_query_points[1]<<" "<<all_query_points[2]<<" f= "<<all_f_cubic[0]<<std::endl;
	//  std::cout<<"procid="<<procid<<" Q="<<all_query_points[3]<<" "<<all_query_points[4]<<" "<<all_query_points[5]<<" f= "<<all_f_cubic[1]<<std::endl;
	//}

	// Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to
	// f_cubic_unordered.
	Real * f_cubic_unordered = (Real*) malloc(N_pts * sizeof(Real) * data_dof); // The reshuffled semi-final interpolated values are stored here
	{
		//PCOUT<<"total_query_points="<<total_query_points<<" N_pts="<<N_pts<<std::endl;
		int dst_r, dst_s;
		MPI_Datatype stype[nprocs], rtype[nprocs];
		for (int i = 0; i < nprocs; ++i) {
			MPI_Type_vector(data_dof, f_index_procs_self_sizes[i], N_pts, MPI_T,
					&rtype[i]);
			MPI_Type_vector(data_dof, f_index_procs_others_sizes[i],
					total_query_points, MPI_T, &stype[i]);
			MPI_Type_commit(&stype[i]);
			MPI_Type_commit(&rtype[i]);
		}
		for (int i = 0; i < nprocs; ++i) {
			dst_r = i; //(procid+i)%nprocs;
			dst_s = i; //(procid-i+nprocs)%nprocs;
			s_request[dst_s] = MPI_REQUEST_NULL;
			request[dst_r] = MPI_REQUEST_NULL;
			// Notice that this is the adjoint of the first comm part
			// because now you are sending others f and receiving your part of f
			int soffset = f_index_procs_others_offset[dst_r];
			int roffset = f_index_procs_self_offset[dst_s];
			//if(procid==0)
			//  std::cout<<"procid="<<procid<<" dst_s= "<<dst_s<<" soffset= "<<soffset<<" s_size="<<f_index_procs_others_sizes[dst_s]<<" dst_r= "<<dst_r<<" roffset="<<roffset<<" r_size="<<f_index_procs_self_sizes[dst_r]<<std::endl;
			//if(f_index_procs_self_sizes[dst_r]!=0)
			//  MPI_Irecv(&f_cubic_unordered[roffset],f_index_procs_self_sizes[dst_r],rtype, dst_r,
			//      0, c_comm, &request[dst_r]);
			//if(f_index_procs_others_sizes[dst_s]!=0)
			//  MPI_Isend(&all_f_cubic[soffset],f_index_procs_others_sizes[dst_s],stype,dst_s,
			//      0, c_comm, &s_request[dst_s]);
			//
			if (f_index_procs_self_sizes[dst_r] != 0)
				MPI_Irecv(&f_cubic_unordered[roffset], 1, rtype[i], dst_r, 0,
						c_comm, &request[dst_r]);
			if (f_index_procs_others_sizes[dst_s] != 0)
				MPI_Isend(&all_f_cubic[soffset], 1, stype[i], dst_s, 0, c_comm,
						&s_request[dst_s]);
		}
		MPI_Status ierr;
		for (int proc = 0; proc < nprocs; ++proc) {
			if (request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&request[proc], &ierr);
			if (s_request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&s_request[proc], &ierr);
		}
		for (int i = 0; i < nprocs; ++i) {
			MPI_Type_free(&stype[i]);
			MPI_Type_free(&rtype[i]);
		}
	}

	// Now copy back f_cubic_unordered to f_cubic in the correct f_index
	for (int dof = 0; dof < data_dof; ++dof) {
		for (int proc = 0; proc < nprocs; ++proc) {
			if (!f_index[proc].empty())
				for (int i = 0; i < f_index[proc].size(); ++i) {
					int ind = f_index[proc][i];
					//f_cubic[ind]=all_f_cubic[f_index_procs_others_offset[proc]+i];
					query_values[ind + dof * N_pts] =
							f_cubic_unordered[f_index_procs_self_offset[proc]
									+ i + dof * N_pts];
				}
		}
	}

	free(query_points);
	free(all_query_points);
	free(all_f_cubic);
	free(f_cubic_unordered);
	delete (s_request);
	delete (request);
	//vector
	for (int proc = 0; proc < nprocs; ++proc) {
		std::vector<int>().swap(f_index[proc]);
		std::vector<Real>().swap(query_outside[proc]);
	}
	return;
}
