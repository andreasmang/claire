

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
#ifdef __unix__
# include <unistd.h>
#elif defined _WIN32
# include <windows.h>
#define sleep(x) Sleep(1000 * x)
#endif

#ifdef INTERP_DEBUG
template <typename T>
void parallel_print(T in, const char* prefix) {

	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  PCOUT << prefix << std::endl;
  for(int proc = 0; proc < nprocs; ++proc) {
    if(procid == proc)
      std::cout << "proc = " << proc << "\t" << in << std::endl;
    else
      sleep(.8);
  }

}
#endif

Interp3_Plan::Interp3_Plan() {
	this->allocate_baked = false;
	this->scatter_baked = false;
  procs_i_recv_from_size_ = 0;
  procs_i_send_to_size_ = 0;
}

void Interp3_Plan::allocate(int N_pts, int* data_dofs, int nplans) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#ifdef INTERP_DEBUG
  PCOUT << "entered allocate\n";
#endif
	query_points = pvfmm::aligned_new<Real>(N_pts * COORD_DIM);
  pvfmm::memset(query_points,0, N_pts*COORD_DIM);

	f_index_procs_others_offset = pvfmm::aligned_new<int>(nprocs); // offset in the all_query_points array
	f_index_procs_self_offset   = pvfmm::aligned_new<int>(nprocs); // offset in the query_outside array
	f_index_procs_self_sizes    = pvfmm::aligned_new<int>(nprocs); // sizes of the number of interpolations that need to be sent to procs
	f_index_procs_others_sizes  = pvfmm::aligned_new<int>(nprocs); // sizes of the number of interpolations that need to be received from procs

	s_request = pvfmm::aligned_new<MPI_Request>(nprocs);
	request = pvfmm::aligned_new<MPI_Request>(2*nprocs);

	f_index = new std::vector<int>[nprocs];
	query_outside = new std::vector<Real>[nprocs];


  this->nplans_ = nplans; // number of reuses of the plan with the same scatter points
  this->data_dofs_ = pvfmm::aligned_new<int>(nplans);

  int max =0;
  for(int i = 0; i < nplans_; ++i) {
    max = std::max(max, data_dofs[i]);
    this->data_dofs_[i] = data_dofs[i];
  }
  this->data_dof_max = max;

	f_cubic_unordered = pvfmm::aligned_new<Real>(N_pts * data_dof_max); // The reshuffled semi-final interpolated values are stored here
  memset(&f_cubic_unordered[0],0, N_pts * sizeof(Real) * data_dof_max);

	stypes = pvfmm::aligned_new<MPI_Datatype>(nprocs*nplans_); // strided for multiple plan calls
	rtypes = pvfmm::aligned_new<MPI_Datatype>(nprocs*nplans_);
	this->allocate_baked = true;
#ifdef INTERP_DEBUG
  PCOUT << "allocate done\n";
#endif
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

#ifdef SORT_QUERIES
static bool ValueCmp(Trip const & a, Trip const & b) {
	return a.z + a.y / a.h[1] * a.N[2] + a.x / a.h[0] * a.N[1] * a.N[2]
			< b.z + b.y / b.h[1] * b.N[2] + b.x / b.h[0] * b.N[1] * b.N[2];
}
#endif

#ifdef SORT_QUERIES
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
		delete[] trip;
	}
	return;
}
#endif

#include "libmorton/libmorton/include/morton.h"

class zTrip {
public:
	zTrip() {
	}
	;
  int mortid_;
  int i_;// bin index

};
inline uint64_t mortonEncode_for(unsigned int x, unsigned int y, unsigned int z){
    uint64_t answer = 0;
    for (uint64_t i = 0; i < (sizeof(uint64_t)* CHAR_BIT)/3; ++i) {
        answer |= ((x & ((uint64_t)1 << i)) << 2*i) | ((y & ((uint64_t)1 << i)) << (2*i + 1)) | ((z & ((uint64_t)1 << i)) << (2*i + 2));
    }
    return answer;
}

#ifdef SORT_QUERIES
static bool zValueCmp(zTrip const & a, zTrip const & b) {
  return (a.mortid_ < b.mortid_);
}
#endif


#ifdef SORT_QUERIES
static void zsort_queries(std::vector<Real>* query_outside,
		std::vector<int>* f_index, int* N_reg, Real* h, MPI_Comm c_comm) {

	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);
  const Real h0 = h[0];
  const Real h1 = h[1];
  const Real h2 = h[2];

#ifdef FAST_INTERP_BINNING
  int bsize_xyz[COORD_DIM];
  const int bsize = 16;
  bsize_xyz[0] = std::ceil(N_reg[0] / (Real)bsize);
  bsize_xyz[1] = std::ceil(N_reg[1] / (Real)bsize);
  bsize_xyz[2] = std::ceil(N_reg[2] / (Real)bsize);
  const size_t total_bsize = bsize_xyz[0] * bsize_xyz[1] * bsize_xyz[2];
  pvfmm::Iterator<zTrip> trip = pvfmm::aligned_new<zTrip>(total_bsize);
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
	for (int proc = 0; proc < nprocs; ++proc) {
		const int qsize = query_outside[proc].size() / COORD_DIM;
    //std::cout << "------------------ total bin size = " << total_bsize << std::endl;

    // double time = -MPI_Wtime();
    std::vector<Real> bins_Q[total_bsize];
    std::vector<int> bins_f[total_bsize];
	  Real* x_ptr = &query_outside[proc][0];

		for (int i = 0; i < qsize; ++i) {
			const int x = (int) std::abs(std::floor(x_ptr[0] / h0) / bsize);
			const int y = (int) std::abs(std::floor(x_ptr[1] / h1) / bsize);
			const int z = (int) std::abs(std::floor(x_ptr[2] / h2) / bsize);
      const int indx = z + y * bsize_xyz[2] + x * bsize_xyz[2] * bsize_xyz[1];
      bins_Q[indx].push_back(x_ptr[0]);
      bins_Q[indx].push_back(x_ptr[1]);
      bins_Q[indx].push_back(x_ptr[2]);
      bins_f[indx].push_back(f_index[proc][i]);
      x_ptr += 3;
		}


		query_outside[proc].clear();
		f_index[proc].clear();
		for (int i = 0; i < (int)total_bsize; ++i) {
      int bindx = trip[i].i_;
      if(!bins_Q[bindx].empty()){
        query_outside[proc].insert(query_outside[proc].end(), bins_Q[bindx].begin(), bins_Q[bindx].end());
        f_index[proc].insert(f_index[proc].end(), bins_f[bindx].begin(), bins_f[bindx].end());
      }

    }
    // time+=MPI_Wtime();
    // std::cout << "*** time = " << time << std::endl;
	}
    // pvfmm::aligned_delete(bins_Q);
    // pvfmm::aligned_delete(bins_f);
  pvfmm::aligned_delete<zTrip>(trip);
	// delete[] trip;
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
		delete[] trip;
	}
#endif
	return;
}
#endif

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
void Interp3_Plan::fast_scatter(int* N_reg, int * isize, int* istart,
		const int N_pts, const int g_size, Real* query_points_in, int* c_dims,
		MPI_Comm c_comm, double * timings) {
	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);

#ifdef INTERP_DEBUG
  PCOUT << "entered fast scatter\n";
#endif
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
    double time = -MPI_Wtime();

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
		memcpy(&query_points[0], query_points_in, N_pts * COORD_DIM * sizeof(Real));
#ifdef INTERP_DEBUG
  PCOUT << "enforcing periodicity\n";
#endif
		// Enforce periodicity
#pragma omp parallel for
		for (int i = 0; i < N_pts; i++) {
      pvfmm::Iterator<Real> Q_ptr = query_points+(i * COORD_DIM);
      //Real* Q_ptr = &query_points[i * COORD_DIM];
			while (Q_ptr[0] <= -h[0]) {
				Q_ptr[0] +=  1;
			}
			while (Q_ptr[1] <= -h[1]) {
				Q_ptr[1] +=  1;
			}
			while (Q_ptr[2] <= -h[2]) {
				Q_ptr[2] +=  1;
			}
			while (Q_ptr[0] >= 1) {
				Q_ptr[0] += - 1;
			}
			while (Q_ptr[1] >= 1) {
				Q_ptr[1] += - 1;
			}
			while (Q_ptr[2] >= 1) {
				Q_ptr[2] += - 1;
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
#ifdef INTERP_DEBUG
  MPI_Barrier(c_comm);
  PCOUT << "sorting\n";
#endif

    //int bsize_xyz[COORD_DIM];
    //const int bsize = 16;
    //bsize_xyz[0] = std::ceil(N_reg[0] / (Real)bsize);
    //bsize_xyz[1] = std::ceil(N_reg[1] / (Real)bsize);
    //bsize_xyz[2] = std::ceil(N_reg[2] / (Real)bsize);
    //const size_t total_bsize = bsize_xyz[0] * bsize_xyz[1] * bsize_xyz[2];
    //pvfmm::Iterator<zTrip> trip = pvfmm::aligned_new<zTrip>(total_bsize);
    //pvfmm::Iterator<std::vector<Real>> bins_Q =
    //  pvfmm::aligned_new<std::vector<Real>>(total_bsize*nprocs);
    //pvfmm::Iterator<std::vector<int>> bins_f =
    //  pvfmm::aligned_new<std::vector<int>>(total_bsize*nprocs);
    //std::vector<Real> bins_Q[nprocs][total_bsize];
    //std::vector<int> bins_f[nprocs][total_bsize];
    //for(int i = 0; i < bsize_xyz[0]; ++i) {
    //for(int j = 0; j < bsize_xyz[1]; ++j) {
    //for(int k = 0; k < bsize_xyz[2]; ++k) {
    //  int indx = k + j * bsize_xyz[2] + i * bsize_xyz[2] * bsize_xyz[1];
    //  trip[indx].mortid_ =  morton3D_32_encode(i, j, k);
    //  trip[indx].i_ = indx;
    //}
    //}
    //}
		//std::sort(trip, trip + total_bsize, zValueCmp);

		for (int i = 0; i < N_pts; i++) {
      // Real* Q_ptr = &query_points[i * COORD_DIM];
      pvfmm::Iterator<Real> Q_ptr = query_points+(i * COORD_DIM);
			//const int x = (int) std::abs(std::floor(Q_ptr[0] / h[0]) / bsize);
			//const int y = (int) std::abs(std::floor(Q_ptr[1] / h[1]) / bsize);
			//const int z = (int) std::abs(std::floor(Q_ptr[2] / h[2]) / bsize);
      //const int indx = z + y * bsize_xyz[2] + x * bsize_xyz[2] * bsize_xyz[1];
			// The if condition checks whether the query points fall into the locally owned domain or not
			//if (iX0[0] - h[0] <= Q_ptr[0]
			//		&& Q_ptr[0] <= iX1[0] + h[0]
			//		&& iX0[1] - h[1] <= Q_ptr[1]
			//		&& Q_ptr[1] <= iX1[1] + h[1]
			//		&& iX0[2] - h[2] <= Q_ptr[2]
			//		&& Q_ptr[2] <= iX1[2] + h[2]) {
			if (iX0[0] - h[0] > Q_ptr[0]
					|| Q_ptr[0] > iX1[0] + h[0]
					|| iX0[1] - h[1] > Q_ptr[1]
					|| Q_ptr[1] > iX1[1] + h[1]
					|| iX0[2] - h[2] > Q_ptr[2]
					|| Q_ptr[2] > iX1[2] + h[2]) {
        // todo create a set for procs_i_communicate
				// If the point does not reside in the processor's domain then we have to
				// first compute which processor owns the point. After computing that
				// we add the query point to the corresponding vector.
				int dproc0 = (int) (Q_ptr[0] / h[0]) / isize0;
				int dproc1 = (int) (Q_ptr[1] / h[1]) / isize1;
				int proc = dproc0 * c_dims[1] + dproc1; // Compute which proc has to do the interpolation
				//PCOUT<<"proc="<<proc<<std::endl;
				query_outside[proc].push_back(Q_ptr[0]);
				query_outside[proc].push_back(Q_ptr[1]);
				query_outside[proc].push_back(Q_ptr[2]);
				f_index[proc].push_back(i);

        //bins_Q[proc*total_bsize+indx].push_back(Q_ptr[0]);
        //bins_Q[proc*total_bsize+indx].push_back(Q_ptr[1]);
        //bins_Q[proc*total_bsize+indx].push_back(Q_ptr[2]);
        //bins_f[proc*total_bsize+indx].push_back(i);
        // bins_Q[proc][indx].push_back(Q_ptr[0]);
        // bins_Q[proc][indx].push_back(Q_ptr[1]);
        // bins_Q[proc][indx].push_back(Q_ptr[2]);
        // bins_f[proc][indx].push_back(i);
				++Q_outside;
				//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
				continue;
			} else {
				query_outside[procid].push_back(Q_ptr[0]);
				query_outside[procid].push_back(Q_ptr[1]);
				query_outside[procid].push_back(Q_ptr[2]);
				f_index[procid].push_back(i);

        //bins_Q[procid*total_bsize+indx].push_back(Q_ptr[0]);
        //bins_Q[procid*total_bsize+indx].push_back(Q_ptr[1]);
        //bins_Q[procid*total_bsize+indx].push_back(Q_ptr[2]);
        //bins_f[procid*total_bsize+indx].push_back(i);
        // bins_Q[procid][indx].push_back(Q_ptr[0]);
        // bins_Q[procid][indx].push_back(Q_ptr[1]);
        // bins_Q[procid][indx].push_back(Q_ptr[2]);
        // bins_f[procid][indx].push_back(i);

				++Q_local;
				//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
				continue;
			}

		}

    //for(int proc = 0; proc < nprocs; ++proc){
		//  query_outside[proc].clear();
		//  f_index[proc].clear();
		//  for (int i = 0; i < total_bsize; ++i) {
    //    int bindx = trip[i].i_;
    //    //if(!bins_Q[proc][bindx].empty()){
    //    if(!bins_Q[proc*total_bsize+bindx].empty()){
    //      query_outside[proc].insert(query_outside[proc].end(),
    //          bins_Q[proc*total_bsize+bindx].begin(),
    //          bins_Q[proc*total_bsize+bindx].end());
    //      f_index[proc].insert(f_index[proc].end(),
    //          bins_f[proc*total_bsize+bindx].begin(),
    //          bins_f[proc*total_bsize+bindx].end());
    //    }

    //  }
    //}
    // pvfmm::aligned_delete<zTrip>(trip);
    // pvfmm::aligned_delete(bins_Q);
    // pvfmm::aligned_delete(bins_f);
    //pvfmm::aligned_delete<Real>(query_points);
		// Now sort the query points in zyx order
#ifdef SORT_QUERIES
		timings[3]+=-MPI_Wtime();
		zsort_queries(query_outside,f_index,N_reg,h,c_comm);
		timings[3]+=+MPI_Wtime();
		//if(procid==0) std::cout<<"Sorting time="<<s_time<<std::endl;;
		//if(procid==0) std::cout<<"Sorting Queries\n";
#endif

		// Now we need to send the query_points that land onto other processor's domain.
		// This is done using a sparse alltoallv.
		// Right now each process knows how much data to send to others, but does not know
		// how much data it should receive. This is a necessary information both for the MPI
		// command as well as memory allocation for received data.
		// So we first do an alltoall to get the f_index[proc].size from all processes.

#ifdef INTERP_DEBUG
  PCOUT << "Communicating sizes\n";
#endif
		for (int proc = 0; proc < nprocs; proc++) {
			if (!f_index[proc].empty()){
				f_index_procs_self_sizes[proc] = f_index[proc].size();
      }
			else
				f_index_procs_self_sizes[proc] = 0;
		}
		timings[0] += -MPI_Wtime();
		MPI_Alltoall(&f_index_procs_self_sizes[0], 1, MPI_INT,
				&f_index_procs_others_sizes[0], 1, MPI_INT, c_comm);
		timings[0] += +MPI_Wtime();

		for (int proc = 0; proc < nprocs; proc++) {
			if (f_index_procs_self_sizes[proc] > 0){
        procs_i_send_to_.push_back(proc);
      }
			if (f_index_procs_others_sizes[proc] > 0 ){
        procs_i_recv_from_.push_back(proc);
      }
		}
    if(!procs_i_recv_from_.empty())
      procs_i_recv_from_size_ = procs_i_recv_from_.size();
    else
      procs_i_recv_from_size_ = 0;
    if(!procs_i_send_to_.empty())
      procs_i_send_to_size_ = procs_i_send_to_.size();
    else
      procs_i_send_to_size_ = 0;

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
      if(proc >0){
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
      pvfmm::aligned_delete<Real>(this->all_query_points);
      pvfmm::aligned_delete<Real>(this->all_f_cubic);
		}
#ifdef INTERP_USE_MORE_MEM_L1
		all_query_points = pvfmm::aligned_new<Real>(
				(total_query_points+16)*(COORD_DIM+1)); // 16 for blocking in interp
#else
		all_query_points = pvfmm::aligned_new<Real>(
				(total_query_points+16)*(COORD_DIM)); // 16 for blocking in interp
#endif
		all_f_cubic = pvfmm::aligned_new<Real>(
				total_query_points * data_dof_max + (16*isize_g[2]*isize_g[1]+16*isize_g[1]+16));

#ifdef INTERP_DEBUG
    //parallel_print(total_query_points, "total_q_points");
    //parallel_print(procs_i_send_to_size_, "procs_i_send_to_size_");
    //parallel_print(procs_i_recv_from_size_, "procs_i_recv_from_size_");
  PCOUT << "communicating query points\n";
#endif
		// Now perform the allotall to send/recv query_points
		timings[0] += -MPI_Wtime();
		{
			int dst_r, dst_s;
			for (int i = 0; i < procs_i_recv_from_size_; ++i) {
				dst_r = procs_i_recv_from_[i];    //(procid+i)%nprocs;
				request[dst_r] = MPI_REQUEST_NULL;
				int roffset = f_index_procs_others_offset[dst_r] * COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
				MPI_Irecv(&all_query_points[roffset],
						f_index_procs_others_sizes[dst_r] * COORD_DIM,
						MPI_T, dst_r, 0, c_comm, &request[dst_r]);
			}
			for (int i = 0; i < procs_i_send_to_size_; ++i) {
				dst_s = procs_i_send_to_[i];    //(procid-i+nprocs)%nprocs;
				s_request[dst_s] = MPI_REQUEST_NULL;
				//int soffset = f_index_procs_self_offset[dst_s] * COORD_DIM;
				MPI_Isend(&query_outside[dst_s][0],
						f_index_procs_self_sizes[dst_s] * COORD_DIM, MPI_T,
						dst_s, 0, c_comm, &s_request[dst_s]);
			}
			for (int i = 0; i < procs_i_recv_from_size_; ++i) {
				int proc = procs_i_recv_from_[i];    //(procid+i)%nprocs;
				if (request[proc] != MPI_REQUEST_NULL)
					MPI_Wait(&request[proc], MPI_STATUS_IGNORE);
     }
			for (int i = 0; i < procs_i_send_to_size_; ++i) {
				int proc = procs_i_send_to_[i];    //(procid+i)%nprocs;
				if (s_request[proc] != MPI_REQUEST_NULL)
					MPI_Wait(&s_request[proc], MPI_STATUS_IGNORE);
     }
			//for (int i = 0; i < nprocs; ++i) {
			//	dst_r = i;    //(procid+i)%nprocs;
			//	dst_s = i;    //(procid-i+nprocs)%nprocs;
			//	s_request[dst_s] = MPI_REQUEST_NULL;
			//	request[dst_r] = MPI_REQUEST_NULL;
			//	int roffset = f_index_procs_others_offset[dst_r] * COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
			//	int soffset = f_index_procs_self_offset[dst_s] * COORD_DIM;
			//	if (f_index_procs_others_sizes[dst_r] != 0)
			//		MPI_Irecv(&all_query_points[roffset],
			//				f_index_procs_others_sizes[dst_r] * COORD_DIM,
			//				MPI_T, dst_r, 0, c_comm, &request[dst_r]);
			//	if (!query_outside[dst_s].empty())
			//		MPI_Isend(&query_outside[dst_s][0],
			//				f_index_procs_self_sizes[dst_s] * COORD_DIM, MPI_T,
			//				dst_s, 0, c_comm, &s_request[dst_s]);
			//}
			//// Wait for all the communication to finish
			//MPI_Status ierr;
			//for (int proc = 0; proc < nprocs; ++proc) {
			//	if (request[proc] != MPI_REQUEST_NULL)
			//		MPI_Wait(&request[proc], &ierr);
			//	if (s_request[proc] != MPI_REQUEST_NULL)
			//		MPI_Wait(&s_request[proc], &ierr);
			//}
		}
		timings[0] += +MPI_Wtime();
    time+=MPI_Wtime();
    //std::cout << "**** time = " << time << std::endl;

	}

#ifdef INTERP_DEBUG
  PCOUT << "done with comm\n";
#endif
  // ParLOG << "nplans_ = " << nplans_ << " data_dof_max = " << data_dof_max << std::endl;
  // ParLOG << "data_dofs[0] = " << data_dofs_[0] << " [1] = " << data_dofs_[1] << std::endl;
  for(int ver = 0; ver < nplans_; ++ver){
	for (int i = 0; i < nprocs; ++i) {
		MPI_Type_vector(data_dofs_[ver], f_index_procs_self_sizes[i], N_pts, MPI_T,
				&rtypes[i+ver*nprocs]);
		MPI_Type_vector(data_dofs_[ver], f_index_procs_others_sizes[i],
				total_query_points, MPI_T, &stypes[i+ver*nprocs]);
		MPI_Type_commit(&stypes[i+ver*nprocs]);
		MPI_Type_commit(&rtypes[i+ver*nprocs]);
	}
  }
#ifdef INTERP_USE_MORE_MEM_L1
  if(total_query_points !=0)
	rescale_xyzgrid(g_size, N_reg, N_reg_g, istart, isize, isize_g, total_query_points,
			all_query_points);
#else
  if(total_query_points !=0)
	rescale_xyz(g_size, N_reg, N_reg_g, istart, isize, isize_g, total_query_points,
			&all_query_points[0]);
#endif

    if (procs_i_recv_from_.size() != 0) procs_i_recv_from_.clear();
    if (procs_i_send_to_.size() != 0) procs_i_send_to_.clear();

	this->scatter_baked = true;
#ifdef INTERP_DEBUG
  PCOUT << "scatter DONE\n";
#endif
	return;
}


/*
 * Phase 2 of the parallel interpolation: This function must be called after the scatter function is called.
 * It performs local interpolation for all the points that the processor has for itself, as well as the interpolations
 * that it has to send to other processors. After the local interpolation is performed, a sparse
 * alltoall is performed so that all the interpolated results are sent/received.
 * version is a number between zero and nplans_ specifying which data_dof to use
 *
 */
void Interp3_Plan::interpolate(Real* __restrict ghost_reg_grid_vals,
		int*__restrict N_reg, int *__restrict isize, int*__restrict istart, const int N_pts, const int g_size,
		Real*__restrict query_values, int*__restrict c_dims, MPI_Comm c_comm, double *__restrict timings, int version) {
	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);
#ifdef INTERP_DEBUG
  PCOUT << "entered interpolate\n";
#endif
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
#ifdef FAST_INTERPV
  const int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];
  const int* N_reg_c = N_reg;
  const int* N_reg_g_c = N_reg_g;
  const int* istart_c = istart;
  const int* isize_g_c = isize_g;
  const int total_query_points_c = total_query_points;
  if(total_query_points!=0)
    for (int k = 0; k < data_dofs_[version]; ++k){
	    vectorized_interp3_ghost_xyz_p(&ghost_reg_grid_vals[k*N_reg3], 1, N_reg_c, N_reg_g_c, isize_g_c,
			  istart_c, total_query_points_c, g_size, &all_query_points[0], &all_f_cubic[k*total_query_points],
			  true);
      //std::cout << "data_dofs_[version ] = " << data_dofs_[version] << std::endl;
      //do{}while(1);
    }
#else
  const int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];
  if(total_query_points!=0)
    for (int k = 0; k < data_dofs_[version]; ++k)
	  optimized_interp3_ghost_xyz_p(&ghost_reg_grid_vals[k*N_reg3], 1, N_reg, N_reg_g, isize_g,
			istart, total_query_points, g_size, &all_query_points[0], &all_f_cubic[k*total_query_points],
			true);
#endif
#else
  if(total_query_points!=0)
	 interp3_ghost_xyz_p(ghost_reg_grid_vals, data_dofs_[version], N_reg, N_reg_g, isize_g,
			istart, total_query_points, g_size, &all_query_points[0], &all_f_cubic[0],
			true);
#endif
	timings[1] += +MPI_Wtime();

	// Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to
	// f_cubic_unordered.
#ifdef INTERP_DEBUG
  PCOUT << "finished interpolation, starting comm\n";
#endif
  double shuffle_time =0;
	timings[0] += -MPI_Wtime();
	{
		int dst_r, dst_s;

		//for (int i = 0; i < procs_i_send_to_size_; ++i) {
		for (int i = procs_i_send_to_size_-1; i >=0; --i) {
			//dst_r = (procid+i)%nprocs;
			dst_r = procs_i_send_to_[i];    //(procid-i+nprocs)%nprocs;
			request[dst_r] = MPI_REQUEST_NULL; //recv
			int roffset = f_index_procs_self_offset[dst_r];

			MPI_Irecv(&f_cubic_unordered[roffset], 1, rtypes[dst_r+version*nprocs], dst_r, 0,
					c_comm, &request[dst_r]);
		}
		for (int i = 0; i < procs_i_recv_from_size_; ++i) {
			// dst_s = (procid-i+nprocs)%nprocs;
			dst_s = procs_i_recv_from_[i];    //(procid+i)%nprocs;
			s_request[dst_s] = MPI_REQUEST_NULL; //send
			int soffset = f_index_procs_others_offset[dst_s];

			MPI_Isend(&all_f_cubic[soffset], 1, stypes[dst_s+version*nprocs], dst_s, 0, c_comm,
					&s_request[dst_s]);
		}
    // wait to receive your part
			for (int i = 0; i < procs_i_send_to_size_; ++i) {
				int proc = procs_i_send_to_[i];    //(procid+i)%nprocs;
				if (request[proc] != MPI_REQUEST_NULL)
					MPI_Wait(&request[proc], MPI_STATUS_IGNORE);
          shuffle_time += -MPI_Wtime();
	        for (int dof = 0; dof < data_dofs_[version]; ++dof) {
            Real* ptr = &f_cubic_unordered[f_index_procs_self_offset[proc]+dof*N_pts];
#pragma omp parallel for
                for (int i = 0; i < (int)f_index[proc].size(); ++i) {
                  int ind = f_index[proc][i];
                  query_values[ind + dof * N_pts] =ptr[i];
                }
          }
          shuffle_time += +MPI_Wtime();
     }
   // wait for send
			for (int i = 0; i < procs_i_recv_from_size_; ++i) {
				int proc = procs_i_recv_from_[i];    //(procid+i)%nprocs;
				if (s_request[proc] != MPI_REQUEST_NULL)
					MPI_Wait(&s_request[proc], MPI_STATUS_IGNORE);
     }

		// for (int i = 0; i < nprocs; ++i) {
		//	dst_r = (procid+i)%nprocs;
		//	dst_s = (procid-i+nprocs)%nprocs;
		//	// s_request[dst_s] = MPI_REQUEST_NULL;
		//	request[2*dst_s] = MPI_REQUEST_NULL; //send
		//	request[2*dst_r+1] = MPI_REQUEST_NULL; //recv
		//	// Notice that this is the adjoint of the first comm part
		//	// because now you are sending others f and receiving your part of f
		//	int soffset = f_index_procs_others_offset[dst_s];
		//	int roffset = f_index_procs_self_offset[dst_r];

		//	if (f_index_procs_self_sizes[dst_r] != 0)
		//		MPI_Irecv(&f_cubic_unordered[roffset], 1, rtype[dst_r], dst_r, 0,
		//				c_comm, &request[2*dst_r+1]);
		//	if (f_index_procs_others_sizes[dst_s] != 0)
		//		MPI_Isend(&all_f_cubic[soffset], 1, stype[dst_s], dst_s, 0, c_comm,
		//				&request[2*dst_s]);
		// }
    // MPI_Waitall(2*nprocs, request, MPI_STATUSES_IGNORE);
	}
	timings[0] += +MPI_Wtime();
  timings[0] -= shuffle_time;

	// Now copy back f_cubic_unordered to f_cubic in the correct f_index
	//for (int dof = 0; dof < data_dof; ++dof) {
	//	for (int proc = 0; proc < nprocs; ++proc) {
	//		if (!f_index[proc].empty())
	//			for (int i = 0; i < f_index[proc].size(); ++i) {
	//				int ind = f_index[proc][i];
	//				//f_cubic[ind]=all_f_cubic[f_index_procs_others_offset[proc]+i];
	//				query_values[ind + dof * N_pts] =
	//						f_cubic_unordered[f_index_procs_self_offset[proc]
	//								+ i + dof * N_pts];
	//			}
	//	}
	//}

#ifdef INTERP_DEBUG
  PCOUT << "interpolation done\n";
#endif
	return;
}

Interp3_Plan::~Interp3_Plan() {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (this->allocate_baked) {
    //if(this->scatter_baked == false)
      pvfmm::aligned_delete<Real>(query_points);

    pvfmm::aligned_delete<int>(f_index_procs_others_offset);
		pvfmm::aligned_delete<int>(f_index_procs_self_offset);
		pvfmm::aligned_delete<int>(f_index_procs_self_sizes);
		pvfmm::aligned_delete<int>(f_index_procs_others_sizes);

    pvfmm::aligned_delete<MPI_Request>(s_request);
		pvfmm::aligned_delete<MPI_Request>(request);
		//vectors
		for (int proc = 0; proc < nprocs; ++proc) {
			std::vector<int>().swap(f_index[proc]);
			std::vector<Real>().swap(query_outside[proc]);
		}
    pvfmm::aligned_delete<Real>(f_cubic_unordered);

	}

	if (this->scatter_baked) {
		for (int ver = 0; ver < nplans_; ++ver)
		for (int i = 0; i < nprocs; ++i) {
			MPI_Type_free(&stypes[i+ver*nprocs]);
			MPI_Type_free(&rtypes[i+ver*nprocs]);
		}
    pvfmm::aligned_delete<Real>(all_query_points);
    pvfmm::aligned_delete<Real>(all_f_cubic);
	}

	if (this->allocate_baked) {
    pvfmm::aligned_delete<MPI_Datatype>(rtypes);
    pvfmm::aligned_delete<MPI_Datatype>(stypes);
    pvfmm::aligned_delete<int>(data_dofs_);
	}

	return;
}

//void Interp3_Plan::high_order_interpolate(Real* ghost_reg_grid_vals, int data_dof,
//		int* N_reg, int * isize, int* istart, const int N_pts, const int g_size,
//		Real* query_values, int* c_dims, MPI_Comm c_comm, double * timings, int interp_order) {
//	int nprocs, procid;
//	MPI_Comm_rank(c_comm, &procid);
//	MPI_Comm_size(c_comm, &nprocs);
//	if (this->allocate_baked == false) {
//		std::cout
//				<< "ERROR Interp3_Plan interpolate called before calling allocate.\n";
//		return;
//	}
//	if (this->scatter_baked == false) {
//		std::cout
//				<< "ERROR Interp3_Plan interpolate called before calling scatter.\n";
//		return;
//	}
//
//	timings[1] += -MPI_Wtime();
//	interp3_ghost_xyz_p(ghost_reg_grid_vals, data_dof, N_reg, N_reg_g, isize_g,
//			istart, total_query_points, g_size, all_query_points, all_f_cubic, interp_order,
//			true);
//	timings[1] += +MPI_Wtime();
//
//	// Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to
//	// f_cubic_unordered.
//	timings[0] += -MPI_Wtime();
//	{
//		int dst_r, dst_s;
//		for (int i = 0; i < nprocs; ++i) {
//			dst_r = i;  //(procid+i)%nprocs;
//			dst_s = i;  //(procid-i+nprocs)%nprocs;
//			s_request[dst_s] = MPI_REQUEST_NULL;
//			request[dst_r] = MPI_REQUEST_NULL;
//			// Notice that this is the adjoint of the first comm part
//			// because now you are sending others f and receiving your part of f
//			int soffset = f_index_procs_others_offset[dst_r];
//			int roffset = f_index_procs_self_offset[dst_s];
//
//			if (f_index_procs_self_sizes[dst_r] != 0)
//				MPI_Irecv(&f_cubic_unordered[roffset], 1, rtype[i], dst_r, 0,
//						c_comm, &request[dst_r]);
//			if (f_index_procs_others_sizes[dst_s] != 0)
//				MPI_Isend(&all_f_cubic[soffset], 1, stype[i], dst_s, 0, c_comm,
//						&s_request[dst_s]);
//		}
//		MPI_Status ierr;
//		for (int proc = 0; proc < nprocs; ++proc) {
//			if (request[proc] != MPI_REQUEST_NULL)
//				MPI_Wait(&request[proc], &ierr);
//			if (s_request[proc] != MPI_REQUEST_NULL)
//				MPI_Wait(&s_request[proc], &ierr);
//		}
//	}
//	timings[0] += +MPI_Wtime();
//
//	// Now copy back f_cubic_unordered to f_cubic in the correct f_index
//	for (int dof = 0; dof < data_dof; ++dof) {
//		for (int proc = 0; proc < nprocs; ++proc) {
//			if (!f_index[proc].empty())
//				for (int i = 0; i < f_index[proc].size(); ++i) {
//					int ind = f_index[proc][i];
//					//f_cubic[ind]=all_f_cubic[f_index_procs_others_offset[proc]+i];
//					query_values[ind + dof * N_pts] =
//							f_cubic_unordered[f_index_procs_self_offset[proc]
//									+ i + dof * N_pts];
//				}
//		}
//	}
//
//	return;
//}
//
//
//
//
void Interp3_Plan::scatter(int* N_reg, int * isize, int* istart,
		const int N_pts, const int g_size, Real* query_points_in, int* c_dims,
		MPI_Comm c_comm, double * timings) {
#ifdef FAST_INTERP
  return fast_scatter(N_reg, isize, istart, N_pts,
      g_size, query_points_in, c_dims, c_comm, timings);
#else
  std::cout << "SCATTER CALLED!\n" << std::endl;
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
      pvfmm::aligned_delete<Real>(this->all_query_points);
      pvfmm::aligned_delete<Real>(this->all_f_cubic);
			all_query_points = pvfmm::aligned_new<Real>(
					all_query_points_allocation);
			all_f_cubic = pvfmm::aligned_new<Real>(
					total_query_points * data_dof_max);
		} else {
			all_query_points = pvfmm::aligned_new<Real>(
					all_query_points_allocation);
			all_f_cubic = pvfmm::aligned_new<Real>(
					total_query_points * data_dof_max);
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
	}

  for(int ver = 0; ver < nplans_; ++ver){
	for (int i = 0; i < nprocs; ++i) {
		MPI_Type_vector(data_dofs_[ver], f_index_procs_self_sizes[i], N_pts, MPI_T,
				&rtypes[i+ver*nprocs]);
		MPI_Type_vector(data_dofs_[ver], f_index_procs_others_sizes[i],
				total_query_points, MPI_T, &stypes[i+ver*nprocs]);
		MPI_Type_commit(&stypes[i+ver*nprocs]);
		MPI_Type_commit(&rtypes[i+ver*nprocs]);
	}
  }

	rescale_xyz(g_size, N_reg, N_reg_g, istart, isize, isize_g, total_query_points,
			all_query_points);
	//rescale_xyz(g_size, N_reg, N_reg_g, istart, isize, total_query_points,
	//		all_query_points);
    if(!procs_i_recv_from_.empty()) procs_i_recv_from_.clear();
    if(!procs_i_send_to_.empty()) procs_i_send_to_.clear();
	this->scatter_baked = true;
	return;
#endif
} // end scatter
