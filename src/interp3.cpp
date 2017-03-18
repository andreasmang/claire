// This function performs a 3D cubic interpolation.

#include <cmath>
#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <vector>

#include <accfft.h>
#include <interp3.hpp>
#define COORD_DIM 3
//#define VERBOSE2
#define sleep(x) ;


void optimized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
	//for (int i = 0; i < 4; i++) {
	//	lagr_denom[i] = 1;
	//	for (int j = 0; j < 4; j++) {
	//		if (i != j)
	//			lagr_denom[i] /= (Real) (i - j);
	//	}
	// }
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;
  //std::cout << " lagr_denom[0] = " << lagr_denom[0]
  //          << " lagr_denom[1] = " << lagr_denom[1]
  //          << " lagr_denom[2] = " << lagr_denom[2]
  //          << " lagr_denom[3] = " << lagr_denom[3] << std::endl;
  //do{}while(1);
	int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];

	for (int i = 0; i < N_pts; i++) {
    {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
  }
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg_g[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
    if(grid_indx[0]> isize_g[0]-3 || grid_indx[1]> isize_g[1]-3 ||grid_indx[2]> isize_g[2] -3)
    {
//#ifdef VERBOSE2
		std::cout<<"***** query point="<<query_points[0]<<" "<<query_points[1]<<" "<<query_points[2]<<std::endl;
		std::cout<<"***** grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;
		std::cout<<"***** point="<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
		std::cout<<"f @grid_index="<<reg_grid_vals[grid_indx[0]*isize_g[1]*isize_g[2]+grid_indx[1]*isize_g[2]+grid_indx[2]]<<std::endl;
		std::cout<<"hp= "<<1./N_reg_g[0]<<std::endl;
		std::cout<<"N_reg_g= "<<N_reg_g[0]<<" "<<N_reg_g[1]<<" "<<N_reg_g[2]<<std::endl;
//#endif
  }
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
      int indx = 0;
			for (int j0 = 0; j0 < 4; j0++) {
				for (int j1 = 0; j1 < 4; j1++) {
            Real M0M1 = M[0][j0]*M[1][j1];
            Real val_ = 0;
			    for (int j2 = 0; j2 < 4; j2++) {
						//int indx = j2
						//		+ isize_g[2] * j1
						//		+ isize_g[2] * isize_g[1] *j0;
						val_ += M[2][j2]
								* reg_grid_vals[indx + indxx + k*N_reg3];
            ++indx;
					}
          val += val_ * M0M1;
          indx += isize_g[2]-4;
				}
        indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];
			}
			query_values[i + k * N_pts] = val;
		}
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p


/*
 * Performs a parallel 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * This function assumes that the input grid values have been padded on all sides
 * by g_size grids.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg: The size of the original grid in each dimension.
 * @param[in] isize: The locally owned sizes that each process owns
 * @param[in] istart: The start index of each process in the global array
 * @param[in] N_pts The number of query points
 * @param[in] g_size The number of ghost points padded around the input array
 *
 * @param[in] query_points_in The coordinates of the query points where the interpolated values are sought
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 * @param[in] c_dims: Size of the cartesian MPI communicator
 * @param[in] c_comm: MPI Communicator
 *
 */

void par_interp3_ghost_xyz_p(Real* ghost_reg_grid_vals, int data_dof,
		int* N_reg, int * isize, int* istart, const int N_pts, const int g_size,
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
} // end of par_interp3_ghost_xyz_p

/*
 * Rescales the query points to [0,1) range for the parallel case. Note that the input query_points are initially
 * in the global range, however, each parallel process needs to rescale it to [0,1) for its local interpolation.
 * Since interpolation is scale invariant, this would not affect the interpolation result.
 * This function assumes that the ghost padding was done only in x, y, and z directions.
 *
 * @param[in] g_size: The ghost padding size
 * @param[in] N_reg: The original (unpadded) global size
 * @param[in] N_reg_g: The padded global size
 * @param[in] istart: The original isize for each process
 * @param[in] N_pts: The number of query points
 * @param[in,out]: query_points: The query points coordinates
 *
 */
void rescale_xyz(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		const int N_pts, Real* query_points) {

	if (g_size == 0)
		return;
	Real hp[3];
	Real h[3];
	hp[0] = 1. / N_reg_g[0]; // New mesh size
	hp[1] = 1. / N_reg_g[1]; // New mesh size
	hp[2] = 1. / N_reg_g[2]; // New mesh size

	h[0] = 1. / (N_reg[0]); // old mesh size
	h[1] = 1. / (N_reg[1]); // old mesh size
	h[2] = 1. / (N_reg[2]); // old mesh size

	Real factor[3];
	factor[0] = (1. - (2. * g_size + 1.) * hp[0]) / (1. - h[0]);
	factor[1] = (1. - (2. * g_size + 1.) * hp[1]) / (1. - h[1]);
	factor[2] = (1. - (2. * g_size + 1.) * hp[2]) / (1. - h[2]);
	for (int i = 0; i < N_pts; i++) {
		query_points[0 + COORD_DIM * i] = (query_points[0 + COORD_DIM * i]
				- istart[0] * h[0]) * factor[0] + g_size * hp[0];
		query_points[1 + COORD_DIM * i] = (query_points[1 + COORD_DIM * i]
				- istart[1] * h[1]) * factor[1] + g_size * hp[1];
		query_points[2 + COORD_DIM * i] = (query_points[2 + COORD_DIM * i]
				- istart[2] * h[2]) * factor[2] + g_size * hp[2];
	}
	return;
} // end of rescale_xyz

/*
 * Rescales the query points to [0,1) range for the parallel case. Note that the input query_points are initially
 * in the global range, however, each parallel process needs to rescale it to [0,1) for its local interpolation.
 * Since interpolation is scale invariant, this would not affect the interpolation result.
 * This function assumes that the ghost padding was done only in x, and y directions and the z direction is not padded.
 *
 * @param[in] g_size: The ghost padding size
 * @param[in] N_reg: The original (unpadded) global size
 * @param[in] N_reg_g: The padded global size
 * @param[in] istart: The original isize for each process
 * @param[in] N_pts: The number of query points
 * @param[in,out]: query_points: The query points coordinates
 *
 */
void rescale(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		const int N_pts, Real* query_points) {

	if (g_size == 0)
		return;
	Real hp[3];
	Real h[3];
	hp[0] = 1. / N_reg_g[0]; // New mesh size
	hp[1] = 1. / N_reg_g[1]; // New mesh size
	hp[2] = 1. / N_reg_g[2]; // New mesh size

	h[0] = 1. / (N_reg[0]); // old mesh size
	h[1] = 1. / (N_reg[1]); // old mesh size
	h[2] = 1. / (N_reg[2]); // old mesh size

	Real factor[3];
	factor[0] = (1. - (2. * g_size + 1.) * hp[0]) / (1. - h[0]);
	factor[1] = (1. - (2. * g_size + 1.) * hp[1]) / (1. - h[1]);
	factor[2] = (1. - (2. * g_size + 1.) * hp[2]) / (1. - h[2]);
	for (int i = 0; i < N_pts; i++) {
		query_points[0 + COORD_DIM * i] = (query_points[0 + COORD_DIM * i]
				- istart[0] * h[0]) * factor[0] + g_size * hp[0];
		query_points[1 + COORD_DIM * i] = (query_points[1 + COORD_DIM * i]
				- istart[1] * h[1]) * factor[1] + g_size * hp[1];
		//query_points[2+COORD_DIM*i]=(query_points[2+COORD_DIM*i]-istart[2]*h[2])*factor[2]+g_size*hp[2];
	}
	return;
} // end of rescale

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * This function assumes that the input grid values have been padded on all sides
 * by g_size grids.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg An integer pointer that specifies the size of the grid in each dimension.
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] g_size The number of ghost points padded around the input array
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 * snafu
 *
 */

void interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
	//for (int i = 0; i < 4; i++) {
	//	lagr_denom[i] = 1;
	//	for (int j = 0; j < 4; j++) {
	//		if (i != j)
	//			lagr_denom[i] /= (Real) (i - j);
	//	}
	// }
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;
  //std::cout << " lagr_denom[0] = " << lagr_denom[0]
  //          << " lagr_denom[1] = " << lagr_denom[1]
  //          << " lagr_denom[2] = " << lagr_denom[2]
  //          << " lagr_denom[3] = " << lagr_denom[3] << std::endl;
  //do{}while(1);
	int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];

	for (int i = 0; i < N_pts; i++) {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg_g[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg_g[j];
		}

#ifdef VERBOSE2
		std::cout<<"***** grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;
		std::cout<<"***** point="<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
		std::cout<<"f @grid_index="<<reg_grid_vals[grid_indx[0]*isize_g[1]*isize_g[2]+grid_indx[1]*isize_g[2]+grid_indx[2]]<<std::endl;
		std::cout<<"hp= "<<1./N_reg_g[0]<<std::endl;
		std::cout<<"N_reg_g= "<<N_reg_g[0]<<" "<<N_reg_g[1]<<" "<<N_reg_g[2]<<std::endl;
#endif

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						int indx = ((grid_indx[2] + j2) % isize_g[2])
								+ isize_g[2]
										* ((grid_indx[1] + j1) % isize_g[1])
								+ isize_g[2] * isize_g[1]
										* ((grid_indx[0] + j0) % isize_g[0]);
						//val += M[0][j0] * M[1][j1] * M[2][j2]
						//		* reg_grid_vals[0];
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			//query_values[0] = val;
			query_values[i + k * N_pts] = val;
		}
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * This function assumes that the input grid values have been padded on all sides
 * by g_size grids.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg An integer pointer that specifies the size of the grid in each dimension.
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] g_size The number of ghost points padded around the input array
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[in] interp_order The order of interpolation (e.g. 3 for cubic)
 * @param[out] query_values The interpolated values
 *
 */

void interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
    int interp_order,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[interp_order + 1];
	for (int i = 0; i < interp_order + 1; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < interp_order + 1; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];

	for (int i = 0; i < N_pts; i++) {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg_g[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg_g[j];
		}

#ifdef VERBOSE2
		std::cout<<"***** grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;
		std::cout<<"***** point="<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
		std::cout<<"f @grid_index="<<reg_grid_vals[grid_indx[0]*isize_g[1]*isize_g[2]+grid_indx[1]*isize_g[2]+grid_indx[2]]<<std::endl;
		std::cout<<"hp= "<<1./N_reg_g[0]<<std::endl;
		std::cout<<"N_reg_g= "<<N_reg_g[0]<<" "<<N_reg_g[1]<<" "<<N_reg_g[2]<<std::endl;
#endif

		Real M[3][interp_order + 1];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < interp_order + 1; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < interp_order + 1; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < interp_order + 1; j2++) {
				for (int j1 = 0; j1 < interp_order + 1; j1++) {
					for (int j0 = 0; j0 < interp_order + 1; j0++) {
						int indx = ((grid_indx[2] + j2) % isize_g[2])
								+ isize_g[2]
										* ((grid_indx[1] + j1) % isize_g[1])
								+ isize_g[2] * isize_g[1]
										* ((grid_indx[0] + j0) % isize_g[0]);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * This function assumes that the input grid values have been padded on x, and y directions
 * by g_size grids.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg An integer pointer that specifies the size of the grid in each dimension.
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] g_size The number of ghost points padded around the input array
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 */

void interp3_ghost_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values) {

	// First we need to rescale the query points to the new padded dimensions
	// To avoid changing the user's input we first copy the query points to a
	// new array
	Real* query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
	memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
	rescale(g_size, N_reg, N_reg_g, istart, N_pts, query_points);

	//std::cout<<"N_reg[0]="<<N_reg[0]<<" N_reg[1]="<<N_reg[1]<<" N_reg[2]="<<N_reg[2]<<std::endl;
	//std::cout<<"N_reg_g[0]="<<N_reg_g[0]<<" N_reg_g[1]="<<N_reg_g[1]<<" N_reg_g[2]="<<N_reg_g[2]<<std::endl;

	Real lagr_denom[4];
	for (int i = 0; i < 4; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < 4; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];
	//int N_pts=query_points.size()/COORD_DIM;

	for (int i = 0; i < N_pts; i++) {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];
		//grid_indx[0]=15;
		//grid_indx[1]=15;
		//grid_indx[2]=14;
		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg_g[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg_g[j];
		}
#ifdef VERBOSE2
		std::cout<<"***** grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;
		std::cout<<"***** point="<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
		std::cout<<"f @grid_index="<<reg_grid_vals[grid_indx[0]*isize_g[1]*isize_g[2]+grid_indx[1]*isize_g[2]+grid_indx[2]]<<std::endl;
		std::cout<<"hp= "<<1./N_reg_g[0]<<std::endl;
		std::cout<<"N_reg_g= "<<N_reg_g[0]<<" "<<N_reg_g[1]<<" "<<N_reg_g[2]<<std::endl;
#endif

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						//int indx = ((grid_indx[2]+j2)%N_reg) + N_reg*((grid_indx[1]+j1)%N_reg) + N_reg*N_reg*((grid_indx[0]+j0)%N_reg);
						//int indx = ((grid_indx[2]+j2)%N_reg_g[2]) + N_reg_g[2]*((grid_indx[1]+j1)%N_reg_g[1]) + N_reg_g[2]*N_reg_g[1]*((grid_indx[0]+j0)%N_reg_g[0]);
						int indx = ((grid_indx[2] + j2) % isize_g[2])
								+ isize_g[2]
										* ((grid_indx[1] + j1) % isize_g[1])
								+ isize_g[2] * isize_g[1]
										* ((grid_indx[0] + j0) % isize_g[0]);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}
	free(query_points);
}            //end of interp3_ghost_p

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg An integer pointer that specifies the size of the grid in each dimension.
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 */

void interp3_p(Real* reg_grid_vals, int data_dof, int* N_reg, const int N_pts,
		Real* query_points, Real* query_values) {

	Real lagr_denom[4];
	for (int i = 0; i < 4; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < 4; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = N_reg[0] * N_reg[1] * N_reg[2];
	//int N_pts=query_points.size()/COORD_DIM;
	//query_values.resize(N_pts*data_dof);

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];
		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg[j];
		}
		//std::cout<<"grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						//int indx = ((grid_indx[2]+j2)%N_reg) + N_reg*((grid_indx[1]+j1)%N_reg) + N_reg*N_reg*((grid_indx[0]+j0)%N_reg);
						int indx = ((grid_indx[2] + j2) % N_reg[2])
								+ N_reg[2] * ((grid_indx[1] + j1) % N_reg[1])
								+ N_reg[2] * N_reg[1]
										* ((grid_indx[0] + j0) % N_reg[0]);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}
} // end of interp3_p

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) ).
 * Limitation: The grid must be cubic, i.e. the number of grid points must be the same
 * in all dimensions.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg The size of the regular grid (The grid must have the same size in all dimensions)
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 */

void interp3_p(Real* reg_grid_vals, int data_dof, int N_reg, const int N_pts,
		Real* query_points, Real* query_values) {

	Real lagr_denom[4];
	for (int i = 0; i < 4; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < 4; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = N_reg * N_reg * N_reg;
	//int N_pts=query_points.size()/COORD_DIM;
	//query_values.resize(N_pts*data_dof);

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];
		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg;
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg;
		}
		//std::cout<<"grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						//int indx = ((grid_indx[0]+j0)%N_reg) + N_reg*((grid_indx[1]+j1)%N_reg) + N_reg*N_reg*((grid_indx[2]+j2)%N_reg);
						int indx = ((grid_indx[2] + j2) % N_reg)
								+ N_reg * ((grid_indx[1] + j1) % N_reg)
								+ N_reg * N_reg * ((grid_indx[0] + j0) % N_reg);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}
} // end of interp3_p

/*
 * Performs a 3D cubic interpolation for a column major periodic input (x \in [0,1) )
 * Limitation: The grid must be cubic, i.e. the number of grid points must be the same
 * in all dimensions.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg The size of the regular grid (The grid must have the same size in all dimensions)
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 */
void interp3_p_col(Real* reg_grid_vals, int data_dof, int N_reg,
		const int N_pts, Real* query_points, Real* query_values) {

	Real lagr_denom[4];
	for (int i = 0; i < 4; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < 4; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = N_reg * N_reg * N_reg;

	Real point[COORD_DIM];
	int grid_indx[COORD_DIM];
	for (int i = 0; i < N_pts; i++) {
		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg;
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg;
		}

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						int indx = ((grid_indx[0] + j0) % N_reg)
								+ N_reg * ((grid_indx[1] + j1) % N_reg)
								+ N_reg * N_reg * ((grid_indx[2] + j2) % N_reg);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}
} // end of interp3_p_col

