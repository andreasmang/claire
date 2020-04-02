// debugging with ib 4 bin/step1 16 16 16

#include <stdlib.h>
#include <math.h> // M_PI
#include <mpi.h>
//#include <accfft.h>
#include <interp3.hpp>
#include <CLAIREUtils.hpp>
//#define VERBOSE2
#include "RegOpt.hpp"

/*
 * Get the left right ghost cells.
 *
 * @param[out] padded_data: The output of the function which pads the input array data, with the ghost
 * cells from left and right neighboring processors.
 * @param[in] data: Input data to be padded
 * @param[in] g_size: The size of the ghost cell padding. Note that it cannot exceed the neighboring processor's
 * local data size
 * @param[in] plan: AccFFT R2C plan
 */
void ghost_left_right(pvfmm::Iterator<Real> padded_data, Real* data, int g_size,
		reg::RegOpt* m_Opt) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//	MPI_Comm c_comm = plan->c_comm;

	/* Get the local pencil size and the allocation size */
	//int isize[3],osize[3],istart[3],ostart[3];
	IntType * isize = m_Opt->m_Domain.isize;
//	int * osize = plan->isize;
//	int * istart = plan->istart;
//	int * ostart = plan->ostart;
//	int alloc_max = plan->alloc_max;

	MPI_Comm row_comm = m_Opt->m_Domain.rowcomm;
	int nprocs_r, procid_r;
	MPI_Comm_rank(row_comm, &procid_r);
	MPI_Comm_size(row_comm, &nprocs_r);
	/* Halo Exchange along y axis
	 * Phase 1: Write local data to be sent to the right process to RS
	 */
#ifdef VERBOSE2
	PCOUT<<"\nGL Row Communication\n";
#endif

	int rs_buf_size = g_size * isize[2] * isize[0];
	Real *RS = (Real*) accfft_alloc(rs_buf_size * sizeof(Real)); // Stores local right ghost data to be sent
	Real *GL = (Real*) accfft_alloc(rs_buf_size * sizeof(Real)); // Left Ghost cells to be received

	for (int x = 0; x < isize[0]; ++x)
		memcpy(&RS[x * g_size * isize[2]],
				&data[x * isize[2] * isize[1] + (isize[1] - g_size) * isize[2]],
				g_size * isize[2] * sizeof(Real));

	/* Phase 2: Send your data to your right process
	 * First question is who is your right process?
	 */
	int dst_s = (procid_r + 1) % nprocs_r;
	int dst_r = (procid_r - 1) % nprocs_r;
	if (procid_r == 0)
		dst_r = nprocs_r - 1;
	MPI_Request rs_s_request, rs_r_request;
	MPI_Status ierr;
	MPI_Isend(RS, rs_buf_size, MPI_T, dst_s, 0, row_comm, &rs_s_request);
	MPI_Irecv(GL, rs_buf_size, MPI_T, dst_r, 0, row_comm, &rs_r_request);
	MPI_Wait(&rs_s_request, &ierr);
	MPI_Wait(&rs_r_request, &ierr);

#ifdef VERBOSE2
	if(procid==0) {
		std::cout<<"procid="<<procid<<" data\n";
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<isize[1];++j)
			std::cout<<data[(i*isize[1]+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<g_size;++j)
			std::cout<<RS[(i*g_size+j)*isize[2]]<<" ";
			//PCOUT<<data[(i*isize[1]+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}
	sleep(1);
	if(procid==1) {
		std::cout<<"procid="<<procid<<" GL data\n";
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<g_size;++j)
			std::cout<<GL[(i*g_size+j)*isize[2]]<<" ";
			//PCOUT<<data[(i*isize[1]+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}
	PCOUT<<"\nGR Row Communication\n";
#endif

	/* Phase 3: Now do the exact same thing for the right ghost side */
	int ls_buf_size = g_size * isize[2] * isize[0];
	Real *LS = (Real*) accfft_alloc(ls_buf_size * sizeof(Real)); // Stores local right ghost data to be sent
	Real *GR = (Real*) accfft_alloc(ls_buf_size * sizeof(Real)); // Left Ghost cells to be received
	for (int x = 0; x < isize[0]; ++x)
		memcpy(&LS[x * g_size * isize[2]], &data[x * isize[2] * isize[1]],
				g_size * isize[2] * sizeof(Real));

	/* Phase 4: Send your data to your right process
	 * First question is who is your right process?
	 */
	dst_s = (procid_r - 1) % nprocs_r;
	dst_r = (procid_r + 1) % nprocs_r;
	if (procid_r == 0)
		dst_s = nprocs_r - 1;
	MPI_Isend(LS, ls_buf_size, MPI_T, dst_s, 0, row_comm, &rs_s_request);
	MPI_Irecv(GR, ls_buf_size, MPI_T, dst_r, 0, row_comm, &rs_r_request);
	MPI_Wait(&rs_s_request, &ierr);
	MPI_Wait(&rs_r_request, &ierr);

#ifdef VERBOSE2
	if(procid==1) {
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<isize[1];++j)
			std::cout<<data[(i*isize[1]+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
		std::cout<<"1 sending to dst_s="<<dst_s<<std::endl;
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<g_size;++j)
			std::cout<<LS[(i*g_size+j)*isize[2]]<<" ";
			//PCOUT<<data[(i*isize[1]+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
		std::cout<<"\n";
	}
	sleep(1);
	if(procid==0) {
		std::cout<<"0 receiving from dst_r="<<dst_r<<std::endl;
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<g_size;++j)
			std::cout<<GR[(i*g_size+j)*isize[2]]<<" ";
			//PCOUT<<data[(i*isize[1]+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}
#endif

	// Phase 5: Pack the data GL+ data + GR
	for (int i = 0; i < isize[0]; ++i) {
		memcpy(&padded_data[i * isize[2] * (isize[1] + 2 * g_size)],
				&GL[i * g_size * isize[2]], g_size * isize[2] * sizeof(Real));
		memcpy(
				&padded_data[i * isize[2] * (isize[1] + 2 * g_size)
						+ g_size * isize[2]], &data[i * isize[2] * isize[1]],
				isize[1] * isize[2] * sizeof(Real));
		memcpy(
				&padded_data[i * isize[2] * (isize[1] + 2 * g_size)
						+ g_size * isize[2] + isize[2] * isize[1]],
				&GR[i * g_size * isize[2]], g_size * isize[2] * sizeof(Real));
	}

#ifdef VERBOSE2
	if(procid==0) {
		std::cout<<"procid= "<<procid<<" padded_array=\n";
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<isize[1]+2*g_size;++j)
			std::cout<<padded_data[(i*(isize[1]+2*g_size)+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}
#endif

	accfft_free(LS);
	accfft_free(GR);
	accfft_free(RS);
	accfft_free(GL);

}

/*
 * Get the top bottom ghost cells AFTER getting the left right ones.
 *
 * @param[out] ghost_data: The output of the function which pads the input array data, with the ghost
 * cells from neighboring processors.
 * @param[in] padded_data: Input data that is already padded with ghost cells from left and right.
 * @param[in] g_size: The size of the ghost cell padding. Note that it cannot exceed the neighboring processor's
 * local data size
 * @param[in] plan: AccFFT R2C plan
 */
void ghost_top_bottom(pvfmm::Iterator<Real> ghost_data, pvfmm::Iterator<Real> padded_data, int g_size,
		reg::RegOpt* m_Opt) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

//	MPI_Comm c_comm = plan->c_comm;

	/* Get the local pencil size and the allocation size */
	//int isize[3],osize[3],istart[3],ostart[3];
	IntType * isize = m_Opt->m_Domain.isize;
//	int * osize = plan->isize;
//	int * istart = plan->istart;
//	int * ostart = plan->ostart;
//	int alloc_max = plan->alloc_max;

	MPI_Comm col_comm = m_Opt->m_Domain.colcomm;
	int nprocs_c, procid_c;
	MPI_Comm_rank(col_comm, &procid_c);
	MPI_Comm_size(col_comm, &nprocs_c);

	/* Halo Exchange along x axis
	 * Phase 1: Write local data to be sent to the bottom process
	 */
#ifdef VERBOSE2
	PCOUT<<"\nGB Col Communication\n";
#endif
	int bs_buf_size = g_size * isize[2] * (isize[1] + 2 * g_size); // isize[1] now includes two side ghost cells
	//Real *BS=(Real*)accfft_alloc(bs_buf_size*sizeof(Real)); // Stores local right ghost data to be sent
  pvfmm::Iterator<Real> GT = pvfmm::aligned_new<Real>(bs_buf_size); // Left Ghost cells to be received
	// snafu: not really necessary to do memcpy, you can simply use padded_data directly
	//memcpy(BS,&padded_data[(isize[0]-g_size)*isize[2]*(isize[1]+2*g_size)],bs_buf_size*sizeof(Real));
  Real* BS = &padded_data[(isize[0] - g_size) * isize[2]
                            * (isize[1] + 2 * g_size)];
	/* Phase 2: Send your data to your bottom process
	 * First question is who is your bottom process?
	 */
	int dst_s = (procid_c + 1) % nprocs_c;
	int dst_r = (procid_c - 1) % nprocs_c;
	if (procid_c == 0)
		dst_r = nprocs_c - 1;
	MPI_Request bs_s_request, bs_r_request;
	MPI_Status ierr;
	MPI_Isend(&BS[0], bs_buf_size, MPI_T, dst_s, 0, col_comm, &bs_s_request);
	MPI_Irecv(&GT[0], bs_buf_size, MPI_T, dst_r, 0, col_comm, &bs_r_request);
	MPI_Wait(&bs_s_request, &ierr);
	MPI_Wait(&bs_r_request, &ierr);

#ifdef VERBOSE2
	if(procid==0) {
		std::cout<<"procid= "<<procid<<" padded_array=\n";
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<isize[1]+2*g_size;++j)
			std::cout<<padded_data[(i*(isize[1]+2*g_size)+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}

		std::cout<<"procid= "<<procid<<" dst_s="<<dst_s<<" BS_array=\n";
		for (int i=0;i<g_size;++i) {
			for (int j=0;j<isize[1]+2*g_size;++j)
			std::cout<<BS[(i*(isize[1]+2*g_size)+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}
	sleep(1);
	if(procid==2) {
		std::cout<<"procid= "<<procid<<" dst_r="<<dst_r<<" GT=\n";
		for (int i=0;i<g_size;++i) {
			for (int j=0;j<isize[1]+2*g_size;++j)
			std::cout<<GT[(i*(isize[1]+2*g_size)+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}

	PCOUT<<"\nGB Col Communication\n";
#endif

	/* Phase 3: Now do the exact same thing for the right ghost side */
	int ts_buf_size = g_size * isize[2] * (isize[1] + 2 * g_size); // isize[1] now includes two side ghost cells
	//Real *TS=(Real*)accfft_alloc(ts_buf_size*sizeof(Real)); // Stores local right ghost data to be sent
  pvfmm::Iterator<Real> GB = pvfmm::aligned_new<Real>(ts_buf_size); // Left Ghost cells to be received
	// snafu: not really necessary to do memcpy, you can simply use padded_data directly
	//memcpy(TS,padded_data,ts_buf_size*sizeof(Real));
	Real *TS = &padded_data[0];

	/* Phase 4: Send your data to your right process
	 * First question is who is your right process?
	 */
	MPI_Request ts_s_request, ts_r_request;
	dst_s = (procid_c - 1) % nprocs_c;
	dst_r = (procid_c + 1) % nprocs_c;
	if (procid_c == 0)
		dst_s = nprocs_c - 1;
	MPI_Isend(&TS[0], ts_buf_size, MPI_T, dst_s, 0, col_comm, &ts_s_request);
	MPI_Irecv(&GB[0], ts_buf_size, MPI_T, dst_r, 0, col_comm, &ts_r_request);
	MPI_Wait(&ts_s_request, &ierr);
	MPI_Wait(&ts_r_request, &ierr);

#ifdef VERBOSE2
	if(procid==0) {
		std::cout<<"procid= "<<procid<<" padded_array=\n";
		for (int i=0;i<isize[0];++i) {
			for (int j=0;j<isize[1]+2*g_size;++j)
			std::cout<<padded_data[(i*(isize[1]+2*g_size)+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}

		std::cout<<"procid= "<<procid<<" dst_s="<<dst_s<<" BS_array=\n";
		for (int i=0;i<g_size;++i) {
			for (int j=0;j<isize[1]+2*g_size;++j)
			std::cout<<TS[(i*(isize[1]+2*g_size)+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}
	sleep(1);
	if(procid==2) {
		std::cout<<"procid= "<<procid<<" dst_r="<<dst_r<<" GB=\n";
		for (int i=0;i<g_size;++i) {
			for (int j=0;j<isize[1]+2*g_size;++j)
			std::cout<<GB[(i*(isize[1]+2*g_size)+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}
#endif

	// Phase 5: Pack the data GT+ padded_data + GB
	memcpy(&ghost_data[0], &GT[0],
			g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(Real));
	memcpy(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)],
			&padded_data[0],
			isize[0] * isize[2] * (isize[1] + 2 * g_size) * sizeof(Real));
	memcpy(
			&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)
					+ isize[0] * isize[2] * (isize[1] + 2 * g_size)], &GB[0],
			g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(Real));

#ifdef VERBOSE2
	if(procid==0) {
		std::cout<<"\n final ghost data\n";
		for (int i=0;i<isize[0]+2*g_size;++i) {
			for (int j=0;j<isize[1]+2*g_size;++j)
			std::cout<<ghost_data[(i*(isize[1]+2*g_size)+j)*isize[2]]<<" ";
			std::cout<<"\n";
		}
	}
#endif

	//accfft_free(TS);
  pvfmm::aligned_delete<Real>(GB);
	//accfft_free(BS);
	//accfft_free(GT);
  pvfmm::aligned_delete<Real>(GT);
}

/*
 * Perform a periodic z padding of size g_size. This function must be called after the ghost_top_bottom.
 * The idea is to have a symmetric periodic padding in all directions not just x, and y.
 *
 * @param[out] ghost_data_z: The output of the function which pads the input array data, with the ghost
 * cells from neighboring processors.
 * @param[in] ghost_data: Input data that is already padded in x, and y direction with ghost cells
 * @param[in] g_size: The size of the ghost cell padding. Note that it cannot exceed the neighboring processor's
 * local data size
 * @param[in] isize_g: An integer array specifying ghost cell padded local sizes.
 * @param[in] plan: AccFFT R2C plan
 */
void ghost_z(Real *ghost_data_z, pvfmm::Iterator<Real> ghost_data, int g_size, int* isize_g,
		reg::RegOpt* m_Opt) {

	IntType * isize = m_Opt->m_Domain.isize;
	for (int i = 0; i < isize_g[0]; ++i)
		for (int j = 0; j < isize_g[1]; ++j) {
			memcpy(&ghost_data_z[(i * isize_g[1] + j) * isize_g[2]],
					&ghost_data[(i * isize_g[1] + j) * isize[2] + isize[2]
							- g_size], g_size * sizeof(Real));
			memcpy(&ghost_data_z[(i * isize_g[1] + j) * isize_g[2] + g_size],
					&ghost_data[(i * isize_g[1] + j) * isize[2]],
					isize[2] * sizeof(Real));
			memcpy(
					&ghost_data_z[(i * isize_g[1] + j) * isize_g[2] + g_size
							+ isize[2]],
					&ghost_data[(i * isize_g[1] + j) * isize[2]],
					g_size * sizeof(Real));
		}
	return;
}

/*
 * Returns the necessary memory allocation in Bytes for the ghost data, as well
 * the local ghost sizes when ghost cell padding is desired only in x and y directions (and not z direction
 * which is locally owned).
 * @param[in] plan: AccFFT plan
 * @param[in] g_size: The number of ghost cells desired. Note that g_size cannot be bigger than
 * the minimum isize in each dimension among all processors. This means that you cannot get ghost
 * cells from a processor that is not a neighbor of the calling processor.
 * @param[out] isize_g: The new local sizes after getting the ghost cells.
 * @param[out] istart_g: Returns the new istart after getting the ghost cells. Note that this is the global
 * istart of the ghost cells. So for example, if processor zero gets 3 ghost cells, the left ghost cells will
 * come from the last process because of the periodicity. Then the istart_g would be the index of those elements
 * (that originally resided in the last processor).
 */

size_t ghost_local_size(reg::RegOpt* m_Opt, int g_size,
		int * isize_g, int* istart_g) {

	size_t alloc_max = m_Opt->m_FFT.nalloc;
	IntType *isize = m_Opt->m_Domain.isize;
	IntType *istart = m_Opt->m_Domain.istart;
	IntType* n = m_Opt->m_Domain.nx;
	istart_g[2] = istart[2];
	isize_g[2] = isize[2];
    
	istart_g[0] = istart[0] - g_size;
	istart_g[1] = istart[1] - g_size;

	if (istart_g[0] < 0)
		istart_g[0] += n[0];
	if (istart_g[1] < 0)
		istart_g[1] += n[1];
    
	isize_g[0] = isize[0] + 2 * g_size;
	isize_g[1] = isize[1] + 2 * g_size;
	return (alloc_max + 2 * g_size * isize[2] * isize[0] * sizeof(Real)
			+ 2 * g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(Real));
}

/*size_t accfft_ghost_local_size_dft_r2c(FFTPlanType* plan, int g_size,
		int * isize_g, int* istart_g) {
  return accfft_ghost_local_size_dft_r2c((FFTPlanType*)plan, g_size, isize_g, istart_g);
}*/
/*
 * Gather the ghost cells for a real input when ghost cell padding is desired only in x and y
 * directions (and not z direction which is locally owned). This function currently has the following limitations:
 *  - The AccFFT plan has to be outplace. Note that for inplace R2C plans, the input array has to be
 *  padded, which would slightly change the pattern of communicating ghost cells.
 *  - The number of ghost cells needed has to be the same in both x and y directions (note that each
 *  processor owns the whole z direction locally, so typically you would not need ghost cells. You can
 *  call accfft_get_ghost_xyz which would actually get ghost cells in z direction as well, but that is
 *  like a periodic padding in z direction).
 *  - The number of ghost cells requested cannot be bigger than the minimum isize of all processors. This
 *  means that a process cannot get a ghost element that does not belong to its immediate neighbors. This
 *  is a limitation that should barely matter, as typically ghost_size is a small integer while the global
 *  array sizes are very large. In mathematical terms g_size< min(isize[0],isize[1]) among isize of all
 *  processors.
 *
 * @param[in] plan: AccFFT plan
 * @param[in] g_size: The number of ghost cells desired. Note that g_size cannot be bigger than
 * the minimum isize in each dimension among all processors. This means that you cannot get ghost
 * cells from a processor that is not a neighbor of the calling processor.
 * @param[in] isize_g: An integer array specifying ghost cell padded local sizes.
 * @param[in] data: The local data whose ghost cells from other processors are sought.
 * @param[out] ghost_data: An array that is the ghost cell padded version of the input data.
 */
void get_ghost(reg::RegOpt* m_Opt, int g_size, int* isize_g, Real* data,
		Real* ghost_data) {
	int nprocs, procid;
	MPI_Comm_rank(m_Opt->m_Domain.mpicomm, &procid);
	MPI_Comm_size(m_Opt->m_Domain.mpicomm, &nprocs);

	/*if (plan->inplace == true) {
		PCOUT << "accfft_get_ghost_r2c does not support inplace transforms."
				<< std::endl;
		return;
	}*/

	if (g_size == 0) {
		memcpy(ghost_data, data, m_Opt->m_FFT.nalloc);
		return;
	}

	IntType *isize = m_Opt->m_Domain.isize;
	//int *istart = plan->istart;
	//int *n = plan->N;
	if (g_size > isize[0] || g_size > isize[1]) {
		std::cout
				<< "accfft_get_ghost_r2c does not support g_size greater than isize."
				<< std::endl;
		return;
	}

  pvfmm::Iterator<Real> padded_data = pvfmm::aligned_new<Real>
    (m_Opt->m_FFT.nalloc + 2 * g_size * isize[2] * isize[0]);
	ghost_left_right(padded_data, data, g_size, m_Opt);
	ghost_top_bottom(ghost_data, padded_data, g_size, m_Opt);
  pvfmm::aligned_delete<Real>(padded_data);
	return;

}

/*void accfft_get_ghost(accfft_plan* plan, int g_size, int* isize_g, Real* data,
		Real* ghost_data) {
  accfft_get_ghost((accfft_plan_t<Real, TC, PL>*)plan, g_size, isize_g, data, ghost_data);
}*/
/*
 * Returns the necessary memory allocation in Bytes for the ghost data, as well
 * the local ghost sizes when padding in all directions (including z direction that is locally owned by each process).
 * @param[in] plan: AccFFT plan
 * @param[in] g_size: The number of ghost cells desired. Note that g_size cannot be bigger than
 * the minimum isize in each dimension among all processors. This means that you cannot get ghost
 * cells from a processor that is not a neighbor of the calling processor.
 * @param[out] isize_g: The new local sizes after getting the ghost cells.
 * @param[out] istart_g: Returns the new istart after getting the ghost cells. Note that this is the global
 * istart of the ghost cells. So for example, if processor zero gets 3 ghost cells, the left ghost cells will
 * come from the last process because of the periodicity. Then the istart_g would be the index of those elements
 * (that originally resided in the last processor).
 */

size_t ghost_xyz_local_size(reg::RegOpt* m_Opt, int g_size,
		int * isize_g, int* istart_g) {

	size_t alloc_max = m_Opt->m_FFT.nalloc;
	IntType *isize = m_Opt->m_Domain.isize;
	IntType *istart = m_Opt->m_Domain.istart;
	IntType* n = m_Opt->m_Domain.nx;

	istart_g[0] = istart[0] - g_size;
	istart_g[1] = istart[1] - g_size;
	istart_g[2] = istart[2] - g_size;

	if (istart_g[0] < 0)
		istart_g[0] += n[0];
	if (istart_g[1] < 0)
		istart_g[1] += n[1];
	if (istart_g[2] < 0)
		istart_g[2] += n[2];

	isize_g[0] = isize[0] + 2 * g_size;
	isize_g[1] = isize[1] + 2 * g_size;
	isize_g[2] = isize[2] + 2 * g_size;

  //printf("\nalloc_max = %zu", alloc_max);
  //printf("\nisize_g[0] = %d", isize_g[0]);
  //printf("\nisize_g[1] = %d", isize_g[1]);
  //printf("\nisize_g[2] = %d", isize_g[2]);

  size_t alloc_max_g = alloc_max + 2 * g_size * isize[2] * isize[0] * sizeof(Real)
			+ 2 * g_size * isize[2] * isize_g[1] * sizeof(Real)
			+ 2 * g_size * isize_g[0] * isize_g[1] * sizeof(Real);
  //printf("\nalloc_max_g before = %zu", alloc_max_g);

  alloc_max_g += (16*isize_g[2]*isize_g[1]+16*isize_g[1]+16)*sizeof(Real); // to account for padding required for peeled loop in interp
  //printf("\nalloc_max_g after = %zu", alloc_max_g); 

  return alloc_max_g;
}

/*size_t accfft_ghost_xyz_local_size_dft_r2c(accfft_planf* plan, int g_size,
		int * isize_g, int* istart_g) {
  return accfft_ghost_xyz_local_size_dft_r2c((accfft_plan_t<Real, TC, PL>*)plan, g_size, isize_g, istart_g);
}

size_t accfft_ghost_xyz_local_size_dft_r2c(accfft_plan* plan, int g_size,
		int * isize_g, int* istart_g) {
  return accfft_ghost_xyz_local_size_dft_r2c((accfft_plan_t<Real, TC, PL>*)plan, g_size, isize_g, istart_g);
}*/
/*
 * Gather the ghost cells for a real input when ghost cell padding is desired in all directions including z direction
 * (which is locally owned). This function currently has the following limitations:
 *  - The AccFFT plan has to be outplace. Note that for inplace R2C plans, the input array has to be
 *  padded, which would slightly change the pattern of communicating ghost cells.
 *  - The number of ghost cells needed has to be the same in both x and y directions (note that each
 *  processor owns the whole z direction locally, so typically you would not need ghost cells. You can
 *  call accfft_get_ghost_xyz which would actually get ghost cells in z direction as well, but that is
 *  like a periodic padding in z direction).
 *  - The number of ghost cells requested cannot be bigger than the minimum isize of all processors. This
 *  means that a process cannot get a ghost element that does not belong to its immediate neighbors. This
 *  is a limitation that should barely matter, as typically ghost_size is a small integer while the global
 *  array sizes are very large. In mathematical terms g_size< min(isize[0],isize[1],isize[2]) among isize of all
 *  processors.
 *
 * @param[in] plan: AccFFT plan
 * @param[in] g_size: The number of ghost cells desired. Note that g_size cannot be bigger than
 * the minimum isize in each dimension among all processors. This means that you cannot get ghost
 * cells from a processor that is not a neighbor of the calling processor.
 * @param[in] isize_g: An integer array specifying ghost cell padded local sizes.
 * @param[in] data: The local data whose ghost cells from other processors are sought.
 * @param[out] ghost_data: An array that is the ghost cell padded version of the input data.
 */
void get_ghost_xyz(reg::RegOpt* m_Opt, int g_size, int* isize_g,
		Real* data, Real* ghost_data) {
	int nprocs, procid;
	MPI_Comm_rank(m_Opt->m_Domain.mpicomm, &procid);
	MPI_Comm_size(m_Opt->m_Domain.mpicomm, &nprocs);

	/*if (plan->inplace == true) {
		PCOUT << "accfft_get_ghost_r2c does not support inplace transforms."
				<< std::endl;
		return;
	}*/

	if (g_size == 0) {
		memcpy(ghost_data, data, m_Opt->m_FFT.nalloc);
		return;
	}

	IntType *isize = m_Opt->m_Domain.isize;
//	int *istart = plan->istart;
//	int *n = plan->N;
	if (g_size > isize[0] || g_size > isize[1]) {
		std::cout
				<< "accfft_get_ghost_r2c does not support g_size greater than isize."
				<< std::endl;
		return;
	}

  pvfmm::Iterator<Real> padded_data = pvfmm::aligned_new<Real>
    (m_Opt->m_FFT.nalloc + 2 * g_size * isize[2] * isize[0]);
  pvfmm::Iterator<Real> ghost_data_xy = pvfmm::aligned_new<Real>(
			m_Opt->m_FFT.nalloc + 2 * g_size * isize[2] * isize[0]
					+ 2 * g_size * isize[2] * isize_g[1]);

	ghost_left_right(padded_data, data, g_size, m_Opt);
	ghost_top_bottom(ghost_data_xy, padded_data, g_size, m_Opt);
	ghost_z(&ghost_data[0], ghost_data_xy, g_size, isize_g, m_Opt);

//#ifdef VERBOSE2
#ifdef VERBOSE3
	if(procid==0) {
		std::cout<<"\n final ghost data\n";
		for (int i=0;i<isize_g[0];++i) {
			for (int j=0;j<isize_g[1];++j) {
			    for (int k=0;k<isize_g[2];++k) {
			        int idx = reg::GetLinearIndex(i,j,k,isize_g);
                    //std::cout<<ghost_data[(i*isize_g[1]+j)*isize_g[2]]<<" ";
                    std::cout<<ghost_data[idx]<<" ";
                }
                std::cout<<"\n";
            }
            std::cout<<"\n\n";
        }
   /* 
		std::cout<<"\n a random z\n";
		int i=3+0*isize_g[0]/2;
		int j=3+0*isize_g[1]/2;
		for(int k=0;k<isize_g[2];++k)
		std::cout<<ghost_data[(i*isize_g[1]+j)*isize_g[2]+k]<<" ";
		std::cout<<"\n";
		*/
	}
#endif

  pvfmm::aligned_delete<Real>(padded_data);
  pvfmm::aligned_delete<Real>(ghost_data_xy);
	return;
}

void share_ghost_layer(reg::RegOpt* m_Opt, int g_size, int* isize_g,
		Real* data, Real* ghost_data, pvfmm::Iterator<Real> padded_data, pvfmm::Iterator<Real> ghost_data_xy) {
	
	int nprocs, procid;
	MPI_Comm_rank(m_Opt->m_Domain.mpicomm, &procid);
	MPI_Comm_size(m_Opt->m_Domain.mpicomm, &nprocs);

	/*if (plan->inplace == true) {
		PCOUT << "accfft_get_ghost_r2c does not support inplace transforms."
				<< std::endl;
		return;
	}*/

	if (g_size == 0) {
		memcpy(ghost_data, data, m_Opt->m_FFT.nalloc);
		return;
	}

	IntType *isize = m_Opt->m_Domain.isize;
//	int *istart = plan->istart;
//	int *n = plan->N;
	if (g_size > isize[0] || g_size > isize[1]) {
		std::cout
				<< "accfft_get_ghost_r2c does not support g_size greater than isize."
				<< std::endl;
		return;
	}

	ghost_left_right(padded_data, data, g_size, m_Opt);
	ghost_top_bottom(ghost_data_xy, padded_data, g_size, m_Opt);
	ghost_z(&ghost_data[0], ghost_data_xy, g_size, isize_g, m_Opt);

//#ifdef VERBOSE2
#ifdef VERBOSE3
	if(procid==0) {
		std::cout<<"\n final ghost data\n";
		for (int i=0;i<isize_g[0];++i) {
			for (int j=0;j<isize_g[1];++j) {
			    for (int k=0;k<isize_g[2];++k) {
			        int idx = reg::GetLinearIndex(i,j,k,isize_g);
                    //std::cout<<ghost_data[(i*isize_g[1]+j)*isize_g[2]]<<" ";
                    std::cout<<ghost_data[idx]<<" ";
                }
                std::cout<<"\n";
            }
            std::cout<<"\n\n";
        }
	}
#endif

	return;
}

/*void accfft_get_ghost_xyz(accfft_plan* plan, int g_size, int* isize_g,
		Real* data, Real* ghost_data) {
  accfft_get_ghost_xyz((accfft_plan_t<Real, TC, PL>*)plan, g_size, isize_g, data, ghost_data);
}*/
