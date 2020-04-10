#include "GhostPlan.hpp"

namespace reg {

GhostPlan::GhostPlan(FFTPlanType* plan, int g_size) {
  this->plan = plan;
  this->g_size = g_size;
}

size_t GhostPlan::get_ghost_local_size_x(int * isize_g, int* istart_g) {
	size_t alloc_max = plan->alloc_max;
	int *isize = plan->isize;
	int *istart = plan->istart;
	int* n = plan->N;
	istart_g[2] = istart[2];
	isize_g[2] = isize[2];
	istart_g[1] = istart[1];
	isize_g[1] = isize[1];
    
	istart_g[0] = istart[0] - g_size;

	if (istart_g[0] < 0)
		istart_g[0] += n[0];
    
	isize_g[0] = isize[0] + 2 * g_size;
	
	for (int i=0; i<3; i++) {
	  this->isize[i] = isize[i];
	  this->isize_g[i] = isize_g[i];
	  this->istart[i] = istart[i];
	  this->istart_g[i] = istart_g[i];
  }
	
	this->g_alloc_max = alloc_max + 2 * g_size * isize[2] * isize[1] * sizeof(ScalarType);	
	return this->g_alloc_max;
}

size_t GhostPlan::get_ghost_local_size_xy(int * isize_g, int* istart_g) {
	size_t alloc_max = plan->alloc_max;
	int *isize = plan->isize;
	int *istart = plan->istart;
	int* n = plan->N;
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
	
	for (int i=0; i<3; i++) {
	  this->isize[i] = isize[i];
	  this->isize_g[i] = isize_g[i];
	  this->istart[i] = istart[i];
	  this->istart_g[i] = istart_g[i];
  }
	
	this->g_alloc_max = alloc_max + 2 * g_size * isize[2] * isize[0] * sizeof(ScalarType)
			+ 2 * g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType);
	
	return this->g_alloc_max;
}

void GhostPlan::allocate() {
  cudaMalloc((void**)&padded_data, sizeof(ScalarType)*( plan->alloc_max + 2*g_size*isize[2]*isize[0] ));

	int rs_buf_size = g_size * isize[2] * isize[0];
  cudaMalloc((void**)&RS, rs_buf_size * sizeof(ScalarType));
  cudaMalloc((void**)&GL, rs_buf_size * sizeof(ScalarType));


	int ls_buf_size = g_size * isize[2] * isize[0];
  cudaMalloc((void**)&LS, ls_buf_size * sizeof(ScalarType));
  cudaMalloc((void**)&GR, ls_buf_size * sizeof(ScalarType));


	int ts_buf_size = g_size * isize[2] * (isize[1] + 2 * g_size);
  cudaMalloc((void**)&GB, sizeof(ScalarType)*ts_buf_size);
  

	int bs_buf_size = g_size * isize[2] * (isize[1] + 2 * g_size);
  cudaMalloc((void**)&GT, sizeof(ScalarType)*bs_buf_size);
  
  stream = new cudaStream_t[num_streams];
  for (int i=0; i<num_streams; i++) {
      cudaStreamCreate(&stream[i]);
  }
}

GhostPlan::~GhostPlan() {
  cudaFree(padded_data);
  cudaFree(RS);
  cudaFree(GL);
  cudaFree(LS);
  cudaFree(GR);
  cudaFree(GB);
  cudaFree(GT);
  for (int i=0; i<num_streams; i++) 
    cudaStreamDestroy(stream[i]);
  delete[] stream;
}

void GhostPlan::share_left_right(ScalarType* data, double* timers) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	/* Get the local pencil size and the allocation size */
	//int* isize = plan->isize;

	MPI_Comm row_comm = plan->row_comm;
	int nprocs_r, procid_r;
	MPI_Comm_rank(row_comm, &procid_r);
	MPI_Comm_size(row_comm, &nprocs_r);
	/* Halo Exchange along y axis
	 * Phase 1: Write local data to be sent to the right process to RS
	 */
#ifdef VERBOSE2
	PCOUT<<"\nGL Row Communication\n";
#endif

  ZeitGeist_define(ghost_cudamemcpy);
  ZeitGeist_tick(ghost_cudamemcpy);
	for (int x = 0; x < isize[0]; ++x) {
	//	reg::gencpy(&RS[x * g_size * isize[2]],
	//			    &data[x * isize[2] * isize[1] + (isize[1] - g_size) * isize[2]],
	//			    g_size * isize[2] * sizeof(ScalarType));
    cudaMemcpyAsync((void*)&RS[x * g_size * isize[2]], 
                    (const void*)&data[x * isize[2] * isize[1] + (isize[1] - g_size) * isize[2]], 
                    g_size * isize[2] * sizeof(ScalarType), cudaMemcpyDeviceToDevice, stream[x%num_streams]);
  }
  for (int i=0; i<num_streams; i++)
    cudaStreamSynchronize(stream[i]);
  ZeitGeist_tock(ghost_cudamemcpy);
  

	/* Phase 2: Send your data to your right process
	 * First question is who is your right process?
	 */
	int dst_s = (procid_r + 1) % nprocs_r;
	int dst_r = (procid_r - 1) % nprocs_r;
	if (procid_r == 0)
		dst_r = nprocs_r - 1;
	
	MPI_Request rs_s_request, rs_r_request;
	MPI_Status ierr;

	int rs_buf_size = g_size * isize[2] * isize[0];
	
	ZeitGeist_define(ghost_comm);
	ZeitGeist_tick(ghost_comm);
  timers[0]+=-MPI_Wtime();
	MPI_Isend(RS, rs_buf_size, MPI_T, dst_s, 0, row_comm, &rs_s_request);
	MPI_Irecv(GL, rs_buf_size, MPI_T, dst_r, 0, row_comm, &rs_r_request);
	MPI_Wait(&rs_s_request, &ierr);
	MPI_Wait(&rs_r_request, &ierr);
	timers[0]+=+MPI_Wtime();
	ZeitGeist_tock(ghost_comm);

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

  ZeitGeist_tick(ghost_cudamemcpy);
	for (int x = 0; x < isize[0]; ++x) {
		//reg::gencpy(&LS[x * g_size * isize[2]], 
		//             &data[x * isize[2] * isize[1]],
		//             g_size * isize[2] * sizeof(ScalarType));
		cudaMemcpyAsync(&LS[x * g_size * isize[2]], 
		                &data[x * isize[2] * isize[1]],
		                g_size * isize[2] * sizeof(ScalarType), cudaMemcpyDeviceToDevice, stream[x%num_streams]);
  }
  for (int i=0; i<num_streams; i++)
    cudaStreamSynchronize (stream[i]);
  ZeitGeist_tock(ghost_cudamemcpy);
  

	/* Phase 4: Send your data to your right process
	 * First question is who is your right process?
	 */
	dst_s = (procid_r - 1) % nprocs_r;
	dst_r = (procid_r + 1) % nprocs_r;
	if (procid_r == 0)
		dst_s = nprocs_r - 1;
  
  ZeitGeist_tick(ghost_comm);
	timers[0]+=-MPI_Wtime();
	MPI_Isend(LS, ls_buf_size, MPI_T, dst_s, 0, row_comm, &rs_s_request);
	MPI_Irecv(GR, ls_buf_size, MPI_T, dst_r, 0, row_comm, &rs_r_request);
	MPI_Wait(&rs_s_request, &ierr);
	MPI_Wait(&rs_r_request, &ierr);
	timers[0]+=+MPI_Wtime();
  ZeitGeist_tock(ghost_comm);


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
  
  int k = 0;
  ZeitGeist_tick(ghost_cudamemcpy);
	// Phase 5: Pack the data GL+ data + GR
	for (int i = 0; i < isize[0]; ++i) {
		//reg::gencpy(&padded_data[i * isize[2] * (isize[1] + 2 * g_size)],
		//		    &GL[i * g_size * isize[2]], 
		//		    g_size * isize[2] * sizeof(ScalarType));
		//reg::gencpy(&padded_data[i * isize[2] * (isize[1] + 2 * g_size) + g_size * isize[2]], 
		//		    &data[i * isize[2] * isize[1]],
		//		    isize[1] * isize[2] * sizeof(ScalarType));
		//reg::gencpy(&padded_data[i * isize[2] * (isize[1] + 2 * g_size) + g_size * isize[2] + isize[2] * isize[1]],
		//		    &GR[i * g_size * isize[2]], 
		//		    g_size * isize[2] * sizeof(ScalarType));
		
		cudaMemcpyAsync(&padded_data[i * isize[2] * (isize[1] + 2 * g_size)],
				    &GL[i * g_size * isize[2]], 
				    g_size * isize[2] * sizeof(ScalarType), cudaMemcpyDeviceToDevice, stream[k%num_streams]);
		k++;
		cudaMemcpyAsync(&padded_data[i * isize[2] * (isize[1] + 2 * g_size) + g_size * isize[2]], 
				    &data[i * isize[2] * isize[1]],
				    isize[1] * isize[2] * sizeof(ScalarType), cudaMemcpyDeviceToDevice, stream[k%num_streams]);
		k++;
		cudaMemcpyAsync(&padded_data[i * isize[2] * (isize[1] + 2 * g_size) + g_size * isize[2] + isize[2] * isize[1]],
				    &GR[i * g_size * isize[2]], 
				    g_size * isize[2] * sizeof(ScalarType), cudaMemcpyDeviceToDevice, stream[k%num_streams]);
	  k++; 	
	}
  for (int i=0; i<num_streams; i++)
    cudaStreamSynchronize (stream[i]);
	ZeitGeist_tock(ghost_cudamemcpy);
  

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

}

void GhostPlan::share_top_bottom(ScalarType* ghost_data,  double* timers) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


	/* Get the local pencil size and the allocation size */
	//int isize[3],osize[3],istart[3],ostart[3];
	//int * isize = plan->isize;

	MPI_Comm col_comm = plan->col_comm;
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
  
	// snafu: not really necessary to do memcpy, you can simply use padded_data directly
  ScalarType* BS = &padded_data[(isize[0] - g_size) * isize[2] * (isize[1] + 2 * g_size)];


  /* Phase 2: Send your data to your bottom process
	 * First question is who is your bottom process?
	 */
	int dst_s = (procid_c + 1) % nprocs_c;
	int dst_r = (procid_c - 1) % nprocs_c;
	if (procid_c == 0)
		dst_r = nprocs_c - 1;
	MPI_Request bs_s_request, bs_r_request;
	MPI_Status ierr;
  
  ZeitGeist_define(ghost_comm);
  ZeitGeist_tick(ghost_comm);
	timers[0]+=-MPI_Wtime();
	MPI_Isend(&BS[0], bs_buf_size, MPI_T, dst_s, 0, col_comm, &bs_s_request);
	MPI_Irecv(&GT[0], bs_buf_size, MPI_T, dst_r, 0, col_comm, &bs_r_request);
	MPI_Wait(&bs_s_request, &ierr);
	MPI_Wait(&bs_r_request, &ierr);
	timers[0]+=+MPI_Wtime();
  ZeitGeist_tock(ghost_comm);

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

	// snafu: not really necessary to do memcpy, you can simply use padded_data directly
	//memcpy(TS,padded_data,ts_buf_size*sizeof(ScalarType));
	ScalarType* TS = &padded_data[0];

	/* Phase 4: Send your data to your right process
	 * First question is who is your right process?
	 */
	MPI_Request ts_s_request, ts_r_request;
	dst_s = (procid_c - 1) % nprocs_c;
	dst_r = (procid_c + 1) % nprocs_c;
	if (procid_c == 0)
		dst_s = nprocs_c - 1;
  
  ZeitGeist_tick(ghost_comm);
	timers[0]+=-MPI_Wtime();
	MPI_Isend(&TS[0], ts_buf_size, MPI_T, dst_s, 0, col_comm, &ts_s_request);
	MPI_Irecv(&GB[0], ts_buf_size, MPI_T, dst_r, 0, col_comm, &ts_r_request);
	MPI_Wait(&ts_s_request, &ierr);
	MPI_Wait(&ts_r_request, &ierr);
	timers[0]+=+MPI_Wtime();
  ZeitGeist_tock(ghost_comm);

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
	//reg::gencpy(&ghost_data[0], 
	//            &GT[0], 
	//            g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType));
	//reg::gencpy(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)], 
	//            &padded_data[0], 
	//            isize[0] * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType));
	//reg::gencpy(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)+ isize[0] * isize[2] * (isize[1] + 2 * g_size)], 
	//            &GB[0],
	//            g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType));
  
  ZeitGeist_define(ghost_cudamemcpy);
  ZeitGeist_tick(ghost_cudamemcpy);
	cudaMemcpyAsync(&ghost_data[0], 
	            &GT[0], 
	            g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	cudaMemcpyAsync(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)], 
	            &padded_data[0], 
	            isize[0] * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	cudaMemcpyAsync(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)+ isize[0] * isize[2] * (isize[1] + 2 * g_size)], 
	            &GB[0],
	            g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	ZeitGeist_tock(ghost_cudamemcpy);

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
}

void GhostPlan::share_ghost_x(ScalarType* data, ScalarType* ghost_data, double* timers) {
	int nprocs, procid;
	MPI_Comm_rank(plan->c_comm, &procid);
	MPI_Comm_size(plan->c_comm, &nprocs);

	if (plan->inplace == true) {
		PCOUT << "accfft_get_ghost_r2c does not support inplace transforms."
				<< std::endl;
		return;
	}

	if (g_size == 0) {
		cudaMemcpy(ghost_data, data, plan->alloc_max, cudaMemcpyDeviceToDevice);
		return;
	}

	if (g_size > isize[0] || g_size > isize[1]) {
		std::cout
				<< "accfft_get_ghost_r2c does not support g_size greater than isize."
				<< std::endl;
		return;
	}

  ZeitGeist_define(ghost_comm_top_bottom);
	ZeitGeist_tick(ghost_comm_top_bottom);
	
	MPI_Comm col_comm = plan->col_comm;
	int nprocs_c, procid_c;
	MPI_Comm_rank(col_comm, &procid_c);
	MPI_Comm_size(col_comm, &nprocs_c);

	/* Halo Exchange along x axis
	 * Phase 1: Write local data to be sent to the bottom process
	 */
#ifdef VERBOSE2
	PCOUT<<"\nGB Col Communication\n";
#endif

	int bs_buf_size = g_size * isize[2] * (isize[1]);
  
	// snafu: not really necessary to do memcpy, you can simply use padded_data directly
  ScalarType* BS = &data[(isize[0] - g_size) * isize[2] * (isize[1])];


  /* Phase 2: Send your data to your bottom process
	 * First question is who is your bottom process?
	 */
	int dst_s = (procid_c + 1) % nprocs_c;
	int dst_r = (procid_c - 1) % nprocs_c;
	if (procid_c == 0)
		dst_r = nprocs_c - 1;
	MPI_Request bs_s_request, bs_r_request;
	MPI_Status ierr;
  
  ZeitGeist_define(ghost_comm);
  ZeitGeist_tick(ghost_comm);
	timers[0]+=-MPI_Wtime();
	MPI_Isend(&BS[0], bs_buf_size, MPI_T, dst_s, 0, col_comm, &bs_s_request);
	MPI_Irecv(&GT[0], bs_buf_size, MPI_T, dst_r, 0, col_comm, &bs_r_request);
	timers[0]+=+MPI_Wtime();
  ZeitGeist_tock(ghost_comm);

  // Do cudamemcpys while waiting for MPI to finish
  // overlap communication and copying
  ZeitGeist_define(ghost_cudamemcpy);
  ZeitGeist_tick(ghost_cudamemcpy);
	cudaMemcpyAsync(&ghost_data[g_size * isize[2] * isize[1]], 
	            &data[0], 
	            isize[0] * isize[2] * isize[1] * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	cudaMemcpyAsync(&ghost_data[0], 
	            &GT[0], 
	            g_size * isize[2] * isize[1] * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	ZeitGeist_tock(ghost_cudamemcpy);

  ZeitGeist_tick(ghost_comm);
	MPI_Wait(&bs_s_request, &ierr);
	MPI_Wait(&bs_r_request, &ierr);
  ZeitGeist_tock(ghost_comm);

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
	int ts_buf_size = g_size * isize[2] * (isize[1]); // isize[1] now includes two side ghost cells

	// snafu: not really necessary to do memcpy, you can simply use padded_data directly
	//memcpy(TS,padded_data,ts_buf_size*sizeof(ScalarType));
	ScalarType* TS = &data[0];

	/* Phase 4: Send your data to your right process
	 * First question is who is your right process?
	 */
	MPI_Request ts_s_request, ts_r_request;
	dst_s = (procid_c - 1) % nprocs_c;
	dst_r = (procid_c + 1) % nprocs_c;
	if (procid_c == 0)
		dst_s = nprocs_c - 1;
  
  ZeitGeist_tick(ghost_comm);
	timers[0]+=-MPI_Wtime();
	MPI_Isend(&TS[0], ts_buf_size, MPI_T, dst_s, 0, col_comm, &ts_s_request);
	MPI_Irecv(&GB[0], ts_buf_size, MPI_T, dst_r, 0, col_comm, &ts_r_request);
	MPI_Wait(&ts_s_request, &ierr);
	MPI_Wait(&ts_r_request, &ierr);
	timers[0]+=+MPI_Wtime();
  ZeitGeist_tock(ghost_comm);

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
	//reg::gencpy(&ghost_data[0], 
	//            &GT[0], 
	//            g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType));
	//reg::gencpy(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)], 
	//            &padded_data[0], 
	//            isize[0] * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType));
	//reg::gencpy(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)+ isize[0] * isize[2] * (isize[1] + 2 * g_size)], 
	//            &GB[0],
	//            g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType));
  
  ZeitGeist_tick(ghost_cudamemcpy);
	cudaMemcpyAsync(&ghost_data[(isize[0] + g_size) * isize[2] * isize[1]], 
	            &GB[0],
	            g_size * isize[2] * isize[1] * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	ZeitGeist_tock(ghost_cudamemcpy);

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

	ZeitGeist_tock(ghost_comm_top_bottom);
	return;
}

void GhostPlan::share_ghost_xy(ScalarType* data, ScalarType* ghost_data, double* timers) {
	int nprocs, procid;
	MPI_Comm_rank(plan->c_comm, &procid);
	MPI_Comm_size(plan->c_comm, &nprocs);

	if (plan->inplace == true) {
		PCOUT << "accfft_get_ghost_r2c does not support inplace transforms."
				<< std::endl;
		return;
	}

	if (g_size == 0) {
		cudaMemcpy(ghost_data, data, plan->alloc_max, cudaMemcpyDeviceToDevice);
		return;
	}

	//int *isize = plan->isize;
	//int *istart = plan->istart;
	//int *n = plan->N;
	if (g_size > isize[0] || g_size > isize[1]) {
		std::cout
				<< "accfft_get_ghost_r2c does not support g_size greater than isize."
				<< std::endl;
		return;
	}

  ZeitGeist_define(ghost_comm_left_right);
  ZeitGeist_define(ghost_comm_top_bottom);

  ZeitGeist_tick(ghost_comm_left_right);
	share_left_right(data,timers);
	ZeitGeist_tock(ghost_comm_left_right);
	ZeitGeist_tick(ghost_comm_top_bottom);
	share_top_bottom(ghost_data, timers);
	ZeitGeist_tock(ghost_comm_top_bottom);
	return;
}

} // namespace reg
