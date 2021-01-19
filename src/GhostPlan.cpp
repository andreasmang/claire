#include "GhostPlan.hpp"

namespace reg {

GhostPlan::GhostPlan(RegOpt* m_Opt, IntType g_size) {
  this->m_Opt = m_Opt;
  this->g_size = static_cast<int>(g_size);
  
  int isize[3], istart[3];
  //int nl = this->m_Opt->m_Domain.nl;
  for (int i=0; i<3; i++) {
    isize[i] = m_Opt->m_Domain.isize[i];
    istart[i] = m_Opt->m_Domain.istart[i];
    this->isize[i] = isize[i];
    this->istart[i] = istart[i];
  }
  
  this->padded_data = nullptr;
  this->RS = nullptr;
  this->GL = nullptr;
  this->LS = nullptr;
  this->GR = nullptr;
  this->GB = nullptr;
  this->GT = nullptr;

/*
  // only needed for xy ghost
  cudaMalloc((void**)&padded_data, sizeof(ScalarType)*( nl + 2*this->g_size*isize[2]*isize[0] ));

	int rs_buf_size = this->g_size * isize[2] * isize[0];
  cudaMalloc((void**)&RS, rs_buf_size * sizeof(ScalarType));
  cudaMalloc((void**)&GL, rs_buf_size * sizeof(ScalarType));


	int ls_buf_size = this->g_size * isize[2] * isize[0];
  cudaMalloc((void**)&LS, ls_buf_size * sizeof(ScalarType));
  cudaMalloc((void**)&GR, ls_buf_size * sizeof(ScalarType));


	//int ts_buf_size = this->g_size * isize[2] * (isize[1] + 2*this->g_size); // for xy ghost
	int ts_buf_size = this->g_size * isize[2] * isize[1];
  cudaMalloc((void**)&GB, sizeof(ScalarType)*ts_buf_size);
*/  

	//int bs_buf_size = this->g_size * isize[2] * (isize[1] + 2 * this->g_size); // for xy ghost
	int bs_buf_size = this->g_size * isize[2] * isize[1]; // only one auxilliary is needed for ghost comm in x
  cudaMalloc((void**)&GT, sizeof(ScalarType)*bs_buf_size);
  
  stream = new cudaStream_t[num_streams];
  for (int i=0; i<num_streams; i++) {
      cudaStreamCreate(&stream[i]);
  }
}

size_t GhostPlan::get_ghost_local_size_x(IntType * isize_g, IntType* istart_g) {
	int nl = this->m_Opt->m_Domain.nl;
	int n[3];
  
  for (int i=0; i<3; i++) {
    n[i] = this->m_Opt->m_Domain.nx[i];
  }

	istart_g[2] = istart[2];
	isize_g[2] = isize[2];
	istart_g[1] = istart[1];
	isize_g[1] = isize[1];
    
	istart_g[0] = istart[0] - g_size;

	if (istart_g[0] < 0)
		istart_g[0] += n[0];
    
	isize_g[0] = isize[0] + 2 * g_size;
	
	for (int i=0; i<3; i++) {
	  this->isize_g[i] = isize_g[i];
	  this->istart_g[i] = istart_g[i];
  }
	
	this->g_alloc_max = (nl + 2 * g_size * isize[2] * isize[1]) * sizeof(ScalarType);	
	return this->g_alloc_max;
}

size_t GhostPlan::get_ghost_local_size_xy(IntType * isize_g, IntType* istart_g) {
	int nl = this->m_Opt->m_Domain.nl;
	int n[3];
  
  for (int i=0; i<3; i++) {
    n[i] = this->m_Opt->m_Domain.nx[i];
  }
	
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
	  this->isize_g[i] = isize_g[i];
	  this->istart_g[i] = istart_g[i];
  }
	
	this->g_alloc_max = (nl + 2 * g_size * isize[2] * isize[0]
			+ 2 * g_size * isize[2] * (isize[1] + 2 * g_size)) * sizeof(ScalarType);
	return this->g_alloc_max;
}

GhostPlan::~GhostPlan() {
  if (padded_data != nullptr) {
    cudaFree(padded_data);
    padded_data = nullptr;
  }

  if (RS != nullptr) {
    cudaFree(RS);
    RS = nullptr;
  }

  if (GL != nullptr) {
    cudaFree(GL);
    GL = nullptr;
  }
  
  if (LS != nullptr) {
    cudaFree(LS);
    LS = nullptr;
  }

  if (GR != nullptr) {
    cudaFree(GR);
    GR = nullptr;
  }
  
  if (GB != nullptr) {
    cudaFree(GB);
    GB = nullptr;
  }

  if (GT != nullptr) {
    cudaFree(GT);
    GT = nullptr;
  }

  for (int i=0; i<num_streams; i++) 
    cudaStreamDestroy(stream[i]);
  delete[] stream;
}

void GhostPlan::share_left_right(const ScalarType* data, double* timers) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	/* Get the local pencil size and the allocation size */
	//int* isize = plan->isize;

	MPI_Comm row_comm = this->m_Opt->m_Domain.rowcomm; 
	int nprocs_r, procid_r;
	MPI_Comm_rank(row_comm, &procid_r);
	MPI_Comm_size(row_comm, &nprocs_r);
	/* Halo Exchange along y axis
	 * Phase 1: Write local data to be sent to the right process to RS
	 */
#ifdef VERBOSE2
	PCOUT<<"\nGL Row Communication\n";
#endif

	for (int x = 0; x < isize[0]; ++x) {
    cudaMemcpyAsync((void*)&RS[x * g_size * isize[2]], 
                    (const void*)&data[x * isize[2] * isize[1] + (isize[1] - g_size) * isize[2]], 
                    g_size * isize[2] * sizeof(ScalarType), cudaMemcpyDeviceToDevice, stream[x%num_streams]);
  }
  for (int i=0; i<num_streams; i++)
    cudaStreamSynchronize(stream[i]);
  

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
	
  timers[0]+=-MPI_Wtime();
	MPI_Isend(RS, rs_buf_size, MPI_T, dst_s, 0, row_comm, &rs_s_request);
	MPI_Irecv(GL, rs_buf_size, MPI_T, dst_r, 0, row_comm, &rs_r_request);
	MPI_Wait(&rs_s_request, &ierr);
	MPI_Wait(&rs_r_request, &ierr);
	timers[0]+=+MPI_Wtime();

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

	for (int x = 0; x < isize[0]; ++x) {
		cudaMemcpyAsync(&LS[x * g_size * isize[2]], 
		                &data[x * isize[2] * isize[1]],
		                g_size * isize[2] * sizeof(ScalarType), cudaMemcpyDeviceToDevice, stream[x%num_streams]);
  }
  for (int i=0; i<num_streams; i++)
    cudaStreamSynchronize (stream[i]);
  

	/* Phase 4: Send your data to your right process
	 * First question is who is your right process?
	 */
	dst_s = (procid_r - 1) % nprocs_r;
	dst_r = (procid_r + 1) % nprocs_r;
	if (procid_r == 0)
		dst_s = nprocs_r - 1;
  
	timers[0]+=-MPI_Wtime();
	MPI_Isend(LS, ls_buf_size, MPI_T, dst_s, 0, row_comm, &rs_s_request);
	MPI_Irecv(GR, ls_buf_size, MPI_T, dst_r, 0, row_comm, &rs_r_request);
	MPI_Wait(&rs_s_request, &ierr);
	MPI_Wait(&rs_r_request, &ierr);
	timers[0]+=+MPI_Wtime();


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
	// Phase 5: Pack the data GL+ data + GR
	for (int i = 0; i < isize[0]; ++i) {
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

	MPI_Comm col_comm = this->m_Opt->m_Domain.colcomm;
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
  
	timers[0]+=-MPI_Wtime();
	MPI_Isend(&BS[0], bs_buf_size, MPI_T, dst_s, 0, col_comm, &bs_s_request);
	MPI_Irecv(&GT[0], bs_buf_size, MPI_T, dst_r, 0, col_comm, &bs_r_request);
	MPI_Wait(&bs_s_request, &ierr);
	MPI_Wait(&bs_r_request, &ierr);
	timers[0]+=+MPI_Wtime();

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
  
	timers[0]+=-MPI_Wtime();
	MPI_Isend(&TS[0], ts_buf_size, MPI_T, dst_s, 0, col_comm, &ts_s_request);
	MPI_Irecv(&GB[0], ts_buf_size, MPI_T, dst_r, 0, col_comm, &ts_r_request);
	MPI_Wait(&ts_s_request, &ierr);
	MPI_Wait(&ts_r_request, &ierr);
	timers[0]+=+MPI_Wtime();

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
  
	cudaMemcpyAsync(&ghost_data[0], 
	            &GT[0], 
	            g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	cudaMemcpyAsync(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)], 
	            &padded_data[0], 
	            isize[0] * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	cudaMemcpyAsync(&ghost_data[g_size * isize[2] * (isize[1] + 2 * g_size)+ isize[0] * isize[2] * (isize[1] + 2 * g_size)], 
	            &GB[0],
	            g_size * isize[2] * (isize[1] + 2 * g_size) * sizeof(ScalarType), cudaMemcpyDeviceToDevice);

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

void GhostPlan::share_ghost_x(const ScalarType* data, ScalarType* ghost_data) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	double timers[4] = {0,0,0,0};

	if (g_size == 0) {
		cudaMemcpy(ghost_data, data, sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToDevice);
		return;
	}

	if (g_size > isize[0]) {
		std::cout
				<< "g_size greater than isize not supported."
				<< std::endl;
		return;
	}

	
	MPI_Comm col_comm = this->m_Opt->m_Domain.colcomm;
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
  const ScalarType* BS = &data[(isize[0] - g_size) * isize[2] * (isize[1])];


  /* Phase 2: Send your data to your bottom process
	 * First question is who is your bottom process?
	 */
	int dst_s = (procid_c + 1) % nprocs_c;
	int dst_r = (procid_c - 1) % nprocs_c;
	if (procid_c == 0)
		dst_r = nprocs_c - 1;
	MPI_Request bs_s_request, bs_r_request;
	MPI_Status ierr;
  
	timers[0]+=-MPI_Wtime();
	MPI_Isend(&BS[0], bs_buf_size, MPI_T, dst_s, 0, col_comm, &bs_s_request);
	MPI_Irecv(&GT[0], bs_buf_size, MPI_T, dst_r, 0, col_comm, &bs_r_request);
	timers[0]+=+MPI_Wtime();

  // Do cudamemcpys while waiting for MPI to finish
  // overlap communication and copying
	MPI_Wait(&bs_s_request, &ierr);
	MPI_Wait(&bs_r_request, &ierr);
  
	cudaMemcpy(&ghost_data[g_size * isize[2] * isize[1]], 
	            &data[0], 
	            isize[0] * isize[2] * isize[1] * sizeof(ScalarType), cudaMemcpyDeviceToDevice);
	cudaMemcpy(&ghost_data[0], 
	            &GT[0], 
	            g_size * isize[2] * isize[1] * sizeof(ScalarType), cudaMemcpyDeviceToDevice);

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
	const ScalarType* TS = &data[0];

	/* Phase 4: Send your data to your right process
	 * First question is who is your right process?
	 */
	MPI_Request ts_s_request, ts_r_request;
	dst_s = (procid_c - 1) % nprocs_c;
	dst_r = (procid_c + 1) % nprocs_c;
	if (procid_c == 0)
		dst_s = nprocs_c - 1;
  
	timers[0]+=-MPI_Wtime();
	MPI_Isend(&TS[0], ts_buf_size, MPI_T, dst_s, 0, col_comm, &ts_s_request);
	MPI_Irecv(&GT[0], ts_buf_size, MPI_T, dst_r, 0, col_comm, &ts_r_request);
	MPI_Wait(&ts_s_request, &ierr);
	MPI_Wait(&ts_r_request, &ierr);
	timers[0]+=+MPI_Wtime();

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

	cudaMemcpy(&ghost_data[(isize[0] + g_size) * isize[2] * isize[1]], 
	            &GT[0],
	            g_size * isize[2] * isize[1] * sizeof(ScalarType), cudaMemcpyDeviceToDevice);

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

  
  if (this->m_Opt->m_Verbosity > 2) reg::DbgMsg("ghost points shared");

	return;
}

void GhostPlan::share_ghost_xy(const ScalarType* data, ScalarType* ghost_data) {
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  double timers[4] = {0,0,0,0};

	if (g_size == 0) {
		cudaMemcpy(ghost_data, data, sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToDevice);
		return;
	}

	if (g_size > isize[0] || g_size > isize[1]) {
		std::cout
				<< "g_size greater than isize not supported"
				<< std::endl;
		return;
	}


	share_left_right(data,timers);
	share_top_bottom(ghost_data, timers);
	return;
}

} // namespace reg
