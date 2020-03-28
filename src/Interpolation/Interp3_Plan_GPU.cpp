#include <interp3_gpu_mpi.hpp>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime_api.h>
#include <cmath>

#include <algorithm>
#include <iostream>

#ifndef ACCFFT_CHECKCUDA_H
#define ACCFFT_CHECKCUDA_H
inline cudaError_t checkCuda_accfft(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
  return result;
}
inline cufftResult checkCuda_accfft(cufftResult result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != CUFFT_SUCCESS) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", result);
    assert(result == CUFFT_SUCCESS);
  }
#endif
  return result;
}
#endif


class Trip_GPU{
  public:
    Trip_GPU() {};
    double x;
    double y;
    double z;
    int ind;
    int N[3];
    double h[3];

};
static bool ValueCmp(Trip_GPU const & a, Trip_GPU const & b)
{
    return a.z + a.y/a.h[1]*a.N[2] + a.x/a.h[0]* a.N[1]*a.N[2]<b.z + b.y/b.h[1]*b.N[2] + b.x/b.h[0]* b.N[1]*b.N[2] ;
}
static void sort_queries(std::vector<Real>* query_outside,std::vector<int>* f_index,int* N_reg,Real* h,MPI_Comm c_comm) {

  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);
  for(int proc=0;proc<nprocs;++proc) {
    int qsize=query_outside[proc].size()/COORD_DIM;
    Trip_GPU* trip=new Trip_GPU[qsize];

    for(int i=0;i<qsize;++i) {
      trip[i].x=query_outside[proc][i*COORD_DIM+0];
      trip[i].y=query_outside[proc][i*COORD_DIM+1];
      trip[i].z=query_outside[proc][i*COORD_DIM+2];
      trip[i].ind=f_index[proc][i];
      trip[i].N[0]=N_reg[0];
      trip[i].N[1]=N_reg[1];
      trip[i].N[2]=N_reg[2];
      trip[i].h[0]=h[0];
      trip[i].h[1]=h[1];
      trip[i].h[2]=h[2];
    }

    std::sort(trip, trip + qsize, ValueCmp);

    query_outside[proc].clear();
    f_index[proc].clear();

    for(int i=0;i<qsize;++i) {
      query_outside[proc].push_back(trip[i].x);
      query_outside[proc].push_back(trip[i].y);
      query_outside[proc].push_back(trip[i].z);
      f_index[proc].push_back(trip[i].ind);
    }
    delete[] trip;
  }
  return;
}



/**********************************************************************************
 * @brief plan constructor
 * *******************************************************************************/
Interp3_Plan_GPU::Interp3_Plan_GPU (size_t g_alloc_max) {
  this->g_alloc_max=g_alloc_max;
  this->allocate_baked=false;
  this->scatter_baked=false;
  int procs_i_recv_from_size_ = 0;
  int procs_i_send_to_size_ = 0;
}

/**********************************************************************************
 * @brief allocate memory for the plan
 * *******************************************************************************/
void Interp3_Plan_GPU::allocate (int N_pts, int* data_dof, int nplans)
{
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  
  // offset in the all_query_points array
  f_index_procs_others_offset = pvfmm::aligned_new<int>(nprocs); 
  // offset in the query_outside array
  f_index_procs_self_offset = pvfmm::aligned_new<int>(nprocs);
  // sizes of the number of interpolations that need to be sent to procs
  f_index_procs_self_sizes = pvfmm::aligned_new<int>(nprocs);
  // sizes of the number of interpolations that need to be received from procs
  f_index_procs_others_sizes = pvfmm::aligned_new<int>(nprocs);

	s_request = pvfmm::aligned_new<MPI_Request>(nprocs);
	request = pvfmm::aligned_new<MPI_Request>(2*nprocs);
  
  // copy of query points
	query_points = pvfmm::aligned_new<Real>(N_pts*COORD_DIM);
  pvfmm::memset(query_points,0, N_pts*COORD_DIM);
  
  f_index = new std::vector<int>[nprocs];
  query_outside = new std::vector<Real>[nprocs];
  
  // number of reuses of the plan with the same scatter points
  this->nplans_ = nplans; 
  this->data_dofs = pvfmm::aligned_new<int>(nplans);
  
  int max = 0;
  for(int i = 0; i < nplans_; ++i) {
    max = std::max(max, data_dof[i]);
    this->data_dofs[i] = data_dof[i];
  }
  this->data_dof_max = max;
  
  // The reshuffled semi-final interpolated values are stored here
  f_cubic_unordered = pvfmm::aligned_new<Real>(N_pts*data_dof_max);
  pvfmm::memset(f_cubic_unordered, 0, N_pts*data_dof_max);
  
#ifdef INTERP_PINNED
  cudaMalloc((void**)&this->ghost_reg_grid_vals_d, g_alloc_max*data_dof_max);
#else
  cudaMalloc((void**)&this->ghost_reg_grid_vals_d, g_alloc_max*data_dof_max);
#endif
  
  // strided for multiple plan calls
	stypes = pvfmm::aligned_new<MPI_Datatype>(nprocs*nplans_); 
	rtypes = pvfmm::aligned_new<MPI_Datatype>(nprocs*nplans_);
  this->allocate_baked = true;
}

/**********************************************************************************
 * @brief interp plan destructor
 * *******************************************************************************/
Interp3_Plan_GPU::~Interp3_Plan_GPU ()
{
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if(this->allocate_baked) {
    pvfmm::aligned_delete<int>(f_index_procs_others_offset);
		pvfmm::aligned_delete<int>(f_index_procs_self_offset);
		pvfmm::aligned_delete<int>(f_index_procs_self_sizes);
		pvfmm::aligned_delete<int>(f_index_procs_others_sizes);
    
    pvfmm::aligned_delete<MPI_Request>(s_request);
		pvfmm::aligned_delete<MPI_Request>(request);
    
    pvfmm::aligned_delete<Real>(query_points);
    pvfmm::aligned_delete<Real>(f_cubic_unordered);
    
    for(int proc=0;proc<nprocs;++proc) {
      std::vector<int>().swap(f_index[proc]);
      std::vector<Real>().swap(query_outside[proc]);
    }
    pvfmm::aligned_delete<MPI_Datatype>(rtypes);
    pvfmm::aligned_delete<MPI_Datatype>(stypes);
    pvfmm::aligned_delete<int>(data_dofs);
  }

  if(this->scatter_baked) {
#ifdef INTERP_PINNED
    cudaFreeHost(ghost_reg_grid_vals_d);
    cudaFreeHost(all_f_cubic_d);
    cudaFreeHost(xq1);
    cudaFreeHost(xq2);
    cudaFreeHost(xq3);
    cudaFreeHost(all_query_points_d);
#else
    cudaFree(ghost_reg_grid_vals_d);
    cudaFree(all_f_cubic_d);
    cudaFree(xq1);
    cudaFree(xq2);
    cudaFree(xq3);
    cudaFree(all_query_points_d);
#endif

		for (int ver = 0; ver < nplans_; ++ver)
		for (int i = 0; i < nprocs; ++i) {
			MPI_Type_free(&stypes[i+ver*nprocs]);
			MPI_Type_free(&rtypes[i+ver*nprocs]);
		}
    pvfmm::aligned_delete<Real>(all_query_points);
    pvfmm::aligned_delete<Real>(all_f_cubic);
	}


  return;
}


/*********************************************************************************************************************
 * Phase 1 of the parallel interpolation: This function computes which query_points needs to be sent to
 * other processors and which ones can be interpolated locally. Then a sparse alltoall is performed and
 * all the necessary information is sent/received including the coordinates for the query_points.
 * At the end, each process has the coordinates for interpolation of its own data and those of the others.
 *
 * IMPORTANT: This function must be called just once for a specific query_points. The reason is because of the
 * optimizations performed which assumes that the query_points do not change. For repeated interpolation you should
 * just call this function once, and instead repeatedly call Interp3_Plan::interpolate function.
 ********************************************************************************************************************/
void Interp3_Plan_GPU::scatter( int* N_reg,  // global grid dimensions
                                int * isize, // local grid dimensions
                                int* istart, // local grid start indices
                                const int N_pts, // local grid point count
                                const int g_size, // ghost layer width
                                Real* query_points_in, // input query points
                                int* c_dims,  // process cartesian grid dimensions
                                MPI_Comm c_comm,  // MPI Comm
                                double * timings) {
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);

  if(this->allocate_baked==false) {
    std::cout<<"ERROR Interp3_Plan_GPU Scatter called before calling allocate.\n";
    return;
  }
  if(this->scatter_baked==true) {
    for(int proc=0;proc<nprocs;++proc) {
      std::vector<int>().swap(f_index[proc]);
      std::vector<Real>().swap(query_outside[proc]);
    }
  }
  all_query_points_allocation=0;

  {
    N_reg_g[0]=N_reg[0]+2*g_size;
    N_reg_g[1]=N_reg[1]+2*g_size;
    N_reg_g[2]=N_reg[2]+2*g_size;

    isize_g[0]=isize[0]+2*g_size;
    isize_g[1]=isize[1]+2*g_size;
    isize_g[2]=isize[2]+2*g_size;

    Real h[3];
    h[0]=1./N_reg[0];
    h[1]=1./N_reg[1];
    h[2]=1./N_reg[2];

    // We copy query_points_in to query_points to aviod overwriting the input coordinates
    memcpy(query_points, query_points_in, N_pts*COORD_DIM*sizeof(Real));
  
    ZeitGeist_define(INTERPOL_MISC);
    ZeitGeist_tick(INTERPOL_MISC);
    // Enforce periodicity
    for(int i=0;i<N_pts;i++) {
      while(query_points[i*COORD_DIM+0]<=-h[0]) {query_points[i*COORD_DIM+0]=query_points[i*COORD_DIM+0]+1;}
      while(query_points[i*COORD_DIM+1]<=-h[1]) {query_points[i*COORD_DIM+1]=query_points[i*COORD_DIM+1]+1;}
      while(query_points[i*COORD_DIM+2]<=-h[2]) {query_points[i*COORD_DIM+2]=query_points[i*COORD_DIM+2]+1;}

      while(query_points[i*COORD_DIM+0]>=1) {query_points[i*COORD_DIM+0]=query_points[i*COORD_DIM+0]-1;}
      while(query_points[i*COORD_DIM+1]>=1) {query_points[i*COORD_DIM+1]=query_points[i*COORD_DIM+1]-1;}
      while(query_points[i*COORD_DIM+2]>=1) {query_points[i*COORD_DIM+2]=query_points[i*COORD_DIM+2]-1;}
    }
    ZeitGeist_tock(INTERPOL_MISC);
    
    // Compute the start and end coordinates that this processor owns
    Real iX0[3],iX1[3];
    for (int j=0;j<3;j++) {
      iX0[j]=istart[j]*h[j];
      iX1[j]=iX0[j]+(isize[j]-1)*h[j];
    }

    // Now march through the query points and split them into nprocs parts.
    // These are stored in query_outside which is an array of vectors of size nprocs.
    // That is query_outside[i] is a vector that contains the query points that need to
    // be sent to process i. Obviously for the case of query_outside[procid], we do not
    // need to send it to any other processor, as we own the necessary information locally,
    // and interpolation can be done locally.
    int Q_local=0, Q_outside=0;

    // This is needed for one-to-one correspondence with output f. This is becaues we are reshuffling
    // the data according to which processor it land onto, and we need to somehow keep the original
    // index to write the interpolation data back to the right location in the output.

    // This is necessary because when we want to compute dproc0 and dproc1 we have to divide by
    // the max isize. If the proc grid is unbalanced, the last proc's isize will be different
    // than others. With this approach we always use the right isize0 for all procs.
    
    ZeitGeist_define(INTERPOL_MISC1);
    ZeitGeist_tick(INTERPOL_MISC1);
    int isize0=std::ceil(N_reg[0]*1./c_dims[0]);
    int isize1=std::ceil(N_reg[1]*1./c_dims[1]);
    for(int i=0;i<N_pts;i++) {
      // The if condition check whether the query points fall into the locally owned domain or not
      if(
          iX0[0]-h[0]<=query_points[i*COORD_DIM+0] && query_points[i*COORD_DIM+0]<=iX1[0]+h[0] &&
          iX0[1]-h[1]<=query_points[i*COORD_DIM+1] && query_points[i*COORD_DIM+1]<=iX1[1]+h[1] &&
          iX0[2]-h[2]<=query_points[i*COORD_DIM+2] && query_points[i*COORD_DIM+2]<=iX1[2]+h[2]
        ) {
        query_outside[procid].push_back(query_points[i*COORD_DIM+0]);
        query_outside[procid].push_back(query_points[i*COORD_DIM+1]);
        query_outside[procid].push_back(query_points[i*COORD_DIM+2]);
        f_index[procid].push_back(i);
        Q_local++;
        continue;
      }
      else{
        // If the point does not reside in the processor's domain then we have to
        // first compute which processor owns the point. After computing that
        // we add the query point to the corresponding vector.
        int dproc0=(int)(query_points[i*COORD_DIM+0]/h[0])/isize0;
        int dproc1=(int)(query_points[i*COORD_DIM+1]/h[1])/isize1;
        // Compute which proc has to do the interpolation
        int proc=dproc0*c_dims[1]+dproc1; 
        
        query_outside[proc].push_back(query_points[i*COORD_DIM+0]);
        query_outside[proc].push_back(query_points[i*COORD_DIM+1]);
        query_outside[proc].push_back(query_points[i*COORD_DIM+2]);
        f_index[proc].push_back(i);
        Q_outside++;
        continue;
      }

    }
    ZeitGeist_tock(INTERPOL_MISC1);
    
    // Now sort the query points in zyx order
#ifdef SORT_QUERIES
    timings[3]+=-MPI_Wtime();
    sort_queries(query_outside,f_index,N_reg,h,c_comm);
    timings[3]+=+MPI_Wtime();
#endif

    // Now we need to send the query_points that land onto other processor's domain.
    // This done using a sparse alltoallv.
    // Right now each process knows how much data to send to others, but does not know
    // how much data it should receive. This is a necessary information both for the MPI
    // command as well as memory allocation for received data.
    // So we first do an alltoall to get the f_index[proc].size from all processes.

    ZeitGeist_define(INTERPOL_COMM);
    ZeitGeist_tick(INTERPOL_COMM);
    for (int proc=0;proc<nprocs;proc++) {
      if(!f_index[proc].empty())
        f_index_procs_self_sizes[proc]=f_index[proc].size();
      else
        f_index_procs_self_sizes[proc]=0;
    }
    timings[0]+=-MPI_Wtime();
    MPI_Alltoall(f_index_procs_self_sizes,1, MPI_INT,
        f_index_procs_others_sizes,1, MPI_INT,
        c_comm);
    timings[0]+=+MPI_Wtime();
    ZeitGeist_tock(INTERPOL_COMM);

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
    f_index_procs_others_offset[0]=0;
    f_index_procs_self_offset[0]=0;
    for (int proc=0;proc<nprocs;++proc) {
      // The reason we multiply by COORD_DIM is that we have three coordinates per interpolation request
      all_query_points_allocation+=f_index_procs_others_sizes[proc]*COORD_DIM;
      if(proc>0) {
        f_index_procs_others_offset[proc]=f_index_procs_others_offset[proc-1]+f_index_procs_others_sizes[proc-1];
        f_index_procs_self_offset[proc]=f_index_procs_self_offset[proc-1]+f_index_procs_self_sizes[proc-1];
      }
    }
    
    total_query_points=all_query_points_allocation/COORD_DIM;
    
    ZeitGeist_define(INTERPOL_MEMALLOC);
    ZeitGeist_tick(INTERPOL_MEMALLOC);
    // This if condition is to allow multiple calls to scatter fucntion with different query points
    // without having to create a new plan
    if(this->scatter_baked==true) {
      pvfmm::aligned_delete<Real>(this->all_query_points);
      pvfmm::aligned_delete<Real>(this->all_f_cubic);
			all_query_points = pvfmm::aligned_new<Real>(all_query_points_allocation);
			all_f_cubic = pvfmm::aligned_new<Real>(total_query_points * data_dof_max);
    } else {
			all_query_points = pvfmm::aligned_new<Real>(all_query_points_allocation);
			all_f_cubic = pvfmm::aligned_new<Real>(total_query_points * data_dof_max);
    }
    ZeitGeist_tock(INTERPOL_MEMALLOC);
    
    ZeitGeist_tick(INTERPOL_COMM);
    // Now perform the allotall to send/recv query_points
    timings[0]+=-MPI_Wtime();
    {
      int dst_r,dst_s;
      for (int i=0;i<nprocs;++i) {
        dst_r=i;
        dst_s=i;
        s_request[dst_s]=MPI_REQUEST_NULL;
        request[dst_r]=MPI_REQUEST_NULL;
        // notice that COORD_DIM is needed because query_points are 3 times f
        int roffset=f_index_procs_others_offset[dst_r]*COORD_DIM; 
        int soffset=f_index_procs_self_offset[dst_s]*COORD_DIM;
        if(f_index_procs_others_sizes[dst_r]!=0)
          MPI_Irecv(&all_query_points[roffset],f_index_procs_others_sizes[dst_r]*COORD_DIM,MPI_T, dst_r,
              0, c_comm, &request[dst_r]);
        if(!query_outside[dst_s].empty())
          MPI_Isend(&query_outside[dst_s][0],f_index_procs_self_sizes[dst_s]*COORD_DIM,MPI_T,dst_s,
              0, c_comm, &s_request[dst_s]);
      }
      // Wait for all the communication to finish
      MPI_Status ierr;
      for (int proc=0;proc<nprocs;++proc) {
        if(request[proc]!=MPI_REQUEST_NULL)
          MPI_Wait(&request[proc], &ierr);
        if(s_request[proc]!=MPI_REQUEST_NULL)
          MPI_Wait(&s_request[proc], &ierr);
      }
    }
    timings[0]+=+MPI_Wtime();
    ZeitGeist_tock(INTERPOL_COMM);

    // Now perform the interpolation on all query points including those that need to
    // be sent to other processors and store them into all_f_cubic
  }

  for(int ver = 0; ver < nplans_; ++ver) {
    for (int i = 0; i < nprocs; ++i) {
      MPI_Type_vector(data_dofs[ver], f_index_procs_self_sizes[i], N_pts, MPI_T, &rtypes[i+ver*nprocs]);
      MPI_Type_vector(data_dofs[ver], f_index_procs_others_sizes[i], total_query_points, MPI_T, &stypes[i+ver*nprocs]);
      MPI_Type_commit(&stypes[i+ver*nprocs]);
      MPI_Type_commit(&rtypes[i+ver*nprocs]);
    }
  }

  // This if condition is to allow multiple calls to scatter fucntion with different query points
  // without having to create a new plan
  ZeitGeist_define(INTERPOL_MEMALLOC);
  ZeitGeist_tick(INTERPOL_MEMALLOC);
  if(this->scatter_baked==true) {
#ifdef INTERP_PINNED
    cudaFreeHost(this->all_query_points_d);
    cudaFreeHost(this->xq1);
    cudaFreeHost(this->xq2);
    cudaFreeHost(this->xq3);
    cudaFreeHost(this->all_f_cubic_d);
    cudaMallocHost((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
    cudaMallocHost((void**)&xq1, total_query_points*sizeof(Real));
    cudaMallocHost((void**)&xq2, total_query_points*sizeof(Real));
    cudaMallocHost((void**)&xq3, total_query_points*sizeof(Real));
    cudaMallocHost((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof_max);
#else
    // freeing the cuda memory is required everytime scatter is called because the distribution of query points might not be uniform across all GPUs
    cudaFree(this->all_query_points_d);
    cudaFree(this->xq1);
    cudaFree(this->xq2);
    cudaFree(this->xq3);
    cudaFree(this->all_f_cubic_d);
    cudaMalloc((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof_max);
    cudaMalloc((void**)&xq1, total_query_points*sizeof(Real));
    cudaMalloc((void**)&xq2, total_query_points*sizeof(Real));
    cudaMalloc((void**)&xq3, total_query_points*sizeof(Real));
    cudaMalloc((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
#endif
  }
  else{
    // if this the first time scatter is being called (scatter_baked = False) then, only allocate the memory on GPU
#ifdef INTERP_PINNED
    cudaMallocHost((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
    cudaMallocHost((void**)&xq1, total_query_points*sizeof(Real));
    cudaMallocHost((void**)&xq2, total_query_points*sizeof(Real));
    cudaMallocHost((void**)&xq3, total_query_points*sizeof(Real));
    cudaMallocHost((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof_max);
#else
    cudaMalloc((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof_max);
    cudaMalloc((void**)&xq1, total_query_points*sizeof(Real));
    cudaMalloc((void**)&xq2, total_query_points*sizeof(Real));
    cudaMalloc((void**)&xq3, total_query_points*sizeof(Real));
    cudaMalloc((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
#endif
  }
  ZeitGeist_tock(INTERPOL_MEMALLOC);


  int proc_coord[2];
  proc_coord[0] = static_cast<int>(istart[0]/isize[0]);
  proc_coord[1] = static_cast<int>(istart[1]/isize[1]);

  ZeitGeist_define(INTERPOL_MEMCPY);
  ZeitGeist_tick(INTERPOL_MEMCPY);
  timings[2]+=-MPI_Wtime();
  cudaMemcpy(all_query_points_d, all_query_points, all_query_points_allocation*sizeof(Real), cudaMemcpyHostToDevice);
  ZeitGeist_tock(INTERPOL_MEMCPY);
  
  ZeitGeist_define(INTERPOL_MISC);
  ZeitGeist_tick(INTERPOL_MISC);
  normalizeQueryPoints(xq1, xq2, xq3, all_query_points_d, total_query_points, isize, N_reg, proc_coord, g_size);
  ZeitGeist_tock(INTERPOL_MISC);
  timings[2]+=+MPI_Wtime();
    
  

  this->scatter_baked=true;
  return;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Phase 2 of the parallel interpolation: This function must be called after the scatter function is called.
 * It performs local interpolation for all the points that the processor has for itself, as well as the interpolations
 * that it has to send to other processors. After the local interpolation is performed, a sparse
 * alltoall is performed so that all the interpolated results are sent/received.
 *
 */

void Interp3_Plan_GPU::interpolate( Real* ghost_reg_grid_vals, // ghost padded regular grid values on CPU
                                    int* N_reg,                // size of global grid points 
                                    int* isize,                // size of the local grid owned by the process
                                    int* istart,               // start point of the local grid owned by the process
                                    int* isize_g,              // size of the local grid (including ghost points)
                                    const int nlghost,         // number of local grid points (including ghost points) owned by process
                                    const int N_pts,           // number of local points owned by the process
                                    const int g_size,          // ghost layer width
                                    Real* query_values,        // interpolation result on CPU
                                    int* c_dims,               // dimensions of the communicator plan
                                    MPI_Comm c_comm,           // MPI communicator
                                    double * timings,          // time variable to store interpolation time
                                    float *tmp1,               // temporary memory for interpolation prefilter
                                    float* tmp2,               // temporary memory for interpolation prefilter
                                    cudaTextureObject_t yi_tex,// texture object for interpolation
                                    int iporder,               // interpolation order
                                    ScalarType* interp_time,   // interpolation time
                                    int version)               // plan version
{
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);
  if(this->allocate_baked==false) {
    std::cout<<"ERROR Interp3_Plan_GPU interpolate called before calling allocate.\n";
    return;
  }
  if(this->scatter_baked==false) {
    std::cout<<"ERROR Interp3_Plan_GPU interpolate called before calling scatter.\n";
    return;
  }

  int data_dof = data_dofs[version];

  ZeitGeist_define(INTERPOL_MEMCPY);
  ZeitGeist_tick(INTERPOL_MEMCPY);
  timings[2]+=-MPI_Wtime();
  cudaMemcpy((void*)ghost_reg_grid_vals_d, (const void*)ghost_reg_grid_vals, nlghost*data_dof*sizeof(Real), cudaMemcpyHostToDevice);
  timings[2]+=+MPI_Wtime();
  ZeitGeist_tock(INTERPOL_MEMCPY);


//#if defined(VERBOSE1) 
//  printf("\ng_alloc_max = %zu", g_alloc_max);
//  printf("\ndata_dof = %d", data_dof);
//  printf("\nisize_g[0] = %d", isize_g[0]);
//  printf("\nisize_g[1] = %d", isize_g[1]);
//  printf("\nisize_g[2] = %d", isize_g[2]);
//  printf("\nipoder = %d", iporder);
//  printf("\ntotal_query_points = %d", total_query_points);
//  printf("\nnlghost = %d", nlghost);
//#endif
  
  ZeitGeist_define(INTERPOL_EXEC);
  ZeitGeist_tick(INTERPOL_EXEC);
  timings[1]+=-MPI_Wtime();
  if (data_dof == 3)
    gpuInterpVec3D(&ghost_reg_grid_vals_d[0*nlghost], 
                   &ghost_reg_grid_vals_d[1*nlghost], 
                   &ghost_reg_grid_vals_d[2*nlghost], 
                   xq1, xq2, xq3, 
                   &all_f_cubic_d[0*total_query_points], 
                   &all_f_cubic_d[1*total_query_points], 
                   &all_f_cubic_d[2*total_query_points], 
                   tmp1, tmp2, isize_g, static_cast<long int>(total_query_points), yi_tex, iporder, interp_time);
  else 
    gpuInterp3D(ghost_reg_grid_vals_d, 
                xq1, xq2, xq3, 
                all_f_cubic_d, 
                tmp1, tmp2, isize_g, static_cast<long int>(total_query_points), yi_tex, 
                iporder, interp_time);
    
  timings[1]+=+MPI_Wtime();
  ZeitGeist_tock(INTERPOL_EXEC);
    
  // copy the interpolated results from the device to the host
  timings[2]+=-MPI_Wtime();
  ZeitGeist_tick(INTERPOL_MEMCPY);
  cudaMemcpy(all_f_cubic, all_f_cubic_d, data_dof*total_query_points*sizeof(Real) ,cudaMemcpyDeviceToHost);
  timings[2]+=+MPI_Wtime();
  ZeitGeist_tock(INTERPOL_MEMCPY);

  ZeitGeist_define(INTERPOL_COMM);
  ZeitGeist_tick(INTERPOL_COMM);
  // Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to f_cubic_unordered.
  timings[0]+=-MPI_Wtime();
  {
    int dst_r,dst_s;
		for (int i = procs_i_send_to_size_-1; i >=0; --i) {
			dst_r = procs_i_send_to_[i];    
			request[dst_r] = MPI_REQUEST_NULL; 
			int roffset = f_index_procs_self_offset[dst_r];
			MPI_Irecv(&f_cubic_unordered[roffset], 1, rtypes[dst_r+version*nprocs], dst_r, 0,
					c_comm, &request[dst_r]);
		}
		for (int i = 0; i < procs_i_recv_from_size_; ++i) {
			dst_s = procs_i_recv_from_[i];    
			s_request[dst_s] = MPI_REQUEST_NULL; 
			int soffset = f_index_procs_others_offset[dst_s];
			MPI_Isend(&all_f_cubic[soffset], 1, stypes[dst_s+version*nprocs], dst_s, 0, c_comm,
					&s_request[dst_s]);
		}
  }
  timings[0]+=+MPI_Wtime();
  ZeitGeist_tock(INTERPOL_COMM);
  
  // reshuffle interpolated values
  ZeitGeist_define(INTERPOL_MISC);
  for (int i = 0; i < procs_i_send_to_size_; ++i) {
    int proc = procs_i_send_to_[i];    
    if (request[proc] != MPI_REQUEST_NULL)
      ZeitGeist_tick(INTERPOL_COMM);
      MPI_Wait(&request[proc], MPI_STATUS_IGNORE);
      ZeitGeist_tock(INTERPOL_COMM);
      for (int dof = 0; dof < data_dofs[version]; ++dof) {
        Real* ptr = &f_cubic_unordered[f_index_procs_self_offset[proc]+dof*N_pts];
        ZeitGeist_tick(INTERPOL_MISC);
#pragma omp parallel for
        for (int i = 0; i < (int)f_index[proc].size(); ++i) {
          int ind = f_index[proc][i];
          query_values[ind + dof * N_pts] =ptr[i];
        }
        ZeitGeist_tock(INTERPOL_MISC);
      }
  }

  // wait for send
  ZeitGeist_tick(INTERPOL_COMM);
  for (int i = 0; i < procs_i_recv_from_size_; ++i) {
    int proc = procs_i_recv_from_[i];    //(procid+i)%nprocs;
    if (s_request[proc] != MPI_REQUEST_NULL)
      MPI_Wait(&s_request[proc], MPI_STATUS_IGNORE);
    }
  ZeitGeist_tock(INTERPOL_COMM);

  return;
}

//    for (int i=0;i<nprocs;++i) {
//      dst_r=i;
//      dst_s=i;
//      s_request[dst_s]=MPI_REQUEST_NULL;
//      request[dst_r]=MPI_REQUEST_NULL;
//      // Notice that this is the adjoint of the first comm part
//      // because now you are sending others f and receiving your part of f
//      int soffset=f_index_procs_others_offset[dst_r];
//      int roffset=f_index_procs_self_offset[dst_s];
//      if(f_index_procs_self_sizes[dst_r]!=0)
//        MPI_Irecv(&f_cubic_unordered[roffset],1,rtype[i], dst_r,
//            0, c_comm, &request[dst_r]);
//      if(f_index_procs_others_sizes[dst_s]!=0)
//        MPI_Isend(&all_f_cubic[soffset],1,stype[i],dst_s,
//            0, c_comm, &s_request[dst_s]);
//    }
//    MPI_Status ierr;
//    for (int proc=0;proc<nprocs;++proc) {
//      if(request[proc]!=MPI_REQUEST_NULL)
//        MPI_Wait(&request[proc], &ierr);
//      if(s_request[proc]!=MPI_REQUEST_NULL)
//        MPI_Wait(&s_request[proc], &ierr);
//    }
//  }
//  
//  ZeitGeist_define(INTERPOL_MISC);
//  ZeitGeist_tick(INTERPOL_MISC);
//  // Now copy back f_cubic_unordered to f_cubic in the correct f_index
//  for(int dof=0;dof<data_dof;++dof) {
//    for(int proc=0;proc<nprocs;++proc) {
//      if(!f_index[proc].empty())
//        for(unsigned long i=0;i<f_index[proc].size();++i) {
//          int ind=f_index[proc][i];
//          query_values[ind+dof*N_pts]=f_cubic_unordered[f_index_procs_self_offset[proc]+i+dof*N_pts];
//        }
//    }
//  }
//  ZeitGeist_tock(INTERPOL_MISC);
//  return;
//}
