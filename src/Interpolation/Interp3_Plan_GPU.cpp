







#include <interp3_gpu_mpi.hpp>
#include <string.h>
#include <stdlib.h>
#include <vector>
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
    Trip_GPU(){};
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
static void sort_queries(std::vector<Real>* query_outside,std::vector<int>* f_index,int* N_reg,Real* h,MPI_Comm c_comm){

  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);
  for(int proc=0;proc<nprocs;++proc){
    int qsize=query_outside[proc].size()/COORD_DIM;
    Trip_GPU* trip=new Trip_GPU[qsize];

    for(int i=0;i<qsize;++i){
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

    for(int i=0;i<qsize;++i){
      query_outside[proc].push_back(trip[i].x);
      query_outside[proc].push_back(trip[i].y);
      query_outside[proc].push_back(trip[i].z);
      f_index[proc].push_back(trip[i].ind);
    }
    delete(trip);
  }
  return;
}



Interp3_Plan_GPU::Interp3_Plan_GPU (size_t g_alloc_max) {
  this->g_alloc_max=g_alloc_max;
  this->allocate_baked=false;
  this->scatter_baked=false;
}

void Interp3_Plan_GPU::allocate (int N_pts, int data_dof)
{
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  query_points=(Real*) malloc(N_pts*COORD_DIM*sizeof(Real));

  f_index_procs_others_offset=(int*)malloc(nprocs*sizeof(int)); // offset in the all_query_points array
  f_index_procs_self_offset  =(int*)malloc(nprocs*sizeof(int)); // offset in the query_outside array
  f_index_procs_self_sizes   =(int*)malloc(nprocs*sizeof(int)); // sizes of the number of interpolations that need to be sent to procs
  f_index_procs_others_sizes =(int*)malloc(nprocs*sizeof(int)); // sizes of the number of interpolations that need to be received from procs

  s_request= new MPI_Request[nprocs];
  request= new MPI_Request[nprocs];

  f_index = new std::vector<int> [nprocs];
  query_outside=new std::vector<Real> [nprocs];

  f_cubic_unordered=(Real*) malloc(N_pts*sizeof(Real)*data_dof); // The reshuffled semi-final interpolated values are stored here

  //double time=0;
  //time=-MPI_Wtime();
#ifdef INTERP_PINNED
  //cudaMallocHost((void**)&this->ghost_reg_grid_vals_d,g_alloc_max*data_dof);
  cudaMalloc((void**)&this->ghost_reg_grid_vals_d, g_alloc_max*data_dof);
#else
  cudaMalloc((void**)&this->ghost_reg_grid_vals_d, g_alloc_max*data_dof);
#endif

  //time+=MPI_Wtime();
  //if(procid==0)
  //  std::cout<<"malloc time="<<time<<std::endl;

  stype= new MPI_Datatype[nprocs];
  rtype= new MPI_Datatype[nprocs];
  this->data_dof=data_dof;
  this->allocate_baked=true;
}

Interp3_Plan_GPU::~Interp3_Plan_GPU ()
{
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if(this->allocate_baked){
    free(query_points);

    free(f_index_procs_others_offset);
    free(f_index_procs_self_offset  );
    free(f_index_procs_self_sizes   );
    free(f_index_procs_others_sizes );

    delete(s_request);
    delete(request);
    //vectors
    for(int proc=0;proc<nprocs;++proc)
    {
      std::vector<int>().swap(f_index[proc]);
      std::vector<Real>().swap(query_outside[proc]);
    }
    free(f_cubic_unordered);

  }

  if(this->scatter_baked){
    free(all_query_points);
    free(all_f_cubic);

#ifdef INTERP_PINNED
    //cudaFreeHost(ghost_reg_grid_vals_d);
    cudaFree(ghost_reg_grid_vals_d);
#else
    cudaFree(ghost_reg_grid_vals_d);
#endif




#ifdef INTERP_PINNED
    cudaFreeHost(all_f_cubic_d);
    cudaFreeHost(all_query_points_d);
#else
    cudaFree(all_f_cubic_d);
    cudaFree(all_query_points_d);
#endif

    for(int i=0;i<nprocs;++i){
      MPI_Type_free(&stype[i]);
      MPI_Type_free(&rtype[i]);
    }

  }

  if(this->allocate_baked){
    delete(stype);
    delete(rtype);
  }
  return;
}

void rescale_xyz(const int g_size,  int* N_reg, int* N_reg_g, int* istart, const int N_pts, Real* query_points);


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
void Interp3_Plan_GPU::scatter( int data_dof,
    int* N_reg, int * isize, int* istart, const int N_pts, const int g_size, Real* query_points_in,
    int* c_dims, MPI_Comm c_comm, double * timings)
{
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);

  if(this->allocate_baked==false){
    std::cout<<"ERROR Interp3_Plan_GPU Scatter called before calling allocate.\n";
    return;
  }
  if(this->scatter_baked==true){
    for(int proc=0;proc<nprocs;++proc)
    {
      std::vector<int>().swap(f_index[proc]);
      std::vector<Real>().swap(query_outside[proc]);
    }
  }
  all_query_points_allocation=0;

  {

    //int N_reg_g[3], isize_g[3];
    N_reg_g[0]=N_reg[0]+2*g_size;
    N_reg_g[1]=N_reg[1]+2*g_size;
    N_reg_g[2]=N_reg[2]+2*g_size;

    isize_g[0]=isize[0]+2*g_size;
    isize_g[1]=isize[1]+2*g_size;
    isize_g[2]=isize[2]+2*g_size;

    Real h[3]; // original grid size along each axis
    h[0]=1./N_reg[0];
    h[1]=1./N_reg[1];
    h[2]=1./N_reg[2];

    // We copy query_points_in to query_points to aviod overwriting the input coordinates
    Real* query_points=(Real*) malloc(N_pts*COORD_DIM*sizeof(Real));
    memcpy(query_points,query_points_in,N_pts*COORD_DIM*sizeof(Real));
    // Enforce periodicity
    for(int i=0;i<N_pts;i++){
      while(query_points[i*COORD_DIM+0]<=-h[0]) {query_points[i*COORD_DIM+0]=query_points[i*COORD_DIM+0]+1;}
      while(query_points[i*COORD_DIM+1]<=-h[1]) {query_points[i*COORD_DIM+1]=query_points[i*COORD_DIM+1]+1;}
      while(query_points[i*COORD_DIM+2]<=-h[2]) {query_points[i*COORD_DIM+2]=query_points[i*COORD_DIM+2]+1;}

      while(query_points[i*COORD_DIM+0]>=1) {query_points[i*COORD_DIM+0]=query_points[i*COORD_DIM+0]-1;}
      while(query_points[i*COORD_DIM+1]>=1) {query_points[i*COORD_DIM+1]=query_points[i*COORD_DIM+1]-1;}
      while(query_points[i*COORD_DIM+2]>=1) {query_points[i*COORD_DIM+2]=query_points[i*COORD_DIM+2]-1;}
    }


    // Compute the start and end coordinates that this processor owns
    Real iX0[3],iX1[3];
    for (int j=0;j<3;j++){
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
    int isize0=std::ceil(N_reg[0]*1./c_dims[0]);
    int isize1=std::ceil(N_reg[1]*1./c_dims[1]);
    for(int i=0;i<N_pts;i++){
      // The if condition check whether the query points fall into the locally owned domain or not
      if(
          iX0[0]-h[0]<=query_points[i*COORD_DIM+0] && query_points[i*COORD_DIM+0]<=iX1[0]+h[0] &&
          iX0[1]-h[1]<=query_points[i*COORD_DIM+1] && query_points[i*COORD_DIM+1]<=iX1[1]+h[1] &&
          iX0[2]-h[2]<=query_points[i*COORD_DIM+2] && query_points[i*COORD_DIM+2]<=iX1[2]+h[2]
        ){
        query_outside[procid].push_back(query_points[i*COORD_DIM+0]);
        query_outside[procid].push_back(query_points[i*COORD_DIM+1]);
        query_outside[procid].push_back(query_points[i*COORD_DIM+2]);
        f_index[procid].push_back(i);
        Q_local++;
        //PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
        continue;
      }
      else{
        // If the point does not reside in the processor's domain then we have to
        // first compute which processor owns the point. After computing that
        // we add the query point to the corresponding vector.
        int dproc0=(int)(query_points[i*COORD_DIM+0]/h[0])/isize0;
        int dproc1=(int)(query_points[i*COORD_DIM+1]/h[1])/isize1;
        int proc=dproc0*c_dims[1]+dproc1; // Compute which proc has to do the interpolation
        //if(proc>=nprocs) std::cout<<"dp0="<<dproc0<<" dp1="<<dproc1<<" proc="<<proc<<" q[0]="<<query_points[i*COORD_DIM+0]<<" h[0]="<<h[0]<< " div="<<(query_points[i*COORD_DIM+0]/h[0])<<" "<<isize[0]<<std::endl;
        //PCOUT<<"proc="<<proc<<std::endl;
        query_outside[proc].push_back(query_points[i*COORD_DIM+0]);
        query_outside[proc].push_back(query_points[i*COORD_DIM+1]);
        query_outside[proc].push_back(query_points[i*COORD_DIM+2]);
        f_index[proc].push_back(i);
        Q_outside++;
        //PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
        continue;
      }

    }


    // Now sort the query points in zyx order
#ifdef SORT_QUERIES
    timings[3]+=-MPI_Wtime();
    sort_queries(query_outside,f_index,N_reg,h,c_comm);
    timings[3]+=+MPI_Wtime();
    //if(procid==0) std::cout<<"Sorting Queries\n";
#endif

    // Now we need to send the query_points that land onto other processor's domain.
    // This done using a sparse alltoallv.
    // Right now each process knows how much data to send to others, but does not know
    // how much data it should receive. This is a necessary information both for the MPI
    // command as well as memory allocation for received data.
    // So we first do an alltoall to get the f_index[proc].size from all processes.


    for (int proc=0;proc<nprocs;proc++){
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


    // Now we need to allocate memory for the receiving buffer of all query
    // points including ours. This is simply done by looping through
    // f_index_procs_others_sizes and adding up all the sizes.
    // Note that we would also need to know the offsets.
    f_index_procs_others_offset[0]=0;
    f_index_procs_self_offset[0]=0;
    for (int proc=0;proc<nprocs;++proc){
      // The reason we multiply by COORD_DIM is that we have three coordinates per interpolation request
      all_query_points_allocation+=f_index_procs_others_sizes[proc]*COORD_DIM;
      if(proc>0){
        f_index_procs_others_offset[proc]=f_index_procs_others_offset[proc-1]+f_index_procs_others_sizes[proc-1];
        f_index_procs_self_offset[proc]=f_index_procs_self_offset[proc-1]+f_index_procs_self_sizes[proc-1];
      }
    }
    total_query_points=all_query_points_allocation/COORD_DIM;

    // This if condition is to allow multiple calls to scatter fucntion with different query points
    // without having to create a new plan
    if(this->scatter_baked==true){
      free(this->all_query_points);
      free(this->all_f_cubic);
      all_query_points=(Real*) malloc(all_query_points_allocation*sizeof(Real));
      all_f_cubic=(Real*)malloc(total_query_points*sizeof(Real)*data_dof);
    }
    else{
      all_query_points=(Real*) malloc(all_query_points_allocation*sizeof(Real));
      all_f_cubic=(Real*)malloc(total_query_points*sizeof(Real)*data_dof);
    }

    // Now perform the allotall to send/recv query_points
    timings[0]+=-MPI_Wtime();
    {
      int dst_r,dst_s;
      for (int i=0;i<nprocs;++i){
        dst_r=i;//(procid+i)%nprocs;
        dst_s=i;//(procid-i+nprocs)%nprocs;
        s_request[dst_s]=MPI_REQUEST_NULL;
        request[dst_r]=MPI_REQUEST_NULL;
        int roffset=f_index_procs_others_offset[dst_r]*COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
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
      for (int proc=0;proc<nprocs;++proc){
        if(request[proc]!=MPI_REQUEST_NULL)
          MPI_Wait(&request[proc], &ierr);
        if(s_request[proc]!=MPI_REQUEST_NULL)
          MPI_Wait(&s_request[proc], &ierr);
      }
    }
    timings[0]+=+MPI_Wtime();

    // Now perform the interpolation on all query points including those that need to
    // be sent to other processors and store them into all_f_cubic
    free(query_points);
  }

  for(int i=0;i<nprocs;++i){
    MPI_Type_vector(data_dof,f_index_procs_self_sizes[i],N_pts, MPI_T, &rtype[i]);
    MPI_Type_vector(data_dof,f_index_procs_others_sizes[i],total_query_points, MPI_T, &stype[i]);
    MPI_Type_commit(&stype[i]);
    MPI_Type_commit(&rtype[i]);
  }

  rescale_xyz(g_size,N_reg,N_reg_g,istart,total_query_points,all_query_points); //SNAFU

  // This if condition is to allow multiple calls to scatter fucntion with different query points
  // without having to create a new plan
  if(this->scatter_baked==true){
#ifdef INTERP_PINNED
    cudaFreeHost(this->all_query_points_d);
    cudaFreeHost(this->all_f_cubic_d);
    cudaMallocHost((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
    cudaMallocHost((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof);
#else
    cudaFree(this->all_query_points_d);
    cudaFree(this->all_f_cubic_d);
    cudaMalloc((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof);
    cudaMalloc((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
#endif
  }
  else{
#ifdef INTERP_PINNED
    cudaMallocHost((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
    cudaMallocHost((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof);
#else
    cudaMalloc((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof);
    cudaMalloc((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
#endif
  }

  timings[2]+=-MPI_Wtime();
  cudaMemcpy(all_query_points_d,all_query_points,all_query_points_allocation*sizeof(Real),cudaMemcpyHostToDevice);
  timings[2]+=+MPI_Wtime();


  this->scatter_baked=true;
  return;
}

/*
 * Phase 2 of the parallel interpolation: This function must be called after the scatter function is called.
 * It performs local interpolation for all the points that the processor has for itself, as well as the interpolations
 * that it has to send to other processors. After the local interpolation is performed, a sparse
 * alltoall is performed so that all the interpolated results are sent/received.
 *
 */
void Interp3_Plan_GPU::interpolate( Real* ghost_reg_grid_vals, int data_dof,
    int* N_reg, int * isize, int* istart, const int N_pts, const int g_size,
    Real* query_values,int* c_dims, MPI_Comm c_comm,double * timings)
{

  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);
  if(this->allocate_baked==false){
    std::cout<<"ERROR Interp3_Plan_GPU interpolate called before calling allocate.\n";
    return;
  }
  if(this->scatter_baked==false){
    std::cout<<"ERROR Interp3_Plan_GPU interpolate called before calling scatter.\n";
    return;
  }


  timings[2]+=-MPI_Wtime();
  cudaMemcpy(ghost_reg_grid_vals_d,ghost_reg_grid_vals,g_alloc_max*data_dof,cudaMemcpyHostToDevice);
  timings[2]+=+MPI_Wtime();

  timings[1]+=-MPI_Wtime();
  gpu_interp3_ghost_xyz_p(ghost_reg_grid_vals_d, data_dof, N_reg, isize,istart,total_query_points,g_size,all_query_points_d, all_f_cubic_d,true);
  timings[1]+=+MPI_Wtime();

  timings[2]+=-MPI_Wtime();
  cudaMemcpy(all_f_cubic,all_f_cubic_d, total_query_points*sizeof(Real)*data_dof ,cudaMemcpyDeviceToHost);
  timings[2]+=+MPI_Wtime();


  // Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to
  // f_cubic_unordered.
  timings[0]+=-MPI_Wtime();
  {
    int dst_r,dst_s;
    for (int i=0;i<nprocs;++i){
      dst_r=i;//(procid+i)%nprocs;
      dst_s=i;//(procid-i+nprocs)%nprocs;
      s_request[dst_s]=MPI_REQUEST_NULL;
      request[dst_r]=MPI_REQUEST_NULL;
      // Notice that this is the adjoint of the first comm part
      // because now you are sending others f and receiving your part of f
      int soffset=f_index_procs_others_offset[dst_r];
      int roffset=f_index_procs_self_offset[dst_s];
      if(f_index_procs_self_sizes[dst_r]!=0)
        MPI_Irecv(&f_cubic_unordered[roffset],1,rtype[i], dst_r,
            0, c_comm, &request[dst_r]);
      if(f_index_procs_others_sizes[dst_s]!=0)
        MPI_Isend(&all_f_cubic[soffset],1,stype[i],dst_s,
            0, c_comm, &s_request[dst_s]);
    }
    MPI_Status ierr;
    for (int proc=0;proc<nprocs;++proc){
      if(request[proc]!=MPI_REQUEST_NULL)
        MPI_Wait(&request[proc], &ierr);
      if(s_request[proc]!=MPI_REQUEST_NULL)
        MPI_Wait(&s_request[proc], &ierr);
    }
  }
  timings[0]+=+MPI_Wtime();

  // Now copy back f_cubic_unordered to f_cubic in the correct f_index
  for(int dof=0;dof<data_dof;++dof){
    for(int proc=0;proc<nprocs;++proc){
      if(!f_index[proc].empty())
        for(int i=0;i<f_index[proc].size();++i){
          int ind=f_index[proc][i];
          query_values[ind+dof*N_pts]=f_cubic_unordered[f_index_procs_self_offset[proc]+i+dof*N_pts];
        }
    }
  }

  return;
}


/*
 * A dummy function that performs the whole process of scatter and interpolation without any planning or
 * prior allocation. Use only for debugging.
 */
void Interp3_Plan_GPU::slow_run( Real* ghost_reg_grid_vals_d, int data_dof,
    int* N_reg, int * isize, int* istart, const int N_pts, const int g_size, Real* query_points_in,
    Real* query_values,int* c_dims, MPI_Comm c_comm)
{
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);

  //printf("====== isize[0]=%d\n",isize[0]);
  int N_reg_g[3], isize_g[3];
  N_reg_g[0]=N_reg[0]+2*g_size;
  N_reg_g[1]=N_reg[1]+2*g_size;
  N_reg_g[2]=N_reg[2]+2*g_size;

  isize_g[0]=isize[0]+2*g_size;
  isize_g[1]=isize[1]+2*g_size;
  isize_g[2]=isize[2]+2*g_size;

  Real h[3]; // original grid size along each axis
  h[0]=1./N_reg[0];
  h[1]=1./N_reg[1];
  h[2]=1./N_reg[2];

  // We copy query_points_in to query_points to aviod overwriting the input coordinates
  Real* query_points=(Real*) malloc(N_pts*COORD_DIM*sizeof(Real));
  memcpy(query_points,query_points_in,N_pts*COORD_DIM*sizeof(Real));
  // Enforce periodicity
  for(int i=0;i<N_pts;i++){
    while(query_points[i*COORD_DIM+0]<=-h[0]) {query_points[i*COORD_DIM+0]=query_points[i*COORD_DIM+0]+1;}
    while(query_points[i*COORD_DIM+1]<=-h[1]) {query_points[i*COORD_DIM+1]=query_points[i*COORD_DIM+1]+1;}
    while(query_points[i*COORD_DIM+2]<=-h[2]) {query_points[i*COORD_DIM+2]=query_points[i*COORD_DIM+2]+1;}

    while(query_points[i*COORD_DIM+0]>=1) {query_points[i*COORD_DIM+0]=query_points[i*COORD_DIM+0]-1;}
    while(query_points[i*COORD_DIM+1]>=1) {query_points[i*COORD_DIM+1]=query_points[i*COORD_DIM+1]-1;}
    while(query_points[i*COORD_DIM+2]>=1) {query_points[i*COORD_DIM+2]=query_points[i*COORD_DIM+2]-1;}
  }


  // Compute the start and end coordinates that this processor owns
  Real iX0[3],iX1[3];
  for (int j=0;j<3;j++){
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
  std::vector <int> f_index[nprocs];
  std::vector <Real> query_outside[nprocs];
  for(int i=0;i<N_pts;i++){
    // The if condition check whether the query points fall into the locally owned domain or not
    if(
        iX0[0]-h[0]<=query_points[i*COORD_DIM+0] && query_points[i*COORD_DIM+0]<=iX1[0]+h[0] &&
        iX0[1]-h[1]<=query_points[i*COORD_DIM+1] && query_points[i*COORD_DIM+1]<=iX1[1]+h[1] &&
        iX0[2]-h[2]<=query_points[i*COORD_DIM+2] && query_points[i*COORD_DIM+2]<=iX1[2]+h[2]
      ){
      query_outside[procid].push_back(query_points[i*COORD_DIM+0]);
      query_outside[procid].push_back(query_points[i*COORD_DIM+1]);
      query_outside[procid].push_back(query_points[i*COORD_DIM+2]);
      f_index[procid].push_back(i);
      Q_local++;
      //PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
      continue;
    }
    else{
      // If the point does not reside in the processor's domain then we have to
      // first compute which processor owns the point. After computing that
      // we add the query point to the corresponding vector.
      int dproc0=(int)(query_points[i*COORD_DIM+0]/h[0])/isize[0];
      int dproc1=(int)(query_points[i*COORD_DIM+1]/h[1])/isize[1];
      int proc=dproc0*c_dims[1]+dproc1; // Compute which proc has to do the interpolation
      //PCOUT<<"proc="<<proc<<std::endl;
      query_outside[proc].push_back(query_points[i*COORD_DIM+0]);
      query_outside[proc].push_back(query_points[i*COORD_DIM+1]);
      query_outside[proc].push_back(query_points[i*COORD_DIM+2]);
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

  for (int proc=0;proc<nprocs;proc++){
    if(!f_index[proc].empty())
      f_index_procs_self_sizes[proc]=f_index[proc].size();
    else
      f_index_procs_self_sizes[proc]=0;
  }
  MPI_Alltoall(f_index_procs_self_sizes,1, MPI_INT,
      f_index_procs_others_sizes,1, MPI_INT,
      c_comm);

#ifdef VERBOSE2
  sleep(1);
  if(procid==0){
    std::cout<<"procid="<<procid<<std::endl;
    std::cout<<"f_index_procs_self[0]="<<f_index_procs_self_sizes[0]<<" [1]= "<<f_index_procs_self_sizes[1]<<std::endl;
    std::cout<<"f_index_procs_others[0]="<<f_index_procs_others_sizes[0]<<" [1]= "<<f_index_procs_others_sizes[1]<<std::endl;
  }
  sleep(1);
  if(procid==1){
    std::cout<<"procid="<<procid<<std::endl;
    std::cout<<"f_index_procs_self[0]="<<f_index_procs_self_sizes[0]<<" [1]= "<<f_index_procs_self_sizes[1]<<std::endl;
    std::cout<<"f_index_procs_others[0]="<<f_index_procs_others_sizes[0]<<" [1]= "<<f_index_procs_others_sizes[1]<<std::endl;
  }
#endif


  // Now we need to allocate memory for the receiving buffer of all query
  // points including ours. This is simply done by looping through
  // f_index_procs_others_sizes and adding up all the sizes.
  // Note that we would also need to know the offsets.
  size_t all_query_points_allocation=0;
  int f_index_procs_others_offset[nprocs]; // offset in the all_query_points array
  int f_index_procs_self_offset[nprocs]; // offset in the query_outside array
  f_index_procs_others_offset[0]=0;
  f_index_procs_self_offset[0]=0;
  for (int proc=0;proc<nprocs;++proc){
    // The reason we multiply by COORD_DIM is that we have three coordinates per interpolation request
    all_query_points_allocation+=f_index_procs_others_sizes[proc]*COORD_DIM;
    if(proc>0){
      f_index_procs_others_offset[proc]=f_index_procs_others_offset[proc-1]+f_index_procs_others_sizes[proc-1];
      f_index_procs_self_offset[proc]=f_index_procs_self_offset[proc-1]+f_index_procs_self_sizes[proc-1];
    }
  }
  int total_query_points=all_query_points_allocation/COORD_DIM;
  Real * all_query_points=(Real*) malloc(all_query_points_allocation*sizeof(Real));
#ifdef VERBOSE2
  if(procid==0){
    std::cout<<"procid="<<procid<<std::endl;
    for (int proc=0;proc<nprocs;++proc)
      std::cout<<"proc= "<<proc<<" others_offset= "<<f_index_procs_others_offset[proc]<<" others_sizes= "<<f_index_procs_others_sizes[proc]<<std::endl;
    for (int proc=0;proc<nprocs;++proc)
      std::cout<<"proc= "<<proc<<" self_offset= "<<f_index_procs_self_offset[proc]<<" self_sizes= "<<f_index_procs_self_sizes[proc]<<std::endl;
  }
#endif

  MPI_Request * s_request= new MPI_Request[nprocs];
  MPI_Request * request= new MPI_Request[nprocs];

  // Now perform the allotall to send/recv query_points
  {
    int dst_r,dst_s;
    for (int i=0;i<nprocs;++i){
      dst_r=i;//(procid+i)%nprocs;
      dst_s=i;//(procid-i+nprocs)%nprocs;
      s_request[dst_s]=MPI_REQUEST_NULL;
      request[dst_r]=MPI_REQUEST_NULL;
      int roffset=f_index_procs_others_offset[dst_r]*COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
      int soffset=f_index_procs_self_offset[dst_s]*COORD_DIM;
      if(f_index_procs_others_sizes[dst_r]!=0)
        MPI_Irecv(&all_query_points[roffset],f_index_procs_others_sizes[dst_r]*COORD_DIM,MPI_T, dst_r,
            0, c_comm, &request[dst_r]);
      if(!query_outside[dst_s].empty())
        MPI_Isend(&query_outside[dst_s][0],f_index_procs_self_sizes[dst_s]*COORD_DIM,MPI_T,dst_s,
            0, c_comm, &s_request[dst_s]);
      //if(procid==1){
      //std::cout<<"soffset="<<soffset<<" roffset="<<roffset<<std::endl;
      //std::cout<<"f_index_procs_self_sizes[0]="<<f_index_procs_self_sizes[0]<<std::endl;
      //std::cout<<"f_index_procs_others_sizes[0]="<<f_index_procs_others_sizes[0]<<std::endl;
      //std::cout<<"q_outside["<<dst_s<<"]="<<query_outside[dst_s][0]<<std::endl;
      //}
    }
    // Wait for all the communication to finish
    MPI_Status ierr;
    for (int proc=0;proc<nprocs;++proc){
      if(request[proc]!=MPI_REQUEST_NULL)
        MPI_Wait(&request[proc], &ierr);
      if(s_request[proc]!=MPI_REQUEST_NULL)
        MPI_Wait(&s_request[proc], &ierr);
    }
  }


  // Now perform the interpolation on all query points including those that need to
  // be sent to other processors and store them into all_f_cubic
  Real* all_f_cubic_d;
  Real* all_query_points_d;
  cudaMalloc((void**)&all_f_cubic_d, total_query_points*sizeof(Real)*data_dof);
  cudaMalloc((void**)&all_query_points_d,all_query_points_allocation*sizeof(Real) );
  cudaMemcpy(all_query_points_d,all_query_points,all_query_points_allocation*sizeof(Real),cudaMemcpyHostToDevice);

  gpu_interp3_ghost_xyz_p(ghost_reg_grid_vals_d, data_dof, N_reg, isize,istart,total_query_points,g_size,all_query_points_d, all_f_cubic_d);
  Real* all_f_cubic=(Real*)malloc(total_query_points*sizeof(Real)*data_dof);
  cudaMemcpy(all_f_cubic,all_f_cubic_d, total_query_points*sizeof(Real)*data_dof ,cudaMemcpyDeviceToHost);


  // Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to
  // f_cubic_unordered.
  Real * f_cubic_unordered=(Real*) malloc(N_pts*sizeof(Real)*data_dof); // The reshuffled semi-final interpolated values are stored here
  {
    //PCOUT<<"total_query_points="<<total_query_points<<" N_pts="<<N_pts<<std::endl;
    int dst_r,dst_s;
    MPI_Datatype stype[nprocs],rtype[nprocs];
    for(int i=0;i<nprocs;++i){
      MPI_Type_vector(data_dof,f_index_procs_self_sizes[i],N_pts, MPI_T, &rtype[i]);
      MPI_Type_vector(data_dof,f_index_procs_others_sizes[i],total_query_points, MPI_T, &stype[i]);
      MPI_Type_commit(&stype[i]);
      MPI_Type_commit(&rtype[i]);
    }
    for (int i=0;i<nprocs;++i){
      dst_r=i;//(procid+i)%nprocs;
      dst_s=i;//(procid-i+nprocs)%nprocs;
      s_request[dst_s]=MPI_REQUEST_NULL;
      request[dst_r]=MPI_REQUEST_NULL;
      // Notice that this is the adjoint of the first comm part
      // because now you are sending others f and receiving your part of f
      int soffset=f_index_procs_others_offset[dst_r];
      int roffset=f_index_procs_self_offset[dst_s];
      //if(procid==0)
      //  std::cout<<"procid="<<procid<<" dst_s= "<<dst_s<<" soffset= "<<soffset<<" s_size="<<f_index_procs_others_sizes[dst_s]<<" dst_r= "<<dst_r<<" roffset="<<roffset<<" r_size="<<f_index_procs_self_sizes[dst_r]<<std::endl;
      //if(f_index_procs_self_sizes[dst_r]!=0)
      //  MPI_Irecv(&f_cubic_unordered[roffset],f_index_procs_self_sizes[dst_r],rtype, dst_r,
      //      0, c_comm, &request[dst_r]);
      //if(f_index_procs_others_sizes[dst_s]!=0)
      //  MPI_Isend(&all_f_cubic[soffset],f_index_procs_others_sizes[dst_s],stype,dst_s,
      //      0, c_comm, &s_request[dst_s]);
      //
      if(f_index_procs_self_sizes[dst_r]!=0)
        MPI_Irecv(&f_cubic_unordered[roffset],1,rtype[i], dst_r,
            0, c_comm, &request[dst_r]);
      if(f_index_procs_others_sizes[dst_s]!=0)
        MPI_Isend(&all_f_cubic[soffset],1,stype[i],dst_s,
            0, c_comm, &s_request[dst_s]);
    }
    MPI_Status ierr;
    for (int proc=0;proc<nprocs;++proc){
      if(request[proc]!=MPI_REQUEST_NULL)
        MPI_Wait(&request[proc], &ierr);
      if(s_request[proc]!=MPI_REQUEST_NULL)
        MPI_Wait(&s_request[proc], &ierr);
    }
    for(int i=0;i<nprocs;++i){
      MPI_Type_free(&stype[i]);
      MPI_Type_free(&rtype[i]);
    }
  }

  // Now copy back f_cubic_unordered to f_cubic in the correct f_index
  for(int dof=0;dof<data_dof;++dof){
    for(int proc=0;proc<nprocs;++proc){
      if(!f_index[proc].empty())
        for(int i=0;i<f_index[proc].size();++i){
          int ind=f_index[proc][i];
          //f_cubic[ind]=all_f_cubic[f_index_procs_others_offset[proc]+i];
          query_values[ind+dof*N_pts]=f_cubic_unordered[f_index_procs_self_offset[proc]+i+dof*N_pts];
        }
    }
  }

  free(query_points);
  free(all_query_points);
  free(all_f_cubic);
  cudaFree(all_f_cubic_d);
  cudaFree(all_query_points_d);
  free(f_cubic_unordered);
  delete(s_request);
  delete(request);
  //vector
  for(int proc=0;proc<nprocs;++proc)
  {
    std::vector<int>().swap(f_index[proc]);
    std::vector<Real>().swap(query_outside[proc]);
  }
  return;
}
