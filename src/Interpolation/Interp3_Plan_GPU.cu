#include <interp3_gpu_mpi.hpp>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime_api.h>
#include <cuda_helper.hpp>
#include <cmath>
#include <time.h>

#include <algorithm>

static void printGPUMemory(int rank) {
    if (rank == 0) {
      size_t free, used;
      cudaMemGetInfo(&free, &used);
      used -= free;
      std::string msg = "Used mem = " + std::to_string(used/1E9) + " GB, Free mem = " + std::to_string(free/1E9) + " GB\n";
      PetscPrintf(PETSC_COMM_WORLD, msg.c_str());
    }
}
//#define VERBOSE1

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

inline size_t get_max_query_allocation(int* isize, int neighbour_query_width) {
  return (isize[0]+2*neighbour_query_width)*isize[1]*isize[2]; 
}

struct is_equal {
    int id;
    is_equal(int comp_id) : id(comp_id) {};
    __host__ __device__ 
    bool operator()(const int &x) {
        return x == id;
    }
};

template <typename Iterator>
class strided_range
{
    public:

    typedef typename thrust::iterator_difference<Iterator>::type difference_type;

    struct stride_functor : public thrust::unary_function<difference_type,difference_type>
    {
        difference_type stride;

        stride_functor(difference_type stride)
            : stride(stride) {}

        __host__ __device__
        difference_type operator()(const difference_type& i) const
        { 
            return stride * i;
        }
    };

    typedef typename thrust::counting_iterator<difference_type>                   CountingIterator;
    typedef typename thrust::transform_iterator<stride_functor, CountingIterator> TransformIterator;
    typedef typename thrust::permutation_iterator<Iterator,TransformIterator>     PermutationIterator;

    // type of the strided_range iterator
    typedef PermutationIterator iterator;

    // construct strided_range for the range [first,last)
    strided_range(Iterator first, Iterator last, difference_type stride)
        : first(first), last(last), stride(stride) {}
   
    iterator begin(void) const
    {
        return PermutationIterator(first, TransformIterator(CountingIterator(0), stride_functor(stride)));
    }

    iterator end(void) const
    {
        return begin() + ((last - first) + (stride - 1)) / stride;
    }
    
    protected:
    Iterator first;
    Iterator last;
    difference_type stride;
};

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

#ifdef SORT_QUERIES
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
    delete[] trip;
  }
  return;
}
#endif


Interp3_Plan_GPU::Interp3_Plan_GPU (size_t g_alloc_max, bool cuda_aware) {
  this->g_alloc_max=g_alloc_max;
  this->allocate_baked=false;
  this->scatter_baked=false;
  this->cuda_aware = cuda_aware;
  this->procs_i_recv_from_size_ = 0;
  this->procs_i_send_to_size_ = 0;
  
  data_dofs = nullptr;
  stypes = nullptr;
  rtypes  = nullptr;

  all_query_points   = nullptr;
  all_f_cubic  = nullptr; 
  f_cubic_unordered   = nullptr;
  query_points  = nullptr;
    
  f_index_procs_others_offset  = nullptr;
  f_index_procs_self_offset    = nullptr;
  f_index_procs_self_sizes     = nullptr;
  f_index_procs_others_sizes   = nullptr;

  s_request   = nullptr;
  request  = nullptr;

  num_query_per_proc  = nullptr;;
  query_outside_offset  = nullptr;;
  f_index_offset  = nullptr;;

  which_proc  = nullptr;;
  query_points_x = nullptr;
  query_points_y = nullptr;
  query_points_z = nullptr;

  all_f_cubic_d  = nullptr;
  all_query_points_d = nullptr;
  f_cubic_unordered_d = nullptr;

  xq1 = nullptr;
  xq2 = nullptr;
  xq3 = nullptr;
}


void Interp3_Plan_GPU::allocate (int N_pts, int* data_dofs, int nplans)
{
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  reg::AllocateArrayOnce<int>(f_index_procs_others_offset, nprocs); // offset in the all_query_points array
  reg::AllocateArrayOnce<int>(f_index_procs_self_offset  , nprocs); // offset in the query_outside array
  reg::AllocateArrayOnce<int>(f_index_procs_self_sizes   , nprocs); // sizes of the number of interpolations that need to be sent to procs
  reg::AllocateArrayOnce<int>(f_index_procs_others_sizes , nprocs); // sizes of the number of interpolations that need to be received from procs

  reg::AllocateArrayOnce<MPI_Request>(s_request, nprocs);
  reg::AllocateArrayOnce<MPI_Request>(request, nprocs);
    
  query_outside = thrust::device_malloc<Real>(COORD_DIM*N_pts);
  f_index = thrust::device_malloc<int>(N_pts);
  reg::AllocateArrayOnce<int>(num_query_per_proc, nprocs);
  reg::AllocateArrayOnce<int>(query_outside_offset, nprocs);
  reg::AllocateArrayOnce<int>(f_index_offset, nprocs);
  
  // number of reuses of the plan with the same scatter points
  this->nplans_ = nplans; 
  reg::AllocateArrayOnce<int>(this->data_dofs, nplans);
  
  int max = 0;
  for(int i = 0; i < nplans_; ++i) {
    max = std::max(max, data_dofs[i]);
    this->data_dofs[i] = data_dofs[i];
  }
  this->data_dof_max = max;

  cudaMalloc((void**)&f_cubic_unordered_d, N_pts*data_dof_max*sizeof(Real));
  cudaMalloc((void**)&query_points_x, sizeof(Real)*N_pts);
  cudaMalloc((void**)&query_points_y, sizeof(Real)*N_pts);
  cudaMalloc((void**)&query_points_z, sizeof(Real)*N_pts);
  cudaMalloc((void**)&which_proc, sizeof(int)*N_pts);
    
  reg::AllocateArrayOnce<MPI_Datatype>(stypes, nprocs*nplans);
  reg::AllocateArrayOnce<MPI_Datatype>(rtypes, nprocs*nplans);
  
  this->allocate_baked=true;
}

Interp3_Plan_GPU::~Interp3_Plan_GPU ()
{
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  thrust::device_free(f_index);
  thrust::device_free(query_outside);
  
  if(this->allocate_baked){
    reg::FreeArray(f_index_procs_others_offset);
    reg::FreeArray(f_index_procs_self_offset  );
    reg::FreeArray(f_index_procs_self_sizes   );
    reg::FreeArray(f_index_procs_others_sizes );

    reg::FreeArray(s_request);
    reg::FreeArray(request);
    
    cudaFree(f_cubic_unordered_d);
    cudaFree(query_points_x);
    cudaFree(query_points_y);
    cudaFree(query_points_z);
    cudaFree(which_proc);
  }

  if(this->scatter_baked) {
    cudaFree(all_f_cubic_d);
    cudaFree(xq1);
    cudaFree(xq2);
    cudaFree(xq3);
    cudaFree(all_query_points_d);

		for (int ver = 0; ver < nplans_; ++ver) {
      for (int i = 0; i < nprocs; ++i) {
        MPI_Type_free(&stypes[i+ver*nprocs]);
        MPI_Type_free(&rtypes[i+ver*nprocs]);
      }
    }
    reg::FreeArray(stypes);
    reg::FreeArray(rtypes);
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
void Interp3_Plan_GPU::scatter( int* N_reg,  // global grid dimensions
                                int * isize, // local grid dimensions
                                int* istart, // local grid start indices
                                const int N_pts, // local grid point count
                                const int g_size, // ghost layer width
                                Real* query_points_in_x, // input query points
                                Real* query_points_in_y, // input query points
                                Real* query_points_in_z, // input query points
                                int* c_dims,  // process cartesian grid dimensions
                                MPI_Comm c_comm,  // MPI Comm
                                double * timings) 
{
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);
  bool verbose = false;

  if(this->allocate_baked==false){
    std::cout<<"ERROR Interp3_Plan_GPU Scatter called before calling allocate.\n";
    return;
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
    
    ZeitGeist_define(scatter_memcpy);
    ZeitGeist_tick(scatter_memcpy);
    cudaMemcpy(query_points_x, query_points_in_x, sizeof(Real)*N_pts, cudaMemcpyDeviceToDevice);
    cudaMemcpy(query_points_y, query_points_in_y, sizeof(Real)*N_pts, cudaMemcpyDeviceToDevice);
    cudaMemcpy(query_points_z, query_points_in_z, sizeof(Real)*N_pts, cudaMemcpyDeviceToDevice);
    ZeitGeist_tock(scatter_memcpy);


    ZeitGeist_define(scatter_create_mpi_buffer);
    ZeitGeist_tick(scatter_create_mpi_buffer);
    // Enforce periodicity // write kernel for this
    timings[3]+=-MPI_Wtime();
    enforcePeriodicity(query_points_x, query_points_y, query_points_z, h, N_pts);
    timings[3]+=+MPI_Wtime();
    ZeitGeist_tock(scatter_create_mpi_buffer);
    
    thrust::device_ptr<ScalarType> query_points_x_ptr = thrust::device_pointer_cast<ScalarType>(query_points_x);
    thrust::device_ptr<ScalarType> query_points_y_ptr = thrust::device_pointer_cast<ScalarType>(query_points_y);
    thrust::device_ptr<ScalarType> query_points_z_ptr = thrust::device_pointer_cast<ScalarType>(query_points_z);

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


    // This is needed for one-to-one correspondence with output f. This is becaues we are reshuffling
    // the data according to which processor it land onto, and we need to somehow keep the original
    // index to write the interpolation data back to the right location in the output.

    // This is necessary because when we want to compute dproc0 and dproc1 we have to divide by
    // the max isize. If the proc grid is unbalanced, the last proc's isize will be different
    // than others. With this approach we always use the right isize0 for all procs.
    int isize0=std::ceil(N_reg[0]*1./c_dims[0]);
    int isize1=std::ceil(N_reg[1]*1./c_dims[1]);
    
    // number of coordinates to be sent to each proc
    int coords_in_proc;
    typedef thrust::device_vector<ScalarType>::iterator Iterator;
    
    ZeitGeist_tick(scatter_create_mpi_buffer);
    timings[3]+=-MPI_Wtime();
    checkDomain(which_proc, query_points_x, query_points_y, query_points_z, iX0, iX1, h, N_pts, procid, isize0, isize1, c_dims[1]);
    ZeitGeist_tock(scatter_create_mpi_buffer);

    thrust::device_ptr<int> which_proc_ptr = thrust::device_pointer_cast<int>(which_proc);
    
    ZeitGeist_define(scatter_memalloc);
    
    // loop over all procs
    f_index_offset[0] = 0;
    query_outside_offset[0] = 0;
    for (int proc=0; proc<nprocs; ++proc) {
        // count how many points belong to proc, will be useful in memory allocation
        ZeitGeist_tick(scatter_create_mpi_buffer);
        get_count(which_proc, N_pts, proc, &coords_in_proc);
        ZeitGeist_tock(scatter_create_mpi_buffer);

#if defined(VERBOSE1) 
        PetscPrintf(PETSC_COMM_WORLD, "proc 0 sending %d points to proc %d\n", coords_in_proc, proc);
#endif
        
        num_query_per_proc[proc] = coords_in_proc;
        if (proc < nprocs-1) {
          f_index_offset[proc+1] = f_index_offset[proc] + coords_in_proc;
          query_outside_offset[proc+1] = query_outside_offset[proc] + coords_in_proc*COORD_DIM;
        }

        if (coords_in_proc > 0) {
            // get indices of coordinates which belong to this proc and store in f_index[proc]
            ZeitGeist_tick(scatter_create_mpi_buffer);
            thrust::copy_if(thrust::device, thrust::make_counting_iterator(0), thrust::make_counting_iterator(N_pts), which_proc_ptr, f_index+f_index_offset[proc], is_equal(proc));
            strided_range<Iterator> strided_x(query_outside+query_outside_offset[proc],   query_outside+query_outside_offset[proc]+coords_in_proc*COORD_DIM, COORD_DIM);
            strided_range<Iterator> strided_y(query_outside+query_outside_offset[proc]+1, query_outside+query_outside_offset[proc]+coords_in_proc*COORD_DIM, COORD_DIM);
            strided_range<Iterator> strided_z(query_outside+query_outside_offset[proc]+2, query_outside+query_outside_offset[proc]+coords_in_proc*COORD_DIM, COORD_DIM);
            thrust::copy_if(thrust::device, 
                            thrust::make_zip_iterator(thrust::make_tuple(query_points_x_ptr, query_points_y_ptr, query_points_z_ptr)), 
                            thrust::make_zip_iterator(thrust::make_tuple(query_points_x_ptr+N_pts, query_points_y_ptr+N_pts, query_points_z_ptr+N_pts)), 
                            which_proc_ptr, 
                            thrust::make_zip_iterator(thrust::make_tuple(strided_x.begin(), strided_y.begin(), strided_z.begin())), 
                            is_equal(proc));
            ZeitGeist_tock(scatter_create_mpi_buffer);
        }
    }
    timings[3]+=+MPI_Wtime();

    // Now we need to send the query_points that land onto other processor's domain.
    // This done using a sparse alltoallv.
    // Right now each process knows how much data to send to others, but does not know
    // how much data it should receive. This is a necessary information both for the MPI
    // command as well as memory allocation for received data.
    // So we first do an alltoall to get the f_index[proc].size from all processes.

    for (int proc=0;proc<nprocs;proc++) {
        f_index_procs_self_sizes[proc]=num_query_per_proc[proc];
    }
    ZeitGeist_define(scatter_comm_query_size);
    ZeitGeist_tick(scatter_comm_query_size);
    timings[0]+=-MPI_Wtime();
    MPI_Alltoall(f_index_procs_self_sizes,1, MPI_INT,
        f_index_procs_others_sizes,1, MPI_INT,
        c_comm);
    timings[0]+=+MPI_Wtime();
    ZeitGeist_tock(scatter_comm_query_size);

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
    
    //if (this->scatter_baked == false) {
    //  total_query_points_prev = total_query_points;
    //} else {
    //  total_query_points_prev = std::max(total_query_points_prev, total_query_points);
    //}
  
  // This if condition is to allow multiple calls to scatter fucntion with different query points
  // without having to create a new plan
  ZeitGeist_tick(scatter_memalloc);
  if(this->scatter_baked==true) {
    
    if (total_query_points > max_query_points_capacity) {
      // Modify the max memory estimate by increasing neighbour query width by 1 until the requrement is reached
      while (total_query_points > max_query_points_capacity) {
        neighbour_query_width++;
        max_query_points_capacity = get_max_query_allocation(isize, neighbour_query_width);
      }
      cudaFree(all_f_cubic_d);
      cudaMalloc((void**)&all_f_cubic_d, max_query_points_capacity*sizeof(Real)*data_dof_max);
      cudaFree(all_query_points_d);
      cudaMalloc((void**)&all_query_points_d, max_query_points_capacity*COORD_DIM*sizeof(Real) );
      cudaFree(xq1);
      cudaMalloc((void**)&xq1, max_query_points_capacity*sizeof(Real));
      cudaFree(xq2);
      cudaMalloc((void**)&xq2, max_query_points_capacity*sizeof(Real));
      cudaFree(xq3);
      cudaMalloc((void**)&xq3, max_query_points_capacity*sizeof(Real));
    }
  }
  else {
    // Make an estimate on the number of query points which the current process expects to receive
    // from neighbouring processes
    neighbour_query_width = g_size;
    max_query_points_capacity = get_max_query_allocation(isize, neighbour_query_width);
    while (total_query_points > max_query_points_capacity) {
      neighbour_query_width++;
      max_query_points_capacity = get_max_query_allocation(isize, neighbour_query_width);
    }

    cudaMalloc((void**)&all_f_cubic_d, max_query_points_capacity*sizeof(Real)*data_dof_max);
    cudaMalloc((void**)&all_query_points_d, max_query_points_capacity*COORD_DIM*sizeof(Real) );
    cudaMalloc((void**)&xq1, max_query_points_capacity*sizeof(Real));
    cudaMalloc((void**)&xq2, max_query_points_capacity*sizeof(Real));
    cudaMalloc((void**)&xq3, max_query_points_capacity*sizeof(Real));
  }
  ZeitGeist_tock(scatter_memalloc);  

    // Now perform the allotall to send/recv query_points
    ZeitGeist_define(scatter_comm_query_points_sendrcv);
    ZeitGeist_tick(scatter_comm_query_points_sendrcv);
    timings[0]+=-MPI_Wtime();
    int dst_r,dst_s;
    for (int i=0;i<nprocs;++i) {
      dst_r=i;
      dst_s=i;
      s_request[dst_s]=MPI_REQUEST_NULL;
      request[dst_r]=MPI_REQUEST_NULL;
      ScalarType* src_ptr = thrust::raw_pointer_cast(query_outside + query_outside_offset[dst_s]);
      int roffset=f_index_procs_others_offset[dst_r]*COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
      if (i != procid) {
        //int soffset=f_index_procs_self_offset[dst_s]*COORD_DIM;
        if(f_index_procs_others_sizes[dst_r]!=0)
          MPI_Irecv(&all_query_points_d[roffset], f_index_procs_others_sizes[dst_r]*COORD_DIM,MPI_T, dst_r, 0, c_comm, &request[dst_r]);

        //if(!query_outside[dst_s].empty())
        if(num_query_per_proc[dst_s] > 0)
          MPI_Isend(src_ptr, f_index_procs_self_sizes[dst_s]*COORD_DIM, MPI_T, dst_s, 0, c_comm, &s_request[dst_s]);

      } else {
        //if (!query_outside[dst_s].empty())
        if (num_query_per_proc[dst_s] > 0)
          reg::gencpy(&all_query_points_d[roffset], src_ptr, f_index_procs_self_sizes[dst_s]*COORD_DIM*sizeof(ScalarType));
      }
    }
    
    // Wait for all the communication to finish
    MPI_Status ierr;
    for (int proc=0;proc<nprocs;++proc) {
      if(request[proc]!=MPI_REQUEST_NULL)
        MPI_Wait(&request[proc], &ierr);
      if(s_request[proc]!=MPI_REQUEST_NULL)
        MPI_Wait(&s_request[proc], &ierr);
    }
    ZeitGeist_tock(scatter_comm_query_points_sendrcv);
    timings[0]+=+MPI_Wtime();
  
    // Now perform the interpolation on all query points including those that need to
    // be sent to other processors and store them into all_f_cubic
    //free(query_points);
  }

  for(int ver = 0; ver < nplans_; ++ver) {
    for (int i = 0; i < nprocs; ++i) {
      MPI_Type_vector(data_dofs[ver], f_index_procs_self_sizes[i], N_pts, MPI_T, &rtypes[i+ver*nprocs]);
      MPI_Type_vector(data_dofs[ver], f_index_procs_others_sizes[i], total_query_points, MPI_T, &stypes[i+ver*nprocs]);
      MPI_Type_commit(&stypes[i+ver*nprocs]);
      MPI_Type_commit(&rtypes[i+ver*nprocs]);
    }
  }

  int proc_coord[2];
  proc_coord[0] = static_cast<int>(istart[0]/isize[0]);
  proc_coord[1] = static_cast<int>(istart[1]/isize[1]);
  
  // transfer query points "all_query_points" from host to device
  ZeitGeist_define(scatter_query_points_normalize_kernel);
  ZeitGeist_tick(scatter_query_points_normalize_kernel);
  timings[3]+=-MPI_Wtime();
  normalizeQueryPoints(xq1, xq2, xq3, all_query_points_d, total_query_points, isize, N_reg, proc_coord, g_size);
  timings[3]+=+MPI_Wtime();
  ZeitGeist_tock(scatter_query_points_normalize_kernel);
  
  this->scatter_baked=true;
  return;
}


void Interp3_Plan_GPU::test_kernel(Real* f, int nq) {
  test(f, nq);
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

void Interp3_Plan_GPU::interpolate( Real* ghost_reg_grid_vals_d, // ghost padded regular grid values on GPU
                                    int* isize_g,              // size of the local grid (including ghost points)
                                    const IntType nlghost,         // number of local grid points (including ghost points) owned by process
                                    const IntType N_pts,           // number of local points owned by the process
                                    Real* query_values_d,      // interpolation result on GPU
                                    MPI_Comm c_comm,           // MPI communicator
                                    float *tmp1,               // temporary memory for interpolation prefilter
                                    float* tmp2,               // temporary memory for interpolation prefilter
                                    cudaTextureObject_t yi_tex,// texture object for interpolation
                                    int iporder,               // interpolation order
                                    ScalarType* interp_time,   // interpolation time
                                    int version)   
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

  bool verbose = false;

  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD, "max ghost value\n");
    print_max(ghost_reg_grid_vals_d, 3*nlghost);
    PetscPrintf(PETSC_COMM_WORLD, "max xq\n");
    print_max(xq1, total_query_points); 
  }

  ZeitGeist_define(interp_kernel);
  ZeitGeist_define(interp_comm_values_sendrcv);
  ZeitGeist_define(interp_values_copy_kernel);

  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD,"\ng_alloc_max = %zu", g_alloc_max);
    PetscPrintf(PETSC_COMM_WORLD,"\ndata_dof = %d", data_dofs[version]);
    PetscPrintf(PETSC_COMM_WORLD,"\nisize_g[0] = %d", isize_g[0]);
    PetscPrintf(PETSC_COMM_WORLD,"\nisize_g[1] = %d", isize_g[1]);
    PetscPrintf(PETSC_COMM_WORLD,"\nisize_g[2] = %d", isize_g[2]);
    PetscPrintf(PETSC_COMM_WORLD,"\nipoder = %d", iporder);
    PetscPrintf(PETSC_COMM_WORLD,"\ntotal_query_points = %d", total_query_points);
    PetscPrintf(PETSC_COMM_WORLD,"\nnlghost = %d\n", nlghost);
  }
  
  // compute the interpolation on the GPU
  ZeitGeist_tick(interp_kernel);
  double interp_kernel_time = -MPI_Wtime();
  if (data_dofs[version] == 3)
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
  ZeitGeist_tock(interp_kernel);
  interp_kernel_time += MPI_Wtime();

  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD, "max interpolated value\n");
    print_max(all_f_cubic_d, 3*total_query_points);
  }

#if defined(VERBOSE1) 
  int device;
  cudaGetDevice(&device);
  double global_interp_time = 0;
  MPI_Barrier(c_comm);
  MPI_Reduce(&interp_kernel_time, &global_interp_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  double all_interp_runtimes[nprocs];
  MPI_Gather(&interp_kernel_time, 1, MPI_DOUBLE, all_interp_runtimes, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "max =  %0.2E\t", global_interp_time);
  if (procid == 0) {
    PetscPrintf(PETSC_COMM_WORLD, "[");
    for (int i=0; i<nprocs; i++) {
      PetscPrintf(PETSC_COMM_WORLD, "%0.2E," , all_interp_runtimes[i]);
    }
    PetscPrintf(PETSC_COMM_WORLD, "]\n");
  }
#endif

  // Now we have to do an alltoall to distribute the interpolated data from all_f_cubic_d to f_cubic_unordered_d
  ZeitGeist_tick(interp_comm_values_sendrcv);
  int dst_r,dst_s;
  for (int i=0;i<nprocs;++i) {
    dst_r=i;
    dst_s=i;
    s_request[dst_s]=MPI_REQUEST_NULL;
    request[dst_r]=MPI_REQUEST_NULL;
    // Notice that this is the adjoint of the first comm part
    // because now you are sending others f and receiving your part of f
    int soffset=f_index_procs_others_offset[dst_r];
    int roffset=f_index_procs_self_offset[dst_s];
    if (true) {
      if(f_index_procs_self_sizes[dst_r]!=0)
        MPI_Irecv(&f_cubic_unordered_d[roffset],1,rtypes[dst_r+version*nprocs], dst_r,
            0, c_comm, &request[dst_r]); 
      if(f_index_procs_others_sizes[dst_s]!=0)
        MPI_Isend(&all_f_cubic_d[soffset],1,stypes[dst_s+version*nprocs],dst_s,
            0, c_comm, &s_request[dst_s]);
    } else {
      reg::gencpy(&f_cubic_unordered_d[roffset], &all_f_cubic_d[soffset], sizeof(ScalarType)*f_index_procs_self_sizes[i]);
    }
  }

  MPI_Status ierr;
  for (int proc=0;proc<nprocs;++proc){
    if(request[proc]!=MPI_REQUEST_NULL)
      MPI_Wait(&request[proc], &ierr);
    if(s_request[proc]!=MPI_REQUEST_NULL)
      MPI_Wait(&s_request[proc], &ierr);
  }
  
  ZeitGeist_tock(interp_comm_values_sendrcv);
  
  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD, "max f_unordered = ");
    print_max(f_cubic_unordered_d, 3*N_pts);
  }
  
  ZeitGeist_tick(interp_values_copy_kernel);
  int* f_index_ptr;
  // Now copy back f_cubic_unordered_d to query_values_d in the correct f_index
  for(int dof=0;dof<data_dofs[version]; ++dof) {
    for(int proc=0;proc<nprocs;++proc) {
      if(num_query_per_proc[proc] > 0) {
          f_index_ptr = thrust::raw_pointer_cast( f_index + f_index_offset[proc] );
          copyQueryValues(&query_values_d[dof*N_pts],
                          &f_cubic_unordered_d[f_index_procs_self_offset[proc]+dof*N_pts], 
                          f_index_ptr, 
                          num_query_per_proc[proc]);
      }
    }
  }
  ZeitGeist_tock(interp_values_copy_kernel);

  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD, "query_values ");
    print_max(query_values_d, 3*N_pts);
  }

/*
  ZeitGeist_tick(interp_comm_values_sendrcv);
  int dst_r,dst_s;
  for (int i = procs_i_send_to_size_-1; i >=0; --i) {
    if (i != procid) {
      dst_r = procs_i_send_to_[i];    
      request[dst_r] = MPI_REQUEST_NULL; 
      int roffset = f_index_procs_self_offset[dst_r];
      MPI_Irecv(&f_cubic_unordered_d[roffset], 1, rtypes[dst_r+version*nprocs], dst_r, 0,
          c_comm, &request[dst_r]);
    }
  }
  for (int i = 0; i < procs_i_recv_from_size_; ++i) {
    if (i != procid) {
      dst_s = procs_i_recv_from_[i];    
      s_request[dst_s] = MPI_REQUEST_NULL; 
      int soffset = f_index_procs_others_offset[dst_s];
      MPI_Isend(&all_f_cubic_d[soffset], 1, stypes[dst_s+version*nprocs], dst_s, 0, c_comm,
          &s_request[dst_s]);
    }
  }
  ZeitGeist_tock(interp_comm_values_sendrcv);
  
  for (int i = 0; i < procs_i_send_to_size_; ++i) {
    int proc = procs_i_send_to_[i];    
    ZeitGeist_tick(interp_comm_values_sendrcv);
    if (request[proc] != MPI_REQUEST_NULL) {
      MPI_Wait(&request[proc], MPI_STATUS_IGNORE);
    }
    ZeitGeist_tock(interp_comm_values_sendrcv);

    ZeitGeist_tick(interp_values_copy_kernel);
    for (int dof = 0; dof < data_dofs[version]; ++dof) {
      int *f_index_ptr = thrust::raw_pointer_cast( f_index + f_index_offset[proc] );
      copyQueryValues(&query_values_d[dof*N_pts],
                      &f_cubic_unordered_d[f_index_procs_self_offset[proc]+dof*N_pts], 
                      f_index_ptr, 
                      num_query_per_proc[proc]);
    }
    ZeitGeist_tock(interp_values_copy_kernel);
  }

  ZeitGeist_tick(interp_comm_values_sendrcv);
  for (int i = 0; i < procs_i_recv_from_size_; ++i) {
    int proc = procs_i_recv_from_[i];    
    if (s_request[proc] != MPI_REQUEST_NULL)
      MPI_Wait(&s_request[proc], MPI_STATUS_IGNORE);
  }
  ZeitGeist_tock(interp_comm_values_sendrcv);
*/
  return;
}
