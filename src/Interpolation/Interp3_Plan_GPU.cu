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

#define verbose false

//static void printGPUMemory(int rank) {
//    if (rank == 0) {
//      size_t free, used;
//      cudaMemGetInfo(&free, &used);
//      used -= free;
//      std::string msg = "Used mem = " + std::to_string(used/1E9) + " GB, Free mem = " + std::to_string(free/1E9) + " GB\n";
//      PetscPrintf(PETSC_COMM_WORLD, msg.c_str());
//    }
//}

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

inline size_t get_query_recv_max_capacity(IntType* isize, int neighbour_query_recv_width) {
  return (isize[0]+2*neighbour_query_recv_width)*isize[1]*isize[2]; 
}
inline size_t get_query_send_max_capacity(IntType* isize, int neighbour_query_send_width) {
  return (2*neighbour_query_send_width)*isize[1]*isize[2]; 
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

static void sort_queries(std::vector<ScalarType>* query_send,std::vector<int>* f_index,int* N_reg,ScalarType* h,MPI_Comm c_comm){

  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);
  for(int proc=0;proc<nprocs;++proc){
    int qsize=query_send[proc].size()/COORD_DIM;
    Trip_GPU* trip=new Trip_GPU[qsize];

    for(int i=0;i<qsize;++i){
      trip[i].x=query_send[proc][i*COORD_DIM+0];
      trip[i].y=query_send[proc][i*COORD_DIM+1];
      trip[i].z=query_send[proc][i*COORD_DIM+2];
      trip[i].ind=f_index[proc][i];
      trip[i].N[0]=N_reg[0];
      trip[i].N[1]=N_reg[1];
      trip[i].N[2]=N_reg[2];
      trip[i].h[0]=h[0];
      trip[i].h[1]=h[1];
      trip[i].h[2]=h[2];
    }

    std::sort(trip, trip + qsize, ValueCmp);

    query_send[proc].clear();
    f_index[proc].clear();

    for(int i=0;i<qsize;++i){
      query_send[proc].push_back(trip[i].x);
      query_send[proc].push_back(trip[i].y);
      query_send[proc].push_back(trip[i].z);
      f_index[proc].push_back(trip[i].ind);
    }
    delete[] trip;
  }
  return;
}
#endif


Interp3_Plan_GPU::Interp3_Plan_GPU (size_t g_alloc_max, bool cuda_aware) {
  this->g_alloc_max=g_alloc_max;
  this->cuda_aware = cuda_aware;
  
  data_dofs = nullptr;
  s_request   = nullptr;
  request  = nullptr;
  
  //query_send = nullptr;
  query_send_offset  = nullptr;

  which_proc  = nullptr;
  f_unordered = nullptr;
  all_f_capacity = 0;
  query_points_send_capacity = 0;
  all_f = nullptr;
  output_baked = false;
  query_baked = false;
  

  for (int i=0; i<numeq; i++) {
    Eq[i].allocate_baked = false;
    Eq[i].scatter_baked = false;
    Eq[i].stypes = nullptr;
    Eq[i].rtypes = nullptr;
    
    //Eq[i].f_index = nullptr;
    Eq[i].f_index_offset = nullptr;
    Eq[i].num_query_per_proc = nullptr;
  
    Eq[i].all_query_points = nullptr;
    Eq[i].xq1 = nullptr;
    Eq[i].xq2 = nullptr;
    Eq[i].xq3 = nullptr;

    Eq[i].f_index_procs_others_offset = nullptr;
    Eq[i].f_index_procs_self_offset = nullptr;
    Eq[i].f_index_procs_self_sizes = nullptr;
    Eq[i].f_index_procs_others_sizes = nullptr;

    Eq[i].total_query_points = 0;
    Eq[i].all_query_points_allocation = 0;
    Eq[i].neighbour_query_recv_width = 0;
    Eq[i].query_points_recv_capacity = 0;
  }
}

/****************************************************************************************************
* Allocate Memory
****************************************************************************************************/
void Interp3_Plan_GPU::allocate (int N_pts, int* data_dofs, int nplans, int gsize, IntType *isize_g)
{
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  // neighbouring procs left to right (excluding self) for slab decomp
  proc_neighbours.push_back(procid-1);
  if (proc_neighbours[0] < 0) proc_neighbours[0] += nprocs;
  proc_neighbours.push_back(procid);
  proc_neighbours.push_back((procid+1)%nprocs);

  // number of reuses of the plan with the same scatter points
  this->nplans_ = nplans; 
  reg::AllocateArrayOnce<int>(this->data_dofs, nplans);
  
  int max = 0;
  for(int i = 0; i < nplans_; ++i) {
    max = std::max(max, data_dofs[i]);
    this->data_dofs[i] = data_dofs[i];
  }
  this->data_dof_max = max;
  
  reg::AllocateArrayOnce<int>(query_send_offset, nprocs);
  reg::AllocateArrayOnce<MPI_Request>(s_request, nprocs);
  reg::AllocateArrayOnce<MPI_Request>(request, nprocs);
  reg::AllocateMemoryOnce(which_proc, sizeof(short)*N_pts);
  
  for (int i=0; i<numeq; i++) {
    reg::AllocateArrayOnce<int>(Eq[i].f_index_procs_others_offset, nprocs); // offset in the all_query_points array
    reg::AllocateArrayOnce<int>(Eq[i].f_index_procs_self_offset  , nprocs); // offset in the query_send array
    reg::AllocateArrayOnce<int>(Eq[i].f_index_procs_self_sizes   , nprocs); // sizes of the number of interpolations that need to be sent to procs
    reg::AllocateArrayOnce<int>(Eq[i].f_index_procs_others_sizes , nprocs); // sizes of the number of interpolations that need to be received from procs

    
    Eq[i].f_index = thrust::device_malloc<int>(N_pts);
    reg::AllocateArrayOnce<int>(Eq[i].f_index_offset, nprocs);
    reg::AllocateArrayOnce<int>(Eq[i].num_query_per_proc, nprocs);
  
    reg::AllocateArrayOnce<MPI_Datatype>(Eq[i].stypes, nprocs*nplans);
    reg::AllocateArrayOnce<MPI_Datatype>(Eq[i].rtypes, nprocs*nplans);
    
    Eq[i].allocate_baked=true;
    
    Eq[i].neighbour_query_recv_width = gsize;
    Eq[i].query_points_recv_capacity = get_query_recv_max_capacity(isize_g, Eq[i].neighbour_query_recv_width);
    
    reg::AllocateMemoryOnce(Eq[i].all_query_points, Eq[i].query_points_recv_capacity*COORD_DIM*sizeof(ScalarType) );
#ifdef BLOCK_COORDINATES
    reg::AllocateMemoryOnce(Eq[i].xq1, Eq[i].query_points_recv_capacity*sizeof(ScalarType));
    reg::AllocateMemoryOnce(Eq[i].xq2, Eq[i].query_points_recv_capacity*sizeof(ScalarType));
    reg::AllocateMemoryOnce(Eq[i].xq3, Eq[i].query_points_recv_capacity*sizeof(ScalarType));
#endif
    Eq[i].scatter_baked=true;
  }
  
  all_f_capacity = get_query_recv_max_capacity(isize_g, gsize);
  reg::AllocateMemoryOnce(all_f, all_f_capacity*sizeof(ScalarType));
  output_baked = true;
  
  neighbour_query_send_width = gsize;
  query_points_send_capacity = get_query_send_max_capacity(isize_g, neighbour_query_send_width);
  query_send = thrust::device_malloc<ScalarType>(COORD_DIM*query_points_send_capacity);
  query_baked = true;
}

/************************************************************
 * Destructor 
 ************************************************************/
Interp3_Plan_GPU::~Interp3_Plan_GPU ()
{
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  reg::FreeArray(this->data_dofs);
  if (!query_send.get()) {
    thrust::device_free(query_send);
    //query_send = nullptr;
  }
  reg::FreeArray(query_send_offset);
  reg::FreeArray(s_request);
  reg::FreeArray(request);
  reg::FreeMemory(which_proc);
  reg::FreeMemory(all_f);
  reg::FreeMemory(f_unordered);
  
  for (int i=0; i<numeq; i++) {
    reg::FreeArray(Eq[i].f_index_procs_others_offset);
    reg::FreeArray(Eq[i].f_index_procs_self_offset  );
    reg::FreeArray(Eq[i].f_index_procs_self_sizes   );
    reg::FreeArray(Eq[i].f_index_procs_others_sizes );
	if (!(Eq[i].f_index).get()) {
      thrust::device_free(Eq[i].f_index);
//      Eq[i].f_index = nullptr;
    }
    reg::FreeArray(Eq[i].f_index_offset);
    reg::FreeArray(Eq[i].num_query_per_proc);
    reg::FreeMemory(Eq[i].all_query_points);
    reg::FreeMemory(Eq[i].xq1);
    reg::FreeMemory(Eq[i].xq2);
    reg::FreeMemory(Eq[i].xq3);
		for (int ver = 0; ver < nplans_; ++ver) {
      for (int j = 0; j < nprocs; ++j) {
        MPI_Type_free(&Eq[i].stypes[j+ver*nprocs]);
        MPI_Type_free(&Eq[i].rtypes[j+ver*nprocs]);
      }
    }
    reg::FreeArray(Eq[i].stypes);
    reg::FreeArray(Eq[i].rtypes);
  }
  return;
}

void rescale_xyz(const int g_size,  int* N_reg, int* N_reg_g, int* istart, const int N_pts, ScalarType* query_points);


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
void Interp3_Plan_GPU::scatter( IntType* N_reg,  // global grid dimensions
                                IntType * isize, // local grid dimensions
                                IntType* istart, // local grid start indices
                                const int N_pts, // local grid point count
                                const int g_size, // ghost layer width
                                ScalarType* query_points_x, // input query points, will be overwritten
                                ScalarType* query_points_y, // input query points, will be overwritten
                                ScalarType* query_points_z, // input query points, will be overwritten
                                int* c_dims,  // process cartesian grid dimensions
                                MPI_Comm c_comm,  // MPI Comm
                                double * timings,
                                std::string flag) 
{
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);
  
  int id = 0;
  if (flag.compare("state") == 0) {
    id = 0;
  } else if (flag.compare("adjoint") == 0) {
    id = 1;
  } else { 
    reg::ThrowError("wrong flag");
  }

  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD, "max xq = ");
    print_max(query_points_x, N_pts);
  }



  if(Eq[id].allocate_baked==false){
    std::cout<<"ERROR Interp3_Plan_GPU Scatter called before calling allocate.\n";
    return;
  }
  
  Eq[id].all_query_points_allocation=0;
 {
    // original grid size along each axis
    ScalarType h[3]; 
    h[0]=1./N_reg[0];
    h[1]=1./N_reg[1];
    h[2]=1./N_reg[2];
    
    // Enforce periodicity // write kernel for this
    enforcePeriodicity(query_points_x, query_points_y, query_points_z, h, N_pts);
    
    // Compute the start and end coordinates that this processor owns
    ScalarType iX0[3],iX1[3];
    for (int j=0;j<3;j++) {
      iX0[j]=istart[j]*h[j];
      iX1[j]=iX0[j]+(isize[j]-1)*h[j];
    }

    // Now march through the query points and split them into nprocs parts.
    // These are stored in query_send which is an array of vectors of size nprocs.
    // That is query_send[i] is a vector that contains the query points that need to
    // be sent to process i. Obviously for the case of query_send[procid], we do not
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
    
    timings[3]+=-MPI_Wtime();
    checkDomain(which_proc, query_points_x, query_points_y, query_points_z, iX0, iX1, h, N_pts, procid, isize0, isize1, c_dims[1]);
    
    thrust::device_ptr<ScalarType> query_points_x_ptr = thrust::device_pointer_cast<ScalarType>(query_points_x);
    thrust::device_ptr<ScalarType> query_points_y_ptr = thrust::device_pointer_cast<ScalarType>(query_points_y);
    thrust::device_ptr<ScalarType> query_points_z_ptr = thrust::device_pointer_cast<ScalarType>(query_points_z);

    thrust::device_ptr<short> which_proc_ptr = thrust::device_pointer_cast<short>(which_proc);
    
    // loop over all procs
    Eq[id].f_index_offset[0] = 0;
    query_send_offset[0] = 0;
    
    size_t query_send_size = 0;
    // specific to slab decomposition
    for (int proc=0; proc<nprocs; ++proc) {
        // count how many points belong to proc
        get_count(which_proc, N_pts, proc, &coords_in_proc);
        
      if (verbose) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "proc %d sending %d points to proc %d\n", procid, coords_in_proc, proc);
        PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
      }
        
      Eq[id].num_query_per_proc[proc] = coords_in_proc;
      if (proc != procid) {
        query_send_size += Eq[id].num_query_per_proc[proc];
      }
      if (proc < nprocs-1) {
        Eq[id].f_index_offset[proc+1] = Eq[id].f_index_offset[proc] + coords_in_proc;
        // only set query offset for neighboring procs
        if (proc != procid) {  
          query_send_offset[proc+1] = query_send_offset[proc] + coords_in_proc*COORD_DIM;
        } else {
          query_send_offset[proc+1] = query_send_offset[proc];
        }
      }

      if (coords_in_proc > 0) {
          // get indices of coordinates which belong to this proc and store in Eq[id].f_index[proc]
          thrust::copy_if(thrust::device, thrust::make_counting_iterator(0), thrust::make_counting_iterator(N_pts), which_proc_ptr, Eq[id].f_index+Eq[id].f_index_offset[proc], is_equal(proc));
      }
    }
    timings[3]+=+MPI_Wtime();
    

    // Now we need to send the query_points that land onto other processor's domain.
    // This done using a sparse alltoallv.
    // Right now each process knows how much data to send to others, but does not know
    // how much data it should receive. This is a necessary information both for the MPI
    // command as well as memory allocation for received data.
    // So we first do an alltoall to get the Eq[id].f_index[proc].size from all processes.

    for (int proc=0;proc<nprocs;proc++) {
        Eq[id].f_index_procs_self_sizes[proc]=Eq[id].num_query_per_proc[proc];
    }
    timings[0]+=-MPI_Wtime();
    MPI_Alltoall(Eq[id].f_index_procs_self_sizes,1, MPI_INT,
        Eq[id].f_index_procs_others_sizes,1, MPI_INT,
        c_comm);
    timings[0]+=+MPI_Wtime();

    // Now we need to allocate memory for the receiving buffer of all query
    // points including ours. This is simply done by looping through
    // Eq[id].f_index_procs_others_sizes and adding up all the sizes.
    // Note that we would also need to know the offsets.
    Eq[id].f_index_procs_others_offset[0]=0;
    Eq[id].f_index_procs_self_offset[0]=0;
    for (int proc=0;proc<nprocs;++proc) {
      // The reason we multiply by COORD_DIM is that we have three coordinates per interpolation request
      Eq[id].all_query_points_allocation+=Eq[id].f_index_procs_others_sizes[proc]*COORD_DIM;
      if(proc>0) {
        Eq[id].f_index_procs_others_offset[proc]=Eq[id].f_index_procs_others_offset[proc-1]+Eq[id].f_index_procs_others_sizes[proc-1];
        Eq[id].f_index_procs_self_offset[proc]=Eq[id].f_index_procs_self_offset[proc-1]+Eq[id].f_index_procs_self_sizes[proc-1];
      }
    }
    Eq[id].total_query_points=Eq[id].all_query_points_allocation/COORD_DIM;
    
  // This if condition is to allow multiple calls to scatter fucntion with different query points
  // without having to create a new plan
  if(Eq[id].scatter_baked && Eq[id].total_query_points > Eq[id].query_points_recv_capacity) {
    // Modify the max memory estimate by increasing neighbour query width by 1 until the requrement is reached
    while (Eq[id].total_query_points > Eq[id].query_points_recv_capacity) {
      Eq[id].neighbour_query_recv_width++;
      Eq[id].query_points_recv_capacity = get_query_recv_max_capacity(isize, Eq[id].neighbour_query_recv_width);
    }

    reg::FreeMemory(Eq[id].all_query_points);
    reg::AllocateMemoryOnce(Eq[id].all_query_points, Eq[id].query_points_recv_capacity*COORD_DIM*sizeof(ScalarType) );

#ifdef BLOCK_COORDINATES
    reg::FreeMemory(Eq[id].xq1);
    reg::AllocateMemoryOnce(Eq[id].xq1, Eq[id].query_points_recv_capacity*sizeof(ScalarType));

    reg::FreeMemory(Eq[id].xq2);
    reg::AllocateMemoryOnce(Eq[id].xq2, Eq[id].query_points_recv_capacity*sizeof(ScalarType));

    reg::FreeMemory(Eq[id].xq3);
    reg::AllocateMemoryOnce(Eq[id].xq3, Eq[id].query_points_recv_capacity*sizeof(ScalarType));
#endif
  }
  
  // this allocates a single vector field for storing output of both adjoint and scalar 
  if (output_baked && Eq[id].query_points_recv_capacity > all_f_capacity) {
    // update maximum
    all_f_capacity = Eq[id].query_points_recv_capacity;
    // free and reallocate
    reg::FreeMemory(all_f);
    reg::AllocateMemoryOnce(all_f, all_f_capacity*sizeof(ScalarType));
  }
  
  // Allocate MPI send buffer here
  if (query_baked && query_send_size > query_points_send_capacity) {
    while (query_send_size > query_points_send_capacity) {
      neighbour_query_send_width++;
      query_points_send_capacity = get_query_send_max_capacity(isize, neighbour_query_send_width);
    }
    thrust::device_free(query_send);
    query_send = thrust::device_malloc<ScalarType>(COORD_DIM*query_points_send_capacity);
  }

  for (int proc=0; proc<nprocs; ++proc) {
    coords_in_proc = Eq[id].num_query_per_proc[proc];
    // only need to copy to query outside buffer for neighbouring procs
    if (coords_in_proc > 0 && proc != procid) {
      strided_range<Iterator> strided_x(query_send+query_send_offset[proc],   query_send+query_send_offset[proc]+coords_in_proc*COORD_DIM, COORD_DIM);
      strided_range<Iterator> strided_y(query_send+query_send_offset[proc]+1, query_send+query_send_offset[proc]+coords_in_proc*COORD_DIM, COORD_DIM);
      strided_range<Iterator> strided_z(query_send+query_send_offset[proc]+2, query_send+query_send_offset[proc]+coords_in_proc*COORD_DIM, COORD_DIM);
      thrust::copy_if(thrust::device, 
                      thrust::make_zip_iterator(thrust::make_tuple(query_points_x_ptr, query_points_y_ptr, query_points_z_ptr)), 
                      thrust::make_zip_iterator(thrust::make_tuple(query_points_x_ptr+N_pts, query_points_y_ptr+N_pts, query_points_z_ptr+N_pts)), 
                      which_proc_ptr, 
                      thrust::make_zip_iterator(thrust::make_tuple(strided_x.begin(), strided_y.begin(), strided_z.begin())), 
                      is_equal(proc));
    } 
    if (coords_in_proc > 0 && proc == procid) {
      int roffset = Eq[id].f_index_procs_others_offset[proc]*COORD_DIM;
      thrust::device_ptr<ScalarType> all_query_points_ptr = thrust::device_pointer_cast<ScalarType>(&Eq[id].all_query_points[roffset]);
      strided_range<Iterator> strided_x(all_query_points_ptr+0, all_query_points_ptr+coords_in_proc*COORD_DIM, COORD_DIM);
      strided_range<Iterator> strided_y(all_query_points_ptr+1, all_query_points_ptr+coords_in_proc*COORD_DIM, COORD_DIM);
      strided_range<Iterator> strided_z(all_query_points_ptr+2, all_query_points_ptr+coords_in_proc*COORD_DIM, COORD_DIM);
      thrust::copy_if(thrust::device, 
                      thrust::make_zip_iterator(thrust::make_tuple(query_points_x_ptr, query_points_y_ptr, query_points_z_ptr)), 
                      thrust::make_zip_iterator(thrust::make_tuple(query_points_x_ptr+N_pts, query_points_y_ptr+N_pts, query_points_z_ptr+N_pts)), 
                      which_proc_ptr, 
                      thrust::make_zip_iterator(thrust::make_tuple(strided_x.begin(), strided_y.begin(), strided_z.begin())), 
                      is_equal(proc));
    }
  }

    // Now perform the allotall to send/recv query_points
    timings[0]+=-MPI_Wtime();
    int dst_r,dst_s;
    for (int i=0;i<nprocs;++i) {
      dst_r=i;
      dst_s=i;
      s_request[dst_s]=MPI_REQUEST_NULL;
      request[dst_r]=MPI_REQUEST_NULL;
      ScalarType* src_ptr = thrust::raw_pointer_cast(query_send + query_send_offset[dst_s]);
      int roffset=Eq[id].f_index_procs_others_offset[dst_r]*COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
      if (i != procid) {
        if(Eq[id].f_index_procs_others_sizes[dst_r]!=0)
          MPI_Irecv(&Eq[id].all_query_points[roffset], Eq[id].f_index_procs_others_sizes[dst_r]*COORD_DIM,MPI_T, dst_r, 0, c_comm, &request[dst_r]);

        if(Eq[id].num_query_per_proc[dst_s] > 0)
          MPI_Isend(src_ptr, Eq[id].f_index_procs_self_sizes[dst_s]*COORD_DIM, MPI_T, dst_s, 0, c_comm, &s_request[dst_s]);

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
    timings[0]+=+MPI_Wtime();
  }

  for(int ver = 0; ver < nplans_; ++ver) {
    for (int i = 0; i < nprocs; ++i) {
      MPI_Type_vector(data_dofs[ver], Eq[id].f_index_procs_self_sizes[i], N_pts, MPI_T, &Eq[id].rtypes[i+ver*nprocs]);
      MPI_Type_vector(data_dofs[ver], Eq[id].f_index_procs_others_sizes[i], Eq[id].total_query_points, MPI_T, &Eq[id].stypes[i+ver*nprocs]);
      MPI_Type_commit(&Eq[id].stypes[i+ver*nprocs]);
      MPI_Type_commit(&Eq[id].rtypes[i+ver*nprocs]);
    }
  }

  int proc_coord[2];
  proc_coord[0] = static_cast<int>(istart[0]/isize[0]);
  proc_coord[1] = static_cast<int>(istart[1]/isize[1]);
  
  // transfer query points "all_query_points" from host to device
  timings[3]+=-MPI_Wtime();
  normalizeQueryPoints(Eq[id].xq1, Eq[id].xq2, Eq[id].xq3, Eq[id].all_query_points, Eq[id].total_query_points, isize, N_reg, proc_coord, g_size);
  timings[3]+=+MPI_Wtime();
  
  Eq[id].scatter_baked=true;
  return;
}


void Interp3_Plan_GPU::test_kernel(ScalarType* f, int nq) {
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

void Interp3_Plan_GPU::interpolate( ScalarType* f_ghost, // ghost padded regular grid values on GPU
                                    IntType* isize_g,              // size of the local grid (including ghost points)
                                    const IntType nlghost,         // number of local grid points (including ghost points) owned by process
                                    const IntType N_pts,           // number of local points owned by the process
                                    ScalarType** query_values,      // interpolation result on GPU
                                    MPI_Comm c_comm,           // MPI communicator
                                    float *tmp1,               // temporary memory for interpolation prefilter
                                    float* tmp2,               // temporary memory for interpolation prefilter
                                    cudaTextureObject_t yi_tex,// texture object for interpolation
                                    int iporder,               // interpolation order
                                    ScalarType* interp_time,   // interpolation time
                                    int version,
                                    std::string flag)   
{
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);
  
  int id = 0;
  if (flag.compare("state") == 0) {
    id = 0;
  } else if (flag.compare("adjoint") == 0) {
    id = 1;
  } else { 
    reg::ThrowError("wrong flag");
  }
  
  
  if(Eq[id].allocate_baked==false){
    std::cout<<"ERROR Interp3_Plan_GPU interpolate called before calling allocate.\n";
    return;
  }
  if(Eq[id].scatter_baked==false){
    std::cout<<"ERROR Interp3_Plan_GPU interpolate called before calling scatter.\n";
    return;
  }


  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD,"\ng_alloc_max = %zu", g_alloc_max);
    PetscPrintf(PETSC_COMM_WORLD,"\ndata_dof = %d", data_dofs[version]);
    PetscPrintf(PETSC_COMM_WORLD,"\nisize_g[0] = %d", isize_g[0]);
    PetscPrintf(PETSC_COMM_WORLD,"\nisize_g[1] = %d", isize_g[1]);
    PetscPrintf(PETSC_COMM_WORLD,"\nisize_g[2] = %d", isize_g[2]);
    PetscPrintf(PETSC_COMM_WORLD,"\nipoder = %d", iporder);
    PetscPrintf(PETSC_COMM_WORLD,"\nEq[id].total_query_points = %d", Eq[id].total_query_points);
    PetscPrintf(PETSC_COMM_WORLD,"\nnlghost = %d\n", nlghost);
    PetscPrintf(PETSC_COMM_WORLD, "max ghost ");
    print_norm(f_ghost, data_dofs[version]*nlghost);
    PetscPrintf(PETSC_COMM_WORLD, "max all query ");
    print_norm(Eq[id].all_query_points, COORD_DIM*Eq[id].total_query_points);
    PetscPrintf(PETSC_COMM_WORLD, "max xq ");
    print_norm(Eq[id].xq1, Eq[id].total_query_points); 
  }

  // change this to have strided or block access to query points
#ifdef BLOCK_COORDINATES
  const ScalarType* xq[3] = {Eq[id].xq1, Eq[id].xq2, Eq[id].xq3};
  reg::Assert(Eq[id].xq1 != nullptr, "Eq[id].xq1 is nullptr");
  reg::Assert(Eq[id].xq2 != nullptr, "Eq[id].xq2 is nullptr");
  reg::Assert(Eq[id].xq3 != nullptr, "Eq[id].xq3 is nullptr");
#else
  const ScalarType* xq[3] = {Eq[id].all_query_points};
#endif

  // compute the interpolation on the GPU
  double interp_kernel_time = -MPI_Wtime();
  if (data_dofs[version] == 3) 
    gpuInterpVec3D(&f_ghost[0*nlghost], 
                   &f_ghost[1*nlghost], 
                   &f_ghost[2*nlghost], 
                   xq, 
                   &all_f[0*Eq[id].total_query_points], 
                   &all_f[1*Eq[id].total_query_points], 
                   &all_f[2*Eq[id].total_query_points], 
                   tmp1, tmp2, isize_g, static_cast<long int>(Eq[id].total_query_points), yi_tex, iporder, interp_time);
  else 
    gpuInterp3D(f_ghost, 
                xq, 
                all_f, 
                tmp1, tmp2, isize_g, static_cast<long int>(Eq[id].total_query_points), yi_tex, 
                iporder, interp_time);
  interp_kernel_time += MPI_Wtime();
  
  if (verbose) {
    reg::DbgMsgCall("interpolation kernel executed");
    PetscPrintf(PETSC_COMM_WORLD, "max interpolated ");
    print_norm(all_f, data_dofs[version]*Eq[id].total_query_points);
  }
  
  // f_unordered needs to be of size nl==N_pts, and f_ghost has nlghost>nl memory, reuse memory
  ScalarType* f_unordered_ptr = f_ghost;

  // Now we have to do an alltoall to distribute the interpolated data from Eq[id].all_f to f_unordered
  int dst_r,dst_s;
  for (int i=0;i<nprocs;++i) {
    dst_r=i;
    dst_s=i;
    s_request[dst_s]=MPI_REQUEST_NULL;
    request[dst_r]=MPI_REQUEST_NULL;
    // Notice that this is the adjoint of the first comm part
    // because now you are sending others f and receiving your part of f
    int soffset=Eq[id].f_index_procs_others_offset[dst_r];
    int roffset=Eq[id].f_index_procs_self_offset[dst_s];
    if (i != procid) {
      if(Eq[id].f_index_procs_self_sizes[dst_r]!=0)
        MPI_Irecv(&f_unordered_ptr[roffset],1,Eq[id].rtypes[dst_r+version*nprocs], dst_r,
            0, c_comm, &request[dst_r]); 
      if(Eq[id].f_index_procs_others_sizes[dst_s]!=0)
        MPI_Isend(&all_f[soffset],1,Eq[id].stypes[dst_s+version*nprocs],dst_s,
            0, c_comm, &s_request[dst_s]);
    } else {
      for (int dof=0; dof<data_dofs[version]; ++dof) {
        // when sending to itself i.e. i==procid, Eq[id].f_index_procs_self_sizes[i] == Eq[id].f_index_procs_others_sizes[i]
        cudaMemcpy(&f_unordered_ptr[roffset+dof*N_pts], &all_f[soffset+dof*Eq[id].total_query_points], sizeof(ScalarType)*Eq[id].f_index_procs_self_sizes[i], cudaMemcpyDeviceToDevice);
      }
    }
  }

  MPI_Status ierr;
  for (int proc=0;proc<nprocs;++proc){
    if(request[proc]!=MPI_REQUEST_NULL)
      MPI_Wait(&request[proc], &ierr);
    if(s_request[proc]!=MPI_REQUEST_NULL)
      MPI_Wait(&s_request[proc], &ierr);
  }
  
  
  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD, "max f_unordered = ");
    //print_norm(f_unordered, data_dofs[version]*N_pts);
    print_norm(f_unordered_ptr, data_dofs[version]*N_pts);
  }
  
  int* f_index_ptr;
  // Now copy back f_unordered to query_values in the correct f_index
  for(int dof=0;dof<data_dofs[version]; ++dof) {
    for(int proc=0;proc<nprocs;++proc) {
      if(Eq[id].num_query_per_proc[proc] > 0) {
          f_index_ptr = thrust::raw_pointer_cast( Eq[id].f_index + Eq[id].f_index_offset[proc] );
          copyQueryValues(&query_values[dof][0],
                          &f_unordered_ptr[Eq[id].f_index_procs_self_offset[proc]+dof*N_pts], 
                          f_index_ptr, 
                          Eq[id].num_query_per_proc[proc]);
      }
    }
  }

  if (verbose) {
    PetscPrintf(PETSC_COMM_WORLD, "query_values ");
    for(int dof=0;dof<data_dofs[version]; ++dof) {
      print_norm(query_values[dof], N_pts);
    }
  }

  return;
}
