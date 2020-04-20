#pragma once

#include <cufft.h>
#include <cuda.h>
#include <mpi.h>
#include <vector>

template<typename T>
class MPIcuFFT {
public:
  MPIcuFFT (MPI_Comm comm=MPI_COMM_WORLD, bool mpi_cuda_aware=false);
  virtual ~MPIcuFFT ();
  
  virtual void initFFT(size_t nx, size_t ny, size_t nz, bool allocate=true);
  virtual void setWorkArea(void *device=nullptr, void *host=nullptr);
  
  virtual inline size_t getDomainSize() const { return domainsize; };
  virtual inline size_t getWorkSizeDevice() const { return worksize_d; };
  virtual inline size_t getWorkSizeHost() const { return worksize_h; };
  
  virtual inline void* getWorkAreaDevice() const { return workarea_d; };
  virtual inline void* getWorkAreaHost() const { return workarea_h; };
  
  virtual void execR2C(void *out, const void *in);
  virtual void execC2R(void *out, const void *in);
  
  virtual inline void getInSize(size_t *isize) { isize[0] = isizex[pidx]; isize[1] = isizey; isize[2] = isizez; };
  virtual inline void getInStart(size_t *istart) { istart[0] = istartx[pidx]; istart[1] = 0; istart[2] = 0; };
  virtual inline void getOutSize(size_t *osize) { osize[0] = osizex; osize[1] = osizey[pidx]; osize[2] = osizez; };
  virtual inline void getOutStart(size_t *ostart) { ostart[0] = 0; ostart[1] = ostarty[pidx]; ostart[2] = 0; };
  
  virtual inline void restrictTo(void *coarse, const void *fine, MPIcuFFT<T>* fft_coarse) {
    changeSize(coarse, fine, fft_coarse, Restrict);
  }
  virtual inline void prolongFrom(void *fine, const void *coarse, MPIcuFFT<T>* fft_coarse) {
    changeSize(fine, coarse, fft_coarse, Prolong);
  }
  virtual inline void prolongFromMerge(void *fine, const void *coarse, MPIcuFFT<T>* fft_coarse) {
    changeSize(fine, coarse, fft_coarse, Prolong, false);
  }
  
  static size_t getSizes(const size_t* nx, 
                         size_t *isize, size_t *istart, 
                         size_t *osize, size_t *ostart,
                         MPI_Comm comm=MPI_COMM_WORLD);
protected:
  enum changeType_e {Prolong, Restrict};
  virtual void changeSize(void *out, const void *in, MPIcuFFT<T>* fft_coarse, changeType_e direction, bool clear_output=true);
  
  cufftHandle planR2C;
  cufftHandle planC2R;
  cufftHandle planC2C;

  std::vector<MPI_Request> send_req;
  std::vector<MPI_Request> recv_req;
  
  MPI_Comm comm;
  int pidx, pcnt;
  
  std::vector<int> comm_order;
  
  std::vector<size_t> isizex;
  std::vector<size_t> istartx;
  std::vector<size_t> osizey;
  std::vector<size_t> ostarty;
  
  size_t isizey, isizez;
  size_t osizex, osizez;
  
  size_t domainsize;
  size_t fft_worksize;
  
  size_t worksize_d;
  size_t worksize_h;
  
  void* workarea_d;
  void* workarea_h;
  
  std::vector<void*> mem_d;
  std::vector<void*> mem_h;
  
  bool allocated_d, allocated_h;
  bool cuda_aware;
  bool initialized;
  bool half_batch;
  bool fft3d;
  
  enum commMode_e {Peer, All2All} comm_mode;
};
