#include <cuda_runtime_api.h>
#include "mpicufft.hpp"
#include <algorithm>

#if (cudaError == 0) && (cufftError == 0)
#include <stdio.h>
#include <stdlib.h>
#define cudaCheck(e) {                                           \
  int err = static_cast<int>(e);                                 \
  if(err) {                                                      \
    printf("CUDA error code %s:%d: %i\n",__FILE__,__LINE__,err); \
    exit(EXIT_FAILURE);                                          \
  }                                                              \
}
#else
#define cudaCheck(e) {e}
#endif

inline int getStartIdx(std::vector<size_t> &idx, size_t begin, size_t size,
                       const std::vector<size_t> &starts, const std::vector<size_t> &sizes, size_t offset,
                       size_t lower, size_t upper, int pcnt) {
  idx.clear();
  size_t pos = std::max(begin, lower);
  int pidx = 0;
  for (int p=0; p<pcnt; ++p) {
    if (pos >= upper) break;
    if (starts[p] + sizes[p] + offset <= std::max(begin, lower)) pidx = p + 1;
    else {
      if (pos >= begin + size) break;
      idx.push_back(pos - begin);
      pos = starts[p] + sizes[p] + offset;
    }
  }
  idx.push_back(std::min(upper - begin, size));
  return pidx;
}

template<typename T> struct cuFFT {
  using C_t = cufftComplex;
  using R_t = cufftReal;
  
  static const cufftType R2Ctype = CUFFT_R2C;
  static const cufftType C2Rtype = CUFFT_C2R;
  static const cufftType C2Ctype = CUFFT_C2C;
  
  static decltype(cufftExecR2C)* execR2C;
  static decltype(cufftExecC2R)* execC2R;
  static decltype(cufftExecC2C)* execC2C;
  
  static inline C_t* complex(void *ptr) { return static_cast<C_t*>(ptr); };
  static inline C_t* complex(const void *ptr) { return static_cast<C_t*>(const_cast<void*>(ptr)); };
  static inline R_t* real(void *ptr) { return static_cast<R_t*>(ptr); };
  static inline R_t* real(const void *ptr) { return static_cast<R_t*>(const_cast<void*>(ptr)); };
};
template<typename T> decltype(cufftExecR2C)* cuFFT<T>::execR2C = cufftExecR2C;
template<typename T> decltype(cufftExecC2R)* cuFFT<T>::execC2R = cufftExecC2R;
template<typename T> decltype(cufftExecC2C)* cuFFT<T>::execC2C = cufftExecC2C;
template<> struct cuFFT<double> {
  using C_t = cufftDoubleComplex;
  using R_t = cufftDoubleReal;
  
  static const cufftType R2Ctype = CUFFT_D2Z;
  static const cufftType C2Rtype = CUFFT_Z2D;
  static const cufftType C2Ctype = CUFFT_Z2Z;
  
  static decltype(cufftExecD2Z)* execR2C;
  static decltype(cufftExecZ2D)* execC2R;
  static decltype(cufftExecZ2Z)* execC2C;
  
  static inline C_t* complex(void *ptr) { return static_cast<C_t*>(ptr); };
  static inline C_t* complex(const void *ptr) { return static_cast<C_t*>(const_cast<void*>(ptr)); };
  static inline R_t* real(void *ptr) { return static_cast<R_t*>(ptr); };
  static inline R_t* real(const void *ptr) { return static_cast<R_t*>(const_cast<void*>(ptr)); };
};
decltype(cufftExecD2Z)* cuFFT<double>::execR2C = cufftExecD2Z;
decltype(cufftExecZ2D)* cuFFT<double>::execC2R = cufftExecZ2D;
decltype(cufftExecZ2Z)* cuFFT<double>::execC2C = cufftExecZ2Z;

template<typename T> MPIcuFFT<T>::MPIcuFFT(MPI_Comm comm, bool cuda_aware)
    : comm(comm), cuda_aware(cuda_aware) {
  domainsize = 0;
  worksize_d = 0;
  worksize_h = 0;
  allocated_d = false;
  allocated_h = false;
  workarea_d = nullptr;
  workarea_h = nullptr;
  
  initialized = false;
  half_batch = false;
  
  planR2C = 0;
  planC2R = 0;
  planC2C = 0;
  
  MPI_Comm_size(comm, &pcnt);
  MPI_Comm_rank(comm, &pidx);
  
  fft3d = (pcnt == 1);
  
  isizex.resize(pcnt, 0);
  istartx.resize(pcnt, 0);
  osizey.resize(pcnt, 0);
  ostarty.resize(pcnt, 0);
  isizez = 0;
  osizez = 0;
  
  send_req.resize(pcnt, MPI_REQUEST_NULL);
  recv_req.resize(pcnt, MPI_REQUEST_NULL);
  
  if (pcnt%2 == 1) {
    for (int i=0; i<pcnt; ++i)
      if ((pcnt + i - pidx)%pcnt != pidx)
        comm_order.push_back((pcnt + i - pidx)%pcnt);
  } else if (((pcnt-1)&pcnt) == 0) {
    for (int i=1; i<pcnt; ++i)
      comm_order.push_back(pidx^i);
  } else {
    for (int i=0; i<pcnt-1;++i) {
      int idle = (pcnt*i/2)%(pcnt-1);
      if (pidx == pcnt - 1) comm_order.push_back(idle);
      else if (pidx == idle) comm_order.push_back(pcnt - 1);
      else comm_order.push_back((pcnt + i - pidx - 1)%(pcnt-1));
    }
  }
}

template<typename T> MPIcuFFT<T>::~MPIcuFFT() {
  if (allocated_d && workarea_d) cudaFree(workarea_d);
  if (allocated_h && workarea_h) cudaCheck(cudaFreeHost(workarea_h));
  if (planR2C) cudaCheck(cufftDestroy(planR2C));
  if (planC2R) cudaCheck(cufftDestroy(planC2R));
  if (planC2C) cudaCheck(cufftDestroy(planC2C));
}

template<typename T> void MPIcuFFT<T>::initFFT(size_t nx, size_t ny, size_t nz, bool allocate) {
  size_t N1    = nx/pcnt;
  size_t N1mod = nx%pcnt;
  for (int p=0; p<pcnt; ++p) {
    isizex[p]  = N1 + (static_cast<size_t>(p) < N1mod?1:0);
    istartx[p] = (p==0?0:istartx[p-1]+isizex[p-1]);
  }
  isizey = ny; isizez = nz;
  half_batch = (isizex[pidx]%2 == 0);
  
  size_t N2    = ny/pcnt;
  size_t N2mod = ny%pcnt;
  for (int p=0; p<pcnt; ++p) {
    osizey[p]  = N2 + (static_cast<size_t>(p) < N2mod?1:0);
    ostarty[p] = (p==0?0:ostarty[p-1]+osizey[p-1]);
  }
  osizex = nx; osizez = (nz / 2) + 1;
  
  domainsize = 2*sizeof(T)*std::max((isizex[pidx]*isizey*isizez/2) + 1, osizex*osizey[pidx]*osizez);
  
  size_t ws_r2c, ws_c2r, ws_c2c;

  cudaCheck(cufftCreate(&planC2R));
  cudaCheck(cufftCreate(&planR2C));
  
  cudaCheck(cufftSetAutoAllocation(planR2C, 0));
  cudaCheck(cufftSetAutoAllocation(planC2R, 0));
  
  if (fft3d) { // combined 3d fft
    cudaCheck(cufftMakePlan3d(planR2C, isizex[pidx], isizey, isizez, cuFFT<T>::R2Ctype, &ws_r2c));
    cudaCheck(cufftMakePlan3d(planC2R, isizex[pidx], isizey, isizez, cuFFT<T>::C2Rtype, &ws_c2r));

    fft_worksize = std::max(ws_r2c, ws_c2r);
  } else { // 2d slab decomposition fft
    size_t batch = (half_batch?isizex[pidx]/2:isizex[pidx]);
    
    cudaCheck(cufftCreate(&planC2C));
    cudaCheck(cufftSetAutoAllocation(planC2C, 0));
    
    long long n[3] = {static_cast<long long>(osizex), static_cast<long long>(isizey), static_cast<long long>(isizez)};
    long long nembed[1] = {1};
    
    cudaCheck(cufftMakePlanMany64(planR2C, 2, &n[1], 0, 0, 0, 0, 0, 0, cuFFT<T>::R2Ctype, batch, &ws_r2c));
    cudaCheck(cufftMakePlanMany64(planC2R, 2, &n[1], 0, 0, 0, 0, 0, 0, cuFFT<T>::C2Rtype, batch, &ws_c2r));
    cudaCheck(cufftMakePlanMany64(planC2C, 1, n, nembed, osizey[pidx]*osizez, 1, nembed, osizey[pidx]*osizez, 1, cuFFT<T>::C2Ctype, osizey[pidx]*osizez, &ws_c2c));
    
    fft_worksize = std::max(ws_r2c, ws_c2r);
    fft_worksize = std::max(fft_worksize, ws_c2c);
  }
  
  if (fft_worksize < domainsize) fft_worksize = domainsize;
  worksize_d = fft_worksize + (fft3d ? 0 : 2*domainsize);
  //worksize_h = (cuda_aware || fft3d ? 0 : 2*domainsize);
  worksize_h = (cuda_aware || fft3d ? 0 : 2*domainsize);
  if (allocate) this->setWorkArea();
  
  cudaCheck(cudaDeviceSynchronize());
}

template<typename T> void MPIcuFFT<T>::setWorkArea(void *device, void *host) {
  if (!domainsize) return;
  if (device && allocated_d) {
    cudaCheck(cudaFree(workarea_d));
    allocated_d = false;
    workarea_d = device;
  } else if (!allocated_d && device) {
    workarea_d = device;
  } else if (!allocated_d && !device) {
    cudaCheck(cudaMalloc(&workarea_d, worksize_d));
    allocated_d = true;
  }
  mem_d.clear();
  for (size_t i=0; i<(worksize_d/domainsize); ++i) mem_d.push_back(&static_cast<char*>(workarea_d)[i*domainsize]);
  
  if (fft3d) {
    cudaCheck(cufftSetWorkArea(planR2C, mem_d[0]));
    cudaCheck(cufftSetWorkArea(planC2R, mem_d[0]));
  } else {
    cudaCheck(cufftSetWorkArea(planR2C, mem_d[2]));
    cudaCheck(cufftSetWorkArea(planC2R, mem_d[2]));
    cudaCheck(cufftSetWorkArea(planC2C, mem_d[2]));
  }
    
  if (host && allocated_h) {
    cudaCheck(cudaFreeHost(workarea_h));
    allocated_h = false;
    workarea_h = host;
  } else if (!allocated_h && host) {
    workarea_h = host;
  } else if (!allocated_h && !host && worksize_h) {
    cudaCheck(cudaMallocHost(&workarea_h, worksize_h));
    allocated_h = true;
  }
  mem_h.clear();
  for (size_t i=0; i<(worksize_h/domainsize); ++i) mem_h.push_back(&static_cast<char*>(workarea_h)[i*domainsize]);
  initialized = true;
}

template<typename T> void MPIcuFFT<T>::execR2C(void *out, const void *in) {
  if (!initialized) return;
  using R_t = typename cuFFT<T>::R_t;
  using C_t = typename cuFFT<T>::C_t;
  R_t *real    = cuFFT<T>::real(in);
  C_t *complex = cuFFT<T>::complex(out);
  if (fft3d) {
    cuFFT<T>::execR2C(planR2C, real, complex);
    cudaCheck(cudaDeviceSynchronize());
  } else {
    C_t *recv_ptr, *send_ptr, *temp_ptr;
    temp_ptr = cuFFT<T>::complex(mem_d[0]);
    if (cuda_aware) {
      recv_ptr = cuFFT<T>::complex(mem_d[0]); // = temp_ptr!
      send_ptr = cuFFT<T>::complex(mem_d[1]);
    } else {
      recv_ptr = cuFFT<T>::complex(mem_h[0]);
      send_ptr = cuFFT<T>::complex(mem_h[1]);
    }
    recv_req[pidx] = MPI_REQUEST_NULL;
    send_req[pidx] = MPI_REQUEST_NULL;
    // compute 2d FFT 
    cudaCheck(cuFFT<T>::execR2C(planR2C, real, complex));
    cudaCheck(cudaDeviceSynchronize());
    size_t batch = 0;
    if (half_batch) { // if batch is halfed, compute second half
      batch = isizex[pidx]/2;
      for (int p=0; p<pcnt; ++p) {
        if (p == pidx) continue;
        size_t oslice = isizex[pidx]*osizez*ostarty[p];
        cudaCheck(cudaMemcpy2DAsync(&send_ptr[oslice], sizeof(C_t)*osizey[p]*osizez,
                                    &complex[ostarty[p]*osizez], sizeof(C_t)*isizey*osizez,
                                    sizeof(C_t)*osizey[p]*osizez, batch,
                                    cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyDeviceToHost));
      }
      // compute missing half of 2d FFT
      cudaCheck(cuFFT<T>::execR2C(planR2C, &real[isizez*isizey*batch], &complex[osizez*isizey*batch]));
      cudaCheck(cudaDeviceSynchronize());
    }
    //MPI_Waitall(pcnt, send_req.data(), MPI_STATUSES_IGNORE);
    for (auto p : comm_order) { // transpose each (missing) block and send it to respective process
    //for (int i=1; i<pcnt; ++i) { // transpose each (missing) block and send it to respective process
      //if (p == pidx) continue;
      //int p = (pcnt + pidx - i)%pcnt;
      MPI_Irecv((&recv_ptr[istartx[p]*osizez*osizey[pidx]]),
        sizeof(C_t)*isizex[p]*osizez*osizey[pidx], MPI_BYTE,
        p, p, comm, &recv_req[p]);
      //p = (pidx + i)%pcnt;
      size_t oslice = isizex[pidx]*osizez*ostarty[p];
      cudaCheck(cudaMemcpy2DAsync(&send_ptr[oslice + batch*osizez*osizey[p]], sizeof(C_t)*osizey[p]*osizez,
                                  &complex[batch*osizez*isizey + ostarty[p]*osizez], sizeof(C_t)*isizey*osizez,
                                  sizeof(C_t)*osizey[p]*osizez, isizex[pidx]-batch,
                                  cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyDeviceToHost));
      cudaCheck(cudaDeviceSynchronize());
      MPI_Isend(&send_ptr[oslice], 
                sizeof(C_t)*isizex[pidx]*osizez*osizey[p], MPI_BYTE, 
                p, pidx, comm, &send_req[p]);
    }
    { // transpose local block
      size_t oslice = osizez*osizey[pidx]*istartx[pidx];
      cudaCheck(cudaMemcpy2DAsync(&temp_ptr[oslice], sizeof(C_t)*osizey[pidx]*osizez,
                                  &complex[ostarty[pidx]*osizez], sizeof(C_t)*isizey*osizez,
                                  sizeof(C_t)*osizey[pidx]*osizez, isizex[pidx],
                                  cudaMemcpyDeviceToDevice));
    }
    if (!cuda_aware) { // copy received blocks to device
      int p;
      do {
        MPI_Waitany(pcnt, recv_req.data(), &p, MPI_STATUSES_IGNORE);
        if (p == MPI_UNDEFINED) break;
        cudaCheck(cudaMemcpyAsync(&temp_ptr[istartx[p]*osizez*osizey[pidx]],
                                  &recv_ptr[istartx[p]*osizez*osizey[pidx]],
                                  isizex[p]*osizez*osizey[pidx]*sizeof(C_t), cudaMemcpyHostToDevice));
      } while(p != MPI_UNDEFINED);
    } else { // just wait for all receives
      MPI_Waitall(pcnt, recv_req.data(), MPI_STATUSES_IGNORE);
    }
    cudaCheck(cudaDeviceSynchronize());
    // compute remaining 1d FFT, for cuda-aware recv and temp buffer are identical
    cudaCheck(cuFFT<T>::execC2C(planC2C, temp_ptr, complex, CUFFT_FORWARD));
    cudaCheck(cudaDeviceSynchronize());
    MPI_Waitall(pcnt, send_req.data(), MPI_STATUSES_IGNORE);
  }
}

template<typename T> void MPIcuFFT<T>::execC2R(void *out, const void *in) {
  if (!initialized) return;
  using R_t = typename cuFFT<T>::R_t;
  using C_t = typename cuFFT<T>::C_t;
  R_t *real    = cuFFT<T>::real(out);
  C_t *complex = cuFFT<T>::complex(in);
  if (fft3d) {
    cuFFT<T>::execC2R(planC2R, complex, real);
    cudaCheck(cudaDeviceSynchronize());
  } else {
    C_t *recv_ptr, *send_ptr, *temp_ptr, *copy_ptr;
    temp_ptr = cuFFT<T>::complex(mem_d[0]);
    copy_ptr = cuFFT<T>::complex(mem_d[1]);
    if (cuda_aware) {
      recv_ptr = cuFFT<T>::complex(mem_d[2]);
      send_ptr = cuFFT<T>::complex(mem_d[0]); // == temp_ptr!
    } else {
      recv_ptr = cuFFT<T>::complex(mem_h[0]);
      send_ptr = cuFFT<T>::complex(mem_h[1]);
    }
    recv_req[pidx] = MPI_REQUEST_NULL;
    send_req[pidx] = MPI_REQUEST_NULL;
    // compute 1d complex to complex FFT in x direction
    cudaCheck(cuFFT<T>::execC2C(planC2C, complex, temp_ptr, CUFFT_INVERSE));
    cudaCheck(cudaDeviceSynchronize());
    //MPI_Waitall(pcnt, send_req.data(), MPI_STATUSES_IGNORE);
    for (auto p : comm_order) { // copy blocks and send to respective process
    //for (int i=1; i<pcnt; ++i) { // copy blocks and send to respective process
      //if (p == pidx) continue;
      //int p = (pcnt + pidx - i)%pcnt;
      MPI_Irecv(&recv_ptr[isizex[pidx]*osizez*ostarty[p]],
                sizeof(C_t)*isizex[pidx]*osizez*osizey[p], MPI_BYTE,
                p, p, comm, &recv_req[p]);
      //p = (pidx + i)%pcnt;
      if (!cuda_aware) {
        cudaCheck(cudaMemcpy(&send_ptr[istartx[p]*osizez*osizey[pidx]],
                             &temp_ptr[istartx[p]*osizez*osizey[pidx]],
                             isizex[p]*osizez*osizey[pidx]*sizeof(C_t), cudaMemcpyDeviceToHost));
        cudaCheck(cudaDeviceSynchronize());
      }
      MPI_Isend(&send_ptr[istartx[p]*osizez*osizey[pidx]], 
                sizeof(C_t)*isizex[p]*osizez*osizey[pidx], MPI_BYTE, 
                p, pidx, comm, &send_req[p]);
    }
    { // transpose local block
      size_t islice = istartx[pidx]*osizez*osizey[pidx];
      cudaCheck(cudaMemcpy2DAsync(&copy_ptr[ostarty[pidx]*osizez], sizeof(C_t)*osizez*isizey,
                                  &temp_ptr[islice], sizeof(C_t)*osizey[pidx]*osizez,
                                  sizeof(C_t)*osizey[pidx]*osizez, isizex[pidx],
                                  cudaMemcpyDeviceToDevice));
    }
    { // transpose received data blocks
      int p;
      do {
        MPI_Waitany(pcnt, recv_req.data(), &p, MPI_STATUSES_IGNORE);
        if (p == MPI_UNDEFINED) break;
        size_t islice = isizex[pidx]*osizez*ostarty[p];
        cudaCheck(cudaMemcpy2DAsync(&copy_ptr[ostarty[p]*osizez], sizeof(C_t)*osizez*isizey,
                                    &recv_ptr[islice], sizeof(C_t)*osizey[p]*osizez,
                                    sizeof(C_t)*osizey[p]*osizez, isizex[pidx],
                                    cudaMemcpyDeviceToDevice));
      } while(p != MPI_UNDEFINED);
    }
    cudaCheck(cudaDeviceSynchronize());
    cudaCheck(cuFFT<T>::execC2R(planC2R, copy_ptr, real));
    cudaCheck(cudaDeviceSynchronize());
    if (half_batch) {
      size_t batch = isizex[pidx]/2;
      cudaCheck(cuFFT<T>::execC2R(planC2R, &copy_ptr[osizez*isizey*batch], &real[isizez*isizey*batch]));
      cudaCheck(cudaDeviceSynchronize());
    }
    MPI_Waitall(pcnt, send_req.data(), MPI_STATUSES_IGNORE);
  }
}

template<typename T> void MPIcuFFT<T>::changeSize(void *out, const void *in, MPIcuFFT<T>* fft_coarse, changeType_e dir) {
  if (!initialized || !fft_coarse->initialized) return;
  MPIcuFFT<T> *fft_o = nullptr, *fft_i = nullptr;
  using C_t = typename cuFFT<T>::C_t;
  C_t *data_o = cuFFT<T>::complex(out);
  C_t *data_i = cuFFT<T>::complex(in);
  switch(dir) {
    case Prolong:
      fft_i = fft_coarse;
      fft_o = this;
      break;
    case Restrict:
      fft_o = fft_coarse;
      fft_i = this;
      break;
  };
  
  cudaCheck(cudaMemset(out, 0, fft_o->domainsize));
  if (pcnt <= 2) { // coarse and fine data is aligned -> only local copies
    size_t pitch_o = fft_o->osizez*sizeof(C_t);
    size_t pitch_i = fft_i->osizez*sizeof(C_t);
    size_t width = (fft_coarse->osizez-1)*sizeof(C_t);
    size_t height_l = fft_coarse->isizey/2;
    size_t height_h = fft_coarse->isizey/2;
    size_t half_x = fft_coarse->osizex/2;
    
    if (fft_coarse->isizey%2 == 0) height_h -= 1; // skip Nyquist
      
    for (size_t x=0; x<half_x; ++x) {
      size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*x;
      size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*x;
      if (pidx == 0) {
        cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                    &data_i[offset_i], pitch_i,
                                    width, height_l,
                                    cudaMemcpyDeviceToDevice));
      }
      if (pidx == 1 || pcnt == 1) {
        offset_o += fft_o->osizez*(fft_o->osizey[pidx] - height_h);
        offset_i += fft_i->osizez*(fft_i->osizey[pidx] - height_h);
        cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                    &data_i[offset_i], pitch_i,
                                    width, height_h,
                                    cudaMemcpyDeviceToDevice));
      }
    }
    if (fft_coarse->osizex%2 == 0) half_x -= 1; // skip Nyquist
    for (size_t x=0; x<half_x; ++x) {
      size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*(fft_o->osizex - half_x + x);
      size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*(fft_i->osizex - half_x + x);
      if (pidx == 0) {
        cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                    &data_i[offset_i], pitch_i,
                                    width, height_l,
                                    cudaMemcpyDeviceToDevice));
      }
      if (pidx == 1 || pcnt == 1) {
        offset_o += fft_o->osizez*(fft_o->osizey[pidx] - height_h);
        offset_i += fft_i->osizez*(fft_i->osizey[pidx] - height_h);
        cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                    &data_i[offset_i], pitch_i,
                                    width, height_h,
                                    cudaMemcpyDeviceToDevice));
      }
    }
  } else {
    C_t *recv_ptr, *send_ptr;
    if (cuda_aware) {
      recv_ptr = cuFFT<T>::complex(mem_d[0]);
      send_ptr = cuFFT<T>::complex(mem_d[1]);
    } else {
      recv_ptr = cuFFT<T>::complex(mem_h[0]);
      send_ptr = cuFFT<T>::complex(mem_h[1]);
    }
    
    size_t pitch_o = fft_o->osizez*sizeof(C_t);
    size_t pitch_i = fft_i->osizez*sizeof(C_t);
    size_t width = (fft_coarse->osizez-1)*sizeof(C_t);
    size_t height_l = fft_coarse->isizey/2;
    size_t height_h = fft_coarse->isizey/2;
    size_t half_x = fft_coarse->osizex/2;
    
    size_t pitch = (fft_coarse->osizex-1)*(fft_coarse->osizez-1);
    
    if (fft_coarse->isizey%2 == 0) height_h -= 1; // skip Nyquist
    
    int psend;
    int precv;
    std::vector<size_t> isend;
    std::vector<size_t> irecv;
    size_t selfrecv = 0;

    psend = getStartIdx(isend, fft_i->ostarty[pidx], fft_i->osizey[pidx],
                        fft_o->ostarty, fft_o->osizey, 0, 0, height_l, pcnt);
    precv = getStartIdx(irecv, fft_o->ostarty[pidx], fft_o->osizey[pidx],
                        fft_i->ostarty, fft_i->osizey, 0, 0, height_l, pcnt);
                        
    for (size_t i=0; i<irecv.size()-1; ++i) {
      int p = precv + i;
      if (p == pidx) {
        selfrecv = irecv[i];
        continue;
      }
      size_t sizey = irecv[i+1] - irecv[i];
      MPI_Irecv(&recv_ptr[irecv[i]*pitch],
                sizeof(C_t)*pitch*sizey, MPI_BYTE,
                p, p, comm, &recv_req[i]);
    }
    for (size_t i=0; i<isend.size()-1; ++i) {
      int p = psend + i;
      size_t sizey = isend[i+1] - isend[i];
      if (p == pidx) {
        half_x = fft_coarse->osizex/2;
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*x;
          size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*x;
          offset_i += isend[i]*fft_i->osizez;
          offset_o += selfrecv*fft_o->osizez;
          cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                      &data_i[offset_i], pitch_i,
                                      width, sizey,
                                      cudaMemcpyDeviceToDevice));
        }
        if (fft_coarse->osizex%2 == 0) half_x -= 1; // skip Nyquist
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*(fft_o->osizex - half_x + x);
          size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*(fft_i->osizex - half_x + x);
          offset_i += isend[i]*fft_i->osizez;
          offset_o += selfrecv*fft_o->osizez;
          cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                      &data_i[offset_i], pitch_i,
                                      width, sizey,
                                      cudaMemcpyDeviceToDevice));
        }
      } else {
        half_x = fft_coarse->osizex/2;
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*x;
          offset_i += isend[i]*fft_i->osizez;
          size_t slice = (fft_coarse->osizez-1)*sizey*x;
          cudaCheck(cudaMemcpy2DAsync(&send_ptr[pitch*isend[i] + slice], width,
                                      &data_i[offset_i], pitch_i,
                                      width, sizey,
                                      cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyDeviceToHost));
        }
        size_t offset = half_x*sizey*(fft_coarse->osizez-1);
        if (fft_coarse->osizex%2 == 0) half_x -= 1; // skip Nyquist
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*(fft_i->osizex - half_x + x);
          offset_i += isend[i]*fft_i->osizez;
          size_t slice = (fft_coarse->osizez-1)*sizey*x;
          cudaCheck(cudaMemcpy2DAsync(&send_ptr[pitch*isend[i] + offset + slice], width,
                                      &data_i[offset_i], pitch_i,
                                      width, sizey,
                                      cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyDeviceToHost));
        }
        cudaDeviceSynchronize();
        MPI_Isend(&send_ptr[isend[i]*pitch],
                  sizeof(C_t)*pitch*sizey, MPI_BYTE,
                  p, pidx, comm, &send_req[i]);
      }
    }
    if (irecv.size() > 1) {
      int i;
      do {
        MPI_Waitany(irecv.size()-1, recv_req.data(), &i, MPI_STATUSES_IGNORE);
        if (i == MPI_UNDEFINED) break;
        int p = psend + i;
        size_t sizey = irecv[i+1] - irecv[i];
        half_x = fft_coarse->osizex/2;
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*x;
          offset_o += irecv[i]*fft_o->osizez;
          size_t slice = (fft_coarse->osizez-1)*sizey*x;
          cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                      &recv_ptr[pitch*irecv[i] + slice], width,
                                      width, sizey,
                                      cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyHostToDevice));
        }
        size_t offset = half_x*sizey*(fft_coarse->osizez-1);
        if (fft_coarse->osizex%2 == 0) half_x -= 1; // skip Nyquist
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*(fft_o->osizex - half_x + x);
          offset_o += irecv[i]*fft_o->osizez;
          size_t slice = (fft_coarse->osizez-1)*sizey*x;
          cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                      &recv_ptr[pitch*irecv[i] + offset + slice], width,
                                      width, sizey,
                                      cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyHostToDevice));
        }
      } while(i != MPI_UNDEFINED);
    }
    if (isend.size() > 1) {
      MPI_Waitall(isend.size()-1, send_req.data(), MPI_STATUSES_IGNORE);
    }
    cudaCheck(cudaDeviceSynchronize());
    
    size_t offset = isizey - fft_coarse->isizey;
    psend = getStartIdx(isend, fft_i->ostarty[pidx]+(dir==Restrict?0:offset), fft_i->osizey[pidx],
                        fft_o->ostarty, fft_o->osizey, (dir==Restrict?offset:0),
                        isizey-height_h, isizey, pcnt);
    precv = getStartIdx(irecv, fft_o->ostarty[pidx]+(dir==Restrict?offset:0), fft_o->osizey[pidx],
                        fft_i->ostarty, fft_i->osizey, (dir==Restrict?0:offset),
                        isizey-height_h, isizey, pcnt);
    
    for (size_t i=0; i<irecv.size()-1; ++i) {
      int p = precv + i;
      if (p == pidx) {
        selfrecv = irecv[i];
        continue;
      }
      size_t sizey = irecv[i+1] - irecv[i];
      MPI_Irecv(&recv_ptr[irecv[i]*pitch],
                sizeof(C_t)*pitch*sizey, MPI_BYTE,
                p, p, comm, &recv_req[i]);
    }
    for (size_t i=0; i<isend.size()-1; ++i) {
      int p = psend + i;
      size_t sizey = isend[i+1] - isend[i];
      if (p == pidx) {
        half_x = fft_coarse->osizex/2;
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*x;
          size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*x;
          offset_i += isend[i]*fft_i->osizez;
          offset_o += selfrecv*fft_o->osizez;
          cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                      &data_i[offset_i], pitch_i,
                                      width, sizey,
                                      cudaMemcpyDeviceToDevice));
        }
        if (fft_coarse->osizex%2 == 0) half_x -= 1; // skip Nyquist
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*(fft_o->osizex - half_x + x);
          size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*(fft_i->osizex - half_x + x);
          offset_i += isend[i]*fft_i->osizez;
          offset_o += selfrecv*fft_o->osizez;
          cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                      &data_i[offset_i], pitch_i,
                                      width, sizey,
                                      cudaMemcpyDeviceToDevice));
        }
      } else {
        half_x = fft_coarse->osizex/2;
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*x;
          offset_i += isend[i]*fft_i->osizez;
          size_t slice = (fft_coarse->osizez-1)*sizey*x;
          cudaCheck(cudaMemcpy2DAsync(&send_ptr[pitch*isend[i] + slice], width,
                                      &data_i[offset_i], pitch_i,
                                      width, sizey,
                                      cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyDeviceToHost));
        }
        size_t offset = half_x*sizey*(fft_coarse->osizez-1);
        if (fft_coarse->osizex%2 == 0) half_x -= 1; // skip Nyquist
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_i = fft_i->osizez*fft_i->osizey[pidx]*(fft_i->osizex - half_x + x);
          offset_i += isend[i]*fft_i->osizez;
          size_t slice = (fft_coarse->osizez-1)*sizey*x;
          cudaCheck(cudaMemcpy2DAsync(&send_ptr[pitch*isend[i] + offset + slice], width,
                                      &data_i[offset_i], pitch_i,
                                      width, sizey,
                                      cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyDeviceToHost));
        }
        cudaDeviceSynchronize();
        MPI_Isend(&send_ptr[isend[i]*pitch],
                  sizeof(C_t)*pitch*sizey, MPI_BYTE,
                  p, pidx, comm, &send_req[i]);
      }
    }
    if (irecv.size() > 1) {
      int i;
      do {
        MPI_Waitany(irecv.size()-1, recv_req.data(), &i, MPI_STATUSES_IGNORE);
        if (i == MPI_UNDEFINED) break;
        int p = psend + i;
        size_t sizey = irecv[i+1] - irecv[i];
        half_x = fft_coarse->osizex/2;
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*x;
          offset_o += irecv[i]*fft_o->osizez;
          size_t slice = (fft_coarse->osizez-1)*sizey*x;
          cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                      &recv_ptr[pitch*irecv[i] + slice], width,
                                      width, sizey,
                                      cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyHostToDevice));
        }
        size_t offset = half_x*sizey*(fft_coarse->osizez-1);
        if (fft_coarse->osizex%2 == 0) half_x -= 1; // skip Nyquist
        for (size_t x=0; x<half_x; ++x) {
          size_t offset_o = fft_o->osizez*fft_o->osizey[pidx]*(fft_o->osizex - half_x + x);
          offset_o += irecv[i]*fft_o->osizez;
          size_t slice = (fft_coarse->osizez-1)*sizey*x;
          cudaCheck(cudaMemcpy2DAsync(&data_o[offset_o], pitch_o,
                                      &recv_ptr[pitch*irecv[i] + offset + slice], width,
                                      width, sizey,
                                      cuda_aware?cudaMemcpyDeviceToDevice:cudaMemcpyHostToDevice));
        }
      } while(i != MPI_UNDEFINED);
    }
    if (isend.size() > 1) {
      MPI_Waitall(isend.size()-1, send_req.data(), MPI_STATUSES_IGNORE);
    }
    cudaCheck(cudaDeviceSynchronize());
  }
}

template class MPIcuFFT<double>;
template class MPIcuFFT<float>;