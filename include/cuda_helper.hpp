#pragma once

#ifndef __CUDA_HELPER_HPP__
#define __CUDA_HELPER_HPP__

#include <cuda_runtime.h>

#define cudaCheckError(ans) cudaAssert((ans), __FILE__, __LINE__,false)
#define cudaCheckKernelError() cudaCheckError(cudaPeekAtLastError())
#define cudaCheckLastError() cudaCheckError(cudaGetLastError())

inline int cudaAssert(cudaError_t code, const char *file, int line, bool abort=true) {
  if (code != cudaSuccess) {
    fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
    return code;
  }
  return 0;
}

inline void cudaPrintDeviceMemory(int dev=0) {
  size_t free_mem;
  size_t total_mem;

  cudaSetDevice(dev);
  cudaMemGetInfo(&free_mem, &total_mem);

  printf("GPU %i memory usage: used = %lf MiB, free = %lf MiB, total = %lf MiB\n",
    dev,
    static_cast<double>(total_mem - free_mem)/1048576.0,
    static_cast<double>(free_mem)/1048576.0,
    static_cast<double>(total_mem)/1048576.0);
}

#endif
