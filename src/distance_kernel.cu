#include "distance_kernel.hpp"
#include "cuda_helper.hpp"

__global__ void sub_kernel(PetscReal *p_l, const PetscReal *p_m, const PetscReal *p_mr, int nl) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx < nl)
    p_l[idx] = p_mr[idx] - p_m[idx];
}

void DistanceMeasureSetFinalGPU(PetscReal *p_l, const PetscReal *p_m, const PetscReal *p_mr, PetscInt nl) {
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  sub_kernel<<<grid,block>>>(p_l, p_m, p_mr, nl);
  cudaCheckKernelError();
}

__global__ void submul_kernel(PetscReal *p_l, const PetscReal *p_m, const PetscReal *p_mr, const PetscReal *p_w, int nl) {
  int idx_nl = threadIdx.x + blockIdx.x*blockDim.x;
  int idx = idx_nl + blockIdx.y*nl;
  if (idx_nl < nl)
    p_l[idx] = p_w[idx_nl]*(p_mr[idx] - p_m[idx]);
}

void DistanceMeasureSetFinalMaskGPU(PetscReal *p_l, const PetscReal *p_m, const PetscReal *p_mr, const PetscReal *p_w, PetscInt nl, PetscInt nc) {
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256, nc);
  submul_kernel<<<grid,block>>>(p_l, p_m, p_mr, p_w, nl);
  cudaCheckKernelError();
}

