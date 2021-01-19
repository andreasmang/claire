#include "TwoLevel.hpp"

template<typename T>
__global__ void RestrictCubicCentralGPU(T* __restrict__ dst,
                                        const T* __restrict__ src,
                                        dim3 nx) {
  T stage0[4];
  T stage1[4];
  
  dim3 gidx;
  
  gidx.x = blockIdx.z;
  gidx.y = blockIdx.y;
  gidx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (gidx.z*2 >= nx.z) return;
  
  unsigned int stride_y = nx.z;
  unsigned int stride_x = nx.z*nx.y;
  
  unsigned int idx_x[4], idx_y[4], idx_z[4];
  
  idx_x[0] = (gidx.x == 0 ? nx.x - 1 : gidx.x*2 - 1);
  idx_x[1] = gidx.x*2;
  idx_x[2] = idx_x[1] + 1;
  idx_x[3] = (idx_x[2] + 1 == nx.x ? 0 : idx_x[2] + 1);
  
  idx_y[0] = (gidx.y == 0 ? nx.y - 1 : gidx.y*2 - 1);
  idx_y[1] = gidx.y*2;
  idx_y[2] = idx_y[1] + 1;
  idx_y[3] = (idx_y[2] + 1 == nx.y ? 0 : idx_y[2] + 1);
  
  idx_z[0] = (gidx.z == 0 ? nx.z - 1 : gidx.z*2 - 1);
  idx_z[1] = gidx.z*2;
  idx_z[2] = idx_z[1] + 1;
  idx_z[3] = (idx_z[2] + 1 == nx.z ? 0 : idx_z[2] + 1);
  
  for (int x=0; x<4; ++x) {
    unsigned int base_x = idx_x[x]*stride_x;
    for (int y=0; y<4; ++y) {
      const T* ptr = &src[base_x + idx_y[y]*stride_y];
      stage0[y] = 0.0625*(- ptr[idx_z[0]] + 9.*ptr[idx_z[1]] + 9.*ptr[idx_z[2]] - ptr[idx_z[3]]);
    }
    stage1[x] = 0.0625*(- stage0[0] + 9.*stage0[1] + 9.*stage0[2] - stage0[3]);
  }
  
  stride_y /= 2;
  stride_x /= 4;
  dst[gidx.z + gidx.y*stride_y + gidx.x*stride_x] = 0.0625*(- stage1[0] + 9.*stage1[1] + 9.*stage1[2] - stage1[3]);
}

template<typename T>
__global__ void ProlongCubicCentralGPU(T* __restrict__ dst,
                                       const T* __restrict__ src,
                                       dim3 nx) {
  T val[8];
  
  dim3 idx;
  
  idx.x = blockIdx.z;
  idx.y = blockIdx.y;
  idx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx.z >= nx.z) return;
  
  unsigned int stride_x = nx.y*nx.z;
  unsigned int stride_y = nx.z;
  
  const T weight[5] = {-3./128., 29./128., 111./128., -9./128., 0.};
  val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.;
  
  unsigned int idx_z[5];
  
  idx_z[0] = (idx.z >= 2 ? idx.z - 2 : nx.z - 2 + idx.z);
  idx_z[1] = (idx.z >= 1 ? idx.z - 1 : nx.z - 1);
  idx_z[2] =  idx.z;
  idx_z[3] = (idx.z + 1 < nx.z ? idx.z + 1 : 0);
  idx_z[4] = (idx.z + 2 < nx.z ? idx.z + 2 : idx.z + 2 - nx.z);

  for (int x = 0; x < 5; ++x) {
    unsigned int idx_x = idx.x + x;
    idx_x = (idx_x >= 2 ? (idx_x - 2 >= nx.x ? idx_x - 2 - nx.x : idx_x - 2) : nx.x - 2 + idx_x);
    idx_x *= stride_x;
    for (int y = 0; y < 5; ++y) {
      unsigned int idx_y = idx.y + y;
      idx_y = (idx_y >= 2 ? (idx_y - 2 >= nx.y ? idx_y - 2 - nx.y : idx_y - 2) : nx.y - 2 + idx_y);
      const T* ptr = &src[idx_x + idx_y*stride_y];
      T cw[4];
      cw[0] = weight[y]  *weight[x];
      cw[1] = weight[4-y]*weight[x];
      cw[2] = weight[y]  *weight[4-x];
      cw[3] = weight[4-y]*weight[4-x];
      for (int z = 0; z < 5; ++z) {
        T inval = ptr[idx_z[z]];
        val[0] += inval*weight[z]  *cw[0];
        val[1] += inval*weight[4-z]*cw[0];
        val[2] += inval*weight[z]  *cw[1];
        val[3] += inval*weight[4-z]*cw[1];
        val[4] += inval*weight[z]  *cw[2];
        val[5] += inval*weight[4-z]*cw[2];
        val[6] += inval*weight[z]  *cw[3];
        val[7] += inval*weight[4-z]*cw[3];
      }
    }
  }
  
  stride_y *= 2;
  stride_x *= 4;
  
  unsigned int base = (idx.z + idx.y*stride_y + idx.x*stride_x)*2;
  
  dst[base + 0 + 0*stride_y + 0*stride_x] = val[0];
  dst[base + 1 + 0*stride_y + 0*stride_x] = val[1];
  dst[base + 0 + 1*stride_y + 0*stride_x] = val[2];
  dst[base + 1 + 1*stride_y + 0*stride_x] = val[3];
  dst[base + 0 + 0*stride_y + 1*stride_x] = val[4];
  dst[base + 1 + 0*stride_y + 1*stride_x] = val[5];
  dst[base + 0 + 1*stride_y + 1*stride_x] = val[6];
  dst[base + 1 + 1*stride_y + 1*stride_x] = val[7];
}

namespace reg {
  
namespace TwoLevelFiniteKernel {

PetscErrorCode InitConstants(IntType* isize, IntType* isize_g, ScalarType* hx, IntType* halo) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode Restrict(ScalarType* dst, const ScalarType *src, IntType *nl_f) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  IntType nl[3];
  nl[0] = nl_f[0]/2; nl[1] = nl_f[1]/2; nl[2] = nl_f[2]/2;
  
  dim3 block, grid;
  
  if (nl[2] <= 1024 && nl[2] >= 32) {
    block.x = (nl[2] + 31)/32;
    block.x *= 32;
    grid.x = 1;
  } else {
    block.x = 128; // 128 threads per block
    grid.x = (nl[2] + 127)/128;  // $\lceil nl_2 / 128 \rceil$
  }
  grid.y = nl[1];
  grid.z = nl[0];
  dim3 nl3;
  nl3.x = nl_f[0]; nl3.y = nl_f[1]; nl3.z = nl_f[2];

  RestrictCubicCentralGPU<<<grid, block>>>(dst, src, nl3);
  ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
  ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);  
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode Prolong(ScalarType* dst, const ScalarType *src, IntType *nl_f) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  IntType nl[3];
  nl[0] = nl_f[0]/2; nl[1] = nl_f[1]/2; nl[2] = nl_f[2]/2;
  
  dim3 block, grid;
  
  if (nl[2] <= 1024 && nl[2] >= 32) {
    block.x = (nl[2] + 31)/32;
    block.x *= 32;
    grid.x = 1;
  } else {
    block.x = 128; // 128 threads per block
    grid.x = (nl[2] + 127)/128;  // $\lceil nl_2 / 128 \rceil$
  }
  grid.y = nl[1];
  grid.z = nl[0];
  dim3 nl3;
  nl3.x = nl[0]; nl3.y = nl[1]; nl3.z = nl[2];

  ProlongCubicCentralGPU<<<grid, block>>>(dst, src, nl3);
  ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
  ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);  
  
  PetscFunctionReturn(ierr);
}

} // namespace TwoLevelFiniteKernel
} // namespace reg
