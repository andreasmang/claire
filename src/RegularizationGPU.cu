#include "RegularizationGPU.hpp"
#include "cuda_helper.hpp"

__device__ inline void ComputeWaveNumber(dim3 &w, dim3 &n) {
  if (w.x >  n.x/2) w.x -= n.x;
  if (w.y >  n.y/2) w.y -= n.y;
  if (w.z >  n.z/2) w.z -= n.z;
}

__device__ inline int GetLinearIndex(int i1, int i2, int i3, dim3 isize) {
  return i1*isize.y*isize.z + i2*isize.z + i3;
}


__global__ void EvaluateGradientH1SN_kernel(ComplexType *v1hat, ComplexType *v2hat, ComplexType *v3hat, dim3 w, dim3 size, dim3 nx, ScalarType opfact) {
  int i1 = threadIdx.x + blockIdx.x*blockDim.x;
  int i2 = blockIdx.y;
  int i3 = blockIdx.z;
  
  if (i1 < size.x) {
    w.x += i1;
    w.y += i2;
    w.z += i3;

    ComputeWaveNumber(w, nx);

    // compute bilaplacian operator and regularization operator
    ScalarType regop = opfact*static_cast<ScalarType>(w.x*w.x + w.y*w.y + w.z*w.z);

    // get linear index
    int i = GetLinearIndex(i1, i2, i3, size);

    // apply to individual components
    v1hat[i][0] *= regop;
    v1hat[i][1] *= regop;

    v2hat[i][0] *= regop;
    v2hat[i][1] *= regop;

    v3hat[i][0] *= regop;
    v3hat[i][1] *= regop;
  }
}

void EvaluateGradientH1SN_GPU(ComplexType *v1hat, ComplexType *v2hat, ComplexType *v3hat, IntType ostart[3], IntType nl[3], IntType nx[3], ScalarType opfact) {
  dim3 block(256,1,1);
  dim3 grid((nl[0] + 255)/256,nl[1],nl[2]);
  dim3 w(ostart[0],ostart[1],ostart[2]);
  dim3 n(nx[0],nx[1],nx[2]);
  dim3 size(nl[0],nl[1],nl[2]);
  EvaluateGradientH1SN_kernel<<<grid,block>>>(v1hat,v2hat,v3hat,w,size,n,opfact);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
}



