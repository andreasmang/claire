#include "Regularization.hpp"
#include "cuda_helper.hpp"

/********************************************************************
 * @brief computes wave number on GPU
 *******************************************************************/
__device__ inline void ComputeWaveNumber(dim3 &w, dim3 &n) {
  if (w.x >  n.x/2) w.x -= n.x;
  if (w.y >  n.y/2) w.y -= n.y;
  if (w.z >  n.z/2) w.z -= n.z;
}

/********************************************************************
 * @brief computes linear array index on GPU
 *******************************************************************/
__device__ inline int GetLinearIndex(int i1, int i2, int i3, dim3 size) {
  return i1*size.y*size.z + i2*size.z + i3;
}

template<int N> __device__ inline ScalarType pow(ScalarType x) {
  return x*pow<N-1>(x);
}
template<> __device__ inline ScalarType pow<0>(ScalarType x) {
  return 1;
}

/********************************************************************
 * @brief computes linear array index on GPU
 *******************************************************************/
template<int N>
__device__ inline ScalarType ComputeLaplaceNumber(dim3 w) {
  ScalarType norm = static_cast<ScalarType>(w.x*w.x + w.y*w.y + w.z*w.z);
  return pow<N>(norm);
}

/********************************************************************
 * @brief computes H[N] kernel in FFT space on GPU
 *******************************************************************/
template<int N>
__global__ void EvaluateGradientKernelGPU(ComplexType *v1, ComplexType *v2, ComplexType *v3,
    dim3 wave, dim3 nx, dim3 nl,
    ScalarType alpha, ScalarType beta) {
  int i1 = threadIdx.x + blockIdx.x*blockDim.x;
  int i2 = blockIdx.y;
  int i3 = blockIdx.z;
  
  if (i1 < nl.x) {
    wave.x += i1;
    wave.y += i2;
    wave.z += i3;

    ComputeWaveNumber(wave, nx);

    // compute bilaplacian operator and regularization operator
    ScalarType lapik = ComputeLaplaceNumber<N>(wave);
    ScalarType regop = alpha*(lapik + beta);

    // get linear index
    int i = GetLinearIndex(i1, i2, i3, nl);

    // apply to individual components
    v1[i][0] *= regop;
    v1[i][1] *= regop;
    v2[i][0] *= regop;
    v2[i][1] *= regop;
    v3[i][0] *= regop;
    v3[i][1] *= regop;
  }
}

namespace reg {

//template<int N>
//PetscErrorCode RegularizationKernel<N>::EvaluateGradient() {
template <int N>
PetscErrorCode EvaluateGradientKernel(ComplexType* v1, ComplexType* v2, ComplexType* v3, IntType nstart[3], IntType nx[3], IntType nl[3], ScalarType beta0, ScalarType beta1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  dim3 block(256,1,1);
  dim3 grid((nl[0] + 255)/256,nl[1],nl[2]);
  dim3 wave(nstart[0],nstart[1],nstart[2]);
  dim3 nx3(nx[0],nx[1],nx[2]);
  dim3 nl3(nl[0],nl[1],nl[2]);
  
  EvaluateGradientKernelGPU<N><<<grid,block>>>(v1,v2,v3,wave,nx3,nl3,beta0,beta1);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

template PetscErrorCode EvaluateGradientKernel<1>(ComplexType*, ComplexType*, ComplexType*, IntType[3], IntType[3], IntType[3], ScalarType, ScalarType);  
template PetscErrorCode EvaluateGradientKernel<2>(ComplexType*, ComplexType*, ComplexType*, IntType[3], IntType[3], IntType[3], ScalarType, ScalarType);  
template PetscErrorCode EvaluateGradientKernel<3>(ComplexType*, ComplexType*, ComplexType*, IntType[3], IntType[3], IntType[3], ScalarType, ScalarType);  

//template PetscErrorCode RegularizationKernel<1>::EvaluateGradientKernel();
//template PetscErrorCode RegularizationKernel<2>::EvaluateGradientKernel();
//template PetscErrorCode RegularizationKernel<3>::EvaluateGradientKernel();

__global__ void ScaleVectorFieldGPU(ScalarType *vR1, ScalarType *vR2, ScalarType *vR3,
    ScalarType *v1, ScalarType *v2, ScalarType *v3,
    IntType nl, ScalarType alpha) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (i < nl) {
    vR1[i] = alpha*v1[i];
    vR2[i] = alpha*v2[i];
    vR3[i] = alpha*v3[i];
  }
}

PetscErrorCode ScaleVectorField(ScalarType *vR1, ScalarType *vR2, ScalarType *vR3,
    ScalarType *v1, ScalarType *v2, ScalarType *v3,
    IntType nl, ScalarType alpha) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  dim3 block(256,1,1);
  dim3 grid((nl + 255)/256,1,1);
  
  ScaleVectorFieldGPU<<<grid,block>>>(vR1,vR2,vR3,v1,v2,v3,nl,alpha);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
  
  PetscFunctionReturn(ierr);
}

} // namespace reg
