#include "adjoint_kernel.hpp"
#include "cuda_helper.hpp"

__global__ void bodyforce_kernel(PetscReal *p_lnext, PetscReal *p_l, PetscReal *p_lx, PetscReal *p_divv, PetscReal *p_divvx, PetscReal *p_vec1, PetscReal *p_vec2, PetscReal *p_vec3, PetscReal *p_b1, PetscReal *p_b2, PetscReal *p_b3, PetscReal ht, PetscReal scale, PetscInt nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    PetscReal lambda  = p_l[i];
    PetscReal lambdax = p_lx[i];

    PetscReal rhs0 = lambdax*p_divvx[i];
    PetscReal rhs1 = (lambdax + ht*rhs0)*p_divv[i];

    // compute \lambda(x,t^{j+1})
    p_lnext[i] = lambdax + 0.5*ht*(rhs0 + rhs1);

    // compute bodyforce
    p_b1[i] += scale*p_vec1[i]*lambda;
    p_b2[i] += scale*p_vec2[i]*lambda;
    p_b3[i] += scale*p_vec3[i]*lambda;
  }
}

void ComputeAdjointBodyForceGPU(PetscReal *p_lnext, PetscReal *p_l, PetscReal *p_lx, PetscReal *p_divv, PetscReal *p_divvx, PetscReal *p_vec1, PetscReal *p_vec2, PetscReal *p_vec3, PetscReal *p_b1, PetscReal *p_b2, PetscReal *p_b3, PetscReal ht, PetscReal scale, PetscInt nl) {
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  bodyforce_kernel<<<grid,block>>>(p_lnext,p_l,p_lx,p_divv,p_divvx,p_vec1,p_vec2,p_vec3,p_b1,p_b2,p_b3,ht,scale,nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
}

__global__ void bodyforce_kernel(PetscReal *p_l, PetscReal *p_vec1, PetscReal *p_vec2, PetscReal *p_vec3, PetscReal *p_b1, PetscReal *p_b2, PetscReal *p_b3, PetscReal scale, PetscInt nl) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < nl) {
    PetscReal lambda  = p_l[i];

    // compute bodyforce
    p_b1[i] += scale*p_vec1[i]*lambda;
    p_b2[i] += scale*p_vec2[i]*lambda;
    p_b3[i] += scale*p_vec3[i]*lambda;
  }
}

void ComputeAdjointBodyForceGPU(PetscReal *p_l, PetscReal *p_vec1, PetscReal *p_vec2, PetscReal *p_vec3, PetscReal *p_b1, PetscReal *p_b2, PetscReal *p_b3, PetscReal scale, PetscInt nl) {
  dim3 block = dim3(256);
  dim3 grid  = dim3((nl + 255)/256);
  bodyforce_kernel<<<grid,block>>>(p_l,p_vec1,p_vec2,p_vec3,p_b1,p_b2,p_b3,scale,nl);
  cudaDeviceSynchronize();
  cudaCheckKernelError();
}
