#ifndef _INTERP3_GPU_HPP_
#define _INTERP3_GPU_HPP_

#define COORD_DIM 3

#include "petsc.h"
#include "petsccuda.h"
#include <cuda.h>
#include <cuda_runtime.h>

void gpuInterp3D(PetscScalar* yi, 
  const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3, 
  PetscScalar* yo, 
  float *tmp1, float* tmp2,
  int* nx, cudaTextureObject_t yi_tex, int iporder, PetscScalar* interp_time);

void gpuInterpVec3D(PetscScalar* yi1, PetscScalar* yi2, PetscScalar* yi3, 
    const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3, 
    PetscScalar* yo1, PetscScalar* yo2, PetscScalar* yo3, 
    float *tmp1, float* tmp2,
    int* nx, cudaTextureObject_t yi_tex, int iporder, PetscScalar* interp_time);

extern "C" cudaTextureObject_t gpuInitEmptyTexture(int* nx);

void interp0(float* m, float* q1, float *q2, float *q3, float *q, int nx[3]);



#endif
