#ifndef _INTERP3_GPU_HPP_
#define _INTERP3_GPU_HPP_

#define COORD_DIM 3

#include "petsc.h"
#include "petsccuda.h"
#include <cuda.h>
#include <cuda_runtime_api.h>


void gpuInterp3D(PetscScalar* yi, const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3, PetscScalar* yo, int* nx, PetscScalar* interp_time);


//__global__ gpuInterp3DKernel(cudaTextureObject_t yi_tex, float* xq1, float* xq2, float* xq3, float* yo ,const float3 inv_nx);


#endif
