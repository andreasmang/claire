#ifndef _INTERP3_GPU_HPP_
#define _INTERP3_GPU_HPP_

#define COORD_DIM 3

#include "petsc.h"
#include "petsccuda.h"
#include <cuda.h>
#include <cuda_runtime.h>


void gpuInterp3D(PetscScalar* yi, const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3, PetscScalar* yo, int* nx, PetscScalar* interp_time);
//void getSemiLagrangianInitialCondition(PetscScalar* x1, PetscScalar* x2, PetscScalar* x3, int* nx, PetscScalar* compute_time);

#endif
