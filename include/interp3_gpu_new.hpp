#ifndef _INTERP3_GPU_HPP_
#define _INTERP3_GPU_HPP_

#define COORD_DIM 3

#include "petsc.h"

#include "TypeDef.hpp"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_helper.hpp>
/*#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/extrema.h>*/

//#define BLOCK_COORDINATES

using ScalarType = PetscScalar;

//void gpuInterp3D(PetscScalar* yi, 
//  const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3, 
//  PetscScalar* yo, 
//  float *tmp1, float* tmp2,
//  int* nx, long int nq, cudaTextureObject_t yi_tex, int iporder, PetscScalar* interp_time);
void gpuInterp3D(PetscScalar* yi, 
  const PetscScalar** query_in, 
  PetscScalar* yo, 
  float *tmp1, float* tmp2,
  IntType* nx, long int nq, cudaTextureObject_t yi_tex, int iporder, PetscScalar* interp_time);

//void gpuInterpVec3D(PetscScalar* yi1, PetscScalar* yi2, PetscScalar* yi3, 
//    const PetscScalar* xq1, const PetscScalar* xq2, const PetscScalar* xq3, 
//    PetscScalar* yo1, PetscScalar* yo2, PetscScalar* yo3, 
//    float *tmp1, float* tmp2,
//    int* nx, long int nq, cudaTextureObject_t yi_tex, int iporder, PetscScalar* interp_time);
void gpuInterpVec3D(PetscScalar* yi1, PetscScalar* yi2, PetscScalar* yi3, 
    const PetscScalar** query_in, 
    PetscScalar* yo1, PetscScalar* yo2, PetscScalar* yo3, 
    float *tmp1, float* tmp2,
    IntType* nx, long int nq, cudaTextureObject_t yi_tex, int iporder, PetscScalar* interp_time);

extern "C" cudaTextureObject_t gpuInitEmptyTexture(IntType* nx);

void interp0(float* m, float* q1, float *q2, float *q3, float *q, int nx[3]);

void normalizeQueryPoints(ScalarType* xq1, ScalarType* xq2, ScalarType* xq3, ScalarType* all_query_points, int nq, IntType* isize, IntType* nx, int* procid, int nghost);

void copyQueryValues(ScalarType* dst, ScalarType* src, int* index, int len);

void enforcePeriodicity(ScalarType* xq, ScalarType* yq, ScalarType* zq, ScalarType* h, int len);

void checkDomain(short* which_proc, ScalarType* xq, ScalarType* yq, ScalarType* zq, ScalarType* iX0, ScalarType* iX1, ScalarType* h, int len, int procid, int isize0, int isize1, int c_dim1);

void initializeGrid(ScalarType* xq, ScalarType* yq, ScalarType* zq, ScalarType* f, ScalarType* ref, ScalarType* h, IntType* isize, IntType* istart, IntType* nx, int caseid);

void test(ScalarType* f, int nq);

#endif
