/*************************************************************************
 *  Copyright (c) 2017.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CLAIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _TESTINTERPOLATION_CPP_
#define _TESTINTERPOLATION_CPP_

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "UnitTestOpt.hpp"
#include "interp3.hpp"
#ifdef REG_HAS_CUDA
#include "interp3_gpu_new.hpp"
#endif

void Reference(ScalarType &val, const ScalarType x1, const ScalarType x2, const ScalarType x3) {
  val = ( PetscSinReal(x1)*PetscSinReal(x1)
        + PetscSinReal(x2)*PetscSinReal(x2)
        + PetscSinReal(x3)*PetscSinReal(x3) )/3.0;
}


namespace reg {
  
PetscErrorCode UnitTestOpt::TestInterpolation() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << "starting interpolation unit test" << std::endl;
  
  srand(time(0));
  
  int nx[3], istart[3], isize[3], nl;
  nl = this->m_Domain.nl;
  ScalarType hx[3], scale;
  const int neval = 1000;
  
  for (int i = 0; i < 3; ++i) {
    nx[i]     = this->m_Domain.nx[i];
    isize[i]  = this->m_Domain.isize[i];
    istart[i] = this->m_Domain.istart[i];
    hx[i]     = this->m_Domain.hx[i];
  }
    
  ScalarType *grid = new ScalarType[isize[0]*isize[1]*isize[2]];
  ScalarType *q    = new ScalarType[neval*3];
  ScalarType *q1   = &q[0];
  ScalarType *q2   = &q[neval];
  ScalarType *q3   = &q[neval*2];
  ScalarType *ref  = new ScalarType[neval];
  ScalarType *eval = new ScalarType[neval];
  
  for (int i1 = 0; i1 < isize[0]; ++i1) {  // x1
    for (int i2 = 0; i2 < isize[1]; ++i2) {  // x2
      for (int i3 = 0; i3 < isize[2]; ++i3) {  // x3
        // compute coordinates (nodal grid)
        ScalarType x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
        ScalarType x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
        ScalarType x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

        // compute linear / flat index
        int i = reg::GetLinearIndex(i1, i2, i3, this->m_Domain.isize);
        
        Reference(grid[i], x1, x2, x3);
      }  // i1
    }  // i2
  }  // i3
  for (int i = 0; i < neval; ++i) {
    ScalarType x1 = static_cast<ScalarType>(rand())*2.*M_PI/static_cast<ScalarType>(RAND_MAX);
    ScalarType x2 = static_cast<ScalarType>(rand())*2.*M_PI/static_cast<ScalarType>(RAND_MAX);
    ScalarType x3 = static_cast<ScalarType>(rand())*2.*M_PI/static_cast<ScalarType>(RAND_MAX);
    
    Reference(ref[i], x1, x2, x3);
    eval[i] = 0.;
    
#ifdef REG_HAS_CUDA
    q1[i] = x1/hx[i];
    q2[i] = x2/hx[i];
    q3[i] = x3/hx[i];
#else
    q[3*i + 0] = x1;
    q[3*i + 1] = x2;
    q[3*i + 2] = x3;
#endif
  }
  
#ifdef REG_HAS_CUDA
  ScalarType *pg, *pq1, *pq2, *pq3, *pe;
  
  cudaMalloc((void**)(&pg), sizeof(ScalarType)*isize[0]*isize[1]*isize[2]);
  cudaMalloc((void**)(&pq1), sizeof(ScalarType)*neval);
  cudaMalloc((void**)(&pq2), sizeof(ScalarType)*neval);
  cudaMalloc((void**)(&pq3), sizeof(ScalarType)*neval);
  cudaMalloc((void**)(&pe), sizeof(ScalarType)*neval);
  
  cudaMemcpy(pg, grid, sizeof(ScalarType)*isize[0]*isize[1]*isize[2], cudaMemcpyHostToDevice);
  cudaMemcpy(pq1, q1, sizeof(ScalarType)*neval, cudaMemcpyHostToDevice);
  cudaMemcpy(pq2, q2, sizeof(ScalarType)*neval, cudaMemcpyHostToDevice);
  cudaMemcpy(pq3, q3, sizeof(ScalarType)*neval, cudaMemcpyHostToDevice);

  cudaTextureObject_t tex = gpuInitEmptyTexture(nx);
  float timer = 0;
  
  gpuInterp3D(pg, pq1, pq2, pq3, pe, nx, tex, &timer);
  
  cudaMemcpy(eval, pe, sizeof(ScalarType)*neval, cudaMemcpyDeviceToHost);
#else
  std::cout << "unit test not implemented" << std::endl;
#endif
  
  double error = 0;
  for (int i = 0; i < neval; ++i) {
    double local = ref[i] - eval[i];
    error += local*local;
  }
  std::cout << "the interpolation has a RMS of " << sqrt(error)/static_cast<double>(neval) << std::endl;
  
  delete[] grid;
  delete[] q;
  delete[] ref;
  delete[] eval;
  
  PetscFunctionReturn(ierr);
}

} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

