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

void TestFunction(ScalarType &val, const ScalarType x1, const ScalarType x2, const ScalarType x3) {
  val = ( PetscSinReal(x1)*PetscSinReal(x1)
        + PetscSinReal(x2)*PetscSinReal(x2)
        + PetscSinReal(x3)*PetscSinReal(x3) )/3.0;
}


namespace reg {
namespace UnitTest {
  
PetscErrorCode TestInterpolation(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << "starting interpolation unit test" << std::endl;
  
  srand(time(0));
  
  int nx[3], istart[3], isize[3], nl;
  ScalarType hx[3], scale;
  
  nl = 1;
  for (int i = 0; i < 3; ++i) {
    nx[i]     = m_Opt->m_Domain.nx[i];
    isize[i]  = m_Opt->m_Domain.isize[i];
    istart[i] = m_Opt->m_Domain.istart[i];
    hx[i]     = m_Opt->m_Domain.hx[i];
    nl *= isize[i];
  }
    
  ScalarType *grid = new ScalarType[nl];
  ScalarType *q    = new ScalarType[nl*3];
  ScalarType *q1   = &q[0];
  ScalarType *q2   = &q[nl];
  ScalarType *q3   = &q[nl*2];
  ScalarType *ref  = new ScalarType[nl];
  ScalarType *eval = new ScalarType[nl];
  
  for (int i1 = 0; i1 < isize[0]; ++i1) {  // x1
    for (int i2 = 0; i2 < isize[1]; ++i2) {  // x2
      for (int i3 = 0; i3 < isize[2]; ++i3) {  // x3
        // compute coordinates (nodal grid)
        double x1 = hx[0]*static_cast<double>(i1 + istart[0]);
        double x2 = hx[1]*static_cast<double>(i2 + istart[1]);
        double x3 = hx[2]*static_cast<double>(i3 + istart[2]);

        // compute linear / flat index
        int i = reg::GetLinearIndex(i1, i2, i3, m_Opt->m_Domain.isize);
        
        TestFunction(grid[i], x1, x2, x3);
    
        x1 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
        x2 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
        x3 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
        
        TestFunction(ref[i], x1, x2, x3);
        eval[i] = 0.;
    
#ifdef REG_HAS_CUDA
        q1[i] = x1*nx[0]/(2.*M_PI);
        q2[i] = x2*nx[1]/(2.*M_PI);
        q3[i] = x3*nx[2]/(2.*M_PI);
#else
        q[3*i + 0] = x1;
        q[3*i + 1] = x2;
        q[3*i + 2] = x3;
#endif
      }  // i1
    }  // i2
  }  // i3


#ifdef REG_HAS_CUDA
  ScalarType *pg, *pq1, *pq2, *pq3, *pe;
  
  cudaMalloc((void**)(&pg), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&pq1), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&pq2), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&pq3), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&pe), sizeof(ScalarType)*nl);
  
  cudaMemcpy(pg, grid, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  cudaMemcpy(pq1, q1, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  cudaMemcpy(pq2, q2, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  cudaMemcpy(pq3, q3, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);

  cudaTextureObject_t tex = gpuInitEmptyTexture(nx);
  float timer = 0;
 
  cudaDeviceSynchronize();
  gpuInterp3D(pg, pq1, pq2, pq3, pe, nx, tex, m_Opt->m_PDESolver.iporder, &timer);
  //interp0(pg,pq1,pq2,pq3,pe,nx);
  cudaDeviceSynchronize();
  
  cudaMemcpy(eval, pe, sizeof(ScalarType)*nl, cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

  cudaFree(pg);
  cudaFree(q1);
  cudaFree(q2);
  cudaFree(q3);
  cudaFree(pe);
  cudaDestroyTextureObject(tex);
#else
  std::cout << "unit test not implemented" << std::endl;
#endif

  
  double error = 0;
  double max = 0;
  for (int i = 0; i < nl; ++i) {
    double local = static_cast<double>(ref[i]) - static_cast<double>(eval[i]);
    double lmax = std::abs(eval[i]);
    if (lmax > max) max = lmax;
    error += local*local;
  }
  std::cout << "the interpolation has a RMS of " << sqrt(error/static_cast<double>(nl)) << std::endl;
  std::cout << "abs max of interpolation is " << max << std::endl;
  
  delete[] grid;
  delete[] q;
  delete[] ref;
  delete[] eval;
  
  PetscFunctionReturn(ierr);
}

}} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

