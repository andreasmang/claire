/*************************************************************************
 *  Copyright (c) 2016.
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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _DIFFERENTIATIONKERNEL_CPP_
#define _DIFFERENTIATIONKERNEL_CPP_

#include "DifferentiationKernel.hpp"
#include "CLAIREUtils.hpp"

template<int N> inline ScalarType pow(ScalarType x) {
  return x*pow<N-1>(x);
}
template<> inline ScalarType pow<0>(ScalarType x) {
  return 1;
}

/********************************************************************
 * @brief computes linear array index on GPU
 *******************************************************************/
template<int N> 
inline ScalarType ComputeNLaplaceNumber(IntType w[3]) {
  ScalarType norm = static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
  return pow<N>(norm);
}

template<int N> struct NLaplacianKernelGPU {
  static inline void call(int i, IntType w[3], 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, ScalarType b0) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w);
    ScalarType regop = b0*lapik;

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NLaplacianRegularizationKernelGPU {
  static inline void call(int i, IntType w[3], 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, ScalarType b0, ScalarType b1) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w) + b1;
    ScalarType regop = b0*lapik;

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NInvLaplacianKernelGPU {
  static inline void call (int i, IntType w[3], 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, ScalarType scale, ScalarType b0) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w);
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NInvLaplacianSqrtKernelGPU {
  static inline void call(int i, IntType w[3], 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, ScalarType scale, ScalarType b0) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w);
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/sqrt(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NInvLaplacianRegularizationKernelGPU {
  static inline void call(int i, IntType w[3], 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0, ScalarType b1) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w) + b1;
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> 
struct NInvLaplacianRegularizationSqrtKernelGPU {
  static inline void call(int i, IntType w[3], 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0, ScalarType b1) {
    ScalarType lapik = ComputeNLaplaceNumber<N>(w) + b1;
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/sqrt(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};

template<typename KernelFn, typename ... Args>
PetscErrorCode SpectralKernelCall(IntType nstart[3], IntType nx[3], IntType nl[3], Args ... args) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
    for (IntType i1 = 0; i1 < nl[0]; ++i1) {
        IntType w[3];
        w[0] = i1 + nstart[0];
        reg::ComputeWaveNumber(w[0], nx[0]);
        for (IntType i2 = 0; i2 < nl[1]; ++i2) {
            w[1] = i2 + nstart[1];
            reg::ComputeWaveNumber(w[1], nx[1]);
            for (IntType i3 = 0; i3 < nl[2]; ++i3) {
                w[2] = i3 + nstart[2];
                reg::ComputeWaveNumber(w[2], nx[2]);
                IntType i = reg::GetLinearIndex(i1, i2, i3, nl);
                
                KernelFn::call(i, w, args...);
            }
        }
    }
}  // pragma omp parallel
  
  PetscFunctionReturn(ierr);
}

namespace reg {
namespace DifferentiationKernel {
  
PetscErrorCode VectorField::Laplacian(ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCall<NLaplacianKernelGPU<1> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], b0*scale); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Laplacian(ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCall<NLaplacianRegularizationKernelGPU<1> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], b0*scale, b1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Bilaplacian(ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCall<NLaplacianKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], b0*scale); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::Bilaplacian(ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCall<NLaplacianRegularizationKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], b0*scale, b1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseBilaplacian(ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCall<NInvLaplacianKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], scale, b0); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseBilaplacianSqrt(ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCall<NInvLaplacianSqrtKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], scale, b0); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseBilaplacian(ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCall<NInvLaplacianRegularizationKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], scale, b0, b1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VectorField::InverseBilaplacianSqrt(ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = SpectralKernelCall<NInvLaplacianRegularizationSqrtKernelGPU<2> >(nstart, nx, nl, pXHat[0], pXHat[1], pXHat[2], scale, b0, b1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

} // namespace DifferentiationKernel
} // namespace reg

#endif
