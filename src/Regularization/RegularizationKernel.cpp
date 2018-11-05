#include "Regularization.hpp"

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
inline ScalarType ComputeLaplaceNumber(IntType w[3]) {
  ScalarType norm = static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
  return pow<N>(norm);
}

namespace reg {

/********************************************************************
 * @brief computes H[N] kernel in FFT space
 *******************************************************************/
template<int N>
PetscErrorCode EvaluateGradientKernel(ComplexType* v1, ComplexType* v2, ComplexType* v3, IntType nstart[3], IntType nx[3], IntType nl[3], ScalarType beta0, ScalarType beta1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
    
#pragma omp parallel
{
#pragma omp for
  for (IntType i1 = 0; i1 < nl[0]; ++i1) {
    for (IntType i2 = 0; i2 < nl[1]; ++i2) {
      for (IntType i3 = 0; i3 < nl[2]; ++i3) {
        IntType w[3];
        w[0] = i1 + nstart[0];
        w[1] = i2 + nstart[1];
        w[2] = i3 + nstart[2];

        ComputeWaveNumber(w, nx);

        // compute N-laplacian operator and regularization operator
        ScalarType lapik = ComputeLaplaceNumber<N>(w);
        ScalarType regop = beta0*(lapik + beta1);

        // get linear index
        IntType i = GetLinearIndex(i1, i2, i3, nl);

        // apply to individual components
        v1[i][0] *= regop;
        v1[i][1] *= regop;
        v2[i][0] *= regop;
        v2[i][1] *= regop;
        v3[i][0] *= regop;
        v3[i][1] *= regop;
      }
    }
  }
}  // pragma omp parallel
  
  PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief computes H[N] kernel inverse in FFT space
 *******************************************************************/
/*template<int N>
PetscErrorCode RegularizationKernel<N>::ApplyInverse(bool usesqrt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
    
#pragma omp parallel
{
#pragma omp for
  for (IntType i1 = 0; i1 < nl[0]; ++i1) {
    for (IntType i2 = 0; i2 < nl[1]; ++i2) {
      for (IntType i3 = 0; i3 < nl[2]; ++i3) {
        IntType w[3];
        w[0] = i1 + nstart[0];
        w[1] = i2 + nstart[1];
        w[2] = i3 + nstart[2];

        ComputeWaveNumber(w, nx);

        // compute N-laplacian operator and regularization operator
        ScalarType lapik = ComputeLaplaceNumber<N>(w);
        ScalarType regop = (lapik + beta1==0?beta0:beta0*(lapik + beta1));
        
        if (usesqrt) regop = std::sqrt(regop);
        
        regop = scale/regop;

        // get linear index
        IntType i = GetLinearIndex(i1, i2, i3, local);

        // apply to individual components
        pXHat[0][i][0] *= regop;
        pXHat[0][i][1] *= regop;
        pXHat[1][i][0] *= regop;
        pXHat[1][i][1] *= regop;
        pXHat[2][i][0] *= regop;
        pXHat[2][i][1] *= regop;
      }
    }
  }
}  // pragma omp parallel
  
  PetscFunctionReturn(ierr);
}*/

template PetscErrorCode EvaluateGradientKernel<1>(ComplexType*, ComplexType*, ComplexType*, IntType[3], IntType[3], IntType[3], ScalarType, ScalarType);  
template PetscErrorCode EvaluateGradientKernel<2>(ComplexType*, ComplexType*, ComplexType*, IntType[3], IntType[3], IntType[3], ScalarType, ScalarType);  
template PetscErrorCode EvaluateGradientKernel<3>(ComplexType*, ComplexType*, ComplexType*, IntType[3], IntType[3], IntType[3], ScalarType, ScalarType); 
    
PetscErrorCode ScaleVectorField(ScalarType *vR1, ScalarType *vR2, ScalarType *vR3,
    ScalarType *v1, ScalarType *v2, ScalarType *v3,
    IntType nl, ScalarType alpha) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  for (IntType i = 0; i < nl; ++i) {
    vR1[i] = alpha*v1[i];
    vR2[i] = alpha*v2[i];
    vR3[i] = alpha*v3[i];
  }
  
  PetscFunctionReturn(ierr);
}

}
