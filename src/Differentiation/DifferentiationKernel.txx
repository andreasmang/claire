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
 
#ifndef _DIFFERENTIATIONKERNEL_TXX_
#define _DIFFERENTIATIONKERNEL_TXX_

#ifndef REG_HAS_CUDA
#define __kernel__
struct {
  int x, y, z;
} dim3;
#else
#define __kernel__ __host__ __device__
#endif

/********************************************************************
 * @brief computes wave number
 *******************************************************************/
__kernel__ inline void ComputeWaveNumber(dim3 &w, dim3 n) {
  if (w.x > n.x/2) w.x -= n.x;
  if (w.y > n.y/2) w.y -= n.y;
  if (w.z > n.z/2) w.z -= n.z;
}

/********************************************************************
 * @brief computes linear array index in Z-Y-X layout
 *******************************************************************/
__kernel__ inline int GetLinearIndex(int i1, int i2, int i3, dim3 size) {
  return i1*size.y*size.z + i2*size.z + i3;
}

/********************************************************************
 * @brief pre-processor evaluated power
 *******************************************************************/
template<int N> __device__ inline ScalarType pow(ScalarType x) {
  return x*pow<N-1>(x);
}
template<> __device__ inline ScalarType pow<0>(ScalarType x) {
  return 1;
}

/********************************************************************
 * @brief computes absolute N-laplacian operator
 *******************************************************************/
/// TODO Note: first cast w, then multiply might be faster. Check accuracy there
template<int N> __kernel__ inline ScalarType AbsNLaplaceNumber(dim3 w) {
  ScalarType norm = static_cast<ScalarType>(w.x*w.x + w.y*w.y + w.z*w.z);
  return pow<N>(norm);
}

template<int N> struct NLaplacianKernel {
  static __kernel__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType b0) {
    ScalarType regop = b0*AbsNLaplaceNumber<N>(w);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct RelaxedNLaplacianKernel {
  static __kernel__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType b0, ScalarType b1) {
    ScalarType regop = b0*(AbsNLaplaceNumber<N>(w) + b1);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct InverseNLaplacianKernel {
  static __kernel__ inline void call (int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0) {
    ScalarType lapik = AbsNLaplaceNumber<N>(w);
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct InverseNLaplacianSqrtKernel {
  static __kernel__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3,
      ScalarType scale, ScalarType b0) {
    ScalarType lapik = AbsNLaplaceNumber<N>(w);
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/sqrt(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct RelaxedInverseNLaplacianKernel {
  static __kernel__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0, ScalarType b1) {
    ScalarType lapik = AbsNLaplaceNumber<N>(w) + b1;
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct RelaxedInverseNLaplacianSqrtKernel {
  static __kernel__ inline void call(int i, dim3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0, ScalarType b1) {
    ScalarType lapik = AbsNLaplaceNumber<N>(w) + b1;
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/sqrt(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};

struct RelaxedTriLaplacianFunctionalKernel {
  static __kernel__ inline void call(int i, dim3 w,
      ComplexType *v1, ComplexType *v2, ComplexType *v3,
      ScalarType b0, ScalarType b1) {
    ScalarType regop[2], gradik, lapik, tmp;
    
    // compute laplacian operator
    lapik = -AbsNLaplaceNumber<1>(w);

    // compute gradient operator
    gradik = static_cast<ScalarType>(w.x);
    // compute regularization operator
    regop[0] = b0*( gradik*lapik + b1);
    regop[1] = b0*(-gradik*lapik + b1);
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v1[i][0]*regop[0];
    v1[i][0] = v1[i][1]*regop[0];
    v1[i][1] = tmp;
    
    // compute gradient operator
    gradik = static_cast<ScalarType>(w.y);
    // compute regularization operator
    regop[0] = b0*( gradik*lapik + b1);
    regop[1] = b0*(-gradik*lapik + b1);
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v2[i][0]*regop[0];
    v2[i][0] = v2[i][1]*regop[0];
    v2[i][1] = tmp;
    
    // compute gradient operator
    gradik = static_cast<ScalarType>(w.z);
    // compute regularization operator
    regop[0] = b0*( gradik*lapik + b1);
    regop[1] = b0*(-gradik*lapik + b1);
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v3[i][0]*regop[0];
    v3[i][0] = v3[i][1]*regop[0];
    v3[i][1] = tmp;
  }
};

struct TriLaplacianFunctionalKernel {
  static __kernel__ inline void call(int i, dim3 w,
      ComplexType *v1, ComplexType *v2, ComplexType *v3,
      ScalarType b0) {
    ScalarType regop, gradik, lapik, tmp;
    
    // compute laplacian operator
    lapik = -AbsNLaplaceNumber<1>(w);

    // compute gradient operator
    gradik = static_cast<ScalarType>(w.x);
    // compute regularization operator
    regop = b0*gradik*lapik;
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v1[i][0]*regop;
    v1[i][0] = -v1[i][1]*regop;
    v1[i][1] = tmp;
    
    // compute gradient operator
    gradik = static_cast<ScalarType>(w.y);
    // compute regularization operator
    regop = b0*gradik*lapik;
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v2[i][0]*regop;
    v2[i][0] = -v2[i][1]*regop;
    v2[i][1] = tmp;
    
    // compute gradient operator
    gradik = static_cast<ScalarType>(w.z);
    // compute regularization operator
    regop = b0*gradik*lapik;
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v3[i][0]*regop;
    v3[i][0] = -v3[i][1]*regop;
    v3[i][1] = tmp;
  }
};

#endif // _DIFFERENTIATIONKERNEL_TXX_
