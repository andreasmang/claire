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
struct int3 { int x, y, z; };
#else
#define __kernel__ __host__ __device__
#endif

/********************************************************************
 * @brief computes wave number
 *******************************************************************/
__kernel__ inline void ComputeWaveNumber(int3 &w, int3 n) {
  if (w.x > n.x/2) w.x -= n.x;
  else if (w.x == n.x/2) w.x = 0;
  if (w.y > n.y/2) w.y -= n.y;
  else if (w.y == n.y/2) w.y = 0;
  if (w.z > n.z/2) w.z -= n.z;
  else if (w.z == n.z/2) w.z = 0;
}

/********************************************************************
 * @brief computes linear array index in Z-Y-X layout
 *******************************************************************/
__kernel__ inline int GetLinearIndex(int i1, int i2, int i3, int3 size) {
  return i1*size.y*size.z + i2*size.z + i3;
}

/********************************************************************
 * @brief pre-processor evaluated power
 *******************************************************************/
template<int N> __kernel__ inline ScalarType pow(ScalarType x) {
  return x*pow<N-1>(x);
}
template<> __kernel__ inline ScalarType pow<0>(ScalarType x) {
  return 1.0;
}

/********************************************************************
 * @brief computes absolute N-laplacian operator
 *******************************************************************/
/// TODO Note: first cast w, then multiply might be faster. Check accuracy there
__kernel__ inline ScalarType LaplaceNumber(int3 w) {
  return -static_cast<ScalarType>(w.x*w.x + w.y*w.y + w.z*w.z);
}

template<int N> struct NLaplacianKernel {
  static __kernel__ inline void call(int i, int3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType b0) {
    ScalarType regop = b0*pow<N>(-LaplaceNumber(w));

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct RelaxedNLaplacianKernel {
  static __kernel__ inline void call(int i, int3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType b0, ScalarType b1) {
    ScalarType regop = b0*pow<N>(-LaplaceNumber(w) + b1);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct InverseNLaplacianKernel {
  static __kernel__ inline void call (int i, int3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0) {
    ScalarType lapik = pow<N>(-LaplaceNumber(w));
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct InverseNLaplacianSqrtKernel {
  static __kernel__ inline void call(int i, int3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3,
      ScalarType scale, ScalarType b0) {
    ScalarType lapik = pow<N>(-LaplaceNumber(w));
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/sqrt(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct RelaxedInverseNLaplacianKernel {
  static __kernel__ inline void call(int i, int3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0, ScalarType b1) {
    ScalarType lapik = pow<N>(-LaplaceNumber(w) + b1);
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct RelaxedInverseNLaplacianSqrtKernel {
  static __kernel__ inline void call(int i, int3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0, ScalarType b1) {
    ScalarType lapik = pow<N>(-LaplaceNumber(w) + b1);
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/sqrt(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};

// TODO: This is not correct! This is not grad(lap(v)), but diag(grad)*lap(v)
struct TriLaplacianFunctionalKernel {
  static __kernel__ inline void call(int i, int3 w,
      ComplexType *v1, ComplexType *v2, ComplexType *v3,
      ScalarType b0) {
    ScalarType regop, gradik, lapik, tmp;
    
    // compute laplacian operator
    lapik = LaplaceNumber(w);

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

struct LerayKernel {
  static __kernel__ inline void call(int i, int3 w,
      ComplexType *v1, ComplexType *v2, ComplexType *v3,
      ScalarType scale, ScalarType b0, ScalarType b1) {
    // compute inverse laplacian operator
    ScalarType lapik = LaplaceNumber(w);

    //lapinvik = round(lapinvik) == 0.0 ? -1.0 : 1.0/lapinvik;
    ScalarType lapinvik = (lapik == 0.0 ? -1.0 : 1.0/lapik);

    // compute gradient operator
    ScalarType gradik1 = static_cast<ScalarType>(w.x);
    ScalarType gradik2 = static_cast<ScalarType>(w.y);
    ScalarType gradik3 = static_cast<ScalarType>(w.z);

    // compute div(v), note v multiplies by i
    ScalarType divVre = -scale*(gradik1*v1[i][1]
                              + gradik2*v2[i][1]
                              + gradik3*v3[i][1]);

    ScalarType divVim =  scale*(gradik1*v1[i][0]
                              + gradik2*v2[i][0]
                              + gradik3*v3[i][0]);

    // compute M^{-1] = (\beta_v (\beta_w(-\ilap + 1))^{-1} + 1)^{-1}
    // TODO Note: could be numerical instable due to lapik +- 1
    ScalarType opik = 1.0/(b1*(-lapik + 1.0));
    opik = 1.0/(b0*opik + 1.0);
    
    // compute lap^{-1} div(b)
    divVre *= opik*lapinvik;
    divVim *= opik*lapinvik;

    // compute x1 gradient of lab^{-1} div(b), note grad multiplies by i
    v1[i][0] = -gradik1*divVim;
    v1[i][1] =  gradik1*divVre;

    // compute x2 gradient of lab^{-1} div(b)
    v2[i][0] = -gradik2*divVim;
    v2[i][1] =  gradik2*divVre;

    // compute x3 gradient of lab^{-1} div(b)
    v3[i][0] = -gradik3*divVim;
    v3[i][1] =  gradik3*divVre;
  }
};

#endif // _DIFFERENTIATIONKERNEL_TXX_
