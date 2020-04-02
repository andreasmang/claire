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

#include "KernelUtils.hpp"

using KernelUtils::real3;
using KernelUtils::pow;
using KernelUtils::array3_t;

template<typename FnOp> struct VecFieldKernel {
  template<typename T, typename ... Args>
  KernelOperator (int i, array3_t<T> w, array3_t<T(*)[2]> v, Args ... args) {
    T op = FnOp::call(w, args...);
    
    v.x[i][0] *= op; v.x[i][1] *= op;
    v.y[i][0] *= op; v.y[i][1] *= op;
    v.z[i][0] *= op; v.z[i][1] *= op;
  }
};

template<int N> struct NLaplacianFilterKernel {
  template<typename T>
  KernelOperator (int i, array3_t<T> w, array3_t<T(*)[2]> v, T b0, T tol) {
    T regop = b0*pow<N>(-LaplaceNumber(w));
    
    if (abs(v.x[i][0]) < tol &&  abs(v.x[i][1]) < tol) {
      v.x[i][0] = 0.0; v.x[i][1] = 0.0;
    } else {
      v.x[i][0] *= regop; v.x[i][1] *= regop;
    }
    if (abs(v.y[i][0]) < tol &&  abs(v.y[i][1]) < tol) {
      v.y[i][0] = 0.0; v.y[i][1] = 0.0;
    } else {
      v.y[i][0] *= regop; v.y[i][1] *= regop;
    }
    if (abs(v.z[i][0]) < tol &&  abs(v.z[i][1]) < tol) {
      v.z[i][0] = 0.0; v.z[i][1] = 0.0;
    } else {
      v.z[i][0] *= regop; v.z[i][1] *= regop;
    }
  }
};

template<int N> struct NLaplacianModKernel {
  KernelOperator (int i, real3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0) {
    ScalarType regop = scale*b0*pow<N>(-LaplaceNumber(w));
    if (regop == 0.0) regop = scale*b0;

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct NLaplacianKernel {
  template<typename T> KernelFunction(T) (array3_t<T> w, T b0) {
    return b0*pow<N>(-LaplaceNumber(w));
  }
  KernelOperator (int i, real3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType b0) {
    ScalarType regop = b0*pow<N>(-LaplaceNumber(w));

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
  KernelOperator (int i, real3 w, 
      ComplexType *v, ScalarType b0) {
    ScalarType regop = b0*pow<N>(-LaplaceNumber(w));

    v[i][0] *= regop; v[i][1] *= regop;
  }
};
template<int N> struct RelaxedNLaplacianKernel {
  template<typename T> KernelFunction(T) (array3_t<T> w, T b0, T b1) {
    return b0*pow<N>(-LaplaceNumber(w) + b1);
  }
  KernelOperator (int i, real3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType b0, ScalarType b1) {
    ScalarType regop = b0*pow<N>(-LaplaceNumber(w) + b1);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct InverseNLaplacianKernel {
  KernelOperator (int i, real3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType scale, ScalarType b0) {
    ScalarType lapik = pow<N>(-LaplaceNumber(w));
    ScalarType regop;
    if (lapik == 0.0) lapik = 1.0;
    regop = scale/(b0*lapik);
    //if (lapik == 0.0) regop = 0.0;
    //else regop = scale/(b0*lapik);

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct InverseNLaplacianSqrtKernel {
  KernelOperator (int i, real3 w, 
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
  KernelOperator (int i, real3 w, 
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
  KernelOperator (int i, real3 w, 
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
  KernelOperator (int i, real3 w,
      ComplexType *v1, ComplexType *v2, ComplexType *v3,
      ScalarType b0) {
    ScalarType regop, gradik, lapik, tmp;
    
    // compute laplacian operator
    lapik = LaplaceNumber(w);

    // compute gradient operator
    gradik = w.x;
    // compute regularization operator
    regop = b0*gradik*lapik;
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v1[i][0]*regop;
    v1[i][0] = -v1[i][1]*regop;
    v1[i][1] = tmp;
    
    // compute gradient operator
    gradik = w.y;
    // compute regularization operator
    regop = b0*gradik*lapik;
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v2[i][0]*regop;
    v2[i][0] = -v2[i][1]*regop;
    v2[i][1] = tmp;
    
    // compute gradient operator
    gradik = w.z;
    // compute regularization operator
    regop = b0*gradik*lapik;
    // apply regularization operator (note: gradient multiplies by i)
    tmp = v3[i][0]*regop;
    v3[i][0] = -v3[i][1]*regop;
    v3[i][1] = tmp;
  }
};

/**
 * TODO Note: The Nyquist is set zero after applying the first gradient, it should
 * be set again after second gradient
 **/
struct LerayKernel {
  KernelOperator (int i, real3 w,
      ComplexType *v1, ComplexType *v2, ComplexType *v3,
      ScalarType scale, ScalarType b0, ScalarType b1) {
    // compute inverse laplacian operator
    ScalarType lapik = LaplaceNumber(w);

    ScalarType lapinvik = (lapik == 0.0 ? -1.0 : 1.0/lapik);

    // compute gradient operator
    ScalarType gradik1 = w.x;
    ScalarType gradik2 = w.y;
    ScalarType gradik3 = w.z;

    // compute div(v), note v multiplies by i
    ScalarType divVre = -scale*(gradik1*v1[i][1]
                              + gradik2*v2[i][1]
                              + gradik3*v3[i][1]);

    ScalarType divVim = scale*(gradik1*v1[i][0]
                             + gradik2*v2[i][0]
                             + gradik3*v3[i][0]);

    // compute M^{-1] = (\beta_v (\beta_w(-\ilap + 1))^{-1} + 1)^{-1}
    // TODO Note: could be numerical instable due to lapik +- 1
    ScalarType opik = 1.0/(b1*(-lapik + 1.0));
    opik = 1.0/(b0*opik + 1.0);
    
    // compute lap^{-1} div(b)
    divVre *= -opik*lapinvik;
    divVim *= -opik*lapinvik;

    // compute x1 gradient of lab^{-1} div(b), note grad multiplies by i
    v1[i][0] = scale*v1[i][0] - gradik1*divVim;
    v1[i][1] = scale*v1[i][1] + gradik1*divVre;

    // compute x2 gradient of lab^{-1} div(b)
    v2[i][0] = scale*v2[i][0] - gradik2*divVim;
    v2[i][1] = scale*v2[i][1] + gradik2*divVre;

    // compute x3 gradient of lab^{-1} div(b)
    v3[i][0] = scale*v3[i][0] - gradik3*divVim;
    v3[i][1] = scale*v3[i][1] + gradik3*divVre;
  }
};

struct GaussianFilterKernel {
  KernelOperator(int i, real3 w,
      ComplexType *x, ScalarType c0, ScalarType c1, ScalarType c2, ScalarType scale) {
    ScalarType sik = 0.5*( (w.x*w.x*c0) + (w.y*w.y*c1) + (w.z*w.z*c2) );
    sik = exp(-sik);

    x[i][0] *= scale*sik;
    x[i][1] *= scale*sik;
  }
};

struct GradientKernel {
  KernelOperator(int i, real3 w,
      ComplexType *x1, ComplexType *x2, ComplexType *x3, ScalarType scale) {
    ScalarType mRe, mIm;
    
    mRe = x1[i][0]*scale;
    mIm = x1[i][1]*scale;
    
    x1[i][0] = -w.x*mIm;
    x1[i][1] =  w.x*mRe;
    x2[i][0] = -w.y*mIm;
    x2[i][1] =  w.y*mRe;
    x3[i][0] = -w.z*mIm;
    x3[i][1] =  w.z*mRe;
  }
};

struct DivergenceKernel {
  KernelOperator(int i, real3 w,
      ComplexType *x1, ComplexType *x2, ComplexType *x3, ScalarType scale) {
    ScalarType mRe1, mIm1;
    ScalarType mRe2, mIm2;
    ScalarType mRe3, mIm3;
    
    mRe1 = x1[i][0]*scale;
    mIm1 = x1[i][1]*scale;
    mRe2 = x2[i][0]*scale;
    mIm2 = x2[i][1]*scale;
    mRe3 = x3[i][0]*scale;
    mIm3 = x3[i][1]*scale;
    
    x1[i][0] = -w.x*mIm1 -w.y*mIm2 -w.z*mIm3;
    x1[i][1] =  w.x*mRe1 +w.y*mRe2 +w.z*mRe3;
  }
};

#endif // _DIFFERENTIATIONKERNEL_TXX_
