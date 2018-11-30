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

template<int N> struct NLaplacianFilterKernel {
  KernelOperator (int i, real3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType b0, ScalarType tol) {
    ScalarType regop = b0*pow<N>(-LaplaceNumber(w));
    
    if (abs(v1[i][0]) < tol &&  abs(v1[i][1]) < tol) {
      v1[i][0] = 0.0; v1[i][1] = 0.0;
    } else {
      v1[i][0] *= regop; v1[i][1] *= regop;
    }
    if (abs(v2[i][0]) < tol &&  abs(v2[i][1]) < tol) {
      v2[i][0] = 0.0; v2[i][1] = 0.0;
    } else {
      v2[i][0] *= regop; v2[i][1] *= regop;
    }
    if (abs(v3[i][0]) < tol &&  abs(v3[i][1]) < tol) {
      v3[i][0] = 0.0; v3[i][1] = 0.0;
    } else {
      v3[i][0] *= regop; v3[i][1] *= regop;
    }
  }
};

template<int N> struct NLaplacianKernel {
  KernelOperator (int i, real3 w, 
      ComplexType *v1, ComplexType *v2, ComplexType *v3, 
      ScalarType b0) {
    ScalarType regop = b0*pow<N>(-LaplaceNumber(w));

    v1[i][0] *= regop; v1[i][1] *= regop;
    v2[i][0] *= regop; v2[i][1] *= regop;
    v3[i][0] *= regop; v3[i][1] *= regop;
  }
};
template<int N> struct RelaxedNLaplacianKernel {
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
    if (lapik == 0.0) lapik = 1.0;
    ScalarType regop = scale/(b0*lapik);

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
