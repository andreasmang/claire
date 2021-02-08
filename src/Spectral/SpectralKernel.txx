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
 
#ifndef _SPECTRALKERNEL_TXX_
#define _SPECTRALKERNEL_TXX_

#include "KernelUtils.hpp"

using KernelUtils::real3;
using KernelUtils::pow;
using KernelUtils::array3_t;

struct LowPassFilterKernel {
  KernelOperator(int i, real3 w, ComplexType *x,
      ScalarType l1, ScalarType l2, ScalarType l3, ScalarType scale) {
    if (!(abs(w.x) < l1 && abs(w.y) < l2 && abs(w.z) < l3))
      scale = 0.;
    
    x[i][0] *= scale;
    x[i][1] *= scale;
  }
};

struct HighPassFilterKernel {
  KernelOperator(int i, real3 w, ComplexType *x,
      ScalarType l1, ScalarType l2, ScalarType l3, ScalarType scale) {
    if (abs(w.x) < l1 && abs(w.y) < l2 && abs(w.z) < l3)
      scale = 0.;
    
    x[i][0] *= scale;
    x[i][1] *= scale;
  }
};

struct ScaleKernel {
  KernelOperator(int i, real3 w, ComplexType *x, ScalarType scale) {
    x[i][0] *= scale;
    x[i][1] *= scale;
  }
};

struct NormKernel {
  ReductionFunctional (int i, real3 w, ComplexType *x, ScalarType l1, ScalarType l2, ScalarType l3) {
    if (abs(w.x) < l1 && abs(w.y) < l2 && abs(w.z) < l3) 
      return x[i][0]*x[i][0] + x[i][1]*x[i][1];
    else 
      return 0.;
  }
};

#endif // _SPECTRAL_KERNEL_
