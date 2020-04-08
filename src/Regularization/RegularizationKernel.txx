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
 
#ifndef _REGULARIZATIONKERNEL_TXX_
#define _REGULARIZATIONKERNEL_TXX_

#include "KernelUtils.hpp"

struct NormKernel {
  template<typename T> ReductionFunctional (int i, const T* x, const T* y, const T* z) {
    return x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
  }
};

#endif // _REGULARIZATIONKERNEL_TXX_
