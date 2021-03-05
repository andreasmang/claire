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
 
#ifndef _PRECONDITIONERKERNEL_TXX_
#define _PRECONDITIONERKERNEL_TXX_

#include "KernelUtils.hpp"

using KernelUtils::array3_t;


struct H0Kernel2 {
  // computes: grad M \otimes grad M
  template<typename T> KernelOperator (int i, T* mvx, T* mvy, T* mvz, 
                                       const T* rx, const T* ry, const T* rz,
                                       const T* gmtx, const T* gmty, const T* gmtz) {
    T gmt00 = gmtx[i]*gmtx[i];
    T gmt01 = gmtx[i]*gmty[i];
    T gmt02 = gmtx[i]*gmtz[i];
    T gmt11 = gmty[i]*gmty[i];
    T gmt12 = gmty[i]*gmtz[i];
    T gmt22 = gmtz[i]*gmtz[i];
    mvx[i] +=  gmt00*rx[i] + gmt01*ry[i] + gmt02*rz[i];
    mvy[i] +=  gmt01*rx[i] + gmt11*ry[i] + gmt12*rz[i];
    mvz[i] +=  gmt02*rx[i] + gmt12*ry[i] + gmt22*rz[i];
  }
  // computes residual
  template<typename T> ReductionFunctional (int i, T* mx, T* my, T* mz,
                                            T* px, T* py, T* pz,
                                            T* rx, T* ry, T* rz,
                                            const T* gmtx, const T* gmty, const T* gmtz,
                                            const T diag) {
    T d0 = gmtx[i]*gmtx[i] + diag;
    T d1 = gmty[i]*gmty[i] + diag;
    T d2 = gmtz[i]*gmtz[i] + diag;
    rx[i] -= mx[i]; ry[i] -= my[i]; rz[i] -= mz[i];
    rx[i] /= d0; ry[i] /= d1; rz[i] /= d2;
    px[i] = rx[i]; py[i] = ry[i]; pz[i] = rz[i];
    
    return rx[i]*rx[i] + ry[i]*ry[i] + rz[i]*rz[i];
  }
  // computes p^T A p
  template<typename T> ReductionFunctional (int i, T* mx, T* my, T* mz,
                                            T* px, T* py, T* pz,
                                            const T* gmtx, const T* gmty, const T* gmtz,
                                            const T diag) {
    T d0 = gmtx[i]*gmtx[i] + diag;
    T d1 = gmty[i]*gmty[i] + diag;
    T d2 = gmtz[i]*gmtz[i] + diag;
    mx[i] /= d0;
    my[i] /= d1;
    mz[i] /= d2;
    return mx[i]*px[i] + my[i]*py[i] + mz[i]*pz[i];
  }
};

struct H0Kernel {
  // computes: grad M \otimes grad M
  template<typename T> KernelOperator (int i, T* mvx, T* mvy, T* mvz, 
                                       const T* rx, const T* ry, const T* rz,
                                       const T* gmtx, const T* gmty, const T* gmtz) {
    T gmt00 = gmtx[i]*gmtx[i];
    T gmt01 = gmtx[i]*gmty[i];
    T gmt02 = gmtx[i]*gmtz[i];
    T gmt11 = gmty[i]*gmty[i];
    T gmt12 = gmty[i]*gmtz[i];
    T gmt22 = gmtz[i]*gmtz[i];
    mvx[i] =  gmt00*rx[i] + gmt01*ry[i] + gmt02*rz[i];
    mvy[i] =  gmt01*rx[i] + gmt11*ry[i] + gmt12*rz[i];
    mvz[i] =  gmt02*rx[i] + gmt12*ry[i] + gmt22*rz[i];
  }
  // computes residual
  template<typename T> ReductionFunctional (int i, T* mx, T* my, T* mz,
                                            T* px, T* py, T* pz,
                                            T* rx, T* ry, T* rz,
                                            const T* vx, const T* vy, const T* vz) {
    T H0x, H0y, H0z;
    H0x = mx[i] + vx[i]; H0y = my[i] + vy[i]; H0z = mz[i] + vz[i];
    rx[i] -= H0x; ry[i] -= H0y; rz[i] -= H0z;
    px[i] = rx[i]; py[i] = ry[i]; pz[i] = rz[i];
    
    return rx[i]*rx[i] + ry[i]*ry[i] + rz[i]*rz[i];
  }
  // computes p^T A p
  template<typename T> ReductionFunctional (int i, T* mx, T* my, T* mz,
                                            T* px, T* py, T* pz) {
    mx[i] += px[i]; my[i] += py[i]; mz[i] += pz[i];
    
    return mx[i]*px[i] + my[i]*py[i] + mz[i]*pz[i];
  }
};
struct H0KernelCG {
  // computes update for x and res
  template<typename T> ReductionFunctional (int i, T* mx, T* my, T* mz,
                                            T* px, T* py, T* pz,
                                            T* rx, T* ry, T* rz,
                                            T* vx, T* vy, T* vz,
                                            const T alpha) {
    vx[i] += alpha*px[i]; vy[i] += alpha*py[i]; vz[i] += alpha*pz[i];
    rx[i] -= alpha*mx[i]; ry[i] -= alpha*my[i]; rz[i] -= alpha*mz[i];
    
    return rx[i]*rx[i] + ry[i]*ry[i] + rz[i]*rz[i];
  }
  // computes update for p
  template<typename T> KernelOperator (int i,
                                       T* px, T* py, T* pz,
                                       T* rx, T* ry, T* rz,
                                       const T alpha) {
    px[i] = alpha*px[i] + rx[i];
    py[i] = alpha*py[i] + ry[i];
    pz[i] = alpha*pz[i] + rz[i];
  }
};

struct CFLKernel {
  template<typename T> ReductionFunctional (int i, const T* v, const T h, const T dt) {
    return ((fabs(v[i])*dt > h) ? 1.0 : 0.0);
  }
};

struct NormKernel {
  template<typename T> ReductionFunctional (int i, const T* v) {
    return v[i]*v[i];
  }
};

#endif // _PRECONDITIONERKERNEL_TXX_
