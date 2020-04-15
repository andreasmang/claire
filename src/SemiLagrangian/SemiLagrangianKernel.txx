/*************************************************************************
 *  Copyright (c) 2018. Malte Brunn
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
 
#ifndef _SEMILAGRANGIANKERNEL_TXX_
#define _SEMILAGRANGIANKERNEL_TXX_

#include "KernelUtils.hpp"

struct RK2Kernel {
  template<typename T> KernelOperator (int i, int3 p, 
                                       T *x, T *y, T *z, 
                                       const T *vx, const T *vy, const T *vz,
                                       const T ix, const T iy, const T iz,
                                       const T hx, const T hy, const T hz) {
    x[i] = ix*static_cast<T>(p.x) - hx*vx[i];
    y[i] = iy*static_cast<T>(p.y) - hy*vy[i];
    z[i] = iz*static_cast<T>(p.z) - hz*vz[i];
  }
  template<typename T> KernelOperator (int i, int3 p, 
                                       T *x, T *y, T *z, 
                                       const T *vx, const T *vy, const T *vz,
                                       const T *vxx, const T *vyx, const T *vzx,
                                       const T ix, const T iy, const T iz,
                                       const T hx, const T hy, const T hz) {
    x[i] = ix*static_cast<T>(p.x) - hx*(vx[i] + vxx[i]);
    y[i] = iy*static_cast<T>(p.y) - hy*(vy[i] + vyx[i]);
    z[i] = iz*static_cast<T>(p.z) - hz*(vz[i] + vzx[i]);
  }
  
  template<typename T> KernelOperator (int i, int3 p, 
                                       const T *vx, const T *vy, const T *vz,
                                       T *vxx, T *vyx, T *vzx,
                                       const T ix, const T iy, const T iz,
                                       const T hx, const T hy, const T hz) {
    vxx[i] = ix*static_cast<T>(p.x) - hx*(vx[i] + vxx[i]);
    vyx[i] = iy*static_cast<T>(p.y) - hy*(vy[i] + vyx[i]);
    vzx[i] = iz*static_cast<T>(p.z) - hz*(vz[i] + vzx[i]);
  }
};

#endif
