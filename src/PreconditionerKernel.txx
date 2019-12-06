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

struct H0Kernel {
  template<typename T>
  KernelOperator (int i, T* vhat[3], const T* gmt[3], const T* rhs[3], const T* reg[3], T omg) {
    array3_t<T> res;
    T gmt00 = gmt[0][i]*gmt[0][i];
    T gmt01 = gmt[0][i]*gmt[1][i];
    T gmt02 = gmt[0][i]*gmt[2][i];
    T gmt11 = gmt[1][i]*gmt[1][i];
    T gmt12 = gmt[1][i]*gmt[2][i];
    T gmt22 = gmt[2][i]*gmt[2][i];
    res.x = rhs[0][i] - reg[0][i] - gmt00*vhat[0][i] - gmt01*vhat[1][i] - gmt02*vhat[2][i];
    res.y = rhs[1][i] - reg[1][i] - gmt01*vhat[0][i] - gmt11*vhat[1][i] - gmt12*vhat[2][i];
    res.z = rhs[2][i] - reg[2][i] - gmt02*vhat[0][i] - gmt12*vhat[1][i] - gmt22*vhat[2][i];
    
    vhat[0][i] += omg*res.x;
    vhat[1][i] += omg*res.y;
    vhat[2][i] += omg*res.z;
  }
  template<typename T> KernelOperator (int i, T* diag[3], const T* gmt[3], T diag_reg) {
    array3_t<T> d;
    d.x = gmt[0][i]*gmt[0][i] + diag_reg;
    //d.x = diag_reg;
    d.y = gmt[1][i]*gmt[1][i] + diag_reg;
    //d.y = diag_reg;
    d.z = gmt[2][i]*gmt[2][i] + diag_reg;
    //d.z = diag_reg;
    if (d.x == 0.) d.x = 1.;
    if (d.y == 0.) d.y = 1.;
    if (d.z == 0.) d.z = 1.;
    diag[0][i] = 1./d.x;
    diag[1][i] = 1./d.y;
    diag[2][i] = 1./d.z;
  }
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
  template<typename T> KernelOperator (int i, T* vhat[3], T* rmv[3], T omg, T beta) {
	  //vhat[0][i] += omg*(rhs[0][i] - (vhat[0][i] + beta*rmv[0][i]));
	  vhat[0][i] += omg*(beta*rmv[0][i] - vhat[0][i]);
	  vhat[1][i] += omg*(beta*rmv[1][i] - vhat[1][i]);
	  vhat[2][i] += omg*(beta*rmv[2][i] - vhat[2][i]);
  }
};

#endif // _PRECONDITIONERKERNEL_TXX_
