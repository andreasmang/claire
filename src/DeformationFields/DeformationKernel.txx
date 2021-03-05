/*************************************************************************
 *  Copyright (c) 2019.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 * 
 *  Author: Malte Brunn
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
 
#ifndef _DEFORMATIONKERNEL_TXX_
#define _DEFORMATIONKERNEL_TXX_

#include "KernelUtils.hpp"

using KernelUtils::real3;
using KernelUtils::array3_t;

template<typename T> using Tp = T*;
template<typename T> using cTp = const T*;
template<typename T> using T3 = array3_t<T>;
template<typename T> using T3p = array3_t<T*>;
template<typename T> using cT3p = array3_t<const T*>;

template<typename FnRHS> struct Euler {
  template<typename T, typename ... Args>
  KernelOperator (int i, Tp<T> pM, T ht, Args ... args) {
    pM[i] += ht*FnRHS::call(i, args...);
  }
  template<typename T, typename ... Args>
  KernelOperator (int i, Tp<T> pMnext, cTp<T> pM, T ht, Args ... args) {
    pMnext[i] = pM[i] + ht*FnRHS::call(i, args...);
  }
};
template<typename FnRHS> struct RK2_1 {
  template<typename T, typename ... Args>
  KernelOperator (int i, Tp<T> pMnext, Tp<T> pRHS, cTp<T> pM, T ht, Args ... args) {
    T rhs = FnRHS::call(i, args...);
    pRHS[i] = rhs;
    pMnext[i] = pM[i] + ht*rhs;
  }
  template<typename T, typename ... Args>
  KernelOperator (int i, Tp<T> pM, Tp<T> pRHS, T ht, Args ... args) {
    T rhs = FnRHS::call(i, args...);
    pRHS[i] = rhs;
    pM[i] += ht*rhs;
  }
};
template<typename FnRHS> struct RK2_2 {
  template<typename T, typename ... Args>
  KernelOperator (int i, Tp<T> pM, cTp<T> pRHS, T ht_2, Args ... args) {
    T rhs = FnRHS::call(i, args...);
    pM[i] += ht_2*(rhs + pRHS[i]);
  }
  template<typename T, typename ... Args>
  KernelOperator (int i, Tp<T> pMnext, cTp<T> pM, cTp<T> pRHS, T ht_2, Args ... args) {
    T rhs = FnRHS::call(i, args...);
    pMnext[i] = pM[i] +  ht_2*(rhs + pRHS[i]);
  }
  
};

struct VecDotFunc {
  template<typename T> KernelFunction(T) (int i, cT3p<T> pG, cT3p<T> pV) {
    return pG.x[i]*pV.x[i] + pG.y[i]*pV.y[i] + pG.z[i]*pV.z[i];
  }
  template<typename T> KernelFunction(T) (int i, cT3p<T> pG, T3<T> V) {
    return pG.x[i]*V.x + pG.y[i]*V.y + pG.z[i]*V.z;
  }
};

struct DetDefGrad_RHS {
  template<typename T>
  KernelFunction(T) (int i,  cTp<T> pJac, cTp<T> pDivV, cT3p<T> pV, cT3p<T> pG, T alpha) {
    return alpha*pJac[i]*pDivV[i] - VecDotFunc::call(i, pV, pG);
  }
};

struct DetDefGradA_RHS {
  template<typename T> 
      KernelFunction(T) (int i, cTp<T> pPhi, cTp<T> pDivV, cTp<T> pDivPhi, cT3p<T> pV, cT3p<T> pGPhi, T alpha) {
    return -0.5*VecDotFunc::call(i, pV, pGPhi) + 
      alpha*pPhi[i]*pDivV[i] - 0.5*pDivPhi[i] + 0.5*pPhi[i]*pDivV[i];
  }
};
struct ScaleVectorKernel {
  template<typename T> KernelOperator (int i, T3p<T> pPhiV, cTp<T> pPhi, cT3p<T> pV) {
    pPhiV.x[i] = pPhi[i]*pV.x[i];
    pPhiV.y[i] = pPhi[i]*pV.y[i];
    pPhiV.z[i] = pPhi[i]*pV.z[i];
  }
};

struct DetDefGradSLKernel {
  template<typename T> KernelOperator (int i, Tp<T> pJ, T v) {
    pJ[i] = v;
  }
  template<typename T> KernelOperator (int i, Tp<T> pJ, 
      cTp<T> pJx, cTp<T> pDivV, cTp<T> pDivVx,
      T alpha, T ht) {
    T jX = pJx[i];
    T rhs0 = alpha*jX*pDivVx[i];
    T rhs1 = alpha*(jX + ht*rhs0)*pDivV[i];
    pJ[i] = jX + 0.5*ht*(rhs0 + rhs1);
  }
};

struct DetDefGradViaDispFieldKernel {
  template<typename T> KernelOperator (int i, Tp<T> pPhi,
      cT3p<T> pGu1, cT3p<T> pGu2, cT3p<T> pGu3) {
    T pGu11 = 1.0 - pGu1.x[i];
    T pGu12 = pGu1.y[i];
    T pGu13 = pGu1.z[i];
    T pGu21 = pGu2.x[i];
    T pGu22 = 1.0 - pGu2.y[i];
    T pGu23 = pGu2.z[i];
    T pGu31 = pGu3.x[i];
    T pGu32 = pGu3.y[i];
    T pGu33 = 1.0 - pGu3.z[i];
    pPhi[i] = pGu11*pGu22*pGu33 + pGu12*pGu23*pGu31 + pGu13*pGu21*pGu32
            - pGu13*pGu22*pGu31 - pGu12*pGu21*pGu33 - pGu11*pGu23*pGu32;
  }
};

struct DefGradSLKernel {
  template<typename T> KernelOperator (int i, 
      T3p<T> pJ1, T3p<T> pJ2, T3p<T> pJ3,
      cT3p<T> pG1T, cT3p<T> pG2T, cT3p<T> pG3T,
      cT3p<T> pGx1T, cT3p<T> pGx2T, cT3p<T> pGx3T,
      cT3p<T> pJx1T, cT3p<T> pJx2T, cT3p<T> pJx3T,
      T ht) {
    T3<T> R1T, R2T, R3T;
    T ht_2 = 0.5*ht;
    R1T.x = VecDotFunc::call(i, pGx1T, pJx1T);
    R1T.y = VecDotFunc::call(i, pGx2T, pJx1T);
    R1T.z = VecDotFunc::call(i, pGx3T, pJx1T);
    R2T.x = VecDotFunc::call(i, pGx1T, pJx2T);
    R2T.y = VecDotFunc::call(i, pGx2T, pJx2T);
    R2T.z = VecDotFunc::call(i, pGx3T, pJx2T);
    R3T.x = VecDotFunc::call(i, pGx1T, pJx3T);
    R3T.y = VecDotFunc::call(i, pGx2T, pJx3T);
    R3T.z = VecDotFunc::call(i, pGx3T, pJx3T);
    pJ1.x[i] = pJx1T.x[i] + ht_2*(R1T.x[i] + pGx1T.x[i] + ht*VecDotFunc::call(i, pG1T, R1T));
    pJ1.y[i] = pJx2T.x[i] + ht_2*(R2T.x[i] + pGx2T.x[i] + ht*VecDotFunc::call(i, pG1T, R2T));
    pJ1.z[i] = pJx3T.x[i] + ht_2*(R3T.x[i] + pGx3T.x[i] + ht*VecDotFunc::call(i, pG1T, R3T));
    pJ2.x[i] = pJx1T.y[i] + ht_2*(R1T.y[i] + pGx1T.y[i] + ht*VecDotFunc::call(i, pG2T, R1T));
    pJ2.y[i] = pJx2T.y[i] + ht_2*(R2T.y[i] + pGx2T.y[i] + ht*VecDotFunc::call(i, pG2T, R2T));
    pJ2.z[i] = pJx3T.y[i] + ht_2*(R3T.y[i] + pGx3T.y[i] + ht*VecDotFunc::call(i, pG2T, R3T));
    pJ3.x[i] = pJx1T.z[i] + ht_2*(R1T.z[i] + pGx1T.z[i] + ht*VecDotFunc::call(i, pG3T, R1T));
    pJ3.y[i] = pJx2T.z[i] + ht_2*(R2T.z[i] + pGx2T.z[i] + ht*VecDotFunc::call(i, pG3T, R2T));
    pJ3.z[i] = pJx3T.z[i] + ht_2*(R3T.z[i] + pGx3T.z[i] + ht*VecDotFunc::call(i, pG3T, R3T));
  }
};

#endif // _DEFORMATIONKERNEL_TXX_
