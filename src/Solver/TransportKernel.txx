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
 
#ifndef _TRANSPORTKERNEL_TXX_
#define _TRANSPORTKERNEL_TXX_

#include "KernelUtils.hpp"

struct AdjointSL {
  template<typename T>
  KernelOperator (int i, 
      T *pB1, T *pB2, T *pB3, T *pLnext,
      const T *pG1, const T *pG2, const T *pG3,
      const T *pL, const T *pLx, const T *pDivV, const T *pDivVx,
      T ht, T ht_2, T ht_c) {
    T lambda  = pL[i];
    T lambdax = pLx[i];

    T rhs0 = lambdax*pDivVx[i];
    T rhs1 = (lambdax + ht*rhs0)*pDivV[i];

    // compute \lambda(x,t^{j+1})
    pLnext[i] = lambdax + ht_2*(rhs0 + rhs1);

    // compute bodyforce
    pB1[i] += ht_c*pG1[i]*lambda;
    pB2[i] += ht_c*pG2[i]*lambda;
    pB3[i] += ht_c*pG3[i]*lambda;
  }
  template<typename T>
  KernelOperator (int i, T *pB1, T *pB2, T *pB3,
      const T *pG1, const T *pG2, const T *pG3,
      const T *pL, T ht_c) {
    T lambda  = pL[i];
    
    // compute bodyforce
    pB1[i] += ht_c*pG1[i]*lambda;
    pB2[i] += ht_c*pG2[i]*lambda;
    pB3[i] += ht_c*pG3[i]*lambda;
  }
};

struct ScalarRHS {
  template<typename T> KernelFunction(T) (int i, const T *pRHS) {
    return pRHS[i];
  }
};
struct VecDotRHS {
  template<typename T> KernelFunction(T) (int i, array3_t<const T*> pG, array3_t<const T*> pV) {
    return pG.x[i]*pV.x[i] + pG.y[i]*pV.y[i] + pG.z[i]*pV.z[i];
  }
};
struct DobuleVecDotRHS {
  template<typename T> KernelFunction(T) (int i,
      array3_t<const T*> pG, array3_t<const T*> pGt,
      array3_t<const T*> pV, array3_t<const T*> pVt) {
    return pGt.x[i]*pV.x[i] + pGt.y[i]*pV.y[i] + pGt.z[i]*pV.z[i]
         + pG.x[i]*pVt.x[i] + pG.y[i]*pVt.y[i] + pG.z[i]*pVt.z[i];
  }
};

template<typename FnRHS> struct Euler {
  template<typename T, typename ... Args>
  KernelOperator (int i, T *pM, T ht, Args ... args) {
    pM[i] += ht*FnRHS::call(i, args...);
  }
  template<typename T, typename ... Args>
  KernelOperator (int i, T *pMnext, const T *pM, T ht, Args ... args) {
    pMnext[i] = pM[i] + ht*FnRHS::call(i, args...);
  }
  template<typename T, typename ... Args>
  KernelOperator (int i, T *pMnext, T *pRHS, const T *pM, T ht, Args ... args) {
    T rhs = FnRHS::call(i, args...);
    pRHS[i] = rhs;
    pMnext[i] = pM[i] + ht*rhs;
  }
};
template<typename FnRHS> struct RK2 {
  template<typename T, typename ... Args>
  KernelOperator (int i, T *pMnext, const T *pM, const T *pRHS, T ht_2, Args ... args) {
    T rhs = FnRHS::call(i, args...);
    pMnext[i] = pM[i] +  ht_2*(rhs + pRHS[i]);
  }
};

#endif // _TRANSPORTKERNEL_TXX_
