/*************************************************************************
 *  Copyright (c) 2018.
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
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _TRANSPORTKERNEL_HPP_
#define _TRANSPORTKERNEL_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"

namespace reg {

struct TransportKernelAdjointSL {
  ScalarType *pB[3];
  ScalarType *pGm[3];
  ScalarType *pDivV;
  ScalarType *pDivVx;
  ScalarType *pL;
  ScalarType *pLnext;
  ScalarType *pLx;
  
  ScalarType scale;
  ScalarType ht;
  
  IntType nl;
  
  PetscErrorCode ComputeBodyForcePart1();
  PetscErrorCode ComputeBodyForcePart1b();
  PetscErrorCode ComputeDiv();
  PetscErrorCode ComputeBodyForcePart2();
  PetscErrorCode ComputeBodyForcePart0();
};

struct TransportKernelIncStateSL {
  ScalarType *pGm[3];
  ScalarType *pGmx[3];
  ScalarType *pMtilde;
  const ScalarType *pVtilde[3];
  const ScalarType *pVtildex[3];
  
  ScalarType hthalf;
  
  IntType nl;
  
  PetscErrorCode TimeIntegrationPart1();
  PetscErrorCode TimeIntegrationPart2();
  PetscErrorCode TimeIntegrationAll();
};

struct TransportKernelAdjoint {
  ScalarType *pL;
  ScalarType *pGm[3];
  ScalarType *pB[3];
  
  ScalarType scale;
  
  IntType nl;
  
  PetscErrorCode ComputeBodyForce();
};

struct TransportKernelAdjointRK2 {
  ScalarType *pVec[3];
  const ScalarType *pV[3];
  ScalarType *pRHS[2];
  ScalarType *pL;
  ScalarType *pLnext;
  ScalarType *pB[3];
  
  ScalarType scale;
  ScalarType ht;
  
  IntType nl;
  
  PetscErrorCode TimeIntegrationPart1();
  PetscErrorCode TimeIntegrationPart2();
  PetscErrorCode TimeIntegrationPart3();
  PetscErrorCode TimeIntegrationPart4();
};

struct TransportKernelIncAdjointRK2 {
  const ScalarType *pVx[3];
  const ScalarType *pVtx[3];
  const ScalarType *pL;
  ScalarType *pLt;
  ScalarType *pLtnext;
  ScalarType *pRHS[2];
  ScalarType *pLtjVx[3];
  
  ScalarType ht;
  
  IntType nl;
  
  PetscErrorCode TimeIntegrationPart1a();
  PetscErrorCode TimeIntegrationPart1b();
  PetscErrorCode TimeIntegrationPart2a();
  PetscErrorCode TimeIntegrationPart2b();
  PetscErrorCode TimeIntegrationPart3b();
};

struct TransportKernelStateRK2 {
  ScalarType *pMbar;
  ScalarType *pGmx[3];
  ScalarType *pM;
  ScalarType *pMnext;
  ScalarType *pRHS;
  const ScalarType *pVx[3];
  
  ScalarType ht;
  
  IntType nl;
  
  PetscErrorCode TimeIntegrationPart1();
  PetscErrorCode TimeIntegrationPart2();
};

struct TransportKernelIncStateRK2 {
  ScalarType *pMtnext;
  ScalarType *pMt;
  ScalarType *pRHS;
  ScalarType *pMtbar;
  ScalarType *pGmx[3];
  ScalarType *pGmtx[3];
  const ScalarType *pVx[3];
  const ScalarType *pVtx[3];
  
  ScalarType ht;
  
  IntType nl;
  
  PetscErrorCode TimeIntegrationEuler();
  PetscErrorCode TimeIntegrationPart1();
  PetscErrorCode TimeIntegrationPart2();
};

struct TransportKernelContinuity {
  ScalarType *pMx;
  ScalarType *pDivV;
  ScalarType *pDivVx;
  ScalarType *pMnext;
  
  ScalarType ht;
  
  IntType nl;
  
  PetscErrorCode TimeIntegration();
};

template<typename T>
PetscErrorCode TransportKernelCopy(T*, T*, IntType);

}

#endif
