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
  PetscErrorCode ComputeBodyForcePart2();
};

struct TransportKernelIncStateSL {
  ScalarType *pGm[3];
  ScalarType *pMtilde;
  const ScalarType *pVtilde[3];
  const ScalarType *pVtildex[3];
  
  ScalarType hthalf;
  
  IntType nl;
  
  PetscErrorCode TimeIntegrationPart1();
  PetscErrorCode TimeIntegrationPart2();
};

struct TransportKernelIncAdjointGN {
  ScalarType *pLtilde;
  ScalarType *pGm[3];
  ScalarType *pBtilde[3];
  
  ScalarType scale;
  
  IntType nl;
  
  PetscErrorCode ComputeBodyForce();
};

}

#endif
