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

#ifndef _DISTANCEMEASUREKERNEL_HPP_
#define _DISTANCEMEASUREKERNEL_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"

namespace reg {
namespace DistanceMeasureKernel {
  
struct EvaluateFunctionalSL2 {
  const ScalarType *pW;
  const ScalarType *pM;
  const ScalarType *pMr;
  const ScalarType *pWts;
  
  IntType nl;
  IntType nc;
  
  ScalarType value;
  
  PetscErrorCode ComputeFunctional();
  PetscErrorCode ComputeFunctionalMask();
};

struct FinalConditionSL2 {
  ScalarType *pL;
  const ScalarType *pM;
  const ScalarType *pMr;
  const ScalarType *pW;
  const ScalarType *pWts;
  
  IntType nl;
  IntType nc;
  
  PetscErrorCode ComputeFinalConditionAE();
  PetscErrorCode ComputeFinalConditionMaskAE();
  PetscErrorCode ComputeFinalConditionIAE();
  PetscErrorCode ComputeFinalConditionMaskIAE();
};

} // namespace DistanceMeasureKernel
} // namespace reg

#endif
