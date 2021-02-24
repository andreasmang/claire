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
  
  ScalarType *res;
  
  IntType nl;
  IntType nc;
  
  ScalarType value;
  
  PetscErrorCode ComputeFunctional();
  PetscErrorCode ComputeFunctionalMask();
};

IntType GetTempResSize(IntType nl, IntType nc);

struct EvaluateFunctionalNCC {
  const ScalarType *pW;
  const ScalarType *pM;
  const ScalarType *pMr;
  const ScalarType *pMt;
  const ScalarType *pWts;
  ScalarType norm_l2_loc;
  ScalarType norm_m1_loc;
  ScalarType norm_mT_loc;
  ScalarType norm_mR_loc;
  ScalarType inpr_m1_mR_loc;
  ScalarType inpr_mT_mR_loc;
  ScalarType sum_m1_loc;
  ScalarType sum_mT_loc;
  ScalarType sum_mR_loc;
  
  IntType nl;
  IntType nc;
  
  ScalarType value;
  
  PetscErrorCode ComputeFunctional();
  PetscErrorCode ComputeScale();
  PetscErrorCode ComputeFunctionalMask();
  PetscErrorCode ComputeScaleMask();
};

struct FinalConditionSL2 {
  ScalarType *pL;
  const ScalarType *pM;
  const ScalarType *pMr;
  const ScalarType *pW;
  const ScalarType *pWts;
  ScalarType norm_mtilde_loc;
  
  IntType nl;
  IntType nc;
  
  PetscErrorCode ComputeFinalConditionAE();
  PetscErrorCode ComputeFinalConditionMaskAE();
  PetscErrorCode ComputeFinalConditionIAE();
  PetscErrorCode ComputeFinalConditionMaskIAE();
};

struct FinalConditionNCC {
  ScalarType *pL;
  ScalarType *pLtilde;
  const ScalarType *pMtilde;
  const ScalarType *pM;
  const ScalarType *pMr;
  const ScalarType *pW;
  const ScalarType *pWts;
  ScalarType const1;
  ScalarType const2;
  ScalarType const3;
  ScalarType const4;
  ScalarType const5;
  ScalarType norm_m1_loc;
  ScalarType norm_mR_loc;
  ScalarType norm_mtilde_loc;
  ScalarType inpr_m1_mR_loc;
  ScalarType inpr_m1_mtilde_loc;
  ScalarType inpr_mR_mtilde_loc;
  ScalarType sum_m1_loc;
  ScalarType sum_mR_loc;
  ScalarType sum_mtilde_loc;
  ScalarType mean_m1;
  ScalarType mean_mR;
  
  
  IntType nl;
  IntType nc;
  
  PetscErrorCode ComputeFinalConditionAE();
  PetscErrorCode ComputeFinalConditionMaskAE();
  PetscErrorCode ComputeFinalConditionIAE();
  PetscErrorCode ComputeFinalConditionMaskIAE();
  PetscErrorCode ComputeInnerProductsFinalConditionAE();
  PetscErrorCode ComputeInnerProductsFinalConditionIAE();
};

} // namespace DistanceMeasureKernel
} // namespace reg

#endif
