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

#include "DistanceMeasureKernel.hpp"

namespace reg {
namespace DistanceMeasureKernel {
  
////////////////////////////////////////////////////////////////////////
//> SSD Distance metric routines 
///////////////////////////////////////////////////////////////////////
PetscErrorCode EvaluateFunctionalSL2::ComputeFunctionalMask() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  value = 0.0;
  for (IntType k = 0; k < nc; ++k) {  // for all image components
      for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
          ScalarType dr = (pMr[k*nl+i] - pM[k*nl+i]);
          value += pW[i]*dr*dr;
      }
  }
  
  PetscFunctionReturn(ierr);
}
  
PetscErrorCode EvaluateFunctionalSL2::ComputeFunctional() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  value = 0.0;
  for (IntType i = 0; i < nc*nl; ++i) {
      ScalarType dr = (pMr[i] - pM[i]);
      value += dr*dr;
  }
  
  PetscFunctionReturn(ierr);
}


PetscErrorCode FinalConditionSL2::ComputeFinalConditionAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nc*nl; ++i) {
    pL[i] = pMr[i] - pM[i];
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionSL2::ComputeFinalConditionMaskAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType k = 0; k < nc; ++k) {  // for all image components
    for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
      pL[k*nl+i] = pW[i]*(pMr[k*nl+i] - pM[k*nl+i]);
    }
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}


PetscErrorCode FinalConditionSL2::ComputeFinalConditionIAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl*nc; ++i) {
    pL[i] = -pM[i]; // / static_cast<ScalarType>(nc);
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}


PetscErrorCode FinalConditionSL2::ComputeFinalConditionMaskIAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType k = 0; k < nc; ++k) {  // for all image components
    for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
      pL[k*nl+i] = -pW[i]*pM[k*nl+i];
    }
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}


////////////////////////////////////////////////////////////////////////
//> NCC Distance metric routines 
///////////////////////////////////////////////////////////////////////
PetscErrorCode EvaluateFunctionalNCC::ComputeFunctionalMask() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
    
  norm_m1_loc = 0;
  norm_mR_loc = 0;
  inpr_m1_mR_loc = 0;

  for (IntType k = 0; k < nc; ++k) {  // for all image components
      for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
          m1i = pM[k*nl+i];
          mRi = pMr[k*nl+i];
          norm_m1_loc    += pW[i]*(m1i*m1i);
          norm_mR_loc    += pW[i]*(mRi*mRi);
          inpr_m1_mR_loc += pW[i]*(m1i*mRi);
      }
  }
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode EvaluateFunctionalNCC::ComputeFunctional() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
    
  norm_m1_loc = 0;
  norm_mR_loc = 0;
  inpr_m1_mR_loc = 0;

  for (IntType i = 0; i < nc*nl; ++i) {
        m1i = pM[i];
        mRi = pMr[i];
        norm_m1_loc    += (m1i*m1i);
        norm_mR_loc    += (mRi*mRi);
        inpr_m1_mR_loc += (m1i*mRi);
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionNCC::ComputeFinalConditionAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nc*nl; ++i) {
        pL[i] = const1*pMr[i] - const2*pM[i];
  }
} // omp
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionNCC::ComputeFinalConditionMaskAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType k = 0; k < nc; ++k) {  // for all image components
    for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
        pL[k*nl+i] = pW[i]*(const1*pMr[k*nl+i] - const2*pM[k*nli]);
  }
} // omp
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionNCC::ComputeFinalConditionIAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl*nc; ++i) {
        ScalarType m1i = pM[i];
        ScalarType mRi = pMr[i];
        ScalarType mtildei = pMtilde[i];
        pLtilde[i] = const1*mRi - const2*mRi + const3*m1i - const4*m1i - const5*mtildei;
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionNCC::ComputeFinalConditionMaskIAE() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType k = 0; k < nc; ++k) {  // for all image components
    for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
        ScalarType m1i = pM[i+k*nl];
        ScalarType mRi = pMr[i+k*nl];
        ScalarType mtildei = pMtilde[i+k*nl];
        pLtilde[i+k*nl] = pW[i]*(const1*mRi - const2*mRi + const3*m1i - const4*m1i - const5*mtilde1i);
    }
  }
}  // omp
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode FinalConditionNCC::ComputeInnerProductsFinalConditionIAE() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    inpr_m1_mtilde_loc = 0;
    inpr_mR_mtilde_loc = 0;
#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i<nl*nl; ++i) {
        m1i = pM[i];
        mRi = pMr[i];
        mtildei = pMtilde[i];
        inpr_m1_mtilde_loc += (m1i*mtilde1i);
        inpr_mR_mtilde_loc += (mRi*mtilde1i);
    }
}
    PetscFunctionReturn(ierr);
}

} // namespace DistanceMeasureKernel
} // namespace reg
