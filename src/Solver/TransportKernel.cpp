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

#include "TransportKernel.hpp"

namespace reg {
  
PetscErrorCode TransportKernelAdjointSL::ComputeBodyForcePart1() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
    ScalarType lambda  = pL[i];
    
    ScalarType lambdax = pLx[i];

    ScalarType rhs0 = lambdax*pDivVx[i];
    ScalarType rhs1 = (lambdax + ht*rhs0)*pDivV[i];

    // compute \lambda(x,t^{j+1})
    pLnext[i] = lambdax + 0.5*ht*(rhs0 + rhs1);

    // compute bodyforce
    pB[0][i] += scale*pGm[0][i]*lambda;
    pB[1][i] += scale*pGm[1][i]*lambda;
    pB[2][i] += scale*pGm[2][i]*lambda;
  }
} // omp

  PetscFunctionReturn(ierr);
}


PetscErrorCode TransportKernelAdjointSL::ComputeBodyForcePart2() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
    ScalarType lambda  = pL[i];
    // compute bodyforce
    pB[0][i] += scale*pGm[0][i]*lambda;
    pB[1][i] += scale*pGm[1][i]*lambda;
    pB[2][i] += scale*pGm[2][i]*lambda;
  }
} // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateSL::TimeIntegrationPart1() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
    pMtilde[i] -= hthalf*(pGm[0][i]*pVtildex[0][i]
                        + pGm[1][i]*pVtildex[1][i]
                        + pGm[2][i]*pVtildex[2][i]);
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateSL::TimeIntegrationPart2() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
    pMtilde[i] -= hthalf*(pGm[0][i]*pVtilde[0][i]
                        + pGm[1][i]*pVtilde[1][i]
                        + pGm[2][i]*pVtilde[2][i]);
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjoint::ComputeBodyForce() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
        // b = \sum_k\int_{\Omega} \lambda_k \grad m_k dt
        for (IntType i = 0; i < nl; ++i) {
            ScalarType lambda = pL[i];
            pB[0][i] += scale*lambda*pGm[0][i];
            pB[1][i] += scale*lambda*pGm[1][i];
            pB[2][i] += scale*lambda*pGm[2][i];
        }
} // omp

  PetscFunctionReturn(ierr);
}

template<typename T>
PetscErrorCode TransportKernelCopy(T* org, T* dest, IntType ne) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  try {
      std::copy(org,org+ne,dest);
  } catch (std::exception& err) {
      ierr = ThrowError(err); CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(ierr);
}
template PetscErrorCode TransportKernelCopy<ScalarType>(ScalarType*, ScalarType*, IntType);


} // namespace reg
