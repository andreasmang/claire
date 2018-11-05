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

PetscErrorCode TransportKernelStateRK2::TimeIntegrationPart1() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
       pRHS[i] = - pGmx[0][i]*pVx[0][i] - pGmx[1][i]*pVx[1][i] - pGmx[2][i]*pVx[2][i];

       // compute intermediate result
       pMbar[i] = pM[i] + ht*pRHS[i];
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelStateRK2::TimeIntegrationPart2() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
      ScalarType rhs1 = - pGmx[0][i]*pVx[0][i] - pGmx[1][i]*pVx[1][i] - pGmx[2][i]*pVx[2][i];

      // m_{j+1} = m_j + 0.5*ht*(RHS0 + RHS1)
      pMnext[i] = pM[i] + 0.5*ht*(pRHS[i] + rhs1);
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointRK2::TimeIntegrationPart1() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i) {  // for all grid points
        ScalarType lambda = pL[i];
        pVec[0][i] = lambda*pV[0][i];
        pVec[1][i] = lambda*pV[1][i];
        pVec[2][i] = lambda*pV[2][i];
    }  // for all grid points
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointRK2::TimeIntegrationPart2() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {  // for all grid points
      ScalarType lambdabar = pL[i] + ht*pRHS[0][i];

      // scale \vect{v} by \bar{\lambda}
      pVec[0][i] = pV[0][i]*lambdabar;
      pVec[1][i] = pV[1][i]*lambdabar;
      pVec[2][i] = pV[2][i]*lambdabar;
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointRK2::TimeIntegrationPart3() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {  // for all grid points
      ScalarType lambda = pL[i];
      // second step of rk2 time integration
      pLnext[i] = lambda + 0.5*ht*(pRHS[0][i] + pRHS[1][i]);

      // compute bodyforce
      pB[0][i] += scale*pVec[0][i]*lambda;
      pB[1][i] += scale*pVec[1][i]*lambda;
      pB[2][i] += scale*pVec[2][i]*lambda;
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelAdjointRK2::TimeIntegrationPart4() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {  // for all grid points
      ScalarType lambda = pL[i];
      // compute bodyforce
      pB[0][i] += 0.5*scale*pVec[0][i]*lambda;
      pB[1][i] += 0.5*scale*pVec[1][i]*lambda;
      pB[2][i] += 0.5*scale*pVec[2][i]*lambda;
  }
} // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart1a() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {  // for all grid points
      ScalarType lambda = pL[i];

      // scale \vect{v} by \lambda
      pLtjVx[0][i] = pVtx[0][i]*lambda;
      pLtjVx[1][i] = pVtx[1][i]*lambda;
      pLtjVx[2][i] = pVtx[2][i]*lambda;
  }  // for all grid points
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart1b() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {  // for all grid points
      ScalarType lambda  = pL[i];
      ScalarType lambdatilde = pLt[i];

      pLtjVx[0][i] = pVx[0][i]*lambdatilde + pVtx[0][i]*lambda;
      pLtjVx[1][i] = pVx[1][i]*lambdatilde + pVtx[1][i]*lambda;
      pLtjVx[2][i] = pVx[2][i]*lambdatilde + pVtx[2][i]*lambda;
  }  // for all grid points
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart2a() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {  // for all grid points
      pLtnext[i] = pLt[i] + ht*pRHS[0][i];
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart2b() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {  // for all grid points
      // \bar{\lambda} = \tilde{\lambda}^j + ht*\idiv(\lambda^j\vect{v})
      ScalarType ltbar = pLt[i] + ht*pRHS[0][i];
      ScalarType lambda = pL[i];

      // v \bar{\lambda} + \vect{\tilde{v}}\lambda^{j+1}
      pLtjVx[0][i] = pVx[0][i]*ltbar + pVtx[0][i]*lambda;
      pLtjVx[1][i] = pVx[1][i]*ltbar + pVtx[1][i]*lambda;
      pLtjVx[2][i] = pVx[2][i]*ltbar + pVtx[2][i]*lambda;
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncAdjointRK2::TimeIntegrationPart3b() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {  // for all grid points
      pLtnext[i] = pLt[i] + 0.5*ht*(pRHS[0][i]+pRHS[1][i]);
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateRK2::TimeIntegrationEuler() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  // the right hand side remains constant;
  // we can reduce the 2 RK2 steps to a single one
  for (IntType i = 0; i < nl; ++i) {
       pMtnext[i] = pMt[i] - ht*(pGmx[0][i]*pVtx[0][i]
                               + pGmx[1][i]*pVtx[1][i]
                               + pGmx[3][i]*pVtx[2][i]);
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateRK2::TimeIntegrationPart1() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
       pRHS[i] = -pGmtx[0][i]*pVx[0][i] - pGmtx[1][i]*pVx[1][i] - pGmtx[2][i]*pVx[2][i]
                 -pGmx[0][i]*pVtx[0][i] - pGmx[1][i]*pVtx[1][i] - pGmx[2][i]*pVtx[2][i];
       // compute intermediate result
       pMtbar[i] = pMt[i] + ht*pRHS[i];
  }
} // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelIncStateRK2::TimeIntegrationPart2() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
      // evaluate right hand side
      ScalarType rhs1 = -pGmtx[0][i]*pVx[0][i] - pGmtx[1][i]*pVx[1][i] - pGmtx[2][i]*pVx[2][i]
                 -pGmx[0][i]*pVtx[0][i] - pGmx[1][i]*pVtx[1][i] - pGmx[2][i]*pVtx[2][i];

      // compute intermediate result
      pMtnext[i] = pMt[i] + 0.5*ht*(rhs1 + pRHS[i]);
  }
}  // omp

  PetscFunctionReturn(ierr);
}

PetscErrorCode TransportKernelContinuity::TimeIntegration() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
      ScalarType mx = pMx[i];

      ScalarType rhs0 = -mx*pDivVx[i];
      ScalarType rhs1 = -(mx + ht*rhs0)*pDivV[i];
      //if (std::abs(p_divv[i]) > 0.1) { std::cout << p_divv[i] << " ";}
      // compute \lambda(x,t^{j+1})
      pMnext[i] = mx + 0.5*ht*(rhs0 + rhs1);
  }
}

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
