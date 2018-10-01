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

#ifndef _TRANSPORTEQUATION_CPP_
#define _TRANSPORTEQUATION_CPP_

#include "TransportEquationRK2.hpp"

namespace reg {

/********************************************************************
 * @brief default constructor
 *******************************************************************/
TransportEquationRK2::TransportEquationRK2() : SuperClass() {
    this->Initialize();
}

/********************************************************************
 * @brief default destructor
 *******************************************************************/
TransportEquationRK2::~TransportEquationRK2() {
    this->ClearMemory();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
TransportEquationRK2::TransportEquationRK2(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode TransportEquationRK2::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode TransportEquationRK2::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the forward problem (i.e., the continuity equation)
 *******************************************************************/
PetscErrorCode TransportEquationRK2::SolveForwardProblem() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
    if (VelocityIsZero) {
        ierr = SuperClass::SolveForwardProblem(); CHKERRQ(ierr);
    } else {
        ierr = this->SolveStateEquation(); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationRK2::SolveAdjointProblem() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
    if (VelocityIsZero) {
        ierr = SuperClass::SolveAdjointProblem(); CHKERRQ(ierr);
    } else {
        ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationRK2::SolveIncForwardProblem() {
    PetscErrorCode ierr = 0;
    IntType nt, nc, lmt, lmtnext;
    bool fullnewton = false;
    bool VelocityIsZero = false;
    const ScalarType *pM = nullptr;
    TransportKernelIncStateRK2 kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {   // gauss newton
        fullnewton = true;
    }

    ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pGmx); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArraysRead(kernel.pVtx); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->m_VelocityField->IsVelocityZero(VelocistyIsZero); CHKERRQ(ierr);
    if (VelocityIsZero) {
        // compute gradient of first time point of image component
        for (IntType k = 0; k < nc; ++k) {
            ierr = this->m_StateVariable->GetArrayRead(pM, k); CHKERRQ(ierr);
            
            // template image is constant in time
            ierr = this->m_Differentiation->Gradient(kernel.pGmx, pM); CHKERRQ(ierr);

            // compute incremental state variable for all time points
            // note: we do not need to store the time history for
            // \tilde{m} if we consider a Gauss--Newton approximation
            for (IntType j = 0; j < nt; ++j) {
                if (fullnewton) {
                    lmt = j; lmtnext = j+1;
                } else {
                    lmt = 0; lmtnext = 0;
                }
                
                ierr = this->m_IncStateVariable->GetArrayReadWrite(kernel.pMt, k, lmt); CHKERRQ(ierr);
                ierr = this->m_IncStateVariable->GetArrayReadWrite(kernel.pMtnext, k, lmtnext); CHKERRQ(ierr);
                
                // the right hand side remains constant;
                // we can reduce the 2 RK2 steps to a single one
                ierr = kernel.TimeIntegrationEuler(); CHKERRQ(ierr);
            }  // for all time points
        }  // for all image components
    } else {  // velocity field is non-zero
        ierr = this->m_WorkScaField[0]->GetArrayWrite(kernel.pMtbar); CHKERRQ(ierr);
        ierr = this->m_WorkScaField[1]->GetArrayWrite(kernel.pRHS); CHKERRQ(ierr);

        ierr = this->m_WorkVecField[1]->GetArraysWrite(kernel.pGmtx); CHKERRQ(ierr);
        ierr = this->m_VelocityField->GetArraysWrite(kernel.pVx); CHKERRQ(ierr);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j) {
            for (IntType k = 0; k < nc; ++k) {
                lm = j*nl*nc; lmnext = (j+1)*nl*nc;
                if (fullnewton) {
                    lmt = j*nl*nc; lmtnext = (j+1)*nl*nc;
                } else {
                    lmt = 0; lmtnext = 0;
                }
                ierr = this->m_StateVariable->GetArrayRead(pM, k, j); CHKERRQ(ierr);
                ierr = this->m_IncStateVariable->GetArrayReadWrite(kernel.pMt, k, lmt); CHKERRQ(ierr);
                ierr = this->m_IncStateVariable->GetArrayReadWrite(kernel.pMtnext, k, lmtnext); CHKERRQ(ierr);

                // compute gradient of m_j
                ierr = this->m_Differentiation->Gradient(kernel.pGmx, pM); CHKERRQ(ierr);

                // compute gradient of \tilde{m}_j
                ierr = this->m_Differentiation->Gradient(kernel.pGmtx, kernel.pMt); CHKERRQ(ierr);

                // compute intermediate result
                ierr = kernel.TimeIntegrationPart1(); CHKERRQ(ierr);
                
                ierr = this->m_StateVariable->GetArrayRead(pM, k, j+1); CHKERRQ(ierr);

                // compute gradient of m_{j+1}
                ierr = this->m_Differentiation->Gradient(kernel.pGmx, pM); CHKERRQ(ierr);

                // compute gradient of \tilde{m}_j
                ierr = this->m_Differentiation->Gradient(kernel.pGmtx, kernel.pMtbar); CHKERRQ(ierr);

                // compute final RK2 step
                ierr = kernel.TimeIntegrationPart2(); CHKERRQ(ierr);
            }  // for all image components
        }  // for all time points

        // copy initial condition to buffer
        ierr = this->m_WorkScaField[0]->RestoreArray(); CHKERRQ(ierr);
        ierr = this->m_WorkScaField[1]->RestoreArray(); CHKERRQ(ierr);

        ierr = this->m_WorkVecField[1]->RestoreArrays(); CHKERRQ(ierr);
        ierr = this->m_VelocityField->RestoreArrays(); CHKERRQ(ierr);
    }  // velzero

    ierr = this->m_IncVelocityField->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(); CHKERRQ(ierr);

    ierr = this->m_IncStateVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationRK2::SolveIncAdjointProblem() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    
    if (this->m_Opt->m_OptPara.method == GAUSSNEWTON) {   // gauss newton
        ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
        if (VelocityIsZero) {
            ierr = SuperClass::SolveIncAdjointProblem(); CHKERRQ(ierr);
        } else {
            ierr = this->SolveIncAdjointEquationGN(); CHKERRQ(ierr);
        }
    } else if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ierr = ThrowError("not tested"); CHKERRQ(ierr);
        ierr = this->SolveIncAdjointEquationFN(); CHKERRQ(ierr);
    } else {
        ierr = ThrowError("update method not defined"); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the incremental adjoint problem with Gauss Newton
 *******************************************************************/
PetscErrorCode TransportEquationRK2::SolveIncAdjointEquationGN() {
    PetscErrorCode ierr = 0;
    IntType = nt, nc;
    const ScalarType *pM = nullptr;
    TransportKernelAdjointRK2 kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    kernel.scale = ht/static_cast<ScalarType>(nc);
    
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariableVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != NULL, "null pointer"); CHKERRQ(ierr);
    
    ierr = this->m_WorkScaField[0]->GetArrayWrite(kernel.pRHS[0]); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->GetArrayWrite(kernel.pRHS[1]); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pVec); CHKERRQ(ierr);
    ierr = this->m_VelocityField->GetArraysRead(kernel.pV); CHKERRQ(ierr);

    // init body force for numerical integration
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArraysReadWrite(kernel.pB); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        if (j == 0) kernel.scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pL, k); CHKERRQ(ierr);
            ierr = this->m_StateVariable->GetArrayRead(pM, k, nt-j); CHKERRQ(ierr);
            kernel.pLnext = kernel.pL;
            
            // scale \vect{v} by \lambda
            ierr = kernel.TimeIntegrationPart1(); CHKERRQ(ierr);

            // compute \idiv(\tilde{\lambda}\vect{v})
            ierr = this->m_Differentiation->Divergence(kernel.pRHS[0], kernel.pVec); CHKERRQ(ierr);

            // compute \bar{\tilde{\lambda}} = \tilde{\lambda}^j + ht*\idiv(\tilde{\lambda}^j\vect{v})
            // scale \vect{v} by \bar{\lambda}
            ierr = kernel.TimeIntegrationPart2(); CHKERRQ(ierr);

            // compute \idiv(\bar{\lambda}\vect{v})
            ierr = this->m_Differentiation->Divergence(kernel.pRHS[1], kernel.pVec); CHKERRQ(ierr);

            // compute gradient of m^j
            ierr = this->m_Differentiation->Gradient(kernel.pVec, pM); CHKERRQ(ierr);

            // compute integration
            // compute incremental body force
            ierr = kernel.TimeIntegrationPart3(); CHKERRQ(ierr);
        }  // for all image components
        if (j == 0) kernel.scale *= 2.0;
    }  // for all time points

    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ierr = this->m_StateVariable->GetArrayRead(pM, k); CHKERRQ(ierr);
        ierr = this->m_IncAdjointVariable(kernel.pL, k); CHKERRQ(ierr);

        // compute gradient of m (for incremental body force)
        ierr = this->m_Differentiation->Gradient(kernel.pVec, pM); CHKERRQ(ierr);

        // compute bodyforce
        ierr = kernel.TimeIntegrationPart4(); CHKERRQ(ierr);
    }

    ierr = this->m_VelocityField->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[0]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_IncAdjointVariable->RestoreArray(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the incremental adjoint problem in full Newton
 *******************************************************************/
PetscErrorCode TransportEquationRK2::SolveIncAdjointEquationFN() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    IntType nt, nc;
    TransportKernelIncAdjointRK2 kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = this->m_WorkScaField[0]->GetArrayWrite(kernel.p_RHS[0]); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pLtjVx); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArraysRead(kernel.pVtx); CHKERRQ(ierr);

    ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
    if (VelocityIsZero) {
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            // lambda and v are constant in time
            ierr = this->m_AdjointVariable->GetArrayRead(kernel.pL, k); CHKERRQ(ierr);
            
            // scale \vect{v} by \lambda
            ierr = kernel.TimeIntegrationPart1a(); CHKERRQ(ierr);

            // compute \idiv(\tilde{\lambda}\vect{v})
            ierr = this->m_Differentiation->Divergence(kernel.pRHS[0],kernel.pLtjVx); CHKERRQ(ierr);

            // compute numerical time integration
            for (IntType j = 0; j < nt; ++j) {  // for all time points
                ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pLt, k, nt-j); CHKERRQ(ierr);
                ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pLtnext, k, nt-(j+1)); CHKERRQ(ierr);
                
                ierr = kernel.TimeIntegrationPart2a(); CHKERRQ(ierr);
            }  // for all time points
        }  // for all image components
    } else {  // velocity is zero
        ierr = this->m_VelocityField->GetArraysRead(kernel.pVx); CHKERRQ(ierr);
        ierr = this->m_WorkScaField[1]->GetArrayWrite(kernel.pRHS[1]); CHKERRQ(ierr);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j) {  // for all time points
            for (IntType k = 0; k < nc; ++k) {  // for all image components
                ierr = this->m_AdjointVariable->GetArrayRead(kernel.pL, k, nt-j); CHKERRQ(ierr);
                ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pLt, k, nt-j); CHKERRQ(ierr);
                ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pLtnext, k, nt-(j+1)); CHKERRQ(ierr);

                ierr = kernel.TimeIntegrationPart1b(); CHKERRQ(ierr);

                // compute \idiv(\tilde{\lambda}\vect{v})
                ierr = this->m_Differentiation->Divergence(kernel.pRHS[0],kernel.pLtjVx); CHKERRQ(ierr);
                
                ierr = this->m_AdjointVariable->GetArrayRead(kernel.pL, k, nt-(j+1)); CHKERRQ(ierr);

                // \bar{\lambda} = \tilde{\lambda}^j + ht*\idiv(\lambda^j\vect{v})
                // v \bar{\lambda} + \vect{\tilde{v}}\lambda^{j+1}
                kernel.TimeIntegrationPart2b(); CHKERRQ(ierr);

                // compute \idiv(\bar{\lambda}\vect{v})
                ierr = this->m_Differentiation->Divergence(kernel.pRHS[1],kernel.pLtjVx); CHKERRQ(ierr);

                ierr = kernel.TimeIntegrationPart3b(); CHKERRQ(ierr);
            }  // for all image components
        }  // for all time points
        ierr = this->m_VelocityField->RestoreArrays(); CHKERRQ(ierr);
        ierr = this->m_WorkScaField[1]->RestoreArray(); CHKERRQ(ierr);
    }  // velzero

    ierr = this->m_AdjointVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_IncAdjointVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[0]->RestoreArray(); CHKERRQ(ierr);

    ierr = this->m_IncVelocityField->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->RestoreArrays(); CHKERRQ(ierr);

    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationRK2::SolveAdjointEquation() {
    PetscErrorCode ierr = 0;
    IntType nt, nc, ll, llnext;
    const ScalarType *pM = nullptr;
    TransportKernelAdjointRK2 kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    kernel.scale = kernel.ht/static_cast<ScalarType>(nc);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != NULL, "null pointer"); CHKERRQ(ierr);

    // for full newton we store $\lambda$
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        fullnewton = true;
    }
    
    ierr = this->m_WorkScaField[0]->GetArrayWrite(kernel.pRHS[0]); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->GetArrayWrite(kernel.pRHS[1]); CHKERRQ(ierr);

    ierr = this->m_VelocityField->GetArraysRead(kernel.pV); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pVec); CHKERRQ(ierr);

    // init body force for numerical integration
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArraysReadWrite(kernel.pB); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        if (fullnewton) {
            ll = nt-j; llnext = nt-(j+1);
        } else {
            ll = 0; llnext = 0;
        }

        // scaling for trapezoidal rule
        if (j == 0) scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL, k, ll); CHKERRQ(ierr);
            ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pLnext, k, llnext); CHKERRQ(ierr);
            ierr = this->m_StateVariable->GetArrayRead(pM, k, nt-j); CHKERRQ(ierr);
            // scale \vect{v} by \lambda
            ierr = kernel.TimeIntegrationPart1(); CHKERRQ(ierr);

            // compute \idiv(\lambda\vect{v})
            ierr = this->m_Differentiation->Divergence(kernel.pRHS[0], kernel.pVec); CHKERRQ(ierr);

            // compute \bar{\lambda} = \lambda_j + ht*\idiv(\lambda\vect{v})
            ierr = kernel.TimeIntegrationPart2();

            // compute \idiv(\bar{\lambda}\vect{v})
            ierr = this->m_Differentiation->Divergence(kernel.pRHS[1], kernel.pVec); CHKERRQ(ierr);

            // grad(m^j)
            ierr = this->m_Differentiation->Gradient(kernel.pVec, pM); CHKERRQ(ierr);

            // second step of rk2 time integration
            ierr = kernel.TimeIntegrationPart3(); CHKERRQ(ierr);
        }  // for all image components
        // trapezoidal rule (revert scaling)
        if (j == 0) scale *= 2.0;
    }  // for all time points

    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ierr = this->m_StateVariable->GetArrayRead(pM, k); CHKERRQ(ierr);
        ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL, k); CHKERRQ(ierr);

        // compute gradient of m (for incremental body force)
        ierr = this->m_Differentiation->Gradient(kernel.pVec, pM); CHKERRQ(ierr);

        ierr = kernel.TimeIntegrationPart4(); CHKERRQ(ierr);
    }

    ierr = this->m_AdjointVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(); CHKERRQ(ierr);

    ierr = this->m_WorkScaField[0]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->RestoreArray(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the state problem
 *******************************************************************/
PetscErrorCode TransportEquationRK2::SolveStateEquation() {
    PetscErrorCode ierr = 0;
    bool store;
    IntType nt, nc, l, lnext;
    TransportKernelStateRK2 kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    // flag to identify if we store the time history
    store = this->m_Opt->m_RegFlags.runinversion;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    kernel.nl = nl;
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = this->m_VelocityField->GetArraysRead(kernel.pVx); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArraysWrite(kernel.pGmx); CHKERRQ(ierr);

    // copy initial condition to buffer
    ierr = this->m_WorkScaField[0]->GetArrayWrite(kernel.pMbar); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->GetArrayWrite(kernel.pRHS); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {
        if (store) {
            l = j; lnext = (j+1);
        } else {
            l = 0; lnext = 0;
        }

        for (IntType k = 0; k < nc; ++k) {
            ierr = this->m_StateVariable->GetArrayReadWrite(kernel.pM, k, l);  CHKERRQ(ierr);
            ierr = this->m_StateVariable->GetArrayReadWrite(kernel.pMnext, k, lnext); CHKERRQ(ierr);
            
            // compute gradient of k-th component of m_j
            ierr = this->m_Differentiation->Gradient(kernel.pGmx, kernel.pM); CHKERRQ(ierr);

            // evaluate right hand side and compute intermediate rk2 step
            ierr = kernel.TimeIntegrationPart1(); CHKERRQ(ierr);

            // compute gradient of \bar{m}
            ierr = this->m_Differentiation->Gradient(kernel.pGmx, kernel.pMbar); CHKERRQ(ierr);

            // evaluate right hand side and wrap up integration
            ierr = kernel.TimeIntegrationPart2(); CHKERRQ(ierr);
        }  // for all components
    }  // for all time points

    // copy initial condition to buffer
    ierr = this->m_WorkScaField[0]->RestoreArray() CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

}  // namespace reg
#endif
