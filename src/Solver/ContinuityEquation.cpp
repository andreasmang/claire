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

#ifndef _CONTINUITYEQUATION_CPP_
#define _CONTINUITYEQUATION_CPP_

#include "ContinuityEquation.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
ContinuityEquation::ContinuityEquation() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
ContinuityEquation::~ContinuityEquation() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
ContinuityEquation::ContinuityEquation(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode ContinuityEquation::Initialize() {
    PetscFunctionBegin;

    this->m_SemiLagrangianMethod = nullptr;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode ContinuityEquation::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    Free(this->m_SemiLagrangianMethod);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem (i.e., the continuity equation)
 *******************************************************************/
PetscErrorCode ContinuityEquation::SolveForwardProblem() {
    PetscErrorCode ierr = 0;
    TransportKernelContinuity kernel;
    ScalarType *pM = nullptr;
    ScalarType nt, nc, l, lnext;
    bool store = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);

    // flag to identify if we store the time history
    store = this->m_Opt->m_RegFlags.runinversion;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    kernel.ht = this->m_Opt->GetTimeStepSize();

    ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);

    ierr = this->m_WorkScaField[0]->GetArrayWrite(kernel.pDivV); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->GetArrayWrite(kernel.pDivVx); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[2]->GetArrayWrite(kernel.pMx); CHKERRQ(ierr);

    // compute divergence of velocity field
    ierr = this->m_Differentiation->Divergence(kernel.pDivV, this->m_VelocityField); CHKERRQ(ierr);

    // evaluate div(v) along characteristic X
    ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pDivVx, kernel.pDivV, "state"); CHKERRQ(ierr);

    // perform numerical time integration for state variable
    for (IntType j = 0; j < nt; ++j) {
        if (store) {
            l = j; lnext = (j+1);
        } else {
            l = 0; lnext = 0;
        }
        // scaling for trapezoidal rule (for body force)
        for (IntType k = 0; k < nc; ++k) {
            ierr = this->m_StateVariable->GetArrayReadWrite(pM, k, l); CHKERRQ(ierr);
            ierr = this->m_StateVariable->GetArrayReadWrite(kernel.pMnext, k, lnext); CHKERRQ(ierr);
            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pMx, pM, "state"); CHKERRQ(ierr);

            ierr = kernel.TimeIntegration(); CHKERRQ(ierr);
        }
    }

    ierr = this->m_WorkScaField[0]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[1]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_WorkScaField[2]->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);
    
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode ContinuityEquation::SolveAdjointProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode ContinuityEquation::SolveIncForwardProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode ContinuityEquation::SolveIncAdjointProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _CONTINUITYEQUATION_CPP_
