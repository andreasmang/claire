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

#include "TransportProblem.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
TransportProblem::TransportProblem() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
TransportProblem::~TransportProblem() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
TransportProblem::TransportProblem(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode TransportProblem::Initialize() {
    PetscFunctionBegin;

    this->m_Opt = nullptr;
    this->m_TemplateImage = nullptr;
    this->m_ReferenceImage = nullptr;
    this->m_StateVariable = nullptr;
    this->m_AdjointVariable = nullptr;
    this->m_IncStateVariable = nullptr;
    this->m_IncAdjointVariable = nullptr;
    
    this->m_VelocityField = nullptr;
    this->m_IncVelocityField = nullptr;
    
    for (int i = 0;i < 5;++i) {
      this->m_WorkScaField[i] = nullptr;
      this->m_WorkVecField[i] = nullptr;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode TransportProblem::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode TransportProblem::SetReferenceImage(Vec mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mR != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_ReferenceImage = mR;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set template image (i.e., the image to be deformed)
 *******************************************************************/
PetscErrorCode TransportProblem::SetTemplateImage(Vec mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mT != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_TemplateImage = mT;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set state variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetStateVariable(Vec m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_StateVariable = m;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set adjoint variable
 *******************************************************************/
PetscErrorCode TransportProblem::SetAdjointVariable(Vec lambda) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_AdjointVariable = lambda;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief set velocity field
 *******************************************************************/
PetscErrorCode TransportProblem::SetControlVariable(VecField *field) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_VelocityField = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set incremental velocity field
 *******************************************************************/
PetscErrorCode TransportProblem::SetIncControlVariable(VecField *field) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_IncVelocityField = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set working scalar field
 ********************************************************************/
PetscErrorCode TransportProblem::SetWorkScaField(Vec field, IntType idx) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(idx >= 0 && idx < 5, "index out of range"); CHKERRQ(ierr);
    
    this->m_WorkScaField[idx] = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set working scalar field
 ********************************************************************/
PetscErrorCode TransportProblem::SetWorkVecField(VecField* field, IntType idx) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(field != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(idx >= 0 && idx < 5, "index out of range"); CHKERRQ(ierr);
    
    this->m_WorkVecField[idx] = field;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set the differentiation method
 *******************************************************************/
PetscErrorCode TransportProblem::SetDifferentiation(Differentiation *diff) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    this->m_Differentiation = diff;
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem for zero velocity fields with Gauss-Newton-scheme
 *******************************************************************/
PetscErrorCode TransportProblem::SolveIncAdjointProblem() {
    PetscErrorCode ierr = 0;
    ScalarType *pM = nullptr, *pLtilde = nullptr;
    IntType nc, nl;
    TransportKernelIncAdjointGN kernel;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    kernel.nl = nl;
    kernel.scale = 1.0/static_cast<ScalarType>(nc);
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);

    // copy terminal condition \tilde{\lambda}_1 = -\tilde{m}_1 to all time points
    ierr = GetRawPointer(this->m_IncAdjointVariable, &pLtilde); CHKERRQ(ierr);

    // m and \lambda are constant in time
    ierr = GetRawPointer(this->m_StateVariable, &pM); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArrays(kernel.pGm); CHKERRQ(ierr);

    // init body force for numerical integration
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArrays(kernel.pBtilde); CHKERRQ(ierr);

    // $m$ and $\tilde{\lambda}$ are constant
    for (IntType k = 0; k < nc; ++k) {  // for all components
        kernel.pLtilde = pLtilde + k*nl;
        // compute gradient of m
        ierr = this->m_Differentiation->Gradient(kernel.pGm,pM + k*nl); CHKERRQ(ierr);

        ierr = kernel.ComputeBodyForce(); CHKERRQ(ierr);
    }
    ierr = this->m_WorkVecField[1]->RestoreArrays(kernel.pBtilde); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->RestoreArrays(kernel.pGm); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &pM); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &kernel.pLtilde); CHKERRQ(ierr);
  
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


}  // namespace reg




#endif  // _CONTINUITYEQUATION_CPP_
