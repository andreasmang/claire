/*************************************************************************
 *  Copyright (c) 2017.
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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _DISTANCEMEASURE_CPP_
#define _DISTANCEMEASURE_CPP_

#include "DistanceMeasure.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DistanceMeasure::DistanceMeasure() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DistanceMeasure::~DistanceMeasure() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DistanceMeasure::DistanceMeasure(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DistanceMeasure::Initialize() {
    PetscFunctionBegin;

    this->m_Opt = nullptr;

    this->m_Mask = nullptr;
    this->m_AuxVar1 = nullptr;
    this->m_AuxVar2 = nullptr;
    this->m_TemplateImage = nullptr;
    this->m_ReferenceImage = nullptr;
    this->m_StateVariable = nullptr;
    this->m_AdjointVariable = nullptr;
    this->m_IncStateVariable = nullptr;
    this->m_IncAdjointVariable = nullptr;
    this->m_WorkVecField1 = nullptr;
    this->m_WorkVecField2 = nullptr;
    this->m_WorkVecField3 = nullptr;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DistanceMeasure::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief set up scale
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetupScale() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetReferenceImage(ScaField* mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mR != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_ReferenceImage = mR;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetTemplateImage(ScaField* mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mT != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_TemplateImage = mT;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief set temporary vector fields 
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetWorkVecField(VecField* v, int id) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != nullptr, "null pointer"); CHKERRQ(ierr);
    switch (id) {
        case 1:
           this->m_WorkVecField1 = v;
           break;
        case 2:
           this->m_WorkVecField2 = v;
           break;
        case 3:
           this->m_WorkVecField3 = v;
           break;
	default:
           ierr = ThrowError("id not defined"); CHKERRQ(ierr);
           break;
  }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief set mask
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetMask(ScaField* mask) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mask != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_Mask = mask;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set state variable
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetStateVariable(ScaField* m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_StateVariable = m;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set incremental state variable
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetIncStateVariable(ScaField* mtilde) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mtilde != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_IncStateVariable = mtilde;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief set adjoint variable
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetAdjointVariable(ScaField* lambda) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_AdjointVariable = lambda;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set adjoint variable
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetIncAdjointVariable(ScaField* lambdatilde) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambdatilde != nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_IncAdjointVariable = lambdatilde;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set auxilary variable
 *******************************************************************/
PetscErrorCode DistanceMeasure::SetAuxVariable(ScaField* x, int id) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(x != nullptr, "null pointer"); CHKERRQ(ierr);
    if (id == 1) {
        this->m_AuxVar1 = x;
    } else if (id == 2) {
        this->m_AuxVar2 = x;
    } else {
        ierr = ThrowError("id not defined"); CHKERRQ(ierr);
    }
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _DISTANCEMEASURE_CPP_
