/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGULARIZATIONREGISTRATIONL2_CPP_
#define _REGULARIZATIONREGISTRATIONL2_CPP_

#include "RegularizationRegistrationL2.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegularizationRegistrationL2::RegularizationRegistrationL2() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegularizationRegistrationL2::~RegularizationRegistrationL2(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
RegularizationRegistrationL2::RegularizationRegistrationL2(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief evaluates the functional
 *******************************************************************/
PetscErrorCode RegularizationRegistrationL2
::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr = 0;
    ScalarType beta, ipxi;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization weight
    beta = this->m_Opt->m_RegNorm.beta[0];

    *R = 0.0;

    // if regularization weight is zero, do noting
    if (beta != 0.0) {
        ierr = VecTDot(v->m_X1, v->m_X1, &ipxi); CHKERRQ(ierr); *R += ipxi;
        ierr = VecTDot(v->m_X2, v->m_X2, &ipxi); CHKERRQ(ierr); *R += ipxi;
        ierr = VecTDot(v->m_X3, v->m_X3, &ipxi); CHKERRQ(ierr); *R += ipxi;
        *R *= 0.5*beta;
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates first variation of regularization norm
 *******************************************************************/
PetscErrorCode RegularizationRegistrationL2
::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr = 0;
    ScalarType beta, *p_dvR1 = NULL, *p_dvR2 = NULL, *p_dvR3 = NULL,
                *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL;
    IntType nl;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization weight
    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = dvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        nl = this->m_Opt->m_Domain.nl;

        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        ierr = dvR->GetArrays(p_dvR1, p_dvR2, p_dvR3); CHKERRQ(ierr);
        for (IntType i = 0; i < nl; ++i) {
            p_dvR1[i] = beta*p_v1[i];
            p_dvR2[i] = beta*p_v2[i];
            p_dvR3[i] = beta*p_v3[i];
        }
        ierr = dvR->RestoreArrays(p_dvR1, p_dvR2, p_dvR3); CHKERRQ(ierr);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies second variation of regularization norm to
 * a vector
 *******************************************************************/
PetscErrorCode RegularizationRegistrationL2
::HessianMatVec(VecField* dvvR, VecField* vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType beta;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(dvvR != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vtilde != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = dvvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        ierr = this->EvaluateGradient(dvvR, vtilde); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
PetscErrorCode RegularizationRegistrationL2
::ApplyInvOp(VecField* Ainvx, VecField* x, bool applysqrt) {
    PetscErrorCode ierr = 0;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL,
                *p_Ainvx1 = NULL, *p_Ainvx2 = NULL, *p_Ainvx3 = NULL, beta;
    IntType nl;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(Ainvx != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = VecCopy(x->m_X1, Ainvx->m_X1); CHKERRQ(ierr);
        ierr = VecCopy(x->m_X2, Ainvx->m_X2); CHKERRQ(ierr);
        ierr = VecCopy(x->m_X3, Ainvx->m_X3); CHKERRQ(ierr);
    } else {
        nl = this->m_Opt->m_Domain.nl;
        ierr = x->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
        ierr = Ainvx->GetArrays(p_Ainvx1, p_Ainvx2, p_Ainvx3); CHKERRQ(ierr);
        for (IntType i = 0; i < nl; ++i) {
            p_Ainvx1[i] = p_x1[i]/beta;
            p_Ainvx2[i] = p_x2[i]/beta;
            p_Ainvx3[i] = p_x3[i]/beta;
        }
        ierr = x->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
        ierr = Ainvx->RestoreArrays(p_Ainvx1, p_Ainvx2, p_Ainvx3); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
PetscErrorCode RegularizationRegistrationL2
::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr = 0;
    ScalarType beta;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // compute largest value for operator
    emin = 1.0/beta;
    emax = 1.0/beta;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _REGULARIZATIONREGISTRATIONL2_CPP_
