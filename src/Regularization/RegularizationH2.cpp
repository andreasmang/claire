/*************************************************************************
 *  Copyright (c) 2016.
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

#ifndef _REGULARIZATIONH2_CPP_
#define _REGULARIZATIONH2_CPP_

#include "RegularizationH2.hpp"

namespace reg {

/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegularizationH2::RegularizationH2() : SuperClass() {
}

/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegularizationH2::~RegularizationH2(void) {
    this->ClearMemory();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
RegularizationH2::RegularizationH2(RegOpt* opt) : SuperClass(opt) {
}

/********************************************************************
 * @brief evaluates the functional
 *******************************************************************/
PetscErrorCode RegularizationH2::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr;
    ScalarType sqrtbeta0, beta1, ipxi, hd;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // get regularization weight
    sqrtbeta0 = sqrt(this->m_Opt->m_RegNorm.beta[0]);
    beta1 = this->m_Opt->m_RegNorm.beta[1];
    hd  = this->m_Opt->GetLebesgueMeasure();   

    *R = 0.0;

    // if regularization weight is zero, do noting
    //if (sqrtbeta[0] != 0.0 && sqrtbeta[1] != 0.0) {
    if (sqrtbeta0 != 0.0) {
        ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_WorkVecField != NULL, "null pointer"); CHKERRQ(ierr);
        
        ierr = this->m_Differentiation->RegLapOp(this->m_WorkVecField, v, sqrtbeta0, beta1); CHKERRQ(ierr);
        
        ierr=VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &ipxi); CHKERRQ(ierr); *R += ipxi;
        ierr=VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &ipxi); CHKERRQ(ierr); *R += ipxi;
        ierr=VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &ipxi); CHKERRQ(ierr); *R += ipxi;

        // multiply with regularization weight
        *R *= 0.5*hd;
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief evaluates first variation of regularization norm
 *******************************************************************/
PetscErrorCode RegularizationH2::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr;
    ScalarType beta[2], hd;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization weight
    beta[0] = this->m_Opt->m_RegNorm.beta[0];
    beta[1] = this->m_Opt->m_RegNorm.beta[1];
    hd  = this->m_Opt->GetLebesgueMeasure();

    // if regularization weight is zero, do noting
    //if ((beta[0] == 0.0) && (beta[1] == 0.0)) {
    if (beta[0] == 0.0) {
        ierr = dvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        ierr = this->m_Differentiation->RegBiLapOp(dvR, v, hd*beta[0], beta[1]); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief applies second variation of regularization norm to
 * a vector
 *******************************************************************/
PetscErrorCode RegularizationH2::HessianMatVec(VecField* dvvR, VecField* vtilde) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = this->EvaluateGradient(dvvR, vtilde); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
PetscErrorCode RegularizationH2::ApplyInverse(VecField* Ainvx, VecField* x, bool applysqrt) {
    PetscErrorCode ierr;
    ScalarType beta[2];
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(Ainvx != NULL, "null pointer"); CHKERRQ(ierr);

    beta[0] = this->m_Opt->m_RegNorm.beta[0];
    beta[1] = this->m_Opt->m_RegNorm.beta[1];

    // if regularization weight is zero, do noting
    //if (beta[0] == 0.0 && beta[1] == 0.0) {
    if (beta[0] == 0.0) {
        ierr = VecCopy(x->m_X1, Ainvx->m_X1); CHKERRQ(ierr);
        ierr = VecCopy(x->m_X2, Ainvx->m_X2); CHKERRQ(ierr);
        ierr = VecCopy(x->m_X3, Ainvx->m_X3); CHKERRQ(ierr);
    } else {
        ierr = this->m_Differentiation->InvRegBiLapOp(Ainvx, x, applysqrt, beta[0], beta[1]); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
PetscErrorCode RegularizationH2::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr = 0;
    ScalarType w[3],beta1,beta2,regop;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    beta1 = this->m_Opt->m_RegNorm.beta[0];
    beta2 = this->m_Opt->m_RegNorm.beta[1];

    // get max value
    w[0] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[0])/2.0;
    w[1] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[1])/2.0;
    w[2] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[2])/2.0;

    // compute largest value for operator
    regop = -(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]); // laplacian
    regop = beta1*(regop*regop) + beta2; // beta * biharmonic
    emin = 1.0/regop;
    emax = 1.0/beta2; // 1/(0*beta_1 + beta_2)

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

} // end of name space

#endif //_REGULARIZATIONH2_CPP_
