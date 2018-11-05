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

#ifndef _REGULARIZATIONH2SN_CPP_
#define _REGULARIZATIONH2SN_CPP_

#include "RegularizationH2SN.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegularizationH2SN::RegularizationH2SN() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegularizationH2SN::~RegularizationH2SN(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
RegularizationH2SN::RegularizationH2SN(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief evaluates the functional (we have to promote everything
 * to double to be able to solve the problem accurately; we loose
 * too many digits here)
 *******************************************************************/
PetscErrorCode RegularizationH2SN::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr = 0;
    ScalarType beta, ipxi, value, hd;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization weight
    beta = static_cast<ScalarType>(this->m_Opt->m_RegNorm.beta[0]);
    hd  = this->m_Opt->GetLebesgueMeasure();   

    *R = 0.0; value = 0.0;

    // if regularization weight is zero, do noting
    if (beta != 0.0) {
        ierr = Assert(this->m_WorkVecField != NULL, "null pointer"); CHKERRQ(ierr);
        
        ierr = this->m_Differentiation->Laplacian(this->m_WorkVecField, v); CHKERRQ(ierr);

        // compute inner product
        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &ipxi); CHKERRQ(ierr); value += ipxi;
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &ipxi); CHKERRQ(ierr); value += ipxi;
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &ipxi); CHKERRQ(ierr); value += ipxi;

        // multiply with regularization weight
        *R = static_cast<ScalarType>(0.5*beta*hd*value);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates first variation of regularization norm
 * @param[in] v velocity field
 * @param[out] dvR gradient of regularization (evaluated)
 *******************************************************************/
PetscErrorCode RegularizationH2SN::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr = 0;
    ScalarType beta, hd;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    beta = static_cast<double>(this->m_Opt->m_RegNorm.beta[0]);
    hd  = this->m_Opt->GetLebesgueMeasure();   

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = VecSet(dvR->m_X1, 0.0); CHKERRQ(ierr);
        ierr = VecSet(dvR->m_X2, 0.0); CHKERRQ(ierr);
        ierr = VecSet(dvR->m_X3, 0.0); CHKERRQ(ierr);
    } else {
        ierr = this->m_Differentiation->Bilaplacian(dvR, v, beta*hd); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies second variation of regularization norm to vector
 * @param dvvR regularization operator applied to vector \tilde{v}
 * @param vtilde incremental vector field \tilde{v}
 *******************************************************************/
PetscErrorCode RegularizationH2SN::HessianMatVec(VecField* dvvR, VecField* vtilde) {
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
PetscErrorCode RegularizationH2SN::ApplyInverse(VecField* ainvv, VecField* v, bool applysqrt) {
    PetscErrorCode ierr = 0;
    ScalarType beta;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(ainvv != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = VecCopy(v->m_X1, ainvv->m_X1); CHKERRQ(ierr);
        ierr = VecCopy(v->m_X2, ainvv->m_X2); CHKERRQ(ierr);
        ierr = VecCopy(v->m_X3, ainvv->m_X3); CHKERRQ(ierr);
    } else {
        if (applysqrt) {
          ierr = this->m_Differentiation->InverseBilaplacianSqrt(ainvv, v, beta); CHKERRQ(ierr);
        } else {
          ierr = this->m_Differentiation->InverseBilaplacian(ainvv, v, beta); CHKERRQ(ierr);
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
PetscErrorCode RegularizationH2SN::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr = 0;
    ScalarType w[3], beta, regop;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // get max value
    w[0] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[0])/2.0;
    w[1] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[1])/2.0;
    w[2] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[2])/2.0;

    // compute largest value for operator
    regop = -(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]); // laplacian
    regop = beta*(regop*regop); // beta * biharmonic
    emin = 1.0/regop;
    emax = 1.0; // by definition; it's 1/0

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _REGULARIZATIONH2SN_CPP_
