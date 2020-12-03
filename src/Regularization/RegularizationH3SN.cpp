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

#ifndef _REGULARIZATIONREGISTRATIONH3SN_CPP_
#define _REGULARIZATIONREGISTRATIONH3SN_CPP_

#include "RegularizationH3SN.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegularizationH3SN::RegularizationH3SN() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegularizationH3SN::~RegularizationH3SN(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
RegularizationH3SN::RegularizationH3SN(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief evaluates the functional
 *******************************************************************/
PetscErrorCode RegularizationH3SN::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr;
    ScalarType beta, hd, value[9];

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // get regularization weight
    beta = this->m_Opt->m_RegNorm.beta[0];
    hd  = this->m_Opt->GetLebesgueMeasure();  

    *R = 0.0;

    // if regularization weight is zero, do noting
    if (beta != 0.0) {
        ierr = Assert(v != NULL,"null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_WorkVecField != NULL, "null pointer"); CHKERRQ(ierr);
        
        ierr = this->m_Differentiation->Laplacian(this->m_WorkVecField, v); CHKERRQ(ierr);
        
        // X1 gradient
        ierr = this->m_Differentiation->Gradient(this->m_WorkVecField, this->m_WorkVecField->m_X1); CHKERRQ(ierr);
        // compute inner products
        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value[0]); CHKERRQ(ierr);
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value[1]); CHKERRQ(ierr);
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value[2]); CHKERRQ(ierr);
        
        ierr = this->m_Differentiation->Laplacian(this->m_WorkVecField, v); CHKERRQ(ierr);
        
        // X2 gradient
        ierr = this->m_Differentiation->Gradient(this->m_WorkVecField, this->m_WorkVecField->m_X2); CHKERRQ(ierr);
        // compute inner products
        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value[3]); CHKERRQ(ierr);
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value[4]); CHKERRQ(ierr);
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value[5]); CHKERRQ(ierr);

        ierr = this->m_Differentiation->Laplacian(this->m_WorkVecField, v); CHKERRQ(ierr);
        
        // X3 gradient
        ierr = this->m_Differentiation->Gradient(this->m_WorkVecField, m_WorkVecField->m_X3); CHKERRQ(ierr);
        // compute inner products
        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value[6]); CHKERRQ(ierr);
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value[7]); CHKERRQ(ierr);
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value[8]); CHKERRQ(ierr);
        
        // multiply with regularization weight
        *R = 0.5*hd*beta*(value[0] + value[1] + value[2] 
                        + value[3] + value[4] + value[5]
                        + value[6] + value[7] + value[8]);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief evaluates first variation of regularization norm
 *******************************************************************/
PetscErrorCode RegularizationH3SN::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr;
    ScalarType beta, hd;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization weight
    beta = this->m_Opt->m_RegNorm.beta[0];
    hd = this->m_Opt->GetLebesgueMeasure();

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = dvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        ierr = this->m_Differentiation->RegTriLapOp(dvR, v, beta*hd); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies second variation of regularization norm to
 * a vector
 *******************************************************************/
PetscErrorCode RegularizationH3SN::HessianMatVec(VecField* dvvR, VecField* vtilde) {
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
PetscErrorCode RegularizationH3SN::ApplyInverse(VecField* Ainvx, VecField* x, bool applysqrt) {
    PetscErrorCode ierr;
    ScalarType beta;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr=Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(Ainvx != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr=VecCopy(x->m_X1, Ainvx->m_X1); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X2, Ainvx->m_X2); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X3, Ainvx->m_X3); CHKERRQ(ierr);
    } else {
        ierr = this->m_Differentiation->InvRegTriLapOp(Ainvx, x, applysqrt, beta); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
PetscErrorCode RegularizationH3SN::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr = 0;
    ScalarType w[3], beta, trihik, regop;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // get max value
    w[0] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[0])/2.0;
    w[1] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[1])/2.0;
    w[2] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[2])/2.0;

    trihik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
    trihik = std::pow(trihik, 3);

    // compute regularization operator
    regop = -beta*trihik;
    emin = 1.0/regop;
    emax = 1.0; // 1/(0*beta_1 + beta_2)

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // end of name space

#endif  // _REGULARIZATIONREGISTRATIONH3SN_CPP_
