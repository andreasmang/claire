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

#ifndef _REGULARIZATIONH1SN_CPP_
#define _REGULARIZATIONH1SN_CPP_

#include "RegularizationH1SN.hpp"
#include "RegularizationKernel.hpp"

namespace reg {

/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegularizationH1SN::RegularizationH1SN() : SuperClass() {
}

/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegularizationH1SN::~RegularizationH1SN(void) {
    this->ClearMemory();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
RegularizationH1SN::RegularizationH1SN(RegOpt* opt) : SuperClass(opt) {
}

/********************************************************************
 * @brief evaluates the functional
 *******************************************************************/
PetscErrorCode RegularizationH1SN::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr = 0;
    ScalarType beta, value, hd;
    IntType rval;
    
    RegularizationKernel kernel;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization parameter
    beta = this->m_Opt->m_RegNorm.beta[0];
    hd  = this->m_Opt->GetLebesgueMeasure();   
    
    kernel.nl = this->m_Opt->m_Domain.nl;

    *R = 0.0;

    // if regularization weight is zero, do noting
    if (beta != 0.0) {
        ierr = Assert(this->m_WorkVecField != NULL, "null pointer"); CHKERRQ(ierr);
        
        /*ierr = this->m_Differentiation->RegLapOp(this->m_WorkVecField, v, 0.5*hd*beta); CHKERRQ(ierr);
        ierr = VecTDot(this->m_WorkVecField->m_X1, v->m_X1, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X2, v->m_X2, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X3, v->m_X3, &value); *R +=value;*/

        // X1 gradient
        ierr = this->m_Differentiation->Gradient(this->m_WorkVecField, v->m_X1); CHKERRQ(ierr);
        // compute inner products
        ierr = this->m_WorkVecField->GetArraysRead(kernel.pX);
        ierr = kernel.LocalNorm(value); CHKERRQ(ierr);
        ierr = this->m_WorkVecField->RestoreArrays();
        *R += value;
        
        //ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value); *R +=value;
        //ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value); *R +=value;
        //ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value); *R +=value;
        
        // X2 gradient
        ierr = this->m_Differentiation->Gradient(this->m_WorkVecField, v->m_X2); CHKERRQ(ierr);
        // compute inner products
        /*ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value); *R +=value;*/
        ierr = this->m_WorkVecField->GetArraysRead(kernel.pX);
        ierr = kernel.LocalNorm(value); CHKERRQ(ierr);
        ierr = this->m_WorkVecField->RestoreArrays();
        *R += value;

        // X3 gradient
        ierr = this->m_Differentiation->Gradient(this->m_WorkVecField, v->m_X3); CHKERRQ(ierr);
        // compute inner products
        /*ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value); *R +=value;*/
        ierr = this->m_WorkVecField->GetArraysRead(kernel.pX);
        ierr = kernel.LocalNorm(value); CHKERRQ(ierr);
        ierr = this->m_WorkVecField->RestoreArrays();
        *R += value;
        
        rval = MPI_Allreduce(MPI_IN_PLACE, R, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
        
        // multiply with regularization weight
        *R = 0.5*hd*beta*(*R);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief evaluates first variation of regularization norm
 *******************************************************************/
PetscErrorCode RegularizationH1SN::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr = 0;
    ScalarType beta, hd;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->m_RegNorm.beta[0];
    hd  = this->m_Opt->GetLebesgueMeasure();

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = dvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        ierr = this->m_Differentiation->RegLapOp(dvR, v, hd*beta); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief applies second variation of regularization norm to
 * a vector
 *******************************************************************/
PetscErrorCode RegularizationH1SN::HessianMatVec(VecField* dvvR, VecField* vtilde) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = this->EvaluateGradient(dvvR, vtilde); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
PetscErrorCode RegularizationH1SN::ApplyInverse(VecField* Ainvx, VecField* x, bool applysqrt) {
    PetscErrorCode ierr;
    ScalarType beta;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(Ainvx != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = Ainvx->Copy(x); CHKERRQ(ierr);
    } else {
        ierr = this->m_Differentiation->InvRegLapOp(Ainvx, x, applysqrt, beta); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
PetscErrorCode RegularizationH1SN::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr = 0;
    ScalarType w[3], beta1, beta2, regop;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    beta1=this->m_Opt->m_RegNorm.beta[0];
    beta2=this->m_Opt->m_RegNorm.beta[1];

    // get max value
    w[0] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[0])/2.0;
    w[1] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[1])/2.0;
    w[2] = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[2])/2.0;

    // compute largest value for operator
    regop = -(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]); // laplacian
    regop = -beta1*regop + beta2; // beta_1*laplacian + beta_2
    emin = 1.0/regop;
    emax = 1.0/beta2; // 1/(\beta_1*0 + \beta_2)

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

} // end of name space

#endif //_REGULARIZATIONH1SN_CPP_
