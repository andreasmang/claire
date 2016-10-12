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

#ifndef _REGULARIZATIONREGISTRATIONH1SN_CPP_
#define _REGULARIZATIONREGISTRATIONH1SN_CPP_

#include "RegularizationRegistrationH1SN.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH1SN"
RegularizationRegistrationH1SN::RegularizationRegistrationH1SN() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegularizationRegistrationH1SN"
RegularizationRegistrationH1SN::~RegularizationRegistrationH1SN(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH1SN"
RegularizationRegistrationH1SN::RegularizationRegistrationH1SN(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief evaluates the functional
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateFunctional"
PetscErrorCode RegularizationRegistrationH1SN::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr = 0;
    ScalarType  *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_gv11 = NULL, *p_gv12 = NULL, *p_gv13 = NULL,
                *p_gv21 = NULL, *p_gv22 = NULL, *p_gv23 = NULL,
                *p_gv31 = NULL, *p_gv32 = NULL, *p_gv33 = NULL;
    double timer[5] = {0, 0, 0, 0, 0};
    std::bitset<3>XYZ = 0; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    ScalarType beta, value;
    PetscFunctionBegin;

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization parameter
    beta = this->m_Opt->GetRegNorm().beta[0];

    *R = 0.0;

    // if regularization weight is zero, do noting
    if (beta != 0.0) {
        if (this->m_WorkVecField == NULL) {
            try{this->m_WorkVecField = new VecField(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        // get arrays
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

        // X1 gradient
        ierr = this->m_WorkVecField->GetArrays(p_gv11, p_gv12, p_gv13); CHKERRQ(ierr);
        accfft_grad(p_gv11, p_gv12, p_gv13, p_v1, this->m_Opt->GetFFT().plan,&XYZ,timer);
        ierr = this->m_WorkVecField->RestoreArrays(p_gv11, p_gv12, p_gv13); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 4);

        // compute inner products
        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value); *R +=value;


        // X2 gradient
        ierr = this->m_WorkVecField->GetArrays(p_gv21, p_gv22, p_gv23); CHKERRQ(ierr);
        accfft_grad(p_gv21, p_gv22, p_gv23, p_v2, this->m_Opt->GetFFT().plan, &XYZ, timer);
        ierr = this->m_WorkVecField->RestoreArrays(p_gv21, p_gv22, p_gv23); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,4);

        // compute inner products
        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value); *R +=value;


        // X3 gradient
        ierr = this->m_WorkVecField->GetArrays(p_gv31, p_gv32, p_gv33); CHKERRQ(ierr);
        accfft_grad(p_gv31, p_gv32, p_gv33, p_v3, this->m_Opt->GetFFT().plan, &XYZ, timer);
        ierr = this->m_WorkVecField->RestoreArrays(p_gv31, p_gv32, p_gv33); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT, 4);

        ierr = VecRestoreArray(this->m_WorkVecField->m_X1, &p_gv31); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_WorkVecField->m_X2, &p_gv32); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_WorkVecField->m_X3, &p_gv33); CHKERRQ(ierr);

        // compute inner products
        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &value); *R +=value;
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &value); *R +=value;

        ierr = v->RestoreArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);

        // multiply with regularization weight
        *R = 0.5*beta*(*R);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timer);
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief evaluates first variation of regularization norm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradient"
PetscErrorCode RegularizationRegistrationH1SN::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_Lv1 = NULL, *p_Lv2 = NULL, *p_Lv3 = NULL;
    ScalarType beta, scale;
    double timer[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegNorm().beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = dvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        ierr = this->Allocate(0); CHKERRQ(ierr);
        ierr = this->Allocate(1); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        // compute forward fft
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v1, this->m_v1hat, timer);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v2, this->m_v2hat, timer);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v3, this->m_v3hat, timer);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,3);

        scale = this->m_Opt->ComputeFFTScale();

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik, regop;
        IntType i, i1, i2, i3;
#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1) {
            for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2) {
                for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3) {
                    w[0] = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                    w[1] = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                    w[2] = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                    CheckWaveNumbers(w,nx);

                    // compute bilaplacian operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = -scale*beta*lapik;

                    // get linear index
                    i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetFFT().osize);

                    // apply to individual components
                    this->m_Lv1hat[i][0] = regop*this->m_v1hat[i][0];
                    this->m_Lv1hat[i][1] = regop*this->m_v1hat[i][1];

                    this->m_Lv2hat[i][0] = regop*this->m_v2hat[i][0];
                    this->m_Lv2hat[i][1] = regop*this->m_v2hat[i][1];

                    this->m_Lv3hat[i][0] = regop*this->m_v3hat[i][0];
                    this->m_Lv3hat[i][1] = regop*this->m_v3hat[i][1];
                }
            }
        }
}  // pragma omp parallel


        // compute inverse fft
        ierr = dvR->GetArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv1hat, p_Lv1, timer);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv2hat, p_Lv2, timer);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv3hat, p_Lv3, timer);
        ierr = dvR->RestoreArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT, 3);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timer);

    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies second variation of regularization norm to
 * a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode RegularizationRegistrationH1SN::HessianMatVec(VecField* dvvR, VecField* vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType beta;
    PetscFunctionBegin;

    ierr = Assert(vtilde != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(dvvR != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegNorm().beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = dvvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        ierr = this->EvaluateGradient(dvvR, vtilde); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvOp"
PetscErrorCode RegularizationRegistrationH1SN::ApplyInvOp(VecField* Ainvx, VecField* x, bool applysqrt) {
    PetscErrorCode ierr;
    int nx[3];
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL,
                *p_Lv1 = NULL, *p_Lv2 = NULL, *p_Lv3 = NULL;
    ScalarType beta, scale;
    double timer[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(Ainvx != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegNorm().beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = VecCopy(x->m_X1, Ainvx->m_X1); CHKERRQ(ierr);
        ierr = VecCopy(x->m_X2, Ainvx->m_X2); CHKERRQ(ierr);
        ierr = VecCopy(x->m_X3, Ainvx->m_X3); CHKERRQ(ierr);
    } else {
        ierr = this->Allocate(0); CHKERRQ(ierr);
        ierr = this->Allocate(1); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        // compute forward fft
        ierr = x->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_x1, this->m_v1hat, timer);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_x2, this->m_v2hat, timer);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_x3, this->m_v3hat, timer);
        ierr = x->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT, 3);
        scale = this->m_Opt->ComputeFFTScale();

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik, regop;
        IntType i, i1, i2, i3;
#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1) {
            for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2) {
                for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3) {
                    w[0] = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                    w[1] = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                    w[2] = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                    CheckWaveNumbersInv(w, nx);

                    // compute bilaplacian operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = (fabs(lapik) == 0.0) ? beta : -beta*lapik;

                    if (applysqrt) regop = sqrt(regop);
                    regop = scale/regop;

                    i=GetLinearIndex(i1, i2, i3, this->m_Opt->GetFFT().osize);

                    // apply to individual components
                    this->m_Lv1hat[i][0] = regop*this->m_v1hat[i][0];
                    this->m_Lv1hat[i][1] = regop*this->m_v1hat[i][1];

                    this->m_Lv2hat[i][0] = regop*this->m_v2hat[i][0];
                    this->m_Lv2hat[i][1] = regop*this->m_v2hat[i][1];

                    this->m_Lv3hat[i][0] = regop*this->m_v3hat[i][0];
                    this->m_Lv3hat[i][1] = regop*this->m_v3hat[i][1];
                }
            }
        }
}  // pragma omp parallel

        // compute inverse fft
        ierr = Ainvx->GetArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv1hat, p_Lv1, timer);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv2hat, p_Lv2, timer);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv3hat, p_Lv3, timer);
        ierr = Ainvx->RestoreArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);

        // increment counter
        this->m_Opt->IncrementCounter(FFT, 3);

        // increase fft timer
        this->m_Opt->IncreaseFFTTimers(timer);
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetExtremeEigValsInvOp"
PetscErrorCode RegularizationRegistrationH1SN::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr = 0;
    ScalarType w[3], beta1, beta2, regop;

    PetscFunctionBegin;

    beta1=this->m_Opt->GetRegNorm().beta[0];
    beta2=this->m_Opt->GetRegNorm().beta[1];

    // get max value
    w[0] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[0])/2.0;
    w[1] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[1])/2.0;
    w[2] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[2])/2.0;

    // compute largest value for operator
    regop = -(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]); // laplacian
    regop = -beta1*regop + beta2; // beta_1*laplacian + beta_2
    emin = 1.0/regop;
    emax = 1.0/beta2; // 1/(\beta_1*0 + \beta_2)

    PetscFunctionReturn(ierr);
}




} // end of name space




#endif //_REGULARIZATIONREGISTRATIONH1SN_CPP_
