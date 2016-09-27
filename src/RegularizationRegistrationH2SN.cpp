/**
 *  Copyright (c) 2015-2016.
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
 *
 */

#ifndef _REGULARIZATIONREGISTRATIONH2SN_CPP_
#define _REGULARIZATIONREGISTRATIONH2SN_CPP_

#include "RegularizationRegistrationH2SN.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH2SN"
RegularizationRegistrationH2SN::RegularizationRegistrationH2SN() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegularizationRegistrationH2SN"
RegularizationRegistrationH2SN::~RegularizationRegistrationH2SN(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH2SN"
RegularizationRegistrationH2SN::RegularizationRegistrationH2SN(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief evaluates the functional
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateFunctional"
PetscErrorCode RegularizationRegistrationH2SN::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_Lv1 = NULL, *p_Lv2 = NULL, *p_Lv3 = NULL;
    ScalarType sqrtbeta, ipxi, scale;
    double timers[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization weight
    sqrtbeta = sqrt(this->m_Opt->GetRegNorm().beta[0]);

    *R=0.0;

    // if regularization weight is zero, do noting
    if (sqrtbeta != 0.0) {
        ierr = this->Allocate(0); CHKERRQ(ierr);
        ierr = this->Allocate(1); CHKERRQ(ierr);

        if (this->m_WorkVecField == NULL){
            try{this->m_WorkVecField = new VecField(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        scale = this->m_Opt->ComputeFFTScale();

        // compute forward fft
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v1, this->m_v1hat, timers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v2, this->m_v2hat, timers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v3, this->m_v3hat, timers);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

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

                    i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetFFT().osize);

                    // compute regularization operator
                    regop = scale*sqrtbeta*lapik;

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
}// pragma omp parallel

        // compute inverse fft
        ierr = this->m_WorkVecField->GetArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv1hat, p_Lv1, timers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv2hat, p_Lv2, timers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv3hat, p_Lv3, timers);
        ierr = this->m_WorkVecField->RestoreArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,3);

        // restore arrays

        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &ipxi); CHKERRQ(ierr); *R += ipxi;
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &ipxi); CHKERRQ(ierr); *R += ipxi;
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &ipxi); CHKERRQ(ierr); *R += ipxi;

        // multiply with regularization weight
        *R *= 0.5;

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timers);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates first variation of regularization norm
 * @param[in] v velocity field
 * @param[out] dvR gradient of regularization (evaluated)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradient"
PetscErrorCode RegularizationRegistrationH2SN::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType beta, scale;
    double ffttimers[5] = {0, 0, 0, 0, 0};
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_Lv1 = NULL, *p_Lv2 = NULL, *p_Lv3 = NULL;

    PetscFunctionBegin;

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegNorm().beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = VecSet(dvR->m_X1, 0.0); CHKERRQ(ierr);
        ierr = VecSet(dvR->m_X2, 0.0); CHKERRQ(ierr);
        ierr = VecSet(dvR->m_X3, 0.0); CHKERRQ(ierr);
    } else {
        ierr = this->Allocate(0); CHKERRQ(ierr);
        ierr = this->Allocate(1); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        scale = this->m_Opt->ComputeFFTScale();

        // compute forward fft
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v1, this->m_v1hat, ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v2, this->m_v2hat, ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v3, this->m_v3hat, ffttimers);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT, 3);

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik,regop;
        IntType i, i1, i2, i3;
#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1) {
            for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2) {
                for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3) {
                    w[0] = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                    w[1] = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                    w[2] = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                    CheckWaveNumbers(w, nx);

                    // compute bilaplacian operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = scale*beta*(lapik*lapik);

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
}// pragma omp parallel

        // compute inverse fft
        ierr = dvR->GetArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv1hat, p_Lv1, ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv2hat, p_Lv2, ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv3hat, p_Lv3, ffttimers);
        ierr = dvR->RestoreArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT, 3);


        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(ffttimers);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies second variation of regularization norm to vector
 * @param dvvR regularization operator applied to vector \tilde{v}
 * @param vtilde incremental vector field \tilde{v}
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode RegularizationRegistrationH2SN::HessianMatVec(VecField* dvvR, VecField* vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType beta;
    PetscFunctionBegin;

    ierr = Assert(dvvR != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vtilde != NULL, "null pointer"); CHKERRQ(ierr);

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
PetscErrorCode RegularizationRegistrationH2SN::ApplyInvOp(VecField* Ainvv, VecField* v, bool applysqrt) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType beta, scale;
    double timers[5] = {0, 0, 0, 0, 0};
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_Lv1 = NULL, *p_Lv2 = NULL, *p_Lv3 = NULL;

    PetscFunctionBegin;

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(Ainvv != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegNorm().beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = VecCopy(v->m_X1, Ainvv->m_X1); CHKERRQ(ierr);
        ierr = VecCopy(v->m_X2, Ainvv->m_X2); CHKERRQ(ierr);
        ierr = VecCopy(v->m_X3, Ainvv->m_X3); CHKERRQ(ierr);
    } else {
        ierr = this->Allocate(0); CHKERRQ(ierr);
        ierr = this->Allocate(1); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        scale = this->m_Opt->ComputeFFTScale();

        // compute forward fft
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v1, this->m_v1hat, timers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v2, this->m_v2hat, timers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_v3, this->m_v3hat, timers);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,3);

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik,regop;
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
                    regop = (std::abs(lapik) == 0.0) ? beta : beta*(lapik*lapik);

                    if (applysqrt) regop = std::sqrt(regop);
                    regop = scale/regop;

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

}// pragma omp parallel

        // compute inverse fft
        ierr = Ainvv->GetArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv1hat, p_Lv1, timers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv2hat, p_Lv2, timers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Lv3hat, p_Lv3, timers);
        ierr = Ainvv->RestoreArrays(p_Lv1, p_Lv2, p_Lv3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timers);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetExtremeEigValsInvOp"
PetscErrorCode RegularizationRegistrationH2SN::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr = 0;
    ScalarType w[3], beta, regop;

    PetscFunctionBegin;

    beta = this->m_Opt->GetRegNorm().beta[0];

    // get max value
    w[0] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[0])/2.0;
    w[1] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[1])/2.0;
    w[2] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[2])/2.0;

    // compute largest value for operator
    regop = -(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]); // laplacian
    regop = beta*(regop*regop); // beta * biharmonic
    emin = 1.0/regop;
    emax = 1.0; // by definition; it's 1/0

    PetscFunctionReturn(ierr);
}




} // end of name space




#endif //_REGULARIZATIONREGISTRATIONH2SN_CPP_
