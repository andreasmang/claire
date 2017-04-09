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

#ifndef _REGULARIZATIONREGISTRATIONH2SN_CPP_
#define _REGULARIZATIONREGISTRATIONH2SN_CPP_

#include "RegularizationRegistrationH2SN.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegularizationRegistrationH2SN::RegularizationRegistrationH2SN() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegularizationRegistrationH2SN::~RegularizationRegistrationH2SN(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
RegularizationRegistrationH2SN::RegularizationRegistrationH2SN(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief evaluates the functional (we have to promote everything
 * to double to be able to solve the problem accurately; we loose
 * too many digits here)
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH2SN::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_bv1 = NULL, *p_bv2 = NULL, *p_bv3 = NULL;
    ScalarType sqrtbeta, ipxi, scale, value;
    double applytime;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v1hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v2hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v3hat != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization weight
    sqrtbeta = static_cast<double>(sqrt(this->m_Opt->GetRegNorm().beta[0]));

    *R = 0.0; value = 0.0;

    // if regularization weight is zero, do noting
    if (sqrtbeta != 0.0) {
        ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_WorkVecField != NULL, "null pointer"); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        scale = static_cast<double>(this->m_Opt->ComputeFFTScale());

        // compute forward fft
        this->m_Opt->StartTimer(FFTSELFEXEC);
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v1, this->m_v1hat, timer);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v2, this->m_v2hat, timer);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v3, this->m_v3hat, timer);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

        applytime = -MPI_Wtime();
#pragma omp parallel
{
        long int w[3];
        double lapik, regop;
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
                    lapik = -static_cast<double>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetFFT().osize);

                    // compute regularization operator
                    regop = scale*sqrtbeta*lapik;

                    // apply to individual components
                    this->m_v1hat[i][0] *= regop;
                    this->m_v1hat[i][1] *= regop;

                    this->m_v2hat[i][0] *= regop;
                    this->m_v2hat[i][1] *= regop;

                    this->m_v3hat[i][0] *= regop;
                    this->m_v3hat[i][1] *= regop;
                }
            }
        }
}// pragma omp parallel
        applytime += MPI_Wtime();
        timer[FFTHADAMARD] += applytime;

        // compute inverse fft
        ierr = this->m_WorkVecField->GetArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v1hat, p_bv1, timer);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v2hat, p_bv2, timer);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v3hat, p_bv3, timer);
        ierr = this->m_WorkVecField->RestoreArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, 3);

        // compute inner product
        ierr = VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &ipxi); CHKERRQ(ierr); value += static_cast<double>(ipxi);
        ierr = VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &ipxi); CHKERRQ(ierr); value += static_cast<double>(ipxi);
        ierr = VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &ipxi); CHKERRQ(ierr); value += static_cast<double>(ipxi);

        // multiply with regularization weight
        *R = static_cast<ScalarType>(0.5*value);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timer);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates first variation of regularization norm
 * @param[in] v velocity field
 * @param[out] dvR gradient of regularization (evaluated)
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH2SN::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType beta, scale;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_bv1 = NULL, *p_bv2 = NULL, *p_bv3 = NULL;
    double applytime;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v1hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v2hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v3hat != NULL, "null pointer"); CHKERRQ(ierr);

    beta = static_cast<double>(this->m_Opt->GetRegNorm().beta[0]);

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = VecSet(dvR->m_X1, 0.0); CHKERRQ(ierr);
        ierr = VecSet(dvR->m_X2, 0.0); CHKERRQ(ierr);
        ierr = VecSet(dvR->m_X3, 0.0); CHKERRQ(ierr);
    } else {
        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        scale = static_cast<double>(this->m_Opt->ComputeFFTScale());

        // compute forward fft
        this->m_Opt->StartTimer(FFTSELFEXEC);
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v1, this->m_v1hat, timer);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v2, this->m_v2hat, timer);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v3, this->m_v3hat, timer);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

        applytime = -MPI_Wtime();
#pragma omp parallel
{
        long int w[3];
        double lapik, regop;
        IntType i;
#pragma omp for
        for (IntType i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1) {
            for (IntType i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2) {
                for (IntType i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3) {
                    w[0] = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                    w[1] = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                    w[2] = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                    CheckWaveNumbers(w, nx);

                    // compute bilaplacian operator
                    lapik = -static_cast<double>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = scale*beta*(lapik*lapik);

                    // get linear index
                    i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetFFT().osize);

                    // apply to individual components
                    this->m_v1hat[i][0] *= regop;
                    this->m_v1hat[i][1] *= regop;

                    this->m_v2hat[i][0] *= regop;
                    this->m_v2hat[i][1] *= regop;

                    this->m_v3hat[i][0] *= regop;
                    this->m_v3hat[i][1] *= regop;
                }
            }
        }
}// pragma omp parallel
        applytime += MPI_Wtime();
        timer[FFTHADAMARD] += applytime;

        // compute inverse fft
        ierr = dvR->GetArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v1hat, p_bv1, timer);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v2hat, p_bv2, timer);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v3hat, p_bv3, timer);
        ierr = dvR->RestoreArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, 3);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timer);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies second variation of regularization norm to vector
 * @param dvvR regularization operator applied to vector \tilde{v}
 * @param vtilde incremental vector field \tilde{v}
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH2SN::HessianMatVec(VecField* dvvR, VecField* vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType beta;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(dvvR != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vtilde != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegNorm().beta[0];

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
PetscErrorCode RegularizationRegistrationH2SN::ApplyInvOp(VecField* Ainvv, VecField* v, bool applysqrt) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType beta, scale;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_bv1 = NULL, *p_bv2 = NULL, *p_bv3 = NULL;
    double applytime;

    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(Ainvv != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v1hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v2hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v3hat != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegNorm().beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = VecCopy(v->m_X1, Ainvv->m_X1); CHKERRQ(ierr);
        ierr = VecCopy(v->m_X2, Ainvv->m_X2); CHKERRQ(ierr);
        ierr = VecCopy(v->m_X3, Ainvv->m_X3); CHKERRQ(ierr);
    } else {
        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        scale = static_cast<double>(this->m_Opt->ComputeFFTScale());

        // compute forward fft
        this->m_Opt->StartTimer(FFTSELFEXEC);
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v1, this->m_v1hat, timer);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v2, this->m_v2hat, timer);
        accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_v3, this->m_v3hat, timer);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

        applytime = -MPI_Wtime();
#pragma omp parallel
{
        long int w[3];
        double lapik, regop;
        IntType i;
#pragma omp for
        for (IntType i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1) {
            for (IntType i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2) {
                for (IntType i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3) {
                    w[0] = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                    w[1] = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                    w[2] = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                    CheckWaveNumbersInv(w, nx);

                    // compute bilaplacian operator
                    lapik = -static_cast<double>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = (std::abs(lapik) == 0.0) ? beta : beta*(lapik*lapik);

                    if (applysqrt) regop = std::sqrt(regop);
                    regop = scale/regop;

                    i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetFFT().osize);

                    // apply to individual components
                    this->m_v1hat[i][0] *= regop;
                    this->m_v1hat[i][1] *= regop;

                    this->m_v2hat[i][0] *= regop;
                    this->m_v2hat[i][1] *= regop;

                    this->m_v3hat[i][0] *= regop;
                    this->m_v3hat[i][1] *= regop;

                }
            }
        }
}  // pragma omp parallel
        applytime += MPI_Wtime();
        timer[FFTHADAMARD] += applytime;

        // compute inverse fft
        ierr = Ainvv->GetArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v1hat, p_bv1, timer);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v2hat, p_bv2, timer);
        accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_v3hat, p_bv3, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        ierr = Ainvv->RestoreArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timer);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH2SN::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr = 0;
    ScalarType w[3], beta, regop;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

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

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _REGULARIZATIONREGISTRATIONH2SN_CPP_
