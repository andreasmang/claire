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

#ifndef _REGULARIZATIONREGISTRATIONH3SN_CPP_
#define _REGULARIZATIONREGISTRATIONH3SN_CPP_

#include "RegularizationRegistrationH3SN.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegularizationRegistrationH3SN::RegularizationRegistrationH3SN() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegularizationRegistrationH3SN::~RegularizationRegistrationH3SN(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
RegularizationRegistrationH3SN::RegularizationRegistrationH3SN(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief evaluates the functional
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH3SN::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_bv1 = NULL, *p_bv2 = NULL, *p_bv3 = NULL;
    ScalarType beta, ipxi, scale;
    IntType nx[3];
    double timer[NFFTTIMERS] = {0}, applytime;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // get regularization weight
    beta = this->m_Opt->m_RegNorm.beta[0];

    ierr = Assert(this->m_v1hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v2hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v3hat != NULL, "null pointer"); CHKERRQ(ierr);

    *R = 0.0;

    // if regularization weight is zero, do noting
    if (beta != 0.0) {
        ierr = Assert(v != NULL,"null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_WorkVecField != NULL, "null pointer"); CHKERRQ(ierr);

        nx[0] = this->m_Opt->m_Domain.nx[0];
        nx[1] = this->m_Opt->m_Domain.nx[1];
        nx[2] = this->m_Opt->m_Domain.nx[2];

        scale = this->m_Opt->ComputeFFTScale();

        // compute forward fft
        this->m_Opt->StartTimer(FFTSELFEXEC);
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_v1, this->m_v1hat, timer);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_v2, this->m_v2hat, timer);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_v3, this->m_v3hat, timer);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

        applytime = -MPI_Wtime();
#pragma omp parallel
{
        ScalarType lapik, regop[6], gradik[3];
        IntType i, i1, i2, i3, w[3];
#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->m_FFT.osize[0]; ++i1){
            for (i2 = 0; i2 < this->m_Opt->m_FFT.osize[1]; ++i2){
                for (i3 = 0; i3 < this->m_Opt->m_FFT.osize[2]; ++i3){
                    w[0] = i1 + this->m_Opt->m_FFT.ostart[0];
                    w[1] = i2 + this->m_Opt->m_FFT.ostart[1];
                    w[2] = i3 + this->m_Opt->m_FFT.ostart[2];

                    ComputeWaveNumber(w, nx);

                    // compute bilaplacian operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute gradient operator
                    gradik[0] = static_cast<ScalarType>(w[0]);
                    gradik[1] = static_cast<ScalarType>(w[1]);
                    gradik[2] = static_cast<ScalarType>(w[2]);

                    // compute regularization operator
                    regop[0] =  scale*gradik[0]*lapik;
                    regop[1] = -scale*gradik[0]*lapik;

                    regop[2] =  scale*gradik[1]*lapik;
                    regop[3] = -scale*gradik[1]*lapik;

                    regop[4] =  scale*gradik[2]*lapik;
                    regop[5] = -scale*gradik[2]*lapik;

                    i = GetLinearIndex(i1, i2, i3, this->m_Opt->m_FFT.osize);

                    // apply to individual components
                    this->m_v1hat[i][0] *= regop[0];
                    this->m_v1hat[i][1] *= regop[1];

                    this->m_v2hat[i][0] *= regop[2];
                    this->m_v2hat[i][1] *= regop[3];

                    this->m_v3hat[i][0] *= regop[4];
                    this->m_v3hat[i][1] *= regop[5];
                }
            }
        }
}  // pragma omp parallel
        applytime += MPI_Wtime();
        timer[FFTHADAMARD] += applytime;

        // compute inverse fft
        ierr = this->m_WorkVecField->GetArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v1hat, p_bv1, timer);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v2hat, p_bv2, timer);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v3hat, p_bv3, timer);
        ierr = this->m_WorkVecField->RestoreArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, 3);

        // compute inner product
        ierr=VecTDot(this->m_WorkVecField->m_X1, this->m_WorkVecField->m_X1, &ipxi); CHKERRQ(ierr); *R += ipxi;
        ierr=VecTDot(this->m_WorkVecField->m_X2, this->m_WorkVecField->m_X2, &ipxi); CHKERRQ(ierr); *R += ipxi;
        ierr=VecTDot(this->m_WorkVecField->m_X3, this->m_WorkVecField->m_X3, &ipxi); CHKERRQ(ierr); *R += ipxi;

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timer);

        // multiply with regularization weight
        *R *= 0.5*beta;
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief evaluates first variation of regularization norm
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH3SN::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr;
    IntType nx[3];
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_bv1 = NULL, *p_bv2 = NULL, *p_bv3 = NULL;
    ScalarType beta, scale;
    double timer[NFFTTIMERS] = {0}, applytime;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = Assert(this->m_v1hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v2hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v3hat != NULL, "null pointer"); CHKERRQ(ierr);

    // get regularization weight
    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0) {
        ierr = dvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        nx[0] = this->m_Opt->m_Domain.nx[0];
        nx[1] = this->m_Opt->m_Domain.nx[1];
        nx[2] = this->m_Opt->m_Domain.nx[2];

        scale = this->m_Opt->ComputeFFTScale();

        // compute forward fft
        this->m_Opt->StartTimer(FFTSELFEXEC);
        ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_v1, this->m_v1hat, timer);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_v2, this->m_v2hat, timer);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_v3, this->m_v3hat, timer);
        ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

        applytime = -MPI_Wtime();
#pragma omp parallel
{
        ScalarType trihik, regop;
        IntType i, i1, i2, i3, w[3];
#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->m_FFT.osize[0]; ++i1){
            for (i2 = 0; i2 < this->m_Opt->m_FFT.osize[1]; ++i2){
                for (i3 = 0; i3 < this->m_Opt->m_FFT.osize[2]; ++i3){
                    w[0] = i1 + this->m_Opt->m_FFT.ostart[0];
                    w[1] = i2 + this->m_Opt->m_FFT.ostart[1];
                    w[2] = i3 + this->m_Opt->m_FFT.ostart[2];

                    ComputeWaveNumber(w, nx);

                    // compute bilaplacian operator
                    trihik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                    trihik = std::pow(trihik, 3);

                    // compute regularization operator
                    regop = -scale*beta*trihik;

                    // get linear index
                    i = GetLinearIndex(i1, i2, i3, this->m_Opt->m_FFT.osize);

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
        ierr = dvR->GetArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v1hat, p_bv1, timer);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v2hat, p_bv2, timer);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v3hat, p_bv3, timer);
        ierr = dvR->RestoreArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, 3);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timer);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies second variation of regularization norm to
 * a vector
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH3SN::HessianMatVec(VecField* dvvR, VecField* vtilde) {
    PetscErrorCode ierr;
    ScalarType beta;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr=Assert(vtilde != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvvR != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0){
        ierr = dvvR->SetValue(0.0); CHKERRQ(ierr);
    } else {
        ierr = this->EvaluateGradient(dvvR, vtilde); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH3SN::ApplyInverse(VecField* Ainvx, VecField* x, bool applysqrt) {
    PetscErrorCode ierr;
    IntType nx[3];
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL,
                *p_bv1 = NULL, *p_bv2 = NULL, *p_bv3 = NULL;
    ScalarType beta, scale;
    double timer[NFFTTIMERS] = {0};
    double applytime;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr=Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(Ainvx != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = Assert(this->m_v1hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v2hat != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_v3hat != NULL, "null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->m_RegNorm.beta[0];

    // if regularization weight is zero, do noting
    if (beta == 0.0){
        ierr=VecCopy(x->m_X1, Ainvx->m_X1); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X2, Ainvx->m_X2); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X3, Ainvx->m_X3); CHKERRQ(ierr);
    } else {
        nx[0] = this->m_Opt->m_Domain.nx[0];
        nx[1] = this->m_Opt->m_Domain.nx[1];
        nx[2] = this->m_Opt->m_Domain.nx[2];

        scale = this->m_Opt->ComputeFFTScale();

        // compute forward fft
        this->m_Opt->StartTimer(FFTSELFEXEC);
        ierr = x->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_x1, this->m_v1hat, timer);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_x2, this->m_v2hat, timer);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_x3, this->m_v3hat, timer);
        ierr = x->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
        this->m_Opt->IncrementCounter(FFT, 3);

        applytime = -MPI_Wtime();
#pragma omp parallel
{
        ScalarType lapik, regop;
        IntType i, i1, i2, i3, w[3];
#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->m_FFT.osize[0]; ++i1) {
            for (i2 = 0; i2 < this->m_Opt->m_FFT.osize[1]; ++i2) {
                for (i3 = 0; i3 < this->m_Opt->m_FFT.osize[2]; ++i3) {
                    w[0] = i1 + this->m_Opt->m_FFT.ostart[0];
                    w[1] = i2 + this->m_Opt->m_FFT.ostart[1];
                    w[2] = i3 + this->m_Opt->m_FFT.ostart[2];

                    ComputeWaveNumber(w, nx);

                    // compute regularization operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                    regop = (std::abs(lapik) == 0.0) ? beta : -beta*std::pow(lapik, 3);

                    if (applysqrt) regop = sqrt(regop);
                    regop = scale/regop;

                    i = GetLinearIndex(i1, i2, i3, this->m_Opt->m_FFT.osize);

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
        ierr = Ainvx->GetArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v1hat, p_bv1, timer);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v2hat, p_bv2, timer);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_v3hat, p_bv3, timer);
        ierr = Ainvx->RestoreArrays(p_bv1, p_bv2, p_bv3); CHKERRQ(ierr);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, 3);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timer);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
PetscErrorCode RegularizationRegistrationH3SN::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
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
