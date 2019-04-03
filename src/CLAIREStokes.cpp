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
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _CLAIRESTOKES_CPP_
#define _CLAIRESTOKES_CPP_

#include <math.h>
#include "CLAIREStokes.hpp"

namespace reg {

/********************************************************************
 * @brief default constructor
 *******************************************************************/
CLAIREStokes::CLAIREStokes() : SuperClass() {
    this->Initialize();
}

/********************************************************************
 * @brief default destructor
 *******************************************************************/
CLAIREStokes::~CLAIREStokes() {
    this->ClearMemory();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
CLAIREStokes::CLAIREStokes(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode CLAIREStokes::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode CLAIREStokes::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIREStokes::SolveAdjointEquationSL() {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt, ll, lm, llnext;
    ScalarType *p_l = NULL,  *p_m=NULL,
                *p_vec1 = NULL, *p_vec2 = NULL, *p_vec3 = NULL,
                *p_b1 = NULL, *p_b2 = NULL, *p_b3 = NULL;
    ScalarType lambda, ht, scale;
    bool fullnewton = false;
    double timer[NFFTTIMERS] = {0};
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;

    PetscFunctionBegin;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ht = this->m_Opt->GetTimeStepSize();
    scale = ht;

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    // compute trajectory
    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    
    ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
    
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        fullnewton = true;
    }

    ierr = VecGetArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(*this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_vec1, p_vec2, p_vec3); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j) {  // for all time points
        lm = (nt-j)*nc*nl;
        if (fullnewton) {
            ll = (nt-j)*nc*nl; llnext = (nt-(j+1))*nc*nl;
        } else {
            ll = 0; llnext = 0;
        }

        // scaling for trapezoidal rule (for body force)
        if (j == 0) scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            // compute gradient of m
            /*this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_vec1, p_vec2, p_vec3, p_m+lm+k*nl, this->m_Opt->m_FFT.plan, &xyz, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);*/
            this->m_Differentiation->Gradient(p_vec1, p_vec2, p_vec3, p_m+lm+k*nl);


            // compute body force
            for (IntType i = 0; i < nl; ++i) {
                lambda = p_l[ll+k*nl+i];
                p_b1[i] += scale*p_vec1[i]*lambda/static_cast<ScalarType>(nc);
                p_b2[i] += scale*p_vec2[i]*lambda/static_cast<ScalarType>(nc);
                p_b3[i] += scale*p_vec3[i]*lambda/static_cast<ScalarType>(nc);
            }

            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_l+llnext+k*nl, p_l+ll+k*nl, "adjoint"); CHKERRQ(ierr);
        }  // for all image components
        // trapezoidal rule (revert scaling; for body force)
        if (j == 0) scale *= 2.0;
    }  // for all time points


    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ll = k*nl; lm = k*nl;

        // compute gradient of m (for incremental body force)
        /*this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_vec1, p_vec2, p_vec3, p_m+lm, this->m_Opt->m_FFT.plan, &xyz, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);*/
        this->m_Differentiation->Gradient(p_vec1, p_vec2, p_vec3, p_m+lm);


        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            lambda = p_l[ll+i];
            // compute bodyforce
            p_b1[i] += 0.5*scale*p_vec1[i]*lambda/static_cast<ScalarType>(nc);
            p_b2[i] += 0.5*scale*p_vec2[i]*lambda/static_cast<ScalarType>(nc);
            p_b3[i] += 0.5*scale*p_vec3[i]*lambda/static_cast<ScalarType>(nc);
        }
    }

    ierr = this->m_WorkVecField1->RestoreArrays(p_vec1, p_vec2, p_vec3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);
    ierr = VecRestoreArray(*this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = VecRestoreArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIREStokes::SolveIncAdjointEquationGNSL() {
    PetscErrorCode ierr = 0;
    IntType nl, nt, nc, lm, ll;
    ScalarType *p_ltilde = NULL, *p_m = NULL,
                *p_btilde1 = NULL, *p_btilde2 = NULL, *p_btilde3 = NULL,
                *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL;
    ScalarType ht, scale, ltilde;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    double timer[NFFTTIMERS] = {0};
    PetscFunctionBegin;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ht = this->m_Opt->GetTimeStepSize();
    scale = ht;

    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);

    if (this->m_SemiLagrangianMethod == NULL) {
        ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }

    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);

    // get variables
    ierr = VecGetArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(*this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_btilde1, p_btilde2, p_btilde3); CHKERRQ(ierr);

    // do numerical time integration
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        lm = (nt-j)*nc*nl;

        if (j == 0) scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            ll = k*nl;

            // compute gradient of m (for incremental body force)
            /*this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m+lm, this->m_Opt->m_FFT.plan, &xyz, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);*/
            this->m_Differentiation->Gradient(p_gradm1, p_gradm2, p_gradm3, p_m+lm);

            // compute incremental bodyforce
            for (IntType i = 0; i < nl; ++i) {
                ltilde = p_ltilde[ll+i];    // get \tilde{\lambda}(x)
                p_btilde1[i] += scale*p_gradm1[i]*ltilde/static_cast<ScalarType>(nc);
                p_btilde2[i] += scale*p_gradm2[i]*ltilde/static_cast<ScalarType>(nc);
                p_btilde3[i] += scale*p_gradm3[i]*ltilde/static_cast<ScalarType>(nc);
            }
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_ltilde+ll, p_ltilde+ll, "adjoint"); CHKERRQ(ierr);
        }  // for all image components
        if (j == 0) scale *= 2.0;
    }  // for all time points


    // incremental compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ll = k*nl; lm = k*nl;

        // compute gradient of m (for incremental body force)
        /*this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m+lm, this->m_Opt->m_FFT.plan, &xyz, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);*/
        this->m_Differentiation->Gradient(p_gradm1, p_gradm2, p_gradm3, p_m+lm);


        // compute incremental bodyforce
        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            ltilde = p_ltilde[ll+i];
            p_btilde1[i] += 0.5*scale*p_gradm1[i]*ltilde/static_cast<ScalarType>(nc);
            p_btilde2[i] += 0.5*scale*p_gradm2[i]*ltilde/static_cast<ScalarType>(nc);
            p_btilde3[i] += 0.5*scale*p_gradm3[i]*ltilde/static_cast<ScalarType>(nc);
        }
    }

    // restore variables
    ierr = this->m_WorkVecField2->RestoreArrays(p_btilde1, p_btilde2, p_btilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = VecRestoreArray(*this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = VecRestoreArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief apply projection to map \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
PetscErrorCode CLAIREStokes::ApplyProjection() {
    PetscErrorCode ierr = 0;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL, scale;
    long int nx[3];
    double applytime;
    ComplexType x1hat, x2hat, x3hat;
    double timer[NFFTTIMERS] = {0};


    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    nx[0] = static_cast<long int>(this->m_Opt->m_Domain.nx[0]);
    nx[1] = static_cast<long int>(this->m_Opt->m_Domain.nx[1]);
    nx[2] = static_cast<long int>(this->m_Opt->m_Domain.nx[2]);

    scale = this->m_Opt->ComputeFFTScale();

    // allocate fields for spectral operations
    ierr = this->SetupSpectralData(); CHKERRQ(ierr);

    // copy input
    ierr = this->m_WorkVecField1->Copy(this->m_WorkVecField2); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    // compute forward fft
    this->m_Opt->StartTimer(FFTSELFEXEC);
    ierr = this->m_Opt->m_FFT.fft->FFT_R2C(p_x1, this->m_x1hat); CHKERRQ(ierr);
    ierr = this->m_Opt->m_FFT.fft->FFT_R2C(p_x2, this->m_x2hat); CHKERRQ(ierr);
    ierr = this->m_Opt->m_FFT.fft->FFT_R2C(p_x3, this->m_x3hat); CHKERRQ(ierr);
    //accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_x1, this->m_x1hat, timer);
    //accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_x2, this->m_x2hat, timer);
    //accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, p_x3, this->m_x3hat, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 3);

    applytime = -MPI_Wtime();
#pragma omp parallel
{
    long int x1, x2, x3, wx1, wx2, wx3;
    ScalarType lapinvik, gradik1, gradik2, gradik3;
    IntType i, i1, i2, i3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->m_FFT.osize[0]; ++i1) {
        for (i2 = 0; i2 < this->m_Opt->m_FFT.osize[1]; ++i2) {
            for (i3 = 0; i3 < this->m_Opt->m_FFT.osize[2]; ++i3) {
                x1 = static_cast<long int>(i1 + this->m_Opt->m_FFT.ostart[0]);
                x2 = static_cast<long int>(i2 + this->m_Opt->m_FFT.ostart[1]);
                x3 = static_cast<long int>(i3 + this->m_Opt->m_FFT.ostart[2]);

                // set wavenumber
                wx1 = x1;
                wx2 = x2;
                wx3 = x3;

                if (x1 > nx[0]/2) wx1 -= nx[0];
                if (x2 > nx[1]/2) wx2 -= nx[1];
                if (x3 > nx[2]/2) wx3 -= nx[2];

                // compute inverse laplacian operator
                lapinvik = static_cast<ScalarType>(wx1*wx1 + wx2*wx2 + wx3*wx3);
                //lapinvik = round(lapinvik) == 0.0 ? -1.0 : 1.0/lapinvik;
                lapinvik = lapinvik == 0.0 ? -1.0 : -1.0/lapinvik;

                if (x1 == nx[0]/2) wx1 = 0;
                if (x2 == nx[1]/2) wx2 = 0;
                if (x3 == nx[2]/2) wx3 = 0;

                // compute gradient operator
                gradik1 = static_cast<ScalarType>(wx1);
                gradik2 = static_cast<ScalarType>(wx2);
                gradik3 = static_cast<ScalarType>(wx3);

                i = GetLinearIndex(i1, i2, i3, this->m_Opt->m_FFT.osize);

                x1hat[0] = this->m_x1hat[i][0];
                x1hat[1] = this->m_x1hat[i][1];

                x2hat[0] = this->m_x2hat[i][0];
                x2hat[1] = this->m_x2hat[i][1];

                x3hat[0] = this->m_x3hat[i][0];
                x3hat[1] = this->m_x3hat[i][1];

                // compute div(b)
                this->m_x1hat[i][0] = -scale*(gradik1*x1hat[0]
                                                     + gradik2*x2hat[0]
                                                     + gradik3*x3hat[0]);

                this->m_x1hat[i][1] =  scale*(gradik1*x1hat[1]
                                                    + gradik2*x2hat[1]
                                                    + gradik3*x3hat[1]);

                // compute lap^{-1} div(b)
                this->m_x1hat[i][0] *= lapinvik;
                this->m_x1hat[i][1] *= lapinvik;

                // compute x2 gradient of lab^{-1} div(b)
                this->m_x2hat[i][0] = -gradik2*this->m_x1hat[i][0];
                this->m_x2hat[i][1] =  gradik2*this->m_x1hat[i][1];

                // compute x3 gradient of lab^{-1} div(b)
                this->m_x3hat[i][0] = -gradik3*this->m_x1hat[i][0];
                this->m_x3hat[i][1] =  gradik3*this->m_x1hat[i][1];

                // compute x1 gradient of lab^{-1} div(b)
                this->m_x1hat[i][0] *= -gradik1;
                this->m_x1hat[i][1] *=  gradik1;
            }
        }
    }
}  // pragma omp parallel
    applytime += MPI_Wtime();
    timer[FFTHADAMARD] += applytime;

    // compute inverse fft
    this->m_Opt->StartTimer(FFTSELFEXEC);
    ierr = this->m_Opt->m_FFT.fft->FFT_C2R(this->m_x1hat, p_x1); CHKERRQ(ierr);
    ierr = this->m_Opt->m_FFT.fft->FFT_C2R(this->m_x2hat, p_x2); CHKERRQ(ierr);
    ierr = this->m_Opt->m_FFT.fft->FFT_C2R(this->m_x3hat, p_x3); CHKERRQ(ierr);
    //accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_x1hat, p_x1, timer);
    //accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_x2hat, p_x2, timer);
    //accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, this->m_x3hat, p_x3, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 3);

    ierr = this->m_WorkVecField1->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->AXPY(1.0, this->m_WorkVecField1); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _CLAIRESTOKES_CPP_
