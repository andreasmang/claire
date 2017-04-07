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
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _OPTIMALCONTROLREGISTRATIONRELAXEDIC_CPP_
#define _OPTIMALCONTROLREGISTRATIONRELAXEDIC_CPP_

#include "OptimalControlRegistrationRelaxedIC.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
OptimalControlRegistrationRelaxedIC::OptimalControlRegistrationRelaxedIC() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
OptimalControlRegistrationRelaxedIC::~OptimalControlRegistrationRelaxedIC(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
OptimalControlRegistrationRelaxedIC::OptimalControlRegistrationRelaxedIC(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationRelaxedIC::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationRelaxedIC::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates the objective value
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationRelaxedIC::EvaluateObjective(ScalarType* J, Vec v) {
    PetscErrorCode ierr = 0;
    ScalarType D = 0.0, Rv = 0.0, Rw = 0.0, hd;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate velocity field
    if (this->m_VelocityField == NULL) {
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("evaluating objective functional"); CHKERRQ(ierr);
    }

    ierr = this->m_Opt->StartTimer(OBJEXEC); CHKERRQ(ierr);

    // get lebesque measure
    hd = this->m_Opt->GetLebesqueMeasure();

    // set components of velocity field
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // evaluate the L2 distance
    ierr = this->EvaluateDistanceMeasure(&D); CHKERRQ(ierr);

    // evaluate the regularization model
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (!this->m_VelocityIsZero) {
        // evaluate the regularization model for v
        if (this->m_WorkVecField1 == NULL) {
            try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        ierr = this->m_Regularization->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_Regularization->EvaluateFunctional(&Rv, this->m_VelocityField); CHKERRQ(ierr);

        // evaluate the regularization model for w = div(v)
        ierr = this->EvaluteRegFunctionalW(&Rw); CHKERRQ(ierr); CHKERRQ(ierr);
    }

    // add up the contributions
    *J = hd*(D + Rv + Rw);

    // store for access
    this->m_Opt->SetJVal(*J);
    this->m_Opt->SetDVal(hd*D);
    this->m_Opt->SetRVal(hd*(Rv + Rw));

    ierr = this->m_Opt->StopTimer(OBJEXEC); CHKERRQ(ierr);

    this->m_Opt->IncrementCounter(OBJEVAL);

    this->m_Opt->Exit(__func__);
    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute the body force
 * b = K[\int_0^1 \igrad m \lambda d t],
 * where K is an operator that projects v onto the manifold of
 * divergence free velocity fields
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationRelaxedIC::EvaluteRegFunctionalW(ScalarType* Rw) {
    PetscErrorCode ierr = 0;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_gdv1 = NULL, *p_gdv2 = NULL, *p_gdv3 = NULL, *p_divv = NULL;
    ScalarType value, regvalue, betaw;
    double timer[NFFTTIMERS] = {0};
    IntType nl, ng;
    std::bitset<3> XYZ = 0; XYZ[0] = 1, XYZ[1] = 1, XYZ[2] = 1;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    if (this->m_WorkVecField1 == NULL) {
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }

    // get regularization weight
    betaw = this->m_Opt->GetRegNorm().beta[2];

    ierr = VecGetArray(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);

    // compute \idiv(\vect{v})
    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(p_divv, p_v1, p_v2, p_v3, this->m_Opt->GetFFT().plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);


    // compute gradient of div(v)
    ierr = this->m_WorkVecField1->GetArrays(p_gdv1, p_gdv2, p_gdv3); CHKERRQ(ierr);
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_grad_t(p_gdv3, p_gdv2, p_gdv1, p_divv, this->m_Opt->GetFFT().plan, &XYZ, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    ierr = this->m_WorkVecField1->RestoreArrays(p_gdv1, p_gdv2, p_gdv3); CHKERRQ(ierr);
    this->m_Opt->IncrementCounter(FFT, FFTGRAD);

    ierr = VecRestoreArray(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);


    // compute inner products ||\igrad w||_L2 + ||w||_L2
    regvalue = 0.0;
    ierr = VecTDot(this->m_WorkVecField1->m_X1, this->m_WorkVecField1->m_X1, &value); CHKERRQ(ierr); regvalue += value;
    ierr = VecTDot(this->m_WorkVecField1->m_X2, this->m_WorkVecField1->m_X2, &value); CHKERRQ(ierr); regvalue += value;
    ierr = VecTDot(this->m_WorkVecField1->m_X3, this->m_WorkVecField1->m_X3, &value); CHKERRQ(ierr); regvalue += value;
    ierr = VecTDot(this->m_WorkScaField1, this->m_WorkScaField1, &value); CHKERRQ(ierr); regvalue += value;

    // add up contributions
    *Rw = 0.5*betaw*regvalue;

    this->m_Opt->IncreaseFFTTimers(timer);
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute the body force
 * b = K[\int_0^1 \igrad m \lambda d t],
 * where K is an operator that projects v onto the manifold of
 * divergence free velocity fields
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationRelaxedIC::ComputeBodyForce() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // assigned to work vec field 2
    ierr = SuperClass::ComputeBodyForce(); CHKERRQ(ierr);

    ierr = this->ApplyProjection(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute the body force
 * b = K[\int_0^1\igrad\tilde{m}\lambda+\igrad m \tilde{\lambda} dt]
 * where K is an operator that projects \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationRelaxedIC::ComputeIncBodyForce() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    // assigned to work vec field 2
    ierr = SuperClass::ComputeIncBodyForce(); CHKERRQ(ierr);

    ierr = this->ApplyProjection(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply projection to map \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationRelaxedIC::ApplyProjection() {
    PetscErrorCode ierr = 0;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL;
    ScalarType beta[3], scale;
    long int nx[3];
    IntType nalloc;
    double applytime;
    ComplexType x1hat, x2hat, x3hat;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    nx[0] = static_cast<long int>(this->m_Opt->GetNumGridPoints(0));
    nx[1] = static_cast<long int>(this->m_Opt->GetNumGridPoints(1));
    nx[2] = static_cast<long int>(this->m_Opt->GetNumGridPoints(2));

    scale = this->m_Opt->ComputeFFTScale();

    // allocate spectral data
    ierr = this->AllocateSpectralData(); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->Copy(this->m_WorkVecField2); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    // compute forward fft
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_x1, this->m_x1hat, timer);
    accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_x2, this->m_x2hat, timer);
    accfft_execute_r2c_t(this->m_Opt->GetFFT().plan, p_x3, this->m_x3hat, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 3);

    beta[0] = this->m_Opt->GetRegNorm().beta[0];
    beta[2] = this->m_Opt->GetRegNorm().beta[2];

    applytime = -MPI_Wtime();
#pragma omp parallel
{
    long int x1, x2, x3, wx1, wx2, wx3;
    ScalarType lapik, lapinvik, gradik1, gradik2, gradik3, opik;
    long int i;
    IntType i1, i2, i3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1) {
        for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2) {
            for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3) {
                x1 = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                x2 = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                x3 = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                // set wavenumber
                wx1 = x1;
                wx2 = x2;
                wx3 = x3;

                if (x1 > nx[0]/2) wx1 -= nx[0];
                if (x2 > nx[1]/2) wx2 -= nx[1];
                if (x3 > nx[2]/2) wx3 -= nx[2];

                // compute inverse laplacian operator
                lapik = -static_cast<ScalarType>(wx1*wx1 + wx2*wx2 + wx3*wx3);

                //lapinvik = round(lapinvik) == 0.0 ? -1.0 : 1.0/lapinvik;
                lapinvik = lapik == 0.0 ? -1.0 : 1.0/lapik;

                if (x1 == nx[0]/2) wx1 = 0;
                if (x2 == nx[1]/2) wx2 = 0;
                if (x3 == nx[2]/2) wx3 = 0;

                // compute gradient operator
                gradik1 = static_cast<ScalarType>(wx1);
                gradik2 = static_cast<ScalarType>(wx2);
                gradik3 = static_cast<ScalarType>(wx3);

                i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetFFT().osize);

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

                // compute M^{-1] = (\beta_v (\beta_w(-\ilap + 1))^{-1} + 1)^{-1}
                opik = 1.0/(beta[2]*(-lapik + 1.0));
                opik = 1.0/(beta[0]*opik + 1.0);

                // compute lap^{-1} div(b)
                this->m_x1hat[i][0] *= opik*lapinvik;
                this->m_x1hat[i][1] *= opik*lapinvik;

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
    accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_x1hat, p_x1, timer);
    accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_x2hat, p_x2, timer);
    accfft_execute_c2r_t(this->m_Opt->GetFFT().plan, this->m_x3hat, p_x3, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 3);

    ierr = this->m_WorkVecField1->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->AXPY(1.0, this->m_WorkVecField1); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




}  // namespace reg




#endif  // _OPTIMALCONTROLREGISTRATIONRELAXEDIC_CPP_
