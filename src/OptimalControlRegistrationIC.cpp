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

#ifndef _OPTIMALCONTROLREGISTRATIONIC_CPP_
#define _OPTIMALCONTROLREGISTRATIONIC_CPP_

#include <math.h>
#include "OptimalControlRegistrationIC.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
OptimalControlRegistrationIC::OptimalControlRegistrationIC() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
OptimalControlRegistrationIC::~OptimalControlRegistrationIC(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
OptimalControlRegistrationIC::OptimalControlRegistrationIC(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationIC::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_x1hat = NULL;
    this->m_x2hat = NULL;
    this->m_x3hat = NULL;

    this->m_Kx1hat = NULL;
    this->m_Kx2hat = NULL;
    this->m_Kx3hat = NULL;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationIC::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_x1hat != NULL) {
        accfft_free(this->m_x1hat);
        this->m_x1hat = NULL;
    }
    if (this->m_x2hat != NULL) {
        accfft_free(this->m_x2hat);
        this->m_x2hat = NULL;
    }
    if (this->m_x3hat != NULL) {
        accfft_free(this->m_x3hat);
        this->m_x3hat = NULL;
    }
    if (this->m_Kx1hat != NULL) {
        accfft_free(this->m_Kx1hat);
        this->m_Kx1hat = NULL;
    }
    if (this->m_Kx2hat != NULL) {
        accfft_free(this->m_Kx2hat);
        this->m_Kx2hat = NULL;
    }
    if (this->m_Kx3hat != NULL) {
        accfft_free(this->m_Kx3hat);
        this->m_Kx3hat = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute the body force
 * b = K[\int_0^1 \igrad m \lambda d t],
 * where K is an operator that projects v onto the manifold of
 * divergence free velocity fields
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationIC::ComputeBodyForce() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_WorkVecField1 == NULL) {
        this->m_WorkVecField1 = new VecField(this->m_Opt);
    }
    if (this->m_WorkVecField2 == NULL) {
        this->m_WorkVecField2 = new VecField(this->m_Opt);
    }

    // assigned to work vec field 2
    ierr = SuperClass::ComputeBodyForce(); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->Copy(this->m_WorkVecField2); CHKERRQ(ierr);
    ierr = this->ApplyProjection(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->AXPY(1.0, this->m_WorkVecField1); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute the body force
 * b = K[\int_0^1\igrad\tilde{m}\lambda+\igrad m \tilde{\lambda} dt]
 * where K is an operator that projects \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationIC::ComputeIncBodyForce() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_WorkVecField1 == NULL) {
        this->m_WorkVecField1 = new VecField(this->m_Opt);
    }
    if (this->m_WorkVecField2 == NULL) {
        this->m_WorkVecField2 = new VecField(this->m_Opt);
    }

    // assigned to work vec field 2
    ierr = SuperClass::ComputeIncBodyForce(); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->Copy(this->m_WorkVecField2); CHKERRQ(ierr);
    ierr = this->ApplyProjection(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->AXPY(1.0, this->m_WorkVecField1); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationIC::SolveAdjointEquationSL() {
    PetscErrorCode ierr = 0;
    ScalarType *p_l = NULL;
    IntType nl, nt, nc, l, lnext;

    PetscFunctionBegin;

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nlocal;

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_WorkVecField1==NULL) {
        this->m_WorkVecField1 = new VecField(this->m_Opt);
    }

    // scale v by -1
    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->Scale(-1.0); CHKERRQ(ierr);

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1, "adjoint"); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            l = (nt-j)*nc*nl + k*nl;
            lnext = (nt-(j+1))*nc*nl + k*nl;

            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_l+lnext, p_l+l, "adjoint"); CHKERRQ(ierr);
        }  // for all image components
    }  // for all time points
    ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);

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
PetscErrorCode OptimalControlRegistrationIC::SolveIncAdjointEquationGNSL(void) {
    PetscErrorCode ierr = 0;
    IntType nl, nt, nc, l, lnext;
    ScalarType *p_ltilde = NULL;

    PetscFunctionBegin;

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nlocal;

    // for all time points
    ierr = VecGetArray(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            l = (nt-j)*nc*nl + k*nl;
            lnext = (nt-(j+1))*nc*nl + k*nl;
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_ltilde+lnext, p_ltilde+l, "adjoint"); CHKERRQ(ierr);
        }  // for all image components
    }  // for all time points
    ierr = VecRestoreArray(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply projection to map \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationIC::ApplyProjection(VecField* x) {
    PetscErrorCode ierr = 0;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL, scale;
    long int nx[3];
    IntType nalloc;
    double timer[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    nx[0] = static_cast<long int>(this->m_Opt->GetNumGridPoints(0));
    nx[1] = static_cast<long int>(this->m_Opt->GetNumGridPoints(1));
    nx[2] = static_cast<long int>(this->m_Opt->GetNumGridPoints(2));

    nalloc = this->m_Opt->GetFFT().nalloc;
    scale = this->m_Opt->ComputeFFTScale();

    if (this->m_x1hat == NULL) {
        this->m_x1hat = reinterpret_cast<FFTScaType*>(accfft_alloc(nalloc));
    }
    if (this->m_x2hat == NULL) {
        this->m_x2hat = reinterpret_cast<FFTScaType*>(accfft_alloc(nalloc));
    }
    if (this->m_x3hat == NULL) {
        this->m_x3hat = reinterpret_cast<FFTScaType*>(accfft_alloc(nalloc));
    }

    if (this->m_Kx1hat == NULL) {
        this->m_Kx1hat = reinterpret_cast<FFTScaType*>(accfft_alloc(nalloc));
    }
    if (this->m_Kx2hat == NULL) {
        this->m_Kx2hat = reinterpret_cast<FFTScaType*>(accfft_alloc(nalloc));
    }
    if (this->m_Kx3hat == NULL) {
        this->m_Kx3hat = reinterpret_cast<FFTScaType*>(accfft_alloc(nalloc));
    }

    ierr = x->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    // compute forward fft
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_x1, this->m_x1hat, timer);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_x2, this->m_x2hat, timer);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan, p_x3, this->m_x3hat, timer);
    this->m_Opt->IncrementCounter(FFT,3);

#pragma omp parallel
{
    long int x1, x2, x3, wx1, wx2, wx3;
    ScalarType lapinvik, gradik1, gradik2, gradik3;
    IntType i, i1, i2, i3;
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

                i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetFFT().osize);

                // compute div(b)
                this->m_Kx1hat[i][0] = -scale*(gradik1*this->m_x1hat[i][0]
                                             + gradik2*this->m_x2hat[i][0]
                                             + gradik3*this->m_x3hat[i][0]);

                this->m_Kx1hat[i][1] =  scale*(gradik1*this->m_x1hat[i][1]
                                             + gradik2*this->m_x2hat[i][1]
                                             + gradik3*this->m_x3hat[i][1]);

                // compute lap^{-1} div(b)
                this->m_Kx1hat[i][0] *= lapinvik;
                this->m_Kx1hat[i][1] *= lapinvik;

                // compute x2 gradient of lab^{-1} div(b)
                this->m_Kx2hat[i][0] = -gradik2*this->m_Kx1hat[i][0];
                this->m_Kx2hat[i][1] =  gradik2*this->m_Kx1hat[i][1];

                // compute x3 gradient of lab^{-1} div(b)
                this->m_Kx3hat[i][0] = -gradik3*this->m_Kx1hat[i][0];
                this->m_Kx3hat[i][1] =  gradik3*this->m_Kx1hat[i][1];

                // compute x1 gradient of lab^{-1} div(b)
                this->m_Kx1hat[i][0] *= -gradik1;
                this->m_Kx1hat[i][1] *=  gradik1;
            }
        }
    }
}  // pragma omp parallel

    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Kx1hat, p_x1, timer);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Kx2hat, p_x2, timer);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan, this->m_Kx3hat, p_x3, timer);
    this->m_Opt->IncrementCounter(FFT,3);

    ierr = x->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _OPTIMALCONTROLREGISTRATIONIC_CPP_
