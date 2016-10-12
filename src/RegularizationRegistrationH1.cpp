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

#ifndef _REGULARIZATIONREGISTRATIONH1_CPP_
#define _REGULARIZATIONREGISTRATIONH1_CPP_

#include "RegularizationRegistrationH1.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH1"
RegularizationRegistrationH1::RegularizationRegistrationH1() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegularizationRegistrationH1"
RegularizationRegistrationH1::~RegularizationRegistrationH1(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH1"
RegularizationRegistrationH1::RegularizationRegistrationH1(RegOpt* opt) : SuperClass(opt) {

}




/********************************************************************
 * @brief evaluates the functional
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateFunctional"
PetscErrorCode RegularizationRegistrationH1::EvaluateFunctional(ScalarType* R, VecField* v) {
    PetscErrorCode ierr;
    ScalarType  *p_v1 = NULL,*p_v2 = NULL, *p_v3 = NULL,
                *p_gv11 = NULL, *p_gv12 = NULL, *p_gv13 = NULL,
                *p_gv21 = NULL, *p_gv22 = NULL, *p_gv23 = NULL,
                *p_gv31 = NULL, *p_gv32 = NULL, *p_gv33 = NULL;
    ScalarType value, beta[2], H1v, L2v;
    std::bitset<3>XYZ = 0; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timers[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    beta[0] = this->m_Opt->GetRegNorm().beta[0];
    beta[1] = this->m_Opt->GetRegNorm().beta[1];

    *R= 0.0;

    if ( (beta[0] != 0.0)  && (beta[1] != 0.0) ) {
        if (this->m_WorkVecField==NULL){
            try{this->m_WorkVecField = new VecField(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        H1v = 0.0;

        // get arrays for velocity field
        ierr=v->GetArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);


        // X1 gradient
        ierr=this->m_WorkVecField->GetArrays(p_gv11,p_gv12,p_gv13); CHKERRQ(ierr);
        accfft_grad(p_gv11,p_gv12,p_gv13,p_v1,this->m_Opt->GetFFT().plan,&XYZ,timers);
        ierr=this->m_WorkVecField->RestoreArrays(p_gv11,p_gv12,p_gv13); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,4);

        // compute inner products
        ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&value); H1v +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&value); H1v +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&value); H1v +=value;


        // X2 gradient
        ierr=this->m_WorkVecField->GetArrays(p_gv21,p_gv22,p_gv23); CHKERRQ(ierr);
        accfft_grad(p_gv21,p_gv22,p_gv23,p_v2,this->m_Opt->GetFFT().plan,&XYZ,timers);
        ierr=this->m_WorkVecField->RestoreArrays(p_gv21,p_gv22,p_gv23); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,4);

        // compute inner products
        ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&value); H1v +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&value); H1v +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&value); H1v +=value;


        // X3 gradient
        ierr=this->m_WorkVecField->GetArrays(p_gv31,p_gv32,p_gv33); CHKERRQ(ierr);
        accfft_grad(p_gv31,p_gv32,p_gv33,p_v3,this->m_Opt->GetFFT().plan,&XYZ,timers);
        ierr=this->m_WorkVecField->RestoreArrays(p_gv31,p_gv32,p_gv33); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,4);

        // compute inner products
        ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&value); H1v +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&value); H1v +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&value); H1v +=value;

        // restore arrays for velocity field
        ierr=v->RestoreArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);

        L2v = 0.0;
        ierr=VecTDot(v->m_X1,v->m_X1,&value); L2v +=value;
        ierr=VecTDot(v->m_X2,v->m_X2,&value); L2v +=value;
        ierr=VecTDot(v->m_X3,v->m_X3,&value); L2v +=value;

        // add up contributions
        *R = 0.5*(beta[0]*H1v + beta[1]*L2v);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timers);
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief evaluates first variation of regularization norm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradient"
PetscErrorCode RegularizationRegistrationH1::EvaluateGradient(VecField* dvR, VecField* v) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta[2],scale;
    double timers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR!=NULL,"null pointer"); CHKERRQ(ierr);

    beta[0] = this->m_Opt->GetRegNorm().beta[0];
    beta[1] = this->m_Opt->GetRegNorm().beta[1];

    // if regularization weight is zero, do noting
    if ( (beta[0] == 0.0)  && (beta[1] == 0.0) ){
        ierr=dvR->SetValue(0.0); CHKERRQ(ierr);
    }
    else{

        ierr=this->Allocate(0); CHKERRQ(ierr);
        ierr=this->Allocate(1); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        // compute forward fft
        ierr=v->GetArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_v1,this->m_v1hat,timers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_v2,this->m_v2hat,timers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_v3,this->m_v3hat,timers);
        ierr=v->RestoreArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,3);

        scale = this->m_Opt->ComputeFFTScale();

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik,regop;
        IntType i,i1,i2,i3;

#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1){
            for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2){
                for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3){

                    w[0] = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                    w[1] = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                    w[2] = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                    CheckWaveNumbers(w,nx);

                    // compute bilaplacian operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = scale*(-beta[0]*lapik + beta[1]);

                    // get linear index
                    i=GetLinearIndex(i1,i2,i3,this->m_Opt->GetFFT().osize);

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
        ierr=dvR->GetArrays(p_Lv1,p_Lv2,p_Lv3); CHKERRQ(ierr);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv1hat,p_Lv1,timers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv2hat,p_Lv2,timers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv3hat,p_Lv3,timers);
        ierr=dvR->RestoreArrays(p_Lv1,p_Lv2,p_Lv3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,3);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timers);

    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies second variation of regularization norm to
 * a vector (note: the first and second variation are identical)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode RegularizationRegistrationH1::HessianMatVec(VecField* dvvR, VecField* vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType beta[2];
    PetscFunctionBegin;

    ierr=Assert(dvvR!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(vtilde!=NULL,"null pointer"); CHKERRQ(ierr);

    beta[0] = this->m_Opt->GetRegNorm().beta[0];
    beta[1] = this->m_Opt->GetRegNorm().beta[1];

    // if regularization weight is zero, do noting
    if ( (beta[0] == 0.0)  && (beta[1] == 0.0) ){
        ierr=dvvR->SetValue(0.0); CHKERRQ(ierr);
    }
    else{ ierr=this->EvaluateGradient(dvvR,vtilde); CHKERRQ(ierr); }


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvOp"
PetscErrorCode RegularizationRegistrationH1::ApplyInvOp(VecField* Ainvx, VecField* x, bool applysqrt) {
    PetscErrorCode ierr = 0;
    int nx[3];
    ScalarType *p_x1=NULL,*p_x2=NULL,*p_x3=NULL,
                *p_Ainvx1=NULL,*p_Ainvx2=NULL,*p_Ainvx3=NULL;
    ScalarType beta[2],scale;
    double timers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(Ainvx != NULL, "null pointer"); CHKERRQ(ierr);

    beta[0] = this->m_Opt->GetRegNorm().beta[0];
    beta[1] = this->m_Opt->GetRegNorm().beta[1];

    // if regularization weight is zero, do noting
    if ( (beta[0] == 0.0)  && (beta[1] == 0.0) ){

        ierr=VecCopy(x->m_X1,Ainvx->m_X1); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X2,Ainvx->m_X2); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X3,Ainvx->m_X3); CHKERRQ(ierr);

    }
    else{

        ierr=this->Allocate(0); CHKERRQ(ierr);
        ierr=this->Allocate(1); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        // compute forward fft
        ierr=x->GetArrays(p_x1,p_x2,p_x3); CHKERRQ(ierr);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x1,this->m_v1hat,timers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x2,this->m_v2hat,timers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x3,this->m_v3hat,timers);
        ierr=x->RestoreArrays(p_x1,p_x2,p_x3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,3);

        scale = this->m_Opt->ComputeFFTScale();

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik,regop;
        IntType i,i1,i2,i3;

#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1){
            for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2){
                for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3){

                    w[0] = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                    w[1] = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                    w[2] = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                    CheckWaveNumbersInv(w,nx);

                    // compute bilaplacian operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = -beta[0]*lapik + beta[1];
                    if (applysqrt) regop = sqrt(regop);
                    regop = scale/regop;

                    i=GetLinearIndex(i1,i2,i3,this->m_Opt->GetFFT().osize);

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
        ierr=Ainvx->GetArrays(p_Ainvx1,p_Ainvx2,p_Ainvx3); CHKERRQ(ierr);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv1hat,p_Ainvx1,timers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv2hat,p_Ainvx2,timers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv3hat,p_Ainvx3,timers);
        ierr=Ainvx->RestoreArrays(p_Ainvx1,p_Ainvx2,p_Ainvx3); CHKERRQ(ierr);

        this->m_Opt->IncrementCounter(FFT,3);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(timers);

    }

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief computes the largest and smallest eigenvalue of
 * the inverse regularization operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetExtremeEigValsInvOp"
PetscErrorCode RegularizationRegistrationH1::GetExtremeEigValsInvOp(ScalarType& emin, ScalarType& emax) {
    PetscErrorCode ierr=0;
    ScalarType w[3],beta1,beta2,regop;

    PetscFunctionBegin;

    beta1=this->m_Opt->GetRegNorm().beta[0];
    beta2=this->m_Opt->GetRegNorm().beta[1];

    // get max value
    w[0] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[0])/2.0;
    w[1] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[1])/2.0;
    w[2] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[2])/2.0;

    // compute largest value for operator
    regop = -(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]); // laplacian
    regop = -beta1*regop + beta2; // -beta_1 * lap + beta_2
    emin = 1.0/regop;
    emax = 1.0/beta2;  // 1.0/(\beta_1*0 + \beta_2)

    PetscFunctionReturn(ierr);
}




} // end of name space

#endif //_REGULARIZATIONREGISTRATIONH1_CPP_
