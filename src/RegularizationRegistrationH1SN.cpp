#ifndef _REGULARIZATIONREGISTRATIONH1SN_CPP_
#define _REGULARIZATIONREGISTRATIONH1SN_CPP_

#include "RegularizationRegistrationH1SN.hpp"




namespace reg
{




/********************************************************************
 * Name: RegularizationRegistrationH1SN
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH1SN"
RegularizationRegistrationH1SN::RegularizationRegistrationH1SN() : SuperClass()
{

}




/********************************************************************
 * Name: ~RegularizationRegistrationH1SN
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegularizationRegistrationH1SN"
RegularizationRegistrationH1SN::~RegularizationRegistrationH1SN(void)
{
    this->ClearMemory();
}




/********************************************************************
 * Name: RegularizationRegistrationH1SN
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH1SN"
RegularizationRegistrationH1SN::RegularizationRegistrationH1SN(RegOpt* opt) : SuperClass(opt)
{

}




/********************************************************************
 * Name: EvaluateFunctional
 * Description: evaluates the functional
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateFunctional"
PetscErrorCode RegularizationRegistrationH1SN::EvaluateFunctional(ScalarType* R, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType  *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_gv11=NULL,*p_gv12=NULL,*p_gv13=NULL,
                *p_gv21=NULL,*p_gv22=NULL,*p_gv23=NULL,
                *p_gv31=NULL,*p_gv32=NULL,*p_gv33=NULL;
    double ffttimers[5]={0,0,0,0,0};
    std::bitset<3>XYZ={111};
    ScalarType beta,hd,value;
    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    // compute hd
    hd = this->m_Opt->GetLebesqueMeasure();

    // get regularization parameter
    beta = this->m_Opt->GetRegularizationWeight();

    *R = 0.0;

    // if regularization weight is zero, do noting
    if (beta != 0.0){

        if (this->m_WorkVecField==NULL){
            try{this->m_WorkVecField = new VecField(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        // get arrays
        ierr=VecGetArray(v->m_X1,&p_v1); CHKERRQ(ierr);
        ierr=VecGetArray(v->m_X2,&p_v2); CHKERRQ(ierr);
        ierr=VecGetArray(v->m_X3,&p_v3); CHKERRQ(ierr);


        // X1
        // compute gradient of first component of vector field
        ierr=VecGetArray(this->m_WorkVecField->m_X1,&p_gv11); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField->m_X2,&p_gv12); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField->m_X3,&p_gv13); CHKERRQ(ierr);

        // compute gradient
        accfft_grad(p_gv11,p_gv12,p_gv13,p_v1,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

        ierr=VecRestoreArray(this->m_WorkVecField->m_X1,&p_gv11); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField->m_X2,&p_gv12); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField->m_X3,&p_gv13); CHKERRQ(ierr);

        // compute inner products
        ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&value); *R +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&value); *R +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&value); *R +=value;


        // X2
        // compute gradient of second component of vector field
        ierr=VecGetArray(this->m_WorkVecField->m_X1,&p_gv21); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField->m_X2,&p_gv22); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField->m_X3,&p_gv23); CHKERRQ(ierr);

        // compute gradient
        accfft_grad(p_gv21,p_gv22,p_gv23,p_v2,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

        ierr=VecRestoreArray(this->m_WorkVecField->m_X1,&p_gv21); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField->m_X2,&p_gv22); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField->m_X3,&p_gv23); CHKERRQ(ierr);

        // compute inner products
        ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&value); *R +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&value); *R +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&value); *R +=value;


        // X3
        // compute gradient of second component of vector field
        ierr=VecGetArray(this->m_WorkVecField->m_X1,&p_gv31); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField->m_X2,&p_gv32); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField->m_X3,&p_gv33); CHKERRQ(ierr);

        // compute gradient
        accfft_grad(p_gv31,p_gv32,p_gv33,p_v3,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

        ierr=VecRestoreArray(this->m_WorkVecField->m_X1,&p_gv31); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField->m_X2,&p_gv32); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField->m_X3,&p_gv33); CHKERRQ(ierr);

        // compute inner products
        ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&value); *R +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&value); *R +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&value); *R +=value;

        ierr=VecRestoreArray(v->m_X1,&p_v1); CHKERRQ(ierr);
        ierr=VecRestoreArray(v->m_X2,&p_v2); CHKERRQ(ierr);
        ierr=VecRestoreArray(v->m_X3,&p_v3); CHKERRQ(ierr);

        // multiply with regularization weight
        *R = 0.5*hd*beta*(*R);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(ffttimers);
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateGradient
 * Description: evaluates first variation of regularization norm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradient"
PetscErrorCode RegularizationRegistrationH1SN::EvaluateGradient(VecField* dvR, VecField* v)
{
    PetscErrorCode ierr;
    unsigned int n[3],iosize[3];
    int isize[3],osize[3],istart[3],ostart[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta,scale,hd;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR != NULL,"null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegularizationWeight();

    // if regularization weight is zero, do noting
    if (beta == 0.0){
        ierr=VecSet(dvR->m_X1,0.0); CHKERRQ(ierr);
        ierr=VecSet(dvR->m_X2,0.0); CHKERRQ(ierr);
        ierr=VecSet(dvR->m_X3,0.0); CHKERRQ(ierr);
    }
    else{

        ierr=this->Allocate(); CHKERRQ(ierr);

        accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                                isize,istart,osize,ostart,
                                                this->m_Opt->m_MiscOpt->c_comm);

        for (unsigned int i=0; i < 3; ++i){
            n[i]      = this->m_Opt->m_MiscOpt->N[i];
            iosize[i] = osize[i];
        }
        scale = this->m_Opt->ComputeFFTScale();

        ierr=VecGetArray(v->m_X1,&p_v1); CHKERRQ(ierr);
        ierr=VecGetArray(v->m_X2,&p_v2); CHKERRQ(ierr);
        ierr=VecGetArray(v->m_X3,&p_v3); CHKERRQ(ierr);

        // compute forward fft
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_v1,this->m_v1hat,ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_v2,this->m_v2hat,ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_v3,this->m_v3hat,ffttimers);
        this->m_Opt->IncrementCounter(FFT,3);

        ierr=VecRestoreArray(v->m_X1,&p_v1); CHKERRQ(ierr);
        ierr=VecRestoreArray(v->m_X2,&p_v2); CHKERRQ(ierr);
        ierr=VecRestoreArray(v->m_X3,&p_v3); CHKERRQ(ierr);

        beta = this->m_Opt->GetRegularizationWeight();

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik,regop;
        IntType i;

#pragma omp for
        for (unsigned int i1 = 0; i1 < iosize[0]; ++i1){
            for (unsigned int i2 = 0; i2 < iosize[1]; i2++){
                for (unsigned int i3 = 0; i3 < iosize[2]; ++i3){

                    w[0] = static_cast<long int>(i1 + ostart[0]);
                    w[1] = static_cast<long int>(i2 + ostart[1]);
                    w[2] = static_cast<long int>(i3 + ostart[2]);

                    CheckWaveNumbers(w,n);

                    // compute bilaplacian operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = -scale*beta*lapik;

                    i=GetLinearIndex(i1,i2,i3,iosize);

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

        ierr=VecGetArray(dvR->m_X1,&p_Lv1); CHKERRQ(ierr);
        ierr=VecGetArray(dvR->m_X2,&p_Lv2); CHKERRQ(ierr);
        ierr=VecGetArray(dvR->m_X3,&p_Lv3); CHKERRQ(ierr);

        // compute inverse fft
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv1hat,p_Lv1,ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv2hat,p_Lv2,ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv3hat,p_Lv3,ffttimers);
        this->m_Opt->IncrementCounter(FFT,3);

        // restore arrays
        ierr=VecRestoreArray(dvR->m_X1,&p_Lv1); CHKERRQ(ierr);
        ierr=VecRestoreArray(dvR->m_X2,&p_Lv2); CHKERRQ(ierr);
        ierr=VecRestoreArray(dvR->m_X3,&p_Lv3); CHKERRQ(ierr);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(ffttimers);

        // compute hd
        hd = this->m_Opt->GetLebesqueMeasure();

        // scale vector field
        ierr=dvR->Scale(hd); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianMatvec
 * Description: applies second variation of regularization norm to
 * a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode RegularizationRegistrationH1SN::HessianMatVec(VecField* dvvR, VecField* vtilde)
{
    PetscErrorCode ierr;
    ScalarType beta,hd;
    PetscFunctionBegin;

    ierr=Assert(vtilde != NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvvR != NULL,"null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegularizationWeight();

    // if regularization weight is zero, do noting
    if (beta == 0.0){
        ierr=VecSet(dvvR->m_X1,0.0); CHKERRQ(ierr);
        ierr=VecSet(dvvR->m_X2,0.0); CHKERRQ(ierr);
        ierr=VecSet(dvvR->m_X3,0.0); CHKERRQ(ierr);
    }
    else{ ierr=this->EvaluateGradient(dvvR,vtilde); CHKERRQ(ierr); }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ApplyInverseOperator
 * Description: apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInverseOperator"
PetscErrorCode RegularizationRegistrationH1SN::ApplyInverseOperator(VecField* Ainvx, VecField* x)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_x1=NULL,*p_x2=NULL,*p_x3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta,scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(Ainvx != NULL,"null pointer"); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegularizationWeight();

    // if regularization weight is zero, do noting
    if (beta == 0.0){
        ierr=VecCopy(x->m_X1,Ainvx->m_X1); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X2,Ainvx->m_X2); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X3,Ainvx->m_X3); CHKERRQ(ierr);
    }
    else{

        ierr=this->Allocate(); CHKERRQ(ierr);

        accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                                isize,istart,osize,ostart,
                                                this->m_Opt->m_MiscOpt->c_comm);

        for (unsigned int i=0; i < 3; ++i){
            n[i] = this->m_Opt->m_MiscOpt->N[i];
            iosize[i] = osize[i];
        }
        scale = this->m_Opt->ComputeFFTScale();

        ierr=VecGetArray(x->m_X1,&p_x1); CHKERRQ(ierr);
        ierr=VecGetArray(x->m_X2,&p_x2); CHKERRQ(ierr);
        ierr=VecGetArray(x->m_X3,&p_x3); CHKERRQ(ierr);

        // compute forward fft
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_x1,this->m_v1hat,ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_x2,this->m_v2hat,ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_x3,this->m_v3hat,ffttimers);
        this->m_Opt->IncrementCounter(FFT,3);

        ierr=VecRestoreArray(x->m_X1,&p_x1); CHKERRQ(ierr);
        ierr=VecRestoreArray(x->m_X2,&p_x2); CHKERRQ(ierr);
        ierr=VecRestoreArray(x->m_X3,&p_x3); CHKERRQ(ierr);

        beta = this->m_Opt->GetRegularizationWeight();

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik,regop;
        IntType i;

#pragma omp for
        for (unsigned int i1 = 0; i1 < iosize[0]; ++i1){
            for (unsigned int i2 = 0; i2 < iosize[1]; i2++){
                for (unsigned int i3 = 0; i3 < iosize[2]; ++i3){

                    w[0] = static_cast<long int>(i1 + ostart[0]);
                    w[1] = static_cast<long int>(i2 + ostart[1]);
                    w[2] = static_cast<long int>(i3 + ostart[2]);

                    CheckWaveNumbersInv(w,n);

                    // compute bilaplacian operator
                    lapik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                    // compute regularization operator
                    regop = (abs(lapik) == 0.0) ? scale/beta : scale/(-beta*lapik);

                    i=GetLinearIndex(i1,i2,i3,iosize);

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

        ierr=VecGetArray(Ainvx->m_X1,&p_Lv1); CHKERRQ(ierr);
        ierr=VecGetArray(Ainvx->m_X2,&p_Lv2); CHKERRQ(ierr);
        ierr=VecGetArray(Ainvx->m_X3,&p_Lv3); CHKERRQ(ierr);

        // compute inverse fft
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv1hat,p_Lv1,ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv2hat,p_Lv2,ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv3hat,p_Lv3,ffttimers);
        this->m_Opt->IncrementCounter(FFT,3);

        // restore arrays
        ierr=VecRestoreArray(Ainvx->m_X1,&p_Lv1); CHKERRQ(ierr);
        ierr=VecRestoreArray(Ainvx->m_X2,&p_Lv2); CHKERRQ(ierr);
        ierr=VecRestoreArray(Ainvx->m_X3,&p_Lv3); CHKERRQ(ierr);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(ffttimers);

    }

    PetscFunctionReturn(0);
}




} // end of name space




#endif //_REGULARIZATIONREGISTRATIONH1SN_CPP_