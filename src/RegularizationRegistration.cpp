#ifndef _REGULARIZATIONREGISTRATION_CPP_
#define _REGULARIZATIONREGISTRATION_CPP_

#include "RegularizationRegistration.hpp"




namespace reg
{




/********************************************************************
 * Name: RegularizationRegistration
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistration"
RegularizationRegistration::RegularizationRegistration()
{
    this->Initialize();
}




/********************************************************************
 * Name: ~RegularizationRegistration
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegularizationRegistration"
RegularizationRegistration::~RegularizationRegistration(void)
{
    this->ClearMemory();
}




/********************************************************************
 * Name: RegularizationRegistration
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistration"
RegularizationRegistration::RegularizationRegistration(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * Name: Initialize
 * Description: init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode RegularizationRegistration::Initialize(void)
{
    PetscFunctionBegin;

    this->m_Opt = NULL;
    this->m_WorkVecField = NULL;

    this->m_v1hat = NULL;
    this->m_v2hat = NULL;
    this->m_v3hat = NULL;

    this->m_Lv1hat = NULL;
    this->m_Lv2hat = NULL;
    this->m_Lv3hat = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode RegularizationRegistration::ClearMemory(void)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_WorkVecField!=NULL){
        delete this->m_WorkVecField;
        this->m_WorkVecField = NULL;
    }

    // clear the fft arrays
    ierr=this->Deallocate(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Allocate
 * Description: allocate arrays for fft (we might have to do this
 * in several functions, so we do it collectively here)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode RegularizationRegistration::Allocate(void)
{
    int isize[3],osize[3],istart[3],ostart[3];
    size_t alloc_max;

    PetscFunctionBegin;

    // get local pencil size and allocation size
    alloc_max=accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                                        isize,istart,osize,ostart,
                                                        this->m_Opt->m_MiscOpt->c_comm);

    if(this->m_v1hat == NULL){
        this->m_v1hat=(FFTScaType*)accfft_alloc(alloc_max);
    }
    if(this->m_v2hat == NULL){
        this->m_v2hat=(FFTScaType*)accfft_alloc(alloc_max);
    }
    if(this->m_v3hat == NULL){
        this->m_v3hat=(FFTScaType*)accfft_alloc(alloc_max);
    }

    if(this->m_Lv1hat == NULL){
        this->m_Lv1hat=(FFTScaType*)accfft_alloc(alloc_max);
    }
    if(this->m_Lv2hat == NULL){
        this->m_Lv2hat=(FFTScaType*)accfft_alloc(alloc_max);
    }
    if(this->m_Lv3hat == NULL){
        this->m_Lv3hat=(FFTScaType*)accfft_alloc(alloc_max);
    }

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: Deallocate
 * Description: deallocate arrays for fft (we might have to do this
 * in several functions, so we do it collectively here)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Deallocate"
PetscErrorCode RegularizationRegistration::Deallocate(void)
{
    PetscFunctionBegin;

    if(this->m_v1hat!=NULL){
        accfft_free(this->m_v1hat);
        this->m_v1hat = NULL;
    }
    if(this->m_v2hat!=NULL){
        accfft_free(this->m_v2hat);
        this->m_v2hat = NULL;
    }
    if(this->m_v3hat!=NULL){
        accfft_free(this->m_v3hat);
        this->m_v3hat = NULL;
    }

    if(this->m_Lv1hat!=NULL){
        accfft_free(this->m_Lv1hat);
        this->m_Lv1hat = NULL;
    }
    if(this->m_Lv2hat!=NULL){
        accfft_free(this->m_Lv2hat);
        this->m_Lv2hat = NULL;
    }
    if(this->m_Lv3hat!=NULL){
        accfft_free(this->m_Lv3hat);
        this->m_Lv3hat = NULL;
    }

    PetscFunctionReturn(0);
}





/********************************************************************
 * Name: EvaluateFunctional
 * Description: evaluates the functional
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateFunctional"
PetscErrorCode RegularizationRegistration::EvaluateFunctional(ScalarType* R, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType beta,hd;
    PetscFunctionBegin;

    beta = this->m_Opt->GetRegularizationWeight();

    // if regularization weight is zero, do noting
    if (beta == 0.0){ *R = 0.0; }
    else{

        switch (this->m_Opt->GetRegNorm()){
            case L2: { ierr=this->EvaluateL2(R,v); CHKERRQ(ierr); break; }
            case H1: { ierr=this->EvaluateH1(R,v); CHKERRQ(ierr); break; }
            case H2: { ierr=this->EvaluateH2(R,v); CHKERRQ(ierr); break; }
            case H1SN: { ierr=this->EvaluateH1S(R,v); CHKERRQ(ierr); break; }
            case H2SN: { ierr=this->EvaluateH2S(R,v); CHKERRQ(ierr); break; }
            default:
            {
                ierr=ThrowError("regularization norm not defined"); CHKERRQ(ierr);
                break;
            }
        }
    }

    // compute hd
    hd = this->m_Opt->GetLebesqueMeasure();

    // multiply with regularization weight
    *R = 0.5*hd*(*R);

    PetscFunctionReturn(0);
}


/********************************************************************
 * Name: EvaluateL2
 * Description: evaluates the L2 norm for a given velocity field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateL2"
PetscErrorCode RegularizationRegistration::EvaluateL2(ScalarType* R, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType value,beta;

    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    *R = 0.0;

    // compute inner products
    ierr=VecTDot(v->m_X1,v->m_X1,&value); *R +=value;
    ierr=VecTDot(v->m_X2,v->m_X2,&value); *R +=value;
    ierr=VecTDot(v->m_X3,v->m_X3,&value); *R +=value;

    beta = this->m_Opt->GetRegularizationWeight();
    *R *= beta;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateH1
 * Description: evaluates the H1 norm for a given velocity field
 * \int_{\Omega} (\igrad + I) \vect{v} \cdot (\igrad + I) \vect{v} dx
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateH1"
PetscErrorCode RegularizationRegistration::EvaluateH1(ScalarType* R, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType value,beta,H1v,L2v;

    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->EvaluateH1S(&H1v,v); CHKERRQ(ierr);

    L2v = 0.0;
    beta = this->m_Opt->GetRegularizationWeight();

    ierr=VecTDot(v->m_X1,v->m_X1,&value); L2v +=value;
    ierr=VecTDot(v->m_X2,v->m_X2,&value); L2v +=value;
    ierr=VecTDot(v->m_X3,v->m_X3,&value); L2v +=value;

    L2v *= beta;

    *R = H1v + L2v;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateH1SN
 * Description: evaluates the H1 seminorm for a given velocity field
 * \int_{\Omega} \igrad \vect{v} \cdot \igrad \vect{v} dx
 *
 * to not allocate to much memory, we'll compute the gradients of
 * each component sepeartely and inbetween evaluate the inner
 * product
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateH1S"
PetscErrorCode RegularizationRegistration::EvaluateH1S(ScalarType* R, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType  *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_gv11=NULL,*p_gv12=NULL,*p_gv13=NULL,
                *p_gv21=NULL,*p_gv22=NULL,*p_gv23=NULL,
                *p_gv31=NULL,*p_gv32=NULL,*p_gv33=NULL,
                value,beta;

    double ffttimers[5]={0,0,0,0,0};
    std::bitset<3>XYZ={111};
    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    if (this->m_WorkVecField==NULL){
        try{this->m_WorkVecField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    *R = 0.0;

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

    beta = this->m_Opt->GetRegularizationWeight();
    *R *= beta;

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateH2
 * Description: evaluates the H2 norm for a given velocity field
 * \int_{\Omega} (\ilap + I) \vect{v} \cdot (\ilap + I)\vect{v} dx
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateH2"
PetscErrorCode RegularizationRegistration::EvaluateH2(ScalarType* R, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType sqrtbeta[2],ipxi,scale;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    double ffttimers[5]={0,0,0,0,0};
    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    if (this->m_WorkVecField==NULL){
        try{this->m_WorkVecField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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

    sqrtbeta[0] = sqrt(this->m_Opt->GetRegularizationWeight());
    sqrtbeta[1] = sqrt(this->m_Opt->GetRegularizationWeight());

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

                i=GetLinearIndex(i1,i2,i3,iosize);

                // compute regularization operator
                regop = scale*(sqrtbeta[0]*lapik + sqrtbeta[1]);

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

    ierr=VecGetArray(this->m_WorkVecField->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField->m_X3,&p_Lv3); CHKERRQ(ierr);

    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv1hat,p_Lv1,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv2hat,p_Lv2,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv3hat,p_Lv3,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

    // restore arrays
    ierr=VecRestoreArray(this->m_WorkVecField->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField->m_X3,&p_Lv3); CHKERRQ(ierr);

    *R=0.0;
    ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&ipxi); CHKERRQ(ierr); *R += ipxi;
    ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&ipxi); CHKERRQ(ierr); *R += ipxi;
    ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&ipxi); CHKERRQ(ierr); *R += ipxi;


    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}






/********************************************************************
 * Name: EvaluateH2S
 * Description: evaluates the H2 seminorm for a given velocity field
 * \int_{\Omega} \ilap \vect{v} \cdot \ilap \vect{v} dx
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateH2S"
PetscErrorCode RegularizationRegistration::EvaluateH2S(ScalarType* R, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType sqrtbeta,ipxi,scale;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    double ffttimers[5]={0,0,0,0,0};
    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    if (this->m_WorkVecField==NULL){
        try{this->m_WorkVecField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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

    sqrtbeta = sqrt(this->m_Opt->GetRegularizationWeight());

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

                i=GetLinearIndex(i1,i2,i3,iosize);

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

    ierr=VecGetArray(this->m_WorkVecField->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField->m_X3,&p_Lv3); CHKERRQ(ierr);

    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv1hat,p_Lv1,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv2hat,p_Lv2,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv3hat,p_Lv3,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

    // restore arrays
    ierr=VecRestoreArray(this->m_WorkVecField->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField->m_X3,&p_Lv3); CHKERRQ(ierr);

    *R=0.0;
    ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&ipxi); CHKERRQ(ierr); *R += ipxi;
    ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&ipxi); CHKERRQ(ierr); *R += ipxi;
    ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&ipxi); CHKERRQ(ierr); *R += ipxi;

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateGradient
 * Description: evaluates first variation of regularization norm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradient"
PetscErrorCode RegularizationRegistration::EvaluateGradient(VecField* dvR, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType beta,hd;
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
        switch (this->m_Opt->GetRegNorm()){
            case L2: { ierr=this->EvaluateGradL2(dvR,v); CHKERRQ(ierr); break; }
            case H1: { ierr=this->EvaluateGradH1(dvR,v); CHKERRQ(ierr); break; }
            case H2: { ierr=this->EvaluateGradH2(dvR,v); CHKERRQ(ierr); break; }
            case H1SN: { ierr=this->EvaluateGradH1S(dvR,v); CHKERRQ(ierr); break; }
            case H2SN: { ierr=this->EvaluateGradH2S(dvR,v); CHKERRQ(ierr); break; }
            default:
            {
                ierr=ThrowError("regularization norm not defined"); CHKERRQ(ierr);
                break;
            }
        }

        // compute hd
        hd = this->m_Opt->GetLebesqueMeasure();

        // scale vector field
        ierr=dvR->Scale(hd); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateGradL2
 * Description: evaluates first variation of L2 norm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradL2"
PetscErrorCode RegularizationRegistration::EvaluateGradL2(VecField* dvR, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType beta;
    PetscFunctionBegin;

    ierr=dvR->Copy(v); CHKERRQ(ierr);

    beta = this->m_Opt->GetRegularizationWeight();
    ierr=dvR->Scale(beta); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateGradH1
 * Description: evaluates first variation of H1 seminorm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradH1"
PetscErrorCode RegularizationRegistration::EvaluateGradH1(VecField* dvR, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta[2],scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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

    beta[0] = this->m_Opt->GetRegularizationWeight();
    beta[1] = this->m_Opt->GetRegularizationWeight();

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
                regop = scale*(-beta[0]*lapik + beta[1]);

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

    PetscFunctionReturn(0);
}






/********************************************************************
 * Name: EvaluateGradH1S
 * Description: evaluates first variation of H1 seminorm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradH1S"
PetscErrorCode RegularizationRegistration::EvaluateGradH1S(VecField* dvR, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta,scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateGradH2
 * Description: evaluates first variation of H2 norm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradH2"
PetscErrorCode RegularizationRegistration::EvaluateGradH2(VecField* dvR, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL,beta[2],scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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

    beta[0] = this->m_Opt->GetRegularizationWeight();
    beta[1] = this->m_Opt->GetRegularizationWeight();

#pragma omp parallel
{
    long int w[3];
    ScalarType biharik,regop;
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
                biharik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                biharik *= biharik;

                // compute regularization operator
                regop = scale*(beta[0]*biharik + beta[1]);

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

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateGradH2S
 * Description: evaluates first variation of H2 seminorm
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradH2S"
PetscErrorCode RegularizationRegistration::EvaluateGradH2S(VecField* dvR, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta,scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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
    ScalarType biharik,regop;
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
                biharik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                biharik *= biharik;

                i=GetLinearIndex(i1,i2,i3,iosize);

                // compute regularization operator
                regop = scale*beta*biharik;

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

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianMatvec
 * Description: applies second variation of regularization norm to
 * a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode RegularizationRegistration::HessianMatVec(VecField* dvvR, VecField* vtilde)
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
    else{
        switch (this->m_Opt->GetRegNorm()){
            case L2: { ierr=this->HessianMatVecL2(dvvR,vtilde); CHKERRQ(ierr); break; }
            case H1: { ierr=this->HessianMatVecH1(dvvR,vtilde); CHKERRQ(ierr); break; }
            case H2: { ierr=this->HessianMatVecH2(dvvR,vtilde); CHKERRQ(ierr); break; }
            case H1SN: { ierr=this->HessianMatVecH1S(dvvR,vtilde); CHKERRQ(ierr); break; }
            case H2SN: { ierr=this->HessianMatVecH2S(dvvR,vtilde); CHKERRQ(ierr); break; }
            default:
            {
                ierr=ThrowError("regularization norm not defined"); CHKERRQ(ierr);
                break;
            }
        }
    }

    // scale with hd
    hd = this->m_Opt->GetLebesqueMeasure();
    ierr=dvvR->Scale(hd); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianMatvecL2
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVecL2"
PetscErrorCode RegularizationRegistration::HessianMatVecL2(VecField* dvvR, VecField* vtilde)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // first and second variation are identical
    ierr=this->EvaluateGradL2(dvvR,vtilde); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianMatvecH1
 * Description: applies the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVecH1"
PetscErrorCode RegularizationRegistration::HessianMatVecH1(VecField* dvvR, VecField* vtilde)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // first and second variation are identical
    ierr=this->EvaluateGradH1(dvvR,vtilde); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianMatvecH1S
 * Description: applies the hessian to a vector (laplacian operator)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVecH1S"
PetscErrorCode RegularizationRegistration::HessianMatVecH1S(VecField* dvvR, VecField* vtilde)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // first and second variation are identical
    ierr=this->EvaluateGradH1S(dvvR,vtilde); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianMatvecH2
 * Description: applies the hessian to a vector (biharmonic operator)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVecH2"
PetscErrorCode RegularizationRegistration::HessianMatVecH2(VecField* dvvR, VecField* vtilde)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // first and second variation are identical
    ierr=this->EvaluateGradH2(dvvR,vtilde); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianMatvecH2S
 * Description: applies the hessian to a vector (biharmonic operator)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVecH2S"
PetscErrorCode RegularizationRegistration::HessianMatVecH2S(VecField* dvvR, VecField* vtilde)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // first and second variation are identical
    ierr=this->EvaluateGradH2S(dvvR,vtilde); CHKERRQ(ierr);

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
PetscErrorCode RegularizationRegistration::ApplyInverseOperator(VecField* Ainvx, VecField* x)
{
    PetscErrorCode ierr;
    ScalarType beta;
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
        switch (this->m_Opt->GetRegNorm()){
            case L2: { ierr=this->ApplyInvOpL2(Ainvx,x); CHKERRQ(ierr); break; }
            case H1: { ierr=this->ApplyInvOpH1(Ainvx,x); CHKERRQ(ierr); break; }
            case H2: { ierr=this->ApplyInvOpH2(Ainvx,x); CHKERRQ(ierr); break; }
            case H1SN: { ierr=this->ApplyInvOpH1S(Ainvx,x); CHKERRQ(ierr); break; }
            case H2SN: { ierr=this->ApplyInvOpH2S(Ainvx,x); CHKERRQ(ierr); break; }
            default:
            {
                ierr=ThrowError("regularization norm not defined"); CHKERRQ(ierr);
                break;
            }
        }
    }

    PetscFunctionReturn(0);
}


/********************************************************************
 * Name: ApplyInvOpL2
 * Description: apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvOpL2"
PetscErrorCode RegularizationRegistration::ApplyInvOpL2(VecField* dvRinvx, VecField* x)
{
    PetscErrorCode ierr;
    ScalarType beta;
    PetscFunctionBegin;

    beta = this->m_Opt->GetRegularizationWeight();

    ierr=dvRinvx->Copy(x); CHKERRQ(ierr);
    ierr=dvRinvx->Scale(1.0/beta); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ApplyInvOpH1
 * Description: apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvOpH1"
PetscErrorCode RegularizationRegistration::ApplyInvOpH1(VecField* dvRinvx, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta[2],scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvRinvx != NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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

    beta[0] = this->m_Opt->GetRegularizationWeight();
    beta[1] = this->m_Opt->GetRegularizationWeight();

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
                regop = scale/(-beta[0]*lapik + beta[1]);

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

    ierr=VecGetArray(dvRinvx->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecGetArray(dvRinvx->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecGetArray(dvRinvx->m_X3,&p_Lv3); CHKERRQ(ierr);

    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv1hat,p_Lv1,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv2hat,p_Lv2,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv3hat,p_Lv3,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

    // restore arrays
    ierr=VecRestoreArray(dvRinvx->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecRestoreArray(dvRinvx->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecRestoreArray(dvRinvx->m_X3,&p_Lv3); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ApplyInvOpH1S
 * Description: apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvOpH1S"
PetscErrorCode RegularizationRegistration::ApplyInvOpH1S(VecField* dvRinvx, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta,scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvRinvx != NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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

    ierr=VecGetArray(dvRinvx->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecGetArray(dvRinvx->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecGetArray(dvRinvx->m_X3,&p_Lv3); CHKERRQ(ierr);

    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv1hat,p_Lv1,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv2hat,p_Lv2,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv3hat,p_Lv3,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

    // restore arrays
    ierr=VecRestoreArray(dvRinvx->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecRestoreArray(dvRinvx->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecRestoreArray(dvRinvx->m_X3,&p_Lv3); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ApplyInvOpH2
 * Description: apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvOpH2"
PetscErrorCode RegularizationRegistration::ApplyInvOpH2(VecField* dvRinvx, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta[2],scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvRinvx != NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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

    beta[0] = this->m_Opt->GetRegularizationWeight();
    beta[1] = this->m_Opt->GetRegularizationWeight();

#pragma omp parallel
{
    long int w[3];
    ScalarType biharik,regop;
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
                biharik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                biharik *= biharik;

                // compute regularization operator
                regop = scale/(beta[0]*biharik + beta[1]);

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

    ierr=VecGetArray(dvRinvx->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecGetArray(dvRinvx->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecGetArray(dvRinvx->m_X3,&p_Lv3); CHKERRQ(ierr);

    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv1hat,p_Lv1,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv2hat,p_Lv2,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv3hat,p_Lv3,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

    // restore arrays
    ierr=VecRestoreArray(dvRinvx->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecRestoreArray(dvRinvx->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecRestoreArray(dvRinvx->m_X3,&p_Lv3); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ApplyInvOpH2S
 * Description: apply the inverse of the regularization operator; we
 * can invert this operator analytically due to the spectral
 * discretization
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvOpH2S"
PetscErrorCode RegularizationRegistration::ApplyInvOpH2S(VecField* dvRinvx, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3];
    unsigned int n[3],iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta,scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvRinvx != NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Allocate(); CHKERRQ(ierr);

    accfft_local_size_dft_r2c_t<ScalarType>(this->m_Opt->m_MiscOpt->N,
                                            isize,istart,osize,ostart,
                                            this->m_Opt->m_MiscOpt->c_comm);

    for (unsigned int i=0; i < 3; ++i){
        n[i] = this->m_Opt->m_MiscOpt->N[i];
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
    ScalarType biharik,regop;
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
                biharik = -static_cast<ScalarType>(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                biharik *= biharik;

                // compute regularization operator
                regop = (abs(biharik) == 0.0) ? scale : scale/(beta*biharik);

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

    ierr=VecGetArray(dvRinvx->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecGetArray(dvRinvx->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecGetArray(dvRinvx->m_X3,&p_Lv3); CHKERRQ(ierr);

    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv1hat,p_Lv1,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv2hat,p_Lv2,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Lv3hat,p_Lv3,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

    // restore arrays
    ierr=VecRestoreArray(dvRinvx->m_X1,&p_Lv1); CHKERRQ(ierr);
    ierr=VecRestoreArray(dvRinvx->m_X2,&p_Lv2); CHKERRQ(ierr);
    ierr=VecRestoreArray(dvRinvx->m_X3,&p_Lv3); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}





} // end of name space

#endif //_REGULARIZATIONREGISTRATION_CPP_
