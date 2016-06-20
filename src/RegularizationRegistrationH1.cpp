#ifndef _REGULARIZATIONREGISTRATIONH1_CPP_
#define _REGULARIZATIONREGISTRATIONH1_CPP_

#include "RegularizationRegistrationH1.hpp"




namespace reg
{




/********************************************************************
 * Name: RegularizationRegistrationH1
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH1"
RegularizationRegistrationH1::RegularizationRegistrationH1() : SuperClass()
{
}




/********************************************************************
 * Name: ~RegularizationRegistrationH1
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegularizationRegistrationH1"
RegularizationRegistrationH1::~RegularizationRegistrationH1(void)
{
    this->ClearMemory();
}




/********************************************************************
 * Name: RegularizationRegistrationH1
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistrationH1"
RegularizationRegistrationH1::RegularizationRegistrationH1(RegOpt* opt) : SuperClass(opt)
{

}






/********************************************************************
 * Name: EvaluateFunctional
 * Description: evaluates the functional
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateFunctional"
PetscErrorCode RegularizationRegistrationH1::EvaluateFunctional(ScalarType* R, VecField* v)
{
    PetscErrorCode ierr;
    ScalarType  *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_gv11=NULL,*p_gv12=NULL,*p_gv13=NULL,
                *p_gv21=NULL,*p_gv22=NULL,*p_gv23=NULL,
                *p_gv31=NULL,*p_gv32=NULL,*p_gv33=NULL;
    ScalarType value,beta[2],H1v,L2v,hd;
    std::bitset<3>XYZ=0;XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    beta[0] = this->m_Opt->GetRegularizationWeight(0);
    beta[1] = this->m_Opt->GetRegularizationWeight(1);

    *R= 0.0;

    if (beta[0] != 0.0 && beta[1] != 0.0 ){

        // compute hd
        hd = this->m_Opt->GetLebesqueMeasure();

        if (this->m_WorkVecField==NULL){
            try{this->m_WorkVecField = new VecField(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        H1v = 0.0;

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
        accfft_grad(p_gv11,p_gv12,p_gv13,p_v1,this->m_Opt->GetFFT().plan,&XYZ,ffttimers);
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
        accfft_grad(p_gv21,p_gv22,p_gv23,p_v2,this->m_Opt->GetFFT().plan,&XYZ,ffttimers);
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
        accfft_grad(p_gv31,p_gv32,p_gv33,p_v3,this->m_Opt->GetFFT().plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

        ierr=VecRestoreArray(this->m_WorkVecField->m_X1,&p_gv31); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField->m_X2,&p_gv32); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField->m_X3,&p_gv33); CHKERRQ(ierr);

        // compute inner products
        ierr=VecTDot(this->m_WorkVecField->m_X1,this->m_WorkVecField->m_X1,&value); H1v +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X2,this->m_WorkVecField->m_X2,&value); H1v +=value;
        ierr=VecTDot(this->m_WorkVecField->m_X3,this->m_WorkVecField->m_X3,&value); H1v +=value;

        ierr=VecRestoreArray(v->m_X1,&p_v1); CHKERRQ(ierr);
        ierr=VecRestoreArray(v->m_X2,&p_v2); CHKERRQ(ierr);
        ierr=VecRestoreArray(v->m_X3,&p_v3); CHKERRQ(ierr);

        L2v = 0.0;
        ierr=VecTDot(v->m_X1,v->m_X1,&value); L2v +=value;
        ierr=VecTDot(v->m_X2,v->m_X2,&value); L2v +=value;
        ierr=VecTDot(v->m_X3,v->m_X3,&value); L2v +=value;

        // add up contributions
        *R = 0.5*hd*(beta[0]*H1v + beta[1]*L2v);

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
PetscErrorCode RegularizationRegistrationH1::EvaluateGradient(VecField* dvR, VecField* v)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3],nx[3];
    IntType iosize[3];
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_Lv1=NULL,*p_Lv2=NULL,*p_Lv3=NULL;
    ScalarType beta[2],scale,hd;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(dvR != NULL, "null pointer"); CHKERRQ(ierr);

    beta[0] = this->m_Opt->GetRegularizationWeight();
    beta[1] = this->m_Opt->GetRegularizationWeight();

    // if regularization weight is zero, do noting
    if ( beta[0]  == 0.0  && beta[1] == 0.0 ){
        ierr=VecSet(dvR->m_X1,0.0); CHKERRQ(ierr);
        ierr=VecSet(dvR->m_X2,0.0); CHKERRQ(ierr);
        ierr=VecSet(dvR->m_X3,0.0); CHKERRQ(ierr);
    }
    else{

        ierr=this->Allocate(); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        accfft_local_size_dft_r2c_t<ScalarType>(nx,isize,istart,osize,ostart,
                                                this->m_Opt->GetFFT().mpicomm);

        for (int i=0; i < 3; ++i){
            iosize[i] = static_cast<IntType>(osize[i]);
        }
        scale = this->m_Opt->ComputeFFTScale();

        ierr=VecGetArray(v->m_X1,&p_v1); CHKERRQ(ierr);
        ierr=VecGetArray(v->m_X2,&p_v2); CHKERRQ(ierr);
        ierr=VecGetArray(v->m_X3,&p_v3); CHKERRQ(ierr);

        // compute forward fft
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_v1,this->m_v1hat,ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_v2,this->m_v2hat,ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_v3,this->m_v3hat,ffttimers);
        this->m_Opt->IncrementCounter(FFT,3);

        ierr=VecRestoreArray(v->m_X1,&p_v1); CHKERRQ(ierr);
        ierr=VecRestoreArray(v->m_X2,&p_v2); CHKERRQ(ierr);
        ierr=VecRestoreArray(v->m_X3,&p_v3); CHKERRQ(ierr);


#pragma omp parallel
{
        long int w[3];
        ScalarType lapik,regop;
        IntType i;

#pragma omp for
        for (IntType i1 = 0; i1 < iosize[0]; ++i1){
            for (IntType i2 = 0; i2 < iosize[1]; i2++){
                for (IntType i3 = 0; i3 < iosize[2]; ++i3){

                    w[0] = static_cast<long int>(i1 + ostart[0]);
                    w[1] = static_cast<long int>(i2 + ostart[1]);
                    w[2] = static_cast<long int>(i3 + ostart[2]);

                    CheckWaveNumbers(w,nx);

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
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv1hat,p_Lv1,ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv2hat,p_Lv2,ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv3hat,p_Lv3,ffttimers);
        this->m_Opt->IncrementCounter(FFT,3);

        // restore arrays
        ierr=VecRestoreArray(dvR->m_X1,&p_Lv1); CHKERRQ(ierr);
        ierr=VecRestoreArray(dvR->m_X2,&p_Lv2); CHKERRQ(ierr);
        ierr=VecRestoreArray(dvR->m_X3,&p_Lv3); CHKERRQ(ierr);

        // compute hd
        hd = this->m_Opt->GetLebesqueMeasure();

        // scale vector field
        ierr=dvR->Scale(hd); CHKERRQ(ierr);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(ffttimers);
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
PetscErrorCode RegularizationRegistrationH1::HessianMatVec(VecField* dvvR, VecField* vtilde)
{
    PetscErrorCode ierr;
    ScalarType beta;
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
        // first and second variation are identical
        ierr=this->EvaluateGradient(dvvR,vtilde); CHKERRQ(ierr);
    }


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
PetscErrorCode RegularizationRegistrationH1::ApplyInverseOperator(VecField* Ainvx, VecField* x)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3],nx[3];
    IntType iosize[3];
    ScalarType *p_x1=NULL,*p_x2=NULL,*p_x3=NULL,
                *p_Ainvx1=NULL,*p_Ainvx2=NULL,*p_Ainvx3=NULL;
    ScalarType beta[2],scale;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(Ainvx != NULL, "null pointer"); CHKERRQ(ierr);

    beta[0] = this->m_Opt->GetRegularizationWeight();
    beta[1] = this->m_Opt->GetRegularizationWeight();

    // if regularization weight is zero, do noting
    if ( beta[0] == 0.0 && beta[1] == 0 ){

        ierr=VecCopy(x->m_X1,Ainvx->m_X1); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X2,Ainvx->m_X2); CHKERRQ(ierr);
        ierr=VecCopy(x->m_X3,Ainvx->m_X3); CHKERRQ(ierr);

    }
    else{

        ierr=this->Allocate(); CHKERRQ(ierr);

        nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
        nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
        nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

        accfft_local_size_dft_r2c_t<ScalarType>(nx,isize,istart,osize,ostart,
                                                this->m_Opt->GetFFT().mpicomm);

        for (unsigned int i=0; i < 3; ++i){
            iosize[i] = static_cast<IntType>(osize[i]);
        }
        scale = this->m_Opt->ComputeFFTScale();

        ierr=VecGetArray(x->m_X1,&p_x1); CHKERRQ(ierr);
        ierr=VecGetArray(x->m_X2,&p_x2); CHKERRQ(ierr);
        ierr=VecGetArray(x->m_X3,&p_x3); CHKERRQ(ierr);

        // compute forward fft
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x1,this->m_v1hat,ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x2,this->m_v2hat,ffttimers);
        accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x3,this->m_v3hat,ffttimers);
        this->m_Opt->IncrementCounter(FFT,3);

        ierr=VecRestoreArray(x->m_X1,&p_x1); CHKERRQ(ierr);
        ierr=VecRestoreArray(x->m_X2,&p_x2); CHKERRQ(ierr);
        ierr=VecRestoreArray(x->m_X3,&p_x3); CHKERRQ(ierr);

#pragma omp parallel
{
        long int w[3];
        ScalarType lapik,regop;
        IntType i;

#pragma omp for
        for (IntType i1 = 0; i1 < iosize[0]; ++i1){
            for (IntType i2 = 0; i2 < iosize[1]; i2++){
                for (IntType i3 = 0; i3 < iosize[2]; ++i3){

                    w[0] = static_cast<long int>(i1 + ostart[0]);
                    w[1] = static_cast<long int>(i2 + ostart[1]);
                    w[2] = static_cast<long int>(i3 + ostart[2]);

                    CheckWaveNumbersInv(w,nx);

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

        ierr=VecGetArray(Ainvx->m_X1,&p_Ainvx1); CHKERRQ(ierr);
        ierr=VecGetArray(Ainvx->m_X2,&p_Ainvx2); CHKERRQ(ierr);
        ierr=VecGetArray(Ainvx->m_X3,&p_Ainvx3); CHKERRQ(ierr);

        // compute inverse fft
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv1hat,p_Ainvx1,ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv2hat,p_Ainvx2,ffttimers);
        accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Lv3hat,p_Ainvx3,ffttimers);
        this->m_Opt->IncrementCounter(FFT,3);

        // restore arrays
        ierr=VecRestoreArray(Ainvx->m_X1,&p_Ainvx1); CHKERRQ(ierr);
        ierr=VecRestoreArray(Ainvx->m_X2,&p_Ainvx2); CHKERRQ(ierr);
        ierr=VecRestoreArray(Ainvx->m_X3,&p_Ainvx3); CHKERRQ(ierr);

        // increment fft timer
        this->m_Opt->IncreaseFFTTimers(ffttimers);

    }

    PetscFunctionReturn(0);
}




} // end of name space

#endif //_REGULARIZATIONREGISTRATIONH1_CPP_
