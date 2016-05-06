#ifndef _OPTIMALCONTROLREGISTRATIONRIC_CPP_
#define _OPTIMALCONTROLREGISTRATIONRIC_CPP_

#include <math.h>
#include "OptimalControlRegistrationRIC.hpp"




namespace reg
{




/********************************************************************
 * Name: OptimalControlRegistrationRIC
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistrationRIC"
OptimalControlRegistrationRIC::OptimalControlRegistrationRIC() : SuperClass()
{
    this->Initialize();
}




/********************************************************************
 * Name: OptimalControlRegistrationRIC
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~OptimalControlRegistrationRIC"
OptimalControlRegistrationRIC::~OptimalControlRegistrationRIC(void)
{
    this->ClearMemory();
}




/********************************************************************
 * Name: OptimalControlRegistrationRIC
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistrationRIC"
OptimalControlRegistrationRIC::OptimalControlRegistrationRIC(RegOpt* opt) : SuperClass(opt)
{
    this->Initialize();
}




/********************************************************************
 * Name: Initialize
 * Description: init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode OptimalControlRegistrationRIC::Initialize(void)
{
    PetscFunctionBegin;

    this->m_x1hat=NULL;
    this->m_x2hat=NULL;
    this->m_x3hat=NULL;

    this->m_Kx1hat=NULL;
    this->m_Kx2hat=NULL;
    this->m_Kx3hat=NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode OptimalControlRegistrationRIC::ClearMemory(void)
{
//    PetscErrorCode ierr;
    PetscFunctionBegin;

    if(this->m_x1hat!=NULL){
        accfft_free(this->m_x1hat);
        this->m_x1hat = NULL;
    }
    if(this->m_x2hat!=NULL){
        accfft_free(this->m_x2hat);
        this->m_x2hat = NULL;
    }
    if(this->m_x3hat!=NULL){
        accfft_free(this->m_x3hat);
        this->m_x3hat = NULL;
    }
    if(this->m_Kx1hat!=NULL){
        accfft_free(this->m_Kx1hat);
        this->m_Kx1hat = NULL;
    }
    if(this->m_Kx2hat!=NULL){
        accfft_free(this->m_Kx2hat);
        this->m_Kx2hat = NULL;
    }
    if(this->m_Kx3hat!=NULL){
        accfft_free(this->m_Kx3hat);
        this->m_Kx3hat = NULL;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeBodyForce
 * Description: compute the body force
 * b = K[\int_0^1 \igrad m \lambda d t],
 * where K is an operator that projects v onto the manifold of
 * divergence free velocity fields
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeBodyForce"
PetscErrorCode OptimalControlRegistrationRIC::ComputeBodyForce()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_WorkVecField1 == NULL){
        this->m_WorkVecField1 = new VecField(this->m_Opt);
    }
    if (this->m_WorkVecField2 == NULL){
        this->m_WorkVecField2 = new VecField(this->m_Opt);
    }

    // assigned to work vec field 2
    ierr=SuperClass::ComputeBodyForce(); CHKERRQ(ierr);

    ierr=this->m_WorkVecField1->Copy(this->m_WorkVecField2); CHKERRQ(ierr);
    ierr=this->ApplyProjection(this->m_WorkVecField1); CHKERRQ(ierr);

    ierr=VecAXPY(this->m_WorkVecField2->m_X1,-1.0,this->m_WorkVecField1->m_X1); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField2->m_X2,-1.0,this->m_WorkVecField1->m_X2); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField2->m_X3,-1.0,this->m_WorkVecField1->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeIncBodyForce
 * Description: compute the body force
 * b = K[\int_0^1\igrad\tilde{m}\lambda+\igrad m \tilde{\lambda} dt]
 * where K is an operator that projects \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeIncBodyForce"
PetscErrorCode OptimalControlRegistrationRIC::ComputeIncBodyForce()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_WorkVecField1 == NULL){
        this->m_WorkVecField1 = new VecField(this->m_Opt);
    }
    if (this->m_WorkVecField2 == NULL){
        this->m_WorkVecField2 = new VecField(this->m_Opt);
    }

    // assigned to work vec field 2
    ierr=SuperClass::ComputeIncBodyForce(); CHKERRQ(ierr);

    ierr=this->m_WorkVecField1->Copy(this->m_WorkVecField2); CHKERRQ(ierr);
    ierr=this->ApplyProjection(this->m_WorkVecField1); CHKERRQ(ierr);

    ierr=VecAXPY(this->m_WorkVecField2->m_X1,-1.0,this->m_WorkVecField1->m_X1); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField2->m_X2,-1.0,this->m_WorkVecField1->m_X2); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField2->m_X3,-1.0,this->m_WorkVecField1->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ApplyProjection
 * Description: apply projection to map \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyProjection"
PetscErrorCode OptimalControlRegistrationRIC::ApplyProjection(VecField* x)
{
    PetscErrorCode ierr;
    ScalarType *p_x1=NULL, *p_x2=NULL, *p_x3=NULL, nx[3], scale;
    int isize[3],osize[3],istart[3],ostart[3];
    IntType alloc_max;
    double ffttimers[5]={0,0,0,0,0};
    int *n;

    PetscFunctionBegin;

    n = this->m_Opt->m_MiscOpt->N;

    scale = 1.0;
    for (unsigned int i=0; i < 3; ++i){
        nx[i] = static_cast<ScalarType>(n[i]);
        scale *= nx[i];
    }
    scale  = 1.0/scale;

    // get local pencil size and allocation size
    alloc_max=accfft_local_size_dft_r2c_t<ScalarType>(n,isize,istart,osize,ostart,this->m_Opt->m_MiscOpt->c_comm);

    if(this->m_x1hat == NULL){
        this->m_x1hat=(FFTScaType*)accfft_alloc(alloc_max);
    }
    if(this->m_x2hat == NULL){
        this->m_x2hat=(FFTScaType*)accfft_alloc(alloc_max);
    }
    if(this->m_x3hat == NULL){
        this->m_x3hat=(FFTScaType*)accfft_alloc(alloc_max);
    }

    if(this->m_Kx1hat == NULL){
        this->m_Kx1hat=(FFTScaType*)accfft_alloc(alloc_max);
    }
    if(this->m_Kx2hat == NULL){
        this->m_Kx2hat=(FFTScaType*)accfft_alloc(alloc_max);
    }
    if(this->m_Kx3hat == NULL){
        this->m_Kx3hat=(FFTScaType*)accfft_alloc(alloc_max);
    }

    ierr=VecGetArray(x->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecGetArray(x->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecGetArray(x->m_X3,&p_x3); CHKERRQ(ierr);

    // compute forward fft
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_x1,this->m_x1hat,ffttimers);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_x2,this->m_x2hat,ffttimers);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->m_MiscOpt->plan,p_x3,this->m_x3hat,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

#pragma omp parallel
{
    long int x1,x2,x3,wx1,wx2,wx3;
    ScalarType lapinvik,gradik1,gradik2,gradik3;
    long int i;
#pragma omp for
    for (unsigned int i1 = 0; i1 < osize[0]; ++i1){
        for (unsigned int i2 = 0; i2 < osize[1]; i2++){
            for (unsigned int i3 = 0; i3 < osize[2]; ++i3){

                x1 = static_cast<long int>(i1 + ostart[0]);
                x2 = static_cast<long int>(i2 + ostart[1]);
                x3 = static_cast<long int>(i3 + ostart[2]);

                // set wavenumber
                wx1 = x1;
                wx2 = x2;
                wx3 = x3;

                if(x1 > n[0]/2) wx1-=n[0];
                if(x2 > n[1]/2) wx2-=n[1];
                if(x3 > n[2]/2) wx3-=n[2];

                // compute inverse laplacian operator
                lapinvik = static_cast<ScalarType>(wx1*wx1 + wx2*wx2 + wx3*wx3);
                //lapinvik = round(lapinvik) == 0.0 ? -1.0 : 1.0/lapinvik;
                lapinvik = lapinvik == 0.0 ? -1.0 : -1.0/lapinvik;

                if(x1 == n[0]/2) wx1 = 0;
                if(x2 == n[1]/2) wx2 = 0;
                if(x3 == n[2]/2) wx3 = 0;

                // compute gradient operator
                gradik1 = static_cast<ScalarType>(wx1);
                gradik2 = static_cast<ScalarType>(wx2);
                gradik3 = static_cast<ScalarType>(wx3);

                i=(i1*osize[1]+i2)*osize[2]+i3;

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
}// pragma omp parallel


    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Kx1hat,p_x1,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Kx2hat,p_x2,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->m_MiscOpt->plan,this->m_Kx3hat,p_x3,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

    ierr=VecRestoreArray(x->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecRestoreArray(x->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecRestoreArray(x->m_X3,&p_x3); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}


} // end of namespace


#endif // _OPTIMALCONTROLREGISTRATIONIC_CPP_
