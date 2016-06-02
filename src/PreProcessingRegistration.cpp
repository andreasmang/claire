#ifndef _PREPROCESSINGREGISTRATION_CPP_
#define _PREPROCESSINGREGISTRATION_CPP_


#include "PreProcessingRegistration.hpp"


namespace reg
{




/********************************************************************
 * Name: PreProcessingRegistration
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PreProcessingRegistration"
PreProcessingRegistration::PreProcessingRegistration()
{
    this->Initialize();
}




/********************************************************************
 * Name: PreProcessingRegistration
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PreProcessingRegistration"
PreProcessingRegistration::PreProcessingRegistration(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * Name: PreProcessingRegistration
 * Description: default deconstructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~PreProcessingRegistration"
PreProcessingRegistration::~PreProcessingRegistration()
{
    this->ClearMemory();
}



/********************************************************************
 * Name: Initialize
 * Description: initialize
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode PreProcessingRegistration::Initialize()
{
    this->m_Opt = NULL;
    this->m_IO = NULL;
    this->m_xhat = NULL;
    this->m_Kxhat = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clear memory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode PreProcessingRegistration::ClearMemory()
{
    //PetscErrorCode ierr;
    if(this->m_xhat!=NULL){
        accfft_free(this->m_xhat);
        this->m_xhat = NULL;
    }
    if(this->m_Kxhat!=NULL){
        accfft_free(this->m_Kxhat);
        this->m_Kxhat = NULL;
    }

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: SetIO
 * Description: set io interface for data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetIO"
PetscErrorCode PreProcessingRegistration::SetIO(PreProcessingRegistration::ReadWriteType* io)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr=Assert(io != NULL,"null pointer"); CHKERRQ(ierr);
    this->m_IO=io;

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: ApplyGaussianSmoothing
 * Description: apply gaussian smoothing operator to input data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyGaussianSmoothing"
PetscErrorCode PreProcessingRegistration::ApplyGaussianSmoothing(Vec y, Vec x)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3],n[3];
    IntType iosize[3];
    size_t alloc_max;
    ScalarType hx[3],nx[3],c[3],scale;
    ScalarType *p_x=NULL, *p_y=NULL;
    double ffttimers[5]={0,0,0,0,0};
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    n[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
    n[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
    n[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

    // get local pencil size and allocation size
    alloc_max=accfft_local_size_dft_r2c_t<ScalarType>(n,isize,istart,osize,ostart,
                                                        this->m_Opt->GetComm());
    if(this->m_xhat == NULL){
        this->m_xhat=(FFTScaType*)accfft_alloc(alloc_max);
    }

    for (int i = 0; i < 3; ++i){

        hx[i] = this->m_Opt->GetSpatialStepSize(i);
        nx[i] = static_cast<ScalarType>(n[i]);

        // sigma is provided by user in # of grid points
        c[i] = this->m_Opt->GetSigma()*hx[i];
        c[i] *= c[i];

        iosize[i] = osize[i];
    }

    scale = this->m_Opt->ComputeFFTScale();

    // compute fft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFTPlan(),p_x,this->m_xhat,ffttimers);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);


#pragma omp parallel
{
#pragma omp for
    for (IntType i1 = 0; i1 < iosize[0]; ++i1){ // x1
        for (IntType i2 = 0; i2 < iosize[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < iosize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                ScalarType k1 = static_cast<ScalarType>(i1 + ostart[0]);
                ScalarType k2 = static_cast<ScalarType>(i2 + ostart[1]);
                ScalarType k3 = static_cast<ScalarType>(i3 + ostart[2]);

                // check if grid index is larger or smaller then
                // half of the total grid size
                bool flagx1 = (k1 <= nx[0]*0.5);
                bool flagx2 = (k2 <= nx[1]*0.5);
                bool flagx3 = (k3 <= nx[2]*0.5);

                k1 = flagx1 ? k1 : -nx[0] + k1;
                k2 = flagx2 ? k2 : -nx[1] + k2;
                k3 = flagx3 ? k3 : -nx[2] + k3;

                ScalarType sik = 0.5*((k1*k1*c[0]) + (k2*k2*c[1]) + (k3*k3*c[2]));
                sik = exp(-sik);

                // compute linear / flat index
                IntType li = GetLinearIndex(i1,i2,i3,iosize);

                this->m_xhat[li][0] *= scale*sik;
                this->m_xhat[li][1] *= scale*sik;

            } // i1
        } // i2
    } // i3

} // pragma omp parallel

    // compute inverse fft
    ierr=VecGetArray(y,&p_y); CHKERRQ(ierr);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFTPlan(),this->m_xhat,p_y,ffttimers);
    ierr=VecRestoreArray(y,&p_y); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);

}




} // end of namespace

#endif // _PREPROCESSINGREGISTRATION_CPP_
