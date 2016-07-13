#ifndef _REGULARIZATIONREGISTRATION_CPP_
#define _REGULARIZATIONREGISTRATION_CPP_

#include "RegularizationRegistration.hpp"




namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistration"
RegularizationRegistration::RegularizationRegistration()
{
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegularizationRegistration"
RegularizationRegistration::~RegularizationRegistration(void)
{
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegularizationRegistration"
RegularizationRegistration::RegularizationRegistration(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
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

    this->m_InvRegOpMinEigval=0.0;
    this->m_InvRegOpMaxEigval=0.0;
    this->m_RegOpMaxEigval=0.0;
    this->m_RegOpMinEigval=0.0;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
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
 * @brief allocate arrays for fft (we might have to do this
 * in several functions, so we do it collectively here)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode RegularizationRegistration::Allocate(void)
{
    int isize[3],osize[3],istart[3],ostart[3],nx[3];
    size_t alloc_max;

    PetscFunctionBegin;

    nx[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
    nx[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
    nx[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

    // get local pencil size and allocation size
    alloc_max=accfft_local_size_dft_r2c_t<ScalarType>(nx,isize,istart,osize,ostart,
                                                        this->m_Opt->GetFFT().mpicomm);

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
 * @brief deallocate arrays for fft (we might have to do this
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



} // end of name space

#endif //_REGULARIZATIONREGISTRATION_CPP_
