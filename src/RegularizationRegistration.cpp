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

    this->m_Dv11hat = NULL;
    this->m_Dv12hat = NULL;
    this->m_Dv13hat = NULL;
    this->m_Dv21hat = NULL;
    this->m_Dv22hat = NULL;
    this->m_Dv23hat = NULL;
    this->m_Dv31hat = NULL;
    this->m_Dv32hat = NULL;
    this->m_Dv33hat = NULL;

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
PetscErrorCode RegularizationRegistration::Allocate(int flag)
{
    PetscFunctionBegin;

    if (flag==0){
        if(this->m_v1hat == NULL){
            this->m_v1hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_v2hat == NULL){
            this->m_v2hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_v3hat == NULL){
            this->m_v3hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
    }
    else if (flag == 1){
        if(this->m_Lv1hat == NULL){
            this->m_Lv1hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_Lv2hat == NULL){
            this->m_Lv2hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_Lv3hat == NULL){
            this->m_Lv3hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
    }
    else if (flag == 2){
        if(this->m_Dv11hat == NULL){
            this->m_Dv11hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_Dv12hat == NULL){
            this->m_Dv12hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_Dv13hat == NULL){
            this->m_Dv13hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }

        if(this->m_Dv21hat == NULL){
            this->m_Dv21hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_Dv22hat == NULL){
            this->m_Dv22hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_Dv23hat == NULL){
            this->m_Dv23hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }

        if(this->m_Dv31hat == NULL){
            this->m_Dv31hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_Dv32hat == NULL){
            this->m_Dv32hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
        if(this->m_Dv33hat == NULL){
            this->m_Dv33hat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
        }
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

    if(this->m_Dv11hat!=NULL){
        accfft_free(this->m_Dv11hat);
        this->m_Dv11hat = NULL;
    }
    if(this->m_Dv12hat!=NULL){
        accfft_free(this->m_Dv12hat);
        this->m_Dv12hat = NULL;
    }
    if(this->m_Dv13hat!=NULL){
        accfft_free(this->m_Dv13hat);
        this->m_Dv13hat = NULL;
    }
    if(this->m_Dv21hat!=NULL){
        accfft_free(this->m_Dv21hat);
        this->m_Dv21hat = NULL;
    }
    if(this->m_Dv22hat!=NULL){
        accfft_free(this->m_Dv22hat);
        this->m_Dv22hat = NULL;
    }
    if(this->m_Dv23hat!=NULL){
        accfft_free(this->m_Dv23hat);
        this->m_Dv23hat = NULL;
    }
    if(this->m_Dv31hat!=NULL){
        accfft_free(this->m_Dv31hat);
        this->m_Dv31hat = NULL;
    }
    if(this->m_Dv32hat!=NULL){
        accfft_free(this->m_Dv32hat);
        this->m_Dv32hat = NULL;
    }
    if(this->m_Dv33hat!=NULL){
        accfft_free(this->m_Dv33hat);
        this->m_Dv33hat = NULL;
    }

    PetscFunctionReturn(0);
}



} // end of name space

#endif //_REGULARIZATIONREGISTRATION_CPP_
