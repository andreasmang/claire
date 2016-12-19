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

#ifndef _REGULARIZATIONREGISTRATION_CPP_
#define _REGULARIZATIONREGISTRATION_CPP_

#include "RegularizationRegistration.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegularizationRegistration::RegularizationRegistration() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
RegularizationRegistration::~RegularizationRegistration(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
RegularizationRegistration::RegularizationRegistration(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode RegularizationRegistration::Initialize(void) {
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
PetscErrorCode RegularizationRegistration::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_WorkVecField!=NULL) {
        delete this->m_WorkVecField;
        this->m_WorkVecField = NULL;
    }

    // clear the fft arrays
    ierr=this->Deallocate(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief allocate arrays for fft (we might have to do this
 * in several functions, so we do it collectively here)
 *******************************************************************/
PetscErrorCode RegularizationRegistration::Allocate(int flag) {
    PetscFunctionBegin;

    if (flag==0) {
        if (this->m_v1hat == NULL) {
            this->m_v1hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_v2hat == NULL) {
            this->m_v2hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_v3hat == NULL) {
            this->m_v3hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
    } else if (flag == 1) {
        if (this->m_Lv1hat == NULL) {
            this->m_Lv1hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_Lv2hat == NULL) {
            this->m_Lv2hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_Lv3hat == NULL) {
            this->m_Lv3hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
    } else if (flag == 2) {
        if (this->m_Dv11hat == NULL) {
            this->m_Dv11hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_Dv12hat == NULL) {
            this->m_Dv12hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_Dv13hat == NULL) {
            this->m_Dv13hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }

        if (this->m_Dv21hat == NULL) {
            this->m_Dv21hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_Dv22hat == NULL) {
            this->m_Dv22hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_Dv23hat == NULL) {
            this->m_Dv23hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }

        if (this->m_Dv31hat == NULL) {
            this->m_Dv31hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_Dv32hat == NULL) {
            this->m_Dv32hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
        if (this->m_Dv33hat == NULL) {
            this->m_Dv33hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
        }
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief deallocate arrays for fft (we might have to do this
 * in several functions, so we do it collectively here)
 *******************************************************************/
PetscErrorCode RegularizationRegistration::Deallocate(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_v1hat != NULL) {
        accfft_free(this->m_v1hat);
        this->m_v1hat = NULL;
    }
    if (this->m_v2hat != NULL) {
        accfft_free(this->m_v2hat);
        this->m_v2hat = NULL;
    }
    if (this->m_v3hat != NULL) {
        accfft_free(this->m_v3hat);
        this->m_v3hat = NULL;
    }

    if (this->m_Lv1hat != NULL) {
        accfft_free(this->m_Lv1hat);
        this->m_Lv1hat = NULL;
    }
    if (this->m_Lv2hat != NULL) {
        accfft_free(this->m_Lv2hat);
        this->m_Lv2hat = NULL;
    }
    if (this->m_Lv3hat != NULL) {
        accfft_free(this->m_Lv3hat);
        this->m_Lv3hat = NULL;
    }

    if (this->m_Dv11hat != NULL) {
        accfft_free(this->m_Dv11hat);
        this->m_Dv11hat = NULL;
    }
    if (this->m_Dv12hat != NULL) {
        accfft_free(this->m_Dv12hat);
        this->m_Dv12hat = NULL;
    }
    if (this->m_Dv13hat != NULL) {
        accfft_free(this->m_Dv13hat);
        this->m_Dv13hat = NULL;
    }
    if (this->m_Dv21hat != NULL) {
        accfft_free(this->m_Dv21hat);
        this->m_Dv21hat = NULL;
    }
    if (this->m_Dv22hat != NULL) {
        accfft_free(this->m_Dv22hat);
        this->m_Dv22hat = NULL;
    }
    if (this->m_Dv23hat != NULL) {
        accfft_free(this->m_Dv23hat);
        this->m_Dv23hat = NULL;
    }
    if (this->m_Dv31hat != NULL) {
        accfft_free(this->m_Dv31hat);
        this->m_Dv31hat = NULL;
    }
    if (this->m_Dv32hat != NULL) {
        accfft_free(this->m_Dv32hat);
        this->m_Dv32hat = NULL;
    }
    if (this->m_Dv33hat != NULL) {
        accfft_free(this->m_Dv33hat);
        this->m_Dv33hat = NULL;
    }

    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _REGULARIZATIONREGISTRATION_CPP_
