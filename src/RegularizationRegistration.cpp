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
    ierr = this->Deallocate(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief allocate arrays for fft (we might have to do this
 * in several functions, so we do it collectively here)
 *******************************************************************/
PetscErrorCode RegularizationRegistration::Allocate() {
    PetscFunctionBegin;

    if (this->m_v1hat == NULL) {
        this->m_v1hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
    }
    if (this->m_v2hat == NULL) {
        this->m_v2hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
    }
    if (this->m_v3hat == NULL) {
        this->m_v3hat = reinterpret_cast<FFTScaType*>(accfft_alloc(this->m_Opt->GetFFT().nalloc));
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

    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _REGULARIZATIONREGISTRATION_CPP_
