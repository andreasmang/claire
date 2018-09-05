/*************************************************************************
 *  Copyright (c) 2018.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CLAIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _DIFFERENTIATIONSM_CPP_
#define _DIFFERENTIATIONSM_CPP_

#include "DifferentiationSM.hpp"

namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DifferentiationSM::DifferentiationSM() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DifferentiationSM::DifferentiationSM(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DifferentiationSM::~DifferentiationSM() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DifferentiationSM::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DifferentiationSM::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Gradient(ScalarType *g1,
                                           ScalarType *g2,
                                           ScalarType *g3,
                                           ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ZeitGeist_define(FFT_GRAD);
    ZeitGeist_tick(FFT_GRAD);
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_grad_t(g1, g2, g3, m, this->m_Opt->m_FFT.plan, &xyz, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTGRAD);
    
    this->m_Opt->IncreaseFFTTimers(timer);
    
    ZeitGeist_tock(FFT_GRAD);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Gradient(VecField *g, ScalarType *m) {
    PetscErrorCode ierr = 0;
    ScalarType *g1 = nullptr, *g2 = nullptr, *g3 = nullptr;
    PetscFunctionBegin;
    
    ierr = g->GetArraysReadWrite(g1, g2, g3); CHKERRQ(ierr);
    
    ierr = this->Gradient(g1, g2, g3, m); CHKERRQ(ierr);
    
    ierr = g->RestoreArraysReadWrite(g1, g2, g3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Gradient(ScalarType **g, ScalarType *m) {
    PetscErrorCode ierr = 0;
    ScalarType *g1 = nullptr, *g2 = nullptr, *g3 = nullptr;
    PetscFunctionBegin;
        
    ierr = this->Gradient(g[0], g[1], g[2], m); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute laplacian of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Laplacian(ScalarType *l,
                                            ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_laplace_t(l, m, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);

    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute laplacian of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Laplacian(ScalarType *l1,
                                            ScalarType *l2,
                                            ScalarType *l3,
                                            ScalarType *v1,
                                            ScalarType *v2,
                                            ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_laplace_t(l1, v1, this->m_Opt->m_FFT.plan, timer);
    accfft_laplace_t(l2, v2, this->m_Opt->m_FFT.plan, timer);
    accfft_laplace_t(l3, v3, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}
 
/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Divergence(ScalarType *l,
                                             ScalarType *v1,
                                             ScalarType *v2,
                                             ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ZeitGeist_define(FFT_DIV);
    ZeitGeist_tick(FFT_DIV);
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(l, v1, v2, v3, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);
    
    this->m_Opt->IncreaseFFTTimers(timer);
    
    ZeitGeist_tock(FFT_DIV);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Divergence(ScalarType *l, VecField *v) {
    PetscErrorCode ierr = 0;
    ScalarType *v1 = nullptr, *v2 = nullptr, *v3 = nullptr;
    PetscFunctionBegin;
    
    ierr = v->GetArraysReadWrite(v1, v2, v3); CHKERRQ(ierr);
    
    ierr = this->Divergence(l, v1, v2, v3); CHKERRQ(ierr);
    
    ierr = v->RestoreArraysReadWrite(v1, v2, v3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Divergence(ScalarType *l, ScalarType **v) {
    PetscErrorCode ierr = 0;
    ScalarType *v1 = nullptr, *v2 = nullptr, *v3 = nullptr;
    PetscFunctionBegin;
        
    ierr = this->Divergence(l, v[0], v[1], v[2]); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief compute biharmonic operator of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Biharmonic(ScalarType *b,
                                             ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_biharmonic_t(b, m, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute biharmonic operator of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Biharmonic(ScalarType *b1,
                                             ScalarType *b2,
                                             ScalarType *b3,
                                             ScalarType *v1,
                                             ScalarType *v2,
                                             ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_biharmonic_t(b1, v1, this->m_Opt->m_FFT.plan, timer);
    accfft_biharmonic_t(b2, v2, this->m_Opt->m_FFT.plan, timer);
    accfft_biharmonic_t(b3, v3, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}


}  // end of name space




#endif  // _DIFFERENTIATION_CPP_
