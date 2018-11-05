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
DifferentiationSM::DifferentiationSM(RegOpt* opt) : SuperClass(opt, Type::Spectral) {
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
    
    this->xyz[0] = 1;
    this->xyz[1] = 1;
    this->xyz[2] = 1;
    
    this->m_XHat[0] = nullptr;
    this->m_XHat[1] = nullptr;
    this->m_XHat[2] = nullptr;
    
    this->m_VecFieldKernel.pXHat[0] = nullptr;
    this->m_VecFieldKernel.pXHat[1] = nullptr;
    this->m_VecFieldKernel.pXHat[2] = nullptr;
    this->m_VecFieldKernel.nx[0] = 0;
    this->m_VecFieldKernel.nx[1] = 0;
    this->m_VecFieldKernel.nx[2] = 0;
    this->m_VecFieldKernel.nl[0] = 0;
    this->m_VecFieldKernel.nl[1] = 0;
    this->m_VecFieldKernel.nl[2] = 0;
    this->m_VecFieldKernel.nstart[0] = 0;
    this->m_VecFieldKernel.nstart[1] = 0;
    this->m_VecFieldKernel.nstart[2] = 0;
    this->m_VecFieldKernel.scale = 0;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DifferentiationSM::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = FreeMemory(this->m_XHat[0]); CHKERRQ(ierr);
    ierr = FreeMemory(this->m_XHat[1]); CHKERRQ(ierr);
    ierr = FreeMemory(this->m_XHat[2]); CHKERRQ(ierr);
    
    this->m_VecFieldKernel.pXHat[0] = nullptr;
    this->m_VecFieldKernel.pXHat[1] = nullptr;
    this->m_VecFieldKernel.pXHat[2] = nullptr;

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::SetupSpectralData(ComplexType *x1, ComplexType *x2, ComplexType *x3) {
    PetscErrorCode ierr = 0;
    IntType nalloc;
    PetscFunctionBegin;
    
    nalloc = this->m_Opt->m_FFT.nalloc;
    
    if (!x1) {
      ierr = AllocateMemoryOnce(this->m_XHat[0], nalloc); CHKERRQ(ierr);
      this->m_VecFieldKernel.pXHat[0] = this->m_XHat[0];
    } else {
      this->m_VecFieldKernel.pXHat[0] = x1;
    }
    if (!x2) {
      ierr = AllocateMemoryOnce(this->m_XHat[1], nalloc); CHKERRQ(ierr);
      this->m_VecFieldKernel.pXHat[1] = this->m_XHat[1];
    } else {
      this->m_VecFieldKernel.pXHat[1] = x2;
    }
    if (!x3) {
      ierr = AllocateMemoryOnce(this->m_XHat[2], nalloc); CHKERRQ(ierr);
      this->m_VecFieldKernel.pXHat[2] = this->m_XHat[2];
    } else {
      this->m_VecFieldKernel.pXHat[2] = x3;
    }
    this->m_VecFieldKernel.nx[0] = this->m_Opt->m_Domain.nx[0];
    this->m_VecFieldKernel.nx[1] = this->m_Opt->m_Domain.nx[1];
    this->m_VecFieldKernel.nx[2] = this->m_Opt->m_Domain.nx[2];
    this->m_VecFieldKernel.nl[0] = this->m_Opt->m_FFT.osize[0];
    this->m_VecFieldKernel.nl[1] = this->m_Opt->m_FFT.osize[1];
    this->m_VecFieldKernel.nl[2] = this->m_Opt->m_FFT.osize[2];
    this->m_VecFieldKernel.nstart[0] = this->m_Opt->m_FFT.ostart[0];
    this->m_VecFieldKernel.nstart[1] = this->m_Opt->m_FFT.ostart[1];
    this->m_VecFieldKernel.nstart[2] = this->m_Opt->m_FFT.ostart[2];
    this->m_VecFieldKernel.scale = this->m_Opt->ComputeFFTScale();
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Gradient(ScalarType *g1,
                                           ScalarType *g2,
                                           ScalarType *g3,
                                           const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Grad");
    
    ZeitGeist_define(FFT_GRAD);
    ZeitGeist_tick(FFT_GRAD);
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_grad_t(g1, g2, g3, const_cast<ScalarType*>(m), this->m_Opt->m_FFT.plan, &xyz, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTGRAD);
    
    this->m_Opt->IncreaseFFTTimers(timer);
    
    ZeitGeist_tock(FFT_GRAD);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Gradient(VecField *g, const ScalarType *m) {
    PetscErrorCode ierr = 0;
    ScalarType *g1 = nullptr, *g2 = nullptr, *g3 = nullptr;
    PetscFunctionBegin;
    
    ierr = g->GetArraysWrite(g1, g2, g3); CHKERRQ(ierr);
    
    ierr = this->Gradient(g1, g2, g3, m); CHKERRQ(ierr);
    
    ierr = g->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Gradient(ScalarType **g, const ScalarType *m) {
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
                                            const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Laplacian");
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_laplace_t(l, const_cast<ScalarType*>(m), this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);

    this->m_Opt->IncreaseFFTTimers(timer);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute laplacian of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Laplacian(ScalarType *l1,
                                            ScalarType *l2,
                                            ScalarType *l3,
                                            const ScalarType *v1,
                                            const ScalarType *v2,
                                            const ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Laplacian Field");
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_laplace_t(l1, const_cast<ScalarType*>(v1), this->m_Opt->m_FFT.plan, timer);
    accfft_laplace_t(l2, const_cast<ScalarType*>(v2), this->m_Opt->m_FFT.plan, timer);
    accfft_laplace_t(l3, const_cast<ScalarType*>(v3), this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    this->m_Opt->IncreaseFFTTimers(timer);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}
 
/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Divergence(ScalarType *l,
                                             const ScalarType *v1,
                                             const ScalarType *v2,
                                             const ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Divergence");
    
    ZeitGeist_define(FFT_DIV);
    ZeitGeist_tick(FFT_DIV);
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(l, 
      const_cast<ScalarType*>(v1), 
      const_cast<ScalarType*>(v2), 
      const_cast<ScalarType*>(v3), this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);
    
    this->m_Opt->IncreaseFFTTimers(timer);
    
    ZeitGeist_tock(FFT_DIV);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Divergence(ScalarType *l, VecField *v) {
    PetscErrorCode ierr = 0;
    const ScalarType *v1 = nullptr, *v2 = nullptr, *v3 = nullptr;
    PetscFunctionBegin;
    
    ierr = v->GetArraysRead(v1, v2, v3); CHKERRQ(ierr);
    
    ierr = this->Divergence(l, v1, v2, v3); CHKERRQ(ierr);
    
    ierr = v->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Divergence(ScalarType *l, const ScalarType *const *v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
        
    ierr = this->Divergence(l, v[0], v[1], v[2]); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute biharmonic operator of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Biharmonic(ScalarType *b,
                                             const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Biharmonic");
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_biharmonic_t(b, const_cast<ScalarType*>(m), this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    this->m_Opt->IncreaseFFTTimers(timer);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute biharmonic operator of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Biharmonic(ScalarType *b1,
                                             ScalarType *b2,
                                             ScalarType *b3,
                                             const ScalarType *v1,
                                             const ScalarType *v2,
                                             const ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Biharmonic Field");
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_biharmonic_t(b1, const_cast<ScalarType*>(v1), this->m_Opt->m_FFT.plan, timer);
    accfft_biharmonic_t(b2, const_cast<ScalarType*>(v2), this->m_Opt->m_FFT.plan, timer);
    accfft_biharmonic_t(b3, const_cast<ScalarType*>(v3), this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    this->m_Opt->IncreaseFFTTimers(timer);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::Laplacian(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    const ScalarType *pV[3] = {nullptr, nullptr, nullptr};
    ScalarType *pBV[3] = {nullptr, nullptr, nullptr};
    ComplexType **pXHat = this->m_VecFieldKernel.pXHat;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Bilaplacian Field");
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = v->GetArraysRead(pV); CHKERRQ(ierr);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[0]), pXHat[0], timer);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[1]), pXHat[1], timer);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[2]), pXHat[2], timer);
    ierr = v->RestoreArrays(); CHKERRQ(ierr);
    
    if (b1 == 0.0) {
      ierr = this->m_VecFieldKernel.Laplacian(b0); CHKERRQ(ierr);
    } else {
      ierr = this->m_VecFieldKernel.Laplacian(b0, b1); CHKERRQ(ierr);
    }
    
    ierr = bv->GetArraysWrite(pBV); CHKERRQ(ierr);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[0], pBV[0], timer);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[1], pBV[1], timer);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[2], pBV[2], timer);
    ierr = bv->RestoreArrays(); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 6);
    this->m_Opt->IncreaseFFTTimers(timer);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::Bilaplacian(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    const ScalarType *pV[3] = {nullptr, nullptr, nullptr};
    ScalarType *pBV[3] = {nullptr, nullptr, nullptr};
    ComplexType **pXHat = this->m_VecFieldKernel.pXHat;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Bilaplacian Field");
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = v->GetArraysRead(pV); CHKERRQ(ierr);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[0]), pXHat[0], timer);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[1]), pXHat[1], timer);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[2]), pXHat[2], timer);
    ierr = v->RestoreArrays(); CHKERRQ(ierr);
    
    if (b1 == 0.0) {
      ierr = this->m_VecFieldKernel.Bilaplacian(b0); CHKERRQ(ierr);
    } else {
      ierr = this->m_VecFieldKernel.Bilaplacian(b0, b1); CHKERRQ(ierr);
    }
    
    ierr = bv->GetArraysWrite(pBV); CHKERRQ(ierr);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[0], pBV[0], timer);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[1], pBV[1], timer);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[2], pBV[2], timer);
    ierr = bv->RestoreArrays(); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 6);
    this->m_Opt->IncreaseFFTTimers(timer);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::InverseBilaplacian(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    const ScalarType *pV[3] = {nullptr, nullptr, nullptr};
    ScalarType *pBV[3] = {nullptr, nullptr, nullptr};
    ComplexType **pXHat = this->m_VecFieldKernel.pXHat;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT inverse Bilaplacian Field");
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = v->GetArraysRead(pV); CHKERRQ(ierr);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[0]), pXHat[0], timer);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[1]), pXHat[1], timer);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[2]), pXHat[2], timer);
    ierr = v->RestoreArrays(); CHKERRQ(ierr);
    
    if (b1 == 0.0) {
      ierr = this->m_VecFieldKernel.InverseBilaplacian(b0); CHKERRQ(ierr);
    } else {
      ierr = this->m_VecFieldKernel.InverseBilaplacian(b0, b1); CHKERRQ(ierr);
    }
    
    ierr = bv->GetArraysWrite(pBV); CHKERRQ(ierr);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[0], pBV[0], timer);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[1], pBV[1], timer);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[2], pBV[2], timer);
    ierr = bv->RestoreArrays(); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 6);
    this->m_Opt->IncreaseFFTTimers(timer);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::InverseBilaplacianSqrt(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    const ScalarType *pV[3] = {nullptr, nullptr, nullptr};
    ScalarType *pBV[3] = {nullptr, nullptr, nullptr};
    ComplexType **pXHat = this->m_VecFieldKernel.pXHat;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT inverse Bilaplacian Field");
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = v->GetArraysRead(pV); CHKERRQ(ierr);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[0]), pXHat[0], timer);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[1]), pXHat[1], timer);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(pV[2]), pXHat[2], timer);
    ierr = v->RestoreArrays(); CHKERRQ(ierr);
    
    if (b1 == 0.0) {
      ierr = this->m_VecFieldKernel.InverseBilaplacianSqrt(b0); CHKERRQ(ierr);
    } else {
      ierr = this->m_VecFieldKernel.InverseBilaplacianSqrt(b0, b1); CHKERRQ(ierr);
    }
    
    ierr = bv->GetArraysWrite(pBV); CHKERRQ(ierr);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[0], pBV[0], timer);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[1], pBV[1], timer);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, pXHat[2], pBV[2], timer);
    ierr = bv->RestoreArrays(); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 6);
    this->m_Opt->IncreaseFFTTimers(timer);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}


}  // end of name space




#endif  // _DIFFERENTIATION_CPP_
