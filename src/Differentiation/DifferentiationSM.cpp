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
    
    if (opt->m_Verbosity > 2) {
      DbgMsg("DifferentiationSM created");
    }
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
    
    this->m_SpectralKernel.pXHat[0] = nullptr;
    this->m_SpectralKernel.pXHat[1] = nullptr;
    this->m_SpectralKernel.pXHat[2] = nullptr;
    this->m_SpectralKernel.nx[0] = 0;
    this->m_SpectralKernel.nx[1] = 0;
    this->m_SpectralKernel.nx[2] = 0;
    this->m_SpectralKernel.nl[0] = 0;
    this->m_SpectralKernel.nl[1] = 0;
    this->m_SpectralKernel.nl[2] = 0;
    this->m_SpectralKernel.nstart[0] = 0;
    this->m_SpectralKernel.nstart[1] = 0;
    this->m_SpectralKernel.nstart[2] = 0;
    this->m_SpectralKernel.scale = 0;

//#ifdef REG_HAS_CUDA
//    this->m_planC2R = nullptr;
//    this->m_planR2C = nullptr;
//#endif

    this->m_FFT = nullptr;
    
    ierr = this->SetupData(); CHKERRQ(ierr);
        
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
    
    this->m_SpectralKernel.pXHat[0] = nullptr;
    this->m_SpectralKernel.pXHat[1] = nullptr;
    this->m_SpectralKernel.pXHat[2] = nullptr;
    
/*#ifdef REG_HAS_CUDA
    cufftDestroy(*this->m_planC2R);
    cufftDestroy(*this->m_planR2C);
    ierr = Free(this->m_planC2R); CHKERRQ(ierr);
    ierr = Free(this->m_planR2C); CHKERRQ(ierr);
#endif*/
    
    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::SetupData(ComplexType *x1, ComplexType *x2, ComplexType *x3) {
    PetscErrorCode ierr = 0;
    IntType nc;
    PetscFunctionBegin;
    
    this->m_FFT = &this->m_Opt->m_FFT;
    
    nc = this->m_FFT->nalloc/(2*sizeof(ScalarType));
    
     this->m_SpectralKernel.pXHat[0] = &this->m_FFT->fft->m_WorkSpace[0*nc];
     this->m_SpectralKernel.pXHat[1] = &this->m_FFT->fft->m_WorkSpace[1*nc];
     this->m_SpectralKernel.pXHat[2] = &this->m_FFT->fft->m_WorkSpace[2*nc];
    /*if (!x1) {
      ierr = AllocateMemoryOnce(this->m_XHat[0], nalloc); CHKERRQ(ierr);
      this->m_SpectralKernel.pXHat[0] = this->m_XHat[0];
    } else {
      this->m_SpectralKernel.pXHat[0] = x1;
    }
    if (!x2) {
      ierr = AllocateMemoryOnce(this->m_XHat[1], nalloc); CHKERRQ(ierr);
      this->m_SpectralKernel.pXHat[1] = this->m_XHat[1];
    } else {
      this->m_SpectralKernel.pXHat[1] = x2;
    }
    if (!x3) {
      ierr = AllocateMemoryOnce(this->m_XHat[2], nalloc); CHKERRQ(ierr);
      this->m_SpectralKernel.pXHat[2] = this->m_XHat[2];
    } else {
      this->m_SpectralKernel.pXHat[2] = x3;
    }*/
    this->m_SpectralKernel.nx[0] = this->m_FFT->nx[0];
    this->m_SpectralKernel.nx[1] = this->m_FFT->nx[1];
    this->m_SpectralKernel.nx[2] = this->m_FFT->nx[2];
    this->m_SpectralKernel.nl[0] = this->m_FFT->osize[0];
    this->m_SpectralKernel.nl[1] = this->m_FFT->osize[1];
    this->m_SpectralKernel.nl[2] = this->m_FFT->osize[2];
    this->m_SpectralKernel.nstart[0] = this->m_FFT->ostart[0];
    this->m_SpectralKernel.nstart[1] = this->m_FFT->ostart[1];
    this->m_SpectralKernel.nstart[2] = this->m_FFT->ostart[2];
    this->m_SpectralKernel.scale = 1./(this->m_FFT->nx[0]*this->m_FFT->nx[1]*this->m_FFT->nx[2]);
        
/*#ifdef REG_HAS_CUDA
    if (this->m_planR2C == nullptr) {
      ierr = AllocateOnce(this->m_planR2C); CHKERRQ(ierr);
      cufftPlan3d(this->m_planR2C, 
        this->m_SpectralKernel.nx[0],  this->m_SpectralKernel.nx[1],  this->m_SpectralKernel.nx[2],
        CUFFT_R2C);
    }
    if (this->m_planC2R == nullptr) {
      ierr = AllocateOnce(this->m_planC2R); CHKERRQ(ierr);
      cufftPlan3d(this->m_planC2R, 
        this->m_SpectralKernel.nx[0],  this->m_SpectralKernel.nx[1],  this->m_SpectralKernel.nx[2],
        CUFFT_C2R);
    }
#endif*/
        
    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::SetFFT(FourierTransform* fft) {
  PetscErrorCode ierr = 0;
  IntType nc;
  PetscFunctionBegin;
  this->m_FFT = fft;
  
  nc = this->m_FFT->nalloc/(2*sizeof(ScalarType));
    
  this->m_SpectralKernel.pXHat[0] = &this->m_FFT->fft->m_WorkSpace[0*nc];
  this->m_SpectralKernel.pXHat[1] = &this->m_FFT->fft->m_WorkSpace[1*nc];
  this->m_SpectralKernel.pXHat[2] = &this->m_FFT->fft->m_WorkSpace[2*nc];
  
  this->m_SpectralKernel.nx[0] = this->m_FFT->nx[0];
  this->m_SpectralKernel.nx[1] = this->m_FFT->nx[1];
  this->m_SpectralKernel.nx[2] = this->m_FFT->nx[2];
  this->m_SpectralKernel.nl[0] = this->m_FFT->osize[0];
  this->m_SpectralKernel.nl[1] = this->m_FFT->osize[1];
  this->m_SpectralKernel.nl[2] = this->m_FFT->osize[2];
  this->m_SpectralKernel.nstart[0] = this->m_FFT->ostart[0];
  this->m_SpectralKernel.nstart[1] = this->m_FFT->ostart[1];
  this->m_SpectralKernel.nstart[2] = this->m_FFT->ostart[2];
  this->m_SpectralKernel.scale = 1./(this->m_FFT->nx[0]*this->m_FFT->nx[1]*this->m_FFT->nx[2]);
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
    
    //ComplexType **pXHat = this->m_SpectralKernel.pXHat;
    
    DebugGPUStartEvent("FFT Grad");
    
    ZeitGeist_define(FFT_GRAD);
    ZeitGeist_tick(FFT_GRAD);
    
//    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    //accfft_grad_t(g1, g2, g3, const_cast<ScalarType*>(m), this->m_Opt->m_FFT->plan, &xyz, timer);
    
    //accfft_execute_r2c_t(this->m_Opt->m_FFT->plan, const_cast<ScalarType*>(m), pXHat[0], timer);    
    ierr = this->ComputeForwardFFT(m);
    ierr = this->m_SpectralKernel.Gradient(); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(g1, g2, g3);
    //accfft_execute_c2r_t(this->m_Opt->m_FFT->plan, pXHat[0], g1, timer);
    //accfft_execute_c2r_t(this->m_Opt->m_FFT->plan, pXHat[1], g2, timer);
    //accfft_execute_c2r_t(this->m_Opt->m_FFT->plan, pXHat[2], g3, timer);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTGRAD);
    
//    this->m_Opt->IncreaseFFTTimers(timer);
    
    ZeitGeist_tock(FFT_GRAD);
    
    DebugGPUStopEvent();

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
    
    ZeitGeist_define(FFT_LAP);
    ZeitGeist_tick(FFT_LAP);
    
//    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    ierr = this->ComputeForwardFFT(m); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.ScalarLaplacian(-1.); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(l); CHKERRQ(ierr);
    //accfft_laplace_t(l, const_cast<ScalarType*>(m), this->m_Opt->m_FFT->plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);

//    this->m_Opt->IncreaseFFTTimers(timer);
    
    ZeitGeist_tock(FFT_LAP);
    
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
    
    ZeitGeist_define(FFT_LAP);
    ZeitGeist_tick(FFT_LAP);
    
//    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    ierr = this->ComputeForwardFFT(v1, v2, v3); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.Laplacian(-1.); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(l1, l2, l3); CHKERRQ(ierr);
    //accfft_laplace_t(l1, const_cast<ScalarType*>(v1), this->m_Opt->m_FFT->plan, timer);
    //accfft_laplace_t(l2, const_cast<ScalarType*>(v2), this->m_Opt->m_FFT->plan, timer);
    //accfft_laplace_t(l3, const_cast<ScalarType*>(v3), this->m_Opt->m_FFT->plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
//    this->m_Opt->IncreaseFFTTimers(timer);
    ZeitGeist_tock(FFT_LAP);
    
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
    
//    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    ierr = this->ComputeForwardFFT(v1, v2, v3); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.Divergence(); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(l);
    /*accfft_divergence_t(l, 
      const_cast<ScalarType*>(v1), 
      const_cast<ScalarType*>(v2), 
      const_cast<ScalarType*>(v3), this->m_Opt->m_FFT->plan, timer);*/
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);
    
//    this->m_Opt->IncreaseFFTTimers(timer);
    
    ZeitGeist_tock(FFT_DIV);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::RegLapModOp(VecField* bv, VecField* v, ScalarType b0) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT modified laplacian regularization operator");
    
    ZeitGeist_define(FFT_REG);
    ZeitGeist_tick(FFT_REG);
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.LaplacianMod(b0); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    ZeitGeist_tock(FFT_REG);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::RegLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT laplacian regularization operator");
    
    ZeitGeist_define(FFT_REG);
    ZeitGeist_tick(FFT_REG);
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    /*if (this->m_FFT->threshold > 0.) {
      ierr = this->m_SpectralKernel.LaplacianTol(b0, b1); CHKERRQ(ierr);
    } else {*/
      ierr = this->m_SpectralKernel.Laplacian(b0, b1); CHKERRQ(ierr);
    //}
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    ZeitGeist_tock(FFT_REG);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::RegBiLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT inverse bilaplacian regularization operator");
    
    ZeitGeist_define(FFT_REG);
    ZeitGeist_tick(FFT_REG);
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.Bilaplacian(b0, b1); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    ZeitGeist_tock(FFT_REG);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::RegTriLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT trilaplacian regularization operator");
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.Trilaplacian(b0, b1); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::RegTriLapFunc(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT trilaplacian regularization functional");
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.TrilaplacianFunctional(b0, b1); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::InvRegLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT inverse laplacian regularization operator");
    
    ZeitGeist_define(FFT_INVREG);
    ZeitGeist_tick(FFT_INVREG);
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.InverseLaplacian(usesqrt, b0, b1); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    ZeitGeist_tock(FFT_INVREG);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::InvRegBiLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT inverse bilaplacian regularization operator");
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.InverseBilaplacian(usesqrt, b0, b1); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::InvRegTriLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT inverse trilaplacian regularization operator");
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.InverseTrilaplacian(usesqrt, b0, b1); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::LerayOperator(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT leray operator");
    
    ZeitGeist_define(FFT_PROJ);
    ZeitGeist_tick(FFT_PROJ);
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.Leray(b0, b1); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    ZeitGeist_tock(FFT_PROJ);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::InvRegLerayOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1, ScalarType b2) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT leray operator");
    
    ZeitGeist_define(FFT_H0);
    ZeitGeist_tick(FFT_H0);
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.InvRegLeray(b0, b1, b2); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    ZeitGeist_tock(FFT_H0);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::GaussianFilter(ScalarType* bv, const ScalarType* v, const ScalarType c[3]) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FFT Gaussian filter");
    
    ZeitGeist_define(FFT_GAUSS);
    ZeitGeist_tick(FFT_GAUSS);
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    
    ierr = this->ComputeForwardFFT(v); CHKERRQ(ierr);
    ierr = this->m_SpectralKernel.GaussianFilter(c); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(bv); CHKERRQ(ierr);
    
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    ZeitGeist_tock(FFT_GAUSS);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::ComputeForwardFFT(VecField* v) {
    PetscErrorCode ierr = 0;
    const ScalarType *pV[3] = {nullptr, nullptr, nullptr};
    //ComplexType **pXHat = this->m_SpectralKernel.pXHat;
    PetscFunctionBegin;
        
//    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    /*if (this->m_FFT->threshold > 0) {
      ScalarType value;
      VecNorm(v->m_X1, NORM_INFINITY, &value);
      this->m_SpectralKernel.tol = value;
      VecNorm(v->m_X2, NORM_INFINITY, &value);
      this->m_SpectralKernel.tol = std::max(this->m_SpectralKernel.tol, std::abs(value));
      VecNorm(v->m_X3, NORM_INFINITY, &value);
      this->m_SpectralKernel.tol = std::max(this->m_SpectralKernel.tol, std::abs(value));
      this->m_SpectralKernel.tol *= this->m_FFT->threshold;
    }*/
        
    ierr = v->GetArraysRead(pV); CHKERRQ(ierr);
    ierr = this->ComputeForwardFFT(pV[0], pV[1], pV[2]); CHKERRQ(ierr);
    ierr = v->RestoreArrays(); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::ComputeForwardFFT(const ScalarType* pV1, const ScalarType *pV2, const ScalarType *pV3) {
    PetscErrorCode ierr = 0;
    ComplexType **pXHat = this->m_SpectralKernel.pXHat;
    PetscFunctionBegin;
    
    this->m_FFT->fft->FFT_R2C(pV1, pXHat[0]);
    this->m_FFT->fft->FFT_R2C(pV2, pXHat[1]);
    this->m_FFT->fft->FFT_R2C(pV3, pXHat[2]);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::ComputeInverseFFT(VecField *v) {
    PetscErrorCode ierr = 0;
    ScalarType *pV[3] = {nullptr, nullptr, nullptr};
    PetscFunctionBegin;
        
    ierr = v->GetArraysWrite(pV); CHKERRQ(ierr);
    ierr = this->ComputeInverseFFT(pV[0], pV[1], pV[2]); CHKERRQ(ierr);
    ierr = v->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationSM::ComputeInverseFFT(ScalarType* pV1, ScalarType *pV2, ScalarType *pV3) {
    PetscErrorCode ierr = 0;
    ComplexType **pXHat = this->m_SpectralKernel.pXHat;
    PetscFunctionBegin;

    this->m_FFT->fft->FFT_C2R(pXHat[0], pV1);
    this->m_FFT->fft->FFT_C2R(pXHat[1], pV2);
    this->m_FFT->fft->FFT_C2R(pXHat[2], pV3);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::ComputeForwardFFT(const ScalarType* v) {
    PetscErrorCode ierr = 0;
    ComplexType **pXHat = this->m_SpectralKernel.pXHat;
    PetscFunctionBegin;
        
    this->m_FFT->fft->FFT_R2C(v, pXHat[0]);
    
    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::ComputeInverseFFT(ScalarType* v) {
    PetscErrorCode ierr = 0;
    ComplexType **pXHat = this->m_SpectralKernel.pXHat;
    PetscFunctionBegin;

    this->m_FFT->fft->FFT_C2R(pXHat[0], v);

    DebugGPUStopEvent();

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::Restrict(ScalarType* vc, const ScalarType* vf, FourierTransform* coarse) {
  PetscErrorCode ierr = 0;
  ComplexType *pXHat[2];
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_RESTRICT);
  ZeitGeist_tick(FFT_RESTRICT);
  ScalarType scale = 1./(this->m_FFT->nx[0]*this->m_FFT->nx[1]*this->m_FFT->nx[2]);
  
  IntType nc = coarse->nalloc/(2*sizeof(ScalarType));
  
  pXHat[0] = &this->m_FFT->fft->m_WorkSpace[0];
  pXHat[1] = &coarse->fft->m_WorkSpace[nc];
  
  this->m_FFT->fft->FFT_R2C(vf, pXHat[0]);
  this->m_FFT->fft->Restrict(pXHat[1], pXHat[0], coarse->fft);
  coarse->fft->Scale(pXHat[1], scale);
  coarse->fft->FFT_C2R(pXHat[1], vc);
  
  ZeitGeist_tock(FFT_RESTRICT);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::Prolong(ScalarType* vf, const ScalarType* vc, FourierTransform* coarse) {
  PetscErrorCode ierr = 0;
  ComplexType *pXHat[2];
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_PROLONG);
  ZeitGeist_tick(FFT_PROLONG);
  
  ScalarType scale = 1./(coarse->nx[0]*coarse->nx[1]*coarse->nx[2]);
  
  IntType nc = coarse->nalloc/(2*sizeof(ScalarType));
  
  pXHat[1] = &this->m_FFT->fft->m_WorkSpace[0];
  pXHat[0] = &coarse->fft->m_WorkSpace[nc];
  
  coarse->fft->FFT_R2C(vc, pXHat[0]);
  //this->m_FFT->fft->FFT_R2C(vf, pXHat[1]);
  coarse->fft->Scale(pXHat[0], scale);
  //this->m_FFT->fft->Scale(pXHat[1], 1./(this->m_FFT->nx[0]*this->m_FFT->nx[1]*this->m_FFT->nx[2]));
  this->m_FFT->fft->Prolong(pXHat[1], pXHat[0], coarse->fft);
  this->m_FFT->fft->FFT_C2R(pXHat[1], vf);
  
  ZeitGeist_tock(FFT_PROLONG);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::RestrictH0(VecField* vc, VecField* vf, FourierTransform* coarse, ScalarType beta) {
  PetscErrorCode ierr = 0;
  ComplexType **pXHat = this->m_SpectralKernel.pXHat;
  ComplexType *pXHat_c[3];
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_H0RESTRICT);
  ZeitGeist_tick(FFT_H0RESTRICT);
  //ScalarType scale = 1./(this->m_FFT->nx[0]*this->m_FFT->nx[1]*this->m_FFT->nx[2]);
  
  const ScalarType *pVf[3] = {nullptr, nullptr, nullptr};
  ScalarType *pVc[3] = {nullptr, nullptr, nullptr};
  
  IntType nc = coarse->nalloc/(2*sizeof(ScalarType));
  
  pXHat_c[0] = &coarse->fft->m_WorkSpace[0*nc];
  pXHat_c[1] = &coarse->fft->m_WorkSpace[1*nc];
  pXHat_c[2] = &coarse->fft->m_WorkSpace[2*nc];
  
  ierr = vf->GetArraysRead(pVf); CHKERRQ(ierr);
  this->m_FFT->fft->FFT_R2C(pVf[0], pXHat[0]);
  this->m_FFT->fft->FFT_R2C(pVf[1], pXHat[1]);
  this->m_FFT->fft->FFT_R2C(pVf[2], pXHat[2]);
  ierr = vf->RestoreArrays(); CHKERRQ(ierr);
  
  ierr = this->m_SpectralKernel.InverseLaplacian(false, beta, 0); CHKERRQ(ierr);
  
  //this->m_FFT->fft->Scale(pXHat[0], scale);
  //this->m_FFT->fft->Scale(pXHat[1], scale);
  //this->m_FFT->fft->Scale(pXHat[2], scale);
  
  this->m_FFT->fft->Restrict(pXHat_c[0], pXHat[0], coarse->fft);
  this->m_FFT->fft->Restrict(pXHat_c[1], pXHat[1], coarse->fft);
  this->m_FFT->fft->Restrict(pXHat_c[2], pXHat[2], coarse->fft);
  
  ierr = vc->GetArraysWrite(pVc); CHKERRQ(ierr);
  coarse->fft->FFT_C2R(pXHat_c[0], pVc[0]);
  coarse->fft->FFT_C2R(pXHat_c[1], pVc[1]);
  coarse->fft->FFT_C2R(pXHat_c[2], pVc[2]);
  ierr = vc->RestoreArrays(); CHKERRQ(ierr);
  
  ZeitGeist_tock(FFT_H0RESTRICT);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationSM::ProlongH0(VecField* vf, VecField* vc, FourierTransform* coarse) {
  PetscErrorCode ierr = 0;
  ComplexType **pXHat = this->m_SpectralKernel.pXHat;
  ComplexType *pXHat_c[3];
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_H0PROLONG);
  ZeitGeist_tick(FFT_H0PROLONG);
  
  ScalarType scale = 1./(coarse->nx[0]*coarse->nx[1]*coarse->nx[2]);
  
  const ScalarType *pVc[3] = {nullptr, nullptr, nullptr};
  ScalarType *pVf[3] = {nullptr, nullptr, nullptr};
  
  IntType nc = coarse->nalloc/(2*sizeof(ScalarType));
  
  pXHat_c[0] = &coarse->fft->m_WorkSpace[0*nc];
  pXHat_c[1] = &coarse->fft->m_WorkSpace[1*nc];
  pXHat_c[2] = &coarse->fft->m_WorkSpace[2*nc];
  
  ierr = vc->GetArraysRead(pVc); CHKERRQ(ierr);
  coarse->fft->FFT_R2C(pVc[0], pXHat_c[0]);
  coarse->fft->FFT_R2C(pVc[1], pXHat_c[1]);
  coarse->fft->FFT_R2C(pVc[2], pXHat_c[2]);
  ierr = vc->RestoreArrays(); CHKERRQ(ierr);
  
  coarse->fft->Scale(pXHat_c[0], scale);
  coarse->fft->Scale(pXHat_c[1], scale);
  coarse->fft->Scale(pXHat_c[2], scale);
  
  this->m_FFT->fft->ProlongMerge(pXHat[0], pXHat_c[0], coarse->fft);
  this->m_FFT->fft->ProlongMerge(pXHat[1], pXHat_c[1], coarse->fft);
  this->m_FFT->fft->ProlongMerge(pXHat[2], pXHat_c[2], coarse->fft);
  
  ierr = vf->GetArraysWrite(pVf); CHKERRQ(ierr);
  this->m_FFT->fft->FFT_C2R(pXHat[0], pVf[0]);
  this->m_FFT->fft->FFT_C2R(pXHat[1], pVf[1]);
  this->m_FFT->fft->FFT_C2R(pXHat[2], pVf[2]);
  ierr = vf->RestoreArrays(); CHKERRQ(ierr);
  
  ZeitGeist_tock(FFT_H0PROLONG);
  
  PetscFunctionReturn(ierr);
}

}  // end of name space




#endif  // _DIFFERENTIATION_CPP_
