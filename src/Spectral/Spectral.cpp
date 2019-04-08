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

#ifndef _SPECTRAL_CPP_
#define _SPECTRAL_CPP_

#include "Spectral.hpp"
#include "RegOpt.hpp"

namespace reg {

/********************************************************************
 * @brief default destructor
 *******************************************************************/
Spectral::~Spectral() {
    this->ClearMemory();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
Spectral::Spectral(RegOpt *opt) {
    this->Initialize();
    this->m_Opt = opt;
    this->SetupFFT();
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Spectral::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    this->m_kernel.nx[0] = 0;
    this->m_kernel.nx[1] = 0;
    this->m_kernel.nx[2] = 0;
    this->m_kernel.nl[0] = 0;
    this->m_kernel.nl[1] = 0;
    this->m_kernel.nl[2] = 0;
    this->m_kernel.nstart[0] = 0;
    this->m_kernel.nstart[1] = 0;
    this->m_kernel.nstart[2] = 0;
    this->m_kernel.scale = 0;

#ifdef REG_HAS_CUDA
    this->m_planC2R = nullptr;
    this->m_planR2C = nullptr;
#endif
    this->m_plan = nullptr;

    this->m_Opt = nullptr;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Spectral::InitFFT() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    this->m_kernel.nx[0] = this->m_Opt->m_Domain.nx[0];
    this->m_kernel.nx[1] = this->m_Opt->m_Domain.nx[1];
    this->m_kernel.nx[2] = this->m_Opt->m_Domain.nx[2];
    this->m_kernel.nl[0] = this->m_Opt->m_FFT.osize[0];
    this->m_kernel.nl[1] = this->m_Opt->m_FFT.osize[1];
    this->m_kernel.nl[2] = this->m_Opt->m_FFT.osize[2];
    this->m_kernel.nstart[0] = this->m_Opt->m_FFT.ostart[0];
    this->m_kernel.nstart[1] = this->m_Opt->m_FFT.ostart[1];
    this->m_kernel.nstart[2] = this->m_Opt->m_FFT.ostart[2];
    this->m_kernel.scale = this->m_Opt->ComputeFFTScale();
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Spectral::SetupFFT() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    int nx[3], nalloc;
    nx[0] = this->m_Opt->m_Domain.nx[0];
    nx[1] = this->m_Opt->m_Domain.nx[1];
    nx[2] = this->m_Opt->m_Domain.nx[2];
    nalloc = this->m_Opt->m_FFT.nalloc;
    
#ifdef REG_HAS_CUDA
    if (this->m_planR2C == nullptr) {
      ierr = AllocateOnce(this->m_planR2C); CHKERRQ(ierr);
      cufftPlan3d(this->m_planR2C, nx[0],  nx[1], nx[2], CUFFT_R2C);
    }
    if (this->m_planC2R == nullptr) {
      ierr = AllocateOnce(this->m_planC2R); CHKERRQ(ierr);
      cufftPlan3d(this->m_planC2R, nx[0], nx[1], nx[2],  CUFFT_C2R);
    }
#else
    std::stringstream ss;
    ScalarType *u = nullptr;
    ComplexType *uk = nullptr;
    
    if (this->m_Opt->m_Verbosity > 2) {
        ss << " >> " << __func__ << ": allocation (size = " << nalloc << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }
    u = reinterpret_cast<ScalarType*>(accfft_alloc(nalloc));
    ierr = Assert(u != nullptr, "allocation failed"); CHKERRQ(ierr);

    // set up the fft
    if (this->m_Opt->m_Verbosity > 2) {
        ss << " >> " << __func__ << ": allocation (size = " << nalloc << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }
    uk = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
    ierr = Assert(uk != nullptr, "allocation failed"); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg("allocating fft plan"); CHKERRQ(ierr);
    }
    
    this->m_plan = accfft_plan_dft_3d_r2c(nx, u, reinterpret_cast<ScalarType*>(uk),
                                              this->m_Opt->m_FFT.mpicomm, ACCFFT_MEASURE);
    ierr = Assert(this->m_plan != nullptr, "allocation failed"); CHKERRQ(ierr);
    
        // clean up
    if (u != NULL) {accfft_free(u); u = NULL;}
    if (uk != NULL) {accfft_free(uk); uk = NULL;}
#endif
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Spectral::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
#ifdef REG_HAS_CUDA
    cufftDestroy(*this->m_planC2R);
    cufftDestroy(*this->m_planR2C);
    ierr = Free(this->m_planC2R); CHKERRQ(ierr);
    ierr = Free(this->m_planR2C); CHKERRQ(ierr);
#endif
    if(this->m_plan) {
      accfft_destroy_plan(this->m_plan);
    }
    this->m_plan = nullptr;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief forward FFT
 *******************************************************************/
PetscErrorCode Spectral::FFT_R2C(const ScalarType *real, ComplexType *complex) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
#ifdef REG_HAS_CUDA
    cufftExecR2C(*this->m_planR2C, const_cast<cufftReal*>(real), reinterpret_cast<cufftComplex*>(complex));
    cudaDeviceSynchronize();
#else
    double timer[NFFTTIMERS] = {0};
    accfft_execute_r2c_t(this->m_plan, const_cast<ScalarType*>(real), complex, timer);
    this->m_Opt->IncrementCounter(FFT, 1);
    this->m_Opt->IncreaseFFTTimers(timer);
#endif

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief inverse FFT
 *******************************************************************/
PetscErrorCode Spectral::FFT_C2R(const ComplexType *complex, ScalarType *real) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
#ifdef REG_HAS_CUDA
    cufftExecC2R(*this->m_planC2R, const_cast<cufftComplex*>(reinterpret_cast<const cufftComplex*>(complex)), reinterpret_cast<cufftReal*>(real));
    cudaDeviceSynchronize();
#else
    double timer[NFFTTIMERS] = {0};
    accfft_execute_c2r_t(this->m_plan, const_cast<ComplexType*>(complex), real, timer);
    this->m_Opt->IncrementCounter(FFT, 1);
    this->m_Opt->IncreaseFFTTimers(timer);
#endif

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Low pass filter
 *******************************************************************/
PetscErrorCode Spectral::LowPassFilter(ComplexType *xHat, ScalarType pct) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = this->m_kernel.LowPassFilter(xHat, pct); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Low pass filter
 *******************************************************************/
PetscErrorCode Spectral::HighPassFilter(ComplexType *xHat, ScalarType pct) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = this->m_kernel.LowPassFilter(xHat, pct); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Restrict to lower Grid
 *******************************************************************/
PetscErrorCode Spectral::Restrict(ComplexType *xc, const ComplexType *xf, const IntType nxc[3]) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    IntType nx_c[3];
    nx_c[0] = nxc[0]; nx_c[1] = nxc[1]; nx_c[2] = nxc[2];
    
    ierr = this->m_kernel.Restrict(xc, xf, nx_c); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Prolong from lower Grid
 *******************************************************************/
PetscErrorCode Spectral::Prolong(ComplexType *xf, const ComplexType *xc, const IntType nxc[3]) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    IntType nx_c[3];
    nx_c[0] = nxc[0]; nx_c[1] = nxc[1]; nx_c[2] = nxc[2];
    
    ierr = this->m_kernel.Prolong(xf, xc, nx_c); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Prolong from lower Grid
 *******************************************************************/
PetscErrorCode Spectral::Scale(ComplexType *x, ScalarType scale) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = this->m_kernel.Scale(x, scale); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


}  // end of namespace

#endif  // _SPECTRAL_CPP_
