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
#include "VecField.hpp"

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
Spectral::Spectral(RegOpt *opt, FourierTransform *fft) {
    this->Initialize();
    this->m_Opt = opt;
    this->m_FFT = fft;
    this->SetupFFT();
    
    if (opt->m_Verbosity > 2) {
      DbgMsg("Spectral created");
    }
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

    this->m_plan = nullptr;
    this->m_WorkSpace = nullptr;
    
    this->m_SharedWorkSpace = false;

    this->m_Opt = nullptr;
    this->m_FFT = nullptr;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Spectral::InitFFT() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    this->m_kernel.nx[0] = this->m_FFT->nx[0];
    this->m_kernel.nx[1] = this->m_FFT->nx[1];
    this->m_kernel.nx[2] = this->m_FFT->nx[2];
    this->m_kernel.nl[0] = this->m_FFT->osize[0];
    this->m_kernel.nl[1] = this->m_FFT->osize[1];
    this->m_kernel.nl[2] = this->m_FFT->osize[2];
    this->m_kernel.nstart[0] = this->m_FFT->ostart[0];
    this->m_kernel.nstart[1] = this->m_FFT->ostart[1];
    this->m_kernel.nstart[2] = this->m_FFT->ostart[2];
    this->m_kernel.scale = 1./(this->m_FFT->nx[0]*this->m_FFT->nx[1]*this->m_FFT->nx[2]);
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Spectral::SetDomain() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
#ifdef REG_HAS_CUDA
    size_t isize[3], istart[3];
    this->m_plan->getInSize(isize);
    this->m_plan->getInStart(istart);
    this->m_Opt->m_Domain.nl = 1;
    this->m_Opt->m_Domain.ng = 1;
    for (int i = 0; i < 3; ++i) {
        this->m_Opt->m_Domain.isize[i] = static_cast<IntType>(isize[i]);
        this->m_Opt->m_Domain.istart[i] = static_cast<IntType>(istart[i]);
        this->m_Opt->m_Domain.nl *= this->m_Opt->m_Domain.isize[i];
        this->m_Opt->m_Domain.ng *= this->m_Opt->m_Domain.nx[i];
    }
    int np[2], periods[2], coord[2];
    MPI_Cart_get(this->m_Opt->m_Domain.mpicomm, 2, np, periods, coord);
    MPI_Comm_split(this->m_Opt->m_Domain.mpicomm, coord[0], coord[1], &this->m_Opt->m_Domain.rowcomm);
    MPI_Comm_split(this->m_Opt->m_Domain.mpicomm, coord[1], coord[0], &this->m_Opt->m_Domain.colcomm);
#else
    int isize[3], istart[3], osize[3], ostart[3];
    // get sizes (n is an integer, so it can overflow)
    accfft_local_size_dft_r2c_t<ScalarType>(nx, isize, istart, osize, ostart, this->m_Domain.mpicomm);
    
    this->m_Opt->m_Domain.nl = 1;
    this->m_Opt->m_Domain.ng = 1;
    for (int i = 0; i < 3; ++i) {
        this->m_Opt->m_Domain.nl *= static_cast<IntType>(isize[i]);
        this->m_Opt->m_Domain.ng *= this->m_Opt->m_Domain.nx[i];
        this->m_Opt->m_Domain.isize[i]  = static_cast<IntType>(isize[i]);
        this->m_Opt->m_Domain.istart[i] = static_cast<IntType>(istart[i]);
    }
    this->m_Opt->m_Domain.rowcomm = this->m_plan->row_comm;
    this->m_Opt->m_Domain.colcomm = this->m_plan->col_comm;
#endif
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Spectral::SetupFFT() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    int nx[3];
    nx[0] = this->m_FFT->nx[0];
    nx[1] = this->m_FFT->nx[1];
    nx[2] = this->m_FFT->nx[2];
    
#ifdef REG_HAS_CUDA
    {
      bool allocmem = true;
      void *sharedmem_d = nullptr, *sharedmem_h = nullptr;
      size_t sharedsize_d = 0, sharedsize_h = 0;
      if (this->m_Opt->m_FFT.fft && this->m_Opt->m_FFT.fft != this && this->m_Opt->m_FFT.fft->m_plan) {
        sharedmem_d = this->m_Opt->m_FFT.fft->m_plan->getWorkAreaDevice();
        sharedmem_h = this->m_Opt->m_FFT.fft->m_plan->getWorkAreaHost();
        sharedsize_d = this->m_Opt->m_FFT.fft->m_plan->getWorkSizeDevice();
        sharedsize_h = this->m_Opt->m_FFT.fft->m_plan->getWorkSizeHost();
        allocmem = false;
      }
#ifdef REG_HAS_MPICUDA
      ierr = AllocateOnce(m_plan, this->m_Opt->m_Domain.mpicomm, true); CHKERRQ(ierr);
#else
      ierr = AllocateOnce(m_plan, this->m_Opt->m_Domain.mpicomm, false); CHKERRQ(ierr);
#endif
      size_t osize[3], ostart[3], isize[3], istart[3];
      this->m_plan->initFFT(nx[0], nx[1], nx[2], allocmem);
      if (sharedsize_d < this->m_plan->getWorkSizeDevice()) {
        sharedmem_d = nullptr;
      }
      if (sharedsize_h < this->m_plan->getWorkSizeHost()) {
        sharedmem_h = nullptr;
      }
      if (!allocmem) {
        this->m_plan->setWorkArea(sharedmem_d, sharedmem_h);
      }
      if (allocmem && this->m_Opt->m_Verbosity > 2) {
        std::stringstream ss;
        ss << "FFT allocated " << this->m_plan->getWorkSizeDevice() << " bytes";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
      }
      this->m_FFT->nalloc = this->m_plan->getDomainSize();
      /*if (sharedalloc > this->m_FFT->nalloc) {
        this->m_WorkSpace = this->m_Opt->m_FFT.fft->m_WorkSpace;
        m_SharedWorkSpace = true;
      } else*/ {
        ierr = AllocateMemoryOnce(this->m_WorkSpace, this->m_FFT->nalloc*3);
        m_SharedWorkSpace = false;
      }
      this->m_plan->getOutSize(osize);
      this->m_plan->getOutStart(ostart);
      this->m_plan->getInSize(isize);
      this->m_plan->getInStart(istart);
      for (int i = 0; i < 3; ++i) {
        this->m_FFT->osize[i] = static_cast<IntType>(osize[i]);
        this->m_FFT->ostart[i] = static_cast<IntType>(ostart[i]);
        this->m_FFT->isize[i] = static_cast<IntType>(isize[i]);
        this->m_FFT->istart[i] = static_cast<IntType>(istart[i]);
      }
    }
#else
    {
      int isize[3], istart[3], osize[3], ostart[3];
      // get sizes (n is an integer, so it can overflow)
      this->m_FFT->nalloc = accfft_local_size_dft_r2c_t<ScalarType>(nx, isize, istart, osize, ostart, this->m_Domain.mpicomm);
      
      for (int i = 0; i < 3; ++i) {
        this->m_FFT.osize[i]  = static_cast<IntType>(osize[i]);
        this->m_FFT.ostart[i] = static_cast<IntType>(ostart[i]);
        this->m_FFT.isize[i]  = static_cast<IntType>(isize[i]);
        this->m_FFT.istart[i] = static_cast<IntType>(istart[i]);
      }
      
      std::stringstream ss;
      ScalarType *u = nullptr;
      ComplexType *uk = nullptr;
      
      if (this->m_Opt->m_Verbosity > 2) {
          ss << " >> " << __func__ << ": allocation (size = " << nalloc << ")";
          ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
          ss.clear(); ss.str(std::string());
      }
      ierr = AllocateMemoryOnce(this->m_WorkSpace, 3*nalloc); CHKERRQ(ierr);
      u = this->m_WorkSpace;
      uk = &this->m_WorkSpace[nalloc];
      //ierr = AllocateMemoryOnce(u, nalloc); CHKERRQ(ierr);
      //ierr = Assert(u != nullptr, "allocation failed"); CHKERRQ(ierr);

      // set up the fft
      if (this->m_Opt->m_Verbosity > 2) {
          ss << " >> " << __func__ << ": allocation (size = " << nalloc << ")";
          ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
          ss.clear(); ss.str(std::string());
      }
      //ierr = AllocateMemoryOnce(uk, nalloc); CHKERRQ(ierr);
      //ierr = Assert(uk != nullptr, "allocation failed"); CHKERRQ(ierr);

      if (this->m_Opt->m_Verbosity > 2) {
          ierr = DbgMsg("allocating fft plan"); CHKERRQ(ierr);
      }
      
      this->m_plan = accfft_plan_dft_3d_r2c(nx, u, reinterpret_cast<ScalarType*>(uk),
                                                this->m_Opt->m_Domain.mpicomm, ACCFFT_MEASURE);
      ierr = Assert(this->m_plan != nullptr, "allocation failed"); CHKERRQ(ierr);
      
          // clean up
      //ierr = FreeMemory(u); CHKERRQ(ierr);
      //ierr = FreeMemory(uk); CHKERRQ(ierr);
    }
#endif
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Spectral::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (!this->m_SharedWorkSpace) {
      ierr = FreeMemory(this->m_WorkSpace); CHKERRQ(ierr);
    }
    
#ifdef REG_HAS_CUDA
    ierr = Free(this->m_plan); CHKERRQ(ierr);
#else
    if(this->m_plan) {
      accfft_destroy_plan(this->m_plan);
    }
#endif
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
    {
      this->m_plan->execR2C(complex, real);
    }
#else
    {
      double timer[NFFTTIMERS] = {0};
      accfft_execute_r2c_t(this->m_plan, const_cast<ScalarType*>(real), complex, timer);
      this->m_Opt->IncrementCounter(FFT, 1);
      this->m_Opt->IncreaseFFTTimers(timer);
    }
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
    {
      this->m_plan->execC2R(real, complex);
    }
#else
    {
      double timer[NFFTTIMERS] = {0};
      accfft_execute_c2r_t(this->m_plan, const_cast<ComplexType*>(complex), real, timer);
      this->m_Opt->IncrementCounter(FFT, 1);
      this->m_Opt->IncreaseFFTTimers(timer);
    }
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

PetscErrorCode Spectral::Norm(ScalarType &norm, ComplexType *x, Spectral* size) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ScalarType *pWS[3];
  
  this->m_WorkVecField->GetArraysWrite(pWS);
  
  this->m_kernel.pWS = pWS[0];
  
  ScalarType lnorm = 0;
  
  IntType w[3];
  
  if (!size) size = this;
  
  w[0] = size->m_FFT->nx[0]/2;
  w[1] = size->m_FFT->nx[1]/2;
  w[2] = size->m_FFT->nx[2]/2;
  
  this->m_kernel.Norm(lnorm, x, w);
  
  MPI_Allreduce(MPI_IN_PLACE, &lnorm, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
  
  norm = lnorm;
  
  this->m_WorkVecField->RestoreArrays();

  PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Restrict to lower Grid
 *******************************************************************/
PetscErrorCode Spectral::Restrict(ComplexType *xc, const ComplexType *xf, Spectral *fft_coarse) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

#if REG_HAS_CUDA
    this->m_plan->restrictTo(xc, xf, fft_coarse->m_plan);
    //ierr = this->m_kernel.Restrict(xc, xf, fft_coarse->m_FFT->nx, fft_coarse->m_FFT->osize, fft_coarse->m_FFT->ostart); CHKERRQ(ierr);
#else
    if (this->m_Opt->rank_cnt == 1) {
      ierr = this->m_kernel.Restrict(xc, xf, fft_coarse->m_FFT->nx, fft_coarse->m_FFT->osize, fft_coarse->m_FFT->ostart); CHKERRQ(ierr);
    } else {
      ierr = ThrowError("Spectral restriction not implemented!"); CHKERRQ(ierr);
    }
#endif

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief Prolong from lower Grid
 *******************************************************************/
PetscErrorCode Spectral::Prolong(ComplexType *xf, const ComplexType *xc, Spectral *fft_coarse) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

#if REG_HAS_CUDA
    this->m_plan->prolongFrom(xf, xc, fft_coarse->m_plan);
    //ierr = this->m_kernel.Prolong(xf, xc, fft_coarse->m_FFT->nx, fft_coarse->m_FFT->osize, fft_coarse->m_FFT->ostart); CHKERRQ(ierr);
#else
    if (this->m_Opt->rank_cnt == 1) {
      ierr = this->m_kernel.Prolong(xf, xc, fft_coarse->m_FFT->nx, fft_coarse->m_FFT->osize, fft_coarse->m_FFT->ostart); CHKERRQ(ierr);
    } else {
      ierr = ThrowError("Spectral restriction not implemented!"); CHKERRQ(ierr);
    }
#endif

    PetscFunctionReturn(ierr);
}
/********************************************************************
 * @brief Prolong from lower Grid
 *******************************************************************/
PetscErrorCode Spectral::ProlongMerge(ComplexType *xf, const ComplexType *xc, Spectral *fft_coarse) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

#if REG_HAS_CUDA
    this->m_plan->prolongFromMerge(xf, xc, fft_coarse->m_plan);
    //ierr = this->m_kernel.Prolong(xf, xc, fft_coarse->m_FFT->nx, fft_coarse->m_FFT->osize, fft_coarse->m_FFT->ostart); CHKERRQ(ierr);
#else
    if (this->m_Opt->rank_cnt == 1) {
      ierr = this->m_kernel.Prolong(xf, xc, fft_coarse->m_FFT->nx, fft_coarse->m_FFT->osize, fft_coarse->m_FFT->ostart); CHKERRQ(ierr);
    } else {
      ierr = ThrowError("Spectral restriction not implemented!"); CHKERRQ(ierr);
    }
#endif

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
