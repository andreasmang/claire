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

#ifndef _DIFFERENTIATIONFD_CPP_
#define _DIFFERENTIATIONFD_CPP_

#include "DifferentiationFD.hpp"

#ifdef REG_HAS_CUDA
#include "TextureDifferentiationKernel.hpp"
#endif



namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DifferentiationFD::DifferentiationFD() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DifferentiationFD::DifferentiationFD(RegOpt* opt) : SuperClass(opt, Type::Finite) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DifferentiationFD::~DifferentiationFD() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DifferentiationFD::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->SetupData(); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DifferentiationFD::SetupData(ScalarType *x1, ScalarType *x2, ScalarType *x3) {
    PetscErrorCode ierr = 0;
    IntType nalloc;
    int isize[3];
    ScalarType hx[3];
    PetscFunctionBegin;
  
    for (int i=0; i<3; i++) {
      isize[i] = this->m_Opt->m_Domain.isize[i];
      hx[i] = this->m_Opt->m_Domain.hx[i];
    }

#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
    this->mtex = gpuInitEmptyGradientTexture(this->m_Opt->m_Domain.nx);
    ierr = initConstants(isize, isize, hx); CHKERRQ(ierr);
#elif defined(REG_HAS_MPICUDA) 
    FFTPlanType *m_plan = this->m_Opt->m_FFT.fft->m_plan;
    this->g_alloc_max = accfft_ghost_xyz_local_size_dft_r2c(m_plan, this->nghost, this->isize_g, this->istart_g);
    this->nlghost = isize_g[0]*isize_g[1]*isize_g[2];
    // work memory allocation for ghost point sharing
    this->m_GhostWork1 = pvfmm::aligned_new<Real> (m_plan->alloc_max + 2 * this->nghost * isize[2] * isize[0]);
    this->m_GhostWork2 = pvfmm::aligned_new<Real> (m_plan->alloc_max + 2 * this->nghost * isize[2] * isize[0] + 2 * this->nghost * isize[2] * this->isize_g[1]);
    // ghost data mem alloc on CPU
    this->m_Ghost = reinterpret_cast<ScalarType*>(accfft_alloc(this->g_alloc_max));
    // ghost data mem alloc on GPU
    cudaMalloc((void**)&this->d_Ghost, this->nlghost*sizeof(ScalarType));
    this->mtex = gpuInitEmptyGradientTexture(this->isize_g); CHKERRQ(ierr);
    ierr = initConstants(isize, this->isize_g, hx); CHKERRQ(ierr);
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DifferentiationFD::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
    if (this->mtex != 0) {
        cudaDestroyTextureObject(this->mtex);
    }
#elif defined(REG_HAS_MPICUDA)
    if (this->mtex != 0) {
        cudaDestroyTextureObject(this->mtex);
    }
    accfft_free(this->m_Ghost);
    cudaFree(this->d_Ghost);
    pvfmm::aligned_delete<Real>(this->m_GhostWork1);
    pvfmm::aligned_delete<Real>(this->m_GhostWork2);
#endif
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Gradient(ScalarType *g1,
                                           ScalarType *g2,
                                           ScalarType *g3,
                                           const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FD Grad");

    ZeitGeist_define(FD_GRAD);
    ZeitGeist_tick(FD_GRAD);
#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
    ierr = computeGradient(g1, g2, g3, m, this->mtex, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
#elif defined(REG_HAS_MPICUDA)
    if (this->m_Opt->m_FFT.fft->m_plan != nullptr) {
      share_ghost_layer(this->m_Opt->m_FFT.fft->m_plan, this->nghost, this->isize_g, (ScalarType*)m, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
    }
    ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
    ierr  = computeGradient(g1, g2, g3, this->d_Ghost, this->mtex, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    ZeitGeist_tock(FD_GRAD);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute laplacian of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Laplacian(ScalarType *l,
                                            const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = DebugNotImplemented(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute laplacian of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Laplacian(ScalarType *l1,
                                            ScalarType *l2,
                                            ScalarType *l3,
                                            const ScalarType *v1,
                                            const ScalarType *v2,
                                            const ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = DebugNotImplemented(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}
 
/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Divergence(ScalarType *l,
                                             const ScalarType *v1,
                                             const ScalarType *v2,
                                             const ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    DebugGPUStartEvent("FD Divergence");
    
    ZeitGeist_define(FD_DIV);
    ZeitGeist_tick(FD_DIV);
#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
    ierr = computeDivergence(l, v1, v2, v3, this->mtex, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
#elif defined(REG_HAS_MPICUDA)
    // Compute divergence in 3 separate function calls to save memory
    // First set the divergence to zero
    cudaMemset((void*)l, 0, this->m_Opt->m_Domain.nl*sizeof(ScalarType));

    share_ghost_layer(this->m_Opt->m_FFT.fft->m_plan, this->nghost, this->isize_g, (ScalarType*)v3, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
    ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
    ierr  = computeDivergenceZ(l, this->d_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  
    share_ghost_layer(this->m_Opt->m_FFT.fft->m_plan, this->nghost, this->isize_g, (ScalarType*)v2, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
    ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
    ierr  = computeDivergenceY(l, this->d_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    
    share_ghost_layer(this->m_Opt->m_FFT.fft->m_plan, this->nghost, this->isize_g, (ScalarType*)v1, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
    ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
    ierr  = computeDivergenceX(l, this->d_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    ZeitGeist_tock(FD_DIV);

    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationFD::RegLapModOp(VecField* bv, VecField* v, ScalarType b0) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::RegLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::RegBiLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::RegTriLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::RegTriLapFunc(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationFD::InvRegLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::InvRegBiLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::InvRegTriLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationFD::LerayOperator(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}




}  // end of namespace




#endif  // _DIFFERENTIATION_CPP_
