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
    
    if (opt->m_Verbosity > 2) {
      DbgMsg("DifferentiationFD created");
    }
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
    
    this->m_tmp = nullptr;
#ifdef REG_HAS_CUDA
    this->mtex = 0;
#endif
    this->m_Ghost = nullptr;
    this->d_Ghost = nullptr;
    //this->m_Work = nullptr;
    this->m_GhostPlan = nullptr;

    ierr = this->SetupData(); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DifferentiationFD::SetupData(ScalarType *x1, ScalarType *x2, ScalarType *x3) {
    PetscErrorCode ierr = 0;
    int isize[3], nx[3];
    ScalarType hx[3];
    PetscFunctionBegin;
     
    for (int i=0; i<3; i++) {
      isize[i] = this->m_Opt->m_Domain.isize[i];
      hx[i] = this->m_Opt->m_Domain.hx[i];
      nx[i] = this->m_Opt->m_Domain.nx[i];
    }

#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)
    if (this->m_Opt->rank_cnt == 1) {
      this->mtex = gpuInitEmptyGradientTexture(nx);
      ierr = initConstants(isize, isize, hx, this->halo); CHKERRQ(ierr);
    } else {
      ierr = AllocateOnce(this->m_GhostPlan, this->m_Opt, this->nghost); CHKERRQ(ierr);
      this->g_alloc_max = this->m_GhostPlan->get_ghost_local_size_x(this->isize_g, this->istart_g);
      this->nlghost = isize_g[0]*isize_g[1]*isize_g[2];

      // ghost data mem alloc on CPU
      //this->m_Ghost = reinterpret_cast<ScalarType*>(accfft_alloc(this->g_alloc_max));
      cudaMalloc((void**)&this->m_Ghost, this->g_alloc_max);
      // work memory on CPU
      //this->m_Work = reinterpret_cast<ScalarType*>(accfft_alloc(sizeof(ScalarType)*this->m_Opt->m_Domain.nl));
      // ghost data mem alloc on GPU
      //cudaMalloc((void**)&this->d_Ghost, this->nlghost*sizeof(ScalarType));
      
      this->mtex = gpuInitEmptyGradientTexture(this->isize_g); CHKERRQ(ierr);
      ierr = initConstants(isize, this->isize_g, hx, this->halo); CHKERRQ(ierr);
    }
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    
    //ierr = AllocateOnce(this->m_tmp, this->m_Opt); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DifferentiationFD::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

#if defined(REG_HAS_CUDA)
    if (this->mtex != 0) {
        cudaDestroyTextureObject(this->mtex);
    }
#endif

    if (this->m_Ghost != nullptr) { 
      //accfft_free(this->m_Ghost);
      cudaFree(this->m_Ghost);
      this->m_Ghost = nullptr;
    }

    /*if (this->m_Work != nullptr) {
      accfft_free(this->m_Work);
      this->m_Work = nullptr;
    }*/

    if (this->d_Ghost != nullptr) {
      cudaFree(this->d_Ghost);
      this->d_Ghost = nullptr;
    }

    if (this->m_GhostPlan != nullptr) {
      ierr = Free(this->m_GhostPlan); CHKERRQ(ierr);
      this->m_GhostPlan = nullptr;
    }
    
    if (this->m_tmp != nullptr) {
      ierr = Free(this->m_tmp); CHKERRQ(ierr);
      this->m_tmp = nullptr;
    }
    
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

#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
    if (this->m_Opt->rank_cnt > 1) {
      //ierr = cudaMemcpy((void*)this->m_Work, (const void*)m, sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost);  CHKERRCUDA(ierr);
      //share_ghost_layer(this->m_Opt, this->nghost, this->isize_g, this->m_Work, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
      //ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
      //ierr = computeGradient(g1, g2, g3, this->d_Ghost, this->mtex, this->m_Opt->m_Domain.isize, true); CHKERRQ(ierr);
      ZeitGeist_define(FD_COMM);
      ZeitGeist_tick(FD_COMM);
      this->m_GhostPlan->share_ghost_x(m, this->m_Ghost); 
      ZeitGeist_tock(FD_COMM);
      ierr = computeGradient(g1, g2, g3, this->m_Ghost, this->mtex, this->m_Opt->m_Domain.isize, true); CHKERRQ(ierr);
    } else {
      ierr = computeGradient(g1, g2, g3, m, this->mtex, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    }
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    ZeitGeist_tock(FD_GRAD);
    DebugGPUStopEvent();
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute laplacian of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Laplacian(ScalarType *l,
                                            const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    DebugGPUStartEvent("FD Laplacian");

    ZeitGeist_define(FD_LAP);
    ZeitGeist_tick(FD_LAP);

#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
    if (this->m_Opt->rank_cnt > 1) {
      //ierr = cudaMemcpy((void*)this->m_Work, (const void*)m, sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
      //share_ghost_layer(this->m_Opt, this->nghost, this->isize_g, this->m_Work, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
      //ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
      //ierr = computeLaplacian(l, this->d_Ghost, this->mtex, this->m_Opt->m_Domain.isize, 1., true); CHKERRQ(ierr);
      ZeitGeist_define(FD_COMM);
      ZeitGeist_tick(FD_COMM);
      this->m_GhostPlan->share_ghost_x(m, this->m_Ghost);
      ZeitGeist_tock(FD_COMM);
      ierr = computeLaplacian(l, this->m_Ghost, this->mtex, this->m_Opt->m_Domain.isize, 1., true); CHKERRQ(ierr);
    } else {
      ierr = computeLaplacian(l, m, this->mtex, this->m_Opt->m_Domain.isize, 1.); CHKERRQ(ierr);
    }
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    ZeitGeist_tock(FD_LAP);
    DebugGPUStopEvent();

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
    const ScalarType *pv[3] = {nullptr, nullptr, nullptr};
    ScalarType *pl[3] = {nullptr, nullptr, nullptr};

    PetscFunctionBegin;
    
    DebugGPUStartEvent("FD Laplacian");
    
    pv[0] = v1;
    pv[1] = v2;
    pv[2] = v3;
    pl[0] = l1;
    pl[1] = l2;
    pl[2] = l3;

    ZeitGeist_define(FD_LAP);
    ZeitGeist_tick(FD_LAP);

    for (int i=0; i<3; i++) {
#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
      if (this->m_Opt->rank_cnt > 1) {
        //ierr = cudaMemcpy((void*)this->m_Work, (const void*)pv[i], sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
        //share_ghost_layer(this->m_Opt, this->nghost, this->isize_g, this->m_Work, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
        //ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
        //ierr = computeLaplacian(pl[i], this->d_Ghost, this->mtex, this->m_Opt->m_Domain.isize, 1., true); CHKERRQ(ierr);
        ierr = DebugNotImplemented(); CHKERRQ(ierr);
      } else {
        ierr = computeLaplacian(pl[i], pv[i], this->mtex, this->m_Opt->m_Domain.isize, 1.); CHKERRQ(ierr);
      }
#else
      ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    }
    ZeitGeist_tock(FD_LAP);
    DebugGPUStopEvent();

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

#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
    if (this->m_Opt->rank_cnt > 1) {
      cudaMemset((void*)l, 0, this->m_Opt->m_Domain.nl*sizeof(ScalarType));
      
      //ierr = cudaMemcpy((void*)this->m_Work, (const void*)v3, sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
      //share_ghost_layer(this->m_Opt, this->nghost, this->isize_g, this->m_Work, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
      //ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
      //ierr  = computeDivergenceZ(l, this->d_Ghost, this->m_Opt->m_Domain.isize, true); CHKERRQ(ierr);
      ZeitGeist_define(FD_COMM);
      ZeitGeist_tick(FD_COMM);
      this->m_GhostPlan->share_ghost_x(v3, this->m_Ghost); 
      ZeitGeist_tock(FD_COMM);
      ierr = computeDivergenceZ(l, this->m_Ghost, this->m_Opt->m_Domain.isize, true); CHKERRQ(ierr);
    
      //ierr = cudaMemcpy((void*)this->m_Work, (const void*)v2, sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
      //share_ghost_layer(this->m_Opt, this->nghost, this->isize_g, this->m_Work, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
      //ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
      //ierr  = computeDivergenceY(l, this->d_Ghost, this->m_Opt->m_Domain.isize, true); CHKERRQ(ierr);
      ZeitGeist_tick(FD_COMM);
      this->m_GhostPlan->share_ghost_x(v2, this->m_Ghost); 
      ZeitGeist_tock(FD_COMM);
      ierr = computeDivergenceY(l, this->m_Ghost, this->m_Opt->m_Domain.isize, true); CHKERRQ(ierr);
      
      //ierr = cudaMemcpy((void*)this->m_Work, (const void*)v1, sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
      //share_ghost_layer(this->m_Opt, this->nghost, this->isize_g, this->m_Work, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
      //ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
      //ierr  = computeDivergenceX(l, this->d_Ghost, this->m_Opt->m_Domain.isize, true); CHKERRQ(ierr);
      ZeitGeist_tick(FD_COMM);
      this->m_GhostPlan->share_ghost_x(v1, this->m_Ghost); 
      ZeitGeist_tock(FD_COMM);
      ierr = computeDivergenceX(l, this->m_Ghost, this->m_Opt->m_Domain.isize, true); CHKERRQ(ierr);
    } else {
      ierr = computeDivergence(l, v1, v2, v3, this->mtex, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    }
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    ZeitGeist_tock(FD_DIV);
    DebugGPUStopEvent();

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
  const ScalarType *pV[3] = {nullptr, nullptr, nullptr};
  ScalarType *pBV[3] = {nullptr, nullptr, nullptr};
  PetscFunctionBegin;
  
  DebugGPUStartEvent("FD Regularization");

  ZeitGeist_define(FD_LAP);
  ZeitGeist_tick(FD_LAP);
  
  ierr = v->GetArraysRead(pV); CHKERRQ(ierr);
  ierr = bv->GetArraysWrite(pBV); CHKERRQ(ierr);
#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
    if (this->m_Opt->rank_cnt > 1) {
      for (int i=0; i<3; i++) {
        //ierr = cudaMemcpy((void*)this->m_Work, (const void*)pV[i], sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
        //share_ghost_layer(this->m_Opt, this->nghost, this->isize_g, this->m_Work, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
        //ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
        ierr = DebugNotImplemented(); CHKERRQ(ierr);
        //ierr = computeLaplacian(pBV[i], this->d_Ghost, this->mtex, this->m_Opt->m_Domain.isize, -b0, true); CHKERRQ(ierr);
      }
    } else {
      for (int i=0; i<3; i++) {
        ierr = computeLaplacian(pBV[i], pV[i], this->mtex, this->m_Opt->m_Domain.isize, -b0); CHKERRQ(ierr);
      }
    }
#else
  ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
  ierr = v->RestoreArrays(); CHKERRQ(ierr);
  ierr = bv->RestoreArrays(); CHKERRQ(ierr);

  ZeitGeist_tock(FD_LAP);
  DebugGPUStopEvent();
  
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
  const ScalarType *pV[3] = {nullptr, nullptr, nullptr};
  ScalarType *pBV[3] = {nullptr, nullptr, nullptr};
  PetscFunctionBegin;

  ZeitGeist_define(FD_INVREG);
  ZeitGeist_tick(FD_INVREG);
  
  ierr = Assert(this->m_tmp != nullptr, "nullptr"); CHKERRQ(ierr);
  
  ScalarType *hx = this->m_Opt->m_Domain.hx;
  ScalarType aii = b0*205./75. * (1./(hx[0]*hx[0]) + 1./(hx[1]*hx[1]) + 1./(hx[2]*hx[2]));
  
  ierr = bv->SetValue(0); CHKERRQ(ierr);
  for (IntType iter=0; iter<2000; iter++) {
    ierr = bv->GetArraysRead(pV); CHKERRQ(ierr);
    ierr = this->m_tmp->GetArraysWrite(pBV); CHKERRQ(ierr);
#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
    if (this->m_Opt->rank_cnt > 1) {
      for (int i=0; i<3; i++) {
        //ierr = cudaMemcpy((void*)this->m_Work, (const void*)pV[i], sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
        //share_ghost_layer(this->m_Opt, this->nghost, this->isize_g, this->m_Work, this->m_Ghost, this->m_GhostWork1, this->m_GhostWork2); 
        //ierr = cudaMemcpy((void*)this->d_Ghost, (const void*)this->m_Ghost, this->nlghost*sizeof(ScalarType), cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
        ierr = DebugNotImplemented(); CHKERRQ(ierr);
        //ierr = computeLaplacian(pBV[i], this->d_Ghost, this->mtex, this->m_Opt->m_Domain.isize, b0, true); CHKERRQ(ierr);
      }
    } else {
      for (int i=0; i<3; i++) {
        ierr = computeLaplacian(pBV[i], pV[i], this->mtex, this->m_Opt->m_Domain.isize, b0); CHKERRQ(ierr);
      }
    }
#else
      ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    ierr = this->m_tmp->RestoreArrays(); CHKERRQ(ierr);
    ierr = bv->RestoreArrays(); CHKERRQ(ierr);
    ierr = this->m_tmp->AXPY(1., v);
  }
    ierr = bv->AXPY(0.1/aii, this->m_tmp);
  ZeitGeist_tock(FD_INVREG);

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
