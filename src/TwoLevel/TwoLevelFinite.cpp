#include "TwoLevel.hpp"

namespace reg {
  
namespace Kernel = TwoLevelFiniteKernel;

TwoLevelFinite::TwoLevelFinite(RegOpt* opt) : TwoLevel(opt) {
    this->m_Ghost = nullptr;
    this->m_GhostPlan = nullptr;

#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)
    if (this->m_Opt->rank_cnt == 1) {
      Kernel::InitConstants(this->m_Opt->m_Domain.isize,
                           this->m_Opt->m_Domain.isize,
                           this->m_Opt->m_Domain.hx,
                           this->halo);
    } else {
      AllocateOnce(this->m_GhostPlan, this->m_Opt, this->nghost);
      this->g_alloc_max = this->m_GhostPlan->get_ghost_local_size_x(this->isize_g, this->istart_g);
      this->nlghost = isize_g[0]*isize_g[1]*isize_g[2];

      cudaMalloc((void**)&this->m_Ghost, this->g_alloc_max);
      
      Kernel::InitConstants(this->m_Opt->m_Domain.isize,
                           this->isize_g,
                           this->m_Opt->m_Domain.hx,
                           this->halo);
    }
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
}

TwoLevelFinite::~TwoLevelFinite() {
    if (this->m_Ghost != nullptr) { 
      cudaFree(this->m_Ghost);
      this->m_Ghost = nullptr;
    }

    if (this->m_GhostPlan != nullptr) {
      Free(this->m_GhostPlan);
      this->m_GhostPlan = nullptr;
    }
}

PetscErrorCode TwoLevelFinite::Restrict(ScalarType* dst, const ScalarType* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  DebugGPUStartEvent("Finite Restrict");

  ZeitGeist_define(FD_RESTRICT);
  ZeitGeist_tick(FD_RESTRICT);

#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
  if (this->m_Opt->rank_cnt > 1) {
    ZeitGeist_define(FD_COMM);
    ZeitGeist_tick(FD_COMM);
    this->m_GhostPlan->share_ghost_x(src, this->m_Ghost); 
    ZeitGeist_tock(FD_COMM);
    
    ierr = Kernel::Restrict(dst, this->m_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  } else {
    ierr = Kernel::Restrict(dst, src, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  }
#else
  ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif

  ZeitGeist_tock(FD_RESTRICT);
  DebugGPUStopEvent();
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelFinite::Restrict(ScaField* dst, ScaField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  const ScalarType *pVf = nullptr;
  ScalarType *pVc = nullptr;
  
  ierr = src->GetArrayRead(pVf); CHKERRQ(ierr);
  ierr = dst->GetArrayWrite(pVc); CHKERRQ(ierr);
  
  ierr = this->Restrict(pVc, pVf);
  
  ierr = src->RestoreArray(); CHKERRQ(ierr);
  ierr = dst->RestoreArray(); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelFinite::Restrict(VecField* dst, VecField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  DebugGPUStartEvent("Finite Restrict");

  ZeitGeist_define(FD_RESTRICT);
  ZeitGeist_tick(FD_RESTRICT);
  
  const ScalarType* pSrc[3] = {nullptr};
  ScalarType* pDst[3] = {nullptr};
  
  src->GetArraysRead(pSrc);
  dst->GetArraysWrite(pDst);

#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
  if (this->m_Opt->rank_cnt > 1) {
    ZeitGeist_define(FD_COMM);
    ZeitGeist_tick(FD_COMM);
    this->m_GhostPlan->share_ghost_x(pSrc[0], this->m_Ghost); 
    ZeitGeist_tock(FD_COMM);
    
    ierr = Kernel::Restrict(pDst[0], this->m_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    
    ZeitGeist_tick(FD_COMM);
    this->m_GhostPlan->share_ghost_x(pSrc[1], this->m_Ghost); 
    ZeitGeist_tock(FD_COMM);
    
    ierr = Kernel::Restrict(pDst[1], this->m_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    
    ZeitGeist_tick(FD_COMM);
    this->m_GhostPlan->share_ghost_x(pSrc[2], this->m_Ghost); 
    ZeitGeist_tock(FD_COMM);
    
    ierr = Kernel::Restrict(pDst[2], this->m_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  } else {
    ierr = Kernel::Restrict(pDst[0], pSrc[0], this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    ierr = Kernel::Restrict(pDst[1], pSrc[1], this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    ierr = Kernel::Restrict(pDst[2], pSrc[2], this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  }
#else
  ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif

  src->RestoreArrays();
  dst->RestoreArrays();

  ZeitGeist_tock(FD_RESTRICT);
  ZeitGeist_inc(FD_RESTRICT);
  ZeitGeist_inc(FD_RESTRICT);
  DebugGPUStopEvent();
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelFinite::Prolong(ScalarType* dst, const ScalarType* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  DebugGPUStartEvent("Finite Prolong");

  ZeitGeist_define(FD_PROLONG);
  ZeitGeist_tick(FD_PROLONG);

#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
  if (this->m_Opt->rank_cnt > 1) {
    ZeitGeist_define(FD_COMM);
    ZeitGeist_tick(FD_COMM);
    this->m_GhostPlan->share_ghost_x(src, this->m_Ghost); 
    ZeitGeist_tock(FD_COMM);
    
    ierr = Kernel::Prolong(dst, this->m_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  } else {
    ierr = Kernel::Prolong(dst, src, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  }
#else
  ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif

  ZeitGeist_tock(FD_PROLONG);
  DebugGPUStopEvent();
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelFinite::Prolong(ScaField* dst, ScaField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  const ScalarType *pVc = nullptr;
  ScalarType *pVf = nullptr;
  
  ierr = src->GetArrayRead(pVc); CHKERRQ(ierr);
  ierr = dst->GetArrayWrite(pVf); CHKERRQ(ierr);
  
  ierr = this->Prolong(pVf, pVc);
  
  ierr = src->RestoreArray(); CHKERRQ(ierr);
  ierr = dst->RestoreArray(); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelFinite::Prolong(VecField* dst, VecField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  DebugGPUStartEvent("Finite Prolong");

  ZeitGeist_define(FD_PROLONG);
  ZeitGeist_tick(FD_PROLONG);
  
  const ScalarType* pSrc[3] = {nullptr};
  ScalarType* pDst[3] = {nullptr};
  
  src->GetArraysRead(pSrc);
  dst->GetArraysWrite(pDst);

#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
  if (this->m_Opt->rank_cnt > 1) {
    ZeitGeist_define(FD_COMM);
    ZeitGeist_tick(FD_COMM);
    this->m_GhostPlan->share_ghost_x(pSrc[0], this->m_Ghost); 
    ZeitGeist_tock(FD_COMM);
    
    ierr = Kernel::Prolong(pDst[0], this->m_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    
    ZeitGeist_tick(FD_COMM);
    this->m_GhostPlan->share_ghost_x(pSrc[1], this->m_Ghost); 
    ZeitGeist_tock(FD_COMM);
    
    ierr = Kernel::Prolong(pDst[1], this->m_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    
    ZeitGeist_tick(FD_COMM);
    this->m_GhostPlan->share_ghost_x(pSrc[2], this->m_Ghost); 
    ZeitGeist_tock(FD_COMM);
    
    ierr = Kernel::Prolong(pDst[2], this->m_Ghost, this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  } else {
    ierr = Kernel::Prolong(pDst[0], pSrc[0], this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    ierr = Kernel::Prolong(pDst[1], pSrc[1], this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
    ierr = Kernel::Prolong(pDst[2], pSrc[2], this->m_Opt->m_Domain.isize); CHKERRQ(ierr);
  }
#else
  ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif

  src->RestoreArrays();
  dst->RestoreArrays();

  ZeitGeist_tock(FD_PROLONG);
  ZeitGeist_inc(FD_PROLONG);
  ZeitGeist_inc(FD_PROLONG);
  DebugGPUStopEvent();
  
  
  PetscFunctionReturn(ierr);
}
  
} // namespace reg
