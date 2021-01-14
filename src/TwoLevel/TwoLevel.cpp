#include "TwoLevel.hpp"
#include "DifferentiationKernel.hpp"

namespace reg {
  
TwoLevel::TwoLevel(RegOpt* opt) : m_Opt(opt) {};
TwoLevel::~TwoLevel() {};

//--------------------------------------------------------------------------------------------------

TwoLevelFFT::TwoLevelFFT(RegOpt* opt) : TwoLevel(opt) {};
TwoLevelFFT::~TwoLevelFFT() {};

PetscErrorCode TwoLevelFFT::Restrict(ScaField* dst, ScaField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_RESTRICT);
  ZeitGeist_tick(FFT_RESTRICT);
  
  FourierTransform *fine   = &this->m_Opt->m_FFT;
  FourierTransform *coarse = &this->m_Opt->m_FFT_coarse;

  ScalarType scale = 1./(fine->nx[0]*fine->nx[1]*fine->nx[2]);
  
  const ScalarType *pVf = nullptr;
  ScalarType *pVc = nullptr;
  
  ComplexType *pXHat_c[1], *pXHat_f[1];
  
  pXHat_f[0] =   &fine->fft->m_WorkSpace[0];  
  pXHat_c[0] = &coarse->fft->m_WorkSpace[0];
  
  ierr = Assert(pXHat_f[0] != pXHat_c[0], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  
  
  ierr = src->GetArrayRead(pVf); CHKERRQ(ierr);
  fine->fft->FFT_R2C(pVf, pXHat_f[0]);
  ierr = src->RestoreArray(); CHKERRQ(ierr);
    
  fine->fft->Restrict(pXHat_c[0], pXHat_f[0], coarse->fft);
  
  coarse->fft->Scale(pXHat_c[0], scale);
  
  ierr = dst->GetArrayWrite(pVc); CHKERRQ(ierr);
  coarse->fft->FFT_C2R(pXHat_c[0], pVc);
  ierr = dst->RestoreArray(); CHKERRQ(ierr);
  
  ZeitGeist_tock(FFT_RESTRICT);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelFFT::Restrict(VecField* dst, VecField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_RESTRICT);
  ZeitGeist_tick(FFT_RESTRICT);
  
  FourierTransform *fine   = &this->m_Opt->m_FFT;
  FourierTransform *coarse = &this->m_Opt->m_FFT_coarse;

  ScalarType scale = 1./(fine->nx[0]*fine->nx[1]*fine->nx[2]);
  
  const ScalarType *pVf[3] = {nullptr, nullptr, nullptr};
  ScalarType *pVc[3] = {nullptr, nullptr, nullptr};
  
  ComplexType *pXHat_c[3], *pXHat_f[3];
  
  IntType nc_f =   fine->nalloc/(2*sizeof(ScalarType));
  IntType nc_c = coarse->nalloc/(2*sizeof(ScalarType));
  
  pXHat_f[0] = &fine->fft->m_WorkSpace[0*nc_f];
  pXHat_f[1] = &fine->fft->m_WorkSpace[1*nc_f];
  pXHat_f[2] = &fine->fft->m_WorkSpace[2*nc_f];
  
  pXHat_c[0] = &coarse->fft->m_WorkSpace[0*nc_c];
  pXHat_c[1] = &coarse->fft->m_WorkSpace[1*nc_c];
  pXHat_c[2] = &coarse->fft->m_WorkSpace[2*nc_c];
  
  ierr = Assert(pXHat_f[0] != pXHat_c[0], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  ierr = Assert(pXHat_f[1] != pXHat_c[1], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  ierr = Assert(pXHat_f[2] != pXHat_c[2], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  
  ierr = src->GetArraysRead(pVf); CHKERRQ(ierr);
  fine->fft->FFT_R2C(pVf[0], pXHat_f[0]);
  fine->fft->FFT_R2C(pVf[1], pXHat_f[1]);
  fine->fft->FFT_R2C(pVf[2], pXHat_f[2]);
  ierr = src->RestoreArrays(); CHKERRQ(ierr);
  
  fine->fft->Restrict(pXHat_c[0], pXHat_f[0], coarse->fft);
  fine->fft->Restrict(pXHat_c[1], pXHat_f[1], coarse->fft);
  fine->fft->Restrict(pXHat_c[2], pXHat_f[2], coarse->fft);
  
  coarse->fft->Scale(pXHat_c[0], scale);
  coarse->fft->Scale(pXHat_c[1], scale);
  coarse->fft->Scale(pXHat_c[2], scale);
  
  ierr = dst->GetArraysWrite(pVc); CHKERRQ(ierr);
  coarse->fft->FFT_C2R(pXHat_c[0], pVc[0]);
  coarse->fft->FFT_C2R(pXHat_c[1], pVc[1]);
  coarse->fft->FFT_C2R(pXHat_c[2], pVc[2]);
  ierr = dst->RestoreArrays(); CHKERRQ(ierr);
  
  ZeitGeist_tock(FFT_RESTRICT);
  ZeitGeist_inc(FFT_RESTRICT);
  ZeitGeist_inc(FFT_RESTRICT);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelFFT::Prolong(ScaField* dst, ScaField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_PROLONG);
  ZeitGeist_tick(FFT_PROLONG);
  
  FourierTransform *fine   = &this->m_Opt->m_FFT;
  FourierTransform *coarse = &this->m_Opt->m_FFT_coarse;
  
  ScalarType scale = 1./(coarse->nx[0]*coarse->nx[1]*coarse->nx[2]);
  
  ScalarType *pVf = nullptr;
  const ScalarType *pVc = nullptr;
  
  ComplexType *pXHat_c[1], *pXHat_f[1];
  
  pXHat_f[0]   = &fine->fft->m_WorkSpace[0];
  pXHat_c[0] = &coarse->fft->m_WorkSpace[0];
  
  ierr = Assert(pXHat_f[0] != pXHat_c[0], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  
  ierr = src->GetArrayRead(pVc); CHKERRQ(ierr);
  coarse->fft->FFT_R2C(pVc, pXHat_c[0]);
  ierr = src->RestoreArray(); CHKERRQ(ierr);
  
  coarse->fft->Scale(pXHat_c[0], scale);
  
  fine->fft->Prolong(pXHat_f[0], pXHat_c[0], coarse->fft);
  
  ierr = dst->GetArrayWrite(pVf); CHKERRQ(ierr);
  fine->fft->FFT_C2R(pXHat_f[0], pVf);
  ierr = dst->RestoreArray(); CHKERRQ(ierr);
  
  ZeitGeist_tock(FFT_PROLONG);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelFFT::Prolong(VecField* dst, VecField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_PROLONG);
  ZeitGeist_tick(FFT_PROLONG);
  
  FourierTransform *fine   = &this->m_Opt->m_FFT;
  FourierTransform *coarse = &this->m_Opt->m_FFT_coarse;
  
  ScalarType scale = 1./(coarse->nx[0]*coarse->nx[1]*coarse->nx[2]);
  
  IntType nc_f =   fine->nalloc/(2*sizeof(ScalarType));
  IntType nc_c = coarse->nalloc/(2*sizeof(ScalarType));
  
  ScalarType *pVf[3] = {nullptr, nullptr, nullptr};
  const ScalarType *pVc[3] = {nullptr, nullptr, nullptr};
  
  ComplexType *pXHat_c[3], *pXHat_f[3];
  
  pXHat_f[0] = &fine->fft->m_WorkSpace[0*nc_f];
  pXHat_f[1] = &fine->fft->m_WorkSpace[1*nc_f];
  pXHat_f[2] = &fine->fft->m_WorkSpace[2*nc_f];
  
  pXHat_c[0] = &coarse->fft->m_WorkSpace[0*nc_c];
  pXHat_c[1] = &coarse->fft->m_WorkSpace[1*nc_c];
  pXHat_c[2] = &coarse->fft->m_WorkSpace[2*nc_c];
  
  ierr = Assert(pXHat_f[0] != pXHat_c[0], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  ierr = Assert(pXHat_f[1] != pXHat_c[1], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  ierr = Assert(pXHat_f[2] != pXHat_c[2], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  
  ierr = src->GetArraysRead(pVc); CHKERRQ(ierr);
  coarse->fft->FFT_R2C(pVc[0], pXHat_c[0]);
  coarse->fft->FFT_R2C(pVc[1], pXHat_c[1]);
  coarse->fft->FFT_R2C(pVc[2], pXHat_c[2]);
  ierr = src->RestoreArrays(); CHKERRQ(ierr);
  
  coarse->fft->Scale(pXHat_c[0], scale);
  coarse->fft->Scale(pXHat_c[1], scale);
  coarse->fft->Scale(pXHat_c[2], scale);
  
  fine->fft->Prolong(pXHat_f[0], pXHat_c[0], coarse->fft);
  fine->fft->Prolong(pXHat_f[1], pXHat_c[1], coarse->fft);
  fine->fft->Prolong(pXHat_f[2], pXHat_c[2], coarse->fft);
  
  ierr = dst->GetArraysWrite(pVf); CHKERRQ(ierr);
  fine->fft->FFT_C2R(pXHat_f[0], pVf[0]);
  fine->fft->FFT_C2R(pXHat_f[1], pVf[1]);
  fine->fft->FFT_C2R(pXHat_f[2], pVf[2]);
  ierr = dst->RestoreArrays(); CHKERRQ(ierr);
  
  ZeitGeist_tock(FFT_PROLONG);
  ZeitGeist_inc(FFT_PROLONG);
  ZeitGeist_inc(FFT_PROLONG);
  
  PetscFunctionReturn(ierr);
}

//--------------------------------------------------------------------------------------------------

TwoLevelRegFFT::TwoLevelRegFFT(RegOpt* opt, ScalarType beta) : TwoLevelFFT(opt), beta(beta) {};
TwoLevelRegFFT::~TwoLevelRegFFT() {};

PetscErrorCode TwoLevelRegFFT::Restrict(VecField* dst, VecField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_RESTRICT);
  ZeitGeist_tick(FFT_RESTRICT);
  
  DifferentiationKernel kernel;
  
  FourierTransform *fine   = &this->m_Opt->m_FFT;
  FourierTransform *coarse = &this->m_Opt->m_FFT_coarse;

  ScalarType scale = 1./(fine->nx[0]*fine->nx[1]*fine->nx[2]);
  
  const ScalarType *pVf[3] = {nullptr, nullptr, nullptr};
  ScalarType *pVc[3] = {nullptr, nullptr, nullptr};
  
  ComplexType *pXHat_c[3], *pXHat_f[3];
  
  IntType nc_f =   fine->nalloc/(2*sizeof(ScalarType));
  IntType nc_c = coarse->nalloc/(2*sizeof(ScalarType));
  
  pXHat_f[0] = &fine->fft->m_WorkSpace[0*nc_f];
  pXHat_f[1] = &fine->fft->m_WorkSpace[1*nc_f];
  pXHat_f[2] = &fine->fft->m_WorkSpace[2*nc_f];
  
  pXHat_c[0] = &coarse->fft->m_WorkSpace[0*nc_c];
  pXHat_c[1] = &coarse->fft->m_WorkSpace[1*nc_c];
  pXHat_c[2] = &coarse->fft->m_WorkSpace[2*nc_c];
  
  kernel.pXHat[0] = pXHat_f[0];
  kernel.pXHat[1] = pXHat_f[1];
  kernel.pXHat[2] = pXHat_f[2];
  kernel.nx[0] = fine->nx[0];
  kernel.nx[1] = fine->nx[1];
  kernel.nx[2] = fine->nx[2];
  kernel.nl[0] = fine->osize[0];
  kernel.nl[1] = fine->osize[1];
  kernel.nl[2] = fine->osize[2];
  kernel.nstart[0] = fine->ostart[0];
  kernel.nstart[1] = fine->ostart[1];
  kernel.nstart[2] = fine->ostart[2];
  kernel.scale = scale;
  
  ierr = Assert(pXHat_f[0] != pXHat_c[0], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  ierr = Assert(pXHat_f[1] != pXHat_c[1], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  ierr = Assert(pXHat_f[2] != pXHat_c[2], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  
  ierr = src->GetArraysRead(pVf); CHKERRQ(ierr);
  fine->fft->FFT_R2C(pVf[0], pXHat_f[0]);
  fine->fft->FFT_R2C(pVf[1], pXHat_f[1]);
  fine->fft->FFT_R2C(pVf[2], pXHat_f[2]);
  ierr = src->RestoreArrays(); CHKERRQ(ierr);
  
  ierr = kernel.InverseLaplacian(false, beta, 0); CHKERRQ(ierr);
  
  fine->fft->Restrict(pXHat_c[0], pXHat_f[0], coarse->fft);
  fine->fft->Restrict(pXHat_c[1], pXHat_f[1], coarse->fft);
  fine->fft->Restrict(pXHat_c[2], pXHat_f[2], coarse->fft);
  
  ierr = dst->GetArraysWrite(pVc); CHKERRQ(ierr);
  coarse->fft->FFT_C2R(pXHat_c[0], pVc[0]);
  coarse->fft->FFT_C2R(pXHat_c[1], pVc[1]);
  coarse->fft->FFT_C2R(pXHat_c[2], pVc[2]);
  ierr = dst->RestoreArrays(); CHKERRQ(ierr);
  
  ZeitGeist_tock(FFT_RESTRICT);
  ZeitGeist_inc(FFT_RESTRICT);
  ZeitGeist_inc(FFT_RESTRICT);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode TwoLevelRegFFT::Prolong(VecField* dst, VecField* src) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ZeitGeist_define(FFT_PROLONG);
  ZeitGeist_tick(FFT_PROLONG);
  
  FourierTransform *fine   = &this->m_Opt->m_FFT;
  FourierTransform *coarse = &this->m_Opt->m_FFT_coarse;
  
  ScalarType scale = 1./(coarse->nx[0]*coarse->nx[1]*coarse->nx[2]);
  
  IntType nc_f =   fine->nalloc/(2*sizeof(ScalarType));
  IntType nc_c = coarse->nalloc/(2*sizeof(ScalarType));
  
  ScalarType *pVf[3] = {nullptr, nullptr, nullptr};
  const ScalarType *pVc[3] = {nullptr, nullptr, nullptr};
  
  ComplexType *pXHat_c[3], *pXHat_f[3];
  
  pXHat_f[0] = &fine->fft->m_WorkSpace[0*nc_f];
  pXHat_f[1] = &fine->fft->m_WorkSpace[1*nc_f];
  pXHat_f[2] = &fine->fft->m_WorkSpace[2*nc_f];
  
  pXHat_c[0] = &coarse->fft->m_WorkSpace[0*nc_c];
  pXHat_c[1] = &coarse->fft->m_WorkSpace[1*nc_c];
  pXHat_c[2] = &coarse->fft->m_WorkSpace[2*nc_c];
  
  ierr = Assert(pXHat_f[0] != pXHat_c[0], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  ierr = Assert(pXHat_f[1] != pXHat_c[1], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  ierr = Assert(pXHat_f[2] != pXHat_c[2], "using same memory for fine and coarse not supported"); CHKERRQ(ierr);
  
  ierr = src->GetArraysRead(pVc); CHKERRQ(ierr);
  coarse->fft->FFT_R2C(pVc[0], pXHat_c[0]);
  coarse->fft->FFT_R2C(pVc[1], pXHat_c[1]);
  coarse->fft->FFT_R2C(pVc[2], pXHat_c[2]);
  ierr = src->RestoreArrays(); CHKERRQ(ierr);
  
  coarse->fft->Scale(pXHat_c[0], scale);
  coarse->fft->Scale(pXHat_c[1], scale);
  coarse->fft->Scale(pXHat_c[2], scale);
  
  fine->fft->ProlongMerge(pXHat_f[0], pXHat_c[0], coarse->fft);
  fine->fft->ProlongMerge(pXHat_f[1], pXHat_c[1], coarse->fft);
  fine->fft->ProlongMerge(pXHat_f[2], pXHat_c[2], coarse->fft);
  
  ierr = dst->GetArraysWrite(pVf); CHKERRQ(ierr);
  fine->fft->FFT_C2R(pXHat_f[0], pVf[0]);
  fine->fft->FFT_C2R(pXHat_f[1], pVf[1]);
  fine->fft->FFT_C2R(pXHat_f[2], pVf[2]);
  ierr = dst->RestoreArrays(); CHKERRQ(ierr);
  
  ZeitGeist_tock(FFT_PROLONG);
  ZeitGeist_inc(FFT_PROLONG);
  ZeitGeist_inc(FFT_PROLONG);
  
  PetscFunctionReturn(ierr);
}
  
} // namespace reg
