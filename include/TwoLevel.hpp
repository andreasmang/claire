#pragma once

#include "RegOpt.hpp"
#include "GhostPlan.hpp"
#include "VecField.hpp"
#include "ScaField.hpp"

namespace reg {

namespace TwoLevelFiniteKernel {
  PetscErrorCode InitConstants(IntType*, IntType*, ScalarType*, IntType*);
  PetscErrorCode Restrict(ScalarType*, const ScalarType*, IntType*);
  PetscErrorCode Prolong(ScalarType*, const ScalarType*, IntType*);
} // namespace TwoLevelFiniteKernel

class TwoLevel {
public:
  TwoLevel(RegOpt* opt) : m_Opt(opt) {}
  virtual ~TwoLevel() {}

  virtual PetscErrorCode Restrict(ScalarType*, const ScalarType*) = 0;
  virtual PetscErrorCode Restrict(ScaField*, ScaField*) = 0;
  virtual PetscErrorCode Restrict(VecField*, VecField*) = 0;
  virtual PetscErrorCode Prolong(ScalarType*, const ScalarType*) = 0;
  virtual PetscErrorCode Prolong(ScaField*, ScaField*) = 0;
  virtual PetscErrorCode Prolong(VecField*, VecField*) = 0;
protected:
  RegOpt *m_Opt;
};

class TwoLevelFFT : public TwoLevel {
public:
  TwoLevelFFT(RegOpt*);
  virtual ~TwoLevelFFT();
  
  virtual PetscErrorCode Restrict(ScalarType*, const ScalarType*);
  virtual PetscErrorCode Restrict(ScaField*, ScaField*);
  virtual PetscErrorCode Restrict(VecField*, VecField*);
  virtual PetscErrorCode Prolong(ScalarType*, const ScalarType*);
  virtual PetscErrorCode Prolong(ScaField*, ScaField*);
  virtual PetscErrorCode Prolong(VecField*, VecField*);
};

class TwoLevelRegFFT : public TwoLevelFFT {
public:
  TwoLevelRegFFT(RegOpt*, ScalarType, ScalarType = 0, bool = false);
  virtual ~TwoLevelRegFFT();
  
  using TwoLevelFFT::Restrict;
  using TwoLevelFFT::Prolong;
  
  virtual PetscErrorCode Restrict(VecField*, VecField*);
  virtual PetscErrorCode Prolong(VecField*, VecField*);
  
  bool restricted;
protected:
  ScalarType beta;
  ScalarType reltol;
  bool sqrtop;
};

class TwoLevelFinite : public TwoLevel {
public:
  TwoLevelFinite(RegOpt*);
  virtual ~TwoLevelFinite();
  
  virtual PetscErrorCode Restrict(ScalarType*, const ScalarType*);
  virtual PetscErrorCode Restrict(ScaField*, ScaField*);
  virtual PetscErrorCode Restrict(VecField*, VecField*);
  virtual PetscErrorCode Prolong(ScalarType*, const ScalarType*);
  virtual PetscErrorCode Prolong(ScaField*, ScaField*);
  virtual PetscErrorCode Prolong(VecField*, VecField*);
protected:
  const IntType nghost = 4;
  IntType halo[3] = {nghost, 0, 0};
  size_t g_alloc_max;
  IntType nlghost, isize_g[3], istart_g[3];
  ScalarType* m_Ghost;
  
  GhostPlan* m_GhostPlan;
};

} // namespace reg
