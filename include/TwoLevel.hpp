#pragma once

#include "RegOpt.hpp"
#include "VecField.hpp"
#include "ScaField.hpp"

namespace reg {

class TwoLevel {
public:
  TwoLevel(RegOpt*);
  virtual ~TwoLevel();

  virtual PetscErrorCode Restrict(ScaField*, ScaField*) = 0;
  virtual PetscErrorCode Restrict(VecField*, VecField*) = 0;
  virtual PetscErrorCode Prolong(ScaField*, ScaField*) = 0;
  virtual PetscErrorCode Prolong(VecField*, VecField*) = 0;
protected:
  RegOpt *m_Opt;
};

class TwoLevelFFT : public TwoLevel {
public:
  TwoLevelFFT(RegOpt*);
  virtual ~TwoLevelFFT();
  
  virtual PetscErrorCode Restrict(ScaField*, ScaField*);
  virtual PetscErrorCode Restrict(VecField*, VecField*);
  virtual PetscErrorCode Prolong(ScaField*, ScaField*);
  virtual PetscErrorCode Prolong(VecField*, VecField*);
};

class TwoLevelRegFFT : public TwoLevelFFT {
public:
  TwoLevelRegFFT(RegOpt*, ScalarType);
  virtual ~TwoLevelRegFFT();
  
  virtual PetscErrorCode Restrict(VecField*, VecField*);
  virtual PetscErrorCode Prolong(VecField*, VecField*);
private:
  ScalarType beta;
};

} // namespace reg
