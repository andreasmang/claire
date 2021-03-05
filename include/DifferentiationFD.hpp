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

#ifndef _DIFFERENTIATIONFD_HPP_
#define _DIFFERENTIATIONFD_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "Differentiation.hpp"
#include "GhostPlan.hpp"


namespace reg {

class DifferentiationFD : public Differentiation {
 public:
    typedef Differentiation SuperClass;
    typedef DifferentiationFD Self;
    
    using SuperClass::Gradient;
    using SuperClass::Divergence;
    using SuperClass::Laplacian;

    DifferentiationFD();
    DifferentiationFD(RegOpt*);
    virtual ~DifferentiationFD();

    virtual PetscErrorCode Gradient(ScalarType*, ScalarType*, ScalarType*, const ScalarType*);
    
    // Laplacian not implemented for FD
    virtual PetscErrorCode Laplacian(ScalarType*, const ScalarType*);
    virtual PetscErrorCode Laplacian(ScalarType*, ScalarType*, ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    
    virtual PetscErrorCode Divergence(ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    
    // Regularization Operators not implemented for FD
    virtual PetscErrorCode RegLapModOp(VecField*, VecField*, ScalarType);
    virtual PetscErrorCode RegLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode RegBiLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode RegTriLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode InvRegLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode InvRegBiLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode InvRegTriLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode RegTriLapFunc(VecField*, VecField*, ScalarType, ScalarType=0.0);
    
    // Leray Operator not implemented for FD
    virtual PetscErrorCode LerayOperator(VecField*, VecField*, ScalarType, ScalarType);
    
 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();
    
    PetscErrorCode SetupData(ScalarType* =nullptr, ScalarType* =nullptr, ScalarType* =nullptr);
    
    VecField *m_tmp;
    
#ifdef REG_HAS_CUDA
  cudaTextureObject_t mtex;
#endif    

  //const int nghost = 4;
  //int halo[3] = {nghost, 0, 0};
  //size_t g_alloc_max;
  //int nlghost, isize_g[3], istart_g[3];
  ScalarType* m_Ghost, *d_Ghost;//, *m_Work;
  IntType isizeg[3];
  IntType halo[3];

  GhostPlan* m_GhostPlan;

};

}  // end of namespace

#endif  // _DIFFERENTIATIONFD_HPP_
