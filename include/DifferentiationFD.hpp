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


namespace reg {

class DifferentiationFD : public Differentiation {
 public:
    typedef Differentiation SuperClass;
    typedef DifferentiationFD Self;
    DifferentiationFD();
    DifferentiationFD(RegOpt*);
    virtual ~DifferentiationFD();

    virtual PetscErrorCode Gradient(ScalarType*, ScalarType*, ScalarType*, const ScalarType*);
    virtual PetscErrorCode Gradient(ScalarType**, const ScalarType*);
    virtual PetscErrorCode Gradient(VecField*, const ScalarType*);
    virtual PetscErrorCode Gradient(VecField*, const Vec);
    
    // Laplacian not implemented for FD
    virtual PetscErrorCode Laplacian(ScalarType*, const ScalarType*);
    virtual PetscErrorCode Laplacian(ScalarType*, ScalarType*, ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    virtual PetscErrorCode Laplacian(Vec, const Vec);
    virtual PetscErrorCode Laplacian(VecField*, VecField*);
    
    virtual PetscErrorCode Divergence(ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    virtual PetscErrorCode Divergence(ScalarType*, const ScalarType*const*);
    virtual PetscErrorCode Divergence(ScalarType*, VecField*);
    virtual PetscErrorCode Divergence(Vec, VecField*);
    
    // Regularization Operators not implemented for FD
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
    
    cudaTextureObject_t mtex;
    
    //IntType nx[3];
    //ScalarType* m_grad[3];
    //ScalarType* m;

};

}  // end of namespace

#endif  // _DIFFERENTIATIONFD_HPP_
