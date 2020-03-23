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

#ifndef _DIFFERENTIATION_HPP_
#define _DIFFERENTIATION_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"
#include "interp3.hpp"


namespace reg {

class Differentiation {
 public:
    enum Type {None, Spectral, Finite};
    typedef Differentiation Self;
    Differentiation();
    Differentiation(RegOpt*, Type);
    virtual ~Differentiation();

    virtual PetscErrorCode Gradient(ScalarType*, ScalarType*, ScalarType*, const ScalarType*) = 0;
    PetscErrorCode Gradient(ScalarType**, const ScalarType*);
    PetscErrorCode Gradient(VecField*, const ScalarType*);
    PetscErrorCode Gradient(VecField*, const Vec);
    
    virtual PetscErrorCode Laplacian(ScalarType*, const ScalarType*) = 0;
    virtual PetscErrorCode Laplacian(ScalarType*, ScalarType*, ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*) = 0;
    PetscErrorCode Laplacian(Vec, const Vec);
    PetscErrorCode Laplacian(VecField*, VecField*);
    
    virtual PetscErrorCode Divergence(ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*) = 0;
    PetscErrorCode Divergence(ScalarType*, const ScalarType*const*);
    PetscErrorCode Divergence(ScalarType*, VecField*);
    PetscErrorCode Divergence(Vec, VecField*);
    
    virtual PetscErrorCode RegLapModOp(VecField*, VecField*, ScalarType) = 0;
    virtual PetscErrorCode RegLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0) = 0;
    virtual PetscErrorCode RegBiLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0) = 0;
    virtual PetscErrorCode RegTriLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0) = 0;
    virtual PetscErrorCode InvRegLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0) = 0;
    virtual PetscErrorCode InvRegBiLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0) = 0;
    virtual PetscErrorCode InvRegTriLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0) = 0;
    virtual PetscErrorCode RegTriLapFunc(VecField*, VecField*, ScalarType, ScalarType=0.0) = 0;
    
    virtual PetscErrorCode LerayOperator(VecField*, VecField*, ScalarType, ScalarType) = 0;
        
    const Type m_Type;

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    RegOpt* m_Opt;
};




}  // end of namespace




#endif  // _DIFFERENTIATION_HPP_
