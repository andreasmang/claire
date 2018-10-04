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

#ifndef _DIFFERENTIATIONSM_HPP_
#define _DIFFERENTIATIONSM_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "Differentiation.hpp"

namespace reg {




class DifferentiationSM : public Differentiation {
 public:
    typedef Differentiation SuperClass;
    typedef DifferentiationSM Self;

    DifferentiationSM();
    DifferentiationSM(RegOpt*);
    ~DifferentiationSM();

    virtual PetscErrorCode Gradient(ScalarType*, ScalarType*, ScalarType*, const ScalarType*);
    virtual PetscErrorCode Gradient(ScalarType**, const ScalarType*);
    virtual PetscErrorCode Gradient(VecField*, const ScalarType*);
    virtual PetscErrorCode Laplacian(ScalarType*, const ScalarType*);
    virtual PetscErrorCode Laplacian(ScalarType*, ScalarType*, ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    virtual PetscErrorCode Divergence(ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    virtual PetscErrorCode Divergence(ScalarType*, const ScalarType*const*);
    virtual PetscErrorCode Divergence(ScalarType*, VecField*);
    virtual PetscErrorCode Biharmonic(ScalarType*, ScalarType*, ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    virtual PetscErrorCode Biharmonic(ScalarType*, const ScalarType*);

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();
    std::bitset<3> xyz;
    double timer[NFFTTIMERS];
    int c_grad;
    int c_div;
};




}  // end of namespace

#endif  // _DIFFERENTIATION_HPP_
