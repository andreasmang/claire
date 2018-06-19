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

#ifndef _DIFFERENTIATIONSMGPU_HPP_
#define _DIFFERENTIATIONSMGPU_HPP_

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

    virtual PetscErrorCode Gradient(ScalarType*, ScalarType*, ScalarType*, ScalarType*);
    virtual PetscErrorCode Laplacian(ScalarType*, ScalarType*, ScalarType*, ScalarType*);

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();
};



}  // end of namespace



#endif  // _DIFFERENTIATIONSMGPU_HPP_
