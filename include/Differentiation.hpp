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




namespace reg {




class Differentiation {
 public:
    typedef Differentiation Self;
    Differentiation();
    Differentiation(RegOpt*);
    virtual ~Differentiation();

    virtual PetscErrorCode Gradient(ScalarType*, ScalarType*, ScalarType*, ScalarType*) = 0;
    virtual PetscErrorCode Laplacian(ScalarType*, ScalarType*) = 0;
    virtual PetscErrorCode Laplacian(ScalarType*, ScalarType*, ScalarType*, ScalarType*) = 0;
    virtual PetscErrorCode Biharmonic(ScalarType*, ScalarType*, ScalarType*, ScalarType*) = 0;

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    RegOpt* m_Opt;
};




}  // end of namespace




#endif  // _DIFFERENTIATION_HPP_
