/**
 *  Description: class (container) for vector field (data type
 *  consists of one PETSc vector for each component)
 *  Copyright (c) 2015-2016.
 *  All rights reserved.
 *  This file is part of PGLISTR library.
 *
 *  PGLISTR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  PGLISTR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PGLISTR.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#ifndef _VECFIELD_H_
#define _VECFIELD_H_

#include "RegUtils.h"
#include "RegOpt.h"

namespace reg
{


class VecField
{

public:
    typedef VecField Self;

    VecField();
    VecField(RegOpt*);
    ~VecField();

    /*! set individual vector components based on a
        flat PETSc vector as an input */
    PetscErrorCode SetComponents(Vec);

    /*! get individual vector components based on a
        flat PETSc vector as an input */
    PetscErrorCode GetComponents(Vec);

    PetscErrorCode SetValue(ScalarType);
    PetscErrorCode Scale(ScalarType);
    PetscErrorCode Scale(Vec);
    PetscErrorCode Copy(VecField*);
    PetscErrorCode Scale(VecField*,Vec);

    // individual components
    Vec m_X1;
    Vec m_X2;
    Vec m_X3;

private:

    PetscErrorCode Allocate(void);
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);

    RegOpt* m_Opt;

};

} // end of name space


#endif // _VECFIELD_H_