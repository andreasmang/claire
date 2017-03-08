/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


#ifndef _VECFIELD_H_
#define _VECFIELD_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"




namespace reg {




class VecField {
 public:
    typedef VecField Self;

    VecField();
    VecField(RegOpt*);
    VecField(RegOpt*,int);
    VecField(IntType,IntType);
    ~VecField();

    PetscErrorCode SetOpt(RegOpt*);

    /*! set individual vector components based on a
        flat PETSc vector as an input */
    PetscErrorCode SetComponents(Vec);

    /*! get individual vector components based on a
        flat PETSc vector as an input */
    PetscErrorCode GetComponents(Vec);

    /*! set all components to a given value*/
    PetscErrorCode SetValue(ScalarType);

    /*! pointwise scaling of individual components of
        vector field by scalar */
    PetscErrorCode Scale(ScalarType);

    /*! pointwise scaling of individual components of
        vector field by scalar field; overwrites vector field */
    PetscErrorCode Scale(Vec);

    /*! pointwise scaling of individual components of
        vector field by scalar field; assigned to input
        vector field */
    PetscErrorCode Scale(VecField*,Vec);

    PetscErrorCode Copy(VecField*);

    PetscErrorCode GetArrays(ScalarType*&,ScalarType*&,ScalarType*&);
    PetscErrorCode RestoreArrays(ScalarType*&,ScalarType*&,ScalarType*&);

    PetscErrorCode GetArraysRead(const ScalarType*&,const ScalarType*&,const ScalarType*&);
    PetscErrorCode RestoreArraysRead(const ScalarType*&,const ScalarType*&,const ScalarType*&);

    PetscErrorCode WAXPY(ScalarType,VecField*,VecField*);
    PetscErrorCode AXPY(ScalarType,VecField*);

    /*! compute norm of vector field */
    PetscErrorCode Norm(Vec);
    PetscErrorCode Norm(ScalarType&);

    // individual components
    Vec m_X1;
    Vec m_X2;
    Vec m_X3;
 private:
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);

    PetscErrorCode Allocate(void);
    PetscErrorCode Allocate(IntType,IntType);
    PetscErrorCode Allocate(int);

    RegOpt* m_Opt;
};




}  // end of name space




#endif  // _VECFIELD_H_
