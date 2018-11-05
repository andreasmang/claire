/*************************************************************************
 *  Copyright (c) 2016.
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



#ifndef _TENFIELD_H_
#define _TENFIELD_H_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"

namespace reg {


class TenField {

public:
    typedef TenField Self;

    TenField();
    TenField(RegOpt*);
    TenField(RegOpt*,int);
    TenField(IntType,IntType);
    virtual ~TenField();

    PetscErrorCode Copy(TenField*);
    PetscErrorCode SetValue(ScalarType);
    PetscErrorCode SetComponents(Vec);
    PetscErrorCode GetComponents(Vec);
    PetscErrorCode Scale(ScalarType);
    PetscErrorCode Scale(Vec);
    PetscErrorCode SetIdentity();

    PetscErrorCode GetArrays(ScalarType*&,ScalarType*&,ScalarType*&,
                             ScalarType*&,ScalarType*&,ScalarType*&,
                             ScalarType*&,ScalarType*&,ScalarType*&);
    PetscErrorCode RestoreArrays(ScalarType*&,ScalarType*&,ScalarType*&,
                                 ScalarType*&,ScalarType*&,ScalarType*&,
                                 ScalarType*&,ScalarType*&,ScalarType*&);
    PetscErrorCode GetArraysRead(const ScalarType*&,const ScalarType*&,const ScalarType*&,
                                 const ScalarType*&,const ScalarType*&,const ScalarType*&,
                                 const ScalarType*&,const ScalarType*&,const ScalarType*&);
    PetscErrorCode RestoreArraysRead(const ScalarType*&,const ScalarType*&,const ScalarType*&,
                                     const ScalarType*&,const ScalarType*&,const ScalarType*&,
                                     const ScalarType*&,const ScalarType*&,const ScalarType*&);

    PetscErrorCode AXPY(ScalarType,TenField*);
    PetscErrorCode WAXPY(ScalarType,TenField*,TenField*);

    // individual components
    Vec m_X11;
    Vec m_X12;
    Vec m_X13;
    Vec m_X21;
    Vec m_X22;
    Vec m_X23;
    Vec m_X31;
    Vec m_X32;
    Vec m_X33;

private:

    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode SetOpt(RegOpt*);

    PetscErrorCode Allocate(void);
    PetscErrorCode Allocate(IntType,IntType);
    PetscErrorCode Allocate(int);

    RegOpt* m_Opt;

};



} // end of name space


#endif // _TENFIELD_H_
