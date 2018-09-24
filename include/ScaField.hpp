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

#ifndef _SCAFIELD_HPP_
#define _SCAFIELD_HPP_

#include "CLAIREUtils.hpp"

namespace reg {

class ScaField {
 public:
    typedef ScaField Self;

    ScaField();
    ScaField(Vec, IntType, IntType=1, IntType=1);
    ScaField(IntType,IntType,IntType=1,IntType=1);
    ~ScaField();
    
    operator Vec& ();

    PetscErrorCode GetArray(ScalarType*&, IntType=0, IntType=0, IntType=0);
    PetscErrorCode GetArrayRead(const ScalarType*&, IntType=0, IntType=0, IntType=0);
    PetscErrorCode GetArrayWrite(ScalarType*&, IntType=0, IntType=0, IntType=0);
    PetscErrorCode GetArrayReadWrite(ScalarType*&, IntType=0, IntType=0, IntType=0);
    
    PetscErrorCode RestoreArray();
    
    PetscErrorCode Copy(Vec);
 private:
    typedef enum {None, Read, Write, ReadWrite} AccessType;
    AccessType m_Type;
    
    Vec m_X;
    
    union {
      ScalarType *m_Ptr;
      const ScalarType *m_ConstPtr;
    };
    
    const IntType m_Dim[3];
    const IntType m_Size[2];
    size_t m_Allocated;
    
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);

    PetscErrorCode Allocate(IntType, IntType);
    PetscErrorCode Assign(Vec, IntType);
};




}  // namespace reg




#endif  // _SCAFIELD_HPP_
