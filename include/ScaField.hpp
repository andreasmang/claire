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
#include "RegOpt.hpp"

namespace reg {
  
class ScaField {
 public:
    typedef ScaField Self;

    ScaField(RegOpt*, bool=false, bool=false);
    ScaField(RegOpt*, ScalarType, bool=false, bool=false);
    ScaField(RegOpt*, Vec, bool=false, bool=false);
    virtual ~ScaField();
    
    PetscErrorCode SetVector(Vec);
    PetscErrorCode SetSize(IntType, IntType, IntType=1, IntType=1);
    
    operator Vec& ();
    
    PetscErrorCode Set(ScalarType);

    PetscErrorCode GetArray(ScalarType*&, IntType=0, IntType=0, IntType=0);
    PetscErrorCode GetArrayRead(const ScalarType*&, IntType=0, IntType=0, IntType=0);
    PetscErrorCode GetArrayWrite(ScalarType*&, IntType=0, IntType=0, IntType=0);
    PetscErrorCode GetArrayReadWrite(ScalarType*&, IntType=0, IntType=0, IntType=0);
    
    PetscErrorCode RestoreArray();
    
    PetscErrorCode DebugInfo(std::string, int, const char*);
    
    PetscErrorCode SetFrame(Vec, IntType);
    PetscErrorCode GetFrame(Vec, IntType);
    
    PetscErrorCode Copy(Vec);
    PetscErrorCode Copy(ScaField*);
    PetscErrorCode CopyFrame(IntType);
    
    size_t GetSize() const;
    
    Vec m_X;
 private:
    typedef enum {None, Read, Write, ReadWrite} AccessType;
    AccessType m_Type;
    
    union {
      ScalarType *m_Ptr;
      const ScalarType *m_ConstPtr;
    };
    
    IntType m_Dim[3];
    IntType m_Size[3];
    size_t m_Allocated;
    
    RegOpt *m_Opt;
    
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);

    PetscErrorCode Allocate(IntType, IntType);
    PetscErrorCode Assign(Vec);
};

inline PetscErrorCode GetRawPointer(ScaField* vec, ScalarType** ptr) {
  return vec->GetArray(*ptr);
}
inline PetscErrorCode GetRawPointerRead(ScaField* vec, const ScalarType** ptr) {
  return vec->GetArrayRead(*ptr);
}
inline PetscErrorCode GetRawPointerReadWrite(ScaField* vec, ScalarType** ptr) {
  return vec->GetArrayReadWrite(*ptr);
}
inline PetscErrorCode GetRawPointerWrite(ScaField* vec, ScalarType** ptr) {
  return vec->GetArrayWrite(*ptr);
}

inline PetscErrorCode RestoreRawPointer(ScaField* vec, ScalarType** ptr) {
  return vec->RestoreArray();
}
inline PetscErrorCode RestoreRawPointerRead(ScaField* vec, const ScalarType** ptr) {
  return vec->RestoreArray();
}
inline PetscErrorCode RestoreRawPointerReadWrite(ScaField* vec, ScalarType** ptr) {
  return vec->RestoreArray();
}
inline PetscErrorCode RestoreRawPointerWrite(ScaField* vec, ScalarType** ptr) {
  return vec->RestoreArray();
}


}  // namespace reg




#endif  // _SCAFIELD_HPP_
