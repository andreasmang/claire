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


#ifndef _VECFIELD_H_
#define _VECFIELD_H_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"




namespace reg {




class VecField {
 public:
    typedef VecField Self;

    VecField();
    VecField(RegOpt*);
    VecField(RegOpt*, Vec);
    VecField(RegOpt*, Vec, Vec, Vec);
    VecField(RegOpt*,int);
    VecField(IntType,IntType);
    virtual ~VecField();

    PetscErrorCode SetOpt(RegOpt*);
    
    PetscErrorCode SetVector(Vec);

    /*! set individual vector components based on a
        flat PETSc vector as an input */
    PetscErrorCode SetComponents(Vec, std::string format="block");
    PetscErrorCode SetComponents(const ScalarType*, std::string format="block");
    PetscErrorCode SetComponents(const ScalarType*, const ScalarType*, const ScalarType*);
    
    //PetscErrorCode SetCopyComponents(Vec, std::string format="block");
    //PetscErrorCode SetCopyComponents(const ScalarType*, std::string format="block");
    //PetscErrorCode SetCopyComponents(const ScalarType*, const ScalarType*, const ScalarType*);

    /*! get individual vector components based on a
        flat PETSc vector as an input */
    PetscErrorCode GetComponents(Vec, std::string format="block");
    PetscErrorCode GetComponents(ScalarType*, std::string format="block");
    PetscErrorCode GetComponents(ScalarType*, ScalarType*, ScalarType*);
    
    //PetscErrorCode GetCopyComponents(Vec, std::string format="block");
    //PetscErrorCode GetCopyComponents(ScalarType*, std::string format="block");
    //PetscErrorCode GetCopyComponents(ScalarType*, ScalarType*, ScalarType*);

    /*! set all components to a given value*/
    PetscErrorCode SetValue(ScalarType);

    /*! pointwise scaling of individual components of
        vector field by scalar */
    PetscErrorCode Scale(ScalarType);
    /*! pointwise scaling of individual components of
        vector field by individual scalars */
    PetscErrorCode Scale(ScalarType*);

    /*! pointwise scaling of individual components of
        vector field by scalar field; overwrites vector field */
    PetscErrorCode Scale(Vec);

    /*! pointwise scaling of individual components of
        vector field by scalar field; assigned to input vector field */
    PetscErrorCode Scale(VecField*,Vec);
    PetscErrorCode Copy(VecField*);

    PetscErrorCode GetSize(IntType&, IntType&);

    PetscErrorCode GetArrays(ScalarType*&, ScalarType*&, ScalarType*&);
    PetscErrorCode GetArrays(ScalarType**);
    PetscErrorCode GetArraysRead(const ScalarType*&, const ScalarType*&, const ScalarType*&);
    PetscErrorCode GetArraysRead(const ScalarType**);
    PetscErrorCode GetArraysReadWrite(ScalarType*&, ScalarType*&, ScalarType*&);
    PetscErrorCode GetArraysReadWrite(ScalarType**);
    PetscErrorCode GetArraysWrite(ScalarType*&, ScalarType*&, ScalarType*&);
    PetscErrorCode GetArraysWrite(ScalarType**);
    
    PetscErrorCode RestoreArrays(ScalarType*&, ScalarType*&, ScalarType*&);
    PetscErrorCode RestoreArrays(ScalarType**);
    PetscErrorCode RestoreArraysRead(const ScalarType*&, const ScalarType*&, const ScalarType*&);
    PetscErrorCode RestoreArraysRead(const ScalarType**);
    PetscErrorCode RestoreArraysReadWrite(ScalarType*&, ScalarType*&, ScalarType*&);
    PetscErrorCode RestoreArraysReadWrite(ScalarType**);
    
    PetscErrorCode RestoreArrays();
    
    PetscErrorCode DebugInfo(std::string, int, const char*);

    PetscErrorCode WAXPY(ScalarType, VecField*, VecField*);
    PetscErrorCode AXPY(ScalarType, VecField*);
    
    /*! compute norm of vector field */
    PetscErrorCode Norm(Vec);
    PetscErrorCode Norm(ScalarType&);
    PetscErrorCode Norm2(ScalarType&);
    PetscErrorCode Norm(ScalarType&, ScalarType&, ScalarType&);
    
    PetscErrorCode IsZero(bool&);

    size_t GetSize() const;


    // individual components
    Vec m_X1;
    Vec m_X2;
    Vec m_X3;
    
    Vec m_X;
    
    PetscErrorCode GetRawVector(ScalarType*&);
 private:
    typedef enum {None, Read, Write, ReadWrite} AccessType;
    AccessType m_Type;
    IntType m_Allocated;
    bool m_Owend;
    bool m_OwendX;
    
    ScalarType *m_RawPtr;
    
    union {
      ScalarType *m_Ptr[3];
      const ScalarType *m_ConstPtr[3];
    };
    
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);

    PetscErrorCode Allocate(void);
    PetscErrorCode Allocate(IntType,IntType);
    PetscErrorCode Allocate(int);

    RegOpt* m_Opt;
};




}  // end of name space




#endif  // _VECFIELD_H_
