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

#ifndef _VECFIELD_CPP_
#define _VECFIELD_CPP_

#include "VecField.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
VecField::VecField() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
VecField::~VecField() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
VecField::VecField(RegOpt* opt) {
    this->Initialize();
    this->SetOpt(opt);
    this->Allocate();
    this->SetValue(0.);
}


/********************************************************************
 * @brief constructor
 *******************************************************************/
VecField::VecField(RegOpt* opt, Vec x1, Vec x2, Vec x3) {
    this->Initialize();
    this->SetOpt(opt);
    this->m_X1 = x1;
    this->m_X2 = x2;
    this->m_X3 = x3;
    this->m_Allocated = this->m_Opt->m_Domain.nl*3;
    this->m_Owend = false;
    this->m_OwendX = false;
    //this->SetValue(0.);
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
VecField::VecField(RegOpt* opt, Vec x) {
    this->Initialize();
    this->SetOpt(opt);
    
    this->m_X = x;
    
    IntType nl = this->m_Opt->m_Domain.nl;
    IntType ng = this->m_Opt->m_Domain.ng;
    
    GetRawPointer(this->m_X, &this->m_RawPtr);
    
    #ifdef REG_HAS_CUDA
      if (this->m_Opt->rank_cnt > 1) {
        VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[0*nl], &this->m_X1);
        VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[1*nl], &this->m_X2);
        VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[2*nl], &this->m_X3);
      } else {
        VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[0*nl], &this->m_X1);
        VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[1*nl], &this->m_X2);
        VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[2*nl], &this->m_X3);
      }
      //ierr = VecSetType(this->m_X1, VECCUDA); CHKERRQ(ierr);
      //ierr = VecSetType(this->m_X2, VECCUDA); CHKERRQ(ierr);
      //ierr = VecSetType(this->m_X3, VECCUDA); CHKERRQ(ierr);
    #else
      if (this->m_Opt->rank_cnt > 1) {
        VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[0*nl], &this->m_X1);
        VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[1*nl], &this->m_X2);
        VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[2*nl], &this->m_X3);
      } else {
        VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[0*nl], &this->m_X1);
        VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[1*nl], &this->m_X2);
        VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[2*nl], &this->m_X3);
      }
    #endif
    
    this->m_Allocated = this->m_Opt->m_Domain.nl*3;
    this->m_Owend = true;
    this->m_OwendX = false;
    //this->SetValue(0.);
}


/********************************************************************
 * @brief constructor
 *******************************************************************/
VecField::VecField(RegOpt* opt, int level) {
    this->Initialize();
    this->SetOpt(opt);
    this->Allocate(level);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
VecField::VecField(IntType nl, IntType ng) {
    this->Initialize();
    this->Allocate(nl, ng);
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode VecField::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;

    this->m_X1 = NULL;
    this->m_X2 = NULL;
    this->m_X3 = NULL;
    this->m_X = NULL;

    this->m_Allocated = 0;
    
    this->m_Type = AccessType::None;
    this->m_Ptr[0] = nullptr;
    this->m_Ptr[1] = nullptr;
    this->m_Ptr[2] = nullptr;
    
    this->m_RawPtr = nullptr;
    
    this->m_Owend = false;
    this->m_OwendX = false;
    
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode VecField::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (this->m_Owend) {
      if (this->m_X1 != NULL) {
          ierr = VecDestroy(&this->m_X1); CHKERRQ(ierr);
          this->m_X1 = NULL;
      }
      if (this->m_X2 != NULL) {
          ierr = VecDestroy(&this->m_X2); CHKERRQ(ierr);
          this->m_X2 = NULL;
      }
      if (this->m_X3 != NULL) {
          ierr = VecDestroy(&this->m_X3); CHKERRQ(ierr);
          this->m_X3 = NULL;
      }
    }
    if (this->m_RawPtr) {
      ierr = RestoreRawPointer(this->m_X, &this->m_RawPtr); CHKERRQ(ierr);
      this->m_RawPtr = nullptr;
    }
    if (this->m_OwendX) {
      if (this->m_X != NULL) {
        ierr = VecDestroy(&this->m_X); CHKERRQ(ierr);
        this->m_X = NULL;
      }
    }
       
    this->m_Type = AccessType::None;
    this->m_Ptr[0] = nullptr;
    this->m_Ptr[1] = nullptr;
    this->m_Ptr[2] = nullptr;

    PetscFunctionReturn(0);
}



/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
PetscErrorCode VecField::SetVector(Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr = this->ClearMemory(); CHKERRQ(ierr);
    
    this->m_X = x;
    
    IntType nl = this->m_Opt->m_Domain.nl;
    IntType ng = this->m_Opt->m_Domain.ng;
    
    GetRawPointer(this->m_X, &this->m_RawPtr);
    
    #ifdef REG_HAS_CUDA
      if (this->m_Opt->rank_cnt > 1) {
        ierr = VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[0*nl], &this->m_X1);
        ierr = VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[1*nl], &this->m_X2);
        ierr = VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[2*nl], &this->m_X3);
      } else {
        ierr = VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[0*nl], &this->m_X1);
        ierr = VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[1*nl], &this->m_X2);
        ierr = VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[2*nl], &this->m_X3);
      }
      //ierr = VecSetType(this->m_X1, VECCUDA); CHKERRQ(ierr);
      //ierr = VecSetType(this->m_X2, VECCUDA); CHKERRQ(ierr);
      //ierr = VecSetType(this->m_X3, VECCUDA); CHKERRQ(ierr);
    #else
      if (this->m_Opt->rank_cnt > 1) {
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[0*nl], &this->m_X1);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[1*nl], &this->m_X2);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[2*nl], &this->m_X3);
      } else {
        ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[0*nl], &this->m_X1);
        ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[1*nl], &this->m_X2);
        ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[2*nl], &this->m_X3);
      }
    #endif
    
    this->m_Allocated = this->m_Opt->m_Domain.nl*3;
    this->m_Owend = true;
    this->m_OwendX = false;

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief set the options
 *******************************************************************/
PetscErrorCode VecField::SetOpt(RegOpt* opt) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(opt != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_Opt = opt;

    PetscFunctionReturn(ierr);
}

size_t VecField::GetSize() const {
  return this->m_Allocated*sizeof(ScalarType);
}



/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
PetscErrorCode VecField::GetSize(IntType& nl, IntType& ng) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    PetscFunctionBegin;

    //  get sizes
    ierr = VecGetSize(this->m_X1, &ng); CHKERRQ(ierr);
    ierr = VecGetLocalSize(this->m_X1, &nl); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
PetscErrorCode VecField::Allocate() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr = this->ClearMemory(); CHKERRQ(ierr);

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    ierr = this->Allocate(nl, ng); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
PetscErrorCode VecField::Allocate(int level) {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr = this->ClearMemory(); CHKERRQ(ierr);

    nl = this->m_Opt->m_GridCont.nl[level];
    ng = this->m_Opt->m_GridCont.ng[level];

    ierr = this->Allocate(nl, ng); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
PetscErrorCode VecField::Allocate(IntType nl, IntType ng) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr = this->ClearMemory(); CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X, 3*nl, 3*ng); CHKERRQ(ierr);
    #ifdef REG_HAS_CUDA
        ierr = VecSetType(this->m_X, VECCUDA); CHKERRQ(ierr);
    #else
        ierr = VecSetFromOptions(this->m_X); CHKERRQ(ierr);
    #endif
    
    ierr = GetRawPointer(this->m_X, &this->m_RawPtr); CHKERRQ(ierr);
    
    #ifdef REG_HAS_CUDA
      if (this->m_Opt->rank_cnt > 1) {
        ierr = VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[0*nl], &this->m_X1);
        ierr = VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[1*nl], &this->m_X2);
        ierr = VecCreateMPICUDAWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[2*nl], &this->m_X3);
      } else {
        ierr = VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[0*nl], &this->m_X1);
        ierr = VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[1*nl], &this->m_X2);
        ierr = VecCreateSeqCUDAWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[2*nl], &this->m_X3);
      }
      //ierr = VecSetType(this->m_X1, VECCUDA); CHKERRQ(ierr);
      //ierr = VecSetType(this->m_X2, VECCUDA); CHKERRQ(ierr);
      //ierr = VecSetType(this->m_X3, VECCUDA); CHKERRQ(ierr);
    #else
      if (this->m_Opt->rank_cnt > 1) {
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[0*nl], &this->m_X1);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[1*nl], &this->m_X2);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nl, ng, &this->m_RawPtr[2*nl], &this->m_X3);
      } else {
        ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[0*nl], &this->m_X1);
        ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[1*nl], &this->m_X2);
        ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nl, &this->m_RawPtr[2*nl], &this->m_X3);
      }
    #endif
    
    /*// allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X1); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X1, nl, ng); CHKERRQ(ierr);
    #ifdef REG_HAS_CUDA
        ierr = VecSetType(this->m_X1, VECCUDA); CHKERRQ(ierr);
    #else
        ierr = VecSetFromOptions(this->m_X1); CHKERRQ(ierr);
    #endif

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X2); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X2, nl, ng); CHKERRQ(ierr);
    #ifdef REG_HAS_CUDA
        ierr = VecSetType(this->m_X2, VECCUDA); CHKERRQ(ierr);
    #else
        ierr = VecSetFromOptions(this->m_X2); CHKERRQ(ierr);
    #endif

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X3); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X3, nl, ng); CHKERRQ(ierr);
    #ifdef REG_HAS_CUDA
        ierr = VecSetType(this->m_X3, VECCUDA); CHKERRQ(ierr);
    #else
        ierr = VecSetFromOptions(this->m_X3); CHKERRQ(ierr);
    #endif*/
    
    this->m_Allocated = 3*nl;
    this->m_Owend = true;
    this->m_OwendX = true;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief copy input vector field
 *******************************************************************/
PetscErrorCode VecField::Copy(VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecCopy(v->m_X1, this->m_X1); CHKERRQ(ierr);
    ierr = VecCopy(v->m_X2, this->m_X2); CHKERRQ(ierr);
    ierr = VecCopy(v->m_X3, this->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief print debug info
 *******************************************************************/
PetscErrorCode VecField::DebugInfo(std::string str, int line, const char *file) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  if (this->m_Opt->m_Verbosity > 2) {
      ScalarType maxval, minval, nvx1, nvx2, nvx3;
      std::stringstream ss;
      
    if (this->m_Opt->m_Verbosity > 3) {
      ScalarType *p_x1, *p_x2, *p_x3;
      IntType nl;
      size_t fingerprint = 0;
      ierr = VecGetArray(this->m_X1, &p_x1); CHKERRQ(ierr);
      ierr = VecGetArray(this->m_X2, &p_x2); CHKERRQ(ierr);
      ierr = VecGetArray(this->m_X3, &p_x3); CHKERRQ(ierr);
      // compute size of each individual component
      nl = this->m_Opt->m_Domain.nl;
  #pragma omp parallel for
      for (IntType i = 0; i < nl; ++i) {
#if defined(PETSC_USE_REAL_SINGLE)
        fingerprint += reinterpret_cast<uint32_t*>(p_x1)[i];
        fingerprint += reinterpret_cast<uint32_t*>(p_x2)[i];
        fingerprint += reinterpret_cast<uint32_t*>(p_x2)[i];
#else
        fingerprint += reinterpret_cast<uint64_t*>(p_x1)[i];
        fingerprint += reinterpret_cast<uint64_t*>(p_x2)[i];
        fingerprint += reinterpret_cast<uint64_t*>(p_x2)[i];
#endif
      }
      ierr = VecRestoreArray(this->m_X1, &p_x1); CHKERRQ(ierr);
      ierr = VecRestoreArray(this->m_X2, &p_x2); CHKERRQ(ierr);
      ierr = VecRestoreArray(this->m_X3, &p_x3); CHKERRQ(ierr);
      
      ss  << str << " hash: " << std::hex << fingerprint;
      ierr = DbgMsgCall(ss.str(), line, file, 3); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();
    }

      ierr = this->Norm(nvx1, nvx2, nvx3); CHKERRQ(ierr);
      ss  << str << " norm: (" << std::scientific
          << nvx1 << "," << nvx2 << "," << nvx3 <<")";
      ierr = DbgMsgCall(ss.str(),line,file,2); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();
      
      ierr = VecMax(this->m_X1, NULL, &maxval); CHKERRQ(ierr);
      ierr = VecMin(this->m_X1, NULL, &minval); CHKERRQ(ierr);
      ss  << str << " min/max: [" << std::scientific
          << minval << "," << maxval << "]";
      ierr = DbgMsgCall(ss.str(),0,0,2); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();
      ierr = VecMax(this->m_X2, NULL, &maxval); CHKERRQ(ierr);
      ierr = VecMin(this->m_X2, NULL, &minval); CHKERRQ(ierr);
      ss  << std::string(str.size(),' ') << "          [" << std::scientific
          << minval << "," << maxval << "] ";
      ierr = DbgMsgCall(ss.str(),0,0,2); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();
      ierr = VecMax(this->m_X3, NULL, &maxval); CHKERRQ(ierr);
      ierr = VecMin(this->m_X3, NULL, &minval); CHKERRQ(ierr);
      ss  << std::string(str.size(), ' ') << "          [" << std::scientific
          << minval << "," << maxval << "]";
      ierr = DbgMsgCall(ss.str(),0,0,2); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();
  }
  
  PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief set value
 *******************************************************************/
PetscErrorCode VecField::SetValue(ScalarType value) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

#ifdef REG_HAS_CUDA
    ScalarType *pW[3];
    ierr = this->GetArraysWrite(pW); CHKERRQ(ierr);
    
    ierr = reg::SetValue(pW[0], value, this->m_Opt->m_Domain.nl); CHKERRQ(ierr);
    ierr = reg::SetValue(pW[1], value, this->m_Opt->m_Domain.nl); CHKERRQ(ierr);
    ierr = reg::SetValue(pW[2], value, this->m_Opt->m_Domain.nl); CHKERRQ(ierr);
    
    ierr = this->RestoreArrays(); CHKERRQ(ierr);
#else
    ierr = VecSet(this->m_X1, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X2, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X3, value); CHKERRQ(ierr);
#endif

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::GetArrays(ScalarType*& p_x1,
                                   ScalarType*& p_x2,
                                   ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Type == AccessType::None) {
      ierr = GetRawPointer(this->m_X1, &this->m_Ptr[0]); CHKERRQ(ierr);
      ierr = GetRawPointer(this->m_X2, &this->m_Ptr[1]); CHKERRQ(ierr);
      ierr = GetRawPointer(this->m_X3, &this->m_Ptr[2]); CHKERRQ(ierr);
      this->m_Type = AccessType::ReadWrite;
    } else {
      ierr = Assert(this->m_Type == AccessType::ReadWrite, "can't access with different types"); CHKERRQ(ierr);
    }
    
    p_x1 = this->m_Ptr[0];
    p_x2 = this->m_Ptr[1];
    p_x3 = this->m_Ptr[2];

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::GetArrays(ScalarType** p_x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = this->GetArrays(p_x[0], p_x[1], p_x[2]); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::GetArraysWrite(ScalarType*& p_x1,
                                   ScalarType*& p_x2,
                                   ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Type == AccessType::None) {
      ierr = GetRawPointerWrite(this->m_X1, &this->m_Ptr[0]); CHKERRQ(ierr);
      ierr = GetRawPointerWrite(this->m_X2, &this->m_Ptr[1]); CHKERRQ(ierr);
      ierr = GetRawPointerWrite(this->m_X3, &this->m_Ptr[2]); CHKERRQ(ierr);
      this->m_Type = AccessType::Write;
    } else {
      ierr = Assert(this->m_Type == AccessType::Write, "can't access with different types"); CHKERRQ(ierr);
    }
    
    p_x1 = this->m_Ptr[0];
    p_x2 = this->m_Ptr[1];
    p_x3 = this->m_Ptr[2];

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::GetArraysWrite(ScalarType** p_x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = this->GetArraysWrite(p_x[0], p_x[1], p_x[2]); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::GetArraysRead(const ScalarType*& p_x1,
                                       const ScalarType*& p_x2,
                                       const ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (this->m_Type == AccessType::None) {
      ierr = GetRawPointerRead(this->m_X1, &this->m_ConstPtr[0]); CHKERRQ(ierr);
      ierr = GetRawPointerRead(this->m_X2, &this->m_ConstPtr[1]); CHKERRQ(ierr);
      ierr = GetRawPointerRead(this->m_X3, &this->m_ConstPtr[2]); CHKERRQ(ierr);
      this->m_Type = AccessType::Read;
    } else {
      ierr = Assert(this->m_Type == AccessType::Read, "can't access with different types"); CHKERRQ(ierr);
    }
    
    p_x1 = this->m_ConstPtr[0];
    p_x2 = this->m_ConstPtr[1];
    p_x3 = this->m_ConstPtr[2];

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::GetArraysRead(const ScalarType** p_x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->GetArraysRead(p_x[0], p_x[1], p_x[2]); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief pointwise scale of a vector field
 *******************************************************************/
PetscErrorCode VecField::RestoreArrays(ScalarType*& p_x1,
                                       ScalarType*& p_x2,
                                       ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief pointwise scale of a vector field
 *******************************************************************/
PetscErrorCode VecField::RestoreArrays(ScalarType** p_x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief restore arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::RestoreArraysRead(const ScalarType*& p_x1,
                                           const ScalarType*& p_x2,
                                           const ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief restore arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::RestoreArraysRead(const ScalarType** p_x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief get arrays of vector field for read/write
 *******************************************************************/
PetscErrorCode VecField::GetArraysReadWrite(ScalarType*& p_x1,
                                            ScalarType*& p_x2,
                                            ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Type == AccessType::None) {
      ierr = GetRawPointer(this->m_X1, &this->m_Ptr[0]); CHKERRQ(ierr);
      ierr = GetRawPointer(this->m_X2, &this->m_Ptr[1]); CHKERRQ(ierr);
      ierr = GetRawPointer(this->m_X3, &this->m_Ptr[2]); CHKERRQ(ierr);
      this->m_Type = AccessType::ReadWrite;
    } else {
      ierr = Assert(this->m_Type == AccessType::ReadWrite, "can't access with different types"); CHKERRQ(ierr);
    }
    
    p_x1 = this->m_Ptr[0];
    p_x2 = this->m_Ptr[1];
    p_x3 = this->m_Ptr[2];
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief get arrays of vector field for read/write
 *******************************************************************/
PetscErrorCode VecField::GetArraysReadWrite(ScalarType** p_x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->GetArraysReadWrite(p_x[0], p_x[1], p_x[2]); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief resotre arrays of vector field for read/write
 *******************************************************************/
PetscErrorCode VecField::RestoreArraysReadWrite(ScalarType*& p_x1,
                                            ScalarType*& p_x2,
                                            ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;

    ierr = this->RestoreArrays(); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief resotre arrays of vector field for read/write
 *******************************************************************/
PetscErrorCode VecField::RestoreArraysReadWrite(ScalarType** p_x) {
    PetscErrorCode ierr = 0;

    ierr = this->RestoreArrays(); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief resotre arrays of vector field
 *******************************************************************/
PetscErrorCode VecField::RestoreArrays() {
    PetscErrorCode ierr = 0;

    ierr = Assert(this->m_Type != AccessType::None, "can't restore without access"); CHKERRQ(ierr);
  
    switch (this->m_Type) {
    case AccessType::Read:
      ierr = RestoreRawPointerRead(this->m_X1, &this->m_ConstPtr[0]); CHKERRQ(ierr);
      ierr = RestoreRawPointerRead(this->m_X2, &this->m_ConstPtr[1]); CHKERRQ(ierr);
      ierr = RestoreRawPointerRead(this->m_X3, &this->m_ConstPtr[2]); CHKERRQ(ierr);
      this->m_Ptr[0] = nullptr;
      this->m_Ptr[1] = nullptr;
      this->m_Ptr[2] = nullptr;
      this->m_Type = AccessType::None;
      break;
    case AccessType::Write:
      ierr = RestoreRawPointerWrite(this->m_X1, &this->m_Ptr[0]); CHKERRQ(ierr);
      ierr = RestoreRawPointerWrite(this->m_X2, &this->m_Ptr[1]); CHKERRQ(ierr);
      ierr = RestoreRawPointerWrite(this->m_X3, &this->m_Ptr[2]); CHKERRQ(ierr);
      this->m_Ptr[0] = nullptr;
      this->m_Ptr[1] = nullptr;
      this->m_Ptr[2] = nullptr;
      this->m_Type = AccessType::None;
      break;
    case AccessType::ReadWrite:
      ierr = RestoreRawPointerReadWrite(this->m_X1, &this->m_Ptr[0]); CHKERRQ(ierr);
      ierr = RestoreRawPointerReadWrite(this->m_X2, &this->m_Ptr[1]); CHKERRQ(ierr);
      ierr = RestoreRawPointerReadWrite(this->m_X3, &this->m_Ptr[2]); CHKERRQ(ierr);
      this->m_Ptr[0] = nullptr;
      this->m_Ptr[1] = nullptr;
      this->m_Ptr[2] = nullptr;
      this->m_Type = AccessType::None;
      break;
    default:
      break;
    };
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief sets the individual components of a vector field;
 * the input is a flat petsc array
 *******************************************************************/
PetscErrorCode VecField::SetComponents(Vec w, std::string format) {
    PetscErrorCode ierr = 0;
    const ScalarType *p_w = NULL;

    PetscFunctionBegin;

    ierr = GetRawPointerRead(w, &p_w); CHKERRQ(ierr);
    ierr = this->SetComponents(p_w, format); CHKERRQ(ierr);
    ierr = RestoreRawPointerRead(w, &p_w); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


PetscErrorCode VecField::SetComponents(const ScalarType *pX, std::string format) {
    PetscErrorCode ierr = 0;
    IntType nl;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL;

    PetscFunctionBegin;
    
    ierr = this->GetArraysWrite(p_x1, p_x2, p_x3);

    // compute size of each individual component
    nl = this->m_Opt->m_Domain.nl;
    
    bool block = true;
    if (format == "block") {
      block = true;
    } else if (format == "stride") {
      block = false;
    } else {
      ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)
    if (block) {
      ierr = cudaMemcpy(static_cast<void*>(p_x1), static_cast<const void*>(pX), 
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
      ierr = cudaMemcpy(static_cast<void*>(p_x2), static_cast<const void*>(&pX[nl]),
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
      ierr = cudaMemcpy(static_cast<void*>(p_x3), static_cast<const void*>(&pX[2*nl]),
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
      ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
    } else {
        //printf("setting strided components\n");
      ierr = CopyStridedFromFlatVec(p_x1, p_x2, p_x3, pX, nl); CHKERRQ(ierr);
    }
#else
    if (block) {
#pragma omp parallel for
      for (IntType i = 0; i < nl; ++i) {
        p_x1[i] = pX[i     ];
        p_x2[i] = pX[i+  nl];
        p_x3[i] = pX[i+2*nl];
      }
    } else {
#pragma omp parallel for
      for (IntType i = 0; i < nl; ++i) {
        p_x1[i] = pX[3*i+0];
        p_x2[i] = pX[3*i+1];
        p_x3[i] = pX[3*i+2];
      }
    }
#endif

    ierr = this->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}
PetscErrorCode VecField::SetComponents(const ScalarType *pX1, const ScalarType *pX2, const ScalarType *pX3) {
    PetscErrorCode ierr = 0;
    IntType nl;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL;

    PetscFunctionBegin;

    ierr = this->GetArraysWrite(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    // compute size of each individual component
    nl = this->m_Opt->m_Domain.nl;

#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)
      ierr = cudaMemcpy(static_cast<void*>(p_x1), static_cast<const void*>(pX1), 
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
      ierr = cudaMemcpy(static_cast<void*>(p_x2), static_cast<const void*>(pX2),
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
      ierr = cudaMemcpy(static_cast<void*>(p_x3), static_cast<const void*>(pX3),
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
#else
#pragma omp parallel for
    for (IntType i = 0; i < nl; ++i) {
      p_x1[i] = pX1[i];
      p_x2[i] = pX2[i];
      p_x3[i] = pX3[i];
    }
#endif

    ierr = this->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VecField::GetRawVector(ScalarType*& ptr) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_OwendX) {
      ptr = this->m_RawPtr;
    } else {
      ptr = nullptr;
    }

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief get components of vector field and store them
 * in a flat vector
 *******************************************************************/
PetscErrorCode VecField::GetComponents(Vec w, std::string format) {
    PetscErrorCode ierr = 0;
    ScalarType *p_w = NULL;

    PetscFunctionBegin;

    ierr = GetRawPointerWrite(w, &p_w); CHKERRQ(ierr);
    ierr = this->GetComponents(p_w, format); CHKERRQ(ierr);
    ierr = RestoreRawPointerWrite(w, &p_w); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


PetscErrorCode VecField::GetComponents(ScalarType *pX, std::string format) {
    PetscErrorCode ierr = 0;
    IntType nl;
    const ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL;

    PetscFunctionBegin;

    ierr = this->GetArraysRead(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    // compute size of each individual component
    nl = this->m_Opt->m_Domain.nl;
    
    bool block = true;
    if (format == "block") {
      block = true;
    } else if (format == "stride") {
      block = false;
    } else {
      ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)
    if (block) {
      ierr = cudaMemcpy(static_cast<void*>(pX), static_cast<const void*>(p_x1), 
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
      ierr = cudaMemcpy(static_cast<void*>(&pX[nl]), static_cast<const void*>(p_x2), 
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
      ierr = cudaMemcpy(static_cast<void*>(&pX[2*nl]), static_cast<const void*>(p_x3), 
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
    } else {
        //printf("getting strided components\n");
      ierr = CopyStridedToFlatVec(pX, p_x1, p_x2, p_x3, nl); CHKERRCUDA(ierr);
    }
#else
    if (block) {
#pragma omp parallel for
      for (IntType i = 0; i < nl; ++i) {
        pX[i     ] = p_x1[i];
        pX[i+  nl] = p_x2[i];
        pX[i+2*nl] = p_x3[i];
      }
    } else {
#pragma omp parallel for
      for (IntType i = 0; i < nl; ++i) {
        pX[3*i+0] = p_x1[i];
        pX[3*i+1] = p_x2[i];
        pX[3*i+2] = p_x3[i];
      }
    }
#endif

    ierr = this->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode VecField::GetComponents(ScalarType *pX1, ScalarType *pX2, ScalarType *pX3) {
    PetscErrorCode ierr = 0;
    IntType nl;
    const ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL;

    PetscFunctionBegin;

    ierr = this->GetArraysRead(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    // compute size of each individual component
    nl = this->m_Opt->m_Domain.nl;

#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)
    ierr = cudaMemcpy(static_cast<void*>(pX1), static_cast<const void*>(p_x1), 
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
    ierr = cudaMemcpy(static_cast<void*>(pX2), static_cast<const void*>(p_x2), 
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
    ierr = cudaMemcpy(static_cast<void*>(pX3), static_cast<const void*>(p_x3), 
                sizeof(ScalarType)*nl, cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
#else
#pragma omp parallel for
    for (IntType i = 0; i < nl; ++i) {
      pX1[i] = p_x1[i];
      pX2[i] = p_x2[i];
      pX3[i] = p_x3[i];
    }
#endif

    ierr = this->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief scale vector by scalar value
 *******************************************************************/
PetscErrorCode VecField::Scale(ScalarType value) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecScale(this->m_X1, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X2, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X3, value); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief scale vector by scalar value
 *******************************************************************/
PetscErrorCode VecField::Scale(ScalarType* value) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecScale(this->m_X1, value[0]); CHKERRQ(ierr);
    ierr = VecScale(this->m_X2, value[1]); CHKERRQ(ierr);
    ierr = VecScale(this->m_X3, value[2]); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief pointwise scale of vector field
 *******************************************************************/
PetscErrorCode VecField::Scale(Vec s) {
    PetscErrorCode ierr = 0;
    IntType nl;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL, *p_s = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(s, &nl); CHKERRQ(ierr);

    // get pointers
    ierr = GetRawPointer(s, &p_s); CHKERRQ(ierr);
    ierr = this->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

#ifdef REG_HAS_CUDA
    ierr = WrngMsg("Not implemented for CUDA"); CHKERRQ(ierr);
#endif

#pragma omp parallel
{
    ScalarType scale;
#pragma omp for
    for (IntType i = 0; i < nl; ++i) {
        scale = p_s[i];
        p_v1[i] = scale*p_v1[i];
        p_v2[i] = scale*p_v2[i];
        p_v3[i] = scale*p_v3[i];
    }
}  // pragma omp parallel

    // get pointers
    ierr = this->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = RestoreRawPointer(s, &p_s); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief pointwise scale of a vector field
 *******************************************************************/
PetscErrorCode VecField::Scale(VecField* v, Vec s) {
    PetscErrorCode ierr = 0;
    IntType nl;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL, *p_s = NULL,
                *p_sv1 = NULL, *p_sv2 = NULL, *p_sv3 = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(s, &nl); CHKERRQ(ierr);

    // get pointers
    ierr = GetRawPointer(s, &p_s); CHKERRQ(ierr);
    ierr = this->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = v->GetArrays(p_sv1, p_sv2, p_sv3); CHKERRQ(ierr);
    
#ifdef REG_HAS_CUDA
    ierr = WrngMsg("Not implemented for CUDA"); CHKERRQ(ierr);
#endif

#pragma omp parallel
{
    ScalarType scale;
#pragma omp for
    for (IntType i = 0; i < nl; ++i) {
        scale = p_s[i];
        p_sv1[i] = scale*p_v1[i];
        p_sv2[i] = scale*p_v2[i];
        p_sv3[i] = scale*p_v3[i];
    }
}  // pragma omp parallel

    // get pointers
    ierr = RestoreRawPointer(s, &p_s); CHKERRQ(ierr);
    ierr = this->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = v->RestoreArrays(p_sv1, p_sv2, p_sv3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief interface for AXPY
 *******************************************************************/
PetscErrorCode VecField::AXPY(ScalarType s, VecField* v) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    ierr = VecAXPY(this->m_X1, s, v->m_X1); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_X2, s, v->m_X2); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_X3, s, v->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief interface for WAXPY
 *******************************************************************/
PetscErrorCode VecField::WAXPY(ScalarType s, VecField* v, VecField* w) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    ierr = VecWAXPY(this->m_X1, s, v->m_X1, w->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X2, s, v->m_X2, w->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X3, s, v->m_X3, w->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief compute pointwise norm of vector field
 *******************************************************************/
PetscErrorCode VecField::Norm(Vec xnorm) {
    PetscErrorCode ierr = 0;
    IntType nl;
    const ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL;
    ScalarType *p_x = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(xnorm, &nl); CHKERRQ(ierr);

    ierr = this->GetArraysRead(p_x1, p_x2, p_x3); CHKERRQ(ierr);
    ierr = GetRawPointer(xnorm, &p_x); CHKERRQ(ierr);
    
#ifdef REG_HAS_CUDA
    ierr = WrngMsg("Not implemented for CUDA"); CHKERRQ(ierr);
#else
#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i) {
        p_x[i] = PetscSqrtReal(p_x1[i]*p_x1[i]
                             + p_x2[i]*p_x2[i]
                             + p_x3[i]*p_x3[i]);
    }
}  // pragma omp parallel
#endif

    ierr = RestoreRawPointer(xnorm, &p_x); CHKERRQ(ierr);
    ierr = this->RestoreArraysRead(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute pointwise norm of vector field
 *******************************************************************/
PetscErrorCode VecField::Norm(ScalarType& value) {
    PetscErrorCode ierr = 0;
    ScalarType nvx1, nvx2, nvx3;

    PetscFunctionBegin;

    ierr = VecNorm(this->m_X1, NORM_2, &nvx1); CHKERRQ(ierr);
    ierr = VecNorm(this->m_X2, NORM_2, &nvx2); CHKERRQ(ierr);
    ierr = VecNorm(this->m_X3, NORM_2, &nvx3); CHKERRQ(ierr);

    value = nvx1 + nvx2 + nvx3;

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute pointwise squared norm of vector field
 *******************************************************************/
PetscErrorCode VecField::Norm2(ScalarType& value) {
    PetscErrorCode ierr = 0;
    ScalarType nvx1, nvx2, nvx3;

    PetscFunctionBegin;

    ierr = VecNorm(this->m_X1, NORM_2, &nvx1); CHKERRQ(ierr);
    ierr = VecNorm(this->m_X2, NORM_2, &nvx2); CHKERRQ(ierr);
    ierr = VecNorm(this->m_X3, NORM_2, &nvx3); CHKERRQ(ierr);

    //value = nvx1 + nvx2 + nvx3;
    value = nvx1*nvx1 + nvx2*nvx2 + nvx3*nvx3;

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief compute norm of individual components of vector field
 *******************************************************************/
PetscErrorCode VecField::Norm(ScalarType& nvx1, ScalarType& nvx2, ScalarType& nvx3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecNorm(this->m_X1, NORM_2, &nvx1); CHKERRQ(ierr);
    ierr = VecNorm(this->m_X2, NORM_2, &nvx2); CHKERRQ(ierr);
    ierr = VecNorm(this->m_X3, NORM_2, &nvx3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief check if field is zero
 *******************************************************************/
PetscErrorCode VecField::IsZero(bool &zero) {
    PetscErrorCode ierr = 0;
    ScalarType normv1 = 0.0, normv2 = 0.0, normv3 = 0.0;
    PetscFunctionBegin;

    ierr = VecNorm(this->m_X1, NORM_INFINITY, &normv1); CHKERRQ(ierr);
    ierr = VecNorm(this->m_X2, NORM_INFINITY, &normv2); CHKERRQ(ierr);
    ierr = VecNorm(this->m_X3, NORM_INFINITY, &normv3); CHKERRQ(ierr);

    zero = (normv1 == 0.0) && (normv2 == 0.0) && (normv3 == 0.0);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _VECFIELD_CPP_
