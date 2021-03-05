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

#ifndef _CLAIREBASE_CPP_
#define _CLAIREBASE_CPP_

#include "CLAIREBase.hpp"
#include "TwoLevel.hpp"
#include "nifti1_io.h"

namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
CLAIREBase::CLAIREBase() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
CLAIREBase::~CLAIREBase() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
CLAIREBase::CLAIREBase(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode CLAIREBase::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // pointer for container of velocity field
    this->m_VelocityField = NULL;
    this->m_IncVelocityField = NULL;

    // pointers to images (set from outside; not to be deleted)
    this->m_TemplateImage = NULL;
    this->m_ReferenceImage = NULL;

    this->m_Mask = NULL;

    // temporary internal variables (all of these have to be deleted)
    this->m_WorkScaField1 = NULL;
    this->m_WorkScaField2 = NULL;
    this->m_WorkScaField3 = NULL;
    this->m_WorkScaField4 = NULL;
    this->m_WorkScaField5 = NULL;

    this->m_WorkScaFieldMC = NULL;

    this->m_WorkVecField1 = NULL;
    this->m_WorkVecField2 = NULL;
    this->m_WorkVecField3 = NULL;
    this->m_WorkVecField4 = NULL;
    this->m_WorkVecField5 = NULL;

    this->m_AuxVariable = NULL;
    this->m_CellDensity = NULL;
    this->m_x1hat = NULL;
    this->m_x2hat = NULL;
    this->m_x3hat = NULL;
    
    this->m_GradientState = nullptr;
    this->m_GradientXState = nullptr;

    // objects
    this->m_ReadWrite = NULL;               ///< read / write object
    this->m_Regularization = NULL;          ///< pointer for regularization class
    this->m_DistanceMeasure = NULL;         ///< distance measure
    this->m_SemiLagrangianMethod = NULL;    ///< semi lagranigan
    this->m_DeformationFields = NULL;       ///< interface for computing deformation field (jacobian; mapping; ...)
    this->m_Differentiation = NULL;         ///< interface for differentiation
    this->m_DifferentiationFD = NULL;         ///< interface for differentiation
    this->m_TransportProblem = NULL;        ///< interface for transport equation

    this->m_VelocityIsZero = false;          ///< flag: is velocity zero
    this->m_StoreTimeHistory = true;         ///< flag: store time history (needed for inversion)

    this->m_DeleteControlVariable = true;    ///< flag: clear memory for control variable
    this->m_DeleteIncControlVariable = true; ///< flag: clear memory for incremental control variable

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode CLAIREBase::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

//    if (this->m_Mask != NULL) {
//        ierr = VecDestroy(&this->m_Mask); CHKERRQ(ierr);
//        this->m_Mask = NULL;
//    }

    if (this->m_Opt->m_Verbosity > 2) {
      std::stringstream ss;
      size_t total = 0;
      if (this->m_WorkScaField1) {
        total += this->m_WorkScaField1->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WSF1 " << this->m_WorkScaField1->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkScaField2) {
        total += this->m_WorkScaField2->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WSF2 " << this->m_WorkScaField2->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkScaField3) {
        total += this->m_WorkScaField3->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WSF3 " << this->m_WorkScaField3->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkScaField4) {
        total += this->m_WorkScaField4->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WSF4 " << this->m_WorkScaField4->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkScaField5) {
        total += this->m_WorkScaField5->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WSF5 " << this->m_WorkScaField5->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkVecField1) {
        total += this->m_WorkVecField1->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WVF1 " << this->m_WorkVecField1->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkVecField2) {
        total += this->m_WorkVecField2->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WVF2 " << this->m_WorkVecField2->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkVecField3) {
        total += this->m_WorkVecField3->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WVF3 " << this->m_WorkVecField3->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkVecField4) {
        total += this->m_WorkVecField4->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WVF4 " << this->m_WorkVecField4->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_WorkVecField5) {
        total += this->m_WorkVecField5->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WVF5 " << this->m_WorkVecField5->GetSize();
        DbgMsg3(ss.str());
      }
      
      if (this->m_WorkScaFieldMC) {
        total += this->m_WorkScaFieldMC->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "WSFMC " << this->m_WorkScaFieldMC->GetSize();
        DbgMsg3(ss.str());
      }
      
      if (this->m_TemplateImage) {
        total += this->m_TemplateImage->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "mT " << this->m_TemplateImage->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_ReferenceImage) {
        total += this->m_ReferenceImage->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "mR " << this->m_ReferenceImage->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_AuxVariable) {
        total += this->m_AuxVariable->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "Aux " << this->m_AuxVariable->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_CellDensity) {
        total += this->m_CellDensity->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "Cell " << this->m_CellDensity->GetSize();
        DbgMsg3(ss.str());
      }
      if (this->m_Mask) {
        total += this->m_Mask->GetSize();
        ss.clear(); ss.str(std::string());
        ss << "mask " << this->m_Mask->GetSize();
        DbgMsg3(ss.str());
      }
      ss.clear(); ss.str(std::string());
      
      if (this->m_GradientState) {
        for (int i=0;i<this->m_Opt->m_Domain.nt;++i)
          total += this->m_GradientState[i]->GetSize();
      }
      
      if (this->m_GradientXState) {
        for (int i=0;i<this->m_Opt->m_Domain.nt;++i)
          total += this->m_GradientXState[i]->GetSize();
      }
      
      if (this->m_DeleteControlVariable && this->m_VelocityField) {
        total += this->m_VelocityField->GetSize();
      }
      if (this->m_DeleteIncControlVariable && this->m_IncVelocityField) {
        total += this->m_IncVelocityField->GetSize();
      }
      
      ss << "memory allocated: "<< std::scientific << total;
      ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
      ss.clear(); ss.str(std::string());
    }
    
    if (this->m_DeleteControlVariable) {
        Free(this->m_VelocityField);
    }
    
    if (this->m_DeleteIncControlVariable) {
        Free(this->m_IncVelocityField);
    }

    Free(this->m_DistanceMeasure);
    Free(this->m_Regularization);
    Free(this->m_SemiLagrangianMethod);
    Free(this->m_DeformationFields);
    Free(this->m_TransportProblem);
    Free(this->m_Differentiation);
    Free(this->m_DifferentiationFD);
    
    Free(this->m_TemplateImage);
    Free(this->m_ReferenceImage);
    Free(this->m_AuxVariable);
    Free(this->m_CellDensity);
    Free(this->m_Mask);
    
    if (this->m_GradientState) {
      for (int i=0;i<this->m_Opt->m_Domain.nc*(this->m_Opt->m_RegFlags.runinversion?this->m_Opt->m_Domain.nt+1:1);++i)
        Free(this->m_GradientState[i]);
    }
    FreeArray(this->m_GradientState);
    if (this->m_GradientXState) {
      for (int i=0;i<this->m_Opt->m_Domain.nt*this->m_Opt->m_Domain.nc;++i)
        Free(this->m_GradientXState[i]);
    }
    FreeArray(this->m_GradientXState);
    
    Free(this->m_WorkScaField1);
    Free(this->m_WorkScaField2);
    Free(this->m_WorkScaField3);
    Free(this->m_WorkScaField4);
    Free(this->m_WorkScaField5);
    Free(this->m_WorkScaFieldMC);
    Free(this->m_WorkVecField1);
    Free(this->m_WorkVecField2);
    Free(this->m_WorkVecField3);
    Free(this->m_WorkVecField4);
    Free(this->m_WorkVecField5);

    FreeMemory(this->m_x1hat);
    FreeMemory(this->m_x2hat);
    FreeMemory(this->m_x3hat);
    /*if (this->m_x1hat != NULL) {
#ifndef REG_HAS_CUDA
        accfft_free(this->m_x1hat);
#else
        cudaFree(this->m_x1hat);
#endif
        this->m_x1hat = NULL;
    }
    if (this->m_x2hat != NULL) {
#ifndef REG_HAS_CUDA
        accfft_free(this->m_x2hat);
#else
        cudaFree(this->m_x2hat);
#endif
        this->m_x2hat = NULL;
    }
    if (this->m_x3hat != NULL) {
#ifndef REG_HAS_CUDA
        accfft_free(this->m_x3hat);
#else
        cudaFree(this->m_x3hat);
#endif
        this->m_x3hat = NULL;
    }*/

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupSpectralData() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
PetscErrorCode CLAIREBase::SetReadWrite(ReadWriteReg* readwrite) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(readwrite != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite = readwrite;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode CLAIREBase::SetReferenceImage(Vec mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);
  
    //ierr = Free(this->m_ReferenceImage); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_ReferenceImage, this->m_Opt, mR, true); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->SetVector(mR); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 3) {
      ierr = this->m_ReferenceImage->DebugInfo("mR", __LINE__, __FILE__);
    }
    
#ifdef BUILD_OPT_TWOLVL_REF
    printf("Beep Boop\n");
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
    
    TwoLevelFinite mgop(this->m_Opt);
    
    this->m_WorkScaField1->Set(0);
    this->m_WorkScaField2->Set(0);
    
    mgop.Restrict(this->m_WorkScaField1, this->m_ReferenceImage);
    mgop.Prolong(this->m_WorkScaField2, this->m_WorkScaField1);
    
    VecAXPY(this->m_WorkScaField2->m_X, -1, this->m_ReferenceImage->m_X);
    
    Vec pvec = this->m_WorkScaField2->m_X;
    ScalarType *data;
    VecGetArray(pvec, &data);
    nifti_image *image = nullptr;
    image = nifti_simple_init_nim();
    image->dim[0] = image->ndim = 3;
    image->dim[1] = image->nx = this->m_Opt->m_Domain.isize[2];
    image->dim[2] = image->ny = this->m_Opt->m_Domain.isize[1];
    image->dim[3] = image->nz = this->m_Opt->m_Domain.isize[0];
    image->nvox = image->nx*image->ny*image->nz;
    image->data = data;
    image->nbyper = sizeof(ScalarType);
    image->datatype = NIFTI_TYPE_FLOAT32;
    image->nifti_type = NIFTI_FTYPE_NIFTI1_1;
    image->fname = nifti_makehdrname("org.nii.gz", image->nifti_type, false, true);
    image->iname = nifti_makeimgname("org.nii.gz", image->nifti_type, false, true);  
    nifti_image_infodump(image);
    nifti_image_write(image);
    image->data = nullptr;
    nifti_image_free(image); 
    VecRestoreArray(pvec, &data);
#endif
    
    // assign pointer
    //this->m_ReferenceImage = mR;
    if (this->m_Opt->m_RegFlags.registerprobmaps) {
        ierr = EnsurePartitionOfUnity(*this->m_ReferenceImage, this->m_Opt->m_Domain.nc); CHKERRQ(ierr);
        ierr = ShowValues(*this->m_ReferenceImage, this->m_Opt->m_Domain.nc); CHKERRQ(ierr);
    }
//    ierr = ShowValues(this->m_ReferenceImage, this->m_Opt->m_Domain.nc); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set template image (i.e., the image to be registered)
 *******************************************************************/
PetscErrorCode CLAIREBase::SetTemplateImage(Vec mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);
    
    //ierr = Free(this->m_TemplateImage); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_TemplateImage, this->m_Opt, mT, true); CHKERRQ(ierr);
    ierr = this->m_TemplateImage->SetVector(mT); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 3) {
      ierr = this->m_TemplateImage->DebugInfo("mT", __LINE__, __FILE__);
    }

    // assign pointer
    //this->m_TemplateImage = mT;
    if (this->m_Opt->m_RegFlags.registerprobmaps) {
        ierr = EnsurePartitionOfUnity(*this->m_TemplateImage, this->m_Opt->m_Domain.nc); CHKERRQ(ierr);
        ierr = ShowValues(*this->m_TemplateImage, this->m_Opt->m_Domain.nc); CHKERRQ(ierr);
    }
//    ierr = ShowValues(this->m_TemplateImage,  this->m_Opt->m_Domain.nc); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode CLAIREBase::GetReferenceImage(Vec& mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
        
    mR = *this->m_ReferenceImage;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get template image (i.e., the image to be registered)
 *******************************************************************/
PetscErrorCode CLAIREBase::GetTemplateImage(Vec& mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    mT = *this->m_TemplateImage;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get template image (i.e., the image to be registered)
 *******************************************************************/
PetscErrorCode CLAIREBase::GetMask(Vec& mask) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // mask can be a null pointer
//    ierr = Assert(this->m_Mask != NULL, "null pointer"); CHKERRQ(ierr);
    if (this->m_Mask) {
      mask = *this->m_Mask;
    } else {
      mask = nullptr;
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief set velocity field
 *******************************************************************/
PetscErrorCode CLAIREBase::SetControlVariable(VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    if (this->m_Opt->m_Verbosity > 2) {
      ierr = DbgMsg2("Setting control variable"); CHKERRQ(ierr);
    }

    if (v == nullptr) {
      ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
      ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
      this->m_DeleteControlVariable = true;
    } else {
      if (this->m_DeleteControlVariable) {
          Free(this->m_VelocityField);
      }
    
      this->m_VelocityField = v;
      this->m_DeleteControlVariable = false;
    }
    
    this->m_VelocityField->DebugInfo("control variable", __LINE__, __FILE__);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set mask (for masking the l2 distance / similarity measure)
 *******************************************************************/
PetscErrorCode CLAIREBase::SetMask(Vec mask) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mask != NULL, "null pointer"); CHKERRQ(ierr);

    //ierr = Free(this->m_Mask); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_Mask, this->m_Opt, mask); CHKERRQ(ierr);
    ierr = this->m_Mask->SetVector(mask); CHKERRQ(ierr);
    // assign pointer
    //this->m_Mask = mask;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set auxilary variable
 *******************************************************************/
PetscErrorCode CLAIREBase::SetAuxVariable(Vec q) {
    PetscErrorCode ierr;
    this->m_Opt->Enter(__func__);

    // switch l2 distance
//    this->m_Opt->m_Distance.type = SL2AUX;

    ierr = Assert(q != NULL, "null pointer"); CHKERRQ(ierr);
    
    //ierr = Free(this->m_AuxVariable); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_AuxVariable, this->m_Opt, q, true); CHKERRQ(ierr);
    ierr = this->m_AuxVariable->SetVector(q); CHKERRQ(ierr);
    
    //this->m_AuxVariable = q;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set auxilary variable
 *******************************************************************/
PetscErrorCode CLAIREBase::SetCellDensity(Vec c) {
    PetscErrorCode ierr;
    this->m_Opt->Enter(__func__);

    // switch l2 distance
//    this->m_Opt->m_Distance.type = SL2AUX;

    ierr = Assert(c != NULL, "null pointer"); CHKERRQ(ierr);
    
    //ierr = Free(this->m_CellDensity); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_CellDensity, this->m_Opt, c, true); CHKERRQ(ierr);
    ierr = this->m_CellDensity->SetVector(c); CHKERRQ(ierr);
    //this->m_CellDensity = c;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
PetscErrorCode CLAIREBase::SetIncControlVariable(VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_DeleteIncControlVariable) {
        Free(this->m_IncVelocityField);
    }
    
    this->m_IncVelocityField = v;
    this->m_DeleteIncControlVariable = false;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
PetscErrorCode CLAIREBase::GetControlVariable(VecField*& v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    // copy buffer
    ierr = v->Copy(this->m_VelocityField); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
PetscErrorCode CLAIREBase::ComputeInitialGuess() {
    PetscErrorCode ierr = 0;
    Vec v = NULL, g = NULL;
    IntType nl, ng;
    PetscFunctionBegin;

    // if velocity field is null pointer, we did not set
    // any initial guess
    ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
    ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);

    if (!this->m_Opt->m_OptPara.usezeroinitialguess) {
        if (this->m_Opt->m_Verbosity > 2) {
          ierr = DbgMsg2("compute nonzero initial guess"); CHKERRQ(ierr);
        }
      
        nl = this->m_Opt->m_Domain.nl;
        ng = this->m_Opt->m_Domain.ng;
        ierr = VecCreate(v, 3*nl, 3*ng); CHKERRQ(ierr);
        ierr = VecCreate(g, 3*nl, 3*ng); CHKERRQ(ierr);

        ierr = this->m_VelocityField->GetComponents(v); CHKERRQ(ierr);

        // TODO: use sobolev gradient
        ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);

        ierr = VecAXPY(v, -1.0, g); CHKERRQ(ierr);

        ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);
    }


    if (v != NULL) {ierr = VecDestroy(&g); CHKERRQ(ierr); v = NULL;}
    if (g != NULL) {ierr = VecDestroy(&v); CHKERRQ(ierr); g = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief check if velocity field is zero
 *******************************************************************/
PetscErrorCode CLAIREBase::IsVelocityZero() {
    PetscErrorCode ierr = 0;
    ScalarType normv1 = 0.0, normv2 = 0.0, normv3 = 0.0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    this->m_VelocityIsZero = false;
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = VecNorm(this->m_VelocityField->m_X1, NORM_INFINITY, &normv1); CHKERRQ(ierr);
    ierr = VecNorm(this->m_VelocityField->m_X2, NORM_INFINITY, &normv2); CHKERRQ(ierr);
    ierr = VecNorm(this->m_VelocityField->m_X3, NORM_INFINITY, &normv3); CHKERRQ(ierr);

    this->m_VelocityIsZero = (normv1 == 0.0) && (normv2 == 0.0) && (normv3 == 0.0);

    if (this->m_Opt->m_Verbosity > 2) {
        if (this->m_VelocityIsZero) {
            ierr = DbgMsg2("zero velocity field"); CHKERRQ(ierr);
        }
    }
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief allocate regularization model
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupDistanceMeasure() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // for the hessian matvec, the template and reference images are not required
//    ierr = Assert(this->m_TemplateImage != NULL, "null pointer (template)"); CHKERRQ(ierr);
//    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer (reference)"); CHKERRQ(ierr);

    // delete regularization if already allocated
    // (should never happen)
    if (this->m_DistanceMeasure != NULL) {
        delete this->m_DistanceMeasure;
        this->m_DistanceMeasure = NULL;
    }

    // switch between regularization norms
    switch (this->m_Opt->m_Distance.type) {
        case SL2:
        {
            ierr = AllocateOnce<DistanceMeasureSL2>(this->m_DistanceMeasure, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case SL2AUX:
        {
            ierr = AllocateOnce<DistanceMeasureSL2aux>(this->m_DistanceMeasure, this->m_Opt); CHKERRQ(ierr);
            // TODO: Fix for 2 level preconditioner (these need to be set in preconditioner)
            ierr = Assert(this->m_CellDensity != NULL, "null pointer (aux 1)"); CHKERRQ(ierr);
            ierr = Assert(this->m_AuxVariable != NULL, "null pointer (aux 2)"); CHKERRQ(ierr);
            ierr = this->m_DistanceMeasure->SetAuxVariable(this->m_CellDensity,1); CHKERRQ(ierr);
            ierr = this->m_DistanceMeasure->SetAuxVariable(this->m_AuxVariable,2); CHKERRQ(ierr);
            break;
        }
        case NCC:
        {
            ierr = AllocateOnce<DistanceMeasureNCC>(this->m_DistanceMeasure, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = reg::ThrowError("distance measure not defined"); CHKERRQ(ierr);
            break;
        }
    }
    // for a pure evaluation of hessian operator, images may not have been set
    if (this->m_TemplateImage != NULL) {
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = DbgMsg2("distance measure: parsing template image"); CHKERRQ(ierr);
        }
        ierr = this->m_DistanceMeasure->SetTemplateImage(this->m_TemplateImage); CHKERRQ(ierr);
    }
    if (this->m_ReferenceImage != NULL) {
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = DbgMsg2("distance measure: parsing reference image"); CHKERRQ(ierr);
        }
       ierr = this->m_DistanceMeasure->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
    }
    if (this->m_Mask != NULL) {
        if (this->m_Opt->m_Verbosity > 1) {
            ierr = DbgMsg1("distance measure: mask enabled"); CHKERRQ(ierr);
        }
        ierr = this->m_DistanceMeasure->SetMask(this->m_Mask); CHKERRQ(ierr);
    }
    
    if (this->m_Opt->m_ObjWts.empty()) {
        for (int i=0; i<this->m_Opt->m_Domain.nc; ++i) this->m_Opt->m_ObjWts.push_back(1); 
    }
    if (this->m_Opt->m_ObjWts.data() != NULL) {
        if (this->m_Opt->m_Verbosity > 1) {
            ierr = DbgMsg("distance measure: objective function weights enabled"); CHKERRQ(ierr);
        }
        ierr = this->m_DistanceMeasure->SetObjectiveFunctionalWeights();
    }

    // Setup scale for distance measure
    // TODO: Do we need to rescale this for the 2level preconditioner? 
    // Currently Temp and Ref image not in scope for restricted problem, so we use old scale value
    if (this->m_TemplateImage != NULL && this->m_ReferenceImage != NULL) { 
    	ierr = this->m_DistanceMeasure->SetupScale(); CHKERRQ(ierr);
    }
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief allocate regularization model
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupRegularization() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // delete regularization if already allocated
    // (should never happen)
    if (this->m_Regularization != NULL) {
        delete this->m_Regularization;
        this->m_Regularization = NULL;
    }

    // switch between regularization norms
    switch (this->m_Opt->m_RegNorm.type) {
        case L2:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate L2 regularization"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<RegularizationL2>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H1:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate H1 regularization"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<RegularizationH1>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H2:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate H2 regularization"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<RegularizationH2>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H3:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate H3 regularization"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<RegularizationH3>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H1SN:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate H1SN regularization"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<RegularizationH1SN>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H2SN:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate H2SN regularization"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<RegularizationH2SN>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H3SN:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate H3SN regularization"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<RegularizationH3SN>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = reg::ThrowError("regularization model not defined"); CHKERRQ(ierr);
        }
    }
    
    ierr = AllocateOnce<DifferentiationSM>(this->m_Differentiation, this->m_Opt); CHKERRQ(ierr);
    ierr = this->m_Regularization->SetDifferentiation(this->m_Differentiation); CHKERRQ(ierr);
    
    //ierr = this->m_Regularization->SetSpectralData(); CHKERRQ(ierr);

    // set the containers for the spectral data
    /*ierr = this->SetupSpectralData(); CHKERRQ(ierr);
    ierr = this->m_Regularization->SetSpectralData(this->m_x1hat,
                                                   this->m_x2hat,
                                                   this->m_x3hat); CHKERRQ(ierr);*/

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief allocate transport problem model
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupTransportProblem() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    // delete regularization if already allocated
    // (should never happen)
    if (this->m_TransportProblem != NULL) {
        delete this->m_TransportProblem;
        this->m_TransportProblem = NULL;
    }

    // switch between regularization norms
    switch (this->m_Opt->m_PDESolver.type) {
        case SL:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate SL transport problem"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<TransportEquationSL>(this->m_TransportProblem, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case RK2:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg1("allocate RK2 transport problem"); CHKERRQ(ierr);
            }
            ierr = AllocateOnce<TransportEquationRK2>(this->m_TransportProblem, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = reg::ThrowError("transport problem model not defined"); CHKERRQ(ierr);
        }
    }

    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField3, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkVecField5, this->m_Opt); CHKERRQ(ierr);

    ierr = this->m_TransportProblem->SetWorkVecField(this->m_WorkVecField1, 1); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkVecField(this->m_WorkVecField2, 2); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkVecField(this->m_WorkVecField3, 3); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkScaField4, this->m_Opt, this->m_WorkVecField5->m_X1); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkScaField5, this->m_Opt, this->m_WorkVecField5->m_X2); CHKERRQ(ierr);
    
    ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField1, 1); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField2, 2); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField3, 3); CHKERRQ(ierr);
    //ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField4, 4); CHKERRQ(ierr);
    //ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField5, 5); CHKERRQ(ierr);
  
    switch (this->m_Opt->m_Diff.diffPDE) {
        case SPECTRAL:
          ierr = AllocateOnce<DifferentiationSM>(this->m_Differentiation, this->m_Opt); CHKERRQ(ierr);
          ierr = this->m_TransportProblem->SetDifferentiation(this->m_Differentiation); CHKERRQ(ierr);
          break;
        case FINITE:
          ierr = AllocateOnce(this->m_DifferentiationFD, this->m_Opt); CHKERRQ(ierr);
          ierr = this->m_TransportProblem->SetDifferentiation(this->m_DifferentiationFD); CHKERRQ(ierr);
          break;
        default:
          ierr = reg::ThrowError("differentiation scheme for transport problem not defined"); CHKERRQ(ierr);
    }
    

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief allocate regularization model
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupDeformationField() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = AllocateOnce(this->m_DeformationFields, this->m_Opt); CHKERRQ(ierr);
    
    ierr = AllocateOnce<DifferentiationSM>(this->m_Differentiation, this->m_Opt); CHKERRQ(ierr);
    if (this->m_Differentiation != NULL) {
      this->m_DeformationFields->SetDifferentiation(this->m_Differentiation);
    }
    
    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField3, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkVecField5, this->m_Opt); CHKERRQ(ierr);
    
    ierr = this->m_DeformationFields->SetWorkVecField(this->m_WorkVecField1, 1); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkVecField(this->m_WorkVecField2, 2); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkVecField(this->m_WorkVecField3, 3); CHKERRQ(ierr);
    
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkScaField4, this->m_Opt, this->m_WorkVecField5->m_X1); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkScaField5, this->m_Opt, this->m_WorkVecField5->m_X2); CHKERRQ(ierr);

    ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField1, 1); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField2, 2); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField3, 3); CHKERRQ(ierr);
    //ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField4, 4); CHKERRQ(ierr);
    //ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField5, 5); CHKERRQ(ierr);

    ierr = this->m_DeformationFields->SetVelocityField(this->m_VelocityField); CHKERRQ(ierr);

    if (this->m_Opt->m_PDESolver.type == SL && this->m_TransportProblem) {
      SemiLagrangian* sl = static_cast<TransportEquationSL*>(this->m_TransportProblem)->GetSemiLagrangian();
      sl->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
      ierr = sl->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
      ierr = this->m_DeformationFields->SetSLM(sl); CHKERRQ(ierr);
    } else {
      ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
      ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
      ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
      ierr = this->m_DeformationFields->SetSLM(this->m_SemiLagrangianMethod); CHKERRQ(ierr);
    }


    ierr = this->m_DeformationFields->SetReadWrite(this->m_ReadWrite); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief evaluate regularization model
 *******************************************************************/
PetscErrorCode CLAIREBase
::EvaluateRegularizationFunctional(ScalarType* value, VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
    ierr = this->m_Regularization->SetWorkScaField(this->m_WorkScaField1); CHKERRQ(ierr);

    ierr = this->m_Regularization->EvaluateFunctional(value, v); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief estimate eigenvalues of hessian
 *******************************************************************/
PetscErrorCode CLAIREBase
::EstimateExtremalHessEigVals(ScalarType &emin,
                              ScalarType &emax) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    ierr = this->m_Regularization->GetExtremeEigValsInvOp(emin, emax); CHKERRQ(ierr);
    emin += 1.0;
    emax += 1.0; // this is crap

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

PetscErrorCode CLAIREBase::PreHessian() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    if (this->m_Opt->m_PDESolver.type == SL) {
      TransportEquationSL *tpeq = static_cast<TransportEquationSL*>(this->m_TransportProblem);
      
      tpeq->GetSemiLagrangian()->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
      tpeq->GetSemiLagrangian()->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
      tpeq->GetSemiLagrangian()->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief pre-processing before the krylov solve
 *******************************************************************/
PetscErrorCode CLAIREBase::PreKrylovSolve(Vec g, Vec x) {
    PetscErrorCode ierr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // switch between hessian operators
    if (this->m_Opt->m_KrylovMethod.matvectype == PRECONDMATVEC) {
        // apply inverse regularization operator
        ierr = this->ApplyInvRegularizationOperator(g, g, false); CHKERRQ(ierr);
    } else if (this->m_Opt->m_KrylovMethod.matvectype == PRECONDMATVECSYM) {
        // apply square root of inverse regularization operator
        ierr = this->ApplyInvRegularizationOperator(g, g, true); CHKERRQ(ierr);
    }

#if 1
    if (this->m_CoarseReg) {
      
      //AllocateOnce(this->m_CoarseReg->m_TemplateImage, this->m_CoarseRegOpt);
      TwoLevelFFT op(this->m_Opt);
      
      //op.Restrict(this->m_CoarseReg->m_TemplateImage, this->m_TemplateImage);
      
     //for (int t=0; t<= this->m_Opt->m_Domain.nt; ++t) {
        const ScalarType *p_m_f;
        ScalarType *p_m_c;
        this->m_StateVariable->GetArrayRead(p_m_f, 0, 0);
        this->m_CoarseReg->m_StateVariable->GetArray(p_m_c, 0, 0);
        
        op.Restrict(p_m_c, p_m_f);
        
        this->m_StateVariable->RestoreArray();
        this->m_CoarseReg->m_StateVariable->RestoreArray();
      //}
      
      op.Restrict(this->m_CoarseReg->m_VelocityField, this->m_VelocityField);
      
      ierr = this->m_CoarseReg->m_TransportProblem->SetStateVariable(this->m_CoarseReg->m_StateVariable); CHKERRQ(ierr);
      ierr = this->m_CoarseReg->m_TransportProblem->SetControlVariable(this->m_CoarseReg->m_VelocityField); CHKERRQ(ierr);
      
      ierr = this->m_CoarseReg->m_TransportProblem->SolveForwardProblem(); CHKERRQ(ierr);
      
      /*{
      Vec pvec = this->m_CoarseReg->m_StateVariable->m_X;
      ScalarType *data;
      VecGetArray(pvec, &data);
      nifti_image *image = nullptr;
      image = nifti_simple_init_nim();
      image->dim[0] = image->ndim = 3;
      image->dim[1] = image->nx = 128;
      image->dim[2] = image->ny = 128;
      image->dim[3] = image->nz = 128;
      image->nvox = image->nx*image->ny*image->nz;
      image->data = &data[0];
      image->nbyper = sizeof(ScalarType);
      image->datatype = NIFTI_TYPE_FLOAT32;
      image->nifti_type = NIFTI_FTYPE_NIFTI1_1;
      image->fname = nifti_makehdrname("m0.nii.gz", image->nifti_type, false, true);
      image->iname = nifti_makeimgname("m0.nii.gz", image->nifti_type, false, true);  
      nifti_image_infodump(image);
      nifti_image_write(image);
      image->data = nullptr;
      nifti_image_free(image); 
      VecRestoreArray(pvec, &data);
      }
      
      {
      Vec pvec = this->m_CoarseReg->m_StateVariable->m_X;
      ScalarType *data;
      VecGetArray(pvec, &data);
      nifti_image *image = nullptr;
      image = nifti_simple_init_nim();
      image->dim[0] = image->ndim = 3;
      image->dim[1] = image->nx = 128;
      image->dim[2] = image->ny = 128;
      image->dim[3] = image->nz = 128;
      image->nvox = image->nx*image->ny*image->nz;
      image->data = &data[128*128*128*4];
      image->nbyper = sizeof(ScalarType);
      image->datatype = NIFTI_TYPE_FLOAT32;
      image->nifti_type = NIFTI_FTYPE_NIFTI1_1;
      image->fname = nifti_makehdrname("m1.nii.gz", image->nifti_type, false, true);
      image->iname = nifti_makeimgname("m1.nii.gz", image->nifti_type, false, true);  
      nifti_image_infodump(image);
      nifti_image_write(image);
      image->data = nullptr;
      nifti_image_free(image); 
      VecRestoreArray(pvec, &data);
      }*/
      
     
      TransportEquationSL *tpeq = static_cast<TransportEquationSL*>(this->m_CoarseReg->m_TransportProblem);
      
      tpeq->GetSemiLagrangian()->SetWorkVecField(this->m_CoarseReg->m_WorkVecField1); CHKERRQ(ierr);
      //tpeq->GetSemiLagrangian()->ComputeTrajectory(this->m_CoarseReg->m_VelocityField, "state"); CHKERRQ(ierr);
      tpeq->GetSemiLagrangian()->ComputeTrajectory(this->m_CoarseReg->m_VelocityField, "adjoint"); CHKERRQ(ierr);
      
      //ierr = AllocateOnce(this->m_CoarseReg->m_SemiLagrangianMethod, this->m_CoarseRegOpt); CHKERRQ(ierr);
      //this->m_CoarseReg->m_SemiLagrangianMethod->SetWorkVecField(this->m_CoarseReg->m_WorkVecField1); CHKERRQ(ierr);
      //this->m_CoarseReg->m_SemiLagrangianMethod->ComputeTrajectory(this->m_CoarseReg->m_VelocityField, "state"); CHKERRQ(ierr);    
      //this->m_CoarseReg->m_SemiLagrangianMethod->ComputeTrajectory(this->m_CoarseReg->m_VelocityField, "adjoint"); CHKERRQ(ierr);
      
      //////////////////////
      /*
      VecField *g_vec = new VecField(m_Opt, g);
      
      
      //this->m_WorkVecField1->SetValue(0.0);
      
      //TwoLevelFFT op(this->m_Opt);
      op.Restrict(this->m_WorkVecField1, g_vec);
      op.Prolong(g_vec, this->m_WorkVecField1);
      delete g_vec;
      */
    }
#endif
    
    /*{
      VecField *g_vec = new VecField(m_Opt, g);
    
      TwoLevelFFT op(this->m_Opt);
      op.Restrict(this->m_WorkVecField1, g_vec);
      op.Prolong(g_vec, this->m_WorkVecField1);
      
      delete g_vec;
    }*/

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief post-processing after the krylov solve
 *******************************************************************/
PetscErrorCode CLAIREBase::PostKrylovSolve(Vec g, Vec x) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    //AllocateOnce(this->m_WorkVecField5, this->m_Opt);
    //ApplyInvHessian(x, g, m_WorkVecField5, true, false);
    //ScalarType hd;
    //hd = this->m_Opt->GetLebesgueMeasure();
    //VecScale(g, 1./hd);
    //SymTwoLevelHessMatVec(x, g);
    //ierr = VecScale(x, hd); CHKERRQ(ierr);
        
    /*{
    ScalarType rTr, pTAp, rel, rel_f;
    
    VecDot(g, g, &rTr);
    
    rel_f = sqrt(rTr);
    
    //if (this->m_Opt->m_KrylovMethod.matvectype == PRECONDMATVEC) {
    //    // apply inverse regularization operator
    //    ierr = this->ApplyInvRegularizationOperator(g, g, false); CHKERRQ(ierr);
    //} else if (this->m_Opt->m_KrylovMethod.matvectype == PRECONDMATVECSYM) {
    //    // apply square root of inverse regularization operator
    //    ierr = this->ApplyInvRegularizationOperator(g, g, true); CHKERRQ(ierr);
    //}
    
    //VecField b_f(this->m_Opt, g);
    //VecField x_f(this->m_Opt, x);
    
    ScalarType beta = this->m_Opt->m_RegNorm.beta[0];
    
    VecField r(this->m_Opt);
    VecField p(this->m_Opt);
    VecField xk(this->m_Opt, x);
    VecField b(this->m_Opt, g);
    VecField m(this->m_Opt);
    
    //TwoLevelRegFFT op(this->m_Opt, beta, 0, true);
    //op.Restrict(&b, &b_f);
        
    r.Copy(&b);
    xk.SetValue(0);
    p.Copy(&r);
    
    //VecDot(r.m_X, r.m_X, &rTr);
    
    printf("%f %f %f ", this->m_Opt->m_KrylovMethod.reltol, rel_f, sqrt(rTr));
    
    rel=sqrt(rTr)*this->m_Opt->m_KrylovMethod.reltol;
    
    for (int k=0;k<50;++k) {
      this->HessianMatVec(m.m_X, p.m_X);
      VecDot(p.m_X, m.m_X, &pTAp);
      ScalarType alpha = rTr/pTAp;
      xk.AXPY(alpha, &p);
      r.AXPY(-alpha, &m);
      ScalarType rTr2;
      VecDot(r.m_X, r.m_X, &rTr2);
      if (sqrt(rTr2) < rel) {
        rTr = rTr2;
        printf("%i ", k);
        break;
      }
      ScalarType beta = rTr2/rTr;
      p.Scale(beta);
      p.AXPY(1., &r);
      rTr = rTr2;
    }
    
    printf("%f\n", sqrt(rTr));
    
    //op.Prolong(&x_f, &xk);
    }*/
    
    if (this->m_Opt->m_KrylovMethod.matvectype == PRECONDMATVECSYM) {
        ierr = this->ApplyInvRegularizationOperator(x, x, true); CHKERRQ(ierr);
    }
    
    /*VecField *g_vec = new VecField(m_Opt, x);
    
    TwoLevelFFT op(this->m_Opt);
    //op.Restrict(this->m_WorkVecField1, g_vec);
    
    this->m_WorkVecField1->Copy(g_vec);
    
    op.Prolong(g_vec, this->m_WorkVecField1);
    
    delete g_vec;*/

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies inverse of regularization operator (used as
 * spectral preconditioner for our problem)
 *******************************************************************/
PetscErrorCode CLAIREBase::ApplyInvRegularizationOperator(Vec ainvx, Vec x, bool flag) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    //ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    //ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    VecField *vinv = nullptr;
    VecField *v = nullptr;
    
    ierr = AllocateOnce(v, this->m_Opt, x); CHKERRQ(ierr);
    ierr = AllocateOnce(vinv, this->m_Opt, ainvx); CHKERRQ(ierr);
    
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    //ierr = this->m_WorkVecField1->SetComponents(x); CHKERRQ(ierr);
    
    /*FILE *handle = fopen("inreg_res.txt", "a");
    
    ScalarType norm_pre, norm_post;
    
    int nit =  this->m_Opt->GetCounter(ITERATIONS) - 1;
    int kit = this->m_Opt->GetCounter(HESSMATVEC);
    
    v->Norm2(norm_pre);
    
    fprintf(handle, "%i,%i,%e", nit, kit, sqrt(norm_pre));*/
    
    ierr = this->m_Regularization->ApplyInverse(vinv, v, flag); CHKERRQ(ierr);
    //ierr = this->SymTwoLevelHessMatVec(ainvx, x); CHKERRQ(ierr);
        
    /*vinv->Norm2(norm_post);
   
    fprintf(handle, ",%e\n", sqrt(norm_post));
    fclose(handle);*/
    
    //ierr = this->m_WorkVecField2->GetComponents(ainvx); CHKERRQ(ierr);
    
    ierr = Free(v); CHKERRQ(ierr);
    ierr = Free(vinv); CHKERRQ(ierr);
  
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief compute a synthetic test problem
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupSyntheticProb(Vec &mR, Vec &mT, VecField* v) {
    PetscErrorCode ierr = 0;
    IntType nx[3], i;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL, *p_mt = NULL;
    ScalarType hx[3], xc1, xc2, xc3, x, sigma, maxval, minval, nvx1, nvx2, nvx3;
    ScalarType x1, x2, x3, v0 = 0.5;
    int vcase = 0, icase = 0;
    bool velocityallocated = false;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    IntType nl = this->m_Opt->m_Domain.nl;
    IntType nc = this->m_Opt->m_Domain.nc;
    IntType ng = this->m_Opt->m_Domain.ng;

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg2("setting up synthetic problem"); CHKERRQ(ierr);
    }
    
    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->m_Domain.hx[i];
        nx[i] = this->m_Opt->m_Domain.nx[i];
    }

    // allocate vector fields
    
    if (this->m_VelocityField == NULL) {
        ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
        velocityallocated = true;
    }


    // allocate reference image
    if (mR == NULL) {
        ierr = VecCreate(mR, nl*nc, ng*nc); CHKERRQ(ierr);
    }
    ierr = VecSet(mR, 0); CHKERRQ(ierr);

    // allocate template image
    if (mT == NULL) {
        ierr = VecCreate(mT, nl*nc, ng*nc); CHKERRQ(ierr);
    }
    ierr = VecSet(mT, 0); CHKERRQ(ierr);

    // compute center coordinates
    xc1 = hx[0]*static_cast<ScalarType>(nx[0])/2.0;
    xc2 = hx[1]*static_cast<ScalarType>(nx[1])/2.0;
    xc3 = hx[2]*static_cast<ScalarType>(nx[2])/2.0;

    // switch use cases
    switch (this->m_Opt->m_RegFlags.synprobid) {
        case 0:
            vcase = 0; icase = 0;
            break;
        case 1:
            vcase = 1; icase = 0;
            break;
        case 2:
            vcase = 2; icase = 0;
            break;
        case 3:
            vcase = 3; icase = 0;
            break;
        case 4:
            vcase = 0; icase = 0; v0 = 1.0;
            break;
        case 5:
            vcase = 4; icase = 0; v0 = 1.0;
            break;
        default:
            vcase = 0; icase = 0;
            break;
    }

    /// for stokes we are going to use an incompressible velocity
    if (this->m_Opt->m_RegModel == STOKES) {vcase = 3;}

    ierr = VecGetArray(this->m_VelocityField->m_X1, &p_vx1); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_VelocityField->m_X2, &p_vx2); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_VelocityField->m_X3, &p_vx3); CHKERRQ(ierr);
    ierr = VecGetArray(mT, &p_mt); CHKERRQ(ierr);
//#pragma omp parallel
//{
    IntType i1, i2, i3;
//#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->m_Domain.isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < this->m_Opt->m_Domain.isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < this->m_Opt->m_Domain.isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + this->m_Opt->m_Domain.istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + this->m_Opt->m_Domain.istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + this->m_Opt->m_Domain.istart[2]);

                // compute linear / flat index
                i = GetLinearIndex(i1, i2, i3, this->m_Opt->m_Domain.isize);

                if (icase == 0) {
                    p_mt[i] =  (PetscSinReal(x1)*PetscSinReal(x1)
                              + PetscSinReal(x2)*PetscSinReal(x2)
                              + PetscSinReal(x3)*PetscSinReal(x3))/3;
                } else if (icase == 1) {
                    sigma = 0.05*2.0*PETSC_PI;
                    x1 -= xc1; x2 -= xc2; x3 -= xc3;
                    x   = PetscSqrtReal(x1*x1 + x2*x2 + x3*x3)/sigma;
                    p_mt[i] = PetscExpReal(-x*x);
                } else if (icase == 2) {
                    p_mt[i] = 0.5;
                }

                // compute the velocity field
                if (vcase == 0) {
                    p_vx1[i] = v0*PetscSinReal(x3)*PetscCosReal(x2)*PetscSinReal(x2);
                    p_vx2[i] = v0*PetscSinReal(x1)*PetscCosReal(x3)*PetscSinReal(x3);
                    p_vx3[i] = v0*PetscSinReal(x2)*PetscCosReal(x1)*PetscSinReal(x1);
                } else if (vcase == 1) {
                    p_vx1[i] = PetscSinReal(2*x1)*PetscCosReal(2*x2)*PetscSinReal(2*x3);
                    p_vx2[i] = PetscSinReal(2*x1)*PetscCosReal(2*x2)*PetscSinReal(2*x3);
                    p_vx3[i] = PetscSinReal(2*x1)*PetscCosReal(2*x2)*PetscSinReal(2*x3);
                } else if (vcase == 2) {
                    p_vx1[i] = PetscCosReal(x1)*PetscSinReal(x2);
                    p_vx2[i] = PetscCosReal(x2)*PetscSinReal(x1);
                    p_vx3[i] = PetscCosReal(x1)*PetscSinReal(x3);
                } else if (vcase == 3) {
                    p_vx1[i] = PetscCosReal(x2)*PetscCosReal(x3);
                    p_vx2[i] = PetscSinReal(x3)*PetscSinReal(x1);
                    p_vx3[i] = PetscCosReal(x1)*PetscCosReal(x2);
                } else if (vcase == 4) {
                    p_vx1[i] = v0;
                    p_vx2[i] = v0;
                    p_vx3[i] = v0;
                }
            }  // i1
        }  // i2
    }  // i3
//}  // pragma omp parallel
    ierr = VecRestoreArray(mT, &p_mt); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_VelocityField->m_X1, &p_vx1); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_VelocityField->m_X2, &p_vx2); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_VelocityField->m_X3, &p_vx3); CHKERRQ(ierr);
    
    if (v && !this->m_Opt->m_RegFlags.zeroinit) {
      ierr = v->Copy(this->m_VelocityField); CHKERRQ(ierr);
    }

    if (this->m_Opt->m_Verbosity > 2) {
        ss  << "synthetic case mt: " <<  icase << ", v: " << vcase;
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
        
        ierr = this->m_VelocityField->Norm(nvx1, nvx2, nvx3); CHKERRQ(ierr);
        ss  << "velocity norm: (" << std::scientific
            << nvx1 << "," << nvx2 << "," << nvx3 <<")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = VecMin(this->m_VelocityField->m_X1, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(this->m_VelocityField->m_X1, NULL, &maxval); CHKERRQ(ierr);
        ss << "velocity x1: (" << minval << "," << maxval <<")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = VecMin(this->m_VelocityField->m_X2, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(this->m_VelocityField->m_X2, NULL, &maxval); CHKERRQ(ierr);
        ss << "velocity x2: (" << minval << "," << maxval <<")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = VecMin(this->m_VelocityField->m_X3, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(this->m_VelocityField->m_X3, NULL, &maxval); CHKERRQ(ierr);
        ss << "velocity x3: (" << minval << "," << maxval <<")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // if the image has more than one component, just copy the
    // content of first image to all other
    if (nc == 2) {
        ierr = VecGetArray(mT, &p_mt); CHKERRQ(ierr);
        for (IntType i = 0; i < nl; ++i) {
            p_mt[nl + i] = 1.0 - p_mt[i];
        }
        ierr = VecRestoreArray(mT, &p_mt); CHKERRQ(ierr);
    }
    if (nc > 2) {
        ierr = VecGetArray(mT, &p_mt); CHKERRQ(ierr);
        for (IntType k = 1; k < nc; ++k) {
            try {std::copy(p_mt, p_mt+nl, p_mt+k*nl);}
            catch (std::exception&) {
                ierr = ThrowError("copy failed"); CHKERRQ(ierr);
            }
        }
        ierr = VecRestoreArray(mT, &p_mt); CHKERRQ(ierr);
    }
    ierr = Rescale(mT, 0, 1, nc); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = VecMin(mT, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(mT, NULL, &maxval); CHKERRQ(ierr);
        ss << "template image: (" << minval << "," << maxval <<")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // solve the forward problem using the computed
    // template image and the computed velocity field as input
    ierr = this->SolveForwardProblem(mR, mT); CHKERRQ(ierr);
    ierr = Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = Rescale(mR, 0, 1, nc); CHKERRQ(ierr);
    if (this->m_Opt->m_Verbosity > 2) {
        ierr = VecMin(mR, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(mR, NULL, &maxval); CHKERRQ(ierr);
        ss << "reference image: (" << minval << "," << maxval <<")";
        ierr = DbgMsg2(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // reset velocity field (to avoid memory leak)
    if (velocityallocated) {
        delete this->m_VelocityField;
        this->m_VelocityField = NULL;
    }

    // reset all the clocks we have used so far
    ierr = this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr = this->m_Opt->ResetCounters(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief copies some input data field to all time points
 *******************************************************************/
PetscErrorCode CLAIREBase::CopyToAllTimePoints(Vec u, Vec uj) {
    PetscErrorCode ierr = 0;
    ScalarType *p_u = NULL, *p_uj = NULL;
    IntType nl, nc, nt;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    // get pointers
    ierr = VecGetArray(u, &p_u); CHKERRQ(ierr);      ///< vec for entire time horizon
    ierr = VecGetArray(uj, &p_uj); CHKERRQ(ierr);    ///< vec at single point in time

    ierr = DebugGPUNotImplemented(); CHKERRQ(ierr);

    // for all time points
    for (IntType j = 0; j <= nt; ++j) {
        // copy data to all time points
        try {std::copy(p_uj, p_uj+nl*nc, p_u+j*nl*nc);}
        catch (std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
    }

    // restore pointers
    ierr = VecRestoreArray(u, &p_u); CHKERRQ(ierr);
    ierr = VecRestoreArray(uj, &p_uj); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief copies some input data field to all time points
 *******************************************************************/
PetscErrorCode CLAIREBase::ComputeCFLCondition() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    ScalarType hx[3], cflnum, vmax, vmaxscaled;
    IntType nt, ntcfl;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    
    
    ierr = AllocateOnce(this->m_WorkVecField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt, this->m_WorkVecField4->m_X1); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt, this->m_WorkVecField4->m_X2); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt, this->m_WorkVecField4->m_X3); CHKERRQ(ierr);
    //if (this->m_WorkScaField1 == NULL) {
    //    ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    //}

    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->m_Domain.hx[i];
    }

    ierr = VecAbs(this->m_WorkVecField1->m_X1); CHKERRQ(ierr);
    ierr = VecAbs(this->m_WorkVecField1->m_X2); CHKERRQ(ierr);
    ierr = VecAbs(this->m_WorkVecField1->m_X3); CHKERRQ(ierr);

    // compute max( |v_1| + |v_2| + |v_3| )
    ierr = VecSet(*this->m_WorkScaField1, 0.0); CHKERRQ(ierr);
    ierr = VecAXPY(*this->m_WorkScaField1, 1.0, this->m_WorkVecField1->m_X1);
    ierr = VecAXPY(*this->m_WorkScaField1, 1.0, this->m_WorkVecField1->m_X2);
    ierr = VecAXPY(*this->m_WorkScaField1, 1.0, this->m_WorkVecField1->m_X3);
    ierr = VecMax(*this->m_WorkScaField1, NULL, &vmax); CHKERRQ(ierr);

    // compute max( |v_1|/hx1 + |v_2|/hx2 + |v_3|/hx3 )
    ierr = VecSet(*this->m_WorkScaField1, 0.0); CHKERRQ(ierr);
    ierr = VecAXPY(*this->m_WorkScaField1, 1.0/hx[0], this->m_WorkVecField1->m_X1);
    ierr = VecAXPY(*this->m_WorkScaField1, 1.0/hx[1], this->m_WorkVecField1->m_X2);
    ierr = VecAXPY(*this->m_WorkScaField1, 1.0/hx[2], this->m_WorkVecField1->m_X3);
    ierr = VecMax(*this->m_WorkScaField1, NULL, &vmaxscaled); CHKERRQ(ierr);

    // if we have a zero velocity field, we do not have to worry
    ntcfl = nt;
    if (vmaxscaled != 0.0) {
        cflnum  = this->m_Opt->m_PDESolver.cflnumber;
        // compute min number of time steps
        ntcfl  = static_cast<IntType>(ceil(vmaxscaled/cflnum));
        cflnum = vmaxscaled/static_cast<ScalarType>(nt);
    }

    if (this->m_Opt->m_PDESolver.monitorcflnumber) {
        ss << "||v||_infty = " << std::scientific
           << vmax << " (cflnum = " << cflnum
           << " -> nt = " << std::setw(3) << ntcfl << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    if (this->m_Opt->m_PDESolver.adapttimestep) {
        if (ntcfl > nt) {
            if (this->m_Opt->m_Verbosity > 1) {
                ss << "changing time step: " << nt << " -> " << ntcfl;
                ierr = DbgMsg1(ss.str()); CHKERRQ(ierr);
                ss.str(std::string()); ss.clear();
            }
            // reset variables
            ierr = this->ClearVariables(); CHKERRQ(ierr);
            this->m_Opt->m_Domain.nt = ntcfl;
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
PetscErrorCode CLAIREBase::CheckBounds(Vec v, bool& boundreached) {
    PetscErrorCode ierr = 0;
    ScalarType detdgradmin, detdgradmax, bound;
    bool minboundreached, maxboundreached;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);

    minboundreached = false;
    maxboundreached = false;

    // parse input velocity field
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // compute determinant of deformation gradient
    ierr = this->ComputeDetDefGrad(); CHKERRQ(ierr);

    detdgradmin = this->m_Opt->m_Monitor.detdgradmin;
    detdgradmax = this->m_Opt->m_Monitor.detdgradmax;
    bound       = this->m_Opt->m_Monitor.detdgradbound;

    // check if jmin < bound and 1/jmax < bound
    minboundreached = detdgradmin     <= bound;
    maxboundreached = 1.0/detdgradmax <= bound;

    boundreached = (minboundreached || maxboundreached) ? true : false;
    if (boundreached) {
        ierr = WrngMsg("jacobian bound reached"); CHKERRQ(ierr);
    }

    // display what's going on
    if (this->m_Opt->m_Verbosity > 1) {
        if (minboundreached) {
            ss << std::scientific
            << "min(det(grad(y^{-1}))) = " << detdgradmin << " <= " << bound;
            ierr = DbgMsg1(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        if (maxboundreached) {
            ss << std::scientific
            << "max(det(grad(y^{-1}))) = " << detdgradmax << " >= " << 1.0/bound
            << " ( = 1/bound )";
            ierr = DbgMsg1(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
PetscErrorCode CLAIREBase::ComputeDetDefGrad(bool write2file, Vec detj) {
    PetscErrorCode ierr = 0;
    std::string filename, detstr;
    std::stringstream ss, ssnum;
    bool inverse = false;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_VelocityField == NULL) {
      ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
      ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    // check cfl condition / update time step
    if (this->m_Opt->m_PDESolver.adapttimestep) {
        ierr = this->ComputeCFLCondition(); CHKERRQ(ierr);
    }

    if (this->m_DeformationFields == NULL) {
        ierr = this->SetupDeformationField(); CHKERRQ(ierr);
    }

    ierr = VecSet(*this->m_WorkScaField1, 1.0); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero == false) {
        ierr = this->m_DeformationFields->ComputeDetDefGrad(); CHKERRQ(ierr);
    }

    if (write2file) {
        filename = inverse ? "inverse-det-deformation-grad" : "det-deformation-grad";
        filename += this->m_Opt->m_FileNames.extension;
        ierr = this->m_ReadWrite->Write(*this->m_WorkScaField1, filename); CHKERRQ(ierr);
    }

    if (detj != NULL) {
        ierr = VecCopy(*this->m_WorkScaField1, detj); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
PetscErrorCode CLAIREBase::ComputeDeformationMaps(bool write2file) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_DeformationFields == NULL) {
        ierr = this->SetupDeformationField(); CHKERRQ(ierr);
    }

    if (this->m_Opt->m_ReadWriteFlags.detdefgrad) {
        ierr = Msg("computing determinant of deformation gradient"); CHKERRQ(ierr);
        ierr = this->ComputeDetDefGrad(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->m_ReadWriteFlags.defgrad) {
        ierr = Msg("computing deformation gradient"); CHKERRQ(ierr);
        ierr = this->m_DeformationFields->ComputeDefGrad(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->m_ReadWriteFlags.defmap) {
        ierr = Msg("computing deformation map"); CHKERRQ(ierr);
        ierr = this->m_DeformationFields->ComputeDeformationMap(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->m_ReadWriteFlags.deffield) {
        ierr = Msg("computing displacement field"); CHKERRQ(ierr);
        ierr = this->m_DeformationFields->ComputeDisplacementField(true); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _CLAIREBASE_CPP_
