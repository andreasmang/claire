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

    // objects
    this->m_ReadWrite = NULL;               ///< read / write object
    this->m_Regularization = NULL;          ///< pointer for regularization class
    this->m_DistanceMeasure = NULL;         ///< distance measure
    this->m_SemiLagrangianMethod = NULL;    ///< semi lagranigan
    this->m_DeformationFields = NULL;       ///< interface for computing deformation field (jacobian; mapping; ...)
    this->m_Differentiation = NULL;         ///< interface for differentiation
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
    Free(this->m_Differentiation);
    Free(this->m_TransportProblem);
    
    Free(this->m_TemplateImage);
    Free(this->m_ReferenceImage);
    Free(this->m_AuxVariable);
    Free(this->m_CellDensity);
    Free(this->m_Mask);
    
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

    if (this->m_x1hat != NULL) {
        //accfft_free(this->m_x1hat);
        cudaFree(this->m_x1hat);
        this->m_x1hat = NULL;
    }
    if (this->m_x2hat != NULL) {
        //accfft_free(this->m_x2hat);
        cudaFree(this->m_x2hat);
        this->m_x2hat = NULL;
    }
    if (this->m_x3hat != NULL) {
        //accfft_free(this->m_x3hat);
        cudaFree(this->m_x3hat);
        this->m_x3hat = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupSpectralData() {
    PetscErrorCode ierr = 0;
    IntType nalloc;
    PetscFunctionBegin;

    nalloc = this->m_Opt->m_FFT.nalloc;
    
    if (this->m_x1hat == NULL) {
        //this->m_x1hat = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
        cudaMalloc(reinterpret_cast<void**>(&this->m_x1hat), nalloc);
    }
    if (this->m_x2hat == NULL) {
        //this->m_x2hat = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
        cudaMalloc(reinterpret_cast<void**>(&this->m_x2hat), nalloc);
    }
    if (this->m_x3hat == NULL) {
        //this->m_x3hat = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
        cudaMalloc(reinterpret_cast<void**>(&this->m_x3hat), nalloc);
    }

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
  
    ierr = AllocateOnce(this->m_ReferenceImage, this->m_Opt, mR, true); CHKERRQ(ierr);

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
    
    ierr = AllocateOnce(this->m_TemplateImage, this->m_Opt, mT, true); CHKERRQ(ierr);

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
    mask = *this->m_Mask;

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

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
//    ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_VelocityField->Copy(v); CHKERRQ(ierr);

//    this->m_VelocityField = v;
//    this->m_DeleteControlVariable = false;

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

    ierr = AllocateOnce(this->m_Mask, this->m_Opt, mask); CHKERRQ(ierr);
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
    
    
    ierr = AllocateOnce(this->m_AuxVariable, this->m_Opt, q, true); CHKERRQ(ierr);
    
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
    
    ierr = AllocateOnce(this->m_CellDensity, this->m_Opt, c, true); CHKERRQ(ierr);
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
        nl = this->m_Opt->m_Domain.nl;
        ng = this->m_Opt->m_Domain.ng;
        ierr = VecCreate(v, nl, ng); CHKERRQ(ierr);
        ierr = VecCreate(g, nl, ng); CHKERRQ(ierr);

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
            ierr = DbgMsg("zero velocity field"); CHKERRQ(ierr);
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
            ierr = Allocate<DistanceMeasureSL2>(this->m_DistanceMeasure, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case SL2AUX:
        {
            ierr = Allocate<DistanceMeasureSL2aux>(this->m_DistanceMeasure, this->m_Opt); CHKERRQ(ierr);
            // TODO: Fix for 2 level preconditioner (these need to be set in preconditioner)
            ierr = Assert(this->m_CellDensity != NULL, "null pointer (aux 1)"); CHKERRQ(ierr);
            ierr = Assert(this->m_AuxVariable != NULL, "null pointer (aux 2)"); CHKERRQ(ierr);
            ierr = this->m_DistanceMeasure->SetAuxVariable(this->m_CellDensity,1); CHKERRQ(ierr);
            ierr = this->m_DistanceMeasure->SetAuxVariable(this->m_AuxVariable,2); CHKERRQ(ierr);
            break;
        }
        case NCC:
        {
            ierr = Allocate<DistanceMeasureNCC>(this->m_DistanceMeasure, this->m_Opt); CHKERRQ(ierr);
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
            ierr = DbgMsg("distance measure: parsing template image"); CHKERRQ(ierr);
        }
        ierr = this->m_DistanceMeasure->SetTemplateImage(this->m_TemplateImage); CHKERRQ(ierr);
    }
    if (this->m_ReferenceImage != NULL) {
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = DbgMsg("distance measure: parsing reference image"); CHKERRQ(ierr);
        }
       ierr = this->m_DistanceMeasure->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
    }
    if (this->m_Mask != NULL) {
        if (this->m_Opt->m_Verbosity > 1) {
            ierr = DbgMsg("distance measure: mask enabled"); CHKERRQ(ierr);
        }
        ierr = this->m_DistanceMeasure->SetMask(this->m_Mask); CHKERRQ(ierr);
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
              ierr = DbgMsg("allocate L2 regularization"); CHKERRQ(ierr);
            }
            ierr = Allocate<RegularizationL2>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H1:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg("allocate H1 regularization"); CHKERRQ(ierr);
            }
            ierr = Allocate<RegularizationH1>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H2:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg("allocate H2 regularization"); CHKERRQ(ierr);
            }
            ierr = Allocate<RegularizationH2>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H3:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg("allocate H3 regularization"); CHKERRQ(ierr);
            }
            ierr = Allocate<RegularizationH3>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H1SN:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg("allocate H1SN regularization"); CHKERRQ(ierr);
            }
            ierr = Allocate<RegularizationH1SN>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H2SN:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg("allocate H2SN regularization"); CHKERRQ(ierr);
            }
            ierr = Allocate<RegularizationH2SN>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case H3SN:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg("allocate H3SN regularization"); CHKERRQ(ierr);
            }
            ierr = Allocate<RegularizationH3SN>(this->m_Regularization, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = reg::ThrowError("regularization model not defined"); CHKERRQ(ierr);
        }
    }

    ierr = this->m_Regularization->SetDifferentiation(Differentiation::Type::Spectral); CHKERRQ(ierr);

    // set the containers for the spectral data
    ierr = this->SetupSpectralData(); CHKERRQ(ierr);
    ierr = this->m_Regularization->SetSpectralData(this->m_x1hat,
                                                   this->m_x2hat,
                                                   this->m_x3hat); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief allocate transport problem model
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupTransportProblem() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

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
              ierr = DbgMsg("allocate SL transport problem"); CHKERRQ(ierr);
            }
            ierr = Allocate<TransportEquationSL>(this->m_TransportProblem, this->m_Opt); CHKERRQ(ierr);
            break;
        }
        case RK2:
        {
            if (this->m_Opt->m_Verbosity > 1) {
              ierr = DbgMsg("allocate RK2 transport problem"); CHKERRQ(ierr);
            }
            ierr = Allocate<TransportEquationRK2>(this->m_TransportProblem, this->m_Opt); CHKERRQ(ierr);
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

    ierr = this->m_TransportProblem->SetWorkVecField(this->m_WorkVecField1, 1); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkVecField(this->m_WorkVecField2, 2); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkVecField(this->m_WorkVecField3, 3); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField5, this->m_Opt); CHKERRQ(ierr);
    
    ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField1, 1); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField2, 2); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField3, 3); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField4, 4); CHKERRQ(ierr);
    ierr = this->m_TransportProblem->SetWorkScaField(this->m_WorkScaField5, 5); CHKERRQ(ierr);
    

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief allocate regularization model
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupDeformationField() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    ierr = AllocateOnce(this->m_DeformationFields, this->m_Opt); CHKERRQ(ierr);
    
    if (this->m_Differentiation != NULL) {
      this->m_DeformationFields->SetDifferentiation(this->m_Differentiation);
    }
    
    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField3, this->m_Opt); CHKERRQ(ierr);
    
    ierr = this->m_DeformationFields->SetWorkVecField(this->m_WorkVecField1, 1); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkVecField(this->m_WorkVecField2, 2); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkVecField(this->m_WorkVecField3, 3); CHKERRQ(ierr);
    
    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField2, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField3, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField4, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkScaField5, this->m_Opt); CHKERRQ(ierr);

    ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField1, 1); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField2, 2); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField3, 3); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField4, 4); CHKERRQ(ierr);
    ierr = this->m_DeformationFields->SetWorkScaField(*this->m_WorkScaField5, 5); CHKERRQ(ierr);

    ierr = this->m_DeformationFields->SetVelocityField(this->m_VelocityField); CHKERRQ(ierr);

    ierr = AllocateOnce(this->m_SemiLagrangianMethod, this->m_Opt); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);

    ierr = this->m_DeformationFields->SetSLM(this->m_SemiLagrangianMethod); CHKERRQ(ierr);

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

    if (this->m_Opt->m_KrylovMethod.matvectype == PRECONDMATVECSYM) {
        ierr = this->ApplyInvRegularizationOperator(x, x, true); CHKERRQ(ierr);
    }

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

    ierr = AllocateOnce(this->m_WorkVecField1, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_WorkVecField2, this->m_Opt); CHKERRQ(ierr);
    
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    ierr = this->m_WorkVecField1->SetComponents(x); CHKERRQ(ierr);
    ierr = this->m_Regularization->ApplyInverse(this->m_WorkVecField2, this->m_WorkVecField1, flag); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetComponents(ainvx); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute a synthetic test problem
 *******************************************************************/
PetscErrorCode CLAIREBase::SetupSyntheticProb(Vec &mR, Vec &mT) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nx[3], i;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL, *p_mt = NULL;
    ScalarType hx[3], xc1, xc2, xc3, x, sigma, maxval, minval, nvx1, nvx2, nvx3;
    ScalarType x1, x2, x3, v0 = 0.5;
    int vcase = 0, icase = 0;
    bool velocityallocated = false;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg("setting up synthetic problem"); CHKERRQ(ierr);
    }
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->m_Domain.hx[i];
        nx[i] = this->m_Opt->m_Domain.nx[i];
    }

    // allocate vector fields
    
    if (this->m_VelocityField == NULL) {
        ierr = Allocate(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
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
        default:
            vcase = 0; icase = 0;
            break;
    }

    /// for stokes we are going to use an incompressible velocity
    if (this->m_Opt->m_RegModel == STOKES) {vcase = 3;}

    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = GetRawPointer(mT, &p_mt); CHKERRQ(ierr);
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
    ierr = RestoreRawPointer(mT, &p_mt); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = this->m_VelocityField->Norm(nvx1, nvx2, nvx3); CHKERRQ(ierr);
        ss  << "velocity norm: (" << std::scientific
            << nvx1 << "," << nvx2 << "," << nvx3 <<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = VecMin(this->m_VelocityField->m_X1, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(this->m_VelocityField->m_X1, NULL, &maxval); CHKERRQ(ierr);
        ss << "velocity x1: (" << minval << "," << maxval <<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = VecMin(this->m_VelocityField->m_X2, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(this->m_VelocityField->m_X2, NULL, &maxval); CHKERRQ(ierr);
        ss << "velocity x2: (" << minval << "," << maxval <<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = VecMin(this->m_VelocityField->m_X3, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(this->m_VelocityField->m_X3, NULL, &maxval); CHKERRQ(ierr);
        ss << "velocity x3: (" << minval << "," << maxval <<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // if the image has more than one component, just copy the
    // content of first image to all other
    if (nc == 2) {
        ierr = GetRawPointer(mT, &p_mt); CHKERRQ(ierr);
        for (IntType i = 0; i < nl; ++i) {
            p_mt[nl + i] = 1.0 - p_mt[i];
        }
        ierr = RestoreRawPointer(mT, &p_mt); CHKERRQ(ierr);
    }
    if (nc > 2) {
        ierr = GetRawPointer(mT, &p_mt); CHKERRQ(ierr);
        for (IntType k = 1; k < nc; ++k) {
            try {std::copy(p_mt, p_mt+nl, p_mt+k*nl);}
            catch (std::exception&) {
                ierr = ThrowError("copy failed"); CHKERRQ(ierr);
            }
        }
        ierr = RestoreRawPointer(mT, &p_mt); CHKERRQ(ierr);
    }
    ierr = Rescale(mT, 0, 1, nc); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = VecMin(mT, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(mT, NULL, &maxval); CHKERRQ(ierr);
        ss << "template image: (" << minval << "," << maxval <<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
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
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
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
    ierr = GetRawPointer(u, &p_u); CHKERRQ(ierr);      ///< vec for entire time horizon
    ierr = GetRawPointer(uj, &p_uj); CHKERRQ(ierr);    ///< vec at single point in time

    // for all time points
    for (IntType j = 0; j <= nt; ++j) {
        // copy data to all time points
        try {std::copy(p_uj, p_uj+nl*nc, p_u+j*nl*nc);}
        catch (std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
    }

    // restore pointers
    ierr = RestoreRawPointer(u, &p_u); CHKERRQ(ierr);
    ierr = RestoreRawPointer(uj, &p_uj); CHKERRQ(ierr);

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
    IntType nl, ng, nt, ntcfl;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    nt = this->m_Opt->m_Domain.nt;

    ierr = AllocateOnce(this->m_WorkScaField1, this->m_Opt); CHKERRQ(ierr);
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
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
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
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        if (maxboundreached) {
            ss << std::scientific
            << "max(det(grad(y^{-1}))) = " << detdgradmax << " >= " << 1.0/bound
            << " ( = 1/bound )";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
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
      ierr = Allocate(this->m_VelocityField, this->m_Opt); CHKERRQ(ierr);
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
