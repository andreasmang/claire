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

#ifndef _OPTIMALCONTROLREGISTRATIONBASE_CPP_
#define _OPTIMALCONTROLREGISTRATIONBASE_CPP_

#include "OptimalControlRegistrationBase.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
OptimalControlRegistrationBase::OptimalControlRegistrationBase() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
OptimalControlRegistrationBase::~OptimalControlRegistrationBase() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
OptimalControlRegistrationBase::OptimalControlRegistrationBase(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // pointer for container of velocity field
    this->m_VelocityField = NULL;
    this->m_IncVelocityField = NULL;

    // pointers to images (set from outside; not to be deleted)
    this->m_TemplateImage = NULL;
    this->m_ReferenceImage = NULL;

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

    this->m_WorkTenField1 = NULL;
    this->m_WorkTenField2 = NULL;
    this->m_WorkTenField3 = NULL;
    this->m_WorkTenField4 = NULL;

    this->m_x1hat = NULL;
    this->m_x2hat = NULL;
    this->m_x3hat = NULL;

    // objects
    this->m_ReadWrite = NULL;               ///< read / write object
    this->m_Regularization = NULL;          ///< pointer for regularization class
    //this->m_DistanceMeasure = NULL;         ///< distance measure
    this->m_SemiLagrangianMethod = NULL;    ///< semi lagranigan

    this->m_VelocityIsZero = false;         ///< flag: is velocity zero
    this->m_StoreTimeHistory = true;        ///< flag: store time history (needed for inversion)
    this->m_ComputeInverseDefMap = false;   ///< flag: compute inverse deformation map

    this->m_DeleteControlVariable = true;   ///< flag: clear memory for control variable
    this->m_DeleteIncControlVariable = true;   ///< flag: clear memory for incremental control variable

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_DeleteControlVariable) {
        if (this->m_VelocityField != NULL) {
            delete this->m_VelocityField;
            this->m_VelocityField = NULL;
        }
    }

    if (this->m_DeleteIncControlVariable) {
        if (this->m_IncVelocityField != NULL) {
            delete this->m_IncVelocityField;
            this->m_IncVelocityField = NULL;
        }
    }

    if (this->m_Regularization != NULL) {
        delete this->m_Regularization;
        this->m_Regularization = NULL;
    }

    if (this->m_SemiLagrangianMethod != NULL) {
        delete this->m_SemiLagrangianMethod;
        this->m_SemiLagrangianMethod = NULL;
    }

    if (this->m_WorkScaField1 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField1); CHKERRQ(ierr);
        this->m_WorkScaField1 = NULL;
    }
    if (this->m_WorkScaField2 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField2); CHKERRQ(ierr);
        this->m_WorkScaField2 = NULL;
    }
    if (this->m_WorkScaField3 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField3); CHKERRQ(ierr);
        this->m_WorkScaField3 = NULL;
    }
    if (this->m_WorkScaField4 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField4); CHKERRQ(ierr);
        this->m_WorkScaField4 = NULL;
    }
    if (this->m_WorkScaField5 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField5); CHKERRQ(ierr);
        this->m_WorkScaField5 = NULL;
    }

    if (this->m_WorkScaFieldMC != NULL) {
        ierr = VecDestroy(&this->m_WorkScaFieldMC); CHKERRQ(ierr);
        this->m_WorkScaFieldMC = NULL;
    }

    if (this->m_WorkVecField1 != NULL) {
        delete this->m_WorkVecField1;
        this->m_WorkVecField1 = NULL;
    }
    if (this->m_WorkVecField2 != NULL) {
        delete this->m_WorkVecField2;
        this->m_WorkVecField2 = NULL;
    }
    if (this->m_WorkVecField3 != NULL) {
        delete this->m_WorkVecField3;
        this->m_WorkVecField3 = NULL;
    }
    if (this->m_WorkVecField4 != NULL) {
        delete this->m_WorkVecField4;
        this->m_WorkVecField4 = NULL;
    }
    if (this->m_WorkVecField5 != NULL) {
        delete this->m_WorkVecField5;
        this->m_WorkVecField5 = NULL;
    }

    if (this->m_WorkTenField1 != NULL) {
        delete this->m_WorkTenField1;
        this->m_WorkTenField1 = NULL;
    }
    if (this->m_WorkTenField2 != NULL) {
        delete this->m_WorkTenField2;
        this->m_WorkTenField2 = NULL;
    }
    if (this->m_WorkTenField3 != NULL) {
        delete this->m_WorkTenField3;
        this->m_WorkTenField3 = NULL;
    }
    if (this->m_WorkTenField4 != NULL) {
        delete this->m_WorkTenField4;
        this->m_WorkTenField4 = NULL;
    }

    if (this->m_x1hat != NULL) {
        accfft_free(this->m_x1hat);
        this->m_x1hat = NULL;
    }
    if (this->m_x2hat != NULL) {
        accfft_free(this->m_x2hat);
        this->m_x2hat = NULL;
    }
    if (this->m_x3hat != NULL) {
        accfft_free(this->m_x3hat);
        this->m_x3hat = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::AllocateSpectralData() {
    PetscErrorCode ierr = 0;
    IntType nalloc;
    PetscFunctionBegin;

    nalloc = this->m_Opt->GetFFT().nalloc;

    if (this->m_x1hat == NULL) {
        this->m_x1hat = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
    }
    if (this->m_x2hat == NULL) {
        this->m_x2hat = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
    }
    if (this->m_x3hat == NULL) {
        this->m_x3hat = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::SetReadWrite(ReadWriteReg* readwrite) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(readwrite != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite = readwrite;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::SetReferenceImage(Vec mR) {
    PetscErrorCode ierr = 0;
    IntType nc;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);
    nc = this->m_Opt->GetDomainPara().nc;

    // by default we rescale the intensity range to [0,1]
    if (this->m_Opt->GetRegFlags().applyrescaling) {
        ierr = Rescale(mR, 0.0, 1.0, nc); CHKERRQ(ierr);
    }

    // assign pointer
    this->m_ReferenceImage = mR;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set template image (i.e., the image to be registered)
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::SetTemplateImage(Vec mT) {
    PetscErrorCode ierr = 0;
    IntType nc;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);
    nc = this->m_Opt->GetDomainPara().nc;

    // by default we rescale the intensity range to [0,1]
    if (this->m_Opt->GetRegFlags().applyrescaling) {
        ierr = Rescale(mT, 0.0, 1.0, nc); CHKERRQ(ierr);
    }

    // assign pointer
    this->m_TemplateImage = mT;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get reference image (i.e., the fixed image)
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::GetReferenceImage(Vec& mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
    mR = this->m_ReferenceImage;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get template image (i.e., the image to be registered)
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::GetTemplateImage(Vec& mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    mT = this->m_TemplateImage;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::SetControlVariable(VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
//    ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_VelocityField->Copy(v); CHKERRQ(ierr);

//    this->m_VelocityField = v;
//    this->m_DeleteControlVariable = false;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::SetIncControlVariable(VecField* v) {
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
PetscErrorCode OptimalControlRegistrationBase::GetControlVariable(VecField*& v) {
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
PetscErrorCode OptimalControlRegistrationBase::ComputeInitialGuess() {
    PetscErrorCode ierr = 0;
    Vec v = NULL, g = NULL;
    IntType nl, ng;
    PetscFunctionBegin;

    // if velocity field is null pointer, we did not set
    // any initial guess
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);

    if (!this->m_Opt->GetOptPara().usezeroinitialguess) {
        nl = this->m_Opt->GetDomainPara().nl;
        ng = this->m_Opt->GetDomainPara().ng;
        ierr = VecCreate(v, nl, ng); CHKERRQ(ierr);
        ierr = VecCreate(g, nl, ng); CHKERRQ(ierr);

        ierr = this->m_VelocityField->GetComponents(v); CHKERRQ(ierr);
        ierr = this->EvaluatePrecondGradient(g, v); CHKERRQ(ierr);

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
PetscErrorCode OptimalControlRegistrationBase::IsVelocityZero() {
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

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief allocate regularization model
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::AllocateRegularization() {
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
    switch (this->m_Opt->GetRegNorm().type) {
        case L2:
        {
            try {this->m_Regularization = new RegularizationRegistrationL2(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H1:
        {
            try {this->m_Regularization = new RegularizationRegistrationH1(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H2:
        {
            try {this->m_Regularization = new RegularizationRegistrationH2(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H3:
        {
            try {this->m_Regularization = new RegularizationRegistrationH3(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H1SN:
        {
            try {this->m_Regularization = new RegularizationRegistrationH1SN(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H2SN:
        {
            try {this->m_Regularization = new RegularizationRegistrationH2SN(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H3SN:
        {
            try {this->m_Regularization = new RegularizationRegistrationH3SN(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        default:
        {
            ierr = reg::ThrowError("regularization model not defined"); CHKERRQ(ierr);
        }
    }


    // set the containers for the spectral data
    ierr = this->AllocateSpectralData(); CHKERRQ(ierr);
    ierr = this->m_Regularization->SetSpectralData(this->m_x1hat,
                                                   this->m_x2hat,
                                                   this->m_x3hat); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute coordinates of regular grid
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase
::ComputeRegularGrid(VecField* x) {
    PetscErrorCode ierr = 0;
    ScalarType hx[3];
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);

    // get grid size
    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->GetDomainPara().hx[i];
    }

    // compute initial condition
    ierr = x->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

#pragma omp parallel
{
    IntType l,i1,i2,i3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->GetDomainPara().isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < this->m_Opt->GetDomainPara().isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < this->m_Opt->GetDomainPara().isize[2]; ++i3) {  // x3
                // compute linear / flat index
                l = GetLinearIndex(i1, i2, i3, this->m_Opt->GetDomainPara().isize);

                // compute coordinates (nodal grid)
                p_x1[l] = hx[0]*static_cast<ScalarType>(i1 + this->m_Opt->GetDomainPara().istart[0]);
                p_x2[l] = hx[1]*static_cast<ScalarType>(i2 + this->m_Opt->GetDomainPara().istart[1]);
                p_x3[l] = hx[2]*static_cast<ScalarType>(i3 + this->m_Opt->GetDomainPara().istart[2]);
            } // i1
        } // i2
    } // i3
}// pragma omp for

    ierr = x->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief evaluate regularization model
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase
::EvaluateRegularizationFunctional(ScalarType* value, VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    ierr = this->m_Regularization->EvaluateFunctional(value, v); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies inverse regularization operator
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ApplyInvRegOp(Vec Ainvx, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_WorkVecField1 == NULL) {
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_WorkVecField2 == NULL) {
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    ierr = this->m_WorkVecField1->SetComponents(x); CHKERRQ(ierr);
    ierr = this->m_Regularization->ApplyInvOp(this->m_WorkVecField2, this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetComponents(Ainvx); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief estimate eigenvalues of hessian
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase
::EstimateExtremalHessEigVals(ScalarType &emin,
                              ScalarType &emax) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
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
PetscErrorCode OptimalControlRegistrationBase
::PreKrylovSolve(Vec g, Vec x) {
    PetscErrorCode ierr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // switch between hessian operators
    switch (this->m_Opt->GetKrylovSolverPara().matvectype) {
        case DEFAULTMATVEC:
        {
            // do nothing
            break;
        }
        case PRECONDMATVEC:
        {
            ierr = this->ApplyInvRegOp(g, g); CHKERRQ(ierr);
            break;
        }
        case PRECONDMATVECSYM:
        {
            ierr = this->ApplyInvRegOpSqrt(g); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("something wrong with setup"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief post-processing after the krylov solve
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::PostKrylovSolve(Vec g, Vec x) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // switch between hessian operators
    switch (this->m_Opt->GetKrylovSolverPara().matvectype) {
        case DEFAULTMATVEC:
        {
            // do nothing
            break;
        }
        case PRECONDMATVEC:
        {
            // do nothing
            break;
        }
        case PRECONDMATVECSYM:
        {
            ierr = this->ApplyInvRegOpSqrt(x); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("something wrong with setup"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply projection operator (this is for the analytically/
 * spectrally preconditioned hessian; it applies the projection
 * operator to the search direction and the gradient)
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ApplyInvRegOpSqrt(Vec x) {
    PetscErrorCode ierr = 0;
    bool usesqrt = true;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_WorkVecField1 == NULL) {
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    // set components
    ierr = this->m_WorkVecField1->SetComponents(x); CHKERRQ(ierr);

    // apply sqrt of inverse regularization operator
    ierr = this->m_Regularization->ApplyInvOp(this->m_WorkVecField2,
                                              this->m_WorkVecField1,
                                              usesqrt); CHKERRQ(ierr);

    ierr = this->m_WorkVecField2->GetComponents(x); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute a synthetic test problem
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::SetupSyntheticProb(Vec &mR, Vec &mT) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nx[3], i;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
               *p_mt = NULL, hx[3], xc1, xc2, xc3, x,
                sigma, maxval, minval, nvx1, nvx2, nvx3;
    ScalarType x1, x2, x3;
    //int vcase = 2;
    int vcase = 0;
    int icase = 0;
    ScalarType v0 = 0.5;
    //ScalarType v0 = PETSC_PI/12.0;
    bool velocityallocated = false;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("setting up synthetic problem"); CHKERRQ(ierr);
    }
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->GetDomainPara().hx[i];
        nx[i] = this->m_Opt->GetDomainPara().nx[i];
    }

    // allocate vector fields
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
        velocityallocated = true;
    }

    if (this->m_Opt->GetRegModel() == STOKES) {vcase = 3;}

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

    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = VecGetArray(mT, &p_mt); CHKERRQ(ierr);
#pragma omp parallel
{
    IntType i1, i2, i3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->GetDomainPara().isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < this->m_Opt->GetDomainPara().isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < this->m_Opt->GetDomainPara().isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + this->m_Opt->GetDomainPara().istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + this->m_Opt->GetDomainPara().istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + this->m_Opt->GetDomainPara().istart[2]);

                // compute linear / flat index
                i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetDomainPara().isize);

                if (icase == 0) {
                    p_mt[i] =  (PetscSinReal(x1)*PetscSinReal(x1)
                              + PetscSinReal(x2)*PetscSinReal(x2)
                              + PetscSinReal(x3)*PetscSinReal(x3))/3;
                } else if (icase == 1) {
                    sigma = 0.05*2.0*PETSC_PI;
                    x1 -= xc1; x2 -= xc2; x3 -= xc3;
                    x   = PetscSqrtReal(x1*x1 + x2*x2 + x3*x3)/sigma;
                    p_mt[i] = PetscExpReal(-x*x);
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
}  // pragma omp parallel
    ierr = VecRestoreArray(mT, &p_mt); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2) {
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

    if (this->m_Opt->GetVerbosity() > 2) {
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
    if (this->m_Opt->GetVerbosity() > 2) {
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
PetscErrorCode OptimalControlRegistrationBase::CopyToAllTimePoints(Vec u, Vec uj) {
    PetscErrorCode ierr = 0;
    ScalarType *p_u = NULL, *p_uj = NULL;
    IntType nl, nc, nt;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;

    // get pointers
    ierr = VecGetArray(u, &p_u); CHKERRQ(ierr);      ///< vec for entire time horizon
    ierr = VecGetArray(uj, &p_uj); CHKERRQ(ierr);    ///< vec at single point in time

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
PetscErrorCode OptimalControlRegistrationBase::ComputeCFLCondition() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    ScalarType hx[3], cflnum, vmax, vmaxscaled;
    IntType nl, ng, nt, ntcfl;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    nt = this->m_Opt->GetDomainPara().nt;

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }

    if (this->m_WorkVecField1 == NULL) {
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->GetDomainPara().hx[i];
    }

    ierr = VecAbs(this->m_WorkVecField1->m_X1); CHKERRQ(ierr);
    ierr = VecAbs(this->m_WorkVecField1->m_X2); CHKERRQ(ierr);
    ierr = VecAbs(this->m_WorkVecField1->m_X3); CHKERRQ(ierr);

    // compute max( |v_1| + |v_2| + |v_3| )
    ierr = VecSet(this->m_WorkScaField1, 0.0); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_WorkScaField1, 1.0, this->m_WorkVecField1->m_X1);
    ierr = VecAXPY(this->m_WorkScaField1, 1.0, this->m_WorkVecField1->m_X2);
    ierr = VecAXPY(this->m_WorkScaField1, 1.0, this->m_WorkVecField1->m_X3);
    ierr = VecMax(this->m_WorkScaField1, NULL, &vmax); CHKERRQ(ierr);

    // compute max( |v_1|/hx1 + |v_2|/hx2 + |v_3|/hx3 )
    ierr = VecSet(this->m_WorkScaField1, 0.0); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_WorkScaField1, 1.0/hx[0], this->m_WorkVecField1->m_X1);
    ierr = VecAXPY(this->m_WorkScaField1, 1.0/hx[1], this->m_WorkVecField1->m_X2);
    ierr = VecAXPY(this->m_WorkScaField1, 1.0/hx[2], this->m_WorkVecField1->m_X3);
    ierr = VecMax(this->m_WorkScaField1, NULL, &vmaxscaled); CHKERRQ(ierr);

    // if we have a zero velocity field, we do not have to worry
    ntcfl = nt;
    if (vmaxscaled != 0.0) {
        cflnum  = this->m_Opt->GetPDESolverPara().cflnumber;
        // compute min number of time steps
        ntcfl  = static_cast<IntType>(ceil(vmaxscaled/cflnum));
        cflnum = vmaxscaled/static_cast<ScalarType>(nt);
    }

    if (this->m_Opt->GetPDESolverPara().monitorcflnumber) {
        ss << "||v||_infty = " << std::scientific
           << vmax << " (cflnum = " << cflnum
           << " -> nt = " << std::setw(3) << ntcfl << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    if (this->m_Opt->GetPDESolverPara().adapttimestep) {
        if (ntcfl > nt) {
            if (this->m_Opt->GetVerbosity() > 1) {
                ss << "changing time step: " << nt << " -> " << ntcfl;
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.str(std::string()); ss.clear();
            }
            // reset variables
            ierr = this->ClearVariables(); CHKERRQ(ierr);
            this->m_Opt->SetNumTimePoints(ntcfl);
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::CheckBounds(Vec v, bool& boundreached) {
    PetscErrorCode ierr = 0;
    ScalarType detdgradmin, detdgradmax, bound;
    bool minboundreached, maxboundreached;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_VelocityField == NULL) {
       try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    minboundreached = false;
    maxboundreached = false;

    // parse input velocity field
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // compute determinant of deformation gradient
    ierr = this->ComputeDetDefGrad(); CHKERRQ(ierr);

    detdgradmin = this->m_Opt->GetRegMonitor().detdgradmin;
    detdgradmax = this->m_Opt->GetRegMonitor().detdgradmax;
    bound       = this->m_Opt->GetRegMonitor().detdgradbound;

    // check if jmin < bound and 1/jmax < bound
    minboundreached = detdgradmin     <= bound;
    maxboundreached = 1.0/detdgradmax <= bound;

    boundreached = (minboundreached || maxboundreached) ? true : false;
    if (boundreached) {
        ierr = WrngMsg("jacobian bound reached"); CHKERRQ(ierr);
    }

    // display what's going on
    if (this->m_Opt->GetVerbosity() > 1) {
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
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGrad(bool write2file, Vec detj) {
    PetscErrorCode ierr = 0;
    ScalarType minddg, maxddg, meanddg;
    IntType nl, ng;
    std::string filename, detstr;
    std::stringstream ss, ssnum;
    bool inverse;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    inverse = this->m_Opt->GetRegFlags().invdefgrad;
    detstr = inverse ? "det(grad(inv(y)))" : "det(grad(y))";

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("computing " + detstr); CHKERRQ(ierr);
    }

    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    if (this->m_VelocityField == NULL) {
       try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    // check cfl condition / update time step
    if (this->m_Opt->GetPDESolverPara().adapttimestep) {
        ierr = this->ComputeCFLCondition(); CHKERRQ(ierr);
    }

    // set initial condition
    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    ierr = VecSet(this->m_WorkScaField1, 1.0); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero == false) {
        // call the solver
        if (this->m_Opt->GetRegFlags().detdefgradfromdeffield) {
            ierr = this->ComputeDetDefGradViaDispField(); CHKERRQ(ierr);
        } else {
            switch (this->m_Opt->GetPDESolverPara().type) {
                case RK2:
                {
                    ierr = this->ComputeDetDefGradRK2(); CHKERRQ(ierr);
                    break;
                }
                case RK2A:
                {
                    ierr = this->ComputeDetDefGradRK2A(); CHKERRQ(ierr);
                    break;
                }
                case SL:
                {
                    ierr = this->ComputeDetDefGradSL(); CHKERRQ(ierr);
                    break;
                }
                default:
                {
                    ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                    break;
                }
            }
        }
    }

    ierr = VecMin(this->m_WorkScaField1, NULL, &minddg); CHKERRQ(ierr);
    ierr = VecMax(this->m_WorkScaField1, NULL, &maxddg); CHKERRQ(ierr);
    ierr = VecSum(this->m_WorkScaField1, &meanddg); CHKERRQ(ierr);
    meanddg /= static_cast<ScalarType>(this->m_Opt->GetDomainPara().ng);

    // remember
    this->m_Opt->SetDetDGradMin(minddg);
    this->m_Opt->SetDetDGradMax(maxddg);
    this->m_Opt->SetDetDGradMean(meanddg);

    if (this->m_Opt->GetVerbosity() > 1 || this->m_Opt->GetRegMonitor().detdgradenabled) {
        ss  << std::scientific << detstr << " : (min, mean, max)="
            << "(" << minddg << ", " << meanddg << ", " << maxddg << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    if (write2file) {
        filename = inverse ? "inverse-det-deformation-grad" : "det-deformation-grad";
        filename += this->m_Opt->GetFileNames().extension;
        ierr = this->m_ReadWrite->Write(this->m_WorkScaField1, filename); CHKERRQ(ierr);
    }

    if (detj != NULL) {
        ierr = VecCopy(this->m_WorkScaField1, detj); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradRK2() {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nt;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL, *p_jbar = NULL,
                *p_gx1 = NULL, *p_gx2 = NULL, *p_gx3 = NULL, *p_divv = NULL,
                *p_jac = NULL,  *p_rhs0 = NULL;
    ScalarType ht, hthalf, alpha, rhs1;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timer[7] = {0};
    bool inverse;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL) {
        ierr = VecCreate(this->m_WorkScaField3, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField4 == NULL) {
        ierr = VecCreate(this->m_WorkScaField4, nl, ng); CHKERRQ(ierr);
    }

    inverse = this->m_Opt->GetRegFlags().invdefgrad;
    alpha = inverse ? -1.0 : 1.0;

    // get pointers
    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gx1, p_gx2, p_gx3); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_WorkScaField1, &p_jac); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2, &p_divv); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField3, &p_jbar); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField4, &p_rhs0); CHKERRQ(ierr);

    // compute div(v)
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(p_divv, p_vx1, p_vx2, p_vx3, this->m_Opt->GetFFT().plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    // for all time points
    for (IntType j = 0; j <= nt; ++j) {
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gx1, p_gx2, p_gx3, p_jac, this->m_Opt->GetFFT().plan, &XYZ, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            // \bar{j} = j (\idiv \vect{v}) - (\vect{v} \cdot \igrad) j
            p_rhs0[i] = alpha*(p_jac[i]*p_divv[i]) - (p_vx1[i]*p_gx1[i] + p_vx2[i]*p_gx2[i] + p_vx3[i]*p_gx3[i]);
            p_jbar[i] = p_jac[i] + ht*p_rhs0[i];
        }

        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gx1, p_gx2, p_gx3, p_jbar, this->m_Opt->GetFFT().plan, &XYZ, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            // \bar{j} = j (\idiv \vect{v}) - (\vect{v} \cdot \igrad) j
            rhs1 = alpha*(p_jbar[i]*p_divv[i]) - (p_vx1[i]*p_gx1[i] + p_vx2[i]*p_gx2[i] + p_vx3[i]*p_gx3[i]);
            p_jac[i] += hthalf*(p_rhs0[i] + rhs1);
        }
    }

    ierr = VecRestoreArray(this->m_WorkScaField4, &p_rhs0); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField3, &p_jbar); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_divv); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_jac); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(p_gx1, p_gx2, p_gx3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradRK2A() {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nt;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL, *p_phibar = NULL,
                *p_gphi1 = NULL, *p_gphi2 = NULL, *p_gphi3 = NULL, *p_divv = NULL,
                *p_phiv1 = NULL, *p_phiv2 = NULL, *p_phiv3 = NULL,
                *p_phi = NULL,  *p_rhs0 = NULL,  *p_divvphi=NULL;
    ScalarType ht, hthalf, alpha;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timer[7] = {0};
    bool inverse;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkVecField1 == NULL) {
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
       try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL) {
        ierr = VecCreate(this->m_WorkScaField3, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField4 == NULL) {
        ierr = VecCreate(this->m_WorkScaField4, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField5 == NULL) {
        ierr = VecCreate(this->m_WorkScaField5, nl, ng); CHKERRQ(ierr);
    }

    // set initial condition
    ierr = VecSet(this->m_WorkScaField1, 1.0); CHKERRQ(ierr);

    // get pointers
    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gphi1, p_gphi2, p_gphi3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_phiv1, p_phiv2, p_phiv3); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_WorkScaField1, &p_phi); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2, &p_divv); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField3, &p_phibar); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField4, &p_rhs0); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField5, &p_divvphi); CHKERRQ(ierr);

    // compute div(v)
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(p_divv, p_vx1, p_vx2, p_vx3, this->m_Opt->GetFFT().plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTGRAD);

    inverse = this->m_Opt->GetRegFlags().invdefgrad;
    alpha = inverse ? -1.0 : 1.0;

#pragma omp parallel
{
#pragma omp  for
    for (IntType i=0; i < nl; ++i) { // for all grid points
        // compute phi \vect{v} = 1 \vect{v}
        p_phiv1[i] = p_vx1[i];
        p_phiv2[i] = p_vx2[i];
        p_phiv3[i] = p_vx3[i];

    }
} // pragma omp


    // for all time points
    for (IntType j = 0; j <= nt; ++j) {
        // compute grad(\phi_j)
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gphi1, p_gphi2, p_gphi3, p_phi, this->m_Opt->GetFFT().plan, &XYZ, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

        // compute div(\vect{v}\phi_j)
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_divergence_t(p_divvphi, p_phiv1, p_phiv2, p_phiv3, this->m_Opt->GetFFT().plan, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTDIV);
#pragma omp parallel
{
#pragma omp  for
        for (IntType i=0; i < nl; ++i) {  // for all grid points
            // \bar{j} = -(\vect{v} \cdot \igrad) \phi + j (\idiv \vect{v})
            //p_rhs0[i] = - alpha*(p_vx1[i]*p_gphi1[i] + p_vx2[i]*p_gphi2[i] + p_vx3[i]*p_gphi3[i])
            //            + 0.5*alpha*p_phi[i]*p_divv[i] + 0.5*alpha*p_divvphi[i]
            //            - 0.5*alpha*(p_gphi1[i]*p_vx1[i] + p_gphi2[i]*p_vx2[i] + p_gphi3[i]*p_vx3[i]);
            p_rhs0[i] = - 0.5*(p_vx1[i]*p_gphi1[i] + p_vx2[i]*p_gphi2[i] + p_vx3[i]*p_gphi3[i])
                        + alpha*p_phi[i]*p_divv[i] - 0.5*p_divvphi[i] + 0.5*p_phi[i]*p_divv[i];

            p_phibar[i] = p_phi[i] + ht*p_rhs0[i];

            // compute \bar{phi} \vect{v}
            p_phiv1[i] = p_phibar[i]*p_vx1[i];
            p_phiv2[i] = p_phibar[i]*p_vx2[i];
            p_phiv3[i] = p_phibar[i]*p_vx3[i];
        }
} // pragma omp

        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gphi1, p_gphi2, p_gphi3, p_phibar, this->m_Opt->GetFFT().plan, &XYZ, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

        // compute div(\vect{v}\bar{\phi}_j)
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_divergence_t(p_divvphi, p_phiv1, p_phiv2, p_phiv3, this->m_Opt->GetFFT().plan, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTDIV);
#pragma omp parallel
{
        ScalarType rhs1;
#pragma omp  for
        for (IntType i=0; i < nl; ++i) {  // for all grid points
            // \bar{j} = -(\vect{v} \cdot \igrad) j + j (\idiv \vect{v})
            //rhs1 = -alpha*(p_vx1[i]*p_gphi1[i] + p_vx2[i]*p_gphi2[i] + p_vx3[i]*p_gphi3[i])
            //        + 0.5*alpha*p_phibar[i]*p_divv[i] + 0.5*alpha*p_divvphi[i]
            //        - 0.5*alpha*(p_gphi1[i]*p_vx1[i] + p_gphi2[i]*p_vx2[i] + p_gphi3[i]*p_vx3[i]);
            rhs1 = - 0.5*(p_vx1[i]*p_gphi1[i] + p_vx2[i]*p_gphi2[i] + p_vx3[i]*p_gphi3[i])
                   + alpha*p_phibar[i]*p_divv[i] - 0.5*p_divvphi[i] + 0.5*p_phibar[i]*p_divv[i];

            p_phi[i] = p_phi[i] + hthalf*(p_rhs0[i] + rhs1);

            // compute \phi \vect{v} for next time step
            p_phiv1[i] = p_phi[i]*p_vx1[i];
            p_phiv2[i] = p_phi[i]*p_vx2[i];
            p_phiv3[i] = p_phi[i]*p_vx3[i];
        }
} // pragma omp
    }

    ierr = VecRestoreArray(this->m_WorkScaField1, &p_phi); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_divv); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField3, &p_phibar); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField4, &p_rhs0); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField5, &p_divvphi); CHKERRQ(ierr);

    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_gphi1, p_gphi2, p_gphi3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_phiv1, p_phiv2, p_phiv3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradSL() {
    PetscErrorCode ierr;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_divv = NULL, *p_divvX = NULL, *p_jac = NULL, *p_jacX=NULL;
    ScalarType ht, hthalf, alpha;
    IntType nl, ng, nt;
    std::stringstream ss;
    std::string ext;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timer[7] = {0};
    bool inverse;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ext = this->m_Opt->GetFileNames().extension;

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_SemiLagrangianMethod == NULL) {
        try{this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL) {
        ierr = VecCreate(this->m_WorkScaField3, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField4 == NULL) {
        ierr = VecCreate(this->m_WorkScaField4, nl, ng); CHKERRQ(ierr);
    }

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);

    // store time series
    if (this->m_Opt->GetReadWriteFlags().timeseries) {
        ss.str(std::string()); ss.clear();
        ss << "det-deformation-grad-j=" << std::setw(3) << std::setfill('0') << 0 << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkScaField1, ss.str()); CHKERRQ(ierr);
    }

    // get pointers
    ierr = VecGetArray(this->m_WorkScaField1, &p_jac); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2, &p_jacX); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField3, &p_divv); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField4, &p_divvX); CHKERRQ(ierr);

    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    inverse = this->m_Opt->GetRegFlags().invdefgrad;
    alpha = inverse ? -1.0 : 1.0;

    // compute div(v)
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(p_divv, p_vx1, p_vx2, p_vx3, this->m_Opt->GetFFT().plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    // compute div(v) at X
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_divvX, p_divv, "state"); CHKERRQ(ierr);


    for (IntType j = 0; j < nt; ++j) {  // for all time points
        // compute J(X,t^j)
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_jacX, p_jac, "state"); CHKERRQ(ierr);

#pragma omp parallel
{
        ScalarType jacX, rhs0, rhs1;
#pragma omp  for
        for (IntType i = 0; i < nl; ++i) { // for all grid points
            jacX = p_jacX[i];
            rhs0 = alpha*jacX*p_divvX[i];
            rhs1 = alpha*(jacX + ht*rhs0)*p_divv[i];
            p_jac[i] = jacX + hthalf*(rhs0 + rhs1);
        }
}  // pragma omp

        // store time series
        if (this->m_Opt->GetReadWriteFlags().timeseries) {
            ierr = VecRestoreArray(this->m_WorkScaField1, &p_jac); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
            ss << "det-deformation-grad-j=" << std::setw(3) << std::setfill('0') << j+1 << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkScaField1, ss.str()); CHKERRQ(ierr);
            ierr = VecGetArray(this->m_WorkScaField1, &p_jac); CHKERRQ(ierr);
        }
    }  // for all time points

    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_WorkScaField4, &p_divvX); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField3, &p_divv); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_jacX); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_jac); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}



/********************************************************************
 * @brief compute determinant of deformation gradient; this
 * implementation first computes the deformation map and then
 * evaluates the determinant of the deformation gradient
 * based on the computed deformation map
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradViaDispField() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    ScalarType  *p_u1 = NULL, *p_u2 = NULL, *p_u3 = NULL, *p_phi = NULL,
                *p_gu11 = NULL, *p_gu12 = NULL, *p_gu13 = NULL,
                *p_gu21 = NULL, *p_gu22 = NULL, *p_gu23 = NULL,
                *p_gu31 = NULL, *p_gu32 = NULL, *p_gu33 = NULL;
    double timer[7] = {0};
    std::bitset<3>XYZ = 0; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    if (this->m_WorkVecField1 == NULL) {
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
       try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3 == NULL) {
       try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField4 == NULL) {
       try{this->m_WorkVecField4 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1,nl,ng); CHKERRQ(ierr);
    }

    // compute deformation map (stored in work vec field one)
    ierr = this->ComputeDisplacementField(); CHKERRQ(ierr);

    // compute the derivatives (jacobian matrix; deformation gradient)
    ierr = this->m_WorkVecField1->GetArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_gu11, p_gu12, p_gu13); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->GetArrays(p_gu21, p_gu22, p_gu23); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->GetArrays(p_gu31, p_gu32, p_gu33); CHKERRQ(ierr);

    // compute gradient of components of displacement field
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_grad_t(p_gu11, p_gu12, p_gu13, p_u1,this->m_Opt->GetFFT().plan,&XYZ,timer);
    accfft_grad_t(p_gu21, p_gu22, p_gu23, p_u2,this->m_Opt->GetFFT().plan,&XYZ,timer);
    accfft_grad_t(p_gu31, p_gu32, p_gu33, p_u3,this->m_Opt->GetFFT().plan,&XYZ,timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 3*FFTGRAD);

    ierr = VecGetArray(this->m_WorkScaField1, &p_phi); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp  for
    for (IntType i = 0; i < nl; ++i) {  // for all grid points
        // compute determinant of deformation gradient
        p_phi[i] = (1.0-p_gu11[i])*(1.0-p_gu22[i])*(1.0-p_gu33[i])
                 + p_gu12[i]*p_gu23[i]*p_gu31[i]
                 + p_gu13[i]*p_gu21[i]*p_gu32[i]
                 - p_gu13[i]*(1.0-p_gu22[i])*p_gu31[i]
                 - p_gu12[i]*p_gu21[i]*(1.0-p_gu33[i])
                 - (1.0-p_gu11[i])*p_gu23[i]*p_gu32[i];
    }
}  // pragma omp

    ierr = VecRestoreArray(this->m_WorkScaField1, &p_phi); CHKERRQ(ierr);

    ierr = this->m_WorkVecField4->RestoreArrays(p_gu31, p_gu32, p_gu33); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_gu21, p_gu22, p_gu23); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_gu11, p_gu12, p_gu13); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief compute deformation gradient
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDefGrad(bool write2file) {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    std::string ext;
    ScalarType *p_phi = NULL, *p_j11 = NULL, *p_j12 = NULL, *p_j13 = NULL,
                *p_j21 = NULL, *p_j22 = NULL, *p_j23 = NULL,
                *p_j31 = NULL, *p_j32 = NULL, *p_j33 = NULL;
    ScalarType minj, meanj, maxj;
    std::stringstream ss, ssnum;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("computing deformation gradient"); CHKERRQ(ierr);
    }

    if (this->m_VelocityField == NULL) {
       try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    // allocate tensor field and set initial condition
    if (this->m_WorkTenField1 == NULL) {
       try {this->m_WorkTenField1 = new TenField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    ierr = this->m_WorkTenField1->SetIdentity(); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1,nl,ng); CHKERRQ(ierr);
    }

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero == false) {
        // call the solver
        if (this->m_Opt->GetRegFlags().detdefgradfromdeffield) {
            ierr = this->ComputeDetDefGradViaDispField(); CHKERRQ(ierr);
        } else {
            switch (this->m_Opt->GetPDESolverPara().type) {
                case RK2:
                {
                    ierr = ThrowError("not implemented");CHKERRQ(ierr);
                    //ierr = this->ComputeDefGradRK2(); CHKERRQ(ierr);
                    break;
                }
                case RK2A:
                {
                    ierr = ThrowError("not implemented");CHKERRQ(ierr);
                    //ierr = this->ComputeDefGradRK2A(); CHKERRQ(ierr);
                    break;
                }
                case SL:
                {
                    ierr = this->ComputeDefGradSL(); CHKERRQ(ierr);
                    break;
                }
                default:
                {
                    ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                    break;
                }
            }
        }
    }

    ierr = VecGetArray(this->m_WorkScaField1,&p_phi); CHKERRQ(ierr);
    ierr = this->m_WorkTenField1->GetArrays(p_j11, p_j12, p_j13,
                                            p_j21, p_j22, p_j23,
                                            p_j31, p_j32, p_j33); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp  for
    for (IntType i=0; i < nl; ++i) { // for all grid points
        // compute determinant of deformation gradient
        p_phi[i] = p_j11[i]*p_j22[i]*p_j33[i]
                 + p_j12[i]*p_j23[i]*p_j31[i]
                 + p_j13[i]*p_j21[i]*p_j32[i]
                 - p_j13[i]*p_j22[i]*p_j31[i]
                 - p_j12[i]*p_j21[i]*p_j33[i]
                 - p_j11[i]*p_j23[i]*p_j32[i];
    }

} // pragma omp

    ierr = this->m_WorkTenField1->RestoreArrays(p_j11, p_j12, p_j13,
                                                p_j21, p_j22, p_j23,
                                                p_j31, p_j32, p_j33); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_phi); CHKERRQ(ierr);

    ierr = VecMin(this->m_WorkScaField1, NULL, &minj); CHKERRQ(ierr);
    ierr = VecMax(this->m_WorkScaField1, NULL, &maxj); CHKERRQ(ierr);
    ierr = VecSum(this->m_WorkScaField1, &meanj); CHKERRQ(ierr);
    meanj /= static_cast<ScalarType>(ng);

    // remember
    this->m_Opt->SetDetDGradMin(minj);
    this->m_Opt->SetDetDGradMax(maxj);
    this->m_Opt->SetDetDGradMean(meanj);

    if (this->m_Opt->GetVerbosity() > 1 || this->m_Opt->GetRegMonitor().detdgradenabled) {
        ss  << std::scientific << "det(grad(y)) : (min, mean, max)="
            << "(" << minj << ", " << meanj << ", " << maxj<<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }


    if (write2file) {
        ext = this->m_Opt->GetFileNames().extension;
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X11, "deformation-grad-x11"+ext); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X12, "deformation-grad-x12"+ext); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X13, "deformation-grad-x13"+ext); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X21, "deformation-grad-x21"+ext); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X22, "deformation-grad-x22"+ext); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X23, "deformation-grad-x23"+ext); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X31, "deformation-grad-x21"+ext); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X32, "deformation-grad-x22"+ext); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkTenField1->m_X33, "deformation-grad-x22"+ext); CHKERRQ(ierr);

        ierr = this->m_ReadWrite->Write(this->m_WorkScaField1,"det-deformation-grad-fulltensor"+ext); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute deformation gradient
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDefGradSL() {
    PetscErrorCode ierr = 0;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    IntType nt, nl;
    ScalarType ht, hthalf;
    ScalarType  *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_j11 = NULL, *p_j12 = NULL, *p_j13 = NULL,
                *p_j21 = NULL, *p_j22 = NULL, *p_j23 = NULL,
                *p_j31 = NULL, *p_j32 = NULL, *p_j33 = NULL,
                *p_gv11 = NULL, *p_gv12 = NULL, *p_gv13 = NULL,
                *p_gv21 = NULL, *p_gv22 = NULL, *p_gv23 = NULL,
                *p_gv31 = NULL, *p_gv32 = NULL, *p_gv33 = NULL,
                *p_j11X = NULL, *p_j12X = NULL, *p_j13X = NULL,
                *p_j21X = NULL, *p_j22X = NULL, *p_j23X = NULL,
                *p_j31X = NULL, *p_j32X = NULL, *p_j33X = NULL,
                *p_gv11X = NULL, *p_gv12X = NULL, *p_gv13X = NULL,
                *p_gv21X = NULL, *p_gv22X = NULL, *p_gv23X = NULL,
                *p_gv31X = NULL, *p_gv32X = NULL, *p_gv33X = NULL;
    double timer[7] = {0};
    PetscFunctionBegin;

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkTenField1 != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_SemiLagrangianMethod == NULL) {
        try{this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    }
    if (this->m_WorkTenField2 == NULL) {
       try {this->m_WorkTenField2 = new TenField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkTenField3 == NULL) {
       try {this->m_WorkTenField3 = new TenField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkTenField4 == NULL) {
       try {this->m_WorkTenField4 = new TenField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkTenField1->GetArrays(p_j11, p_j12, p_j13,
                                            p_j21, p_j22, p_j23,
                                            p_j31, p_j32, p_j33); CHKERRQ(ierr);

    ierr = this->m_WorkTenField2->GetArrays(p_j11X, p_j12X, p_j13X,
                                            p_j21X, p_j22X, p_j23X,
                                            p_j31X, p_j32X, p_j33X); CHKERRQ(ierr);

    ierr = this->m_WorkTenField3->GetArrays(p_gv11, p_gv12, p_gv13,
                                            p_gv21, p_gv22, p_gv23,
                                            p_gv31, p_gv32, p_gv33); CHKERRQ(ierr);

    ierr = this->m_WorkTenField4->GetArrays(p_gv11X, p_gv12X, p_gv13X,
                                            p_gv21X, p_gv22X, p_gv23X,
                                            p_gv31X, p_gv32X, p_gv33X); CHKERRQ(ierr);

    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_grad_t(p_gv11, p_gv12, p_gv13, p_v1, this->m_Opt->GetFFT().plan, &XYZ, timer);  ///< X1 gradient
    accfft_grad_t(p_gv21, p_gv22, p_gv23, p_v2, this->m_Opt->GetFFT().plan, &XYZ, timer);  ///< X2 gradient
    accfft_grad_t(p_gv31, p_gv32, p_gv33, p_v3, this->m_Opt->GetFFT().plan, &XYZ, timer);  ///< X3 gradient
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, 3*FFTGRAD);

    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    // interpolate gradient of velocity field
    // (TODO: write interpolation for tensor field)
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_gv11X, p_gv12X, p_gv13X,
                                                     p_gv11, p_gv12, p_gv13, "state"); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_gv21X, p_gv22X, p_gv23X,
                                                     p_gv21, p_gv22, p_gv23, "state"); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_gv31X, p_gv32X, p_gv33X,
                                                     p_gv31, p_gv32, p_gv33, "state"); CHKERRQ(ierr);

    // for all time points
    for (IntType j = 0; j < nt; ++j) {
        // evaluate j at X
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_j11X, p_j12X, p_j13X,
                                                         p_j11, p_j12, p_j13, "state"); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_j21X, p_j22X, p_j23X,
                                                         p_j21, p_j22, p_j23, "state"); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_j31X, p_j32X, p_j33X,
                                                         p_j31, p_j32, p_j33, "state"); CHKERRQ(ierr);

#pragma omp parallel
{
        ScalarType rhs11, rhs12, rhs13,
                   rhs21, rhs22, rhs23,
                   rhs31, rhs32, rhs33;
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {
            // evaluate right hand side at t^j
            rhs11 = p_gv11X[i]*p_j11X[i] + p_gv21X[i]*p_j21X[i] + p_gv31X[i]*p_j31X[i];
            rhs12 = p_gv11X[i]*p_j12X[i] + p_gv21X[i]*p_j22X[i] + p_gv31X[i]*p_j32X[i];
            rhs13 = p_gv11X[i]*p_j13X[i] + p_gv21X[i]*p_j23X[i] + p_gv31X[i]*p_j33X[i];

            rhs21 = p_gv12X[i]*p_j11X[i] + p_gv22X[i]*p_j21X[i] + p_gv32X[i]*p_j31X[i];
            rhs22 = p_gv12X[i]*p_j12X[i] + p_gv22X[i]*p_j22X[i] + p_gv32X[i]*p_j32X[i];
            rhs23 = p_gv12X[i]*p_j13X[i] + p_gv22X[i]*p_j23X[i] + p_gv32X[i]*p_j33X[i];

            rhs31 = p_gv13X[i]*p_j11X[i] + p_gv23X[i]*p_j21X[i] + p_gv33X[i]*p_j31X[i];
            rhs32 = p_gv13X[i]*p_j12X[i] + p_gv23X[i]*p_j22X[i] + p_gv33X[i]*p_j32X[i];
            rhs33 = p_gv13X[i]*p_j13X[i] + p_gv23X[i]*p_j23X[i] + p_gv33X[i]*p_j33X[i];

            // second stage of rk2 scheme
            p_j11[i] = p_j11X[i] + hthalf*( rhs11 + ( p_gv11X[i] + ht*(p_gv11[i]*rhs11 + p_gv21[i]*rhs21 + p_gv31[i]*rhs31) ) );
            p_j12[i] = p_j12X[i] + hthalf*( rhs12 + ( p_gv12X[i] + ht*(p_gv11[i]*rhs12 + p_gv21[i]*rhs22 + p_gv31[i]*rhs32) ) );
            p_j13[i] = p_j13X[i] + hthalf*( rhs13 + ( p_gv13X[i] + ht*(p_gv11[i]*rhs13 + p_gv21[i]*rhs23 + p_gv31[i]*rhs33) ) );

            p_j21[i] = p_j21X[i] + hthalf*( rhs21 + ( p_gv21X[i] + ht*(p_gv12[i]*rhs11 + p_gv22[i]*rhs21 + p_gv32[i]*rhs31) ) );
            p_j22[i] = p_j22X[i] + hthalf*( rhs22 + ( p_gv22X[i] + ht*(p_gv12[i]*rhs12 + p_gv22[i]*rhs22 + p_gv32[i]*rhs32) ) );
            p_j23[i] = p_j23X[i] + hthalf*( rhs23 + ( p_gv23X[i] + ht*(p_gv12[i]*rhs13 + p_gv22[i]*rhs23 + p_gv32[i]*rhs33) ) );

            p_j31[i] = p_j31X[i] + hthalf*( rhs31 + ( p_gv31X[i] + ht*(p_gv13[i]*rhs11 + p_gv23[i]*rhs21 + p_gv33[i]*rhs31) ) );
            p_j32[i] = p_j32X[i] + hthalf*( rhs32 + ( p_gv32X[i] + ht*(p_gv13[i]*rhs12 + p_gv23[i]*rhs22 + p_gv33[i]*rhs32) ) );
            p_j33[i] = p_j33X[i] + hthalf*( rhs33 + ( p_gv33X[i] + ht*(p_gv13[i]*rhs13 + p_gv23[i]*rhs23 + p_gv33[i]*rhs33) ) );
        }  // for all grid points
}  // pragma omp parallel
    }  // for all time points

    ierr = this->m_WorkTenField4->RestoreArrays(p_gv11X, p_gv12X, p_gv13X,
                                                p_gv21X, p_gv22X, p_gv23X,
                                                p_gv31X, p_gv32X, p_gv33X); CHKERRQ(ierr);

    ierr = this->m_WorkTenField3->RestoreArrays(p_gv11, p_gv12, p_gv13,
                                                p_gv21, p_gv22, p_gv23,
                                                p_gv31, p_gv32, p_gv33); CHKERRQ(ierr);

    ierr = this->m_WorkTenField2->RestoreArrays(p_j11X, p_j12X, p_j13X,
                                                p_j21X, p_j22X, p_j23X,
                                                p_j31X, p_j32X, p_j33X); CHKERRQ(ierr);

    ierr = this->m_WorkTenField1->RestoreArrays(p_j11, p_j12, p_j13,
                                                p_j21, p_j22, p_j23,
                                                p_j31, p_j32, p_j33); CHKERRQ(ierr);
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::CheckDefMapConsistency() {
    PetscErrorCode ierr = 0;
    VecField *y = NULL, *x = NULL;
    ScalarType value, normx;
    std::stringstream ss;
    bool flag;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate velocity field
    if (this->m_VelocityField == NULL) {
       try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    try {y = new VecField(this->m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {x = new VecField(this->m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // remember
    flag = this->m_ComputeInverseDefMap;

    // compute initial condition
    ierr = this->ComputeRegularGrid(x); CHKERRQ(ierr);

    ierr = y->Copy(x); CHKERRQ(ierr);
    this->m_ComputeInverseDefMap = false;
    ierr = this->ComputeDeformationMap(false, y); CHKERRQ(ierr);

    this->m_ComputeInverseDefMap = true;
    ierr = this->ComputeDeformationMap(false, y); CHKERRQ(ierr);

    // reset
    this->m_ComputeInverseDefMap = flag;

    ierr = y->AXPY(-1.0, x); CHKERRQ(ierr);

    ierr = VecNorm(y->m_X1, NORM_2, &value); CHKERRQ(ierr);
    ierr = VecNorm(x->m_X1, NORM_2, &normx); CHKERRQ(ierr);
    ss  << "error x1 " << std::scientific << value/normx << " (" << value << ")";
    ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();

    ierr = VecNorm(y->m_X2, NORM_2, &value); CHKERRQ(ierr);
    ierr = VecNorm(x->m_X2, NORM_2, &normx); CHKERRQ(ierr);
    ss  << "error x2 " << std::scientific << value/normx << " (" << value << ")";
    ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();

    ierr = VecNorm(y->m_X3, NORM_2, &value); CHKERRQ(ierr);
    ierr = VecNorm(x->m_X3, NORM_2, &normx); CHKERRQ(ierr);
    ss  << "error x3 " << std::scientific << value/normx << " (" << value << ")";
    ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();

    if (y != NULL) {delete y; y = NULL;}
    if (x != NULL) {delete x; x = NULL;}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMap(bool write2file, VecField* y) {
    PetscErrorCode ierr = 0;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate velocity field
    if (this->m_VelocityField == NULL) {
       try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }
    if(this->m_WorkVecField1 == NULL) {
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("computing deformation map"); CHKERRQ(ierr);
    }

    // check cfl condition / update time step
    if (this->m_Opt->GetPDESolverPara().adapttimestep) {
        ierr = this->ComputeCFLCondition(); CHKERRQ(ierr);
    }

    if (y == NULL) {
        // compute initial condition
        ierr = this->ComputeRegularGrid(this->m_WorkVecField1); CHKERRQ(ierr);
    } else {
        ierr = this->m_WorkVecField1->Copy(y); CHKERRQ(ierr);
    }

    // call the solver
    switch (this->m_Opt->GetPDESolverPara().type) {
        case RK2:
        {
            ierr = this->ComputeDeformationMapRK2(); CHKERRQ(ierr);
            break;
        }
        case SL:
        {
            switch (this->m_Opt->GetPDESolverPara().rungekuttaorder) {
                case 2:
                {
                    ierr = this->ComputeDeformationMapSLRK2(); CHKERRQ(ierr);
                    break;
                }
                case 4:
                {
                    ierr = this->ComputeDeformationMapSLRK4(); CHKERRQ(ierr);
                    break;
                }
                default:
                {
                    ierr = ThrowError("order not available"); CHKERRQ(ierr);
                }
            }
            break;
        }
        default:
        {
            ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
            break;
        }
    }

    if (write2file) {
        ext = this->m_Opt->GetFileNames().extension;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1, "deformation-map"+ext); CHKERRQ(ierr);
    }

    if (y != NULL) {
        ierr = y->Copy(this->m_WorkVecField1); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapRK2() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = ThrowError("not implemented"); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapRK2A() {
    PetscErrorCode ierr = 0;
/*
    IntType nl,ng,nt;
    ScalarType *p_u1 = NULL, *p_u2 = NULL, *p_u3 = NULL,
                *p_rhs01 = NULL, *p_rhs02 = NULL, *p_rhs03 = NULL,
                *p_gu11 = NULL, *p_gu12 = NULL, *p_gu13 = NULL,
                *p_gu21 = NULL, *p_gu22 = NULL, *p_gu23 = NULL,
                *p_gu31 = NULL, *p_gu32 = NULL, *p_gu33 = NULL,
                *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_divv=NULL;
    ScalarType ht=0.0,hthalf=0.0;
    double timer[7] = {0};
    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if(this->m_WorkVecField1 == NULL) {
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if(this->m_WorkVecField2 == NULL) {
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if(this->m_WorkVecField3 == NULL) {
        try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if(this->m_WorkVecField4 == NULL) {
        try{this->m_WorkVecField4 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1,nl,ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2,nl,ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL) {
        ierr = VecCreate(this->m_WorkScaField3,nl,ng); CHKERRQ(ierr);
    }

    ierr = VecGetArray(this->m_WorkScaField1,&p_divv); CHKERRQ(ierr);

    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField2->GetArrays(p_gu11, p_gu12, p_gu13); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->GetArrays(p_gu21, p_gu22, p_gu23); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->GetArrays(p_gu31, p_gu32, p_gu33); CHKERRQ(ierr);

    ierr = this->m_WorkVecField5->GetArrays(p_rhs01, p_rhs02, p_rhs03); CHKERRQ(ierr);
    ierr = this->m_WorkVecField6->GetArrays(p_rhs01, p_rhs02, p_rhs03); CHKERRQ(ierr);

    // compute div(v)
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(p_divv, p_v1, p_v2, p_v3,this->m_Opt->GetFFT().plan,timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, DIV);

    // copy initial condition to buffer
    ierr = VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2,&p_mj); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField1,&p_mbar); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);

    // copy memory (m_0 to m_j)
    try{ std::copy(p_m, p_m+nl, p_mj); }
    catch(std::exception&) {
        ierr = ThrowError("copy failed"); CHKERRQ(ierr);
    }

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {

        // compute gradient of m_j
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gu11, p_gu12, p_gu13, p_u1,this->m_Opt->GetFFT().plan,&XYZ,timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gu21, p_gu22, p_gu23, p_u2,this->m_Opt->GetFFT().plan,&XYZ,timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gu31, p_gu32, p_gu33, p_u3,this->m_Opt->GetFFT().plan,&XYZ,timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
        for (IntType i=0; i < nl; ++i) {

             p_rhs01[i] = -p_gu11[i]*p_v1[i]-p_gu12[i]*p_v2[i]-p_gu13[i]*p_v3[i];
             p_rhs02[i] = -p_gu21[i]*p_v1[i]-p_gu22[i]*p_v2[i]-p_gu23[i]*p_v3[i];
             p_rhs03[i] = -p_gu31[i]*p_v1[i]-p_gu32[i]*p_v2[i]-p_gu33[i]*p_v3[i];

             // compute intermediate result
             p_mbar[i] = p_mj[i] + ht*p_rhs0[i];

        }
} // pragma omp parallel

        // compute div(v)
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_divergence_t(p_divv, p_v1, p_v2, p_v3,this->m_Opt->GetFFT().plan,timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTDIV);


        // compute gradient of \bar{m}
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gmx1, p_gmx2, p_gmx3, p_mbar,this->m_Opt->GetFFT().plan,&XYZ,timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
        for (IntType i=0; i < nl; ++i) {

            ScalarType rhs1 = -p_gmx1[i]*p_vx1[i]
                              -p_gmx2[i]*p_vx2[i]
                              -p_gmx3[i]*p_vx3[i];

            // we have overwritten m_j with intermediate result
            // m_j = m_{j-1} + 0.5*ht*(RHS0 + RHS1)
            p_mj[i] = p_mj[i] + hthalf*(p_rhs0[i] + rhs1);
        }
} // parallel

        // copy to buffer
        try{ std::copy(p_mj, p_mj+nl, p_m+(j+1)*nl); }
        catch(std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }

    } // for all time points

    // copy initial condition to buffer
    ierr = VecGetArray(this->m_WorkScaField1,&p_divv); CHKERRQ(ierr);

    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField2->RestoreArrays(p_gu11, p_gu12, p_gu13); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_gu21, p_gu22, p_gu23); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->RestoreArrays(p_gu31, p_gu32, p_gu33); CHKERRQ(ierr);

    ierr = this->m_WorkVecField5->RestoreArrays(p_rhs01, p_rhs02, p_rhs03); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);
*/

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map if we consider a semi-lagrangian
 * time integrator; the scheme is full lagrangian; we use an
 * rk2 scheme to compute the characteristic;
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapSLRK2() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    std::string ext;
    IntType nl, nt;
    ScalarType ht, hthalf;
    ScalarType *p_y1 = NULL, *p_y2 = NULL, *p_y3 = NULL,
                *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_vy1 = NULL, *p_vy2 = NULL, *p_vy3 = NULL,
                *p_vytilde1 = NULL, *p_vytilde2 = NULL, *p_vytilde3 = NULL,
                *p_ytilde1 = NULL, *p_ytilde2 = NULL, *p_ytilde3 = NULL;

    PetscFunctionBegin;

    // TODO: fix this
    ierr = ThrowError("not implemented"); CHKERRQ(ierr);

    this->m_Opt->Enter(__func__);

    ext = this->m_Opt->GetFileNames().extension;

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate vector fields
    if (this->m_WorkVecField1 == NULL) {
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3 == NULL) {
        try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField4 == NULL) {
        try{this->m_WorkVecField4 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField5 == NULL) {
        try{this->m_WorkVecField5 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate semi-lagrangian solver
    if(this->m_SemiLagrangianMethod == NULL) {
        try{this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField5); CHKERRQ(ierr);


    // store time series
    if (this->m_Opt->GetReadWriteFlags().timeseries ) {
        ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
        ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << 0 << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1,ss.str()); CHKERRQ(ierr);
    }


    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    // if we request the inverse deformation map
    if (this->m_ComputeInverseDefMap) {ht *= -1.0; hthalf *= -1.0;}

    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->GetArrays(p_vy1, p_vy2, p_vy3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->GetArrays(p_vytilde1, p_vytilde2, p_vytilde3); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {
        // evaluate v(y)
//        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3,
//                                                         p_v1, p_v2, p_v3,
//                                                         p_y1, p_y2, p_y3 ); CHKERRQ(ierr);

        // compute intermediate variable (fist stage of RK2)
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {
            p_ytilde1[i] = p_y1[i] - ht*p_vy1[i];
            p_ytilde2[i] = p_y2[i] - ht*p_vy2[i];
            p_ytilde3[i] = p_y3[i] - ht*p_vy3[i];
        }
}// end of pragma omp parallel

        // evaluate v(ytilde)
//        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vytilde1, p_vytilde2, p_vytilde3,
//                                                         p_v1, p_v2, p_v3,
//                                                         p_ytilde1, p_ytilde2, p_ytilde3 ); CHKERRQ(ierr);

        // update deformation map (second stage of RK2)
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {
            p_y1[i] = p_y1[i] - hthalf*(p_vytilde1[i] + p_vy1[i]);
            p_y2[i] = p_y2[i] - hthalf*(p_vytilde2[i] + p_vy2[i]);
            p_y3[i] = p_y3[i] - hthalf*(p_vytilde3[i] + p_vy3[i]);
        }
}// end of pragma omp parallel

        // store time series
        if (this->m_Opt->GetReadWriteFlags().timeseries ) {
            ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);
            ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
            ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << j+1 << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkVecField1, ss.str()); CHKERRQ(ierr);
            ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
        }
    } // for all time points

    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_vy1, p_vy2, p_vy3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->RestoreArrays(p_vytilde1, p_vytilde2, p_vytilde3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map if we consider a semi-lagrangian
 * time integrator; the scheme is full lagrangian; we use an
 * rk4 scheme to compute the characteristic;
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapSLRK4() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    std::string ext;
    IntType nl, nt;
    ScalarType ht, hthalf, htby6;
    ScalarType *p_y1 = NULL, *p_y2 = NULL, *p_y3 = NULL,
                *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_vy1 = NULL, *p_vy2 = NULL, *p_vy3 = NULL,
                *p_dy1 = NULL, *p_dy2 = NULL, *p_dy3 = NULL,
                *p_ytilde1 = NULL, *p_ytilde2 = NULL, *p_ytilde3 = NULL;

    PetscFunctionBegin;

    // TODO: fix this
    ierr = ThrowError("not implemented"); CHKERRQ(ierr);

    this->m_Opt->Enter(__func__);

    ext = this->m_Opt->GetFileNames().extension;

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate vector fields
    if (this->m_WorkVecField1 == NULL) {
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3 == NULL) {
        try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField4 == NULL) {
        try{this->m_WorkVecField4 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate semi-lagrangian solver
    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // store time series
    if (this->m_Opt->GetReadWriteFlags().timeseries ) {
        ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
        ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << 0 << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1,ss.str()); CHKERRQ(ierr);
    }

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;
    htby6  = ht/6.0;

    // if we request the inverse deformation map
    if (this->m_ComputeInverseDefMap) {ht *= -1.0; hthalf *= -1.0; htby6 *= -1.0;}

    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->GetArrays(p_vy1, p_vy2, p_vy3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->GetArrays(p_dy1, p_dy2, p_dy3); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {
        // evaluate right hand side v(y) (i.e., F0)
//        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3, p_v1, p_v2, p_v3,
//                                                         p_y1, p_y2, p_y3); CHKERRQ(ierr);

        // compute intermediate variable (fist stage of RK4); F0
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {
            p_ytilde1[i] = p_y1[i] - hthalf*p_vy1[i];
            p_ytilde2[i] = p_y2[i] - hthalf*p_vy2[i];
            p_ytilde3[i] = p_y3[i] - hthalf*p_vy3[i];

            // F0
            p_dy1[i] = p_vy1[i];
            p_dy2[i] = p_vy2[i];
            p_dy3[i] = p_vy3[i];
        }
}  // end of pragma omp parallel

        // evaluate right hand side v(ytilde) (i.e., F1)
//        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3, p_v1, p_v2, p_v3,
//                                                         p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);

        // compute intermediate variable (sedond stage of RK4)
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {
            p_ytilde1[i] = p_y1[i] - hthalf*p_vy1[i];
            p_ytilde2[i] = p_y2[i] - hthalf*p_vy2[i];
            p_ytilde3[i] = p_y3[i] - hthalf*p_vy3[i];

            // F0 + 2.0*F1
            p_dy1[i] += 2.0*p_vy1[i];
            p_dy2[i] += 2.0*p_vy2[i];
            p_dy3[i] += 2.0*p_vy3[i];
        }
}  // end of pragma omp parallel

        // evaluate right hand side v(ytilde) (i.e., F2)
//        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3, p_v1, p_v2, p_v3,
//                                                         p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);

        // compute intermediate variable (sedond stage of RK4)
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {
            p_ytilde1[i] = p_y1[i] - ht*p_vy1[i];
            p_ytilde2[i] = p_y2[i] - ht*p_vy2[i];
            p_ytilde3[i] = p_y3[i] - ht*p_vy3[i];

            // F0 + 2.0*F1 + 2.0*F2
            p_dy1[i] += 2.0*p_vy1[i];
            p_dy2[i] += 2.0*p_vy2[i];
            p_dy3[i] += 2.0*p_vy3[i];
        }
}  // end of pragma omp parallel

        // evaluate right hand side v(ytilde) (i.e., F3)
//        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3, p_v1, p_v2, p_v3,
//                                                         p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);

        // compute intermediate variable (sedond stage of RK4)
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {
            // y_{k+1} = y_k - ht/6*(F0 + 2.0*F1 + 2.0*F2 + F3)
            p_y1[i] = p_y1[i] - htby6*(p_dy1[i]+p_vy1[i]);
            p_y2[i] = p_y2[i] - htby6*(p_dy2[i]+p_vy2[i]);
            p_y3[i] = p_y3[i] - htby6*(p_dy3[i]+p_vy3[i]);
        }
}  // end of pragma omp parallel

        // store time series
        if (this->m_Opt->GetReadWriteFlags().timeseries) {
            ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);
            ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
            ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << j+1 << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkVecField1, ss.str()); CHKERRQ(ierr);
            ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
        }
    }  // for all time points

    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_vy1, p_vy2, p_vy3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->RestoreArrays(p_dy1, p_dy2, p_dy3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute displacement field
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDisplacementField(bool write2file) {
    PetscErrorCode ierr = 0;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate velocity field
    if (this->m_VelocityField == NULL) {
       try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("computing displacement field"); CHKERRQ(ierr);
    }

    // call the solver
    switch (this->m_Opt->GetPDESolverPara().type) {
        case RK2:
        {
            // compute displacement field using rk2 time integrator
            ierr = this->ComputeDisplacementFieldRK2(); CHKERRQ(ierr);
            break;
        }
        case SL:
        {
            // compute displacement field using sl time integrator
            ierr = this->ComputeDisplacementFieldSL(); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
            break;
        }
    }

    if (write2file) {
        ext = this->m_Opt->GetFileNames().extension;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1, "displacement-field"+ext); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute displacement field
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDisplacementFieldRK2() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = ThrowError("not implemented"); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute displacement field
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDisplacementFieldSL() {
    PetscErrorCode ierr = 0;
    IntType nl, nt;
    ScalarType ht, hthalf;
    std::stringstream ss;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_vX1 = NULL, *p_vX2 = NULL, *p_vX3 = NULL,
                *p_u1 = NULL, *p_u2 = NULL, *p_u3 = NULL,
                *p_uX1 = NULL, *p_uX2 = NULL, *p_uX3 = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_WorkVecField1 == NULL) {
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3 == NULL) {
        try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate semi-lagrangian solver
    if(this->m_SemiLagrangianMethod == NULL) {
        try{this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1, "state"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    // set initial condition
    ierr = this->m_WorkVecField1->SetValue(0.0); CHKERRQ(ierr);

    // evaluate v(y)
    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField2, this->m_VelocityField, "state"); CHKERRQ(ierr);

    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {
        // interpolate u^j at X
        ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField3, this->m_WorkVecField1, "state"); CHKERRQ(ierr);

        ierr = this->m_WorkVecField1->GetArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);
        ierr = this->m_WorkVecField3->GetArrays(p_uX1, p_uX2, p_uX3); CHKERRQ(ierr);
        // update deformation field (RK2)
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {
            p_u1[i] = p_uX1[i] + hthalf*(p_vX1[i] + p_v1[i]);
            p_u2[i] = p_uX2[i] + hthalf*(p_vX2[i] + p_v2[i]);
            p_u3[i] = p_uX3[i] + hthalf*(p_vX3[i] + p_v3[i]);
        }
}  // end of pragma omp parallel
        ierr = this->m_WorkVecField3->RestoreArrays(p_uX1, p_uX2, p_uX3); CHKERRQ(ierr);
        ierr = this->m_WorkVecField1->RestoreArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);
    }  // for all time points

    ierr = this->m_WorkVecField2->RestoreArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map from a displacement field
 *******************************************************************/
PetscErrorCode OptimalControlRegistrationBase::ComputeDefMapFromDisplacement() {
    PetscErrorCode ierr = 0;
    ScalarType hx[3];
    ScalarType *p_u1 = NULL,*p_u2 = NULL, *p_u3 = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate velocity field
    if (this->m_VelocityField == NULL) {
       try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("computing displacement field"); CHKERRQ(ierr);
    }

    // call the solver
    switch (this->m_Opt->GetPDESolverPara().type) {
        case RK2:
        {
            // compute displacement field using rk2 time integrator
            ierr = this->ComputeDisplacementFieldRK2(); CHKERRQ(ierr);
            break;
        }
        case SL:
        {
            // compute displacement field using sl time integrator
            ierr = this->ComputeDisplacementFieldSL(); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
            break;
        }
    }

    // get spatial step size
    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->GetDomainPara().hx[i];
    }

    ierr = this->m_WorkVecField1->GetArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);
#pragma omp parallel
{
    IntType i, i1, i2, i3;
    ScalarType x1, x2, x3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->GetDomainPara().isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < this->m_Opt->GetDomainPara().isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < this->m_Opt->GetDomainPara().isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + this->m_Opt->GetDomainPara().istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + this->m_Opt->GetDomainPara().istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + this->m_Opt->GetDomainPara().istart[2]);

                // compute linear / flat index
                i = GetLinearIndex(i1, i2, i3, this->m_Opt->GetDomainPara().isize);

                // assign values
                p_u1[i] = x1 + p_u1[i];
                p_u2[i] = x2 + p_u2[i];
                p_u3[i] = x3 + p_u3[i];

            } // i1
        } // i2
    } // i3
}// pragma omp for
    ierr = this->m_WorkVecField1->RestoreArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _OPTIMALCONTROLREGISTRATIONBASE_CPP_
