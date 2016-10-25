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
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistrationBase"
OptimalControlRegistrationBase::OptimalControlRegistrationBase() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~OptimalControlRegistrationBase"
OptimalControlRegistrationBase::~OptimalControlRegistrationBase(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistrationBase"
OptimalControlRegistrationBase::OptimalControlRegistrationBase(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode OptimalControlRegistrationBase::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // pointer for container of velocity field
    this->m_VelocityField = NULL;

    // pointer for container of incremental velocity field
    this->m_IncVelocityField = NULL;

    // pointers to images
    this->m_TemplateImage = NULL;
    this->m_ReferenceImage = NULL;

    this->m_ReadWrite = NULL;             ///< read / write object
    this->m_SemiLagrangianMethod = NULL;  ///< semi lagranigan

    this->m_WorkScaField1 = NULL;
    this->m_WorkScaField2 = NULL;
    this->m_WorkScaField3 = NULL;
    this->m_WorkScaField4 = NULL;
    this->m_WorkScaField5 = NULL;

    this->m_WorkVecField1 = NULL;
    this->m_WorkVecField2 = NULL;
    this->m_WorkVecField3 = NULL;
    this->m_WorkVecField4 = NULL;
    this->m_WorkVecField5 = NULL;

    this->m_WorkTenField1 = NULL;
    this->m_WorkTenField2 = NULL;
    this->m_WorkTenField3 = NULL;
    this->m_WorkTenField4 = NULL;

    this->m_VelocityIsZero = false;  ///< flag for velocity field

    this->m_Regularization = NULL;   ///< pointer for regularization class

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode OptimalControlRegistrationBase::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_VelocityField != NULL) {
        delete this->m_VelocityField;
        this->m_VelocityField = NULL;
    }

    if (this->m_IncVelocityField != NULL) {
        delete this->m_IncVelocityField;
        this->m_IncVelocityField = NULL;
    }

    if (this->m_Regularization != NULL) {
        delete this->m_Regularization;
        this->m_Regularization = NULL;
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

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReadWrite"
PetscErrorCode OptimalControlRegistrationBase::SetReadWrite(ReadWriteReg* rw) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(rw != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite = rw;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReferenceImage"
PetscErrorCode OptimalControlRegistrationBase::SetReferenceImage(Vec mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Rescale(mR, 0.0, 1.0); CHKERRQ(ierr);
    this->m_ReferenceImage = mR;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set template image (i.e., the image to be registered)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetTemplateImage"
PetscErrorCode OptimalControlRegistrationBase::SetTemplateImage(Vec mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Rescale(mT, 0.0, 1.0); CHKERRQ(ierr);
    this->m_TemplateImage = mT;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get reference image (i.e., the fixed image)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetReferenceImage"
PetscErrorCode OptimalControlRegistrationBase::GetReferenceImage(Vec& mR) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
    mR = this->m_ReferenceImage;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get template image (i.e., the image to be registered)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetTemplateImage"
PetscErrorCode OptimalControlRegistrationBase::GetTemplateImage(Vec& mT) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    mT = this->m_TemplateImage;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetControlVariable"
PetscErrorCode OptimalControlRegistrationBase::SetControlVariable(VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate velocity field
    if(this->m_VelocityField == NULL) {
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // copy buffer/values
    ierr = this->m_VelocityField->Copy(v); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetControlVariable"
PetscErrorCode OptimalControlRegistrationBase::GetControlVariable(VecField*& v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    // copy buffer
    ierr = v->Copy(this->m_VelocityField); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetControlVariable"
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
        nl = this->m_Opt->GetDomainPara().nlocal;
        ng = this->m_Opt->GetDomainPara().nglobal;
        ierr = VecCreate(v, nl, ng); CHKERRQ(ierr);
        ierr = VecCreate(g, nl, ng); CHKERRQ(ierr);

        ierr = this->m_VelocityField->GetComponents(v); CHKERRQ(ierr);
        ierr = this->EvaluatePrecondGradient(g, v); CHKERRQ(ierr);

        ierr = VecAXPY(v, -1.0, g); CHKERRQ(ierr);

        ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}






/********************************************************************
 * @brief check if velocity field is zero
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "IsVelocityZero"
PetscErrorCode OptimalControlRegistrationBase::IsVelocityZero() {
    PetscErrorCode ierr = 0;
    ScalarType normv1, normv2, normv3;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    this->m_VelocityIsZero = false;
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = VecNorm(this->m_VelocityField->m_X1, NORM_INFINITY, &normv1); CHKERRQ(ierr);
    ierr = VecNorm(this->m_VelocityField->m_X2, NORM_INFINITY, &normv2); CHKERRQ(ierr);
    ierr = VecNorm(this->m_VelocityField->m_X3, NORM_INFINITY, &normv3); CHKERRQ(ierr);

    this->m_VelocityIsZero = (normv1 == 0.0) && (normv2 == 0.0) && (normv3 == 0.0);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief allocate regularization model
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "AllocateRegularization"
PetscErrorCode OptimalControlRegistrationBase::AllocateRegularization() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // delete regularization if already allocated
    // (should never happen)
    if (this->m_Regularization != NULL) {
        delete this->m_Regularization;
        this->m_Regularization = NULL;
    }

    // switch between regularization norms
    switch (this->m_Opt->GetRegNorm().type) {
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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies inverse regularization operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvRegOp"
PetscErrorCode OptimalControlRegistrationBase::ApplyInvRegOp(Vec Ainvx, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief estimate eigenvalues of hessian
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EstimateExtremalHessEigVals"
PetscErrorCode OptimalControlRegistrationBase
::EstimateExtremalHessEigVals(ScalarType &emin,
                              ScalarType &emax) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    ierr = this->m_Regularization->GetExtremeEigValsInvOp(emin, emax); CHKERRQ(ierr);
    emin += 1.0;
    emax += 1.0; // this is crap

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief pre-processing before the krylov solve
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PreKrylovSolve"
PetscErrorCode OptimalControlRegistrationBase
::PreKrylovSolve(Vec g, Vec x) {
    PetscErrorCode ierr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief post-processing after the krylov solve
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PostKrylovSolve"
PetscErrorCode OptimalControlRegistrationBase::PostKrylovSolve(Vec g, Vec x) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply projection operator (this is for the analytically/
 * spectrally preconditioned hessian; it applies the projection
 * operator to the search direction and the gradient)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvRegOpSqrt"
PetscErrorCode OptimalControlRegistrationBase::ApplyInvRegOpSqrt(Vec x) {
    PetscErrorCode ierr = 0;
    bool usesqrt = true;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute a synthetic test problem
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupSyntheticProb"
PetscErrorCode OptimalControlRegistrationBase::SetupSyntheticProb(Vec &mR, Vec &mT) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nx[3];
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
               *p_mt = NULL, hx[3], xc1, xc2, xc3, x, sigma;
    int vcase = 3;
    int icase = 0;
    ScalarType v0 = 0.2;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("setting up synthetic test problem"); CHKERRQ(ierr);
    }
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

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
    }

    if (this->m_Opt->GetRegModel() == STOKES) {vcase = 4;}

    // allocate reference image
    if (mR == NULL) {ierr = VecCreate(mR, nl, ng); CHKERRQ(ierr);}
    ierr = VecSet(mR, 0.0); CHKERRQ(ierr);

    // allocate template image
    if (mT == NULL) {ierr = VecCreate(mT, nl, ng); CHKERRQ(ierr);}
    ierr = VecSet(mT, 0.0); CHKERRQ(ierr);


    // compute coordinates (nodal grid)
    xc1 = hx[0]*static_cast<ScalarType>(nx[0])/2.0;
    xc2 = hx[1]*static_cast<ScalarType>(nx[1])/2.0;
    xc3 = hx[2]*static_cast<ScalarType>(nx[2])/2.0;

    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = VecGetArray(mT, &p_mt); CHKERRQ(ierr);
#pragma omp parallel
{
    ScalarType x1, x2, x3;
    IntType i1, i2, i3, i;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->GetDomainPara().isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < this->m_Opt->GetDomainPara().isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < this->m_Opt->GetDomainPara().isize[2]; ++i3) {  // x3

                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + this->m_Opt->GetDomainPara().istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + this->m_Opt->GetDomainPara().istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + this->m_Opt->GetDomainPara().istart[2]);

                // compute linear / flat index
                i = GetLinearIndex(i1,i2,i3,this->m_Opt->GetDomainPara().isize);

                if (icase == 0) {
                    p_mt[i] = (sin(x1)*sin(x1) + sin(x2)*sin(x2) + sin(x3)*sin(x3))/3.0;
                } else if (icase == 1) {
                    sigma = 0.05*2.0*PETSC_PI;
                    x1 -= xc1; x2 -= xc2; x3 -= xc3;
                    x   = sqrt(x1*x1 + x2*x2 + x3*x3)/sigma;
                    p_mt[i] = exp(-x*x);
                }

                if (vcase == 0) {
                    // compute the velocity field
                    p_vx1[i] = v0*sin(x3)*cos(x2)*sin(x2);
                    p_vx2[i] = v0*sin(x1)*cos(x3)*sin(x3);
                    p_vx3[i] = v0*sin(x2)*cos(x1)*sin(x1);
                } else if (vcase == 1) {
                    // compute the velocity field
                    p_vx1[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                    p_vx2[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                    p_vx3[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                } else if (vcase == 3) {
                    p_vx1[i] = cos(x1)*sin(x2);
                    p_vx2[i] = cos(x2)*sin(x1);
                    p_vx3[i] = cos(x1)*sin(x3);
                } else if (vcase == 4) {
                    p_vx1[i] = cos(x2)*cos(x3);
                    p_vx2[i] = sin(x3)*sin(x1);
                    p_vx3[i] = cos(x1)*cos(x2);
                } else if (vcase == 5) {
                    p_vx1[i] = v0;
                    p_vx2[i] = v0;
                    p_vx3[i] = v0;
                }
            } // i1
        } // i2
    } // i3
} // pragma omp parallel

    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = VecRestoreArray(mT, &p_mt); CHKERRQ(ierr);

    ierr = Rescale(mT, 0.0, 1.0); CHKERRQ(ierr);

    // solve the forward problem using the computed
    // template image and the computed velocity field as input
    ierr = this->SolveForwardProblem(mR, mT); CHKERRQ(ierr);
    ierr = Rescale(mR, 0.0, 1.0); CHKERRQ(ierr);

    // reset velocity field
    ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief copies some input data field to all time points
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CopyToAllTimePoints"
PetscErrorCode OptimalControlRegistrationBase::CopyToAllTimePoints(Vec u, Vec uj) {
    PetscErrorCode ierr = 0;
    ScalarType *p_u = NULL, *p_uj = NULL;
    IntType nl, nt;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;

    // get pointers
    ierr = VecGetArray(u, &p_u); CHKERRQ(ierr);      ///< vec for entire time horizon
    ierr = VecGetArray(uj, &p_uj); CHKERRQ(ierr);    ///< vec at single point in time

    // for all time points
    for (IntType j = 0; j <= nt; ++j) {
        // copy data to all time points
        try{ std::copy(p_uj, p_uj+nl, p_u+j*nl); }
        catch(std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
    }

    // restore pointers
    ierr = VecRestoreArray(u, &p_u); CHKERRQ(ierr);
    ierr = VecRestoreArray(uj, &p_uj); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief copies some input data field to all time points
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeCFLCondition"
PetscErrorCode OptimalControlRegistrationBase::ComputeCFLCondition() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    ScalarType hx[3], c, vmax, vmaxscaled;
    IntType nl, ng, ntcfl;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

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
    ntcfl = this->m_Opt->GetDomainPara().nt;
    if ( vmaxscaled != 0.0 ) {
        c  = 1.0;//this->m_opt->cflnum;
        // compute min number of time steps
        ntcfl = static_cast<IntType>(ceil(vmaxscaled/c));
    }

    ss << "||v||_infty = " << std::scientific << std::fixed << vmax << " nt_CFL = " << ntcfl;
    ierr = DbgMsg(ss.str()); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CheckBounds"
PetscErrorCode OptimalControlRegistrationBase::CheckBounds(Vec v, bool& boundreached) {
    PetscErrorCode ierr = 0;
    ScalarType jmin, jmax, jbound;
    bool minboundreached, maxboundreached;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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

    jmin   = this->m_Opt->GetRegMonitor().jacmin;
    jmax   = this->m_Opt->GetRegMonitor().jacmax;
    jbound = this->m_Opt->GetRegMonitor().jacbound;

    // check if jmin < bound and 1/jmax < bound
    minboundreached = jmin <= jbound;
    maxboundreached = 1.0/jmax <= jbound;

    boundreached = (minboundreached || maxboundreached) ? true : false;
    if (boundreached) { ierr = WrngMsg("jacobian bound reached"); CHKERRQ(ierr); }

    // display what's going on
    if (this->m_Opt->GetVerbosity() > 1) {
        if(minboundreached) {
            ss << std::scientific
            << "min(det(grad(y))) = "<< jmin << " <= " << jbound;
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();
        }
        if(maxboundreached) {
            ss << std::scientific
            << "max(det(grad(y))) = "<< jmax << " >= " << 1.0/jbound
            << " ( = 1/bound )";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();
        }
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDetDefGrad"
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGrad(bool write2file) {
    PetscErrorCode ierr = 0;
    ScalarType minddg, maxddg, meanddg;
    IntType nl, ng;
    std::string ext;
    std::stringstream ss, ssnum;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("computing determinant of deformation gradient"); CHKERRQ(ierr);
    }

    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    if (this->m_VelocityField == NULL) {
       try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
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
            switch (this->m_Opt->GetPDESolver().type) {
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
    meanddg /= static_cast<ScalarType>(this->m_Opt->GetDomainPara().nglobal);

    // remember
    this->m_Opt->SetJacMin(minddg);
    this->m_Opt->SetJacMax(maxddg);
    this->m_Opt->SetJacMean(meanddg);

    if (this->m_Opt->GetVerbosity() > 1 || this->m_Opt->GetRegMonitor().JAC) {
        ss  << std::scientific << "det(grad(y)) : (min, mean, max)="
            << "(" << minddg << ", " << meanddg << ", " << maxddg<<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    if (write2file) {
        ext = this->m_Opt->GetReadWriteFlags().extension;
        ierr = this->m_ReadWrite->Write(this->m_WorkScaField1, "det-deformation-grad"+ext); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDefGrad"
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

    this->m_Opt->Enter(__FUNCT__);

    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

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
            switch (this->m_Opt->GetPDESolver().type) {
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
    this->m_Opt->SetJacMin(minj);
    this->m_Opt->SetJacMax(maxj);
    this->m_Opt->SetJacMean(meanj);

    if (this->m_Opt->GetVerbosity() > 1 || this->m_Opt->GetRegMonitor().JAC) {
        ss  << std::scientific << "det(grad(y)) : (min, mean, max)="
            << "(" << minj << ", " << meanj << ", " << maxj<<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }


    if (write2file) {
        ext = this->m_Opt->GetReadWriteFlags().extension;
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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief compute deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDefGradSL"
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
    double timer[5] = {0, 0, 0, 0, 0};
    PetscFunctionBegin;

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkTenField1 != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_SemiLagrangianMethod == NULL) {
        try{this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
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

    accfft_grad(p_gv11, p_gv12, p_gv13, p_v1, this->m_Opt->GetFFT().plan, &XYZ, timer);  ///< X1 gradient
    accfft_grad(p_gv21, p_gv22, p_gv23, p_v2, this->m_Opt->GetFFT().plan, &XYZ, timer);  ///< X2 gradient
    accfft_grad(p_gv31, p_gv32, p_gv33, p_v3, this->m_Opt->GetFFT().plan, &XYZ, timer);  ///< X3 gradient

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
        for (IntType i=0; i < nl; ++i) {
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
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDetDefGradRK2"
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradRK2() {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nt;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL, *p_jbar = NULL,
                *p_gx1 = NULL, *p_gx2 = NULL, *p_gx3 = NULL, *p_divv = NULL,
                *p_jac = NULL,  *p_rhs0 = NULL;
    ScalarType ht, hthalf, velsign;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timings[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_Opt->GetRegFlags().invdefgrad) {
        velsign = -1.0;
    } else {
        velsign = 1.0;
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

    // set initial condition
    ierr = VecSet(this->m_WorkScaField1, 1.0); CHKERRQ(ierr);

    // get pointers
    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gx1, p_gx2, p_gx3); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_WorkScaField1, &p_jac); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2, &p_divv); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField3, &p_jbar); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField4, &p_rhs0); CHKERRQ(ierr);

    // compute div(v)
    accfft_divergence(p_divv, p_vx1, p_vx2, p_vx3,this->m_Opt->GetFFT().plan,timings);

    // for all time points
    for (IntType j = 0; j <= nt; ++j) {
        accfft_grad(p_gx1, p_gx2, p_gx3, p_jac,this->m_Opt->GetFFT().plan,&XYZ,timings);

#pragma omp parallel
{
#pragma omp  for
        for (IntType i=0; i < nl; ++i) {  // for all grid points
            // \bar{j} = -(\vect{v} \cdot \igrad) j + j (\idiv \vect{v})
            p_rhs0[i] = -( velsign*p_vx1[i]*p_gx1[i]
                         + velsign*p_vx2[i]*p_gx2[i]
                         + velsign*p_vx3[i]*p_gx3[i] )
                         + p_jac[i]*velsign*p_divv[i];

            p_jbar[i] = p_jac[i] + ht*p_rhs0[i];
        }
} // pragma omp

        accfft_grad(p_gx1, p_gx2, p_gx3, p_jbar,this->m_Opt->GetFFT().plan,&XYZ,timings);

#pragma omp parallel
{
#pragma omp  for
        for (IntType i=0; i < nl; ++i) {  // for all grid points
            // \bar{j} = -(\vect{v} \cdot \igrad) j + j (\idiv \vect{v})
            ScalarType rhs1 = -( velsign*p_vx1[i]*p_gx1[i]
                               + velsign*p_vx2[i]*p_gx2[i]
                               + velsign*p_vx3[i]*p_gx3[i] )
                               + p_jbar[i]*velsign*p_divv[i];

            p_jac[i] = p_jac[i] + hthalf*(p_rhs0[i] + rhs1);
        }
}  // pragma omp
    }

    ierr = VecRestoreArray(this->m_WorkScaField1, &p_jac); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_divv); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField3, &p_jbar); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField4, &p_rhs0); CHKERRQ(ierr);

    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_gx1, p_gx2, p_gx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDetDefGradRK2A"
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradRK2A() {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nt;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL, *p_phibar = NULL,
                *p_gphi1 = NULL, *p_gphi2 = NULL, *p_gphi3 = NULL, *p_divv = NULL,
                *p_phiv1 = NULL, *p_phiv2 = NULL, *p_phiv3 = NULL,
                *p_phi = NULL,  *p_rhs0 = NULL,  *p_divvphi=NULL;
    ScalarType ht, hthalf;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timings[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;
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
    accfft_divergence(p_divv, p_vx1, p_vx2, p_vx3,this->m_Opt->GetFFT().plan,timings);


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
        accfft_grad(p_gphi1, p_gphi2, p_gphi3, p_phi,this->m_Opt->GetFFT().plan,&XYZ,timings);

        // compute div(\vect{v}\phi_j)
        accfft_divergence(p_divvphi, p_phiv1, p_phiv2, p_phiv3,this->m_Opt->GetFFT().plan,timings);

#pragma omp parallel
{
#pragma omp  for
        for (IntType i=0; i < nl; ++i) {  // for all grid points
            // \bar{j} = -(\vect{v} \cdot \igrad) \phi + j (\idiv \vect{v})
            p_rhs0[i] = - ( p_vx1[i]*p_gphi1[i]
                          + p_vx2[i]*p_gphi2[i]
                          + p_vx3[i]*p_gphi3[i] )
                        + 0.5*p_phi[i]*p_divv[i]
                        + 0.5*p_divvphi[i]
                        - 0.5*( p_gphi1[i]*p_vx1[i]
                              + p_gphi2[i]*p_vx2[i]
                              + p_gphi3[i]*p_vx3[i] );

            p_phibar[i] = p_phi[i] + ht*p_rhs0[i];

            // compute \bar{phi} \vect{v}
            p_phiv1[i] = p_phibar[i]*p_vx1[i];
            p_phiv2[i] = p_phibar[i]*p_vx2[i];
            p_phiv3[i] = p_phibar[i]*p_vx3[i];
        }
} // pragma omp

        accfft_grad(p_gphi1, p_gphi2, p_gphi3, p_phibar, this->m_Opt->GetFFT().plan, &XYZ, timings);

        // compute div(\vect{v}\bar{\phi}_j)
        accfft_divergence(p_divvphi, p_phiv1, p_phiv2, p_phiv3, this->m_Opt->GetFFT().plan, timings);

#pragma omp parallel
{
        ScalarType rhs1;
#pragma omp  for
        for (IntType i=0; i < nl; ++i) {  // for all grid points
            // \bar{j} = -(\vect{v} \cdot \igrad) j + j (\idiv \vect{v})
            rhs1 = -( p_vx1[i]*p_gphi1[i]
                    + p_vx2[i]*p_gphi2[i]
                    + p_vx3[i]*p_gphi3[i] )
                    + 0.5*p_phibar[i]*p_divv[i]
                    + 0.5*p_divvphi[i]
                    - 0.5*( p_gphi1[i]*p_vx1[i]
                          + p_gphi2[i]*p_vx2[i]
                          + p_gphi3[i]*p_vx3[i] );

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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute determinant of deformation gradient; this
 * implementation first computes the deformation map and then
 * evaluates the determinant of the deformation gradient
 * based on the computed deformation map
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDetDefGradViaDispField"
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradViaDispField() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    ScalarType  *p_u1 = NULL, *p_u2 = NULL, *p_u3 = NULL, *p_phi = NULL,
                *p_gu11 = NULL, *p_gu12 = NULL, *p_gu13 = NULL,
                *p_gu21 = NULL, *p_gu22 = NULL, *p_gu23 = NULL,
                *p_gu31 = NULL, *p_gu32 = NULL, *p_gu33 = NULL;
    double timer[5] = {0, 0, 0, 0, 0};
    std::bitset<3>XYZ = 0; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

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

    // X1 gradient
    accfft_grad(p_gu11, p_gu12, p_gu13, p_u1,this->m_Opt->GetFFT().plan,&XYZ,timer);

    // X2 gradient
    accfft_grad(p_gu21, p_gu22, p_gu23, p_u2,this->m_Opt->GetFFT().plan,&XYZ,timer);

    // X3 gradient
    accfft_grad(p_gu31, p_gu32, p_gu33, p_u3,this->m_Opt->GetFFT().plan,&XYZ,timer);

    ierr = VecGetArray(this->m_WorkScaField1, &p_phi); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp  for
    for (IntType i=0; i < nl; ++i) {  // for all grid points
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

    ierr = this->m_WorkVecField1->RestoreArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_gu11, p_gu12, p_gu13); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_gu21, p_gu22, p_gu23); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->RestoreArrays(p_gu31, p_gu32, p_gu33); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDetDefGradSL"
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradSL() {
    PetscErrorCode ierr;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_divv = NULL, *p_divvX = NULL, *p_j = NULL, *p_jX=NULL;
    ScalarType ht, hthalf, velsign;
    IntType nl, ng, nt;
    std::stringstream ss;
    std::string ext;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timings[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ext = this->m_Opt->GetReadWriteFlags().extension;

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_Opt->GetRegFlags().invdefgrad) {
        velsign = -1.0;
    }
    else velsign = 1.0;

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


    // compute trajectory
    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1, "state"); CHKERRQ(ierr);

    // store time series
    if ( this->m_Opt->GetReadWriteFlags().timeseries ) {
        ss.str(std::string()); ss.clear();
        ss << "det-deformation-grad-j=" << std::setw(3) << std::setfill('0') << 0 << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkScaField1,ss.str()); CHKERRQ(ierr);
    }

    // get pointers
    ierr = VecGetArray(this->m_WorkScaField1, &p_j); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2, &p_jX); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField3, &p_divv); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField4, &p_divvX); CHKERRQ(ierr);

    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    // compute div(v)
    accfft_divergence(p_divv, p_vx1, p_vx2, p_vx3, this->m_Opt->GetFFT().plan, timings);

    // compute div(v) at X
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_divvX, p_divv, "state"); CHKERRQ(ierr);


    for (IntType j = 0; j < nt; ++j) {  // for all time points
        // compute J(X,t^j)
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_jX, p_j, "state"); CHKERRQ(ierr);

#pragma omp parallel
{
        ScalarType jX, rhs0, rhs1;
#pragma omp  for
        for (IntType i = 0; i < nl; ++i) { // for all grid points
            jX = p_jX[i];
            rhs0 = jX*velsign*p_divvX[i];
            rhs1 = (jX + ht*rhs0)*velsign*p_divv[i];
            p_j[i] = jX + hthalf*(rhs0+rhs1);
        }
}  // pragma omp

        // store time series
        if (this->m_Opt->GetReadWriteFlags().timeseries) {
            ierr = VecRestoreArray(this->m_WorkScaField1, &p_j); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
            ss << "det-deformation-grad-j=" << std::setw(3) << std::setfill('0') << j+1 << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkScaField1, ss.str()); CHKERRQ(ierr);
            ierr = VecGetArray(this->m_WorkScaField1, &p_j); CHKERRQ(ierr);
        }
    }  // for all time points

    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_WorkScaField4, &p_divvX); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField3, &p_divv); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_jX); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_j); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMap"
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMap(bool write2file) {
    PetscErrorCode ierr = 0;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // allocate velocity field
    if (this->m_VelocityField == NULL) {
       try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("computing deformation map"); CHKERRQ(ierr);
    }

    // call the solver
    switch (this->m_Opt->GetPDESolver().type) {
        case RK2:
        {
            ierr = this->ComputeDeformationMapRK2(); CHKERRQ(ierr);
            break;
        }
        case SL:
        {
            ierr = this->ComputeDeformationMapSL(); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
            break;
        }
    }

    if (write2file) {
        ext = this->m_Opt->GetReadWriteFlags().extension;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1, "deformation-map-x1"+ext,
                                                               "deformation-map-x2"+ext,
                                                               "deformation-map-x3"+ext); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMapRK2"
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapRK2() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = ThrowError("not implemented"); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMapRK2A"
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
    double timers[5]={0,0,0,0,0};
    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;
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
    accfft_divergence(p_divv, p_v1, p_v2, p_v3,this->m_Opt->GetFFT().plan,timers);


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
        accfft_grad(p_gu11, p_gu12, p_gu13, p_u1,this->m_Opt->GetFFT().plan,&XYZ,timers);
        this->m_Opt->IncrementCounter(FFT,4);

        accfft_grad(p_gu21, p_gu22, p_gu23, p_u2,this->m_Opt->GetFFT().plan,&XYZ,timers);
        this->m_Opt->IncrementCounter(FFT,4);

        accfft_grad(p_gu31, p_gu32, p_gu33, p_u3,this->m_Opt->GetFFT().plan,&XYZ,timers);
        this->m_Opt->IncrementCounter(FFT,4);

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
        accfft_divergence(p_divv, p_v1, p_v2, p_v3,this->m_Opt->GetFFT().plan,timers);


        // compute gradient of \bar{m}
        accfft_grad(p_gmx1, p_gmx2, p_gmx3, p_mbar,this->m_Opt->GetFFT().plan,&XYZ,timers);
        this->m_Opt->IncrementCounter(FFT,4);

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

    this->m_Opt->IncreaseFFTTimers(timers);

    this->m_Opt->Exit(__FUNCT__);
*/

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMapSL"
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapSL() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    switch (this->m_Opt->GetPDESolver().order) {
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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map if we consider a semi-lagrangian
 * time integrator; the scheme is full lagrangian; we use an
 * rk2 scheme to compute the characteristic;
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMapSLRK2"
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapSLRK2() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    std::string ext;
    IntType nl, nt;
    ScalarType hx[3], ht, hthalf;
    ScalarType *p_y1 = NULL, *p_y2 = NULL, *p_y3 = NULL,
                *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_vy1 = NULL, *p_vy2 = NULL, *p_vy3 = NULL,
                *p_vytilde1 = NULL, *p_vytilde2 = NULL, *p_vytilde3 = NULL,
                *p_ytilde1 = NULL, *p_ytilde2 = NULL, *p_ytilde3 = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ext = this->m_Opt->GetReadWriteFlags().extension;

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate vector fields
    if (this->m_WorkVecField1 == NULL) {
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2==NULL) {
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3==NULL) {
        try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField4==NULL) {
        try{this->m_WorkVecField4 = new VecField(this->m_Opt);}
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
    ierr = this->m_SemiLagrangianMethod->SetReadWrite(this->m_ReadWrite); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->GetDomainPara().hx[i];
    }

    // compute initial condition
    ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);

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
                p_y1[l] = hx[0]*static_cast<ScalarType>(i1 + this->m_Opt->GetDomainPara().istart[0]);
                p_y2[l] = hx[1]*static_cast<ScalarType>(i2 + this->m_Opt->GetDomainPara().istart[1]);
                p_y3[l] = hx[2]*static_cast<ScalarType>(i3 + this->m_Opt->GetDomainPara().istart[2]);
            } // i1
        } // i2
    } // i3
}// pragma omp for

    ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);

    // store time series
    if (this->m_Opt->GetReadWriteFlags().timeseries ) {
        ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

        // write out y1
        ss.str(std::string()); ss.clear();
        ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << 0 << "-x1" << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X1,ss.str()); CHKERRQ(ierr);

        ss.str(std::string()); ss.clear();
        ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << 0 << "-x2" << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X2,ss.str()); CHKERRQ(ierr);

        ss.str(std::string()); ss.clear();
        ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << 0 << "-x3" << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X3,ss.str()); CHKERRQ(ierr);
    }


    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;
    ht = this->m_Opt->GetTimeStepSize();
//    if (inverse) { ht *= -1.0; }
    hthalf = 0.5*ht;

    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->GetArrays(p_vy1, p_vy2, p_vy3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->GetArrays(p_vytilde1, p_vytilde2, p_vytilde3); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {
        // evaluate v(y)
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3,
                                                         p_v1, p_v2, p_v3,
                                                         p_y1, p_y2, p_y3 ); CHKERRQ(ierr);

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
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vytilde1, p_vytilde2, p_vytilde3,
                                                         p_v1, p_v2, p_v3,
                                                         p_ytilde1, p_ytilde2, p_ytilde3 ); CHKERRQ(ierr);

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

            // write out y1
            ss.str(std::string()); ss.clear();
            ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << j+1 << "-x1" << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X1,ss.str()); CHKERRQ(ierr);

            ss.str(std::string()); ss.clear();
            ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << j+1 << "-x2" << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X2,ss.str()); CHKERRQ(ierr);

            ss.str(std::string()); ss.clear();
            ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << j+1 << "-x3" << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X3,ss.str()); CHKERRQ(ierr);

            ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
        }
    } // for all time points

    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_vy1, p_vy2, p_vy3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->RestoreArrays(p_vytilde1, p_vytilde2, p_vytilde3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute deformation map if we consider a semi-lagrangian
 * time integrator; the scheme is full lagrangian; we use an
 * rk4 scheme to compute the characteristic;
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMapSLRK4"
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapSLRK4() {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    std::string ext;
    IntType nl,nt;
    ScalarType hx[3],ht,hthalf,htby6;
    ScalarType *p_y1 = NULL, *p_y2 = NULL, *p_y3 = NULL,
                *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_vy1 = NULL, *p_vy2 = NULL, *p_vy3 = NULL,
                *p_dy1 = NULL, *p_dy2 = NULL, *p_dy3 = NULL,
                *p_ytilde1 = NULL, *p_ytilde2 = NULL, *p_ytilde3 = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ext = this->m_Opt->GetReadWriteFlags().extension;

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
    if(this->m_SemiLagrangianMethod == NULL) {
        try{this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    ierr = this->m_SemiLagrangianMethod->SetReadWrite(this->m_ReadWrite); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->GetDomainPara().hx[i];
    }

    // compute initial condition y = x
    ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
#pragma omp parallel
{
    IntType l, i1, i2, i3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->GetDomainPara().isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < this->m_Opt->GetDomainPara().isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < this->m_Opt->GetDomainPara().isize[2]; ++i3) {  // x3
                // compute linear / flat index
                l = GetLinearIndex(i1,i2,i3,this->m_Opt->GetDomainPara().isize);

                // compute coordinates (nodal grid)
                p_y1[l] = hx[0]*static_cast<ScalarType>(i1 + this->m_Opt->GetDomainPara().istart[0]);
                p_y2[l] = hx[1]*static_cast<ScalarType>(i2 + this->m_Opt->GetDomainPara().istart[1]);
                p_y3[l] = hx[2]*static_cast<ScalarType>(i3 + this->m_Opt->GetDomainPara().istart[2]);
            }  // i1
        }  // i2
    }  // i3
}  // pragma omp for
    ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);

    // store time series
    if (this->m_Opt->GetReadWriteFlags().timeseries ) {
        ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

        // write out y1
        ss.str(std::string()); ss.clear();
        ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << 0 << "-x1" << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X1,ss.str()); CHKERRQ(ierr);

        ss.str(std::string()); ss.clear();
        ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << 0 << "-x2" << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X2,ss.str()); CHKERRQ(ierr);

        ss.str(std::string()); ss.clear();
        ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << 0 << "-x3" << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X3,ss.str()); CHKERRQ(ierr);
    }

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;
    ht = this->m_Opt->GetTimeStepSize();
//    if (inverse) { ht *= -1.0; }
    hthalf = 0.5*ht;
    htby6  = ht/6.0;

    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->GetArrays(p_vy1, p_vy2, p_vy3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->GetArrays(p_dy1, p_dy2, p_dy3); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {

        // evaluate right hand side v(y) (i.e., F0)
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3, p_v1, p_v2, p_v3,
                                                         p_y1, p_y2, p_y3); CHKERRQ(ierr);

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
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3, p_v1, p_v2, p_v3,
                                                         p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);

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
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3, p_v1, p_v2, p_v3,
                                                         p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);

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
        ierr = this->m_SemiLagrangianMethod->Interpolate(p_vy1, p_vy2, p_vy3, p_v1, p_v2, p_v3,
                                                         p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);

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
        if (this->m_Opt->GetReadWriteFlags().timeseries ) {
            ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

            ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);

            // write out y1
            ss.str(std::string()); ss.clear();
            ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << j+1 << "-x1" << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X1,ss.str()); CHKERRQ(ierr);

            ss.str(std::string()); ss.clear();
            ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << j+1 << "-x2" << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X2,ss.str()); CHKERRQ(ierr);

            ss.str(std::string()); ss.clear();
            ss << "deformation-map-j=" << std::setw(3) << std::setfill('0') << j+1 << "-x3" << ext;
            ierr = this->m_ReadWrite->Write(this->m_WorkVecField1->m_X3,ss.str()); CHKERRQ(ierr);

            ierr = this->m_WorkVecField1->GetArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
        }
    }  // for all time points

    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(p_y1, p_y2, p_y3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_ytilde1, p_ytilde2, p_ytilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_vy1, p_vy2, p_vy3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->RestoreArrays(p_dy1, p_dy2, p_dy3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute displacement field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDisplacementField"
PetscErrorCode OptimalControlRegistrationBase::ComputeDisplacementField(bool write2file) {
    PetscErrorCode ierr = 0;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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
    switch (this->m_Opt->GetPDESolver().type) {
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
        ext = this->m_Opt->GetReadWriteFlags().extension;
        ierr = this->m_ReadWrite->Write(this->m_WorkVecField1,"displacement-field-x1"+ext,
                                                              "displacement-field-x2"+ext,
                                                              "displacement-field-x3"+ext); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute displacement field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDisplacementFieldRK2"
PetscErrorCode OptimalControlRegistrationBase::ComputeDisplacementFieldRK2() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = ThrowError("not implemented"); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute displacement field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDisplacementFieldSL"
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

    this->m_Opt->Enter(__FUNCT__);

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
    ierr = this->m_SemiLagrangianMethod->SetReadWrite(this->m_ReadWrite); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1, "state"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;
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
}// end of pragma omp parallel
        ierr = this->m_WorkVecField3->RestoreArrays(p_uX1, p_uX2, p_uX3); CHKERRQ(ierr);
        ierr = this->m_WorkVecField1->RestoreArrays(p_u1, p_u2, p_u3); CHKERRQ(ierr);

    } // for all time points
    ierr = this->m_WorkVecField2->RestoreArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief compute deformation map from a displacement field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDefMapFromDisplacement"
PetscErrorCode OptimalControlRegistrationBase::ComputeDefMapFromDisplacement() {
    PetscErrorCode ierr = 0;
    ScalarType hx[3];
    ScalarType *p_u1 = NULL,*p_u2 = NULL, *p_u3 = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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
    switch (this->m_Opt->GetPDESolver().type) {
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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif // _OPTIMALCONTROLREGISTRATIONBASE_CPP_
