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

#ifndef _VECFIELD_CPP_
#define _VECFIELD_CPP_

#include "VecField.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecField"
VecField::VecField() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~VecField"
VecField::~VecField() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecField"
VecField::VecField(RegOpt* opt) {
    this->Initialize();
    this->SetOpt(opt);
    this->Allocate();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecField"
VecField::VecField(RegOpt* opt, int level) {
    this->Initialize();
    this->SetOpt(opt);
    this->Allocate(level);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecField"
VecField::VecField(IntType nl, IntType ng) {
    this->Initialize();
    this->Allocate(nl, ng);
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode VecField::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;

    this->m_X1 = NULL;
    this->m_X2 = NULL;
    this->m_X3 = NULL;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode VecField::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

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

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set the options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetOpt"
PetscErrorCode VecField::SetOpt(RegOpt* opt) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(opt != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_Opt = opt;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode VecField::Allocate() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr = this->ClearMemory(); CHKERRQ(ierr);

    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    ierr = this->Allocate(nl, ng); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode VecField::Allocate(int level) {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr = this->ClearMemory(); CHKERRQ(ierr);

    nl = this->m_Opt->GetGridContPara().nlocal[level];
    ng = this->m_Opt->GetGridContPara().nglobal[level];

    ierr = this->Allocate(nl, ng); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode VecField::Allocate(IntType nl, IntType ng) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr = this->ClearMemory(); CHKERRQ(ierr);

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X1); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X1, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X1); CHKERRQ(ierr);

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X2); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X2, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X2); CHKERRQ(ierr);

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X3); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X3, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief Copy
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Copy"
PetscErrorCode VecField::Copy(VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecCopy(v->m_X1, this->m_X1); CHKERRQ(ierr);
    ierr = VecCopy(v->m_X2, this->m_X2); CHKERRQ(ierr);
    ierr = VecCopy(v->m_X3, this->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set value
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetValue"
PetscErrorCode VecField::SetValue(ScalarType value) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecSet(this->m_X1, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X2, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X3, value); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetArrays"
PetscErrorCode VecField::GetArrays(ScalarType*& p_x1,
                                   ScalarType*& p_x2,
                                   ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;

    ierr = VecGetArray(this->m_X1, &p_x1); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_X2, &p_x2); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_X3, &p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetArrays"
PetscErrorCode VecField::GetArraysRead(const ScalarType*& p_x1,
                                       const ScalarType*& p_x2,
                                       const ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;

    ierr = VecGetArrayRead(this->m_X1, &p_x1); CHKERRQ(ierr);
    ierr = VecGetArrayRead(this->m_X2, &p_x2); CHKERRQ(ierr);
    ierr = VecGetArrayRead(this->m_X3, &p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief pointwise scale of a vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RestoreArrays"
PetscErrorCode VecField::RestoreArrays(ScalarType*& p_x1,
                                       ScalarType*& p_x2,
                                       ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;

    ierr = VecRestoreArray(this->m_X1, &p_x1); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_X2, &p_x2); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_X3, &p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief restore arrays of vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetArrays"
PetscErrorCode VecField::RestoreArraysRead(const ScalarType*& p_x1,
                                           const ScalarType*& p_x2,
                                           const ScalarType*& p_x3) {
    PetscErrorCode ierr = 0;

    ierr = VecRestoreArrayRead(this->m_X1, &p_x1); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(this->m_X2, &p_x2); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(this->m_X3, &p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief sets the individual components of a vector field;
 * the input is a flat petsc array
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetComponents"
PetscErrorCode VecField::SetComponents(Vec w) {
    PetscErrorCode ierr = 0;
    IntType nl, n;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL;
    const ScalarType *p_w = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(w, &n); CHKERRQ(ierr);

    ierr = Assert(this->m_X1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X3 != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = VecGetArrayRead(w, &p_w); CHKERRQ(ierr);
    ierr = this->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    // compute size of each individual component
    nl = n / 3;
//    ierr = Assert(nl==this->m_Opt->GetDomainPara().nlocal,"dimension mismatch"); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i) {
        p_x1[i] = p_w[i     ];
        p_x2[i] = p_w[i+  nl];
        p_x3[i] = p_w[i+2*nl];
    }
}  // pragma omp parallel

    ierr = VecRestoreArrayRead(w, &p_w); CHKERRQ(ierr);
    ierr = this->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get components of vector field and store them
 * in a flat vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetComponents"
PetscErrorCode VecField::GetComponents(Vec w) {
    PetscErrorCode ierr = 0;
    IntType nl, n;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL, *p_w = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(w, &n); CHKERRQ(ierr);

    ierr = VecGetArray(w, &p_w); CHKERRQ(ierr);
    ierr = this->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    // compute size of each individual component
    nl = n / 3;

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i) {
        p_w[i     ] = p_x1[i];
        p_w[i+  nl] = p_x2[i];
        p_w[i+2*nl] = p_x3[i];
    }
}  // pragma omp parallel

    ierr = this->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
    ierr = VecRestoreArray(w, &p_w); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief scale vector by scalar value
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Scale"
PetscErrorCode VecField::Scale(ScalarType value) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecScale(this->m_X1, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X2, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X3, value); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief pointwise scale of vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Scale"
PetscErrorCode VecField::Scale(Vec s) {
    PetscErrorCode ierr = 0;
    IntType nl;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL, *p_s = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(s, &nl); CHKERRQ(ierr);

    // get pointers
    ierr = VecGetArray(s, &p_s); CHKERRQ(ierr);
    ierr = this->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

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
    ierr = VecRestoreArray(s, &p_s); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief pointwise scale of a vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Scale"
PetscErrorCode VecField::Scale(VecField* v, Vec s) {
    PetscErrorCode ierr = 0;
    IntType nl;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL, *p_s = NULL,
                *p_sv1 = NULL, *p_sv2 = NULL, *p_sv3 = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(s, &nl); CHKERRQ(ierr);

    // get pointers
    ierr = VecGetArray(s, &p_s); CHKERRQ(ierr);
    ierr = this->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = v->GetArrays(p_sv1, p_sv2, p_sv3); CHKERRQ(ierr);

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
    ierr = VecRestoreArray(s, &p_s); CHKERRQ(ierr);
    ierr = this->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = v->RestoreArrays(p_sv1, p_sv2, p_sv3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief interface for AXPY
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "AXPY"
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
#undef __FUNCT__
#define __FUNCT__ "AXPY"
PetscErrorCode VecField::WAXPY(ScalarType s, VecField* v, VecField* w) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    ierr = VecWAXPY(this->m_X1, s, v->m_X1, w->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X2, s, v->m_X2, w->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X3, s, v->m_X3, w->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief compute norm of vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Scale"
PetscErrorCode VecField::Norm(Vec xnorm) {
    PetscErrorCode ierr = 0;
    IntType i, nl;
    ScalarType *p_x1 = NULL, *p_x2 = NULL, *p_x3 = NULL, *p_x = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(xnorm, &nl); CHKERRQ(ierr);

    ierr = this->GetArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);
    ierr = VecGetArray(xnorm, &p_x); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (i = 0; i < nl; ++i) {
        p_x[i] = PetscSqrtReal(p_x1[i]*p_x1[i]
                             + p_x2[i]*p_x2[i]
                             + p_x3[i]*p_x3[i]);
    }
}  // pragma omp parallel

    ierr = VecRestoreArray(xnorm, &p_x); CHKERRQ(ierr);
    ierr = this->RestoreArrays(p_x1, p_x2, p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _VECFIELD_CPP_
