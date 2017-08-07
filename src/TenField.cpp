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

#ifndef _TENFIELD_CPP_
#define _TENFIELD_CPP_

#include "TenField.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
TenField::TenField() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
TenField::~TenField() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
TenField::TenField(RegOpt* opt) {
    this->Initialize();
    this->SetOpt(opt);
    this->Allocate();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
TenField::TenField(RegOpt* opt, int level) {
    this->Initialize();
    this->SetOpt(opt);
    this->Allocate(level);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
TenField::TenField(IntType nl, IntType ng) {
    this->Initialize();
    this->Allocate(nl, ng);

}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode TenField::Initialize(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;

    this->m_X11 = NULL;
    this->m_X12 = NULL;
    this->m_X13 = NULL;
    this->m_X21 = NULL;
    this->m_X22 = NULL;
    this->m_X23 = NULL;
    this->m_X31 = NULL;
    this->m_X32 = NULL;
    this->m_X33 = NULL;

    PetscFunctionReturn(ierr);

}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode TenField::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_X11 != NULL) {ierr = VecDestroy(&this->m_X11); CHKERRQ(ierr); this->m_X11 = NULL;}
    if (this->m_X12 != NULL) {ierr = VecDestroy(&this->m_X12); CHKERRQ(ierr); this->m_X12 = NULL;}
    if (this->m_X13 != NULL) {ierr = VecDestroy(&this->m_X13); CHKERRQ(ierr); this->m_X13 = NULL;}

    if (this->m_X21 != NULL) {ierr = VecDestroy(&this->m_X21); CHKERRQ(ierr); this->m_X21 = NULL;}
    if (this->m_X22 != NULL) {ierr = VecDestroy(&this->m_X22); CHKERRQ(ierr); this->m_X22 = NULL;}
    if (this->m_X23 != NULL) {ierr = VecDestroy(&this->m_X23); CHKERRQ(ierr); this->m_X23 = NULL;}

    if (this->m_X31 != NULL) {ierr = VecDestroy(&this->m_X31); CHKERRQ(ierr); this->m_X31 = NULL;}
    if (this->m_X32 != NULL) {ierr = VecDestroy(&this->m_X32); CHKERRQ(ierr); this->m_X32 = NULL;}
    if (this->m_X33 != NULL) {ierr = VecDestroy(&this->m_X33); CHKERRQ(ierr); this->m_X33 = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
PetscErrorCode TenField::SetOpt(RegOpt* opt) {
    PetscErrorCode ierr = 0;

    ierr = Assert(opt != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_Opt = opt;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
PetscErrorCode TenField::Allocate() {
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
PetscErrorCode TenField::Allocate(int level) {
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
PetscErrorCode TenField::Allocate(IntType nl, IntType ng) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr = this->ClearMemory(); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X11); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X11, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X11); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X12); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X12, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X12); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X13); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X13, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X13); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X21); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X21, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X21); CHKERRQ(ierr);

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X22); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X22, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X22); CHKERRQ(ierr);

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X23); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X23, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X23); CHKERRQ(ierr);

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X31); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X31, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X31); CHKERRQ(ierr);

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X32); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X32, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X32); CHKERRQ(ierr);

    // allocate vector field
    ierr = VecCreate(PETSC_COMM_WORLD, &this->m_X33); CHKERRQ(ierr);
    ierr = VecSetSizes(this->m_X33, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(this->m_X33); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}






/********************************************************************
 * @brief Copy
 *******************************************************************/
PetscErrorCode TenField::Copy(TenField* t) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecCopy(t->m_X11, this->m_X11); CHKERRQ(ierr);
    ierr = VecCopy(t->m_X12, this->m_X12); CHKERRQ(ierr);
    ierr = VecCopy(t->m_X13, this->m_X13); CHKERRQ(ierr);

    ierr = VecCopy(t->m_X21, this->m_X21); CHKERRQ(ierr);
    ierr = VecCopy(t->m_X22, this->m_X22); CHKERRQ(ierr);
    ierr = VecCopy(t->m_X23, this->m_X23); CHKERRQ(ierr);

    ierr = VecCopy(t->m_X31, this->m_X31); CHKERRQ(ierr);
    ierr = VecCopy(t->m_X32, this->m_X32); CHKERRQ(ierr);
    ierr = VecCopy(t->m_X33, this->m_X33); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set value
 *******************************************************************/
PetscErrorCode TenField::SetValue(ScalarType value) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecSet(this->m_X11, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X12, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X13, value); CHKERRQ(ierr);

    ierr = VecSet(this->m_X21, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X22, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X23, value); CHKERRQ(ierr);

    ierr = VecSet(this->m_X31, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X32, value); CHKERRQ(ierr);
    ierr = VecSet(this->m_X33, value); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
PetscErrorCode TenField::GetArrays(ScalarType*& p_x11, ScalarType*& p_x12, ScalarType*& p_x13,
                                   ScalarType*& p_x21, ScalarType*& p_x22, ScalarType*& p_x23,
                                   ScalarType*& p_x31, ScalarType*& p_x32, ScalarType*& p_x33) {
    PetscErrorCode ierr = 0;

    ierr = VecGetArray(this->m_X11, &p_x11); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_X12, &p_x12); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_X13, &p_x13); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_X21, &p_x21); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_X22, &p_x22); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_X23, &p_x23); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_X31, &p_x31); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_X32, &p_x32); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_X33, &p_x33); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get arrays of vector field
 *******************************************************************/
PetscErrorCode TenField::GetArraysRead(const ScalarType*& p_x11, const ScalarType*& p_x12, const ScalarType*& p_x13,
                                       const ScalarType*& p_x21, const ScalarType*& p_x22, const ScalarType*& p_x23,
                                       const ScalarType*& p_x31, const ScalarType*& p_x32, const ScalarType*& p_x33) {
    PetscErrorCode ierr = 0;

    ierr = VecGetArrayRead(this->m_X11, &p_x11); CHKERRQ(ierr);
    ierr = VecGetArrayRead(this->m_X12, &p_x12); CHKERRQ(ierr);
    ierr = VecGetArrayRead(this->m_X13, &p_x13); CHKERRQ(ierr);

    ierr = VecGetArrayRead(this->m_X21, &p_x21); CHKERRQ(ierr);
    ierr = VecGetArrayRead(this->m_X22, &p_x22); CHKERRQ(ierr);
    ierr = VecGetArrayRead(this->m_X23, &p_x23); CHKERRQ(ierr);

    ierr = VecGetArrayRead(this->m_X31, &p_x31); CHKERRQ(ierr);
    ierr = VecGetArrayRead(this->m_X32, &p_x32); CHKERRQ(ierr);
    ierr = VecGetArrayRead(this->m_X33, &p_x33); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief pointwise scale of a vector field
 *******************************************************************/
PetscErrorCode TenField::RestoreArrays(ScalarType*& p_x11,ScalarType*& p_x12,ScalarType*& p_x13,
                                       ScalarType*& p_x21,ScalarType*& p_x22,ScalarType*& p_x23,
                                       ScalarType*& p_x31,ScalarType*& p_x32,ScalarType*& p_x33) {
    PetscErrorCode ierr = 0;

    ierr = VecRestoreArray(this->m_X11, &p_x11); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_X12, &p_x12); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_X13, &p_x13); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_X21, &p_x21); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_X22, &p_x22); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_X23, &p_x23); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_X31, &p_x31); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_X32, &p_x32); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_X33, &p_x33); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief restore arrays of vector field
 *******************************************************************/
PetscErrorCode TenField::RestoreArraysRead(const ScalarType*& p_x11, const ScalarType*& p_x12, const ScalarType*& p_x13,
                                           const ScalarType*& p_x21, const ScalarType*& p_x22, const ScalarType*& p_x23,
                                           const ScalarType*& p_x31, const ScalarType*& p_x32, const ScalarType*& p_x33) {
    PetscErrorCode ierr = 0;

    ierr = VecRestoreArrayRead(this->m_X11, &p_x11); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(this->m_X12, &p_x12); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(this->m_X13, &p_x13); CHKERRQ(ierr);

    ierr = VecRestoreArrayRead(this->m_X21,&p_x21); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(this->m_X22,&p_x22); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(this->m_X23,&p_x23); CHKERRQ(ierr);

    ierr = VecRestoreArrayRead(this->m_X31,&p_x31); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(this->m_X32,&p_x32); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(this->m_X33,&p_x33); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief sets the individual components of a vector field;
 * the input is a flat petsc array
 *******************************************************************/
PetscErrorCode TenField::SetComponents(Vec w) {
    PetscErrorCode ierr = 0;
    IntType nl, n;
    ScalarType *p_x11 = NULL, *p_x12 = NULL, *p_x13 = NULL,
                *p_x21 = NULL, *p_x22 = NULL, *p_x23 = NULL,
                *p_x31 = NULL, *p_x32 = NULL, *p_x33 = NULL;
    const ScalarType *p_w = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(w, &n); CHKERRQ(ierr);

    ierr = Assert(this->m_X11 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X12 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X13 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X21 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X22 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X23 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X31 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X32 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X33 != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = VecGetArrayRead(w, &p_w); CHKERRQ(ierr);
    ierr = this->GetArrays(p_x11, p_x12, p_x13,
                           p_x21, p_x22, p_x23,
                           p_x31, p_x32, p_x33); CHKERRQ(ierr);

    //compute size of each individual component
    nl = n / 9;
    //ierr = Assert(nl==this->m_Opt->m_Domain.nl,"dimension mismatch"); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i){
        p_x11[i] = p_w[i+0*nl];
        p_x12[i] = p_w[i+1*nl];
        p_x13[i] = p_w[i+2*nl];

        p_x21[i] = p_w[i+3*nl];
        p_x22[i] = p_w[i+4*nl];
        p_x23[i] = p_w[i+5*nl];

        p_x31[i] = p_w[i+6*nl];
        p_x32[i] = p_w[i+7*nl];
        p_x33[i] = p_w[i+8*nl];

    }
} // pragma omp parallel


    ierr = VecRestoreArrayRead(w, &p_w); CHKERRQ(ierr);
    ierr = this->RestoreArrays(p_x11, p_x12, p_x13,
                               p_x21, p_x22, p_x23,
                               p_x31, p_x32, p_x33); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief sets the individual components of a vector field;
 * the input is a flat petsc array
 *******************************************************************/
PetscErrorCode TenField::GetComponents(Vec w) {
    PetscErrorCode ierr = 0;
    IntType nl, n;
    const ScalarType *p_x11 = NULL, *p_x12 = NULL, *p_x13 = NULL,
                     *p_x21 = NULL, *p_x22 = NULL, *p_x23 = NULL,
                     *p_x31 = NULL, *p_x32 = NULL, *p_x33 = NULL;
    ScalarType *p_w = NULL;

    PetscFunctionBegin;

    // get local size of vector field
    ierr = VecGetLocalSize(w, &n); CHKERRQ(ierr);

    ierr = Assert(this->m_X11 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X12 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X13 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X21 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X22 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X23 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X31 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X32 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_X33 != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = VecGetArray(w, &p_w); CHKERRQ(ierr);
    ierr = this->GetArraysRead(p_x11, p_x12, p_x13,
                               p_x21, p_x22, p_x23,
                               p_x31, p_x32, p_x33); CHKERRQ(ierr);

    //compute size of each individual component
    nl = n / 9;
    //ierr = Assert(nl==this->m_Opt->m_Domain.nl,"dimension mismatch"); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i){
       p_w[i+0*nl] = p_x11[i];
       p_w[i+1*nl] = p_x12[i];
       p_w[i+2*nl] = p_x13[i];

       p_w[i+3*nl] = p_x21[i];
       p_w[i+4*nl] = p_x22[i];
       p_w[i+5*nl] = p_x23[i];

       p_w[i+6*nl] = p_x31[i];
       p_w[i+7*nl] = p_x32[i];
       p_w[i+8*nl] = p_x33[i];

    }
} // pragma omp parallel

    ierr = VecRestoreArray(w, &p_w); CHKERRQ(ierr);
    ierr = this->RestoreArraysRead(p_x11, p_x12, p_x13,
                                   p_x21, p_x22, p_x23,
                                   p_x31, p_x32, p_x33); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief scale tensor field by scalar value
 *******************************************************************/
PetscErrorCode TenField::Scale(ScalarType value) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecScale(this->m_X11, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X12, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X13, value); CHKERRQ(ierr);

    ierr = VecScale(this->m_X21, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X22, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X23, value); CHKERRQ(ierr);

    ierr = VecScale(this->m_X31, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X32, value); CHKERRQ(ierr);
    ierr = VecScale(this->m_X33, value); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief pointwise scale of vector field
 *******************************************************************/
PetscErrorCode TenField::SetIdentity() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecSet(this->m_X11, 1.0); CHKERRQ(ierr);
    ierr = VecSet(this->m_X12, 0.0); CHKERRQ(ierr);
    ierr = VecSet(this->m_X13, 0.0); CHKERRQ(ierr);

    ierr = VecSet(this->m_X21, 0.0); CHKERRQ(ierr);
    ierr = VecSet(this->m_X22, 1.0); CHKERRQ(ierr);
    ierr = VecSet(this->m_X23, 0.0); CHKERRQ(ierr);

    ierr = VecSet(this->m_X31, 0.0); CHKERRQ(ierr);
    ierr = VecSet(this->m_X32, 0.0); CHKERRQ(ierr);
    ierr = VecSet(this->m_X33, 1.0); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief pointwise scale of vector field
 *******************************************************************/
PetscErrorCode TenField::Scale(Vec s) {
    PetscErrorCode ierr;
    IntType nl;
    ScalarType *p_x11 = NULL, *p_x12 = NULL, *p_x13 = NULL,
               *p_x21 = NULL, *p_x22 = NULL, *p_x23 = NULL,
               *p_x31 = NULL, *p_x32 = NULL, *p_x33 = NULL, *p_s = NULL;

    PetscFunctionBegin;

    // get pointers
    ierr = VecGetArray(s, &p_s); CHKERRQ(ierr);

    // get local size of vector field
    ierr = VecGetLocalSize(s, &nl); CHKERRQ(ierr);
    ierr = this->GetArrays(p_x11, p_x12, p_x13,
                           p_x21, p_x22, p_x23,
                           p_x31, p_x32, p_x33); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i) {
        ScalarType scale = p_s[i];
        p_x11[i] *= scale;
        p_x12[i] *= scale;
        p_x13[i] *= scale;
        p_x21[i] *= scale;
        p_x22[i] *= scale;
        p_x23[i] *= scale;
        p_x31[i] *= scale;
        p_x32[i] *= scale;
        p_x33[i] *= scale;
    }
} // pragma omp parallel

    // get pointers
    ierr = VecRestoreArray(s, &p_s); CHKERRQ(ierr);
    ierr = this->RestoreArrays(p_x11, p_x12, p_x13,
                               p_x21, p_x22, p_x23,
                               p_x31, p_x32, p_x33); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief interface for AXPY
 *******************************************************************/
PetscErrorCode TenField::AXPY(ScalarType s,TenField* t) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecAXPY(this->m_X11, s, t->m_X11); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_X12, s, t->m_X12); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_X13, s, t->m_X13); CHKERRQ(ierr);

    ierr = VecAXPY(this->m_X21, s, t->m_X21); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_X22, s, t->m_X22); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_X23, s, t->m_X23); CHKERRQ(ierr);

    ierr = VecAXPY(this->m_X31, s, t->m_X31); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_X32, s, t->m_X32); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_X33, s, t->m_X33); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief interface for WAXPY
 *******************************************************************/
PetscErrorCode TenField::WAXPY(ScalarType s, TenField* tv, TenField* tw) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = VecWAXPY(this->m_X11, s,tv->m_X11, tw->m_X11); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X12, s,tv->m_X12, tw->m_X12); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X13, s,tv->m_X13, tw->m_X13); CHKERRQ(ierr);

    ierr = VecWAXPY(this->m_X21, s, tv->m_X21, tw->m_X21); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X22, s, tv->m_X22, tw->m_X22); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X23, s, tv->m_X23, tw->m_X23); CHKERRQ(ierr);

    ierr = VecWAXPY(this->m_X31, s, tv->m_X31, tw->m_X31); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X32, s, tv->m_X32, tw->m_X32); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X33, s, tv->m_X33, tw->m_X33); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _TENFIELD_CPP_
