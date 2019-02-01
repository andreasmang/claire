/*************************************************************************
 *  Copyright (c) 2018.
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

#ifndef _DIFFERENTIATIONFD_CPP_
#define _DIFFERENTIATIONFD_CPP_

#include "DifferentiationFD.hpp"
#include "TextureDifferentiationKernel.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DifferentiationFD::DifferentiationFD() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DifferentiationFD::DifferentiationFD(RegOpt* opt) : SuperClass(opt, Type::Finite) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DifferentiationFD::~DifferentiationFD() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DifferentiationFD::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->SetupData(); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DifferentiationFD::SetupData(ScalarType *x1, ScalarType *x2, ScalarType *x3) {
    PetscErrorCode ierr = 0;
    IntType nalloc;
    PetscFunctionBegin;

    this->mtex = gpuInitEmptyGradientTexture(this->m_Opt->m_Domain.nx);
    ierr = initConstants(this->m_Opt->m_Domain.nx); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DifferentiationFD::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->mtex != 0) {
        cudaDestroyTextureObject(this->mtex);
    }
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Gradient(ScalarType *g1,
                                           ScalarType *g2,
                                           ScalarType *g3,
                                           const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = computeTextureGradient(g1, g2, g3, m, this->mtex, this->m_Opt->m_Domain.nx); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Gradient(VecField *g, const ScalarType *m) {
    PetscErrorCode ierr = 0;
    ScalarType *g1 = nullptr, *g2 = nullptr, *g3 = nullptr;
    PetscFunctionBegin;
    
    ierr = g->GetArraysWrite(g1, g2, g3); CHKERRQ(ierr);
    
    ierr = this->Gradient(g1, g2, g3, m); CHKERRQ(ierr);
    
    ierr = g->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Gradient(ScalarType **g, const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = this->Gradient(g[0], g[1], g[2], m); CHKERRQ(ierr);
            
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Gradient(VecField *g, const Vec m) {
    PetscErrorCode ierr = 0;
    ScalarType *g1 = nullptr, *g2 = nullptr, *g3 = nullptr;
    const ScalarType *pm = nullptr;
    PetscFunctionBegin;
    
    ierr = g->GetArraysWrite(g1, g2, g3); CHKERRQ(ierr);
    ierr = GetRawPointerRead(m, &pm); CHKERRQ(ierr);
    
    ierr = this->Gradient(g1, g2, g3, pm); CHKERRQ(ierr);
    
    ierr = RestoreRawPointerRead(m, &pm); CHKERRQ(ierr);
    ierr = g->RestoreArrays(); CHKERRQ(ierr);
            
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute laplacian of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Laplacian(ScalarType *l,
                                            const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = DebugNotImplemented(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute laplacian of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Laplacian(VecField *l,
                                            VecField *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = DebugNotImplemented(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute laplacian of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Laplacian(Vec l,
                                            const Vec m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = DebugNotImplemented(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute laplacian of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Laplacian(ScalarType *l1,
                                            ScalarType *l2,
                                            ScalarType *l3,
                                            const ScalarType *v1,
                                            const ScalarType *v2,
                                            const ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = DebugNotImplemented(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}
 
/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Divergence(ScalarType *l,
                                             const ScalarType *v1,
                                             const ScalarType *v2,
                                             const ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = computeTextureDivergence(l, v1, v2, v3, this->mtex, this->m_Opt->m_Domain.nx); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Divergence(ScalarType *l, VecField *v) {
    PetscErrorCode ierr = 0;
    const ScalarType *v1 = nullptr, *v2 = nullptr, *v3 = nullptr;
    PetscFunctionBegin;
    
    ierr = v->GetArraysRead(v1, v2, v3); CHKERRQ(ierr);
    
    ierr = this->Divergence(l, v1, v2, v3); CHKERRQ(ierr);
    
    ierr = v->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Divergence(Vec l, VecField *v) {
    PetscErrorCode ierr = 0;
    const ScalarType *v1 = nullptr, *v2 = nullptr, *v3 = nullptr;
    ScalarType *pl;
    PetscFunctionBegin;
    
    ierr = v->GetArraysRead(v1, v2, v3); CHKERRQ(ierr);
    ierr = GetRawPointerWrite(l, &pl); CHKERRQ(ierr);
    
    ierr = this->Divergence(pl, v1, v2, v3); CHKERRQ(ierr);
    
    ierr = RestoreRawPointerWrite(l, &pl); CHKERRQ(ierr);
    ierr = v->RestoreArrays(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode DifferentiationFD::Divergence(ScalarType *l, const ScalarType *const *v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
        
    ierr = this->Divergence(l, v[0], v[1], v[2]); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationFD::RegLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::RegBiLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::RegTriLapOp(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::RegTriLapFunc(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationFD::InvRegLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::InvRegBiLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
PetscErrorCode DifferentiationFD::InvRegTriLapOp(VecField* bv, VecField* v, bool usesqrt, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode DifferentiationFD::LerayOperator(VecField* bv, VecField* v, ScalarType b0, ScalarType b1) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  ierr = DebugNotImplemented(); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}




}  // end of namespace




#endif  // _DIFFERENTIATION_CPP_
