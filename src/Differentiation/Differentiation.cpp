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

#ifndef _DIFFERENTIATION_CPP_
#define _DIFFERENTIATION_CPP_

#include "Differentiation.hpp"




namespace reg {





/********************************************************************
 * @brief default constructor
 *******************************************************************/
Differentiation::Differentiation() : m_Type(Type::None) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
Differentiation::~Differentiation() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
Differentiation::Differentiation(RegOpt* opt, Type type) : m_Type(type) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode Differentiation::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode Differentiation::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode Differentiation::Gradient(VecField *g, const ScalarType *m) {
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
PetscErrorCode Differentiation::Gradient(ScalarType **g, const ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    ierr = this->Gradient(g[0], g[1], g[2], m); CHKERRQ(ierr);
            
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode Differentiation::Gradient(VecField *g, const Vec m) {
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
 * @brief compute laplacian of a vector field
 *******************************************************************/
PetscErrorCode Differentiation::Laplacian(VecField *l,
                                            VecField *m) {
    PetscErrorCode ierr = 0;
    ScalarType *l1 = nullptr, *l2 = nullptr, *l3 = nullptr;
    const ScalarType *m1 = nullptr, *m2 = nullptr, *m3 = nullptr;
    PetscFunctionBegin;
    
    ierr = m->GetArraysRead(m1, m2, m3); CHKERRQ(ierr);
    ierr = l->GetArraysWrite(l1, l2, l3); CHKERRQ(ierr);
    
    ierr = this->Laplacian(l1, l2, l3, m1, m2, m3); CHKERRQ(ierr);

    ierr = l->RestoreArrays(); CHKERRQ(ierr);
    ierr = m->RestoreArrays(); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute divergence of a scalar field
 *******************************************************************/
PetscErrorCode Differentiation::Laplacian(Vec l, const Vec v) {
    PetscErrorCode ierr = 0;
    const ScalarType *pv = nullptr;
    ScalarType *pl;
    PetscFunctionBegin;

    ierr = GetRawPointerRead(v, &pv); CHKERRQ(ierr);
    ierr = GetRawPointerWrite(l, &pl); CHKERRQ(ierr);
    
    ierr = this->Laplacian(pl, pv); CHKERRQ(ierr);
    
    ierr = RestoreRawPointerWrite(l, &pl); CHKERRQ(ierr);
    ierr = RestoreRawPointerRead(v, &pv); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute divergence of a vector field
 *******************************************************************/
PetscErrorCode Differentiation::Divergence(ScalarType *l, VecField *v) {
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
PetscErrorCode Differentiation::Divergence(Vec l, VecField *v) {
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
PetscErrorCode Differentiation::Divergence(ScalarType *l, const ScalarType *const *v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
        
    ierr = this->Divergence(l, v[0], v[1], v[2]); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}



}  // end of namespace




#endif  // _DIFFERENTIATION_CPP_
