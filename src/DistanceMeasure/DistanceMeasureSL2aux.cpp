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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _DISTANCEMEASURESL2AUX_CPP_
#define _DISTANCEMEASURESL2AUX_CPP_



#include "DistanceMeasureSL2aux.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DistanceMeasureSL2aux::DistanceMeasureSL2aux() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DistanceMeasureSL2aux::~DistanceMeasureSL2aux() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DistanceMeasureSL2aux::DistanceMeasureSL2aux(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2aux::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluate the functional (i.e., the distance measure)
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2aux::EvaluateFunctional(ScalarType* D) {
    PetscErrorCode ierr = 0;
    ScalarType *p_mr = NULL, *p_m = NULL, *p_q = NULL, *p_c = NULL;
    IntType nt, nc, nl, l;
    int rval;
    ScalarType dr, value, val1, val2, l2distance, hx;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_AuxVar1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AuxVar2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    hx  = this->m_Opt->GetLebesgueMeasure();   

    ierr = VecGetArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(*this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);

    l = nt*nl*nc;
    value = 0.0, val1 = 0.0, val2 = 0.0;
    ierr = VecGetArray(*this->m_AuxVar1, &p_c); CHKERRQ(ierr);
    ierr = VecGetArray(*this->m_AuxVar2, &p_q); CHKERRQ(ierr);
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        for (IntType i = 0; i < nl; ++i) {
            // mismatch: mr - m1*(1-c1)
            dr    = (p_mr[k*nl+i] - p_m[l+k*nl+i]*(1.0 - p_c[i]));
            val1 += dr*dr;
            // q^k*m^k
            val2 += p_q[k*nl+i]*p_m[l+k*nl+i];
        }
    }
    ierr = VecRestoreArray(*this->m_AuxVar2, &p_q); CHKERRQ(ierr);
    ierr = VecRestoreArray(*this->m_AuxVar1, &p_c); CHKERRQ(ierr);

    // all reduce
    rval = MPI_Allreduce(&val1, &value, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    l2distance  = value;

    rval = MPI_Allreduce(&val2, &value, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    l2distance += value;
    // parse value to registration monitor for display
    this->m_Opt->m_Monitor.qmval = value;

    ierr = VecRestoreArray(*this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = VecRestoreArray(*this->m_StateVariable, &p_m); CHKERRQ(ierr);

    // objective value
    *D = 0.5*hx*l2distance/static_cast<ScalarType>(nc);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for adjoint equaiton (varies for
 * different distance measres)
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2aux::SetFinalConditionAE() {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt, l, ll;
    ScalarType *p_mr = NULL, *p_m = NULL, *p_l = NULL, *p_c = NULL, scale;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_AuxVar1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AuxVar2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    // index for final condition
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ll = nt*nc*nl;
    } else {
        ll = 0;
    }

    ierr = GetRawPointer(this->m_AuxVar1, &p_c); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);

    l = nt*nc*nl;
#pragma omp parallel
{
#pragma omp for
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            scale = (1.0 - p_c[i]);
            // compute initial condition
            p_l[ll+k*nl+i] = (p_mr[k*nl+i] - p_m[l+k*nl+i]*scale)*scale;
        }
    }
}  // omp

    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_AuxVar1, &p_c); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for incremental adjoint equation
 * (varies for different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2aux::SetFinalConditionIAE() {
    PetscErrorCode ierr = 0;
    IntType nt, nc, nl, ll;
    ScalarType *p_mtilde = NULL, *p_ltilde = NULL, *p_c = NULL;
    ScalarType scale = 0.0;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_AuxVar1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    // index for final condition
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ll = nt*nc*nl;
    } else {
        ll = 0;
    }

    ierr = GetRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AuxVar1, &p_c); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
            scale = (1.0 - p_c[i]);
            scale *= scale;
            p_ltilde[ll+k*nl+i] = -scale*p_mtilde[ll+k*nl+i]; // / static_cast<ScalarType>(nc);
        }
    }
}  // omp
    ierr = RestoreRawPointer(this->m_AuxVar1, &p_c); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _DISTANCEMEASURESL2AUX_CPP_
