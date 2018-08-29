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

#ifndef _DISTANCEMEASURESL2_CPP_
#define _DISTANCEMEASURESL2_CPP_

#include "DistanceMeasureSL2.hpp"

#ifdef REG_HAS_CUDA
#include "distance_kernel.hpp"
#endif



namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DistanceMeasureSL2::DistanceMeasureSL2() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DistanceMeasureSL2::~DistanceMeasureSL2() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DistanceMeasureSL2::DistanceMeasureSL2(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief set up scale
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2::SetupScale() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief evaluate the functional (i.e., the distance measure)
 * D = (1/2)*||m1 - mR||_L2
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2::EvaluateFunctional(ScalarType* D) {
    PetscErrorCode ierr = 0;
    ScalarType *p_mr = NULL, *p_m = NULL, *p_w = NULL;
    IntType nt, nc, nl, l;
    int rval;
    ScalarType dr, value, l2distance, hx;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    hx  = this->m_Opt->GetLebesgueMeasure();   

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);

    l = nt*nl*nc;
    value = 0.0;
    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = GetRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
                dr = (p_mr[k*nl+i] - p_m[l+k*nl+i]);
                value += p_w[i]*dr*dr;
            }
        }
        ierr = RestoreRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
        for (IntType i = 0; i < nc*nl; ++i) {
            dr = (p_mr[i] - p_m[l+i]);
            value += dr*dr;
        }
    }
    // all reduce
    rval = MPI_Allreduce(&value, &l2distance, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    ierr = RestoreRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    // objective value
    *D = 0.5*hx*l2distance/static_cast<ScalarType>(nc);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for adjoint equation (varies for
 * different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2::SetFinalConditionAE() {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt, l, ll;
    //ScalarType *p_mr = NULL, *p_m = NULL, *p_l = NULL, *p_w = NULL;
    ScalarType *p_l = NULL;
    const ScalarType *p_mr = NULL, *p_m = NULL, *p_w = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    // index for final condition
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ll = nt*nc*nl;
    } else {
        ll = 0;
    }

    ierr = GetRawPointerRead(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointerRead(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    
    l = nt*nc*nl;
    // compute terminal condition \lambda_1 = -(m_1 - m_R) = m_R - m_1
    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = GetRawPointerRead(this->m_Mask, &p_w); CHKERRQ(ierr);
#ifdef REG_HAS_CUDA
        DistanceMeasureSetFinalMaskGPU(&p_l[ll],&p_m[l],p_mr,p_w,nl,nc);
#else
#pragma omp parallel
{
#pragma omp for
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
                p_l[ll+k*nl+i] = p_w[i]*(p_mr[k*nl+i] - p_m[l+k*nl+i]);
            }
        }
}  // omp
#endif
        ierr = RestoreRawPointerRead(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
#ifdef REG_HAS_CUDA
        DistanceMeasureSetFinalGPU(&p_l[ll],&p_m[l],p_mr,nc*nl);
#else
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nc*nl; ++i) {
            p_l[ll+i] = p_mr[i] - p_m[l+i];
        }
}  // omp
#endif
    }
        
    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = RestoreRawPointerRead(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = RestoreRawPointerRead(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for incremental adjoint equation
 * (varies for different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2::SetFinalConditionIAE() {
    PetscErrorCode ierr = 0;
    IntType nt, nc, nl, l;
    ScalarType *p_mtilde = NULL, *p_ltilde = NULL, *p_w = NULL;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    // index for final condition
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        l = nt*nc*nl;
    } else {
        l = 0;
    }
    ierr = GetRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    // compute terminal condition \tilde{\lambda}_1 = -\tilde{m}_1
    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = GetRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
                p_ltilde[l+k*nl+i] = -p_w[i]*p_mtilde[l+k*nl+i];
            }
        }
}  // omp
        ierr = RestoreRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl*nc; ++i) {
            p_ltilde[l+i] = -p_mtilde[l+i]; // / static_cast<ScalarType>(nc);
        }
}  // omp
    }
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _DISTANCEMEASURESL2_CPP_
