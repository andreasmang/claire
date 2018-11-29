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
#include "DistanceMeasureKernel.hpp"


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
    IntType nt;
    int rval;
    DistanceMeasureKernel::EvaluateFunctionalSL2 kernel;
    ScalarType l2distance, hd;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->m_Domain.nt;
    kernel.nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    hd  = this->m_Opt->GetLebesgueMeasure();

    ierr = this->m_StateVariable->GetArrayRead(kernel.pM, 0, nt); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->GetArrayRead(kernel.pMr); CHKERRQ(ierr);

    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = this->m_Mask->GetArrayRead(kernel.pW); CHKERRQ(ierr);
        
        ierr = kernel.ComputeFunctionalMask(); CHKERRQ(ierr);
    
        ierr = this->m_Mask->RestoreArray(); CHKERRQ(ierr);
    } else {
      ierr = kernel.ComputeFunctional(); CHKERRQ(ierr);
    }
    // all reduce
    rval = MPI_Allreduce(&kernel.value, &l2distance, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    ierr = this->m_ReferenceImage->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);

    // objective value
    *D = 0.5*hd*l2distance/static_cast<ScalarType>(kernel.nc);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for adjoint equation (varies for
 * different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2::SetFinalConditionAE() {
    PetscErrorCode ierr = 0;
    IntType nt;
    DistanceMeasureKernel::FinalConditionSL2 kernel;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    kernel.nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;

    ierr = this->m_ReferenceImage->GetArrayRead(kernel.pMr); CHKERRQ(ierr);
    // index for final condition
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
      ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL, 0, nt); CHKERRQ(ierr);
    } else {
      ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL); CHKERRQ(ierr);
    }
    ierr = this->m_StateVariable->GetArrayRead(kernel.pM, 0, nt); CHKERRQ(ierr);
    
    // compute terminal condition \lambda_1 = -(m_1 - m_R) = m_R - m_1
    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = this->m_Mask->GetArrayRead(kernel.pW); CHKERRQ(ierr);
        
        ierr = kernel.ComputeFinalConditionMaskAE(); CHKERRQ(ierr);
/*#ifdef REG_HAS_CUDA
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
#endif*/
        ierr = this->m_Mask->RestoreArray(); CHKERRQ(ierr);
    } else {
        kernel.ComputeFinalConditionAE();
/*#ifdef REG_HAS_CUDA
        DistanceMeasureSetFinalGPU(&p_l[ll],&p_m[l],p_mr,nc*nl);
#else
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nc*nl; ++i) {
            p_l[ll+i] = p_mr[i] - p_m[l+i];
        }
}  // omp
#endif*/
    }
        
    ierr = this->m_AdjointVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for incremental adjoint equation
 * (varies for different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureSL2::SetFinalConditionIAE() {
    PetscErrorCode ierr = 0;
    IntType nt;
    DistanceMeasureKernel::FinalConditionSL2 kernel;
    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    kernel.nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;

    // index for final condition
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pL, 0, nt); CHKERRQ(ierr);
        ierr = this->m_IncStateVariable->GetArrayRead(kernel.pM, 0, nt); CHKERRQ(ierr);
    } else {
      ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pL); CHKERRQ(ierr);
      ierr = this->m_IncStateVariable->GetArrayRead(kernel.pM); CHKERRQ(ierr);
    }
    
    // compute terminal condition \tilde{\lambda}_1 = -\tilde{m}_1
    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = this->m_Mask->GetArrayRead(kernel.pW); CHKERRQ(ierr);
        
        ierr = kernel.ComputeFinalConditionMaskIAE(); CHKERRQ(ierr);
/*
#pragma omp parallel
{
#pragma omp for
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
                p_ltilde[l+k*nl+i] = -p_w[i]*p_mtilde[l+k*nl+i];
            }
        }
}  // omp*/
        ierr = this->m_Mask->RestoreArray(); CHKERRQ(ierr);
    } else {
        ierr = kernel.ComputeFinalConditionIAE(); CHKERRQ(ierr);
/*
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl*nc; ++i) {
            p_ltilde[l+i] = -p_mtilde[l+i]; // / static_cast<ScalarType>(nc);
        }
}  // omp*/
    }
    ierr = this->m_IncAdjointVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_IncStateVariable->RestoreArray(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



}  // namespace reg




#endif  // _DISTANCEMEASURESL2_CPP_
