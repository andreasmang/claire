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

#ifndef _DISTANCEMEASURENCC_CPP_
#define _DISTANCEMEASURENCC_CPP_

#include "DistanceMeasureNCC.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DistanceMeasureNCC::DistanceMeasureNCC() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DistanceMeasureNCC::~DistanceMeasureNCC() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DistanceMeasureNCC::DistanceMeasureNCC(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DistanceMeasureNCC::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluate the functional (i.e., the distance measure)
 * D = 1.0 - <m1,mR>_L2/(||m1||_L2 * ||mR||_L2)
 *******************************************************************/
PetscErrorCode DistanceMeasureNCC::EvaluateFunctional(ScalarType* D) {
    PetscErrorCode ierr = 0;
    ScalarType *p_mr = NULL, *p_m = NULL, *p_w = NULL;
    IntType nt, nc, nl, l;
    int rval;
    ScalarType norm_m1_loc, norm_mR_loc, inpr_m1_mR_loc, norm_m1, norm_mR, inpr_m1_mR;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // Get sizes
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);

    l = nt*nl*nc;
    norm_m1_loc = 0.0;
    norm_mR_loc = 0.0;
    inpr_m1_mR_loc = 0.0;
    if (this->m_Mask != NULL) {
        // Mask objective functional
        ierr = GetRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
                norm_m1_loc += p_w[i]*p_m[l+k*nl+i]*p_m[l+k*nl+i];
		norm_mR_loc += p_w[i]*p_mr[k*nl+i]*p_mr[k*nl+i];
		inpr_m1_mR_loc += p_w[i]*p_m[l+k*nl+i]*p_mr[k*nl+i];
            }
        }
        ierr = RestoreRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
        for (IntType i = 0; i < nc*nl; ++i) {
	    norm_m1_loc += p_m[l+i]*p_m[l+i];
            norm_mR_loc += p_mr[i]*p_mr[i];
            inpr_m1_mR_loc += p_m[l+i]*p_mr[i];
        }
    }
    // All reduce various pieces
    rval = MPI_Allreduce(&norm_m1_loc, &norm_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&norm_mR_loc, &norm_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&inpr_m1_mR_loc, &inpr_m1_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    ierr = RestoreRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    // Objective value
    *D = 0.5 - 0.5*(inpr_m1_mR*inpr_m1_mR)/(norm_m1*norm_mR);
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for adjoint equation 
 * (varies for different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureNCC::SetFinalConditionAE() {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt, l, ll;
    int rval;
    ScalarType *p_mr = NULL, *p_m = NULL, *p_l = NULL, *p_w = NULL;
    ScalarType norm_m1_loc, norm_mR_loc, inpr_m1_mR_loc, norm_m1, norm_mR, inpr_m1_mR; 
    ScalarType const1, const2, hx;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    hx  = this->m_Opt->GetLebesgueMeasure();   

    // Index for final condition
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ll = nt*nc*nl;
    } else {
        ll = 0;
    }

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);

    l = nt*nc*nl;
    norm_m1_loc = 0.0;
    norm_mR_loc = 0.0;
    inpr_m1_mR_loc = 0.0;
    
    /* compute terminal condition  
	lambda = (<m1,mR>/<m1,m1><mR,mR>) * (mR - (<m1,mR>/<m1,m1>)*m1)
    */

    // First, calculate necessary inner products locally
    if (this->m_Mask != NULL) {
        // Mask objective functional
        ierr = GetRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
                norm_m1_loc += p_w[i]*p_m[l+k*nl+i]*p_m[l+k*nl+i];
		norm_mR_loc += p_w[i]*p_mr[k*nl+i]*p_mr[k*nl+i];
		inpr_m1_mR_loc += p_w[i]*p_m[l+k*nl+i]*p_mr[k*nl+i];
            }
        }
        ierr = RestoreRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
        for (IntType i = 0; i < nc*nl; ++i) {
	    norm_m1_loc += p_m[l+i]*p_m[l+i];
            norm_mR_loc += p_mr[i]*p_mr[i];
            inpr_m1_mR_loc += p_m[l+i]*p_mr[i];
        }
    }
    
    // All reduce for full inner products
    rval = MPI_Allreduce(&norm_m1_loc, &norm_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&norm_mR_loc, &norm_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&inpr_m1_mR_loc, &inpr_m1_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    // Now, write the terminal condition to lambda
    const1 = inpr_m1_mR/(hx*norm_m1*norm_mR); 
    const2 = (inpr_m1_mR*inpr_m1_mR)/(hx*norm_m1*norm_m1*norm_mR);

   if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = GetRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            for (IntType i = 0; i < nl; ++i) {  // for all grid nodes
                p_l[ll+k*nl+i] = p_w[i]*(const1*p_mr[k*nl+i] - const2*p_m[l+k*nl+i]);
            }
        }
}  // omp
        ierr = RestoreRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nc*nl; ++i) {
	    p_l[ll+i] = const1*p_mr[i] - const2*p_m[l+i];
        }
}  // omp
    }
    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for incremental adjoint equation
 * (varies for different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureNCC::SetFinalConditionIAE() {
    PetscErrorCode ierr = 0;
    IntType nt, nc, nl, ll, l;
    int rval;
    ScalarType *p_m = NULL, *p_mr = NULL, *p_mtilde = NULL, *p_ltilde = NULL, *p_w = NULL;
    ScalarType const1, const2, const3, const4, const5, hx; 
    ScalarType norm_m1_loc, norm_mR_loc, inpr_m1_mR_loc, inpr_m1_mtilde_loc, inpr_mR_mtilde_loc;
    ScalarType norm_m1, norm_mR, inpr_m1_mR, inpr_m1_mtilde, inpr_mR_mtilde;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    hx  = this->m_Opt->GetLebesgueMeasure();   

    // Index for final condition
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ll = nt*nc*nl;
    } else {
        ll = 0;
    }
    ierr = GetRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    
    // Compute terminal condition 

    // First, calculate necessary inner products locally
    norm_m1_loc = 0.0;
    norm_mR_loc = 0.0;
    inpr_m1_mR_loc = 0.0;
    inpr_m1_mtilde_loc = 0.0;
    inpr_mR_mtilde_loc = 0.0;
    l = nt*nc*nl;

    if (this->m_Mask != NULL) {
        // Mask objective functional
        ierr = GetRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            for (IntType i = 0; i < nl; ++i) {
	        norm_m1_loc += p_w[i]*p_m[l+k*nl+i]*p_m[l+k*nl+i];
                norm_mR_loc += p_w[i]* p_mr[k*nl+i]*p_mr[k*nl+i];
                inpr_m1_mR_loc += p_w[i]*p_m[l+k*nl+i]*p_mr[k*nl+i];
                inpr_m1_mtilde_loc += p_w[i]*p_m[l+k*nl+i]*p_mtilde[l+k*nl+i];
                inpr_mR_mtilde_loc += p_w[i]*p_mtilde[l+k*nl+i]*p_mr[k*nl+i];
            }
        }
        ierr = RestoreRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
        for (IntType i = 0; i < nc*nl; ++i) {
	    norm_m1_loc += p_m[l+i]*p_m[l+i];
            norm_mR_loc += p_mr[i]*p_mr[i];
            inpr_m1_mR_loc += p_m[l+i]*p_mr[i];
            inpr_m1_mtilde_loc += p_m[l+i]*p_mtilde[l+i];
            inpr_mR_mtilde_loc += p_mr[i]*p_mtilde[l+i];
        }
    }
    
    // All reduce for full inner products
    rval = MPI_Allreduce(&norm_m1_loc, &norm_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&norm_mR_loc, &norm_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&inpr_m1_mR_loc, &inpr_m1_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&inpr_m1_mtilde_loc, &inpr_m1_mtilde, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&inpr_mR_mtilde_loc, &inpr_mR_mtilde, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    // Now, write the terminal condition to lambda tilde
    const1 = inpr_mR_mtilde/(hx*norm_m1*norm_mR); 
    const2 = 2.0*(inpr_m1_mR*inpr_m1_mtilde)/(hx*norm_m1*norm_m1*norm_mR);
    const3 = 4.0*(inpr_m1_mR*inpr_m1_mR*inpr_m1_mtilde)/(hx*norm_m1*norm_m1*norm_m1*norm_mR);
    const4 = 2.0*(inpr_m1_mR*inpr_mR_mtilde)/(hx*norm_m1*norm_m1*norm_mR);
    const5 = (inpr_m1_mR*inpr_m1_mR)/(hx*norm_m1*norm_m1*norm_mR);
 
    if (this->m_Mask != NULL) {
        // Mask objective functional
        ierr = GetRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
        for (IntType k = 0; k < nc; ++k) {  // For all image components
            for (IntType i = 0; i < nl; ++i) {  // For all grid nodes
                p_ltilde[ll+k*nl+i] = p_w[i]*(const1*p_mr[k*nl+i] - const2*p_mr[k*nl+i] + const3*p_m[l+k*nl+i] - const4*p_m[l+k*nl+i] - const5*p_mtilde[l+k*nl+i]);
            }
        }
}  // omp
        ierr = RestoreRawPointer(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl*nc; ++i) {
            p_ltilde[ll+i] = const1*p_mr[i] - const2*p_mr[i] + const3*p_m[l+i] - const4*p_m[l+i] - const5*p_mtilde[l+i];
        }
}  // omp
    }
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

}  // namespace reg




#endif  // _DISTANCEMEASURENCC_CPP_
