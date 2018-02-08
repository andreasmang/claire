/*************************************************************************
 *  Copyright (c) 2017.
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

#ifndef _L2DISTANCE_CPP_
#define _L2DISTANCE_CPP_

#include "DistanceMeasure.hpp"
#include "L2Distance.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
L2Distance::L2Distance() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
L2Distance::~L2Distance() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
L2Distance::L2Distance(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode L2Distance::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluate the functional (i.e., the distance measure)
 * D = (1/2)*||m1 - mR||_L2
 *******************************************************************/
PetscErrorCode L2Distance::EvaluateFunctional(ScalarType* D) {
    PetscErrorCode ierr = 0;
    ScalarType *p_mr = NULL, *p_m = NULL, *p_w = NULL;
    IntType nt, nc, nl, l;
    int rval;
    ScalarType dr, value, l2distance;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);

    l = nt*nl*nc;
    value = 0.0;
    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = VecGetArray(this->m_Mask, &p_w); CHKERRQ(ierr);
        for (IntType i = 0; i < nc*nl; ++i) {
            dr = (p_mr[i] - p_m[l+i]);
            value += p_w[i]*dr*dr;
        }
        ierr = VecRestoreArray(this->m_Mask, &p_w); CHKERRQ(ierr);
    } else {
        for (IntType i = 0; i < nc*nl; ++i) {
            dr = (p_mr[i] - p_m[l+i]);
            value += dr*dr;
        }
    }
    // all reduce
    rval = MPI_Allreduce(&value, &l2distance, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    // objective value
    *D = 0.5*l2distance/static_cast<ScalarType>(nc);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for adjoint equaiton (varies for
 * different distance measres)
 *******************************************************************/
PetscErrorCode L2Distance::SetFinalCondition(Vec lambda) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}


}  // namespace reg







#endif  // _L2DISTANCE_CPP_
