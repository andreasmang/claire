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




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
L2Distance::L2Distance() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
L2Distance::~L2Distance(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
L2Distance::L2Distance(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
L2Distance::EvaluateFunctional(ScalarType* D, Vec m) {
    PetscErrorCode ierr = 0;
    ScalarType *p_mr = NULL, *p_m = NULL;
    IntType nt, nc, nl, l;
    int rval;
    ScalarType dr, value, l2distance;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;

    ierr = VecGetArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);

    l = nt*nl*nc;
    value = 0.0;
#pragma omp parallel for private(dr) reduction(+:value)
    for (IntType i = 0; i < nc*nl; ++i) {
        dr = (p_mr[i] - p_m[l+i]);
        value += dr*dr;
    }

    ierr = VecRestoreArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);

    rval = MPI_Allreduce(&value, &l2distance, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    // objective value
    *D = 0.5*l2distance/static_cast<ScalarType>(nc);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


}  // namespace reg




#endif  // _L2DISTANCE_CPP_
