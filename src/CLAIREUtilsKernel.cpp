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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _CLAIREUTILSKERNEL_CPP_
#define _CLAIREUTILSKERNEL_CPP_

#include "CLAIREUtils.hpp"




namespace reg {


/********************************************************************
 * @brief interface to create vector
 *******************************************************************/
PetscErrorCode VecCreate(Vec& x, IntType nl, IntType ng) {
    PetscErrorCode ierr = 0;

    if (x != NULL) {
        ierr = VecDestroy(&x); CHKERRQ(ierr);
        x = NULL;
    }

    ierr = VecCreate(PETSC_COMM_WORLD, &x); CHKERRQ(ierr);
    ierr = VecSetSizes(x, nl, ng); CHKERRQ(ierr);
#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)
    ierr = VecSetType(x, VECCUDA); CHKERRQ(ierr);
#else
    ierr = VecSetType(x, VECSTANDARD); CHKERRQ(ierr);
#endif
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief clip image to be in [0,1]
 *******************************************************************/
PetscErrorCode Clip(Vec x, IntType nc) {
    PetscErrorCode ierr = 0;
    ScalarType *p_x = NULL;
    IntType nl;

    PetscFunctionBegin;

    ierr = VecGetLocalSize(x, &nl); CHKERRQ(ierr);
    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
    for (IntType i = 0; i < nl; ++i) {
        if (p_x[i] < 0.0) p_x[i] = 0.0;
        if (p_x[i] > 1.0) p_x[i] = 1.0;
    }  // for all components
    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief ensure that sum accross all image components is in [0,1]
 * for a particular location x
 *******************************************************************/
PetscErrorCode EnsurePartitionOfUnity(Vec x, IntType nc) {
    PetscErrorCode ierr = 0;
    ScalarType *p_x = NULL, sum, background;
    IntType nl, l;

    PetscFunctionBegin;

    ierr = VecGetLocalSize(x, &nl); CHKERRQ(ierr);
    nl /= nc;

    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
    for (IntType i = 0; i < nl; ++i) {
        sum = 0.0;
        for (IntType k = 0; k < nc; ++k) {
            l = k*nl + i;

            if (p_x[l] < 0.0) p_x[l] = 0.0;
            if (p_x[l] > 1.0) p_x[l] = 1.0;

            sum += p_x[l];
        }
        background = 1.0 - sum;
        if (background <= 0.0) background = 0.0;
        if (background >= 1.0) background = 1.0;

        sum += background;

        for (IntType k = 0; k < nc; ++k) {
            p_x[k*nl + i] /= sum;
        }

    }  // for all components

    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}
  
}

#endif
