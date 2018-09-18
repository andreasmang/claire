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
 * @brief view vector entries (transpose output)
 *******************************************************************/
PetscErrorCode VecView(Vec x) {
    PetscErrorCode ierr = 0;
    ScalarType *p_x = NULL;
    IntType nl;
    int rank;
    PetscFunctionBegin;

    ierr = VecGetLocalSize(x, &nl); CHKERRQ(ierr);
    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (rank == 0) {
        std::cout << " VEC VIEW" << std::endl;
        std::cout << " ";
        for (IntType i = 0; i < nl; ++i) {
            std::cout << p_x[i] << " ";
        }
        std::cout << std::endl;
    }

    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


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
    ierr = VecSetFromOptions(x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief rescale data to [0,1]
 *******************************************************************/
PetscErrorCode Normalize(Vec x, IntType nc) {
    PetscErrorCode ierr = 0;
    ScalarType xmin, xmax, xmin_g, xmax_g, *p_x = NULL;
    IntType nl, l;
    int rval;
    std::stringstream ss;

    PetscFunctionBegin;

    if (nc == 1) {
        // get max and min values
        ierr = VecMin(x, NULL, &xmin); CHKERRQ(ierr);
        ierr = VecMax(x, NULL, &xmax); CHKERRQ(ierr);

        if (xmin < 0.0) {
            ss << "negative values in input data detected "
               << xmin << " (setting to zero)";
            ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());

            // compute local size from input vector
            ierr = VecGetLocalSize(x, &nl); CHKERRQ(ierr);

            xmin = 0.0; // resetting
            ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
            for (IntType i = 0; i < nc*nl; ++i) {
                if (p_x[i] < 0.0) p_x[i] = 0.0;
            }
            ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);
        }

        ierr = VecShift(x, -xmin); CHKERRQ(ierr);
        ierr = VecScale(x, 1.0/xmax); CHKERRQ(ierr);
    } else {
        // compute local size from input vector
        ierr = VecGetLocalSize(x, &nl); CHKERRQ(ierr);
        nl /= nc;
        ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {
            xmin = std::numeric_limits<ScalarType>::max();
            xmax = std::numeric_limits<ScalarType>::min();

            // get min and max values
            for (IntType i = 0; i < nl; ++i) {
                l = k*nl + i;
                if (p_x[l] < xmin) {xmin = p_x[l];}
                if (p_x[l] > xmax) {xmax = p_x[l];}
            }

            // get min accross all procs
            rval = MPI_Allreduce(&xmin, &xmin_g, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD);
            ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

            // get max accross all procs
            rval = MPI_Allreduce(&xmax, &xmax_g, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD);
            ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

            if (xmin_g < 0.0) {
                ss << "negative values in input data detected "
                   << xmin << " (setting to zero)";
                ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
                ss.clear(); ss.str(std::string());

                xmin_g = 0.0; // resetting
                for (IntType i = 0; i < nl; ++i) {
                    if (p_x[k*nl + i] < 0.0) p_x[k*nl + i] = 0.0;
                }
            }

            // make sure we do not devide by zero
            xmax_g = (xmax_g != 0.0) ? xmax_g : 1.0;

            // apply shift and scale
            for (IntType i = 0; i < nl; ++i) {
                p_x[k*nl + i] = (p_x[k*nl + i] - xmin_g) / xmax_g;
            }
        }  // for all components
        ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);
    }  // if else

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
