/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGUTILS_CPP_
#define _REGUTILS_CPP_


#include "RegUtils.hpp"
#include <time.h>
#include <limits>

#ifdef REG_HAS_PNETCDF
#include "pnetcdf.h"
#endif




namespace reg {




/********************************************************************
 * @brief error handling: check if condition is valid, and if
 * not throw an error PETSc style
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Assert"
PetscErrorCode Assert(bool condition,std::string msg) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if(condition == false) {
        ierr = ThrowError(msg); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get the filename of an image
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetFileName"
PetscErrorCode GetFileName(std::string& filename, std::string file) {
    PetscErrorCode ierr = 0;
    std::string path;
    size_t sep;
    PetscFunctionBegin;

    sep = file.find_last_of("\\/");

    if (sep != std::string::npos) {
        path=file.substr(0,sep);
        filename=file.substr(sep + 1);
    }

    if (filename == "") { filename = file; }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get the filename, path, and extension
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetFileName"
PetscErrorCode GetFileName(std::string& path, std::string& filename,
                           std::string& extension, std::string file) {
    PetscErrorCode ierr = 0;
    std::string::size_type idx;

    PetscFunctionBegin;

    // get path
    idx = file.find_last_of("\\/");
    if (idx != std::string::npos) {
        path = file.substr(0,idx);
        filename = file.substr(idx + 1);
    }
    if (filename == "") {
        filename = file;
    }

    // get extension
    idx = filename.rfind(".");
    if (idx != std::string::npos) {
        extension = filename.substr(idx+1);

        // handle zipped files
        if (strcmp(extension.c_str(),"gz") == 0) {
            filename = filename.substr(0,idx);
            idx = filename.rfind(".");
            if(idx != std::string::npos) {
                extension = filename.substr(idx+1);
                extension = extension + ".gz";
            }
        }
        extension = "." + extension;
        filename  = filename.substr(0,idx);

    } else {
        ierr = ThrowError("no extension found"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief check if file exists
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "FileExists"
bool FileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}




/********************************************************************
 * @brief print msg (interfaces petsc)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Msg"
PetscErrorCode Msg(std::string msg) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;

    PetscFunctionBegin;

    ss << std::left << msg;
    msg = " "  + ss.str() + "\n";

    // display message
    ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief print msg (interfaces petsc)
 * Author: Andreas Mang
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DbgMsg"
PetscErrorCode DbgMsg(std::string msg) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;

    PetscFunctionBegin;

    ss << std::left << std::setw(98) << msg;
    msg = "\x001b[90m[ "  + ss.str() + "]\x1b[0m\n";

    // display message
    ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief print warning msg (interfaces petsc)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WrngMsg"
PetscErrorCode WrngMsg(std::string msg) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;

    PetscFunctionBegin;

    ss << std::left << std::setw(98) << msg;
    msg = "\x1b[33m[ " + ss.str() + "]\x1b[0m\n";

    // display error
    ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief throw error
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ThrowError"
PetscErrorCode ThrowError(std::string msg) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    std::string errmsg = "\x1b[31mERROR: " + msg + "\x1b[0m";
    ierr = PetscError(PETSC_COMM_WORLD,__LINE__,
                    PETSC_FUNCTION_NAME,__FILE__,
                    1,PETSC_ERROR_INITIAL,
                    errmsg.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief mpi error handling
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MPIERRQ"
PetscErrorCode MPIERRQ(int cerr) {
    int rank;
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if (cerr != MPI_SUCCESS) {
        char error_string[BUFSIZ];
        int length_of_error_string, error_class;

        MPI_Error_class(cerr, &error_class);
        MPI_Error_string(error_class, error_string, &length_of_error_string);
        ierr = ThrowError(error_string); CHKERRQ(ierr);

    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief mpi error handling
 *******************************************************************/
#ifdef REG_HAS_PNETCDF
#undef __FUNCT__
#define __FUNCT__ "NCERRQ"
PetscErrorCode NCERRQ(int cerr) {
    int rank;
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (cerr != NC_NOERR) {
        ss << ncmpi_strerror(cerr);
        ierr = ThrowError(ss.str()); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}
#endif




/********************************************************************
 * @brief function to slow down code
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "isleep"
void isleep( unsigned int nanosec ) {
    clock_t wait = (clock_t) nanosec;
    clock_t start_time = clock();
    while( clock() != start_time + wait ) {};
    return;
}




/********************************************************************
 * @brief setup library
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "InitializeDataDistribution"
PetscErrorCode InitializeDataDistribution(int nthreads,
                                          int *c_grid,
                                          MPI_Comm& c_comm) {
    PetscErrorCode ierr = 0;
    int nprocs, ompthreads, np;
    std::stringstream ss;

    PetscFunctionBegin;

    omp_set_dynamic(0);
    omp_set_num_threads(nthreads);

    // check if number of threads is consistent with user options
    ompthreads=omp_get_max_threads();
    ss << "max number of openmp threads is not a match (user,set)=("
       << nthreads <<"," << ompthreads <<")\n";
    ierr = Assert(ompthreads == nthreads,ss.str().c_str()); CHKERRQ(ierr);
    ss.str( std::string() );
    ss.clear();

    // set up MPI/cartesian grid
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    np = c_grid[0]*c_grid[1];

    // check number of procs
    if (np != nprocs) {
        // update cartesian grid layout
        c_grid[0]=0;
        c_grid[1]=0;
        MPI_Dims_create(nprocs, 2, c_grid);
    }

    if (c_comm != NULL) {MPI_Comm_free(&c_comm); c_comm = NULL;}

    // initialize accft
    accfft_create_comm(PETSC_COMM_WORLD, c_grid, &c_comm);
    accfft_init(nthreads);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief view vector entries (transpose output)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Finalize"
PetscErrorCode Finalize() {
    PetscErrorCode ierr = 0;

    accfft_cleanup();

    // clean up petsc
    ierr = PetscFinalize(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief view vector entries (transpose output)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecView"
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
#undef __FUNCT__
#define __FUNCT__ "VecCreate"
PetscErrorCode VecCreate(Vec& x, IntType nl, IntType ng) {
    PetscErrorCode ierr = 0;

    if(x != NULL) {ierr = VecDestroy(&x); CHKERRQ(ierr); x = NULL;}

    ierr = VecCreate(PETSC_COMM_WORLD, &x); CHKERRQ(ierr);
    ierr = VecSetSizes(x, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief rescale data to [xminout,xmaxout]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Rescale"
PetscErrorCode Rescale(Vec x, ScalarType xminout, ScalarType xmaxout, IntType nc) {
    PetscErrorCode ierr = 0;
    ScalarType xmin, xmax, xmin_g, xmax_g, xshift, xscale, *p_x = NULL;
    IntType nl, l;
    int rval;
    std::stringstream ss;

    PetscFunctionBegin;

    if (nc == 1) {
        // get max and min values
        ierr = VecMin(x, NULL, &xmin); CHKERRQ(ierr);
        ierr = VecMax(x, NULL, &xmax); CHKERRQ(ierr);

        xshift = xminout - xmin;
        ierr = VecShift(x, xshift); CHKERRQ(ierr);

        xmax = (xmax != 0.0) ? xmax : 1.0;
        xscale = (xmaxout == 0.0) ? 1.0 : xmaxout / xmax;

        ierr = VecScale(x, xscale); CHKERRQ(ierr);
    } else {
        // compute local size from input vector
        ierr = VecGetSize(x, &nl); CHKERRQ(ierr);
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
            rval = MPI_Allreduce(&xmin, &xmin_g, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
            ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

            // get max accross all procs
            rval = MPI_Allreduce(&xmax, &xmax_g, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
            ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

            // compute shift and scale
            xshift = xminout - xmin_g;
            xmax_g = (xmax_g != 0.0) ? xmax_g : 1.0;
            xscale = (xmaxout == 0.0) ? 1.0 : xmaxout / xmax_g;

            // apply shift and scale
            for (IntType i = 0; i < nl; ++i) {
                p_x[k*nl + i] = xscale*(p_x[k*nl + i] + xshift);
            }
        }  // for all components
        ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);
    }  // if else
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief parse string of NUMxNUMxNUM into a vector
 *******************************************************************/
std::vector<unsigned int> String2Vec(const std::string & str) {
    std::vector<unsigned int> vect;
    int ival;
    std::string::size_type xpos = str.find('x',0);

    if (xpos == std::string::npos) {
        // only one uint
        vect.push_back( static_cast<unsigned int>( atoi(str.c_str()) ));
        return vect;
    }

    // first uint$
    ival = atoi((str.substr(0, xpos)).c_str());
    vect.push_back(static_cast<unsigned int>(ival));

    while (true) {
        std::string::size_type newxpos = xpos;
        xpos = str.find('x', newxpos+1);

        if (xpos == std::string::npos) {
            ival = atoi((str.substr(newxpos+1, str.length()-newxpos-1)).c_str());
            vect.push_back(static_cast<unsigned int>(ival));
            return vect;
        }
        ival = atoi((str.substr(newxpos+1, xpos-newxpos-1)).c_str() );
        vect.push_back(static_cast<unsigned int>(ival));
    }
}




}  //  namespace reg




#endif   // _REGUTILS_CPP_
