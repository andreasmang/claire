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

#ifndef _CLAIREUTILS_CPP_
#define _CLAIREUTILS_CPP_

#include "CLAIREUtils.hpp"

//#include "cuda_helper.hpp"

namespace reg {
  
const char * basename(const char* path) {
  const char *pos = path;
  const char *ptr = path;
  while (*ptr) {
    if (*ptr == '/') pos = ptr+1;
    ptr++;
  }
  return pos;
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
 * @brief compute pointwise norm of vector field
 *******************************************************************/
PetscErrorCode VecFieldPointWiseNorm(Vec norm, Vec m_X1, Vec m_X2, Vec m_X3) {
    PetscErrorCode ierr = 0;
    ScalarType *p_m = NULL;
    const ScalarType *p_X1 = NULL, *p_X2 = NULL, *p_X3 = NULL;
    IntType nl;
    
    PetscFunctionBegin;
    
    ierr = GetRawPointer(norm, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointerRead(m_X1, &p_X1); CHKERRQ(ierr);
    ierr = GetRawPointerRead(m_X2, &p_X2); CHKERRQ(ierr);
    ierr = GetRawPointerRead(m_X3, &p_X3); CHKERRQ(ierr);
    
    ierr = VecGetLocalSize(norm, &nl); CHKERRQ(ierr);
#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)
    ierr = reg::VecFieldPointWiseNormGPU(p_m, p_X1, p_X2, p_X3, nl); CHKERRQ(ierr);    
#else
#pragma omp parallel
{
#pragma omp for
  for (IntType i = 0; i < nl; ++i) {
    p_m[i] = sqrtf(p_X1[i]*p_X1[i] + p_X2[i]*p_X2[i] + p_X3[i]*p_X3[i]);
  }
}
#endif

    ierr = RestoreRawPointer(norm, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointerRead(m_X1, &p_X1); CHKERRQ(ierr);
    ierr = RestoreRawPointerRead(m_X2, &p_X2); CHKERRQ(ierr);
    ierr = RestoreRawPointerRead(m_X3, &p_X3); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);

}
  
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
 * @brief error handling: check if condition is valid, and if
 * not throw an error PETSc style
 *******************************************************************/
PetscErrorCode Assert(bool condition, std::string msg) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (condition == false) {
        ierr = ThrowError(msg); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get the filename of an image
 ********************************************************************/
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
bool FileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}




/********************************************************************
 * @brief show extremal values and norm of vector
 *******************************************************************/
PetscErrorCode ShowValues(Vec x, IntType nc) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    IntType nl;
    ScalarType value, maxval, minval, minval_g, maxval_g;
    ScalarType *p_x = NULL;
    int rval;
    PetscFunctionBegin;

    if (nc == 1) {
        ierr = VecNorm(x, NORM_2, &value); CHKERRQ(ierr);
        ierr = VecMax(x, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(x, NULL, &minval); CHKERRQ(ierr);

        ss << "(norm,min,max) = (" << std::scientific << value
           << "," << minval << "," << maxval << ")";

        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    } else {
        ierr = VecGetLocalSize(x, &nl); CHKERRQ(ierr);
        nl /= nc;
        ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {
            minval = std::numeric_limits<ScalarType>::max();
            maxval = std::numeric_limits<ScalarType>::min();

            // get min and max values
            for (IntType i = 0; i < nl; ++i) {
                value = p_x[k*nl + i];
                if (value < minval) {minval = value;}
                if (value > maxval) {maxval = value;}
            }

            // get min accross all procs
            rval = MPI_Allreduce(&minval, &minval_g, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD);
            ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

            // get max accross all procs
            rval = MPI_Allreduce(&maxval, &maxval_g, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD);
            ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

            // TODO: add norm
            ss  << "component " << k << " (min,max) = ("
                << std::scientific <<  minval_g << "," << maxval_g << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
        ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief print msg (interfaces petsc)
 *******************************************************************/
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
/*PetscErrorCode DbgMsgCall(std::string msg) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;

    PetscFunctionBegin;

    ss << std::left << std::setw(98) << msg;
    msg = "\x001b[90m[ "  + ss.str() + "]\x1b[0m\n";

    // display message
    ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
    
    PetscFunctionReturn(ierr);
}*/
PetscErrorCode DbgMsgCall(std::string msg, int line, const char *file, int level) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    std::stringstream ss2;

    PetscFunctionBegin;
    
    std::string color = "\x001b[90m"; // dark grey
    switch (level) {
    case 0:
      if (file)
        color = "\x1b[35m"; // magenta
      else
        color = "\x1b[90m"; // darkgray
      break;
    case 1:
      if (file)
        color = "\x1b[34m"; // blue
      else
        color = "\x1b[90m"; // darkgray
      break;
    case 2:
      if (file)
        color = "\x1b[36m"; // cyan
      else
        color = "\x1b[90m"; // darkgray
      break;
    case 3:
      if (file)
        color = "\x1b[32m"; // green
      else
        color = "\x1b[90m"; // darkgray
      break;
    };

    if (file)
      ss2 << basename(file) << ":" << line;
    ss << std::setw(98-ss2.str().size()) << std::left << msg << std::right << ss2.str();
    //ss << std::left << msg;
    msg = color + "[ "  + ss.str() + "]\x1b[0m\n";

    // display message
    ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
    
    //cudaPrintDeviceMemory();
    
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief print warning msg (interfaces petsc)
 *******************************************************************/
PetscErrorCode WrngMsgCall(std::string msg, int line, const char* file) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    std::stringstream ss2;

    PetscFunctionBegin;

    ss2 << basename(file) << ":" << line;
    ss << std::setw(98-ss2.str().size()) << std::left << msg << std::right << ss2.str();
    msg = "\x1b[33m[ " + ss.str() + "]\x1b[0m\n";

    // display error
    ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief throw error
 *******************************************************************/
PetscErrorCode ThrowErrorMsg(std::bad_alloc& err, int line, const char *file) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;

    PetscFunctionBegin;

    ss << "allocation error " << err.what();
    ierr = ThrowErrorMsg(ss.str(), line, file); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief throw error
 *******************************************************************/
PetscErrorCode ThrowErrorMsg(std::exception& err, int line, const char *file) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    PetscFunctionBegin;

    ss << "exception caught: " << err.what();
    ierr = ThrowErrorMsg(ss.str(), line, file); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief throw error
 *******************************************************************/
PetscErrorCode ThrowErrorMsg(std::string msg, int line, const char *file) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;
    std::stringstream ss2;

    PetscFunctionBegin;

    ss2 << basename(file) << ":" << line;
    ss << std::setw(98-ss2.str().size()) << std::left << msg << std::right << ss2.str();
    std::string errmsg = "\x1b[31mERROR: " + ss.str() + "\x1b[0m";
    ierr = PetscError(PETSC_COMM_WORLD, __LINE__, PETSC_FUNCTION_NAME, __FILE__, 1, PETSC_ERROR_INITIAL, errmsg.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief mpi error handling
 *******************************************************************/
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
void isleep(unsigned int nanosec) {
    clock_t wait = (clock_t) nanosec;
    clock_t start_time = clock();
    while( clock() != start_time + wait ) {};
    return;
}


/********************************************************************
 * @brief function to copy memory
 ********************************************************************/
PetscErrorCode gencpy(ScalarType* dst, ScalarType* src, size_t size) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
#if defined(REG_HAS_MPICUDA)
    cudaMemcpy((void*)dst , (const void*)src, size, cudaMemcpyDeviceToDevice);
#else
    memcpy((void*)dst, (void*)src, size);
#endif
    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief setup library
 *******************************************************************/
PetscErrorCode InitializeDataDistribution(int nthreads, int *c_grid, MPI_Comm& c_comm, bool c_exists) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;

    PetscFunctionBegin;

    omp_set_dynamic(0);
    omp_set_num_threads(nthreads);
    // check if number of threads is consistent with user options
    int ompthreads = omp_get_max_threads();
    ss << "openmp threads (user,set)=("
       << nthreads <<"," << ompthreads << ")\n";
    ierr = Assert(ompthreads == nthreads, ss.str().c_str()); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();

    int nprocs;
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
#if defined(REG_HAS_CUDA)
    c_grid[0] = nprocs; c_grid[1] = 1;
    int period[2] = {0, 0};
    int reorder = 0;
    if (c_comm != MPI_COMM_NULL)  MPI_Comm_free(&c_comm);
    MPI_Cart_create(PETSC_COMM_WORLD, 2, c_grid, period, reorder, &c_comm);
#else
    int val;
    // set up MPI/cartesian grid
    int np = c_grid[0]*c_grid[1];
    if (np != nprocs) {
        // update cartesian grid layout
        c_grid[0] = 0; c_grid[1] = 0;
        rval = MPI_Dims_create(nprocs, 2, c_grid);
        ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    }

    if (c_exists) {
        rval = MPI_Comm_free(&c_comm);
        ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    }
    
    // initialize accfft
    accfft_create_comm(PETSC_COMM_WORLD, c_grid, &c_comm);
    accfft_init(nthreads);
    accfft_init();
#endif

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief view vector entries (transpose output)
 *******************************************************************/
PetscErrorCode Finalize() {
    PetscErrorCode ierr = 0;

#if !defined(REG_HAS_CUDA)
    accfft_cleanup();
#endif

    // clean up petsc
    ierr = PetscFinalize(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief ensure that sum accross all image components is in [0,1]
 * for a particular location x
 *******************************************************************/
PetscErrorCode ComputeBackGround(Vec background, Vec x, IntType nc) {
    PetscErrorCode ierr = 0;
    ScalarType *p_x = NULL, *p_b = NULL, sum, bval;
    IntType nl, l;

    PetscFunctionBegin;

    ierr = VecGetLocalSize(x, &nl); CHKERRQ(ierr);
    nl /= nc;

    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
    ierr = VecGetArray(background, &p_b); CHKERRQ(ierr);
    for (IntType i = 0; i < nl; ++i) {
        sum = 0.0; // reset sum
        for (IntType k = 0; k < nc; ++k) {
            l = k*nl + i;
            sum += p_x[l]; // sum all components
        }

        // compute background value
        bval = 1.0 - sum;
        if (bval <= 0.0) bval = 0.0;
        if (bval >= 1.0) bval = 1.0;

        // assign background value
        p_b[i] = bval;

    }  // for all components

    ierr = VecRestoreArray(background, &p_b); CHKERRQ(ierr);
    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief rescale data to [xminout,xmaxout]
 *******************************************************************/
PetscErrorCode Rescale(Vec x, ScalarType xminout, ScalarType xmaxout, IntType nc) {
    PetscErrorCode ierr = 0;
    ScalarType xmin, xmax, xmin_g, xmax_g, xscale, xshift, *p_x = NULL;
    IntType nl, l;
    int rval;
    std::stringstream ss;

    PetscFunctionBegin;


    if (nc == 1) {
        // get max and min values
        ierr = VecMin(x, NULL, &xmin); CHKERRQ(ierr);
        ierr = VecMax(x, NULL, &xmax); CHKERRQ(ierr);

//        ss << "rescaling intensities: [" << xmin << "," << xmax << "]"
//           << " -> [" << xminout << "," << xmaxout << "]" ;
//        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
//        ss.clear(); ss.str(std::string());

        xshift = -xmin;

        if (xshift > 1e-6) {
            ierr = VecShift(x, xshift); CHKERRQ(ierr);
        }

        if (std::abs(xmax - xmin) > 1e-6) {
            xscale = xmaxout / (xmax - xmin);
        } else {
            xscale = 1.0;
        }
        ierr = VecScale(x, xscale); CHKERRQ(ierr);
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

            // compute shift and scale
            xshift = xminout - xmin_g;
            xmax_g = (xmax_g  < 1e-5) ? 1.0 : xmax_g;
            xscale = (xmaxout < 1e-5) ? 1.0 : xmaxout / xmax_g;

//            ss << "rescaling intensities: [" << xmin << "," << xmax << "]"
//               << " -> [" << xminout << "," << xmaxout << "]" ;
//            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
//            ss.clear(); ss.str(std::string());

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
 * @brief rescale data to [xminout,xmaxout]
 *******************************************************************/
PetscErrorCode VecNorm(Vec x, IntType nc) {
    PetscErrorCode ierr = 0;
    ScalarType minval, maxval, value, xmin_g, xmax_g, *p_x = NULL;
    IntType nl, l;
    int rval;
    std::stringstream ss;

    PetscFunctionBegin;

    if (nc == 1) {
        // get max and min values
        ierr = VecMin(x, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(x, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecNorm(x, NORM_2, &value); CHKERRQ(ierr);
        ss << "(norm,min,max) = (" << std::scientific << value
           << "," << minval << "," << maxval << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    } else {
        // compute local size from input vector
        ierr = VecGetLocalSize(x, &nl); CHKERRQ(ierr);
        nl /= nc;
        ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {
            minval = std::numeric_limits<ScalarType>::max();
            maxval = std::numeric_limits<ScalarType>::min();

            // get min and max values
            for (IntType i = 0; i < nl; ++i) {
                l = k*nl + i;
                if (p_x[l] < minval) {minval = p_x[l];}
                if (p_x[l] > maxval) {maxval = p_x[l];}
            }

            // get min accross all procs
            rval = MPI_Allreduce(&minval, &xmin_g, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD);
            ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

            // get max accross all procs
            rval = MPI_Allreduce(&maxval, &xmax_g, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD);
            ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

            // TODO: compute norm
            ss << "(norm,min,max) = (" << std::scientific << 0.0
               << "," << xmin_g << "," << xmax_g << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());

        }  // for all components
        ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);
    }  // if else

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief parse string of NUMxNUMxNUM into a vector
 *******************************************************************/
std::vector<ScalarType> String2VecScalarType(const std::string & str, std::string sep) {
    std::vector<ScalarType> vec;
    std::stringstream ss(str);
    std::string num;

    while (std::getline(ss, num, ','))
        vec.push_back(std::stod(num));
    return vec;
}


/********************************************************************
 * @brief parse string of NUMxNUMxNUM into a vector
 *******************************************************************/
std::vector<int> String2Vec(const std::string & str) {
    std::vector<int> vect;
    int ival;
    std::string::size_type xpos = str.find('x',0);

    if (xpos == std::string::npos) {
        // only one uint
        vect.push_back(static_cast<int>(atoi(str.c_str())));
        return vect;
    }

    // first uint$
    ival = atoi((str.substr(0, xpos)).c_str());
    vect.push_back(static_cast<int>(ival));

    while (true) {
        std::string::size_type newxpos = xpos;
        xpos = str.find('x', newxpos+1);

        if (xpos == std::string::npos) {
            ival = atoi((str.substr(newxpos+1, str.length()-newxpos-1)).c_str());
            vect.push_back(static_cast<int>(ival));
            return vect;
        }
        ival = atoi((str.substr(newxpos+1, xpos-newxpos-1)).c_str() );
        vect.push_back(static_cast<int>(ival));
    }
}



/********************************************************************
 * @brief parse string of NUMxNUMxNUM into a vector
 *******************************************************************/
std::vector<int> String2Vec(const std::string & str, std::string sep) {
    std::vector<int> vect;
    int ival;
    std::string::size_type xpos = str.find(sep,0);

    if (xpos == std::string::npos) {
        // only one uint
        vect.push_back(static_cast<int>(atoi(str.c_str())));
        return vect;
    }

    // first uint$
    ival = atoi((str.substr(0, xpos)).c_str());
    vect.push_back(static_cast<int>(ival));

    while (true) {
        std::string::size_type newxpos = xpos;
        xpos = str.find(sep, newxpos+1);

        if (xpos == std::string::npos) {
            ival = atoi((str.substr(newxpos+1, str.length()-newxpos-1)).c_str());
            vect.push_back(static_cast<int>(ival));
            return vect;
        }
        ival = atoi((str.substr(newxpos+1, xpos-newxpos-1)).c_str() );
        vect.push_back(static_cast<int>(ival));
    }
}


/********************************************************************
 * @brief return/restore raw pointers for vector for write/read purpose
 *******************************************************************/
PetscErrorCode GetRawPointer(Vec v, ScalarType** a) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    #if defined(REG_HAS_CUDA)
        ierr = VecCUDAGetArray(v, a); CHKERRQ(ierr);
    #else
        ierr = VecGetArray(v, a); CHKERRQ(ierr);
    #endif

    PetscFunctionReturn(ierr);
}

PetscErrorCode RestoreRawPointer(Vec v, ScalarType** a) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    #if defined(REG_HAS_CUDA)
        ierr = VecCUDARestoreArray(v, a); CHKERRQ(ierr);
    #else
        ierr = VecRestoreArray(v, a); CHKERRQ(ierr);
    #endif

    PetscFunctionReturn(ierr);
}

PetscErrorCode GetRawPointerRead(Vec v, const ScalarType** a) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;
    
    #if defined(REG_HAS_CUDA)
        ierr = VecCUDAGetArrayRead(v, a); CHKERRQ(ierr);
    #else
        ierr = VecGetArrayRead(v, a); CHKERRQ(ierr);
    #endif

    PetscFunctionReturn(ierr);
}


PetscErrorCode RestoreRawPointerRead(Vec v, const ScalarType** a) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;
    
    #if defined(REG_HAS_CUDA)
        ierr = VecCUDARestoreArrayRead(v, a); CHKERRQ(ierr);
    #else
        ierr = VecRestoreArrayRead(v, a); CHKERRQ(ierr);
    #endif

    PetscFunctionReturn(ierr);
}

PetscErrorCode GetRawPointerReadWrite(Vec v, ScalarType** a) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;
    
    ierr = GetRawPointer(v, a); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode RestoreRawPointerReadWrite(Vec v, ScalarType** a) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;
    
    ierr = RestoreRawPointer(v, a); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

PetscErrorCode GetRawPointerWrite(Vec v, ScalarType** a) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;
    
    #ifdef REG_HAS_CUDA
        ierr = VecCUDAGetArrayWrite(v, a); CHKERRQ(ierr);
    #else
        ierr = VecGetArray(v, a); CHKERRQ(ierr);
    #endif

    PetscFunctionReturn(ierr);
}

PetscErrorCode RestoreRawPointerWrite(Vec v, ScalarType** a) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;
  
    #ifdef REG_HAS_CUDA
        ierr = VecCUDARestoreArrayWrite(v, a); CHKERRQ(ierr);
    #else
        ierr = VecRestoreArray(v, a); CHKERRQ(ierr);
    #endif
    
    PetscFunctionReturn(ierr);
}


PetscErrorCode PrintVectorMemoryLocation(Vec v, std::string msg) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}

}  //  namespace reg




#endif   // _CLAIREUTILS_CPP_
