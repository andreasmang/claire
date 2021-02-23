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
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


#ifndef _CLAIREUTILS_HPP_
#define _CLAIREUTILS_HPP_

//#define _REG_DEBUG_

// global includes
#include <fstream>
#include <iomanip>
#include <sstream>
#include <omp.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <ctime>
#include <limits>

#include <sys/types.h>
#include <sys/stat.h>

#include "zeitgeist.hpp"


#ifdef REG_HAS_PNETCDF
#include "pnetcdf.h"
#endif


// local includes
#include "TypeDef.hpp"
#include "MemoryUtils.hpp"

#define DBGCHK() printf("dbg %s:[%i] in %s\n",__FILE__,__LINE__,__FUNCTION__)

namespace reg {


/*! assert (PETSc interface) */
PetscErrorCode Assert(bool, std::string);

/*! throw error (PETSc interface) */
PetscErrorCode ThrowErrorMsg(std::string, int, const char*);
/*! throw error (PETSc interface) */
PetscErrorCode ThrowErrorMsg(std::bad_alloc&, int, const char*);
PetscErrorCode ThrowErrorMsg(std::exception&, int, const char*);


/*! mpi error handling */
PetscErrorCode MPIERRQ(int);
#ifdef REG_HAS_PNETCDF
PetscErrorCode NCERRQ(int);
#endif

/*! check if file exists */
void isleep(unsigned int);

/*! check if file exists */
bool FileExists(const std::string&);

/*! display message (PETSc interface) */
PetscErrorCode Msg(std::string);


/*! display warning message (PETSc interface) */
PetscErrorCode WrngMsgCall(std::string, int, const char*);

/*! display dgb message (PETSc interface) */
//PetscErrorCode DbgMsgCall(std::string);
PetscErrorCode DbgMsgCall(std::string, int=0, const char(*)=0, int=0);

/*! interface to create a vector (essentially simplifies
 * the petsc vector creation) */
PetscErrorCode VecCreate(Vec&, IntType, IntType);

/*! display scalar field */
PetscErrorCode VecView(Vec);

/*! compute norm of vector field */
PetscErrorCode VecNorm(Vec, IntType);

/*! compute norm of vector field */
PetscErrorCode ShowValues(Vec, IntType nc = 1);

/*! rescale vector field to given bounds [xmin,xmax] */
PetscErrorCode Rescale(Vec, ScalarType, ScalarType, IntType nc = 1);

/*! normalize field to [0,1] */
PetscErrorCode Normalize(Vec, IntType nc = 1);

/*! compte pointwise norm of vector field */
PetscErrorCode VecFieldPointWiseNorm(Vec, Vec, Vec, Vec);
PetscErrorCode VecFieldPointWiseNormGPU(ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*, IntType);

/*! clip field to [0,1] */
PetscErrorCode Clip(Vec, IntType nc = 1);

/* Copy vector field to a flat array in strided fasion */
PetscErrorCode CopyStridedToFlatVec(ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*, IntType);

/* Copy from a strided flat array to vecfield strided fasion */
PetscErrorCode CopyStridedFromFlatVec(ScalarType*, ScalarType*, ScalarType*, const ScalarType*, IntType);

/*! ensure partition of unity */
PetscErrorCode EnsurePartitionOfUnity(Vec, IntType);

/*! generic copy function based CPU or GPU implementation */
PetscErrorCode gencpy(ScalarType*, ScalarType*, size_t);

PetscErrorCode SetValue(ScalarType*, ScalarType, IntType);

/*! ensure partition of unity */
PetscErrorCode ComputeBackGround(Vec, Vec, IntType);

PetscErrorCode GetFileName(std::string&, std::string);

PetscErrorCode GetFileName(std::string&, std::string&, std::string&, std::string);

std::vector<int> String2Vec(const std::string&);
std::vector<int> String2Vec(const std::string&, std::string);
std::vector<ScalarType> String2VecScalarType(const std::string&, std::string);

PetscErrorCode InitializeDataDistribution(int, int*, MPI_Comm&, bool);

PetscErrorCode Finalize();

/* get raw pointer to write, read and read,write */
PetscErrorCode GetRawPointer(Vec, ScalarType**);
PetscErrorCode RestoreRawPointer(Vec, ScalarType**);
PetscErrorCode GetRawPointerRead(Vec, const ScalarType**);
PetscErrorCode RestoreRawPointerRead(Vec, const ScalarType**);
PetscErrorCode GetRawPointerReadWrite(Vec, ScalarType**);
PetscErrorCode RestoreRawPointerReadWrite(Vec, ScalarType**);
PetscErrorCode GetRawPointerWrite(Vec, ScalarType**);
PetscErrorCode RestoreRawPointerWrite(Vec, ScalarType**);
PetscErrorCode PrintVectorMemoryLocation(Vec, std::string);


inline PetscErrorCode DebugInfo(Vec vec, std::string str, int line, const char* file) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
      ScalarType maxval, minval, nvx1;
      std::stringstream ss;
      
      ScalarType *p_x1;
      IntType nl;
      size_t fingerprint = 0;
      ierr = VecGetArray(vec, &p_x1); CHKERRQ(ierr);
      // compute size of each individual component
      ierr = VecGetLocalSize(vec, &nl); CHKERRQ(ierr);
  #pragma omp parallel for
      for (IntType i = 0; i < nl; ++i) {
#if defined(PETSC_USE_REAL_SINGLE)
        fingerprint += reinterpret_cast<uint32_t*>(p_x1)[i];
#else
        fingerprint += reinterpret_cast<uint64_t*>(p_x1)[i];
#endif
      }
      ierr = VecRestoreArray(vec, &p_x1); CHKERRQ(ierr);
      
      ss  << str << " hash: " << std::hex << fingerprint;
      ierr = DbgMsgCall(ss.str(), line, file); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();

      ierr = VecNorm(vec, NORM_2, &nvx1); CHKERRQ(ierr);
      ss  << str << " 2-norm: " << std::scientific
          << nvx1;
      ierr = DbgMsgCall(ss.str()); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();
      
      ierr = VecMax(vec, NULL, &maxval); CHKERRQ(ierr);
      ierr = VecMin(vec, NULL, &minval); CHKERRQ(ierr);
      ss  << str << " min/max: [" << std::scientific
          << minval << "," << maxval << "]";
      ierr = DbgMsgCall(ss.str()); CHKERRQ(ierr);
      ss.str(std::string()); ss.clear();
  PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief map 3d index to linear index (accfft style)
 *******************************************************************/
inline IntType GetLinearIndex(IntType i[3], IntType isize[3]) {
    // row major order (ACCFFT)$
    return i[0]*isize[1]*isize[2] + i[1]*isize[2] + i[2];
}




/********************************************************************
 * @brief map 3d index to linear index (accfft style)
 *******************************************************************/
inline IntType GetLinearIndex(IntType i, IntType j, IntType k, IntType isize[3]) {
    // row major order (ACCFFT)$
    return i*isize[1]*isize[2] + j*isize[2] + k;
}




/********************************************************************
 * @brief check wave numbers
 *******************************************************************/
inline void CheckWaveNumbersInv(long int w[3], int n[3]) {
    if (w[0] > n[0]/2) w[0] -= n[0];
    else if (w[0] == n[0]/2) w[0]  = 0;
    if (w[1] > n[1]/2) w[1] -= n[1];
    else if (w[1] == n[1]/2) w[1]  = 0;
    if (w[2] > n[2]/2) w[2] -= n[2];
    else if (w[2] == n[2]/2) w[2]  = 0;
}




/********************************************************************
 * @brief check wave numbers
 *******************************************************************/
inline void CheckWaveNumbers(long int w[3], int n[3]) {
    if      (w[0] >  n[0]/2) w[0] -= n[0];
    else if (w[0] == n[0]/2) w[0]  = 0;
    if      (w[1] >  n[1]/2) w[1] -= n[1];
    else if (w[1] == n[1]/2) w[1]  = 0;
    if      (w[2] >  n[2]/2) w[2] -= n[2];
    else if (w[2] == n[2]/2) w[2]  = 0;
};


/********************************************************************
 * @brief check wave numbers
 *******************************************************************/
inline void ComputeWaveNumber(IntType w[3], IntType n[3]) {
    if      (w[0] >  n[0]/2) w[0] -= n[0];
    else if (w[0] == n[0]/2) w[0]  = 0;
    if      (w[1] >  n[1]/2) w[1] -= n[1];
    else if (w[1] == n[1]/2) w[1]  = 0;
    if      (w[2] >  n[2]/2) w[2] -= n[2];
    else if (w[2] == n[2]/2) w[2]  = 0;
};

/********************************************************************
 * @brief check wave number
 *******************************************************************/
inline void ComputeWaveNumber(IntType &w, IntType n) {
    if (w > n/2) w -= n;
};




/********************************************************************
 * @brief map 3d index to linear index (accfft style)
 *******************************************************************/
inline bool OnMaster() {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    return (rank == 0);
};


}   // namespace reg




#endif  // _CLAIREUTILS_HPP_
