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

#ifndef _READWRITEREG_H_
#define _READWRITEREG_H_

// library includes
#ifdef REG_HAS_NIFTI
#include "nifti1_io.h"
#endif

#if defined(PETSC_HAVE_HDF5)
#include "petscviewerhdf5.h"
#endif

#ifdef REG_HAS_PNETCDF
#include "pnetcdf.h"
#endif

#include "RegOpt.hpp"
#include "VecField.hpp"




namespace reg {




class ReadWriteReg {
 public:
    typedef ReadWriteReg Self;

    ReadWriteReg(void);
    ReadWriteReg(RegOpt*);
    ~ReadWriteReg(void);

    PetscErrorCode Read(Vec*, std::string);
    PetscErrorCode Read(VecField*, std::string, std::string, std::string);
    PetscErrorCode ReadBlock(Vec, int*, std::string);
    PetscErrorCode ReadTimeSeries(Vec, std::string);

    PetscErrorCode Write(Vec, std::string);
    PetscErrorCode Write(VecField*, std::string, std::string, std::string);
    PetscErrorCode WriteBlock(Vec, int*, std::string);
    PetscErrorCode WriteTimeSeries(Vec, std::string);

 private:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

#if defined(PETSC_HAVE_HDF5)
    PetscErrorCode ReadHDF5(Vec*, std::string);
    PetscErrorCode WriteHDF5(Vec, std::string);
#endif

#ifdef REG_HAS_PNETCDF
    PetscErrorCode ReadNC(Vec*, std::string);
    PetscErrorCode WriteNC(Vec, std::string);
#endif

    PetscErrorCode ReadBIN(Vec*, std::string);
    PetscErrorCode WriteBIN(Vec, std::string);

    PetscErrorCode ReadNetCDF(Vec, std::string);
    PetscErrorCode ReadTimeSeriesNetCDF(Vec, std::string);
    PetscErrorCode ReadBlockNetCDF(Vec, int*, std::string);

    PetscErrorCode WriteNetCDF(Vec, std::string);
    PetscErrorCode WriteTimeSeriesNetCDF(Vec, std::string);
    PetscErrorCode WriteBlockNetCDF(Vec, int*, std::string);

#ifdef REG_HAS_NIFTI
    PetscErrorCode ReadNII(Vec*, std::string);
    PetscErrorCode ReadNII(VecField*, std::string, std::string, std::string);
    PetscErrorCode ReadNII(nifti_image*, std::string);
    template <typename T> PetscErrorCode ReadNII(nifti_image*, std::string);

    PetscErrorCode WriteNII(Vec, std::string);
    PetscErrorCode WriteNII(nifti_image**, Vec, std::string);
    template <typename T> PetscErrorCode WriteNII(nifti_image**, Vec, std::string);

    PetscErrorCode GetComponentTypeNII(nifti_image*);;
    PetscErrorCode AllocateNII(nifti_image**, Vec);
#endif

    enum VoxelType{CHAR, UCHAR, SHORT, USHORT, INT, UINT, FLOAT, DOUBLE, UNDEF};
    VoxelType m_ComponentType;

    RegOpt* m_Opt;
    IntType* m_iSizeC;
    IntType* m_iStartC;
    int* m_nOffset;
    int* m_nSend;

    ScalarType* m_Data;
    IntType m_nx[3];
};




}  // namespace reg




#endif  // _READWRITEREG_H_
