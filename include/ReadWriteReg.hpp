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

#ifndef _READWRITEREG_H_
#define _READWRITEREG_H_

// library includes
#ifdef REG_HAS_NIFTI
#include "nifti1_io.h"
#endif

#ifdef REG_HAS_PNETCDF
#include "pnetcdf.h"
#endif

#include "RegOpt.hpp"
#include "VecField.hpp"




namespace reg {




enum DataType {CHAR, UCHAR, SHORT, USHORT, INT, UINT, FLOAT, DOUBLE, UNDEF};

struct ImageType {
#ifdef REG_HAS_NIFTI
    nifti_image* data;
#endif
    DataType datatype;

    ScalarType minval;
    ScalarType maxval;
    IntType nx[3];
    bool read;
    bool write;
};




class ReadWriteReg {
 public:
    typedef ReadWriteReg Self;

    ReadWriteReg(void);
    ReadWriteReg(RegOpt*);
    virtual ~ReadWriteReg(void);

    /*! read reference image */
    PetscErrorCode ReadR(Vec*, std::vector < std::string >);

    /*! read template image */
    PetscErrorCode ReadT(Vec*, std::vector < std::string >);

    PetscErrorCode Read(Vec*, std::vector < std::string >);
    PetscErrorCode Read(Vec*, std::string);
    PetscErrorCode Read(VecField*, std::string, std::string, std::string);

    /*! write reference image */
    PetscErrorCode WriteR(Vec, std::string, bool multicomponent = false);

    /*! write template image */
    PetscErrorCode WriteT(Vec, std::string, bool multicomponent = false);

    PetscErrorCode Write(Vec, std::string, bool multicomponent = false);
    PetscErrorCode Write(VecField*, std::string);

 private:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    PetscErrorCode Read(Vec*);

    PetscErrorCode Write(Vec);
    PetscErrorCode Write(VecField*);

#ifdef REG_HAS_PNETCDF
    PetscErrorCode ReadNC(Vec*);
    PetscErrorCode WriteNC(Vec);
#endif

    PetscErrorCode ReadBIN(Vec*);
    PetscErrorCode WriteBIN(Vec);

    PetscErrorCode ReadNetCDF(Vec);
    PetscErrorCode ReadTimeSeriesNetCDF(Vec);
    PetscErrorCode ReadBlockNetCDF(Vec, int*);

    PetscErrorCode WriteNetCDF(Vec);
    PetscErrorCode WriteTimeSeriesNetCDF(Vec);
    PetscErrorCode WriteBlockNetCDF(Vec, int*);

    PetscErrorCode CollectSizes();

#ifdef REG_HAS_NIFTI
    PetscErrorCode ReadNII(Vec*);
    PetscErrorCode ReadNII(VecField*);
    PetscErrorCode ReadNII(nifti_image*);
    template <typename T> PetscErrorCode ReadNII(nifti_image*);

    PetscErrorCode WriteNII(Vec);
    PetscErrorCode WriteNII(nifti_image**);
    template <typename T> PetscErrorCode WriteNII(nifti_image**, Vec);

    PetscErrorCode GetComponentType(nifti_image*, DataType&);;
    PetscErrorCode AllocateImage(nifti_image**, Vec);

    nifti_image* m_ImageData;
#endif

    ImageType m_TemplateImage;
    ImageType m_ReferenceImage;

    RegOpt* m_Opt;
    IntType* m_iSizeC;
    IntType* m_iStartC;
    int* m_nOffset;
    int* m_nSend;
    int m_NumProcs;

    ScalarType* m_Data;
    IntType m_nx[3];

    std::string m_FileName;
};




}  // namespace reg




#endif  // _READWRITEREG_H_
