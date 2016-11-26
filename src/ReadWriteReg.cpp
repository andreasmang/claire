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
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _READWRITEREG_CPP_
#define _READWRITEREG_CPP_




#include "ReadWriteReg.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadWriteReg"
ReadWriteReg::ReadWriteReg() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadWriteReg"
ReadWriteReg::ReadWriteReg(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~ReadWriteReg"
ReadWriteReg::~ReadWriteReg() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode ReadWriteReg::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;
    this->m_Data = NULL;

    this->m_iSizeC = NULL;
    this->m_iStartC = NULL;

    this->m_nSend = NULL;
    this->m_nOffset = NULL;

    this->m_nx[0] = -1;
    this->m_nx[1] = -1;
    this->m_nx[2] = -1;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clear class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode ReadWriteReg::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Data != NULL) {
        delete [] this->m_Data; this->m_Data = NULL;
    }
    if (this->m_iSizeC != NULL) {
        delete [] this->m_iSizeC; this->m_iSizeC = NULL;
    }
    if (this->m_iStartC != NULL) {
        delete [] this->m_iStartC; this->m_iStartC = NULL;
    }
    if (this->m_nOffset != NULL) {
        delete [] this->m_nOffset; this->m_nOffset = NULL;
    }
    if (this->m_nSend != NULL) {
        delete [] this->m_nSend; this->m_nSend = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Read"
PetscErrorCode ReadWriteReg::Read(Vec* x, std::string filename) {
    PetscErrorCode ierr = 0;
    std::string file, msg;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // get file name without path
    ierr = GetFileName(file, filename); CHKERRQ(ierr);

    // check if file exists
    msg = "file " + file + " does not exist";
    ierr = Assert(FileExists(filename), msg); CHKERRQ(ierr);

    // display what we are doing
    if (this->m_Opt->GetVerbosity() > 2) {
        msg = "reading " + file;
        ierr = DbgMsg(msg); CHKERRQ(ierr);
    }

    if (filename.find(".nii") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->ReadNII(x,filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (filename.find(".nii.gz") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->ReadNII(x,filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (filename.find(".hdr") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->ReadNII(x,filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (filename.find(".bin") != std::string::npos) {
        ierr = this->ReadBIN(x,filename); CHKERRQ(ierr);
    } else if (filename.find(".hdf5") != std::string::npos) {
#if defined(PETSC_HAVE_HDF5)
        ierr = this->ReadHDF5(x,filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install hdf5 library (petsc)/enable hdf5 support"); CHKERRQ(ierr);
#endif
    } else if (filename.find(".nc") != std::string::npos) {
#ifdef REG_HAS_PNETCDF
        ierr = this->ReadNC(x,filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install pnetcdf library/enable pnetcdf support"); CHKERRQ(ierr);
#endif
    } else {
        ierr = ThrowError("could not read: data type not supported"); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Read"
PetscErrorCode ReadWriteReg::Read(VecField* v, std::string fnx1,
                                    std::string fnx2, std::string fnx3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = this->Read(&v->m_X1, fnx1); CHKERRQ(ierr);
    ierr = this->Read(&v->m_X1, fnx2); CHKERRQ(ierr);
    ierr = this->Read(&v->m_X1, fnx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write time series data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteTimeSeries"
PetscErrorCode ReadWriteReg::WriteTimeSeries(Vec x, std::string filename) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);

    filename = this->m_Opt->GetReadWriteFlags().xfolder + filename;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read temporal data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadTimeSeries"
PetscErrorCode ReadWriteReg::ReadTimeSeries(Vec x, std::string filename) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadBlock"
PetscErrorCode ReadWriteReg::ReadBlock(Vec x, int isize[3], std::string filename) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteBlock"
PetscErrorCode ReadWriteReg::WriteBlock(Vec x, int isize[3], std::string filename) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    filename = this->m_Opt->GetReadWriteFlags().xfolder + filename;
    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write multi component image to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteMC"
PetscErrorCode ReadWriteReg::WriteMC(Vec x, std::string filename) {
    PetscErrorCode ierr;
    IntType nc, nl, ng;
    Vec xk = NULL;
    ScalarType *p_xk = NULL, *p_x = NULL;
    std::stringstream ss;
    std::string msg, path, file, ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);

    // get file name without path
    ierr = GetFileName(file, filename); CHKERRQ(ierr);

    // display what we are doing
    if (this->m_Opt->GetVerbosity() > 2) {
        msg = "writing " + file;
        ierr = DbgMsg(msg); CHKERRQ(ierr);
    }

    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    // allocate data
    ierr = VecCreate(xk, nl, ng); CHKERRQ(ierr);

    filename = this->m_Opt->GetReadWriteFlags().xfolder + "/" + filename;
    ierr = GetFileName(path, file, ext, filename); CHKERRQ(ierr);
    filename = file;

    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
    for (IntType k = 0; k < nc; ++k) {
        ierr = VecGetArray(xk, &p_xk); CHKERRQ(ierr);
        try {std::copy(p_x+k*nl, p_x+(k+1)*nl, p_xk);}
        catch(std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(xk, &p_xk); CHKERRQ(ierr);

        // construct file names for velocity field components
        ss  << filename << "-c=" << std::setw(3) << std::setfill('0') << k << ext;
        ierr = this->Write(xk, ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }
    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);


    if (xk != NULL) {ierr = VecDestroy(&xk); CHKERRQ(ierr); xk = NULL;}

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Write"
PetscErrorCode ReadWriteReg::Write(Vec x, std::string filename) {
    PetscErrorCode ierr;
    std::string file, msg;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);

    // get file name without path
    ierr = GetFileName(file, filename); CHKERRQ(ierr);

    // display what we are doing
    if (this->m_Opt->GetVerbosity() > 2) {
        msg = "writing " + file;
        ierr = DbgMsg(msg); CHKERRQ(ierr);
    }

    filename = this->m_Opt->GetReadWriteFlags().xfolder + filename;
    if (filename.find(".nii") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->WriteNII(x, filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (filename.find(".nii.gz") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->WriteNII(x, filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (filename.find(".hdr") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->WriteNII(x, filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (filename.find(".bin") != std::string::npos) {
        ierr = this->WriteBIN(x, filename); CHKERRQ(ierr);
    } else if (filename.find(".hdf5") != std::string::npos) {
#if defined(PETSC_HAVE_HDF5)
        ierr = this->WriteHDF5(x, filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install hdf library/enable hdf5 support"); CHKERRQ(ierr);
#endif
    } else if (filename.find(".nc") != std::string::npos) {
#ifdef REG_HAS_PNETCDF
        ierr = this->WriteNC(x, filename); CHKERRQ(ierr);
#else
        ierr = ThrowError("install pnetcdf library/enable pnetcdf support"); CHKERRQ(ierr);
#endif
    } else {
        ierr = ThrowError("could not write: data type not supported"); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Write"
PetscErrorCode ReadWriteReg::Write(VecField* v,
                                   std::string fnx1,
                                   std::string fnx2,
                                   std::string fnx3) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = this->Write(v->m_X1, fnx1); CHKERRQ(ierr);
    ierr = this->Write(v->m_X2, fnx2); CHKERRQ(ierr);
    ierr = this->Write(v->m_X3, fnx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get component type of NII images
 *******************************************************************/
#ifdef REG_HAS_NIFTI
#undef __FUNCT__
#define __FUNCT__ "GetComponentTypeNII"
PetscErrorCode ReadWriteReg::GetComponentTypeNII(nifti_image* niiimage) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    switch (niiimage->datatype) {
        case NIFTI_TYPE_UINT8:
        {
            this->m_ComponentType = UCHAR;
            break;
        }
        case NIFTI_TYPE_INT8:
        {
            this->m_ComponentType = CHAR;
            break;
        }
        case NIFTI_TYPE_UINT16:
        {
            this->m_ComponentType = USHORT;
            break;
        }
        case NIFTI_TYPE_INT16:
        {
            this->m_ComponentType = SHORT;
            break;
        }
        case NIFTI_TYPE_UINT32:
        {
            this->m_ComponentType = UINT;
            break;
        }
        case NIFTI_TYPE_INT32:
        {
            this->m_ComponentType = INT;
            break;
        }
        case NIFTI_TYPE_FLOAT32:
        {
            this->m_ComponentType = FLOAT;
            break;
        }
        case NIFTI_TYPE_FLOAT64:
        {
            this->m_ComponentType = DOUBLE;
            break;
        }
        default:
        {
            ierr = ThrowError("image data not supported"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief read nifty image
 *******************************************************************/
#ifdef REG_HAS_NIFTI
#undef __FUNCT__
#define __FUNCT__ "ReadNII"
PetscErrorCode ReadWriteReg::ReadNII(Vec* x, std::string filename) {
    PetscErrorCode ierr = 0;
    std::string msg, file;
    std::stringstream ss;
    int nprocs, rank, rval;
    IntType nx[3], isize[3], istart[3];
    IntType ng, ngx, nl, i1, i2, i3, j1, j2, j3, l, k;
    ScalarType *p_x = NULL, *comdata = NULL;
    nifti_image *image = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    // get file name without path
    ierr = GetFileName(file, filename); CHKERRQ(ierr);

    // read header file
    image = nifti_image_read(filename.c_str(), false);
    msg = "could not read nifti image " + file;
    ierr = Assert(image != NULL, msg); CHKERRQ(ierr);

    // get number of grid points
    nx[0] = static_cast<IntType>(image->nx);
    nx[1] = static_cast<IntType>(image->ny);
    nx[2] = static_cast<IntType>(image->nz);

    for (int i = 0; i < 3; ++i) {
        if (nx[i] % 2 != 0) {
            std::cout << "grid size is not odd " << nx[i] << std::endl;
//            nx[i]++;
        }
    }

    // if we read images, we want to make sure that they have the same size
    if (   (this->m_nx[0] == -1)
        && (this->m_nx[1] == -1)
        && (this->m_nx[2] == -1) ) {

        for (int i = 0; i < 3; ++i) {
            this->m_nx[i] = nx[i];
        }

        if(this->m_Opt->GetVerbosity() > 2) {
            ss << "grid size (" << nx[0] << "," << nx[1] << "," << nx[2] << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        }
    } else {
        msg = "grid size of images varies: perform affine registration first";
        for (int i = 0; i < 3; ++i) {
            ierr = Assert(this->m_nx[i] == nx[i], msg); CHKERRQ(ierr);
        }
    }

    // pass number of grid points to options
    for (int i = 0; i < 3; ++i) {
        this->m_Opt->SetNumGridPoints(i, nx[i]);
    }

    // do the setup before running the code (this essentially
    // concerns the memory distribution/the setup of accfft
    if ( !this->m_Opt->SetupDone() ) {
        ierr = this->m_Opt->DoSetup(); CHKERRQ(ierr);
    }

    // get local size
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;
    for (int i = 0; i < 3; ++i) {
        isize[i] = this->m_Opt->GetDomainPara().isize[i];
        istart[i] = this->m_Opt->GetDomainPara().istart[i];
    }

    //check global size
    ngx = 1;
    for (int i = 0; i < 3; ++i) {
        ngx *= nx[i];
    }
    ierr = Assert(ng == ngx, "global size mismatch"); CHKERRQ(ierr);

    // allocate vector
    if (*x != NULL) {ierr = VecDestroy(x); CHKERRQ(ierr);}
    ierr = VecCreate(PETSC_COMM_WORLD, x); CHKERRQ(ierr);
    ierr = VecSetSizes(*x, nl, ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(*x); CHKERRQ(ierr);


    // read data only on master rank
    if (rank == 0) {
        // allocate data buffer
        if (this->m_Data != NULL) {
            delete this->m_Data;
            this->m_Data = NULL;
        }

        // allocate data buffer
        try {this->m_Data = new ScalarType[ng];}
        catch (std::bad_alloc&) {
            ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        ierr = this->ReadNII(image, filename); CHKERRQ(ierr);

        // get all the sizes to read and assign data correctly
        if (this->m_iSizeC == NULL) {
            try {this->m_iSizeC = new IntType[3*nprocs];}
            catch(std::bad_alloc&) {
                ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        if (this->m_iStartC == NULL) {
            try {this->m_iStartC = new IntType[3*nprocs];}
            catch(std::bad_alloc&) {
                ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
    }

    if (this->m_nSend == NULL) {
        try{ this->m_nSend = new int[3*nprocs]; }
        catch(std::bad_alloc&) {
             ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_nOffset == NULL) {
        try{ this->m_nOffset = new int[3*nprocs]; }
        catch(std::bad_alloc&) {
             ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // gather isize and istart on master rank
    rval = MPI_Gather(isize, 3, MPIU_INT, this->m_iSizeC, 3, MPIU_INT, 0, PETSC_COMM_WORLD);
    ierr = MPIERRQ(rval); CHKERRQ(ierr);
    rval = MPI_Gather(istart, 3, MPIU_INT, this->m_iStartC, 3, MPIU_INT, 0, PETSC_COMM_WORLD);
    ierr = MPIERRQ(rval); CHKERRQ(ierr);

    // compute offset and number of entries to send
    if (rank == 0) {
        IntType offset = 0;
        for (int p = 0; p < nprocs; ++p) {
            IntType nsend = 1;
            for (int i = 0; i < 3; ++i) {
               nsend *= this->m_iSizeC[p*3+i];
            }
            this->m_nSend[p] = static_cast<int>(nsend);
            this->m_nOffset[p] = offset;
            offset += nsend;
        }

        // allocate data buffer
        try {comdata = new ScalarType[ng];}
        catch (std::bad_alloc&) {
            ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        k = 0;
        for (int p = 0; p < nprocs; ++p) {
            for (i1 = 0; i1 < this->m_iSizeC[3*p+0]; ++i1) {  // x1
                for (i2 = 0; i2 < this->m_iSizeC[3*p+1]; ++i2) {  // x2
                    for (i3 = 0; i3 < this->m_iSizeC[3*p+2]; ++i3) {  // x3
                        j1 = i1 + this->m_iStartC[3*p+0];
                        j2 = i2 + this->m_iStartC[3*p+1];
                        j3 = i3 + this->m_iStartC[3*p+2];

                        l = GetLinearIndex(j1, j2, j3, nx);
                        comdata[k++] = this->m_Data[l];
                    }  // for i1
                }  // for i2
            }  // for i3
        }  // for all procs
    }  // on master rank

    int nrecv = static_cast<int>(nl);

    ierr = VecGetArray(*x, &p_x); CHKERRQ(ierr);
    rval = MPI_Scatterv(comdata, this->m_nSend, this->m_nOffset, MPI_DOUBLE, p_x, nrecv, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    ierr = VecRestoreArray(*x, &p_x); CHKERRQ(ierr);

    if (comdata != NULL) {delete [] comdata; comdata = NULL;}
    if (image != NULL) {nifti_image_free(image); image = NULL;}
    if (this->m_Data != NULL) {delete [] this->m_Data; this->m_Data = NULL;}

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief read nifty image with right component type
 *******************************************************************/
#ifdef REG_HAS_NIFTI
#undef __FUNCT__
#define __FUNCT__ "ReadNII"
PetscErrorCode ReadWriteReg::ReadNII(nifti_image* niiimage, std::string filename) {
    PetscErrorCode ierr;
    std::string msg;
    int rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    msg = "should only be called on master/root rank";
    ierr = Assert(rank==0,msg); CHKERRQ(ierr);

    switch (niiimage->datatype) {
        case NIFTI_TYPE_UINT8:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("reading data of type uint8 (uchar)"); CHKERRQ(ierr);
            }
            this->m_ComponentType = UCHAR;
            ierr = this->ReadNII<unsigned char>(niiimage, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT8:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("reading data of type int8 (char)"); CHKERRQ(ierr);
            }
            this->m_ComponentType = CHAR;
            ierr = this->ReadNII<char>(niiimage, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT16:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("reading data of type uint16 (unsigned short)"); CHKERRQ(ierr);
            }
            this->m_ComponentType = USHORT;
            ierr = this->ReadNII<unsigned short>(niiimage, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT16:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("reading data of type int16 (short)"); CHKERRQ(ierr);
            }
            this->m_ComponentType = SHORT;
            ierr = this->ReadNII<short>(niiimage, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT32:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("reading data of type uint32 (unsigned int)"); CHKERRQ(ierr);
            }
            this->m_ComponentType = UINT;
            ierr = this->ReadNII<unsigned int>(niiimage, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT32:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("reading data of type int32 (int)"); CHKERRQ(ierr);
            }
            this->m_ComponentType = INT;
            ierr = this->ReadNII<int>(niiimage, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT32:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("reading data of type float32 (float)"); CHKERRQ(ierr);
            }
            this->m_ComponentType = FLOAT;
            ierr = this->ReadNII<float>(niiimage, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT64:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("reading data of type float64 (double)"); CHKERRQ(ierr);
            }
            this->m_ComponentType = DOUBLE;
            ierr = this->ReadNII<double>(niiimage, filename); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("image data not supported"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief get component type of NII images
 *******************************************************************/
#ifdef REG_HAS_NIFTI
#undef __FUNCT__
#define __FUNCT__ "ReadNII"
template <typename T> PetscErrorCode ReadWriteReg::ReadNII(nifti_image* niiimage, std::string filename) {
    PetscErrorCode ierr;
    T *data = NULL;
    std::string msg;
    IntType ng;
    int rank;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(this->m_Data != NULL, "null pointer"); CHKERRQ(ierr);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    msg = "should only be called on master/root rank";
    ierr = Assert(rank == 0, msg); CHKERRQ(ierr);

    if (nifti_image_load(niiimage) == -1) {
        msg="could not read image " + filename;
        ierr = ThrowError(msg); CHKERRQ(ierr);
    }

    // assign data
    data = static_cast<T*>(niiimage->data);
    ierr = Assert(data != NULL, "null pointer"); CHKERRQ(ierr);

    // get global number of points
    ng = static_cast<IntType>(this->m_Opt->GetDomainPara().nglobal);

    // copy buffer
    for (IntType i = 0; i < ng; ++i) {
        this->m_Data[i] = static_cast<ScalarType>(data[i]);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief write buffer to nii files
 *******************************************************************/
#ifdef REG_HAS_NIFTI
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
PetscErrorCode ReadWriteReg::WriteNII(Vec x, std::string filename) {
    PetscErrorCode ierr;
    int rank;
    nifti_image* image = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // TODO
    // we want to rescale the image to whatever
    // the input intensity range was
    //rescale = this->m_RescaleImage;
    //this->m_RescaleImage = true;

    // write x buffer to nifti image
    //ierr = this->WriteNII(&this->m_Image,x,filename); CHKERRQ(ierr);
    ierr = this->WriteNII(&image, x, filename); CHKERRQ(ierr);

    //this->m_RescaleImage = rescale;

    // at rank zero write out
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(rank == 0) { nifti_image_write(image); }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief write buffer to nii files
 *******************************************************************/
#ifdef REG_HAS_NIFTI
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
PetscErrorCode ReadWriteReg::WriteNII(nifti_image** niiimage, Vec x, std::string filename) {
    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // if nifty image is null pointer default to double
    if ((*niiimage) == NULL) {
        ierr = this->WriteNII<ScalarType>(niiimage, x, filename); CHKERRQ(ierr);
        this->m_Opt->Exit(__FUNCT__);
        PetscFunctionReturn(0);
    }

    switch ((*niiimage)->datatype) {
        case NIFTI_TYPE_UINT8:
        {
            ierr = this->WriteNII<unsigned char>(niiimage, x, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT8:
        {
            ierr = this->WriteNII<char>(niiimage, x, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT16:
        {
            ierr = this->WriteNII<unsigned short>(niiimage, x, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT16:
        {
            ierr = this->WriteNII<short>(niiimage, x, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT32:
        {
            ierr = this->WriteNII<unsigned int>(niiimage, x, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT32:
        {
            ierr = this->WriteNII<int>(niiimage, x, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT32:
        {
            ierr = this->WriteNII<float>(niiimage, x, filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT64:
        {
            ierr = this->WriteNII<double>(niiimage, x, filename); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("image data not supported"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}
#endif



/********************************************************************
 * @brief write buffer to nii files
 *******************************************************************/
#ifdef REG_HAS_NIFTI
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
template <typename T>
PetscErrorCode ReadWriteReg::WriteNII(nifti_image** image, Vec x, std::string filename) {
    PetscErrorCode ierr;
    IntType ng;
    T* data = NULL;
    ScalarType *p_xc = NULL;
    int nprocs,rank,rval;
    IntType istart[3], isize[3], nx[3];
    IntType i1, i2, i3, j1, j2, j3, l;
    Vec xcollect = NULL;
    VecScatter scatterctx = NULL;
    std::string msg;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // get number of ranks
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    // get local size
    ng = static_cast<IntType>(this->m_Opt->GetDomainPara().nglobal);

    // allocate the index buffers on master rank
    if (rank == 0) {
        // we need to allocate the image if it's a zero pointer; this
        // will also create a standard header file; not tested (might need
        // to parse the dimensions of the data)
        if ((*image) == NULL) {
            if (this->m_Opt->GetVerbosity() >= 4) {
                msg = "allocating buffer for nifti image";
                ierr = DbgMsg(msg); CHKERRQ(ierr);
            }
            ierr = this->AllocateNII(image, x); CHKERRQ(ierr);
        }

        msg="nifty image is null pointer";
        ierr = Assert((*image) != NULL, msg); CHKERRQ(ierr);

        // construct file name
        std::string file(filename);

        // get temp extension
        const char* exttemp = nifti_find_file_extension(file.c_str());
        if (exttemp == NULL) {exttemp = ".nii";}

        // set extension
        const std::string ext(exttemp);

        // get base file name
        char* bnametemp = nifti_makebasename(file.c_str());
        const std::string bname(bnametemp);
        free(bnametemp);

        // is file compressed
        const std::string::size_type sep = ext.rfind(".gz");
        const bool iscompressed = (sep == std::string::npos) ? false : true;

        if ((ext == ".nii") || (ext == ".nii.gz")) {
            (*image)->nifti_type = NIFTI_FTYPE_NIFTI1_1;
        } else if (ext == ".nia") {
            (*image)->nifti_type = NIFTI_FTYPE_ASCII;
        } else if ((ext == ".hdr") || (ext == ".img") || (ext == ".hdr.gz") || (ext == ".img.gz")) {
            (*image)->nifti_type = NIFTI_FTYPE_NIFTI1_2;
        } else {
            ierr = ThrowError("file extension not supported"); CHKERRQ(ierr);
        }

        (*image)->fname = nifti_makehdrname(bname.c_str(), (*image)->nifti_type, false, iscompressed);
        (*image)->iname = nifti_makeimgname(bname.c_str(), (*image)->nifti_type, false, iscompressed);
    }

    // get array
    if (nprocs == 1) {
        ierr = VecGetArray(x, &p_xc); CHKERRQ(ierr);

        // cast pointer of nifti image data
        data = static_cast<T*>((*image)->data);
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < ng; ++i) {
            data[i] = static_cast<T>(p_xc[i]);
        }
}  /// pragma omp parallel
        ierr = VecRestoreArray(x, &p_xc); CHKERRQ(ierr);
    } else {
        // allocate the index buffers on master rank
        if (rank == 0) {
            // get all the sizes to read and assign data correctly
            if (this->m_iSizeC == NULL) {
                try{ this->m_iSizeC = new IntType[3*nprocs]; }
                catch(std::bad_alloc&) {
                    ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
                }
            }
            if (this->m_iStartC == NULL) {
                try{ this->m_iStartC = new IntType[3*nprocs]; }
                catch(std::bad_alloc&) {
                    ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
                }
            }
        }

        isize[0] = this->m_Opt->GetDomainPara().isize[0];
        isize[1] = this->m_Opt->GetDomainPara().isize[1];
        isize[2] = this->m_Opt->GetDomainPara().isize[2];

        istart[0] = this->m_Opt->GetDomainPara().istart[0];
        istart[1] = this->m_Opt->GetDomainPara().istart[1];
        istart[2] = this->m_Opt->GetDomainPara().istart[2];

        // gather the indices
        rval = MPI_Gather(istart, 3, MPIU_INT, this->m_iStartC, 3, MPIU_INT, 0, PETSC_COMM_WORLD);
        ierr = MPIERRQ(rval); CHKERRQ(ierr);

        rval = MPI_Gather(isize, 3, MPIU_INT, this->m_iSizeC, 3, MPIU_INT, 0, PETSC_COMM_WORLD);
        ierr = MPIERRQ(rval); CHKERRQ(ierr);

        // create scatter object
        ierr = VecScatterCreateToZero(x, &scatterctx, &xcollect); CHKERRQ(ierr);

        // gather the data
        ierr = VecScatterBegin(scatterctx, x, xcollect, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(scatterctx, x, xcollect, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

        nx[0] = this->m_Opt->GetNumGridPoints(0);
        nx[1] = this->m_Opt->GetNumGridPoints(1);
        nx[2] = this->m_Opt->GetNumGridPoints(2);

        // if we are on master rank
        if (rank == 0) {
            data = static_cast<T*>((*image)->data);
            ierr = VecGetArray(xcollect, &p_xc); CHKERRQ(ierr);
            IntType k = 0;
            for(int p = 0; p < nprocs; ++p) {
                for (i1 = 0; i1 < this->m_iSizeC[3*p+0]; ++i1) {  // x1
                    for (i2 = 0; i2 < this->m_iSizeC[3*p+1]; ++i2) {  // x2
                        for (i3 = 0; i3 < this->m_iSizeC[3*p+2]; ++i3) {  // x3
                            j1 = i1 + this->m_iStartC[3*p+0];
                            j2 = i2 + this->m_iStartC[3*p+1];
                            j3 = i3 + this->m_iStartC[3*p+2];

                            l = GetLinearIndex(j1, j2, j3, nx);
                            data[l] = static_cast<T>(p_xc[k++]);
                        }  // for i1
                    }  // for i2
                }  // for i3
            }  // for all procs
            ierr = VecRestoreArray(xcollect, &p_xc); CHKERRQ(ierr);
        }   // if on master
    }  // else

    // clear memory
    if (xcollect != NULL) {ierr = VecDestroy(&xcollect); CHKERRQ(ierr);}
    if (scatterctx != NULL) {ierr = VecScatterDestroy(&scatterctx); CHKERRQ(ierr);}

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}
#endif



/********************************************************************
 * @brief allocate buffer for nifty image
 *******************************************************************/
#ifdef REG_HAS_NIFTI
#undef __FUNCT__
#define __FUNCT__ "AllocateNII"
PetscErrorCode ReadWriteReg::AllocateNII(nifti_image** image, Vec x) {
    PetscErrorCode ierr = 0;
    IntType n;
    int rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // init nifty image
    *image = nifti_simple_init_nim();

    // dimensionalty of data: default is 5 (space, time, components)
    (*image)->dim[0] = (*image)->ndim = 5;

    ierr = VecGetLocalSize(x, &n); CHKERRQ(ierr);
    (*image)->dim[1] = (*image)->nx = this->m_Opt->GetDomainPara().nx[0];
    (*image)->dim[2] = (*image)->ny = this->m_Opt->GetDomainPara().nx[1];
    (*image)->dim[3] = (*image)->nz = this->m_Opt->GetDomainPara().nx[2];

    (*image)->pixdim[1] = static_cast<float>(this->m_Opt->GetDomainPara().hx[0]);  // x direction
    (*image)->pixdim[2] = static_cast<float>(this->m_Opt->GetDomainPara().hx[1]);  // y direction
    (*image)->pixdim[3] = static_cast<float>(this->m_Opt->GetDomainPara().hx[2]);  // z direction

    // TODO: add temporal support
    if (n == this->m_Opt->GetDomainPara().nlocal) {  // scalar field
        (*image)->dim[4] = (*image)->nt = 1;
        (*image)->dim[5] = (*image)->nu = 1;

        // temporal step size
        (*image)->pixdim[4] = 1.0;
    } else if (n == 2*this->m_Opt->GetDomainPara().nlocal) {  // 2D vector field
        (*image)->dim[4] = (*image)->nt = 1;
        (*image)->dim[5] = (*image)->nu = 2;

        // temporal step size
        (*image)->pixdim[4] = 1.0;

        // step size (vector field)
        (*image)->pixdim[5] = (*image)->du = static_cast<float>(this->m_Opt->GetDomainPara().hx[0]);
        (*image)->pixdim[6] = (*image)->dv = static_cast<float>(this->m_Opt->GetDomainPara().hx[1]);
    } else if (n == 3*this->m_Opt->GetDomainPara().nlocal) {  // 3D vector field
        (*image)->dim[4] = (*image)->nt = 1;
        (*image)->dim[5] = (*image)->nu = 3;

        // temporal step size
        (*image)->pixdim[4] = 1.0;

        // step size (vector field)
        (*image)->pixdim[5] = (*image)->du = static_cast<float>(this->m_Opt->GetDomainPara().hx[0]);
        (*image)->pixdim[6] = (*image)->dv = static_cast<float>(this->m_Opt->GetDomainPara().hx[1]);
        (*image)->pixdim[7] = (*image)->dw = static_cast<float>(this->m_Opt->GetDomainPara().hx[2]);
    }

    // currently fixed to double precision: TODO flexible...
    (*image)->datatype = NIFTI_TYPE_FLOAT64;
    (*image)->nbyper = sizeof(ScalarType);

    (*image)->nvox = 1;

    // compute number of voxels (space and time)
    for (int i = 1; i <= 4; ++i) {
        (*image)->nvox *= (*image)->dim[i];
    }
    // add components
    (*image)->nvox *= (*image)->nu;

    // only allocate image buffer on root
    if (rank == 0) {
        try {(*image)->data = new ScalarType[(*image)->nvox];}
        catch (std::bad_alloc&) {
            ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief write binary data set to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadBIN"
PetscErrorCode ReadWriteReg::ReadBIN(Vec* x, std::string filename) {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    PetscViewer viewer=NULL;
    PetscFunctionBegin;

    if (*x != NULL) {ierr = VecDestroy(x); CHKERRQ(ierr); *x = NULL;}

    if (!this->m_Opt->SetupDone()) {
        ierr = this->m_Opt->DoSetup(); CHKERRQ(ierr);
    }
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;
    ierr = VecCreate(*x, nl, ng); CHKERRQ(ierr);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = Assert(viewer != NULL, "could not read binary file"); CHKERRQ(ierr);
    ierr = PetscViewerBinarySetFlowControl(viewer, 2); CHKERRQ(ierr);
    ierr = VecLoad(*x, viewer); CHKERRQ(ierr);

    // clean up
    if (viewer!=NULL) {
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        viewer=NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write bin to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteBIN"
PetscErrorCode ReadWriteReg::WriteBIN(Vec x, std::string filename) {
    PetscErrorCode ierr = 0;
    PetscViewer viewer = NULL;
    PetscFunctionBegin;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = Assert(viewer != NULL, "could not write binary file"); CHKERRQ(ierr);
    ierr = VecView(x, viewer); CHKERRQ(ierr);

    // clean up
    if (viewer != NULL) {
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        viewer = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read hdf5 image
 *******************************************************************/
#if defined(PETSC_HAVE_HDF5)
#undef __FUNCT__
#define __FUNCT__ "ReadHDF5"
PetscErrorCode ReadWriteReg::ReadHDF5(Vec* x, std::string filename) {
    PetscErrorCode ierr = 0;
    PetscViewer viewer = NULL;
    PetscFunctionBegin;

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = Assert(viewer!=NULL,"could not read hdf5 file"); CHKERRQ(ierr);
    ierr = VecLoad(*x, viewer); CHKERRQ(ierr);

    // clean up
    if (viewer != NULL) {
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        viewer = NULL;
    }

    PetscFunctionReturn(ierr);
}
#endif




/********************************************************************
 * @brief read hdf5 image
 *******************************************************************/
#if defined(PETSC_HAVE_HDF5)
#undef __FUNCT__
#define __FUNCT__ "WriteHDF5"
PetscErrorCode ReadWriteReg::WriteHDF5(Vec x, std::string filename) {
    PetscErrorCode ierr = 0;
    PetscViewer viewer=NULL;
    PetscFunctionBegin;

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = Assert(viewer!=NULL,"could not write hdf5 file"); CHKERRQ(ierr);
//    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);
    ierr = VecView(x, viewer); CHKERRQ(ierr);

    // clean up
    if (viewer != NULL) {
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        viewer = NULL;
    }

    PetscFunctionReturn(ierr);
}
#endif




/********************************************************************
 * @brief write netcdf to file
 *******************************************************************/
#ifdef REG_HAS_PNETCDF
#undef __FUNCT__
#define __FUNCT__ "ReadNC"
PetscErrorCode ReadWriteReg::ReadNC(Vec* x, std::string filename) {
    PetscErrorCode ierr = 0;
    int rank, ncerr, fileid, ndims, nvars, ngatts, unlimited, varid[1];
    IntType nl, ng;
    ScalarType *p_x = NULL;
//    double *data=NULL;
    MPI_Offset istart[3], isize[3];
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if ((*x) != NULL) {ierr = VecDestroy(x); CHKERRQ(ierr); *x = NULL;}

    if (!this->m_Opt->SetupDone()) {
        ierr = this->m_Opt->DoSetup(); CHKERRQ(ierr);
    }
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;
    ierr = VecCreate(*x, nl, ng); CHKERRQ(ierr);

    // open file
    ncerr = ncmpi_open(PETSC_COMM_WORLD,filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &fileid);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    // query info about field named "data"
    ncerr=ncmpi_inq(fileid, &ndims, &nvars, &ngatts, &unlimited);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ncerr=ncmpi_inq_varid(fileid, "data", &varid[0]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    istart[0] = this->m_Opt->GetDomainPara().istart[0];
    istart[1] = this->m_Opt->GetDomainPara().istart[1];
    istart[2] = this->m_Opt->GetDomainPara().istart[2];

    isize[0] = this->m_Opt->GetDomainPara().isize[0];
    isize[1] = this->m_Opt->GetDomainPara().isize[1];
    isize[2] = this->m_Opt->GetDomainPara().isize[2];

    ierr = VecGetArray(*x, &p_x); CHKERRQ(ierr);
    //data=static_cast<double*>(p_x);
    //ncerr=ncmpi_get_vara_all(fileid,varid[0],istart,isize,data,nl,MPI_DOUBLE);
    ncerr=ncmpi_get_vara_all(fileid, varid[0], istart, isize, p_x, nl, MPI_DOUBLE);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ierr = VecRestoreArray(*x, &p_x); CHKERRQ(ierr);


    ncerr=ncmpi_close(fileid);

    PetscFunctionReturn(ierr);
}
#endif




/********************************************************************
 * @brief write netcdf to file
 *******************************************************************/
#ifdef REG_HAS_PNETCDF
#undef __FUNCT__
#define __FUNCT__ "WriteNC"
PetscErrorCode ReadWriteReg::WriteNC(Vec x, std::string filename) {
    PetscErrorCode ierr = 0;
    int ncerr, mode, dims[3], varid[1], nx[3], iscdf5, fileid;
    int nl;
    MPI_Offset istart[3], isize[3];
    MPI_Comm c_comm;
    ScalarType *p_x = NULL;
    bool usecdf5 = false;   // CDF-5 is mandatory for large files (>= 2x10^9 cells)

    PetscFunctionBegin;

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);

    // file creation mode
    mode=NC_CLOBBER;
    if (usecdf5) {
        mode = NC_CLOBBER | NC_64BIT_DATA;
    } else {
        mode = NC_CLOBBER | NC_64BIT_OFFSET;
    }

    c_comm = this->m_Opt->GetFFT().mpicomm;

    // create netcdf file
    ncerr = ncmpi_create(c_comm, filename.c_str(), mode, MPI_INFO_NULL, &fileid);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    nx[0] = static_cast<int>(this->m_Opt->GetDomainPara().nx[0]);
    nx[1] = static_cast<int>(this->m_Opt->GetDomainPara().nx[1]);
    nx[2] = static_cast<int>(this->m_Opt->GetDomainPara().nx[2]);
    nl = static_cast<int>(this->m_Opt->GetDomainPara().nlocal);

    // set size
    ncerr = ncmpi_def_dim(fileid, "x", nx[0], &dims[0]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ncerr = ncmpi_def_dim(fileid, "y", nx[1], &dims[1]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ncerr = ncmpi_def_dim(fileid, "z", nx[2], &dims[2]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    // define name for output field
    ncerr = ncmpi_def_var(fileid, "data", NC_DOUBLE, 3, dims, &varid[0]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    iscdf5 = usecdf5 ? 1 : 0;
    ncerr = ncmpi_put_att_int(fileid, NC_GLOBAL, "CDF-5 mode", NC_INT, 1, &iscdf5);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ncerr = ncmpi_enddef(fileid);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    // get local sizes
    istart[0] = static_cast<MPI_Offset>( this->m_Opt->GetDomainPara().istart[0] );
    istart[1] = static_cast<MPI_Offset>( this->m_Opt->GetDomainPara().istart[1] );
    istart[2] = static_cast<MPI_Offset>( this->m_Opt->GetDomainPara().istart[2] );

//    std::cout<< istart[0] << " " << istart[1] << " " << istart[2] << std::endl;


    isize[0] = static_cast<MPI_Offset>( this->m_Opt->GetDomainPara().isize[0] );
    isize[1] = static_cast<MPI_Offset>( this->m_Opt->GetDomainPara().isize[1] );
    isize[2] = static_cast<MPI_Offset>( this->m_Opt->GetDomainPara().isize[2] );
//    std::cout<< isize[0] << " " << isize[1] << " " << isize[2] << std::endl;
    ierr = Assert(nl == isize[0]*isize[1]*isize[2], "size error"); CHKERRQ(ierr);

    // write data to file
    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
    ncerr = ncmpi_put_vara_all(fileid, varid[0], istart, isize, p_x, nl, MPI_DOUBLE);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);

    // close file
    ncerr = ncmpi_close(fileid);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}
#endif




}  // namespace reg




#endif  // _READWRITEREG_CPP_
