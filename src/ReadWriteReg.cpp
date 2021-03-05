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

#ifndef _READWRITEREG_CPP_
#define _READWRITEREG_CPP_




#include "ReadWriteReg.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
ReadWriteReg::ReadWriteReg() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
ReadWriteReg::ReadWriteReg(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
ReadWriteReg::~ReadWriteReg() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init class variables
 *******************************************************************/
PetscErrorCode ReadWriteReg::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;
    this->m_Data = NULL;

#ifdef REG_HAS_NIFTI
    this->m_ReferenceImage.data = NULL;
#endif
    this->m_ReferenceImage.datatype = DOUBLE;
    this->m_ReferenceImage.nx[0] = -1;
    this->m_ReferenceImage.nx[1] = -1;
    this->m_ReferenceImage.nx[2] = -1;
    this->m_ReferenceImage.read = false;
    this->m_ReferenceImage.write = false;
    this->m_ReferenceImage.minval = -1.0;
    this->m_ReferenceImage.maxval = -1.0;

#ifdef REG_HAS_NIFTI
    this->m_TemplateImage.data = NULL;
#endif
    this->m_TemplateImage.datatype = DOUBLE;
    this->m_TemplateImage.nx[0] = -1;
    this->m_TemplateImage.nx[1] = -1;
    this->m_TemplateImage.nx[2] = -1;
    this->m_TemplateImage.read = false;
    this->m_TemplateImage.write = false;
    this->m_TemplateImage.minval = -1.0;
    this->m_TemplateImage.maxval = -1.0;

#ifdef REG_HAS_NIFTI
    this->m_ImageData = NULL;
#endif

    this->m_NumProcs = 0;
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
PetscErrorCode ReadWriteReg::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Data != NULL) {
        delete [] this->m_Data;
        this->m_Data = NULL;
    }
    if (this->m_iSizeC != NULL) {
        delete [] this->m_iSizeC;
        this->m_iSizeC = NULL;
    }
    if (this->m_iStartC != NULL) {
        delete [] this->m_iStartC;
        this->m_iStartC = NULL;
    }
    if (this->m_nOffset != NULL) {
        delete [] this->m_nOffset;
        this->m_nOffset = NULL;
    }
    if (this->m_nSend != NULL) {
        delete [] this->m_nSend;
        this->m_nSend = NULL;
    }

#ifdef REG_HAS_NIFTI
    if (this->m_ReferenceImage.data != NULL) {
        nifti_image_free(this->m_ReferenceImage.data);
        this->m_ReferenceImage.data = NULL;
    }
    if (this->m_TemplateImage.data != NULL) {
        nifti_image_free(this->m_TemplateImage.data);
        this->m_TemplateImage.data = NULL;
    }

    if (this->m_ImageData != NULL) {
        nifti_image_free(this->m_ImageData);
        this->m_ImageData = NULL;
    }
#endif

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief collect data distribution sizes and collect them on master
 *******************************************************************/
PetscErrorCode ReadWriteReg::CollectSizes() {
    PetscErrorCode ierr = 0;
    int nprocs, rank, rval;
    IntType isize[3], istart[3], offset, nsend;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    this->m_NumProcs = nprocs;

    isize[0] = this->m_Opt->m_Domain.isize[0];
    isize[1] = this->m_Opt->m_Domain.isize[1];
    isize[2] = this->m_Opt->m_Domain.isize[2];

    istart[0] = this->m_Opt->m_Domain.istart[0];
    istart[1] = this->m_Opt->m_Domain.istart[1];
    istart[2] = this->m_Opt->m_Domain.istart[2];

    // read data only on master rank
    if (rank == 0) {
        // get all the sizes to read and assign data correctly
        if (this->m_iSizeC == NULL) {
            try {this->m_iSizeC = new IntType[3*nprocs];}
            catch (std::bad_alloc&) {
                ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        if (this->m_iStartC == NULL) {
            try {this->m_iStartC = new IntType[3*nprocs];}
            catch (std::bad_alloc&) {
                ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
    }

    if (this->m_nSend == NULL) {
        try {this->m_nSend = new int[3*nprocs];}
        catch (std::bad_alloc&) {
             ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_nOffset == NULL) {
        try {this->m_nOffset = new int[3*nprocs];}
        catch (std::bad_alloc&) {
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
        offset = 0;
        for (int p = 0; p < nprocs; ++p) {
            nsend = 1;
            for (int i = 0; i < 3; ++i) {
               nsend *= this->m_iSizeC[p*3+i];
            }
            this->m_nSend[p] = static_cast<int>(nsend);
            this->m_nOffset[p] = offset;
            offset += nsend;
        }
    }  // on master rank

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief read reference image from filename
 *******************************************************************/
PetscErrorCode ReadWriteReg::ReadR(Vec* x, std::vector< std::string > filenames) {
    PetscErrorCode ierr = 0;
    ScalarType maxval, minval;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

#ifdef REG_HAS_NIFTI
    if (this->m_ReferenceImage.data != NULL) {
        nifti_image_free(this->m_ReferenceImage.data);
        this->m_ReferenceImage.data = NULL;
    }
#endif
    this->m_ReferenceImage.read = true;

    ierr = this->Read(x, filenames); CHKERRQ(ierr);

    this->m_ReferenceImage.read = false;

    ierr = VecMin(*x, NULL, &minval); CHKERRQ(ierr);
    ierr = VecMax(*x, NULL, &maxval); CHKERRQ(ierr);
    this->m_ReferenceImage.minval = minval;
    this->m_ReferenceImage.maxval = maxval;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read reference image from filename
 *******************************************************************/
PetscErrorCode ReadWriteReg::ReadT(Vec* x, std::vector< std::string > filenames) {
    PetscErrorCode ierr = 0;
    ScalarType maxval, minval;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

#ifdef REG_HAS_NIFTI
    if (this->m_TemplateImage.data != NULL) {
        nifti_image_free(this->m_TemplateImage.data);
        this->m_TemplateImage.data = NULL;
    }
#endif
    this->m_TemplateImage.read = true;

    ierr = this->Read(x, filenames); CHKERRQ(ierr);

    this->m_TemplateImage.read = false;

    ierr = VecMin(*x, NULL, &minval); CHKERRQ(ierr);
    ierr = VecMax(*x, NULL, &maxval); CHKERRQ(ierr);
    this->m_TemplateImage.minval = minval;
    this->m_TemplateImage.maxval = maxval;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
PetscErrorCode ReadWriteReg::Read(Vec* x, std::vector< std::string > filenames) {
    PetscErrorCode ierr = 0;
    std::string file, filename;
    IntType nc, nl, ng;
    std::stringstream ss;
    Vec xk = NULL;
    ScalarType *p_x = NULL, *p_xk = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    //ierr = Assert(!filename.empty(), "filename not set"); CHKERRQ(ierr);

    nc = this->m_Opt->m_Domain.nc;
    ierr = Assert(filenames.size() == static_cast<unsigned int>(nc), "size mismatch"); CHKERRQ(ierr);

    for (IntType k = 0; k < nc; ++k) {
        filename = filenames[k];

        // get file name without path
        ierr = GetFileName(file, filename); CHKERRQ(ierr);

        // check if file exists
        ss << "file " << file << " does not exist";
        ierr = Assert(FileExists(filename), ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());

        // display what we are doing
        if (this->m_Opt->m_Verbosity > 2) {
            ss << "reading " << file;
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }

        // read component
        this->m_FileName = filename;
        ierr = this->Read(&xk); CHKERRQ(ierr);

        // display how we are doing
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = ShowValues(xk); CHKERRQ(ierr);
        }

        ierr = Assert(this->m_Opt->m_SetupDone, "error in setup"); CHKERRQ(ierr);
        nl = this->m_Opt->m_Domain.nl;
        ng = this->m_Opt->m_Domain.ng;

        if (*x == NULL) {
            ierr = VecCreate(PETSC_COMM_WORLD, x); CHKERRQ(ierr);
            ierr = VecSetSizes(*(x), nl*nc, ng*nc); CHKERRQ(ierr);
            ierr = VecSetType(*(x), VECSTANDARD); CHKERRQ(ierr);
            //ierr = VecCreate(*x, nc*nl, nc*ng); CHKERRQ(ierr);
        }

        ierr = VecGetArray(*x, &p_x); CHKERRQ(ierr);
        // extract individual components
        ierr = VecGetArray(xk, &p_xk); CHKERRQ(ierr);
        try {std::copy(p_xk, p_xk+nl, p_x+k*nl);}
        catch (std::exception& err) {
            ierr = ThrowError(err); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(xk, &p_xk); CHKERRQ(ierr);
        ierr = VecRestoreArray(*x, &p_x); CHKERRQ(ierr);

        // delete temporary variable
        if (xk != NULL) {
            ierr = VecDestroy(&xk); CHKERRQ(ierr); xk = NULL;
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
PetscErrorCode ReadWriteReg::Read(Vec* x, std::string filename) {
    PetscErrorCode ierr = 0;
    std::string file, msg;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(!filename.empty(), "filename not set"); CHKERRQ(ierr);

    // get file name without path
    ierr = GetFileName(file, filename); CHKERRQ(ierr);

    // check if file exists
    msg = "file " + file + " does not exist";
    ierr = Assert(FileExists(filename), msg); CHKERRQ(ierr);

    // display what we are doing
    if (this->m_Opt->m_Verbosity > 2) {
        msg = "reading " + file;
        ierr = DbgMsg(msg); CHKERRQ(ierr);
    }
    this->m_FileName = filename;
    ierr = this->Read(x); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
PetscErrorCode ReadWriteReg::Read(Vec* x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    ZeitGeist_define(IO_READ);
    ZeitGeist_tick(IO_READ);

    ierr = Assert(!this->m_FileName.empty(), "filename not set"); CHKERRQ(ierr);

    if (this->m_FileName.find(".nii") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->ReadNII(x); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (this->m_FileName.find(".nii.gz") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->ReadNII(x); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (this->m_FileName.find(".hdr") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->ReadNII(x); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (this->m_FileName.find(".bin") != std::string::npos) {
        ierr = this->ReadBIN(x); CHKERRQ(ierr);
    } else if (this->m_FileName.find(".nc") != std::string::npos) {
#ifdef REG_HAS_PNETCDF
        ierr = this->ReadNC(x); CHKERRQ(ierr);
#else
        ierr = ThrowError("install pnetcdf library/enable pnetcdf support"); CHKERRQ(ierr);
#endif
    } else {
        ierr = ThrowError("could not read: data type not supported"); CHKERRQ(ierr);
    }
    
    ZeitGeist_tock(IO_READ);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
PetscErrorCode ReadWriteReg::Read(VecField* v, std::string fnx1,
                                               std::string fnx2,
                                               std::string fnx3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(!fnx1.empty(), "filename not set"); CHKERRQ(ierr);
    ierr = Assert(!fnx2.empty(), "filename not set"); CHKERRQ(ierr);
    ierr = Assert(!fnx3.empty(), "filename not set"); CHKERRQ(ierr);

    ierr = this->Read(&v->m_X1, fnx1); CHKERRQ(ierr);
    ierr = this->Read(&v->m_X2, fnx2); CHKERRQ(ierr);
    ierr = this->Read(&v->m_X3, fnx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write template image to file
 *******************************************************************/
PetscErrorCode ReadWriteReg::WriteR(Vec x, std::string filename, bool multicomponent) {
    PetscErrorCode ierr = 0;
    ScalarType maxval, minval;
    IntType nc;
    std::stringstream ss;
    bool rescaled = false;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nc = this->m_Opt->m_Domain.nc;

    this->m_ReferenceImage.write = true;
    if (this->m_Opt->m_RegFlags.applyrescaling) {
        if (this->m_ReferenceImage.minval != -1.0 && this->m_ReferenceImage.maxval != -1.0) {
            minval = this->m_ReferenceImage.minval;
            maxval = this->m_ReferenceImage.maxval;
            if (this->m_Opt->m_Verbosity > 1) {
                ss << "rescaling output [0,1] -> [" << minval << "," << maxval << "]";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            }
            if (multicomponent) {
                ierr = Rescale(x, minval, maxval, nc); CHKERRQ(ierr);
            } else {
                ierr = Rescale(x, minval, maxval); CHKERRQ(ierr);
            }
            rescaled = true;
        }
    }

    ierr = this->Write(x, filename, multicomponent); CHKERRQ(ierr);

    this->m_ReferenceImage.write = false;
    if (rescaled) {
        if (multicomponent) {
            ierr = Normalize(x, nc); CHKERRQ(ierr);
        } else {
            ierr = Normalize(x); CHKERRQ(ierr);
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write template image to file
 *******************************************************************/
PetscErrorCode ReadWriteReg::WriteT(Vec x, std::string filename, bool multicomponent) {
    PetscErrorCode ierr = 0;
    ScalarType maxval, minval;
    IntType nc;
    std::stringstream ss;
    bool rescaled = false;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);


    this->m_TemplateImage.write = true;

    if (this->m_Opt->m_RegFlags.applyrescaling) {
        if (this->m_TemplateImage.minval != -1.0 && this->m_TemplateImage.maxval != -1.0) {
            nc = this->m_Opt->m_Domain.nc;
            minval = this->m_TemplateImage.minval;
            maxval = this->m_TemplateImage.maxval;
            if (this->m_Opt->m_Verbosity > 1) {
                ss << "rescaling output [0,1] -> [" << minval << "," << maxval << "]";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            }

            if (multicomponent) {
                ierr = Rescale(x, minval, maxval, nc); CHKERRQ(ierr);
            } else {
                ierr = Rescale(x, minval, maxval); CHKERRQ(ierr);
            }
            rescaled = true;
        }
    }

    ierr = this->Write(x, filename, multicomponent); CHKERRQ(ierr);

    this->m_TemplateImage.write = false;
    if (rescaled) {
        if (multicomponent) {
            ierr = Normalize(x, nc); CHKERRQ(ierr);
        } else {
            ierr = Normalize(x); CHKERRQ(ierr);
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
PetscErrorCode ReadWriteReg::Write(Vec x, std::string filename, bool multicomponent) {
    PetscErrorCode ierr = 0;
    IntType nc, nl, ng;
    Vec xk = NULL;
    std::string msg, path, file, ext;
    std::stringstream ss;
    ScalarType *p_xk = NULL, *p_x = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(!filename.empty(), "filename not set"); CHKERRQ(ierr);
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);

    // get file name without path
    ierr = GetFileName(file, filename); CHKERRQ(ierr);

    // display what we are doing
    if (this->m_Opt->m_Verbosity > 2) {
        msg = "writing " + file;
        ierr = DbgMsg(msg); CHKERRQ(ierr);
        ierr = ShowValues(x); CHKERRQ(ierr);
    }

    if (multicomponent == false) {
        this->m_FileName = this->m_Opt->m_FileNames.xfolder + filename;
        ierr = this->Write(x); CHKERRQ(ierr);
    } else {
        ierr = GetFileName(path, file, ext, filename); CHKERRQ(ierr);
        if (path.empty()) {
            filename = file;
        } else {
            filename = path + "/" + file;
        }

        nc = this->m_Opt->m_Domain.nc;
        nl = this->m_Opt->m_Domain.nl;
        ng = this->m_Opt->m_Domain.ng;

        // allocate data
        //ierr = VecCreate(xk, nl, ng); CHKERRQ(ierr);
        ierr = VecCreate(PETSC_COMM_WORLD, &xk); CHKERRQ(ierr);
        ierr = VecSetSizes(xk, nl, ng); CHKERRQ(ierr);
        ierr = VecSetType(xk, VECSTANDARD); CHKERRQ(ierr);

        ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {
            if (this->m_Opt->m_Verbosity > 2) {
                //msg = "writing component " + std::to_string(static_cast<int>(k));
                msg = "writing component " + std::to_string(static_cast<long long>(k));
                ierr = DbgMsg(msg); CHKERRQ(ierr);
            }

            // extract individual components
            ierr = VecGetArray(xk, &p_xk); CHKERRQ(ierr);
            try {std::copy(p_x+k*nl, p_x+(k+1)*nl, p_xk);}
            catch(std::exception&) {
                ierr = ThrowError("copy failed"); CHKERRQ(ierr);
            }
            ierr = VecRestoreArray(xk, &p_xk); CHKERRQ(ierr);

            // construct file name and write out component
            ss  << filename << "-" << std::setw(3) << std::setfill('0') << k << ext;
            this->m_FileName = this->m_Opt->m_FileNames.xfolder + ss.str();

            ierr = this->Write(xk); CHKERRQ(ierr);

            ss.str(std::string()); ss.clear();
        }
        ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);
    }

    if (xk != NULL) {
        ierr = VecDestroy(&xk); CHKERRQ(ierr); xk = NULL;
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
PetscErrorCode ReadWriteReg::Write(Vec x) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    ZeitGeist_define(IO_WRITE);
    ZeitGeist_tick(IO_WRITE);

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(!this->m_FileName.empty(), "filename not set"); CHKERRQ(ierr);

    if (this->m_FileName.find(".nii") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->WriteNII(x); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (this->m_FileName.find(".nii.gz") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->WriteNII(x); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (this->m_FileName.find(".hdr") != std::string::npos) {
#ifdef REG_HAS_NIFTI
        ierr = this->WriteNII(x); CHKERRQ(ierr);
#else
        ierr = ThrowError("install nifit library/enable nifti support"); CHKERRQ(ierr);
#endif
    } else if (this->m_FileName.find(".bin") != std::string::npos) {
        ierr = this->WriteBIN(x); CHKERRQ(ierr);
    } else if (this->m_FileName.find(".nc") != std::string::npos) {
#ifdef REG_HAS_PNETCDF
        ierr = this->WriteNC(x); CHKERRQ(ierr);
#else
        ierr = ThrowError("install pnetcdf library/enable pnetcdf support"); CHKERRQ(ierr);
#endif
    } else {
        ierr = ThrowError("could not write: data type not supported"); CHKERRQ(ierr);
    }
    
    ZeitGeist_tock(IO_WRITE);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
PetscErrorCode ReadWriteReg::Write(VecField* v, std::string filename) {
    PetscErrorCode ierr = 0;
    std::string fnx1, fnx2, fnx3, path, file, ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(!filename.empty(), "filename not set"); CHKERRQ(ierr);

    ierr = GetFileName(path, file, ext, filename); CHKERRQ(ierr);
    if (path.empty()) {
        fnx1 = file + "-x1" + ext;
        fnx2 = file + "-x2" + ext;
        fnx3 = file + "-x3" + ext;
    } else {
        fnx1 = path + "/" + file + "-x1" + ext;
        fnx2 = path + "/" + file + "-x2" + ext;
        fnx3 = path + "/" + file + "-x3" + ext;
    }

    ierr = this->Write(v->m_X1, fnx1); CHKERRQ(ierr);
    ierr = this->Write(v->m_X2, fnx2); CHKERRQ(ierr);
    ierr = this->Write(v->m_X3, fnx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief read nifty image
 *******************************************************************/
#ifdef REG_HAS_NIFTI
PetscErrorCode ReadWriteReg::ReadNII(Vec* x) {
    PetscErrorCode ierr = 0;
    std::string file;
    std::stringstream ss;
    int rank, rval;
    IntType ng, nl, nglobal, nx[3];
    ScalarType *p_x = NULL;
    nifti_image *image = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // get file name without path
    ierr = GetFileName(file, this->m_FileName); CHKERRQ(ierr);

    // read header file
    image = nifti_image_read(this->m_FileName.c_str(), false);
    ss << "could not read image " + file;
    ierr = Assert(image != NULL, ss.str()); CHKERRQ(ierr);
    ss.clear(); ss.str(std::string());

    // get number of grid points
//    nx[0] = static_cast<IntType>(image->nx);
//    nx[1] = static_cast<IntType>(image->ny);
//    nx[2] = static_cast<IntType>(image->nz);
    nx[2] = static_cast<IntType>(image->nx);
    nx[1] = static_cast<IntType>(image->ny);
    nx[0] = static_cast<IntType>(image->nz);

    if (this->m_ReferenceImage.read) {
        if (this->m_Opt->m_Verbosity > 2) {
            ss << "reading reference image " + file;
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
        if (this->m_ReferenceImage.data == NULL) {
            this->m_ReferenceImage.data = image;
            this->m_ReferenceImage.nx[0] = nx[0];
            this->m_ReferenceImage.nx[1] = nx[1];
            this->m_ReferenceImage.nx[2] = nx[2];
        } else {
            ierr = Assert(this->m_ReferenceImage.nx[0] == nx[0], "dimension mismatch"); CHKERRQ(ierr);
            ierr = Assert(this->m_ReferenceImage.nx[1] == nx[1], "dimension mismatch"); CHKERRQ(ierr);
            ierr = Assert(this->m_ReferenceImage.nx[2] == nx[2], "dimension mismatch"); CHKERRQ(ierr);
        }
    } else if (this->m_TemplateImage.read) {
        if (this->m_Opt->m_Verbosity > 2) {
            ss << "reading template image " + file;
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }

        if (this->m_TemplateImage.data == NULL) {
            this->m_TemplateImage.data = image;
            this->m_TemplateImage.nx[0] = nx[0];
            this->m_TemplateImage.nx[1] = nx[1];
            this->m_TemplateImage.nx[2] = nx[2];
        } else {
            ierr = Assert(this->m_TemplateImage.nx[0] == nx[0], "dimension mismatch"); CHKERRQ(ierr);
            ierr = Assert(this->m_TemplateImage.nx[1] == nx[1], "dimension mismatch"); CHKERRQ(ierr);
            ierr = Assert(this->m_TemplateImage.nx[2] == nx[2], "dimension mismatch"); CHKERRQ(ierr);
        }
    } else {
        // do nothing
    }


    // if we read images, we want to make sure that they have the same size
    if ((this->m_nx[0] == -1) && (this->m_nx[1] == -1) && (this->m_nx[2] == -1)) {
        for (int i = 0; i < 3; ++i) {
            this->m_nx[i] = nx[i];
        }
        if (this->m_Opt->m_Verbosity > 2) {
            ss << "reading image with grid size (" << nx[0] << "," << nx[1] << "," << nx[2] << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
    } else {
        ss << "grid size of input images varies: perform affine registration first";
        for (int i = 0; i < 3; ++i) {
            ierr = Assert(this->m_nx[i] == nx[i], ss.str()); CHKERRQ(ierr);
        }
        ss.clear(); ss.str(std::string());
    }

    // pass number of grid points to options
    for (int i = 0; i < 3; ++i) {
        this->m_Opt->m_Domain.nx[i] = nx[i];
    }

    // do the setup before running the code (this essentially
    // concerns the memory distribution/the setup of accfft
    if (!this->m_Opt->m_SetupDone) {
        ierr = this->m_Opt->DoSetup(); CHKERRQ(ierr);
    }

    // get local size
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    //check global size
    nglobal = 1;
    for (int i = 0; i < 3; ++i) {
        nglobal *= nx[i];
    }
    ierr = Assert(ng == nglobal, "problem in setup"); CHKERRQ(ierr);

    // allocate vector
    //if (*x != NULL) {
    //    ierr = VecDestroy(x); CHKERRQ(ierr); *x = NULL;
    //}
    if (*x == NULL) {
      //ierr = VecCreate(*x, nl, ng); CHKERRQ(ierr);
      ierr = VecCreate(PETSC_COMM_WORLD, x); CHKERRQ(ierr);
      ierr = VecSetSizes(*(x), nl, ng); CHKERRQ(ierr);
      ierr = VecSetType(*(x), VECSTANDARD); CHKERRQ(ierr);
    }

    // compute offset and number of entries to send
    ierr = this->CollectSizes(); CHKERRQ(ierr);

    // read the image data
    if (rank == 0) {
        ierr = this->ReadNII(image); CHKERRQ(ierr);
    }

    ierr = VecGetArray(*x, &p_x); CHKERRQ(ierr);
    rval = MPI_Scatterv(this->m_Data, this->m_nSend, this->m_nOffset, MPIU_SCALAR, p_x, nl, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
    ierr = MPIERRQ(rval); CHKERRQ(ierr);
    ierr = VecRestoreArray(*x, &p_x); CHKERRQ(ierr);

    if (!this->m_ReferenceImage.read && !this->m_TemplateImage.read) {
        if (image != NULL) {
            nifti_image_free(image); image = NULL;
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief read nifty image with right component type
 *******************************************************************/
#ifdef REG_HAS_NIFTI
PetscErrorCode ReadWriteReg::ReadNII(nifti_image* image) {
    PetscErrorCode ierr;
    DataType datatype = DOUBLE;
    std::string msg;
    int rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    msg = "should only be called on master/root rank";
    ierr = Assert(rank == 0, msg); CHKERRQ(ierr);

    switch (image->datatype) {
        case NIFTI_TYPE_UINT8:
        {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("reading data of type uint8 (uchar)"); CHKERRQ(ierr);
            }
            datatype = UCHAR;
            ierr = this->ReadNII<unsigned char>(image); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT8:
        {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("reading data of type int8 (char)"); CHKERRQ(ierr);
            }
            datatype = CHAR;
            ierr = this->ReadNII<char>(image); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT16:
        {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("reading data of type uint16 (unsigned short)"); CHKERRQ(ierr);
            }
            datatype = USHORT;
            ierr = this->ReadNII<unsigned short>(image); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT16:
        {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("reading data of type int16 (short)"); CHKERRQ(ierr);
            }
            datatype = SHORT;
            ierr = this->ReadNII<short>(image); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT32:
        {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("reading data of type uint32 (unsigned int)"); CHKERRQ(ierr);
            }
            datatype = UINT;
            ierr = this->ReadNII<unsigned int>(image); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT32:
        {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("reading data of type int32 (int)"); CHKERRQ(ierr);
            }
            datatype = INT;
            ierr = this->ReadNII<int>(image); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT32:
        {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("reading data of type float32 (float)"); CHKERRQ(ierr);
            }
            datatype = FLOAT;
            ierr = this->ReadNII<float>(image); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT64:
        {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("reading data of type float64 (double)"); CHKERRQ(ierr);
            }
            datatype = DOUBLE;
            ierr = this->ReadNII<double>(image); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("image data not supported"); CHKERRQ(ierr);
            break;
        }
    }

    // if we read the reference image and the template
    // image we have to remember the data type
    if (this->m_ReferenceImage.read) {
        this->m_ReferenceImage.datatype = datatype;
    }

    if (this->m_TemplateImage.read) {
        this->m_TemplateImage.datatype = datatype;
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief get component type of NII images
 *******************************************************************/
#ifdef REG_HAS_NIFTI
template <typename T> PetscErrorCode ReadWriteReg::ReadNII(nifti_image* image) {
    PetscErrorCode ierr = 0;
    T *data = NULL;
    std::string msg;
    IntType ng, nx[3];
    int rank, master = 0;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    msg = "should only be called on master/root rank";
    ierr = Assert(rank == master, msg); CHKERRQ(ierr);

    // get number of grid points
    ng = this->m_Opt->m_Domain.ng;

    // allocate data buffer
    if (this->m_Data == NULL) {
        try {this->m_Data = new ScalarType[ng];}
        catch (std::bad_alloc&) {
            ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // load the image
    if (nifti_image_load(image) == -1) {
        msg = "could not read image " + this->m_FileName;
        ierr = ThrowError(msg); CHKERRQ(ierr);
    }

    // assign data
    data = static_cast<T*>(image->data);
    ierr = Assert(data != NULL, "null pointer"); CHKERRQ(ierr);

    // get global number of points
    ng = this->m_Opt->m_Domain.ng;
    nx[0] = this->m_Opt->m_Domain.nx[0];
    nx[1] = this->m_Opt->m_Domain.nx[1];
    nx[2] = this->m_Opt->m_Domain.nx[2];

    IntType k = 0;
    for (int p = 0; p < this->m_NumProcs; ++p) {
        for (IntType i1 = 0; i1 < this->m_iSizeC[3*p+0]; ++i1) {  // x1
            for (IntType i2 = 0; i2 < this->m_iSizeC[3*p+1]; ++i2) {  // x2
                for (IntType i3 = 0; i3 < this->m_iSizeC[3*p+2]; ++i3) {  // x3
                    IntType j1 = i1 + this->m_iStartC[3*p+0];
                    IntType j2 = i2 + this->m_iStartC[3*p+1];
                    IntType j3 = i3 + this->m_iStartC[3*p+2];
                    IntType l = GetLinearIndex(j1, j2, j3, nx);
                    this->m_Data[k++] = static_cast<ScalarType>(data[l]);
                }  // for i1
            }  // for i2
        }  // for i3
    }  // for all procs

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief write buffer to nii files
 *******************************************************************/
#ifdef REG_HAS_NIFTI
PetscErrorCode ReadWriteReg::WriteNII(Vec x) {
    PetscErrorCode ierr;
    nifti_image* image = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // if we write the template or refernence image to file,
    // we'll use the data we have read in; if this data has not
    // been read, we'll pass a NULL pointer resulting in an image
    // being allocated
    if (this->m_TemplateImage.write) {
        image = this->m_TemplateImage.data;
    } else if (this->m_ReferenceImage.write) {
        image = this->m_ReferenceImage.data;
    } else {
        // now we have to deal with the case that a reference
        // image might have been read, so we want to use
        // all the information within its header (including the
        // orientation), but write out the data using the
        // datatype we have used for our computations
        if (this->m_ImageData == NULL) {   // only do this once
            // if we have read the reference image
            // note: the following will only be true
            // on the master rank, because we allocate
            // image data only there
            if (this->m_ReferenceImage.data != NULL) {
                this->m_ImageData = nifti_copy_nim_info(this->m_ReferenceImage.data);
            } else {
                if (this->m_TemplateImage.data != NULL) {
                    this->m_ImageData = nifti_copy_nim_info(this->m_TemplateImage.data);
                }
            }
            if (this->m_ImageData != NULL) {
                // switch precision in io
#if defined(PETSC_USE_REAL_SINGLE)
                this->m_ImageData->datatype = NIFTI_TYPE_FLOAT32; // single precision
#else
                this->m_ImageData->datatype = NIFTI_TYPE_FLOAT64; // double precision
#endif
                this->m_ImageData->nbyper = sizeof(ScalarType);
                // allocate image buffer
                try {this->m_ImageData->data = new ScalarType[this->m_ImageData->nvox];}
                catch (std::bad_alloc&) {
                    ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
                }
            }
        }
        image = this->m_ImageData;
    }

    // if nifty image is null pointer default to double
    if (image == NULL) {
        ierr = this->WriteNII<ScalarType>(&image, x); CHKERRQ(ierr);
    } else {
        switch (image->datatype) {
            case NIFTI_TYPE_UINT8:
            {
                ierr = this->WriteNII<unsigned char>(&image, x); CHKERRQ(ierr);
                break;
            }
            case NIFTI_TYPE_INT8:
            {
                ierr = this->WriteNII<char>(&image, x); CHKERRQ(ierr);
                break;
            }
            case NIFTI_TYPE_UINT16:
            {
                ierr = this->WriteNII<unsigned short>(&image, x); CHKERRQ(ierr);
                break;
            }
            case NIFTI_TYPE_INT16:
            {
                ierr = this->WriteNII<short>(&image, x); CHKERRQ(ierr);
                break;
            }
            case NIFTI_TYPE_UINT32:
            {
                ierr = this->WriteNII<unsigned int>(&image, x); CHKERRQ(ierr);
                break;
            }
            case NIFTI_TYPE_INT32:
            {
                ierr = this->WriteNII<int>(&image, x); CHKERRQ(ierr);
                break;
            }
            case NIFTI_TYPE_FLOAT32:
            {
                ierr = this->WriteNII<float>(&image, x); CHKERRQ(ierr);
                break;
            }
            case NIFTI_TYPE_FLOAT64:
            {
                ierr = this->WriteNII<double>(&image, x); CHKERRQ(ierr);
                break;
            }
            default:
            {
                ierr = ThrowError("image data not supported"); CHKERRQ(ierr);
                break;
            }
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief write buffer to nii files
 *******************************************************************/
#ifdef REG_HAS_NIFTI
template <typename T>
PetscErrorCode ReadWriteReg::WriteNII(nifti_image** image, Vec x) {
    PetscErrorCode ierr;
    T* data = NULL;
    ScalarType *p_xc = NULL;
    int nprocs, rank, rval, master = 0;
    IntType nx[3], ng, nl;
    bool deleteimage = false;
    std::string msg;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // get number of ranks
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    ng = this->m_Opt->m_Domain.ng;
    nl = this->m_Opt->m_Domain.nl;

    // allocate the index buffers on master rank
    if (rank == master) {
        // we need to allocate the image if it's a zero pointer; this
        // will also create a standard header file; not tested (might need
        // to parse the dimensions of the data)
        if ((*image) == NULL) {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("allocating nifti image buffer"); CHKERRQ(ierr);
            }
            ierr = this->AllocateImage(image, x); CHKERRQ(ierr);
            deleteimage = true;
        }
        ierr = Assert((*image) != NULL, "null pointer"); CHKERRQ(ierr);

        // construct file name
        std::string file(this->m_FileName);

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

    // allocate data buffer
    if (this->m_Data == NULL) {
        try {this->m_Data = new ScalarType[ng];}
        catch (std::bad_alloc&) {
            ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // collect sizes and compute number of data to send
    ierr = this->CollectSizes(); CHKERRQ(ierr);

    // gather data on master rank
    ierr = VecGetArray(x, &p_xc); CHKERRQ(ierr);
    rval = MPI_Gatherv(p_xc, nl, MPIU_SCALAR, this->m_Data, this->m_nSend, this->m_nOffset, MPIU_SCALAR, master, PETSC_COMM_WORLD);
    ierr = MPIERRQ(rval); CHKERRQ(ierr);
    ierr = VecRestoreArray(x, &p_xc); CHKERRQ(ierr);

    nx[0] = this->m_Opt->m_Domain.nx[0];
    nx[1] = this->m_Opt->m_Domain.nx[1];
    nx[2] = this->m_Opt->m_Domain.nx[2];

    if (rank == master) {
        // cast pointer of nifti image data
        data = reinterpret_cast<T*>((*image)->data);

        IntType k = 0;
        for (int p = 0; p < nprocs; ++p) {
            for (IntType i1 = 0; i1 < this->m_iSizeC[3*p+0]; ++i1) {  // x1
                for (IntType i2 = 0; i2 < this->m_iSizeC[3*p+1]; ++i2) {  // x2
                    for (IntType i3 = 0; i3 < this->m_iSizeC[3*p+2]; ++i3) {  // x3
                        IntType j1 = i1 + this->m_iStartC[3*p+0];
                        IntType j2 = i2 + this->m_iStartC[3*p+1];
                        IntType j3 = i3 + this->m_iStartC[3*p+2];
                        IntType l = GetLinearIndex(j1, j2, j3, nx);
                        data[l] = static_cast<T>(this->m_Data[k++]);
                    }  // for i1
                }  // for i2
            }  // for i3
        }  // for all procs

        // write image to file
        nifti_image_write(*image);
    }  // if on master


    if (deleteimage) {
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = DbgMsg("deleting nifti image buffer"); CHKERRQ(ierr);
        }
        nifti_image_free(*image); *image = NULL;
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}
#endif



/********************************************************************
 * @brief allocate buffer for nifty image
 *******************************************************************/
#ifdef REG_HAS_NIFTI
PetscErrorCode ReadWriteReg::AllocateImage(nifti_image** image, Vec x) {
    PetscErrorCode ierr = 0;
    IntType n, nl;
    int rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    nl = this->m_Opt->m_Domain.nl;

    // init nifty image
    *image = nifti_simple_init_nim();

    // dimensionalty of data: default is 5 (space, time, components)
    (*image)->dim[0] = (*image)->ndim = 5;

    ierr = VecGetLocalSize(x, &n); CHKERRQ(ierr);
    (*image)->dim[1] = (*image)->nx = this->m_Opt->m_Domain.nx[2];
    (*image)->dim[2] = (*image)->ny = this->m_Opt->m_Domain.nx[1];
    (*image)->dim[3] = (*image)->nz = this->m_Opt->m_Domain.nx[0];

//    (*image)->dim[1] = (*image)->nx = this->m_Opt->m_Domain.nx[0];
//    (*image)->dim[2] = (*image)->ny = this->m_Opt->m_Domain.nx[1];
//    (*image)->dim[3] = (*image)->nz = this->m_Opt->m_Domain.nx[2];

    (*image)->pixdim[1] = static_cast<float>(this->m_Opt->m_Domain.hx[0]);  // x direction
    (*image)->pixdim[2] = static_cast<float>(this->m_Opt->m_Domain.hx[1]);  // y direction
    (*image)->pixdim[3] = static_cast<float>(this->m_Opt->m_Domain.hx[2]);  // z direction

    // TODO: add temporal support
    if (n == nl) {  // scalar field
        (*image)->dim[4] = (*image)->nt = 1;
        (*image)->dim[5] = (*image)->nu = 1;

        // temporal step size
        (*image)->pixdim[4] = 1.0;
    } else if (n == 2*nl) {  // 2D vector field
        (*image)->dim[4] = (*image)->nt = 1;
        (*image)->dim[5] = (*image)->nu = 2;

        // temporal step size
        (*image)->pixdim[4] = 1.0;

        // step size (vector field)
        (*image)->pixdim[5] = (*image)->du = static_cast<float>(this->m_Opt->m_Domain.hx[0]);
        (*image)->pixdim[6] = (*image)->dv = static_cast<float>(this->m_Opt->m_Domain.hx[1]);
    } else if (n == 3*nl) {  // 3D vector field
        (*image)->dim[4] = (*image)->nt = 1;
        (*image)->dim[5] = (*image)->nu = 3;

        // temporal step size
        (*image)->pixdim[4] = 1.0;

        // step size (vector field)
        (*image)->pixdim[5] = (*image)->du = static_cast<float>(this->m_Opt->m_Domain.hx[0]);
        (*image)->pixdim[6] = (*image)->dv = static_cast<float>(this->m_Opt->m_Domain.hx[1]);
        (*image)->pixdim[7] = (*image)->dw = static_cast<float>(this->m_Opt->m_Domain.hx[2]);
    }

    // switch precision in io
#if defined(PETSC_USE_REAL_SINGLE)
    (*image)->datatype = NIFTI_TYPE_FLOAT32; // single precision
#else
    (*image)->datatype = NIFTI_TYPE_FLOAT64; // double precision
#endif
    (*image)->nbyper = sizeof(ScalarType);
    // compute number of voxels (space and time)
    (*image)->nvox = 1;
    for (int i = 1; i <= 4; ++i) {
        (*image)->nvox *= (*image)->dim[i];
    }
    (*image)->nvox *= (*image)->nu;

    // only allocate image buffer on root
    if (rank == 0) {
        try {(*image)->data = new ScalarType[(*image)->nvox];}
        catch (std::bad_alloc&) {
            ierr = ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}
#endif




/********************************************************************
 * @brief write binary data set to file
 *******************************************************************/
PetscErrorCode ReadWriteReg::ReadBIN(Vec* x) {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    PetscViewer viewer = NULL;
    PetscFunctionBegin;

    //if (*x != NULL) {
    //   ierr = VecDestroy(x); CHKERRQ(ierr);
    //    *x = NULL;
    //}

    if (!this->m_Opt->m_SetupDone) {
        ierr = this->m_Opt->DoSetup(); CHKERRQ(ierr);
    }
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    if (*x == NULL) {
      //ierr = VecCreate(*x, nl, ng); CHKERRQ(ierr);
      ierr = VecCreate(PETSC_COMM_WORLD, x); CHKERRQ(ierr);
      ierr = VecSetSizes(*(x), nl, ng); CHKERRQ(ierr);
      ierr = VecSetType(*(x), VECSTANDARD); CHKERRQ(ierr);
    }

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, this->m_FileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = Assert(viewer != NULL, "could not read binary file"); CHKERRQ(ierr);
    ierr = PetscViewerBinarySetFlowControl(viewer, 2); CHKERRQ(ierr);
    ierr = VecLoad(*x, viewer); CHKERRQ(ierr);

    // clean up
    if (viewer != NULL) {
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        viewer = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write bin to file
 *******************************************************************/
PetscErrorCode ReadWriteReg::WriteBIN(Vec x) {
    PetscErrorCode ierr = 0;
    PetscViewer viewer = NULL;
    PetscFunctionBegin;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, this->m_FileName.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
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
 * @brief write netcdf to file
 *******************************************************************/
#ifdef REG_HAS_PNETCDF
PetscErrorCode ReadWriteReg::ReadNC(Vec* x) {
    PetscErrorCode ierr = 0;
    int rank, ncerr, fileid, ndims, nvars, ngatts, unlimited, varid[1];
    IntType nl, ng;
    ScalarType *p_x = NULL;
    MPI_Offset istart[3], isize[3];
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //if ((*x) != NULL) {ierr = VecDestroy(x); CHKERRQ(ierr); *x = NULL;}

    if (!this->m_Opt->m_SetupDone) {
        ierr = this->m_Opt->DoSetup(); CHKERRQ(ierr);
    }
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    if (*x == NULL) {
      //ierr = VecCreate(*x, nl, ng); CHKERRQ(ierr);
      ierr = VecCreate(PETSC_COMM_WORLD, x); CHKERRQ(ierr);
      ierr = VecSetSizes(*(x), nl, ng); CHKERRQ(ierr);
      ierr = VecSetType(*(x), VECSTANDARD); CHKERRQ(ierr);
    }

    // open file
    ncerr = ncmpi_open(PETSC_COMM_WORLD,this->m_FileName.c_str(), NC_NOWRITE, MPI_INFO_NULL, &fileid);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    // query info about field named "data"
    ncerr = ncmpi_inq(fileid, &ndims, &nvars, &ngatts, &unlimited);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ncerr = ncmpi_inq_varid(fileid, "data", &varid[0]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    istart[0] = this->m_Opt->m_Domain.istart[0];
    istart[1] = this->m_Opt->m_Domain.istart[1];
    istart[2] = this->m_Opt->m_Domain.istart[2];

    isize[0] = this->m_Opt->m_Domain.isize[0];
    isize[1] = this->m_Opt->m_Domain.isize[1];
    isize[2] = this->m_Opt->m_Domain.isize[2];

    ierr = VecGetArray(*x, &p_x); CHKERRQ(ierr);
    ncerr = ncmpi_get_vara_all(fileid, varid[0], istart, isize, p_x, nl, MPIU_SCALAR);
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
PetscErrorCode ReadWriteReg::WriteNC(Vec x) {
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

    c_comm = this->m_Opt->m_Domain.mpicomm;

    // create netcdf file
    ncerr = ncmpi_create(c_comm, this->m_FileName.c_str(), mode, MPI_INFO_NULL, &fileid);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    nx[0] = static_cast<int>(this->m_Opt->m_Domain.nx[0]);
    nx[1] = static_cast<int>(this->m_Opt->m_Domain.nx[1]);
    nx[2] = static_cast<int>(this->m_Opt->m_Domain.nx[2]);
    nl = static_cast<int>(this->m_Opt->m_Domain.nl);

    // set size
    ncerr = ncmpi_def_dim(fileid, "x", nx[0], &dims[0]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ncerr = ncmpi_def_dim(fileid, "y", nx[1], &dims[1]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ncerr = ncmpi_def_dim(fileid, "z", nx[2], &dims[2]);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    // define name for output field
#if defined(PETSC_USE_REAL_SINGLE)
    ncerr = ncmpi_def_var(fileid, "data", NC_FLOAT, 3, dims, &varid[0]);
#else
    ncerr = ncmpi_def_var(fileid, "data", NC_DOUBLE, 3, dims, &varid[0]);
#endif
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    iscdf5 = usecdf5 ? 1 : 0;
    ncerr = ncmpi_put_att_int(fileid, NC_GLOBAL, "CDF-5 mode", NC_INT, 1, &iscdf5);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);
    ncerr = ncmpi_enddef(fileid);
    ierr = NCERRQ(ncerr); CHKERRQ(ierr);

    // get local sizes
    istart[0] = static_cast<MPI_Offset>(this->m_Opt->m_Domain.istart[0]);
    istart[1] = static_cast<MPI_Offset>(this->m_Opt->m_Domain.istart[1]);
    istart[2] = static_cast<MPI_Offset>(this->m_Opt->m_Domain.istart[2]);

    isize[0] = static_cast<MPI_Offset>(this->m_Opt->m_Domain.isize[0]);
    isize[1] = static_cast<MPI_Offset>(this->m_Opt->m_Domain.isize[1]);
    isize[2] = static_cast<MPI_Offset>(this->m_Opt->m_Domain.isize[2]);

    ierr = Assert(nl == isize[0]*isize[1]*isize[2], "size error"); CHKERRQ(ierr);

    // write data to file
    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
    ncerr = ncmpi_put_vara_all(fileid, varid[0], istart, isize, p_x, nl, MPIU_SCALAR);
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
