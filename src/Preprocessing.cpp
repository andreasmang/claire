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

#ifndef _PREPROCESSING_CPP_
#define _PREPROCESSING_CPP_

#include "Preprocessing.hpp"
#include "DifferentiationSM.hpp"
#include <time.h>

namespace reg {

/********************************************************************
 * @brief default constructor
 *******************************************************************/
Preprocessing::Preprocessing() {
    this->Initialize();
}

/********************************************************************
 * @brief constructor
 *******************************************************************/
Preprocessing::Preprocessing(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}

/********************************************************************
 * @brief default deconstructor
 *******************************************************************/
Preprocessing::~Preprocessing() {
    this->ClearMemory();
}

/********************************************************************
 * @brief initialize
 *******************************************************************/
PetscErrorCode Preprocessing::Initialize() {
    PetscErrorCode ierr = 0;
    this->m_Opt = nullptr;

    this->m_GridChangeOpsSet = false;
    this->m_ResetGridChangeOps = false;
    this->m_IndicesCommunicated = false;
    this->m_GridChangeIndicesComputed = false;

    this->m_ReadWrite = nullptr;

    //this->m_XHatFine = nullptr;
    //this->m_XHatCoarse = nullptr;
    this->m_FFTFinePlan = nullptr;
    this->m_FFTCoarsePlan = nullptr;

    this->m_FourierCoeffSendF = nullptr;
    this->m_FourierCoeffSendC = nullptr;

    this->m_FourierCoeffRecvF = nullptr;
    this->m_FourierCoeffRecvC = nullptr;

    this->m_FourierIndicesRecvF = nullptr;
    this->m_FourierIndicesRecvC = nullptr;

    this->m_FourierIndicesSendF = nullptr;
    this->m_FourierIndicesSendC = nullptr;

    this->m_NumSend = nullptr;
    this->m_NumRecv = nullptr;

    this->m_OffsetSend = nullptr;
    this->m_OffsetRecv = nullptr;

    this->m_nAllocSend = 0;
    this->m_nAllocRecv = 0;

    this->m_SendRequest = nullptr;
    this->m_RecvRequest = nullptr;

    this->m_OverlapMeasures = nullptr;
//    this->m_LabelValues = nullptr;
//    this->m_NoLabel = -99;

    //this->m_xhat = nullptr;
    //this->m_yhat = nullptr;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clear memory
 *******************************************************************/
PetscErrorCode Preprocessing::ClearMemory() {
    PetscErrorCode ierr = 0;

    /*if (this->m_xhat != nullptr) {
        accfft_free(this->m_xhat);
        this->m_xhat = nullptr;
    }
    if (this->m_yhat != nullptr) {
        accfft_free(this->m_yhat);
        this->m_yhat = nullptr;
    }

    if (this->m_XHatFine != nullptr) {
        accfft_free(this->m_XHatFine);
        this->m_XHatFine = nullptr;
    }
    if (this->m_XHatCoarse != nullptr) {
        accfft_free(this->m_XHatCoarse);
        this->m_XHatCoarse = nullptr;
    }*/

    if (this->m_FFTFinePlan != nullptr) {
        accfft_destroy_plan(this->m_FFTFinePlan);
        this->m_FFTFinePlan = nullptr;
    }
    if (this->m_FFTCoarsePlan != nullptr) {
        accfft_destroy_plan(this->m_FFTCoarsePlan);
        this->m_FFTCoarsePlan = nullptr;
    }

    FreeArray(this->m_FourierCoeffSendF);
    FreeArray(this->m_FourierCoeffSendC);
    FreeArray(this->m_FourierCoeffRecvF);
    FreeArray(this->m_FourierCoeffRecvC);
    FreeArray(this->m_FourierIndicesRecvF);
    FreeArray(this->m_FourierIndicesRecvC);
    FreeArray(this->m_FourierIndicesSendF);
    FreeArray(this->m_FourierIndicesSendC);
    FreeArray(this->m_NumSend);
    FreeArray(this->m_NumRecv);
    FreeArray(this->m_OffsetRecv);
    FreeArray(this->m_OffsetSend);
    FreeArray(this->m_SendRequest);
    FreeArray(this->m_RecvRequest);
    FreeArray(this->m_OverlapMeasures);
/*
    if (this->m_LabelValues != nullptr) {
        delete [] this->m_LabelValues;
        this->m_LabelValues = nullptr;
    }
*/
    PetscFunctionReturn(ierr);

}

/********************************************************************
 * @brief set read/write object for data
 *******************************************************************/
PetscErrorCode Preprocessing::SetReadWrite(Preprocessing::ReadWriteType* readwrite) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(readwrite !=  nullptr, "null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite = readwrite;

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief do setup for applying prolongation and restriction
 * operators
 * @param nx_f grid size on fine grid
 * @param nx_c grid size on coarse grid
 *******************************************************************/
PetscErrorCode Preprocessing::SetupGridChangeOps(IntType* nx_f, IntType* nx_c) {
    PetscErrorCode ierr = 0;
    IntType nalloc_c, nalloc_f;
    int _nx_f[3], _ostart_f[3], _osize_f[3], _isize_f[3], _istart_f[3],
        _nx_c[3], _ostart_c[3], _osize_c[3], _isize_c[3], _istart_c[3];
    ScalarType *p_xfd = nullptr, *p_xcd = nullptr;
    ComplexType *p_xfdhat = nullptr, *p_xcdhat = nullptr;
    std::stringstream ss;
    MPI_Comm mpicomm;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Opt->m_Verbosity > 2) {
        ss  << "setup gridchange operator ( (" << nx_c[0]
            << "," << nx_c[1] << "," << nx_c[2]
            << ") <=> (" << nx_f[0] << "," << nx_f[1]
            << "," << nx_f[2] << ") )";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    /*if (this->m_XHatCoarse !=  nullptr) {
        accfft_free(this->m_XHatCoarse);
        this->m_XHatCoarse = nullptr;
    }
    if (this->m_XHatFine !=  nullptr) {
        accfft_free(this->m_XHatFine);
        this->m_XHatFine = nullptr;
    }*/
    if (this->m_FFTFinePlan !=  nullptr) {
        accfft_destroy_plan(this->m_FFTFinePlan);
        this->m_FFTFinePlan = nullptr;
    }
    if (this->m_FFTCoarsePlan !=  nullptr) {
        accfft_destroy_plan(this->m_FFTCoarsePlan);
        this->m_FFTCoarsePlan = nullptr;
    }

    // parse input sizes
    for (int i = 0; i < 3; ++i) {
        _nx_f[i] = static_cast<int>(nx_f[i]);
        _nx_c[i] = static_cast<int>(nx_c[i]);
    }

    this->m_FFTFineScale = 1.0;
    this->m_FFTCoarseScale = 1.0;
    for (int i = 0; i < 3; ++i) {
        this->m_FFTFineScale *= static_cast<ScalarType>(nx_f[i]);
        this->m_FFTCoarseScale *= static_cast<ScalarType>(nx_c[i]);
    }
    this->m_FFTFineScale = 1.0/this->m_FFTFineScale;
    this->m_FFTCoarseScale = 1.0/this->m_FFTCoarseScale;

    // get communicator
    mpicomm = this->m_Opt->m_FFT.mpicomm;

    nalloc_c = accfft_local_size_dft_r2c_t<ScalarType>(_nx_c, _isize_c, _istart_c, _osize_c, _ostart_c, mpicomm);
    nalloc_f = accfft_local_size_dft_r2c_t<ScalarType>(_nx_f, _isize_f, _istart_f, _osize_f, _ostart_f, mpicomm);
    if (this->m_Opt->m_Verbosity > 2) {
        ss << "sizes on coarse grid: isize=("
           << _isize_c[0] << "," << _isize_c[1] << "," << _isize_c[2]
           << ") istart=(" << _istart_c[0] << "," << _istart_c[1]
           << "," << _istart_c[2] << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }
    //ierr = Assert(nalloc_c > 0, "alloc problems 1"); CHKERRQ(ierr);
    //ierr = Assert(nalloc_f > 0, "alloc problems 2"); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i) {
        this->m_osizeC[i] = static_cast<IntType>(_osize_c[i]);
        this->m_osizeF[i] = static_cast<IntType>(_osize_f[i]);
        this->m_ostartC[i] = static_cast<IntType>(_ostart_c[i]);
        this->m_ostartF[i] = static_cast<IntType>(_ostart_f[i]);
    }

    if (this->m_Opt->m_Verbosity > 2) {
        ss  << "coarse: osize=(" << this->m_osizeC[0]
            << "," << this->m_osizeC[1]
            << "," << this->m_osizeC[2]
            << "); ostart=(" << this->m_ostartC[0]
            << "," << this->m_ostartC[1]
            << "," << this->m_ostartC[2]
            << "); n=" << nalloc_c;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
        ss  << "fine: osize=(" << this->m_osizeF[0]
            << "," << this->m_osizeF[1]
            << "," << this->m_osizeF[2]
            << "); ostart=(" << this->m_ostartF[0]
            << "," << this->m_ostartF[1]
            << "," << this->m_ostartF[2]
            << "); n=" << nalloc_f;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    /*if (this->m_XHatCoarse == nullptr) {
        this->m_XHatCoarse = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc_c));
    }
    ierr = Assert(this->m_XHatCoarse != nullptr,"allocation failed"); CHKERRQ(ierr);

    if (this->m_XHatFine == nullptr) {
        this->m_XHatFine = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc_f));
    }
    ierr = Assert(this->m_XHatFine != nullptr,"allocation failed"); CHKERRQ(ierr);
    */
    ierr = this->m_XHatCoarse.Resize(nalloc_c); CHKERRQ(ierr);
    ierr = this->m_XHatFine.Resize(nalloc_f); CHKERRQ(ierr);

    // allocate plan for fine grid
    if (this->m_FFTFinePlan == nullptr) {
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = DbgMsg("initializing fft plan (fine grid)"); CHKERRQ(ierr);
        }

        // TODO: this is not nice! p_xfd and p_xfdhat are only used to tell accfft that it's not working inplace... accfft stores the pointer, but the memory is freed afterwards
        ierr = AllocateMemoryOnce(p_xfd, nalloc_f); CHKERRQ(ierr);
        ierr = AllocateMemoryOnce(p_xfdhat, nalloc_f); CHKERRQ(ierr);

        this->m_FFTFinePlan = accfft_plan_dft_3d_r2c(_nx_f, p_xfd, reinterpret_cast<ScalarType*>(p_xfdhat),
                                                     this->m_Opt->m_FFT.mpicomm, ACCFFT_MEASURE);
        ierr = Assert(this->m_FFTFinePlan != nullptr, "malloc failed"); CHKERRQ(ierr);

        ierr = FreeMemory(p_xfd); CHKERRQ(ierr);
        ierr = FreeMemory(p_xfdhat); CHKERRQ(ierr);
    }

    // allocate plan for coarse grid
    if (this->m_FFTCoarsePlan == nullptr) {
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = DbgMsg("initializing fft plan (coarse grid)"); CHKERRQ(ierr);
        }
        
        ierr = AllocateMemoryOnce(p_xcd, nalloc_c); CHKERRQ(ierr);
        ierr = AllocateMemoryOnce(p_xcdhat, nalloc_c); CHKERRQ(ierr);

        this->m_FFTCoarsePlan = accfft_plan_dft_3d_r2c(_nx_c, p_xcd, reinterpret_cast<ScalarType*>(p_xcdhat),
                                                       this->m_Opt->m_FFT.mpicomm, ACCFFT_MEASURE);
        ierr = Assert(this->m_FFTCoarsePlan != nullptr, "malloc failed"); CHKERRQ(ierr);

        ierr = FreeMemory(p_xcd); CHKERRQ(ierr);
        ierr = FreeMemory(p_xcdhat); CHKERRQ(ierr);
    }

    // set flag
    this->m_GridChangeOpsSet = true;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief convert sharp label image to smooth multi-component image
 * @param labelim label image
 * @param multi-component image
 *******************************************************************/
PetscErrorCode Preprocessing::Labels2MultiCompImage(Vec m, Vec labelmap) {
    PetscErrorCode ierr = 0;
    IntType nl, nc;
    const ScalarType *p_labelmap = nullptr;
    ScalarType *p_m = nullptr;
    int label;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    ierr = Assert(this->m_Opt->m_LabelIDs.size() == static_cast<unsigned int>(nc), "size mismatch"); CHKERRQ(ierr);

    // now assign the individual labels to the
    // individual components
    ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
    ierr = VecGetArrayRead(labelmap, &p_labelmap); CHKERRQ(ierr);
    for (int l = 0; l < nc; ++l) {
        label = this->m_Opt->m_LabelIDs[l];
        for (IntType i = 0; i < nl; ++i) {
            // get current label
            if (label == p_labelmap[i]) {
                p_m[l*nl + i] = 1.0;
            } else {
                p_m[l*nl + i] = 0.0;
            }
        }
    }
    ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(labelmap, &p_labelmap); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/*
PetscErrorCode Preprocessing::Labels2MultiCompImage(Vec m, Vec labelmap) {
    PetscErrorCode ierr = 0;
    IntType nl, nc;
    const ScalarType *p_labelmap = nullptr;
    ScalarType *p_m = nullptr;
    int label, numlabelfound = 0;
    bool labelfound, labelset;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    // TODO: scheme does not work in parallel
    // allocate array for labels (account for background)
    if (this->m_LabelValues == nullptr) {
        try{this->m_LabelValues = new int[nc];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // initialize; we assume the labels are not negative
    for (int l = 0; l < nc; ++l){
        this->m_LabelValues[l] = this->m_NoLabel;
    }

    ierr = VecGetArrayRead(labelmap, &p_labelmap); CHKERRQ(ierr);

    // for all image points
    for (IntType i = 0; i < nl; ++i) {
        // get current label
        label = p_labelmap[i];
        if (label > 0) {
            labelfound = false;
            // figure out which labels exist
            for (int l = 0; l < nc; ++l) {
                if (this->m_LabelValues[l] == label) {
                    labelfound = true;
                }
            }
            // if we have not yet identified this label
            // remember it in the list
            if (labelfound == false) {
                labelset = false;
                for (int l = 0; l < nc; ++l){
                    if (this->m_LabelValues[l] == this->m_NoLabel && labelset == false) {
                        this->m_LabelValues[l] = label;
                        numlabelfound++;
                        labelset = true;
                    }
                }
            }
        }
        if (numlabelfound == nc) {
            // stop serach if we have identified
            // all labels in the label map
            break;
        }
    }

    // now assign the individual labels to the
    // individual components
    ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
    for (int l = 0; l < nc; ++l) {
        label = this->m_LabelValues[l];
        for (IntType i = 0; i < nl; ++i) {
            // get current label
            if (label == p_labelmap[i]) {
                p_m[l*nl + i] = 1.0;
            } else {
                p_m[l*nl + i] = 0.0;
            }
        }
    }
    ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);

    ierr = VecRestoreArrayRead(labelmap, &p_labelmap); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}
*/

PetscErrorCode Preprocessing::EnsurePatitionOfUnity(Vec m) {
    PetscErrorCode ierr = 0;
    IntType nl, nc;
    ScalarType *p_labelprobs = NULL;
    ScalarType *p_m = NULL;
    ScalarType labelsum;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    ierr = Assert(this->m_Opt->m_LabelIDs.size() == static_cast<unsigned int>(nc), "size mismatch"); CHKERRQ(ierr);

    try {p_labelprobs = new ScalarType[nc+1];}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }

    // set dummy values
    ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
    for (IntType i = 0; i < nl; ++i) {

        // compute value for background
        labelsum = 0.0; p_labelprobs[nc] = 1.0;
        for (int l = 0; l < nc; ++l){
            // get label probability
            p_labelprobs[l] = p_m[l*nl + i];

            // compute background
            p_labelprobs[nc] -= p_labelprobs[l];

            // accumulate label probabilities
            labelsum += p_labelprobs[l];
        }

        // set to zero, if negative
        if (p_labelprobs[nc] < 0.0) p_labelprobs[nc] = 0.0;
        labelsum += p_labelprobs[nc];
//        labelsum = labelsum > 1E-1 ? labelsum : 1.0;

        // normalize (partition of unity)
        for (int l = 0; l < nc+1; ++l) {
            p_labelprobs[l] /= labelsum;
            p_m[l*nl + i] = p_labelprobs[l];
        }
    }

    ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);


    if (p_labelprobs != NULL) {delete [] p_labelprobs;}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief convert smooth multi-component image to sharp label image
 * @param labelim label image
 * @param multi-component image
 *******************************************************************/
PetscErrorCode Preprocessing::MultiCompImage2Labels(Vec labelim, Vec m) {
    PetscErrorCode ierr = 0;
    IntType nl, nc;
    int majoritylabel;
    ScalarType *p_labels = nullptr, *p_labelprobs = nullptr;
    const ScalarType *p_m = nullptr;
    ScalarType value, majorityvote, labelsum;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    ierr = Assert(this->m_Opt->m_LabelIDs.size() == static_cast<unsigned int>(nc), "size mismatch"); CHKERRQ(ierr);

    ierr = AllocateArrayOnce(p_labelprobs, nc+1); CHKERRQ(ierr);

    // set dummy values
    ierr = VecGetArrayRead(m, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(labelim, &p_labels); CHKERRQ(ierr);
    for (IntType i = 0; i < nl; ++i) {
        majorityvote  = 0.0;
        majoritylabel = -1;

        // compute value for background
        labelsum = 0.0; p_labelprobs[nc] = 1.0;
        for (int l = 0; l < nc; ++l){
            // get label probability
            p_labelprobs[l] = p_m[l*nl + i];

            // compute background
            p_labelprobs[nc] -= p_labelprobs[l];

            // accumulate label probabilities
            labelsum += p_labelprobs[l];
        }

        // set to zero, if negative
        if (p_labelprobs[nc] < 0.0) p_labelprobs[nc] = 0.0;
        labelsum += p_labelprobs[nc];
//        labelsum = labelsum > 1E-1 ? labelsum : 1.0;

        // normalize (partition of unity)
        for (int l = 0; l < nc+1; ++l) {
            p_labelprobs[l] /= labelsum;
        }

        for (int l = 0; l < nc + 1; ++l) {
            value = p_labelprobs[l];

            // get largest value
            if (value > majorityvote) {
                majoritylabel = l;
                majorityvote  = value;
            }

            // get current label
            if (majoritylabel == nc) {
                p_labels[i] = 0;
            } else {
                p_labels[i] = this->m_Opt->m_LabelIDs[majoritylabel];
            }
        }
    }
    ierr = VecRestoreArray(labelim, &p_labels); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(m, &p_m); CHKERRQ(ierr);

    ierr = FreeArray(p_labelprobs); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}
/*
PetscErrorCode Preprocessing::MultiCompImage2Labels(Vec labelim, Vec m) {
    PetscErrorCode ierr = 0;
    IntType nl, nc;
    int majoritylabel;
    ScalarType *p_labels = nullptr, *p_labelprobs = nullptr;
    const ScalarType *p_m = nullptr;
    ScalarType value, majorityvote, labelsum;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    if (this->m_LabelValues == nullptr) {
        try {this->m_LabelValues = new int[nc];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        // set dummy values
        for (int l = 0; l < nc; ++l){
            this->m_LabelValues[l] = (l+1)*10;
        }
    }
    try {p_labelprobs = new double[nc+1];}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }

    // set dummy values
    ierr = VecGetArrayRead(m, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(labelim, &p_labels); CHKERRQ(ierr);
    for (IntType i = 0; i < nl; ++i) {
        majorityvote  = 0.0;
        majoritylabel = -1;

        // compute value for background
        labelsum = 0.0; p_labelprobs[nc] = 1.0;
        for (int l = 0; l < nc; ++l){
            // get label probability
            p_labelprobs[l] = p_m[l*nl + i];

            // compute background
            p_labelprobs[nc] -= p_labelprobs[l];

            // accumulate label probabilities
            labelsum += p_labelprobs[l];
        }

        // set to zero, if negative
        if (p_labelprobs[nc] < 0.0) p_labelprobs[nc] = 0.0;
        labelsum += p_labelprobs[nc];
//        labelsum = labelsum > 1E-1 ? labelsum : 1.0;

        // normalize (partition of unity)
        for (int l = 0; l < nc+1; ++l) {
            p_labelprobs[l] /= labelsum;
        }

        for (int l = 0; l < nc + 1; ++l) {
            value = p_labelprobs[l];

            // get largest value
            if (value > majorityvote) {
                majoritylabel = l;
                majorityvote  = value;
            }

            // get current label
            if (majoritylabel == nc) {
                p_labels[i] = 0;
            } else {
                p_labels[i] = this->m_LabelValues[majoritylabel];
            }
        }
    }
    ierr = VecRestoreArray(labelim, &p_labels); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(m, &p_m); CHKERRQ(ierr);

    if (p_labelprobs != nullptr) {delete [] p_labelprobs;}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}
*/

/********************************************************************
 * @brief restrict vector field
 * @param v input vector field
 * @param vcoarse output vector field v_c = R[v]
 *******************************************************************/
PetscErrorCode Preprocessing::Restrict(VecField* vcoarse, VecField* vfine, IntType* nx_c, IntType* nx_f) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vfine != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vcoarse != nullptr, "null pointer"); CHKERRQ(ierr);

    ierr = this->Restrict(&vcoarse->m_X1, vfine->m_X1, nx_c, nx_f); CHKERRQ(ierr);
    ierr = this->Restrict(&vcoarse->m_X2, vfine->m_X2, nx_c, nx_f); CHKERRQ(ierr);
    ierr = this->Restrict(&vcoarse->m_X3, vfine->m_X3, nx_c, nx_f); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief restrict data
 * @param x input vector
 * @param xcoarse output vector xcoarse = R[x]
 * @param nx_c number of grid points on coarse grid
 * @param dosetup flag to identify if we have to do the setup step;
 * this is essentially to prevent an additional setup step if we
 * apply this function to each component of a vector field, or a
 * time dependend field; if the parameter is not set, it is true
 *******************************************************************/
PetscErrorCode Preprocessing::Restrict(Vec* x_c, Vec x_f, IntType* nx_c, IntType* nx_f) {
    PetscErrorCode ierr = 0;
    ScalarType scale, coeff[2], value, *p_xc = nullptr;
    const ScalarType *p_xf = nullptr;
    IntType n, l, k_c[3], i_c[3], nr, os_recv, nyqfreqid[3];
    std::stringstream ss;
    int rank, nprocs;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "applying restriction operator ["
           << nx_f[0] << "," << nx_f[1] << "," << nx_f[2] << "]"
           << " -> [" << nx_c[0] << "," << nx_c[1] << "," << nx_c[2] << "]";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    ierr = Assert(x_f != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(x_c != nullptr, "null pointer"); CHKERRQ(ierr);

    if ((nx_c[0] == nx_f[0]) && (nx_c[1] == nx_f[1]) && (nx_c[2] == nx_f[2])) {
        ierr = VecCopy(x_f, *x_c); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    for (int i = 0; i < 3; ++i) {
        this->m_nxC[i] = nx_c[i];
        this->m_nxF[i] = nx_f[i];
        value = static_cast<ScalarType>(nx_c[i])/2.0;
        nyqfreqid[i] = static_cast<IntType>(std::ceil(value));
    }
//    for(int i = 0; i < 3; ++i) {
//        value = static_cast<ScalarType>(nx_c[i])/2.0;
//        nxhalf_c[i] = static_cast<IntType>(std::ceil(value));
//    }

    // set up fft operators
    if (this->m_ResetGridChangeOps) {
        this->m_GridChangeOpsSet = false;
        this->m_GridChangeIndicesComputed = false;
    }
    if (!this->m_GridChangeOpsSet) {
        ierr = this->SetupGridChangeOps(nx_f, nx_c); CHKERRQ(ierr);
    }

    ierr = Assert(this->m_FFTFinePlan != nullptr, "null pointer"); CHKERRQ(ierr);

    n  = this->m_osizeC[0];
    n *= this->m_osizeC[1];
    n *= this->m_osizeC[2];
    
    ierr = this->m_XHatCoarse.AllocateHost(); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    // set freqencies to zero
    for (l = 0; l < n; ++l) {
        this->m_XHatCoarse[l][0] = 0.0;
        this->m_XHatCoarse[l][1] = 0.0;
    }
} // #pragma omp parallel

    // compute fft of data on fine grid
    ierr = GetRawPointerRead(x_f, &p_xf); CHKERRQ(ierr);
    accfft_execute_r2c_t(this->m_FFTFinePlan, const_cast<ScalarType*>(p_xf), this->m_XHatFine.WriteDevice(), timer);
    ierr = RestoreRawPointerRead(x_f, &p_xf); CHKERRQ(ierr);
    ierr = this->m_XHatFine.CopyDeviceToHost(); CHKERRQ(ierr);

    // compute indices
    if (!this->m_GridChangeIndicesComputed) {
        ierr = this->ComputeGridChangeIndices(nx_f, nx_c); CHKERRQ(ierr);
    }
    ierr = this->GridChangeCommDataRestrict(); CHKERRQ(ierr);

    // get grid sizes/fft scales
    scale = this->m_FFTFineScale;

    // get number of entries we are going to assign
    for (int p = 0; p < nprocs; ++p) {
        nr = this->m_NumRecv[p];
        os_recv = this->m_OffsetRecv[p];

        for (IntType j = 0; j < nr; ++j) {
            bool outofbounds = false;

            for (int i = 0; i < 3; ++i) {
//                k_f[i] = this->m_FourierIndicesRecvF[3*j + i + 3*os_recv];
                k_c[i] = this->m_FourierIndicesRecvC[3*j + i + 3*os_recv] ;

                // get wave number index on coarse grid from index on fine grid
//                k_c[i] = k_f[i] <= nxhalf_c[i] ? k_f[i] : nx_c[i] - nx_f[i] + k_f[i];
                i_c[i] = k_c[i] - this->m_ostartC[i];

                if ( (k_c[i] < this->m_ostartC[i]) || (k_c[i] > this->m_ostartC[i] + this->m_osizeC[i]) ) {
                    outofbounds = true;
                }
                if ( (i_c[i] < 0) || (i_c[i] > this->m_osizeC[i]) ) {
                    outofbounds = true;
                }
            }
            if (outofbounds) {
                std::cout << i_c[0] << " " << i_c[1] << " " << i_c[2] << std::endl;
            }
            if (!outofbounds) {
                // compute flat index
                l = GetLinearIndex(i_c[0], i_c[1], i_c[2], this->m_osizeC);
                bool setvalue = true;

                // check for nyquist frequency
                for (int i = 0; i < 3; ++i) {
                    if (i_c[i] == nyqfreqid[i]) {
                        setvalue = false;
                    }
                }

                if (setvalue) {
                    // get fourier coefficients
                    coeff[0] = this->m_FourierCoeffRecvF[2*j + 0 + 2*os_recv];
                    coeff[1] = this->m_FourierCoeffRecvF[2*j + 1 + 2*os_recv];

                    // assign values to coarse grid
                    this->m_XHatCoarse[l][0] = scale*coeff[0];
                    this->m_XHatCoarse[l][1] = scale*coeff[1];
                }
            }
        }
    }

    ierr = GetRawPointerWrite(*x_c, &p_xc); CHKERRQ(ierr);
    accfft_execute_c2r_t(this->m_FFTCoarsePlan, const_cast<ComplexType*>(this->m_XHatCoarse.ReadDevice()), p_xc, timer);
    ierr = RestoreRawPointerWrite(*x_c, &p_xc); CHKERRQ(ierr);

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(timer);

    // increment counter
    this->m_Opt->IncrementCounter(FFT, 2);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief do setup for applying restriction operator
 * @param nx_c grid size on coarse grid
 *******************************************************************/
PetscErrorCode Preprocessing::ComputeGridChangeIndices(IntType* nx_f, IntType* nx_c) {
    PetscErrorCode ierr = 0;
    int rank, nprocs, nowned, ncommunicate, nprocessed, xrank, cart_grid[2], p1, p2;
    IntType oend_c[3], osc_x2, osc_x3, i_f[3], k_f[3], k_c[3], nxhalf_c[3];
    ScalarType nc[2];
    bool locallyowned,oncoarsegrid;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr = Assert(nx_c[0] <= nx_f[0], "grid size in restriction wrong"); CHKERRQ(ierr);
    ierr = Assert(nx_c[1] <= nx_f[1], "grid size in restriction wrong"); CHKERRQ(ierr);
    ierr = Assert(nx_c[2] <= nx_f[2], "grid size in restriction wrong"); CHKERRQ(ierr);

    // allocate if necessary
    if (this->m_IndicesC.empty()) {
        this->m_IndicesC.resize(nprocs);
    }
    if (this->m_IndicesF.empty()) {
        this->m_IndicesF.resize(nprocs);
    }

    for (int i = 0; i < nprocs; ++i) {
        if (!this->m_IndicesC[i].empty()) {
            this->m_IndicesC[i].clear();
        }
        if (!this->m_IndicesF[i].empty()) {
            this->m_IndicesF[i].clear();
        }
    }

    for(int i = 0; i < 3; ++i) {
        nxhalf_c[i] = static_cast<IntType>(std::ceil(static_cast<ScalarType>(nx_c[i])/2.0));
        oend_c[i] = this->m_ostartC[i] + this->m_osizeC[i];
    }

    // get cartesian grid (MPI)
    cart_grid[0] = this->m_Opt->m_CartGridDims[0];
    cart_grid[1] = this->m_Opt->m_CartGridDims[1];

    nc[0] = static_cast<ScalarType>(nx_c[1]);
    nc[1] = static_cast<ScalarType>(nx_c[2])/2.0 + 1.0;
    osc_x2 = static_cast<IntType>(std::ceil(nc[0]/static_cast<ScalarType>(cart_grid[0])));
    osc_x3 = static_cast<IntType>(std::ceil(nc[1]/static_cast<ScalarType>(cart_grid[1])));

    // for all points on fine grid
    nowned = 0; ncommunicate = 0; nprocessed = 0;
    for (i_f[0] = 0; i_f[0] < this->m_osizeF[0]; ++i_f[0]) { // x1
        for (i_f[1] = 0; i_f[1] < this->m_osizeF[1]; ++i_f[1]) { // x2
            for (i_f[2] = 0; i_f[2] < this->m_osizeF[2]; ++i_f[2]) { // x3
                oncoarsegrid = true;
                for (int i = 0; i < 3; ++i) {
                    // compute wave number index on fine grid
                    k_f[i] = i_f[i] + this->m_ostartF[i];

                    // only if current fourier entry is represented in spectral
                    // domain of coarse grid; we ignore the nyquist frequency nx_i/2
                    // because it's not informative
                    if (k_f[i] >= nxhalf_c[i] && k_f[i] <= (nx_f[i]-nxhalf_c[i])) {
                        oncoarsegrid = false;
                    }
                }

                if (oncoarsegrid) {
                    ++nprocessed;
                    locallyowned = true;
                    for (int i = 0; i < 3; ++i) {
                        // get wave number index on coarse grid from index on fine grid
                        k_c[i] = k_f[i] <= nxhalf_c[i] ? k_f[i] : nx_c[i] - nx_f[i] + k_f[i];

                        // sanity checks
                        if ( (k_c[i] < 0.0) || (k_c[i] > nx_c[i]) ) {
                            std::cout << "index out of bounds" << std::endl;
                        }

                        // check if represented on current grid
                        if ( (k_c[i] < this->m_ostartC[i]) || (k_c[i] >= oend_c[i]) ) {
                            locallyowned = false;
                        }
                    }

                    // compute processor id (we do this outside, to check if
                    // we indeed land on the current processor if the points
                    // are owned; sanity check)
                    p1 = static_cast<int>(k_c[1]/osc_x2);
                    p2 = static_cast<int>(k_c[2]/osc_x3);

                    // compute rank
                    xrank = p1*cart_grid[1] + p2;

                    if (locallyowned) { // if owned by local processor
                        // assign computed indices to array (for given rank)
                        for (int i = 0; i < 3; ++i) {
                            this->m_IndicesC[rank].push_back(k_c[i]);
                            this->m_IndicesF[rank].push_back(k_f[i]);
                        }

                        // check if woned is really owned
                        if (rank != xrank) {
                            std::cout << "rank not owned: " << rank << " " << xrank << std::endl;
                        }
                        ++nowned;
                    } else {
                        // assign computed indices to array (for given rank)
                        for (int i = 0; i < 3; ++i) {
                            this->m_IndicesC[xrank].push_back(k_c[i]);
                            this->m_IndicesF[xrank].push_back(k_f[i]);
                        }

                        if (rank == xrank) {
                            std::cout << "rank owned: " << rank << " " << xrank << std::endl;
                        }
                        ++ncommunicate;
                    }
                }
            }  // i1
        }  // i2
    }  // i3

    MPI_Barrier(PETSC_COMM_WORLD);
    // do the communication
    ierr = this->GridChangeCommIndices(); CHKERRQ(ierr);

    this->m_GridChangeIndicesComputed = true;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief communicate indices
 *******************************************************************/
PetscErrorCode Preprocessing::GridChangeCommIndices() {
    PetscErrorCode ierr = 0;
    int merr, nprocs, rank, i_recv, i_send;
    IntType n, k_c[3], k_f[3], os_send, os_recv, nr, ns, n_c, n_f;
    MPI_Status status;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    ierr = AllocateArrayOnce(this->m_OffsetSend, nprocs); CHKERRQ(ierr);
    ierr = AllocateArrayOnce(this->m_OffsetRecv, nprocs); CHKERRQ(ierr);
    ierr = AllocateArrayOnce(this->m_NumSend, nprocs); CHKERRQ(ierr);
    ierr = AllocateArrayOnce(this->m_NumRecv, nprocs); CHKERRQ(ierr);
    ierr = AllocateArrayOnce(this->m_SendRequest, nprocs); CHKERRQ(ierr);
    ierr = AllocateArrayOnce(this->m_RecvRequest, nprocs); CHKERRQ(ierr);

    // compute size to be allocated
    for (int p = 0; p < nprocs; ++p) {
        this->m_NumSend[p] = 0;
        this->m_NumRecv[p] = 0;
        if (!this->m_IndicesF[p].empty()) {
            n_f=this->m_IndicesF[p].size();
            n_c=this->m_IndicesC[p].size();
            ierr = Assert(n_f == n_c,"error in setup"); CHKERRQ(ierr);
            this->m_NumSend[p] = n_f/3;
        }

    }
    // communicate the amount of data we will send from one
    // processor to another (all to all)
    merr = MPI_Alltoall(&this->m_NumSend[0], 1, MPIU_INT, &this->m_NumRecv[0], 1, MPIU_INT, PETSC_COMM_WORLD);
    ierr = MPIERRQ(merr); CHKERRQ(ierr);

    ierr = Assert(this->m_NumSend[rank] == this->m_NumRecv[rank], "alltoall error"); CHKERRQ(ierr);

    // now we compute the size of the arrays, we have to allocate locally to
    // send and recv all the data
    this->m_nAllocSend = this->m_NumSend[0];
    this->m_nAllocRecv = this->m_NumRecv[0];
    this->m_OffsetSend[0] = 0;
    this->m_OffsetRecv[0] = 0;
    for (int p = 1; p < nprocs; ++p) {
        this->m_OffsetSend[p] = this->m_OffsetSend[p-1] + this->m_NumSend[p-1];
        this->m_OffsetRecv[p] = this->m_OffsetRecv[p-1] + this->m_NumRecv[p-1];
        this->m_nAllocSend += this->m_NumSend[p];
        this->m_nAllocRecv += this->m_NumRecv[p];
    }

    // if we actually need to allocate something
    if (this->m_nAllocSend > 0) {
        if (this->m_ResetGridChangeOps) {
            ierr = FreeArray(this->m_FourierIndicesSendC); CHKERRQ(ierr);
            ierr = FreeArray(this->m_FourierIndicesSendF); CHKERRQ(ierr);
        }

        ierr = AllocateArrayOnce(this->m_FourierIndicesSendC, this->m_nAllocSend*3); CHKERRQ(ierr);
        ierr = AllocateArrayOnce(this->m_FourierIndicesSendF, this->m_nAllocSend*3); CHKERRQ(ierr);

        // assign indices accross all procs
        for (int p = 0; p < nprocs; ++p) {
            if (this->m_NumSend[p] != 0) {
                // sanity check
                n = this->m_IndicesF[p].size()/3;

                if (n != this->m_NumSend[p]) {
                    ss << "size mismatch " << n << "!=" << this->m_NumSend[p];
                    ierr = ThrowError(ss.str()); CHKERRQ(ierr);
                }

                // do the assignment
                for (IntType j = 0; j < n; ++j) {
                    // get index
                    for (int i = 0; i < 3; ++i) {
                        k_f[i] = this->m_IndicesF[p][j*3 + i];
                        k_c[i] = this->m_IndicesC[p][j*3 + i];
                    }

                    // get offset
                    os_send = this->m_OffsetSend[p];

                    // assign to flat array
                    this->m_FourierIndicesSendF[3*j+0+3*os_send] = k_f[0];
                    this->m_FourierIndicesSendF[3*j+1+3*os_send] = k_f[1];
                    this->m_FourierIndicesSendF[3*j+2+3*os_send] = k_f[2];

                    // assign to flat array
                    this->m_FourierIndicesSendC[3*j+0+3*os_send] = k_c[0];
                    this->m_FourierIndicesSendC[3*j+1+3*os_send] = k_c[1];
                    this->m_FourierIndicesSendC[3*j+2+3*os_send] = k_c[2];
                }  // for all points
            }  // if indices are not empty
        }  // for all procs
    }  // alloc

    // allocate receiving array
    if (this->m_nAllocRecv > 0) {
        if (this->m_ResetGridChangeOps) {
            ierr = FreeArray(this->m_FourierIndicesRecvF); CHKERRQ(ierr);
            ierr = FreeArray(this->m_FourierIndicesRecvC); CHKERRQ(ierr);
        }

        ierr = AllocateArrayOnce(this->m_FourierIndicesRecvF, this->m_nAllocRecv*3); CHKERRQ(ierr);
        ierr = AllocateArrayOnce(this->m_FourierIndicesRecvC, this->m_nAllocRecv*3); CHKERRQ(ierr);
    } // alloc

    // for all procs, send indices
    for (int i = 0; i < nprocs; ++i) {
        i_send = i; i_recv = i;
        this->m_SendRequest[i_send] = MPI_REQUEST_NULL;
        this->m_RecvRequest[i_recv] = MPI_REQUEST_NULL;

        ns = this->m_NumSend[i];
        os_send = this->m_OffsetSend[i];
        if (ns > 0) {
            ierr = Assert(&this->m_FourierIndicesSendF[3*os_send] != nullptr, "null pointer"); CHKERRQ(ierr);
            merr = MPI_Isend(&this->m_FourierIndicesSendF[3*os_send],
                             3*ns, MPIU_INT, i_send, 0, PETSC_COMM_WORLD,
                             &this->m_SendRequest[i_send]);
            ierr = MPIERRQ(merr); CHKERRQ(ierr);
        }

        nr = this->m_NumRecv[i];
        os_recv = this->m_OffsetRecv[i];
        if (nr > 0) {
            ierr = Assert(&this->m_FourierIndicesRecvF[3*os_recv] != nullptr, "null pointer"); CHKERRQ(ierr);
            merr = MPI_Irecv(&this->m_FourierIndicesRecvF[3*os_recv],
                             3*nr, MPIU_INT, i_recv, 0, PETSC_COMM_WORLD,
                             &this->m_RecvRequest[i_recv]);
            ierr = MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    // we have to wait until all communication is
    // finished before we proceed
    for (int i = 0; i < nprocs; ++i) {
        if (this->m_SendRequest[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&this->m_SendRequest[i], &status);
        }
        if (this->m_RecvRequest[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&this->m_RecvRequest[i], &status);
        }
    }

    // for all procs, send indices
    for (int i = 0; i < nprocs; ++i) {
        i_send = i; i_recv = i;
        this->m_SendRequest[i_send] = MPI_REQUEST_NULL;
        this->m_RecvRequest[i_recv] = MPI_REQUEST_NULL;

        ns = this->m_NumSend[i];
        os_send = this->m_OffsetSend[i];
        if (ns > 0) {
            merr = MPI_Isend(&this->m_FourierIndicesSendC[3*os_send],
                             3*ns, MPIU_INT, i_send, 0, PETSC_COMM_WORLD,
                             &this->m_SendRequest[i_send]);
            ierr = MPIERRQ(merr); CHKERRQ(ierr);
        }

        nr = this->m_NumRecv[i];
        os_recv = this->m_OffsetRecv[i];
        if (nr > 0) {
            merr = MPI_Irecv(&this->m_FourierIndicesRecvC[3*os_recv],
                             3*nr, MPIU_INT, i_recv, 0, PETSC_COMM_WORLD,
                             &this->m_RecvRequest[i_recv]);
            ierr = MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    // we have to wait until all communication is
    // finished before we proceed
    for (int i = 0; i < nprocs; ++i) {
        if (this->m_SendRequest[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&this->m_SendRequest[i], &status);
        }
        if (this->m_RecvRequest[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&this->m_RecvRequest[i], &status);
        }
    }

    // we only have to communicate these indices once
    this->m_IndicesCommunicated = true;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief do setup for applying restriction operator
 * @param nx_c grid size on coarse grid
 *******************************************************************/
PetscErrorCode Preprocessing::GridChangeCommDataRestrict() {
    PetscErrorCode ierr = 0;
    int merr,nprocs,rank,i_recv,i_send;
    IntType n,l,i_f[3],os_send,os_recv,nr,ns;
    MPI_Status status;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr = Assert(this->m_OffsetSend != nullptr, "error in setup"); CHKERRQ(ierr);
    ierr = Assert(this->m_OffsetRecv != nullptr, "error in setup"); CHKERRQ(ierr);
    ierr = Assert(this->m_NumSend != nullptr, "error in setup"); CHKERRQ(ierr);
    ierr = Assert(this->m_NumRecv != nullptr, "error in setup"); CHKERRQ(ierr);

    ierr = AllocateArrayOnce(this->m_SendRequest, nprocs); CHKERRQ(ierr);
    ierr = AllocateArrayOnce(this->m_RecvRequest, nprocs); CHKERRQ(ierr);

    // if we actually need to allocate something
    if (this->m_nAllocSend > 0) {
        if (this->m_ResetGridChangeOps) {
            ierr = FreeArray(this->m_FourierCoeffSendF); CHKERRQ(ierr);
        }

        ierr = AllocateArrayOnce(this->m_FourierCoeffSendF, this->m_nAllocSend*2); CHKERRQ(ierr);

        for (int p = 0; p < nprocs; ++p) {
            n = this->m_NumSend[p];
            if (n != 0) {
                os_send = this->m_OffsetSend[p];

                for (IntType j = 0; j < n; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        // get index (in local space)
                        i_f[i] = this->m_FourierIndicesSendF[3*j + i + 3*os_send] - this->m_ostartF[i];

                        // check if we're inside expected range
                        if ( (i_f[i] >= this->m_osizeF[i]) || (i_f[i] < 0) ) {
                            std::cout << " r "<< rank << " " << i_f[i]
                                      << ">=" << this->m_osizeF[i] << "   " << i_f[i]
                                      << "<0" << std::endl;
                        }
                    }

                    // compute flat index
                    l = GetLinearIndex(i_f[0],i_f[1],i_f[2],this->m_osizeF);

                    // assign values to coarse grid
                    this->m_FourierCoeffSendF[2*j+0+2*os_send] = this->m_XHatFine[l][0];
                    this->m_FourierCoeffSendF[2*j+1+2*os_send] = this->m_XHatFine[l][1];
                }  // for all points
            }  // if indices are not empty
        }  // for all procs
    }

    if (this->m_nAllocRecv > 0) {
        if (this->m_ResetGridChangeOps) {
            ierr = FreeArray(this->m_FourierCoeffRecvF); CHKERRQ(ierr);
        }
        
        ierr = AllocateArrayOnce(this->m_FourierCoeffRecvF, this->m_nAllocRecv*2); CHKERRQ(ierr);
    }

    // send and recv fourier coefficients on fine grid
    for (int i = 0; i < nprocs; ++i) {
        i_send = i; i_recv = i;
        this->m_SendRequest[i_send] = MPI_REQUEST_NULL;
        this->m_RecvRequest[i_recv] = MPI_REQUEST_NULL;

        os_send = this->m_OffsetSend[i];
        ns = this->m_NumSend[i];
        if (ns > 0) {
            merr = MPI_Isend(&this->m_FourierCoeffSendF[2*os_send],
                             2*ns, MPIU_REAL, i_send, 0, PETSC_COMM_WORLD,
                             &this->m_SendRequest[i_send]);
            ierr = MPIERRQ(merr); CHKERRQ(ierr);
        }

        os_recv = this->m_OffsetRecv[i];
        nr      = this->m_NumRecv[i];
        if (nr > 0) {
            merr = MPI_Irecv(&this->m_FourierCoeffRecvF[2*os_recv],
                             2*nr, MPIU_REAL, i_recv, 0, PETSC_COMM_WORLD,
                             &this->m_RecvRequest[i_recv]);
            ierr = MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    for (int i = 0; i < nprocs; ++i) {
        if (this->m_SendRequest[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&this->m_SendRequest[i], &status);
        }
        if (this->m_RecvRequest[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&this->m_RecvRequest[i], &status);
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief communicate data for prolongation (this is the transpose
 * of the restriction operator); we define the restriction as the
 * forward operation and this as the adjoint operation; therefore,
 * we send here, what has been received on the coarse grid
 * @param nx_c grid size on coarse grid
 *******************************************************************/
PetscErrorCode Preprocessing::GridChangeCommDataProlong() {
    PetscErrorCode ierr = 0;
    int merr,nprocs,rank,i_recv,i_send;
    IntType n,l,i_c[3],os_send,os_recv,nr,ns;
    MPI_Status status;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    ierr = Assert(this->m_NumSend != nullptr, "error in setup"); CHKERRQ(ierr);
    ierr = Assert(this->m_NumRecv != nullptr, "error in setup"); CHKERRQ(ierr);
    ierr = Assert(this->m_OffsetSend != nullptr, "error in setup"); CHKERRQ(ierr);
    ierr = Assert(this->m_OffsetRecv != nullptr, "error in setup"); CHKERRQ(ierr);

    ierr = AllocateArrayOnce(this->m_SendRequest, nprocs); CHKERRQ(ierr);
    ierr = AllocateArrayOnce(this->m_RecvRequest, nprocs); CHKERRQ(ierr);

    // if we actually need to allocate something
    if (this->m_nAllocRecv > 0) {
        if (this->m_ResetGridChangeOps) {
            Free(this->m_FourierCoeffSendC);
        }
        ierr = AllocateArrayOnce(this->m_FourierCoeffSendC, this->m_nAllocRecv*2); CHKERRQ(ierr);

        for (int p = 0; p < nprocs; ++p) {
            n = this->m_NumRecv[p];

            if (n != 0) {
                os_recv = this->m_OffsetRecv[p];
                for (IntType j = 0; j < n; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        // get index (in local space)
                        i_c[i] = this->m_FourierIndicesRecvC[3*j + i + 3*os_recv] - this->m_ostartC[i];

                        // check if we're inside expected range
                        if ( (i_c[i] >= this->m_osizeC[i]) || (i_c[i] < 0) ) {
                            std::cout << " r " << rank << " " << i_c[i]
                                      << ">=" << this->m_osizeC[i]
                                      << "   " << i_c[i] << "<0" << std::endl;
                        }
                    }

                    // compute flat index
                    l = GetLinearIndex(i_c[0], i_c[1], i_c[2], this->m_osizeC);

                    // assign values to coarse grid
                    this->m_FourierCoeffSendC[2*j+0+2*os_recv] = this->m_XHatCoarse[l][0];
                    this->m_FourierCoeffSendC[2*j+1+2*os_recv] = this->m_XHatCoarse[l][1];
                }  // for all points
            }  // if indices are not empty
        }  // for all procs
    }

    if (this->m_nAllocSend > 0) {
        if (this->m_ResetGridChangeOps) {
            ierr = FreeArray(this->m_FourierCoeffRecvC); CHKERRQ(ierr);
        }

        ierr = AllocateArrayOnce(this->m_FourierCoeffRecvC, this->m_nAllocSend*2); CHKERRQ(ierr);
    }

    // send and recv fourier coefficients on fine grid
    for (int i = 0; i < nprocs; ++i) {
        i_send = i; i_recv = i;
        this->m_SendRequest[i_send] = MPI_REQUEST_NULL;
        this->m_RecvRequest[i_recv] = MPI_REQUEST_NULL;

        os_send = this->m_OffsetRecv[i];
        ns = this->m_NumRecv[i];
        if (ns > 0) {
            merr = MPI_Isend(&this->m_FourierCoeffSendC[2*os_send],
                             2*ns, MPIU_REAL, i_send, 0, PETSC_COMM_WORLD,
                             &this->m_SendRequest[i_send]);
            ierr = MPIERRQ(merr); CHKERRQ(ierr);
        }

        os_recv = this->m_OffsetSend[i];
        nr      = this->m_NumSend[i];
        if (nr > 0) {
            merr = MPI_Irecv(&this->m_FourierCoeffRecvC[2*os_recv],
                             2*nr, MPIU_REAL, i_recv, 0, PETSC_COMM_WORLD,
                             &this->m_RecvRequest[i_recv]);
            ierr = MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    for (int i = 0; i < nprocs; ++i) {
        if (this->m_SendRequest[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&this->m_SendRequest[i], &status);
        }
        if (this->m_RecvRequest[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&this->m_RecvRequest[i], &status);
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief prolong vector field
 * @param vcoarse input vector field
 * @param vfine output vector field vfine = P[vcoarse]
 *******************************************************************/
PetscErrorCode Preprocessing::Prolong(VecField* v_f, VecField* v_c, IntType* nx_f, IntType* nx_c) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v_f != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(v_c != nullptr, "null pointer"); CHKERRQ(ierr);

    ierr = this->Prolong(&v_f->m_X1, v_c->m_X1, nx_f, nx_c); CHKERRQ(ierr);
    ierr = this->Prolong(&v_f->m_X2, v_c->m_X2, nx_f, nx_c); CHKERRQ(ierr);
    ierr = this->Prolong(&v_f->m_X3, v_c->m_X3, nx_f, nx_c); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief prolong scalar field
 * @param x input vector
 * @param xfine output vector xfine = P[x]
 * @param nx_f number of grid points on fine grid
 * @param dosetup flag to identify if we have to do the setup step;
 * this is essentially to prevent an additional setup step if we
 * apply this function to each component of a vector field, or a
 * time dependend field; if the parameter is not set, it is true
 *******************************************************************/
PetscErrorCode Preprocessing::Prolong(Vec* x_f, Vec x_c, IntType* nx_f, IntType* nx_c) {
    PetscErrorCode ierr = 0;
    int rank, nprocs;
    IntType l, n, ns, os_send, k_f[3], i_f[3], nyqfreqid[3];
    ScalarType *p_xf = nullptr, scale, coeff[2], value;
    const ScalarType *p_xc = nullptr;
    std::stringstream ss;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    ierr = Assert(x_c != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(x_f != nullptr, "null pointer"); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "applying prolongation operator [" << nx_c[0] << "," << nx_c[1] << "," << nx_c[2] << "]"
           << " -> [" << nx_f[0] << "," << nx_f[1] << "," << nx_f[2] << "]";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
    }

    if ( (nx_c[0] == nx_f[0]) && (nx_c[1] == nx_f[1]) && (nx_c[2] == nx_f[2]) ) {
        ierr = VecCopy(x_c, *x_f); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    for (int i = 0; i < 3; ++i) {
        this->m_nxC[i] = nx_c[i];
        this->m_nxF[i] = nx_f[i];
        value = static_cast<ScalarType>(nx_c[i])/2.0;
        nyqfreqid[i] = static_cast<IntType>(std::ceil(value));
    }

    if (this->m_ResetGridChangeOps) {
        this->m_GridChangeOpsSet = false;
        this->m_GridChangeIndicesComputed = false;
    }

    // set up fft operators
    if (!this->m_GridChangeOpsSet) {
        ierr = this->SetupGridChangeOps(nx_f, nx_c); CHKERRQ(ierr);
    }

    n  = this->m_osizeF[0];
    n *= this->m_osizeF[1];
    n *= this->m_osizeF[2];
    
    ierr = this->m_XHatFine.AllocateHost(); CHKERRQ(ierr);

#pragma omp parallel
{
    IntType l;
#pragma omp for
    // set freqencies to zero
    for (l = 0; l < n; ++l) {
        this->m_XHatFine[l][0] = 0.0;
        this->m_XHatFine[l][1] = 0.0;
    }
} // pragma omp parallel

    // compute fft of data on fine grid
    ierr = GetRawPointerRead(x_c, &p_xc); CHKERRQ(ierr);
    accfft_execute_r2c_t(this->m_FFTCoarsePlan, const_cast<ScalarType*>(p_xc), this->m_XHatCoarse.WriteDevice(), timer);
    ierr = RestoreRawPointerRead(x_c, &p_xc); CHKERRQ(ierr);
    ierr = this->m_XHatCoarse.CopyDeviceToHost(); CHKERRQ(ierr);

    // compute indices for mapping from coarse grid to fine grid
    if (!this->m_GridChangeIndicesComputed) {
        ierr = this->ComputeGridChangeIndices(nx_f, nx_c); CHKERRQ(ierr);
    }
    ierr = this->GridChangeCommDataProlong(); CHKERRQ(ierr);

    // get grid sizes/fft scales
    scale = this->m_FFTCoarseScale;

    // get number of entries we are going to assign
    for (int p = 0; p < nprocs; ++p) {
        ns = this->m_NumSend[p];
        os_send = this->m_OffsetSend[p];

        for (IntType j = 0; j < ns; ++j) {
            bool outofbounds = false;
            for (int i = 0; i < 3; ++i) {
                k_f[i] = this->m_FourierIndicesSendF[3*j + i + 3*os_send] ;

                // get wave number index on coarse grid from index on fine grid
//                k_c[i] = k_f[i] <= nxhalf_c[i] ? k_f[i] : nx_c[i] - nx_f[i] + k_f[i];
                i_f[i] = k_f[i] - this->m_ostartF[i];

                if ((k_f[i] < this->m_ostartF[i]) || (k_f[i] > this->m_ostartF[i] + this->m_osizeF[i])) {
                    outofbounds = true;
                }
                if ((i_f[i] < 0) || (i_f[i] > this->m_osizeF[i])) {
                    outofbounds = true;
                }
            }
            if (outofbounds) {
                std::cout << i_f[0] << " " << i_f[1] << " " << i_f[2] << std::endl;
            }
            if (!outofbounds) {

                bool setvalue = true;
                for (int i = 0; i < 3; ++i) {
                    if (i_f[i] == nyqfreqid[i]) {
                        setvalue = false;
                    }
                }

                if (setvalue) {
                    // compute flat index
                    l = GetLinearIndex(i_f[0], i_f[1], i_f[2], this->m_osizeF);

                    // get fourier coefficients
                    coeff[0] = this->m_FourierCoeffRecvC[2*j + 0 + 2*os_send];
                    coeff[1] = this->m_FourierCoeffRecvC[2*j + 1 + 2*os_send];

                    // assign values to coarse grid
                    this->m_XHatFine[l][0] = scale*coeff[0];
                    this->m_XHatFine[l][1] = scale*coeff[1];
                }
            }
        }
    }

    ierr = GetRawPointerWrite(*x_f, &p_xf); CHKERRQ(ierr);
    accfft_execute_c2r_t(this->m_FFTFinePlan, const_cast<ComplexType*>(this->m_XHatFine.ReadDevice()), p_xf, timer);
    ierr = RestoreRawPointerWrite(*x_f, &p_xf); CHKERRQ(ierr);

    // set fft timeri
    this->m_Opt->IncreaseFFTTimers(timer);

    // increment counter
    this->m_Opt->IncrementCounter(FFT, 2);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief apply cutoff frequency filter
 * @param xflt output/filtered x
 * @param x input
 * @param pct cut off precentage (provide 0.5 for 50%)
 * @param lowpass flag to switch on low pass filter; default is true
 *******************************************************************/
PetscErrorCode Preprocessing::ApplyRectFreqFilter(VecField* vflt, VecField* v, ScalarType pct, bool lowpass) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = this->ApplyRectFreqFilter(vflt->m_X1, v->m_X1, pct, lowpass); CHKERRQ(ierr);
    ierr = this->ApplyRectFreqFilter(vflt->m_X2, v->m_X2, pct, lowpass); CHKERRQ(ierr);
    ierr = this->ApplyRectFreqFilter(vflt->m_X3, v->m_X3, pct, lowpass); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief apply cutoff frequency filter
 * @param xflt output/filtered x
 * @param x input
 * @param pct cut off precentage (provide 0.5 for 50%)
 * @param lowpass flag to switch on low pass filter; default is true
 *******************************************************************/
PetscErrorCode Preprocessing::ApplyRectFreqFilter(Vec xflt, Vec x, ScalarType pct, bool lowpass) {
    PetscErrorCode ierr = 0;
    IntType nalloc;
    ScalarType *p_xflt = nullptr;
    const ScalarType *p_x = nullptr;
    ScalarType nxhalf[3], scale, cfreq[3][2], indicator[2], indic;
    int nx[3];
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(x != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xflt != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(pct >= 0.0 && pct <= 1.0, "parameter error"); CHKERRQ(ierr);

    if (pct == 1.0 && lowpass) {
        ierr = VecCopy(x, xflt); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }
    if (pct == 1.0 && !lowpass) {
        ierr = VecSet(xflt, 0.0); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    indicator[0] = 1;
    indicator[1] = 0;

    if (!lowpass) {
        indicator[0] = 0;
        indicator[1] = 1;
    }

    // get local pencil size and allocation size
    nalloc = this->m_Opt->m_FFT.nalloc;

    // allocate
    /*if (this->m_xhat == nullptr) {
        this->m_xhat = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
    }*/
    ierr = this->m_XHat.Resize(nalloc); CHKERRQ(ierr);
    ierr = this->m_XHat.AllocateHost(); CHKERRQ(ierr);

    // get parameters
    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        nxhalf[i] = static_cast<ScalarType>(nx[i]/2.0);
    }
    scale = this->m_Opt->ComputeFFTScale();

    // compute fft
    ierr = GetRawPointerRead(x, &p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(p_x), this->m_XHat.WriteDevice(), timer);
    ierr = RestoreRawPointerRead(x, &p_x); CHKERRQ(ierr);
    ierr = this->m_XHat.CopyDeviceToHost(); CHKERRQ(ierr);

    // compute cutoff frequency
    cfreq[0][0] = pct*(nxhalf[0]-1);
    cfreq[1][0] = pct*(nxhalf[1]-1);
    cfreq[2][0] = pct*(nxhalf[2]-1);

    cfreq[0][1] = static_cast<ScalarType>(nx[0]) - pct*(nxhalf[0]);
    cfreq[1][1] = static_cast<ScalarType>(nx[1]) - pct*(nxhalf[1]);
    cfreq[2][1] = static_cast<ScalarType>(nx[2]) - pct*(nxhalf[2]);

#pragma omp parallel
{
    long int w[3];
    IntType li,i1,i2,i3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->m_FFT.osize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < this->m_Opt->m_FFT.osize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < this->m_Opt->m_FFT.osize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                w[0] = static_cast<ScalarType>(i1 + this->m_Opt->m_FFT.ostart[0]);
                w[1] = static_cast<ScalarType>(i2 + this->m_Opt->m_FFT.ostart[1]);
                w[2] = static_cast<ScalarType>(i3 + this->m_Opt->m_FFT.ostart[2]);

                bool inside = true;
                inside = ( ( (w[0] < cfreq[0][0]) || (w[0] > cfreq[0][1]) ) && inside ) ? true : false;
                inside = ( ( (w[1] < cfreq[1][0]) || (w[1] > cfreq[1][1]) ) && inside ) ? true : false;
                inside = ( ( (w[2] < cfreq[2][0]) || (w[2] > cfreq[2][1]) ) && inside ) ? true : false;

                indic = inside ? indicator[0] : indicator[1];

                // compute linear / flat index
                li = GetLinearIndex(i1, i2, i3, this->m_Opt->m_FFT.osize);

                if (indic == 0) {
                    this->m_XHat[li][0] = 0.0;
                    this->m_XHat[li][1] = 0.0;
                } else {
                    this->m_XHat[li][0] *= scale;
                    this->m_XHat[li][1] *= scale;
                }
            } // i1
        } // i2
    } // i3

} // pragma omp parallel

    // compute inverse fft
    ierr = GetRawPointerWrite(xflt, &p_xflt); CHKERRQ(ierr);
    accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, const_cast<ComplexType*>(this->m_XHat.ReadDevice()), p_xflt, timer);
    ierr = RestoreRawPointerWrite(xflt, &p_xflt); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(timer);

    // increment counter
    this->m_Opt->IncrementCounter(FFT, 2);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief apply gaussian smoothing operator to input data
 *******************************************************************/
PetscErrorCode Preprocessing::Smooth(Vec xs, Vec x, IntType nc) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(x != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xs != nullptr, "null pointer"); CHKERRQ(ierr);

    ierr = this->GaussianSmoothing(xs, x, nc); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief apply gaussian smoothing operator to input data
 *******************************************************************/
PetscErrorCode Preprocessing::GaussianSmoothing(Vec xs, Vec x, IntType nc) {
    PetscErrorCode ierr = 0;
    IntType nalloc, nl;
    std::stringstream ss;
    ScalarType *p_xs = nullptr, c[3], scale; //, nx[3];
    const ScalarType *p_x = nullptr;
    int nx[3];
    double timer[NFFTTIMERS] = {0};
    DifferentiationSM *spectral = nullptr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(x != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xs != nullptr, "null pointer"); CHKERRQ(ierr);

    // get local pencil size and allocation size
    nl     = this->m_Opt->m_Domain.nl;
    nalloc = this->m_Opt->m_FFT.nalloc;
    scale  = this->m_Opt->ComputeFFTScale();

    ierr = AllocateOnce(spectral, this->m_Opt);
    //ierr = spectral->SetupSpectralData(); CHKERRQ(ierr);
    //if (this->m_xhat == nullptr) {
    //    this->m_xhat = reinterpret_cast<ComplexType*>(accfft_alloc(nalloc));
    //}
    //ierr = this->m_XHat.Resize(nalloc/sizeof(ComplexType)); CHKERRQ(ierr);
    //ierr = this->m_XHat.AllocateHost(); CHKERRQ(ierr);
    //ierr = this->m_XHat.AllocateDevice(); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 1) {
        ss << "applying smoothing: ("
           << this->m_Opt->m_Sigma[0]
           << ", " << this->m_Opt->m_Sigma[1]
           << ", " << this->m_Opt->m_Sigma[2] << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }
    
    if (this->m_Opt->m_Verbosity > 3) {
        ScalarType norm;
        ierr = VecNorm(x, NORM_2, &norm);
        ss << "norm: " << norm;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    // get parameters
    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        // sigma is provided by user in # of grid points
        c[i] = this->m_Opt->m_Sigma[i]*this->m_Opt->m_Domain.hx[i];
        c[i] *= c[i];
    }
    
    ierr = GetRawPointerRead(x, &p_x); CHKERRQ(ierr);
    ierr = GetRawPointerWrite(xs, &p_xs); CHKERRQ(ierr);

    for (IntType k = 0; k < nc; ++k) {
        // compute fft
        //ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
       
        ierr = spectral->GaussianFilter(p_xs + k*nl, p_x + k*nl, c); CHKERRQ(ierr);
        
        /*ierr = GetRawPointerRead(x, &p_x); CHKERRQ(ierr);
        accfft_execute_r2c_t(this->m_Opt->m_FFT.plan, const_cast<ScalarType*>(p_x + k*nl), this->m_XHat.WriteDevice(), timer);
        ierr = RestoreRawPointerRead(x, &p_x); CHKERRQ(ierr);
        ierr = this->m_XHat.CopyDeviceToHost(); CHKERRQ(ierr);
        //ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);
#pragma omp parallel
{
        IntType i1, i2, i3, li;
        ScalarType sik;
        long int k1, k2, k3;
#pragma omp for
        for (i1 = 0; i1 < this->m_Opt->m_FFT.osize[0]; ++i1) {  // x1
            for (i2 = 0; i2 < this->m_Opt->m_FFT.osize[1]; ++i2) {  // x2
                for (i3 = 0; i3 < this->m_Opt->m_FFT.osize[2]; ++i3) {  // x3
                    // compute coordinates (nodal grid)
                    k1 = static_cast<long int>(i1 + this->m_Opt->m_FFT.ostart[0]);
                    k2 = static_cast<long int>(i2 + this->m_Opt->m_FFT.ostart[1]);
                    k3 = static_cast<long int>(i3 + this->m_Opt->m_FFT.ostart[2]);

                    // check if grid index is larger or smaller then
                    // half of the total grid size
                    if (k1 > nx[0]/2) k1 -= nx[0];
                    if (k2 > nx[1]/2) k2 -= nx[1];
                    if (k3 > nx[2]/2) k3 -= nx[2];

                    sik = 0.5*( (k1*k1*c[0]) + (k2*k2*c[1]) + (k3*k3*c[2]) );
                    sik = exp(-sik);

                    // compute linear / flat index
                    li = GetLinearIndex(i1, i2, i3, this->m_Opt->m_FFT.osize);

                    this->m_XHat[li][0] *= scale*sik;
                    this->m_XHat[li][1] *= scale*sik;
                }  // i1
            }  // i2
        }  // i3
}  // pragma omp parallel
        ierr = this->m_XHat.CopyHostToDevice(); CHKERRQ(ierr);
        //ierr = VecGetArray(xs, &p_xs); CHKERRQ(ierr);
        ierr = GetRawPointerReadWrite(xs, &p_xs); CHKERRQ(ierr);
        accfft_execute_c2r_t(this->m_Opt->m_FFT.plan, const_cast<ComplexType*>(this->m_XHat.ReadDevice()), p_xs + k*nl, timer);
        //ierr = VecRestoreArray(xs, &p_xs); CHKERRQ(ierr);
        ierr = RestoreRawPointerReadWrite(xs, &p_xs); CHKERRQ(ierr);*/
    }
    
    ierr = RestoreRawPointerRead(x, &p_x); CHKERRQ(ierr);
    ierr = RestoreRawPointerWrite(xs, &p_xs); CHKERRQ(ierr);
    
    if (this->m_Opt->m_Verbosity > 3) {
        ScalarType norm;
        ierr = VecNorm(xs, NORM_2, &norm);
        ss << "norm: " << norm;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }
    
    // increment fft timer
    //this->m_Opt->IncreaseFFTTimers(timer);
    //this->m_Opt->IncrementCounter(FFT, 2);
    
    ierr = Free(spectral); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

}  // namespace reg

#endif   // _PREPROCESSING_CPP_
