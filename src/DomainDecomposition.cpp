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

#ifndef _DOMAINDECOMPOSITION_CPP_
#define _DOMAINDECOMPOSITION_CPP_

#include "DomainDecomposition.hpp"




namespace reg {




/********************************************************************
 * @brief DomainDecomposition
 *******************************************************************/
DomainDecomposition::DomainDecomposition() {
    this->Initialize();
}




/********************************************************************
 * @brief DomainDecomposition
 *******************************************************************/
DomainDecomposition::DomainDecomposition(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief DomainDecomposition
 *******************************************************************/
DomainDecomposition::~DomainDecomposition() {
    this->ClearMemory();
}



/********************************************************************
 * @brief Initialize
 *******************************************************************/
PetscErrorCode DomainDecomposition::Initialize() {
    PetscErrorCode ierr = 0;

    this->m_Opt = NULL;
    this->m_IO = NULL;
    this->m_xhat = NULL;
    this->m_Kxhat = NULL;

    this->m_DDPara.nglobal = 0;
    this->m_DDPara.nsubdom = 0;
    this->m_DDPara.nshared = 0;
    this->m_DDPara.nzeropad = 0;

    this->m_DDPara.nlocal = NULL;
    this->m_DDPara.isize = NULL;
    this->m_DDPara.istart = NULL;
    this->m_DDPara.iend = NULL;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief ClearMemory
 *******************************************************************/
PetscErrorCode DomainDecomposition::ClearMemory() {
    //PetscErrorCode ierr;
    if(this->m_xhat != NULL){
        accfft_free(this->m_xhat);
        this->m_xhat = NULL;
    }
    if(this->m_Kxhat != NULL){
        accfft_free(this->m_Kxhat);
        this->m_Kxhat = NULL;
    }

    if(this->m_DDPara.nlocal != NULL){
        delete this->m_DDPara.nlocal;
        this->m_DDPara.nlocal = NULL;
    }
    if(this->m_DDPara.isize != NULL){
        delete this->m_DDPara.isize;
        this->m_DDPara.isize = NULL;
    }
    if(this->m_DDPara.istart != NULL){
        delete this->m_DDPara.istart;
        this->m_DDPara.istart = NULL;
    }
    if(this->m_DDPara.iend != NULL){
        delete this->m_DDPara.iend;
        this->m_DDPara.iend = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief SetIO
 *******************************************************************/
PetscErrorCode DomainDecomposition::SetIO(DomainDecomposition::ReadWriteType* io) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;
    ierr = Assert(io != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_IO = io;

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief ResetDDData
 *******************************************************************/
PetscErrorCode DomainDecomposition::ResetDDData(int np) {
    PetscErrorCode ierr = 0;

    if(this->m_DDPara.nlocal!=NULL){
        delete this->m_DDPara.nlocal;
        this->m_DDPara.nlocal = NULL;
    }
    if(this->m_DDPara.isize!=NULL){
        delete this->m_DDPara.isize;
        this->m_DDPara.isize = NULL;
    }
    if(this->m_DDPara.istart!=NULL){
        delete this->m_DDPara.istart;
        this->m_DDPara.istart = NULL;
    }
    if(this->m_DDPara.iend!=NULL){
        delete this->m_DDPara.iend;
        this->m_DDPara.iend = NULL;
    }

    try{this->m_DDPara.nlocal = new unsigned long[np];}
    catch (std::bad_alloc&){
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{this->m_DDPara.isize = new unsigned int[3*np];}
    catch (std::bad_alloc&){
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{this->m_DDPara.istart = new int[3*np];}
    catch (std::bad_alloc&){
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{this->m_DDPara.iend = new unsigned int[3*np];}
    catch (std::bad_alloc&){
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief ApplyGaussianSmoothing
 *******************************************************************/
PetscErrorCode DomainDecomposition::ApplyGaussianSmoothing(Vec y, Vec x) {
    PetscErrorCode ierr = 0;
    int isize[3], osize[3], istart[3], ostart[3], n[3];
    IntType iosize[3];
    size_t alloc_max;
    ScalarType hx[3], nx[3], c[3], scale;
    ScalarType *p_x = NULL, *p_y = NULL;
	double timer[7] = {0};
    PetscFunctionBegin;

    ierr = Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    n[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
    n[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
    n[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

    // get local pencil size and allocation size
    alloc_max = accfft_local_size_dft_r2c_t<ScalarType>(n,isize,istart,osize,ostart,
                                                        this->m_Opt->GetComm());
    if(this->m_xhat == NULL){
        this->m_xhat=(FFTScaType*)accfft_alloc(alloc_max);
    }

    for (int i = 0; i < 3; ++i){
        hx[i] = this->m_Opt->GetSpatialStepSize(i);
        nx[i] = static_cast<ScalarType>(n[i]);

        // sigma is provided by user in # of grid points
        c[i] = this->m_Opt->GetSigma()*hx[i];
        c[i] *= c[i];

        iosize[i] = osize[i];
    }

    scale = this->m_Opt->ComputeFFTScale();

    // compute fft
    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFTPlan(), p_x, this->m_xhat, timer);
    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);

#pragma omp parallel {
    IntType i1, i2, i3, li;
    ScalarType k1, k2, k3, sik;
    bool flagx1, flagx2, flagx3;
#pragma omp for
    for (i1 = 0; i1 < iosize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < iosize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < iosize[2]; ++i3) {  // x3

                // compute coordinates (nodal grid)
                k1 = static_cast<ScalarType>(i1 + ostart[0]);
                k2 = static_cast<ScalarType>(i2 + ostart[1]);
                k3 = static_cast<ScalarType>(i3 + ostart[2]);

                // check if grid index is larger or smaller then
                // half of the total grid size
                flagx1 = (k1 <= nx[0]*0.5);
                flagx2 = (k2 <= nx[1]*0.5);
                flagx3 = (k3 <= nx[2]*0.5);

                k1 = flagx1 ? k1 : -nx[0] + k1;
                k2 = flagx2 ? k2 : -nx[1] + k2;
                k3 = flagx3 ? k3 : -nx[2] + k3;

                sik = 0.5*((k1*k1*c[0]) + (k2*k2*c[1]) + (k3*k3*c[2]));
                sik = exp(-sik);

                // compute linear / flat index
                li = GetLinearIndex(i1, i2, i3, iosize);

                this->m_xhat[li][0] *= scale*sik;
                this->m_xhat[li][1] *= scale*sik;

            } // i1
        } // i2
    } // i3

} // pragma omp parallel

    // compute inverse fft
    ierr = VecGetArray(y, &p_y); CHKERRQ(ierr);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFTPlan(), this->m_xhat, p_y, timer);
    ierr = VecRestoreArray(y, &p_y); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief SetupDomainComposition
 *******************************************************************/
PetscErrorCode DomainDecomposition::SetupDomainComposition(unsigned int n) {
    PetscErrorCode ierr = 0;
    ScalarType nx[3],nd;
    IntType isize[3], nsub[3], np, p1, p2, p3, j;
    PetscFunctionBegin;

    this->m_DDPara.nglobal = this->m_Opt->GetNGlobal();
    this->m_DDPara.nsubdom = n;
    this->m_DDPara.nshared = 4;
    this->m_DDPara.nzeropad = 4;

    np = static_cast<IntType>(pow(n,3));
    nd = static_cast<ScalarType>(n);

    ierr = this->ResetDDData(np); CHKERRQ(ierr);

    // compute identifiers for domain decomposition
    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<ScalarType>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<IntType>(ceil(nx[i]/nd));
        nsub[i]  = n;
    }

    for (p1 = 0; p1 < nsub[0]; ++p1) {
        for (p2 = 0; p2 < nsub[1]; ++p2) {
            for (p3 = 0; p3 < nsub[2]; ++p3) {
                j = GetLinearIndex(p1, p2, p3, nsub);

                this->m_DDPara.isize[j*3+0] = isize[0];
                this->m_DDPara.isize[j*3+1] = isize[1];
                this->m_DDPara.isize[j*3+2] = isize[2];

                this->m_DDPara.istart[j*3+0] = p1*isize[0];
                this->m_DDPara.istart[j*3+1] = p2*isize[1];
                this->m_DDPara.istart[j*3+2] = p3*isize[2];

                this->m_DDPara.iend[j*3+0] = (p1+1)*isize[0];
                this->m_DDPara.iend[j*3+1] = (p2+1)*isize[1];
                this->m_DDPara.iend[j*3+2] = (p3+1)*isize[2];

                this->m_DDPara.nlocal[j] = this->m_DDPara.isize[j*3+0]
                                          *this->m_DDPara.isize[j*3+1]
                                          *this->m_DDPara.isize[j*3+2];
            }
        }
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief SetupDomainDecomposition
 *******************************************************************/
PetscErrorCode DomainDecomposition::SetupDomainDecomposition(unsigned int n) {
    PetscErrorCode ierr = 0;
    ScalarType nx[3],nd;
    IntType isize[3],nsub[3],np;
    PetscFunctionBegin;

    this->m_DDPara.nglobal = this->m_Opt->GetNGlobal();
    this->m_DDPara.nsubdom = n;
    this->m_DDPara.nshared = 4;
    this->m_DDPara.nzeropad = 4;

    np = static_cast<IntType>(pow(n,3));
    nd = static_cast<ScalarType>(n);

    ierr = this->ResetDDData(np); CHKERRQ(ierr);

    // compute identifiers for domain decomposition
    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<ScalarType>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<IntType>(ceil(nx[i]/nd));
        nsub[i]  = n;
    }

    unsigned int nshared = this->m_DDPara.nshared;

    for (IntType p1 = 0; p1 < nsub[0]; ++p1){
        for (IntType p2 = 0; p2 < nsub[1]; ++p2){
            for (IntType p3 = 0; p3 < nsub[2]; ++p3){

                IntType j = GetLinearIndex(p1,p2,p3,nsub);

                this->m_DDPara.isize[j*3+0] = isize[0] + 2*nshared;
                this->m_DDPara.isize[j*3+1] = isize[1] + 2*nshared;
                this->m_DDPara.isize[j*3+2] = isize[2] + 2*nshared;

                this->m_DDPara.istart[j*3+0] = p1*isize[0] - nshared;
                this->m_DDPara.istart[j*3+1] = p2*isize[1] - nshared;
                this->m_DDPara.istart[j*3+2] = p3*isize[2] - nshared;

                this->m_DDPara.iend[j*3+0] = (p1+1)*isize[0] + nshared;
                this->m_DDPara.iend[j*3+1] = (p2+1)*isize[1] + nshared;
                this->m_DDPara.iend[j*3+2] = (p3+1)*isize[2] + nshared;

                this->m_DDPara.nlocal[j] = this->m_DDPara.isize[j*3+0]
                                          *this->m_DDPara.isize[j*3+1]
                                          *this->m_DDPara.isize[j*3+2];
            }
        }
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief DecompositionData
 *******************************************************************/
PetscErrorCode DomainDecomposition::DecompositionData(Vec x, unsigned int n, std::string filename) {
    PetscErrorCode ierr = 0;
    Vec yj;
    ScalarType *p_yj=NULL,*p_x=NULL;
    std::string::size_type pos;
    std::ostringstream ss;;
    std::string fn;
    IntType nsub[3],nxblock[3],nzeropad;
    int is[3],ie[3],nxj[3];
    IntType nlj;

    PetscFunctionBegin;

    ierr = Assert(x != NULL,"null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IO != NULL,"null pointer"); CHKERRQ(ierr);

    // remove extension from file
    pos = filename.find_last_of(".");
    filename = (std::string::npos == pos) ? filename : filename.substr(0, pos);

    ierr = this->SetupDomainDecomposition(n); CHKERRQ(ierr);

    for(int i=0; i < 3; ++i){
        this->m_nx[i] = this->m_Opt->GetNumGridPoints(i);
        nsub[i] = n;
    }

    nzeropad = this->m_DDPara.nzeropad;

    ierr = VecGetArray(x,&p_x); CHKERRQ(ierr);

    // for all domains
    for (IntType p1 = 0; p1 < nsub[0]; ++p1){
        for (IntType p2 = 0; p2 < nsub[1]; ++p2){
            for (IntType p3 = 0; p3 < nsub[2]; ++p3){

                IntType j = GetLinearIndex(p1,p2,p3,nsub);

                nlj = 1;

               for(int i=0; i < 3; ++i){

                    ie[i] = this->m_DDPara.iend[j*3+i];
                    is[i] = this->m_DDPara.istart[j*3+i];

                    nxj[i]     = this->m_DDPara.isize[j*3+i] + 2*nzeropad;
                    nxblock[i] = this->m_DDPara.isize[j*3+i] + 2*nzeropad;

                    nlj *= nxblock[i];
                }

                ierr = VecCreate(PETSC_COMM_WORLD,&yj); CHKERRQ(ierr);
                ierr = VecSetSizes(yj,nlj,nlj); CHKERRQ(ierr);
                ierr = VecSetFromOptions(yj); CHKERRQ(ierr);
                ierr = VecSet(yj,0.0); CHKERRQ(ierr);
                ierr = VecGetArray(yj,&p_yj); CHKERRQ(ierr);

                unsigned int k1 = nzeropad;
                for (int i1=is[0]; i1 < ie[0]; ++i1){
                    unsigned int k2 = nzeropad;
                    for (int i2=is[1]; i2 < ie[1]; ++i2){
                        unsigned int k3 = nzeropad;
                        for (int i3=is[2]; i3 < ie[2]; ++i3){

                            IntType lk = GetLinearIndex(k1,k2,k3,nxblock);

                            IntType li = this->GetIndex(i1,i2,i3);

                            p_yj[lk] = p_x[li];
                            ++k3;
                        }
                        ++k2;
                    }
                    ++k1;
                }

                ierr = VecRestoreArray(yj,&p_yj); CHKERRQ(ierr);

                ss << j;
                fn = filename + "-dd-" + ss.str() + ".nc";
                ierr = this->m_IO->WriteBlock(yj,nxj,fn); CHKERRQ(ierr);
                ss.str("");
                ierr = VecDestroy(&yj); CHKERRQ(ierr);
            }
        }
    } // for all domains

    ierr = VecRestoreArray(x,&p_x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief DecompositionData
 *******************************************************************/
PetscErrorCode DomainDecomposition::CompositionData(Vec x, unsigned int n, std::string filename) {
    PetscErrorCode ierr = 0;
    Vec yj;
    ScalarType *p_yj=NULL,*p_x=NULL;
    std::string::size_type pos;
    std::ostringstream ss;;
    std::string fn;
    int is[3],ie[3],nxj[3];
    IntType nlj,nsub[3],nxblock[3];

    PetscFunctionBegin;

    ierr = Assert(x != NULL,"null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IO != NULL,"null pointer"); CHKERRQ(ierr);

    ierr = this->SetupDomainComposition(n); CHKERRQ(ierr);

    for(unsigned int i=0; i < 3; ++i){
        this->m_nx[i] = this->m_Opt->GetNumGridPoints(i);
        nsub[i] = n;
    }

    // remove extension from file
    pos = filename.find_last_of(".");
    filename = (std::string::npos == pos) ? filename : filename.substr(0, pos);

    unsigned int nshared = this->m_DDPara.nshared;

    ierr = VecGetArray(x,&p_x); CHKERRQ(ierr);

    // for all domains
    for (IntType p1 = 0; p1 < nsub[0]; ++p1){
        for (IntType p2 = 0; p2 < nsub[1]; ++p2){
            for (IntType p3 = 0; p3 < nsub[2]; ++p3){

                IntType j = GetLinearIndex(p1,p2,p3,nsub);

                nlj = 1;
                for(int i=0; i < 3; ++i){
                    ie[i] = this->m_DDPara.iend[j*3+i];
                    is[i] = this->m_DDPara.istart[j*3+i];
                    nxj[i] = this->m_DDPara.isize[j*3+i] + 2*nshared;
                    nxblock[i] = this->m_DDPara.isize[j*3+i] + 2*nshared;
                    nlj *= nxj[i];
                }

                ierr = VecCreate(PETSC_COMM_WORLD,&yj); CHKERRQ(ierr);
                ierr = VecSetSizes(yj,nlj,nlj); CHKERRQ(ierr);
                ierr = VecSetFromOptions(yj); CHKERRQ(ierr);
                ierr = VecSet(yj,0.0); CHKERRQ(ierr);

                ss << j;
                fn = this->m_Opt->GetIFolder() + filename + "-dd-" + ss.str() + ".nc";
                ierr = this->m_IO->ReadBlock(yj,nxj,fn); CHKERRQ(ierr);
                ss.str("");

                ierr = VecGetArray(yj,&p_yj); CHKERRQ(ierr);

                IntType j1 = nshared;
                for (IntType i1=is[0]; i1 < ie[0]; ++i1){

                    IntType j2 = nshared;
                    for (IntType i2=is[1]; i2 < ie[1]; ++i2){

                        IntType j3 = nshared;
                        for (IntType i3=is[2]; i3 < ie[2]; ++i3){

                            IntType lk = GetLinearIndex(j1,j2,j3,nxblock);
                            IntType li = this->GetIndex(i1,i2,i3);

                            p_x[li] = p_yj[lk];

                            ++j3;
                        }
                        ++j2;
                    }
                    ++j1;
                }

                ierr = VecRestoreArray(yj,&p_yj); CHKERRQ(ierr);
                ierr = VecDestroy(&yj); CHKERRQ(ierr);

            }
        }
    } // for all domains

    ierr = VecRestoreArray(x,&p_x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief CompositionTimeDependentData
 *******************************************************************/
PetscErrorCode DomainDecomposition::CompositionTimeDependentData(Vec x, unsigned int n, std::string filename) {
    PetscErrorCode ierr;
    ScalarType *p_x=NULL, *p_xj=NULL;
    std::string::size_type pos;
    std::ostringstream ss;
    std::string fn;
    Vec xj;
    IntType nt,nl,ng;

    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ng = this->m_Opt->GetNGlobal();

    ierr = VecCreate(PETSC_COMM_WORLD,&xj); CHKERRQ(ierr);
    ierr = VecSetSizes(xj,nl,ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(xj); CHKERRQ(ierr);

    // remove extension from file
    pos = filename.find_last_of(".");
    filename = (std::string::npos == pos) ? filename : filename.substr(0, pos);

    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);

    for (IntType j = 0; j <= nt; ++j){

        // construct file name
        ss << std::setw(3) << std::setfill('0') << j;
        fn = filename + "-j-" + ss.str() + ".nc"; ss.str("");

        ierr = this->CompositionData(xj,n,fn); CHKERRQ(ierr);

        ierr = VecGetArray(xj,&p_xj); CHKERRQ(ierr);
        try{ std::copy(p_xj,p_xj+nl,p_x+j*nl); }
        catch(std::exception&){
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(xj,&p_xj); CHKERRQ(ierr);

    }

    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief DecompositionTimeDependentData
 *******************************************************************/
PetscErrorCode DomainDecomposition::DecompositionTimeDependentData(Vec x, unsigned int n, std::string filename) {
    PetscErrorCode ierr = 0;
    ScalarType *p_x=NULL, *p_xj=NULL;
    std::string::size_type pos;
    std::ostringstream ss;
    std::string fn;
    Vec xj;
    IntType nl,ng,nt;

    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ng = this->m_Opt->GetNGlobal();

    ierr = VecCreate(PETSC_COMM_WORLD,&xj); CHKERRQ(ierr);
    ierr = VecSetSizes(xj,nl,ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(xj); CHKERRQ(ierr);

    // remove extension from file
    pos = filename.find_last_of(".");
    filename = (std::string::npos == pos) ? filename : filename.substr(0, pos);

    ierr = VecGetArray(x,&p_x); CHKERRQ(ierr);

    for (IntType j = 0; j <= nt; ++j){

        // construct file name
        ss << std::setw(3) << std::setfill('0') << j;
        fn = filename + "-j-" + ss.str() + ".nc"; ss.str("");

        ierr = VecGetArray(xj,&p_xj); CHKERRQ(ierr);
        try{ std::copy(p_x+j*nl,p_x+(j+1)*nl,p_xj); }
        catch(std::exception&){
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(xj,&p_xj); CHKERRQ(ierr);

        ierr = this->DecompositionData(xj,n,fn); CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(x,&p_x); CHKERRQ(ierr);
    ierr = VecDestroy(&xj); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _DOMAINDECOMPOSITION_CPP_
