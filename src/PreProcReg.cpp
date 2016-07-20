#ifndef _PREPROCREG_CPP_
#define _PREPROCREG_CPP_


#include "PreProcReg.hpp"
#include <unistd.h>

namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PreProcReg"
PreProcReg::PreProcReg()
{
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PreProcReg"
PreProcReg::PreProcReg(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief default deconstructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~PreProcReg"
PreProcReg::~PreProcReg()
{
    this->ClearMemory();
}




/********************************************************************
 * @brief initialize
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode PreProcReg::Initialize()
{
    this->m_Opt = NULL;

    this->m_GridChangeOpsSet=false;
    this->m_ResetGridChangeOps=false;
    this->m_GridChangeIndicesComputed=false;
    this->m_IndicesCommunicated=false;

    this->m_ReadWrite = NULL;

    this->m_XHatFine = NULL;
    this->m_XHatCoarse = NULL;
    this->m_FFTFinePlan = NULL;
    this->m_FFTCoarsePlan = NULL;

    this->m_FourierCoeffSendF = NULL;
    this->m_FourierCoeffSendC = NULL;

    this->m_FourierCoeffRecvF = NULL;
    this->m_FourierCoeffRecvC = NULL;

    this->m_FourierIndicesRecvF = NULL;
    this->m_FourierIndicesRecvC = NULL;

    this->m_FourierIndicesSendF = NULL;
    this->m_FourierIndicesSendC = NULL;

    this->m_NumSend=NULL;
    this->m_NumRecv=NULL;

    this->m_OffsetSend=NULL;
    this->m_OffsetRecv=NULL;

    this->m_nAllocSend=0;
    this->m_nAllocRecv=0;

    this->m_xhat = NULL;
    this->m_yhat = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clear memory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode PreProcReg::ClearMemory()
{
    //PetscErrorCode ierr;
    if(this->m_xhat!=NULL){
        accfft_free(this->m_xhat);
        this->m_xhat = NULL;
    }
    if(this->m_yhat!=NULL){
        accfft_free(this->m_yhat);
        this->m_yhat = NULL;
    }

    if(this->m_XHatFine != NULL){
        accfft_free(this->m_XHatFine);
        this->m_XHatFine=NULL;
    }
    if(this->m_XHatCoarse != NULL){
        accfft_free(this->m_XHatCoarse);
        this->m_XHatCoarse=NULL;
    }

    if(this->m_FFTFinePlan != NULL){
        accfft_destroy_plan(this->m_FFTFinePlan);
        this->m_FFTFinePlan=NULL;
    }
    if(this->m_FFTCoarsePlan != NULL){
        accfft_destroy_plan(this->m_FFTCoarsePlan);
        this->m_FFTCoarsePlan=NULL;
    }


    if (this->m_FourierCoeffSendF == NULL){
        delete [] this->m_FourierCoeffSendF;
        this->m_FourierCoeffSendF = NULL;
    }

    if (this->m_FourierCoeffSendC == NULL){
        delete [] this->m_FourierCoeffSendC;
        this->m_FourierCoeffSendC = NULL;
    }

    if (this->m_FourierCoeffRecvF == NULL){
        delete [] this->m_FourierCoeffRecvF;
        this->m_FourierCoeffRecvF = NULL;
    }

    if (this->m_FourierCoeffRecvC == NULL){
        delete [] this->m_FourierCoeffRecvC;
        this->m_FourierCoeffRecvC = NULL;
    }

    if (this->m_FourierIndicesRecvF == NULL){
        delete [] this->m_FourierIndicesRecvF;
        this->m_FourierIndicesRecvF = NULL;
    }

    if (this->m_FourierIndicesRecvC == NULL){
        delete [] this->m_FourierIndicesRecvC;
        this->m_FourierIndicesRecvC = NULL;
    }

    if (this->m_FourierIndicesSendF == NULL){
        delete [] this->m_FourierIndicesSendF;
        this->m_FourierIndicesSendF = NULL;
    }

    if (this->m_FourierIndicesSendC == NULL){
        delete [] this->m_FourierIndicesSendC;
        this->m_FourierIndicesSendC = NULL;
    }

    if (this->m_NumSend==NULL){
        delete [] this->m_NumSend;
        this->m_NumSend=NULL;
    }
    if (this->m_NumRecv==NULL){
        delete [] this->m_NumRecv;
        this->m_NumRecv=NULL;
    }

    if (this->m_OffsetSend==NULL){
        delete [] this->m_OffsetSend;
        this->m_OffsetSend=NULL;
    }
    if (this->m_OffsetRecv==NULL){
        delete [] this->m_OffsetRecv;
        this->m_OffsetRecv=NULL;
    }


    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set read/write object for data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReadWrite"
PetscErrorCode PreProcReg::SetReadWrite(PreProcReg::ReadWriteType* readwrite)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(readwrite!= NULL,"null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite=readwrite;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief do setup for applying prolongation and restriction
 * operators
 * @param nx_f grid size on fine grid
 * @param nx_c grid size on coarse grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupGridChangeOps"
PetscErrorCode PreProcReg::SetupGridChangeOps(IntType* nx_f,IntType* nx_c)
{
    PetscErrorCode ierr;
    IntType nalloc_c,nalloc_f;
    int _nx_f[3],_ostart_f[3],_osize_f[3],_isize_f[3],_istart_f[3],
        _nx_c[3],_ostart_c[3],_osize_c[3],_isize_c[3],_istart_c[3];
    ScalarType *p_xfd=NULL,*p_xcd=NULL;
    Complex *p_xfdhat=NULL,*p_xcdhat=NULL;
    MPI_Comm mpicomm;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_Opt->GetVerbosity() > 1){
        ierr=DbgMsg("setting up grid change operators"); CHKERRQ(ierr);
    }

    if (this->m_XHatCoarse != NULL){
        accfft_free(this->m_XHatCoarse);
        this->m_XHatCoarse=NULL;
    }
    if (this->m_XHatFine != NULL){
        accfft_free(this->m_XHatFine);
        this->m_XHatFine=NULL;
    }
    if(this->m_FFTFinePlan != NULL){
        accfft_destroy_plan(this->m_FFTFinePlan);
        this->m_FFTFinePlan=NULL;
    }
    if(this->m_FFTCoarsePlan != NULL){
        accfft_destroy_plan(this->m_FFTCoarsePlan);
        this->m_FFTCoarsePlan=NULL;
    }

    // parse input sizes
    for (int i=0; i< 3; ++i){
        _nx_f[i] = static_cast<int>(nx_f[i]);
        _nx_c[i] = static_cast<int>(nx_c[i]);
    }

    this->m_FFTFineScale = 1.0;
    this->m_FFTCoarseScale = 1.0;
    for (int i=0; i < 3; ++i){
        this->m_FFTFineScale *= static_cast<ScalarType>(nx_f[i]);
        this->m_FFTCoarseScale *= static_cast<ScalarType>(nx_c[i]);
    }
    this->m_FFTFineScale = 1.0/this->m_FFTFineScale;
    this->m_FFTCoarseScale = 1.0/this->m_FFTCoarseScale;

    // get communicator
    mpicomm = this->m_Opt->GetFFT().mpicomm;


    nalloc_c = accfft_local_size_dft_r2c(_nx_c,_isize_c,_istart_c,_osize_c,_ostart_c,mpicomm);
    nalloc_f = accfft_local_size_dft_r2c(_nx_f,_isize_f,_istart_f,_osize_f,_ostart_f,mpicomm);

    for (int i = 0; i < 3; ++i){
        this->m_osizeC[i] = static_cast<IntType>(_osize_c[i]);
        this->m_osizeF[i] = static_cast<IntType>(_osize_f[i]);

        this->m_ostartC[i] = static_cast<IntType>(_ostart_c[i]);
        this->m_ostartF[i] = static_cast<IntType>(_ostart_f[i]);
    }

    if (this->m_XHatCoarse == NULL){
        this->m_XHatCoarse = (ScalarTypeFD*)accfft_alloc(nalloc_c);
    }

    if (this->m_XHatFine == NULL){
        this->m_XHatFine = (ScalarTypeFD*)accfft_alloc(nalloc_f);
    }

    // allocate plan for fine grid
    if (this->m_FFTFinePlan == NULL){

        p_xfd = (ScalarType*)accfft_alloc(nalloc_f);
        p_xfdhat = (Complex*)accfft_alloc(nalloc_f);

        this->m_FFTFinePlan = accfft_plan_dft_3d_r2c(_nx_f,p_xfd,(double*)p_xfdhat,mpicomm,ACCFFT_MEASURE);

        accfft_free(p_xfd); p_xfd=NULL;
        accfft_free(p_xfdhat); p_xfdhat=NULL;
    }

    // allocate plan for coarse grid
    if (this->m_FFTCoarsePlan == NULL){

        p_xcd = (ScalarType*)accfft_alloc(nalloc_c);
        p_xcdhat = (Complex*)accfft_alloc(nalloc_c);

        this->m_FFTCoarsePlan = accfft_plan_dft_3d_r2c(_nx_c,p_xcd,(double*)p_xcdhat,mpicomm,ACCFFT_MEASURE);

        accfft_free(p_xcd); p_xcd=NULL;
        accfft_free(p_xcdhat); p_xcdhat=NULL;
    }

    // set flag
    this->m_GridChangeOpsSet=true;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief restrict vector field
 * @param v input vector field
 * @param vcoarse output vector field v_c = R[v]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Restrict"
PetscErrorCode PreProcReg::Restrict(VecField* vcoarse, VecField* vfine, IntType* nx_c, IntType* nx_f)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(vfine!=NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(vcoarse!=NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Restrict(&vcoarse->m_X1,vfine->m_X1,nx_c,nx_f); CHKERRQ(ierr);
    ierr=this->Restrict(&vcoarse->m_X2,vfine->m_X2,nx_c,nx_f); CHKERRQ(ierr);
    ierr=this->Restrict(&vcoarse->m_X3,vfine->m_X3,nx_c,nx_f); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

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
#undef __FUNCT__
#define __FUNCT__ "Restrict"
PetscErrorCode PreProcReg::Restrict(Vec* x_c, Vec x_f, IntType* nx_c, IntType* nx_f)
{
    PetscErrorCode ierr=0;
    ScalarType *p_xf=NULL,*p_xc=NULL,scale,coeff[2];//,value
    IntType n,l,k_c[3],i_c[3],nr,os_recv;//,k_f[3],nxhalf_c[3];
    int rank,nprocs;
    double timer[5]={0,0,0,0,0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("applying restriction operator"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr=Assert(x_f!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(x_c!=NULL,"null pointer"); CHKERRQ(ierr);

    if ( (nx_c[0] == nx_f[0]) && (nx_c[1] == nx_f[1]) && (nx_c[2] == nx_f[2]) ){
        ierr=VecCopy(x_f,*x_c); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    for (int i = 0; i < 3; ++i){
        this->m_nxC[i] = nx_c[i];
        this->m_nxF[i] = nx_f[i];
    }
//    for(int i = 0; i < 3; ++i){
//        value = static_cast<ScalarType>(nx_c[i])/2.0;
//        nxhalf_c[i] = static_cast<IntType>(std::ceil(value));
//    }

    // set up fft operators
    if (this->m_ResetGridChangeOps){ this->m_GridChangeOpsSet=false; }
    if (!this->m_GridChangeOpsSet){
        ierr=this->SetupGridChangeOps(nx_f,nx_c); CHKERRQ(ierr);
    }

    n  = this->m_osizeC[0];
    n *= this->m_osizeC[1];
    n *= this->m_osizeC[2];

#pragma omp parallel
{
#pragma omp for
    // set freqencies to zero
    for (l = 0; l < n; ++l){
        this->m_XHatCoarse[l][0] = 0.0;
        this->m_XHatCoarse[l][1] = 0.0;
    }
} // #pragma omp parallel

    // compute fft of data on fine grid
    ierr=VecGetArray(x_f,&p_xf); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_FFTFinePlan,p_xf,this->m_XHatFine,timer);
    ierr=VecRestoreArray(x_f,&p_xf); CHKERRQ(ierr);

    // compute indices
    if (!this->m_GridChangeIndicesComputed){
        ierr=this->ComputeGridChangeIndices(nx_f,nx_c); CHKERRQ(ierr);
    }
    ierr=this->GridChangeCommDataRestrict(); CHKERRQ(ierr);

    // get grid sizes/fft scales
    scale = this->m_FFTFineScale;

    // get number of entries we are going to assign
    for (int p = 0; p < nprocs; ++p){

        nr = this->m_NumRecv[p];
        os_recv = this->m_OffsetRecv[p];

        for (int j = 0; j < nr; ++j){

            bool outofbounds=false;

            for ( int i = 0; i < 3; ++i ){

//                k_f[i] = this->m_FourierIndicesRecvF[3*j + i + 3*os_recv];
                k_c[i] = this->m_FourierIndicesRecvC[3*j + i + 3*os_recv] ;

                // get wave number index on coarse grid from index on fine grid
//                k_c[i] = k_f[i] <= nxhalf_c[i] ? k_f[i] : nx_c[i] - nx_f[i] + k_f[i];
                i_c[i] = k_c[i] - this->m_ostartC[i];

                if ( (k_c[i] < this->m_ostartC[i]) || (k_c[i] > this->m_ostartC[i] + this->m_osizeC[i]) ){
                    outofbounds=true;
                }
                if ( (i_c[i] < 0) || (i_c[i] > this->m_osizeC[i]) ){
                    outofbounds=true;
                }
            }
            if (outofbounds) std::cout << i_c[0] << " " << i_c[1] << " " << i_c[2] << std::endl;
            if (!outofbounds){

                // compute flat index
                l = GetLinearIndex(i_c[0],i_c[1],i_c[2],this->m_osizeC);

                // get fourier coefficients
                coeff[0] = this->m_FourierCoeffRecvF[2*j + 0 + 2*os_recv];
                coeff[1] = this->m_FourierCoeffRecvF[2*j + 1 + 2*os_recv];

                // assign values to coarse grid
                this->m_XHatCoarse[l][0] = scale*coeff[0];
                this->m_XHatCoarse[l][1] = scale*coeff[1];
            }
        }
    }

/*
    n = this->m_IndicesC[rank].size()/3;
#pragma omp parallel
{
    IntType l_f,l_c,k_c,k_f,i_c[3],i_f[3];
#pragma omp for
    for (IntType j = 0; j < n; ++j){

        // compute local index
        for (int i = 0; i < 3; ++i){

            k_c = this->m_IndicesC[rank][j*3 + i];
            k_f = this->m_IndicesF[rank][j*3 + i];

            // convert global index into local index
            i_c[i] = k_c - this->m_ostartC[i];
            i_f[i] = k_f - this->m_ostartF[i];

            // check if we're inside expected range
            if ( (i_c[i] >= this->m_osizeC[i]) || (i_c[i] < 0) ){
                std::cout<<" r "<<rank<<" "<<i_c[i]<<">="<<this->m_osizeC[i]<<"   "<<i_c[i]<<"<0"<<std::endl;
            }
            if ( (i_f[i] >= this->m_osizeF[i]) || (i_f[i] < 0) ){
                std::cout<<" r "<<rank<<" "<<i_f[i]<<">="<<this->m_osizeF[i]<<"   "<<i_f[i]<<"<0"<< std::endl;
            }

        }

        // compute flat index
        l_c = GetLinearIndex(i_c[0],i_c[1],i_c[2],this->m_osizeC);
        l_f = GetLinearIndex(i_f[0],i_f[1],i_f[2],this->m_osizeF);

        // assign values to coarse grid
        this->m_XHatCoarse[l_c][0] = scale*this->m_XHatFine[l_f][0];
        this->m_XHatCoarse[l_c][1] = scale*this->m_XHatFine[l_f][1];

    }

} // pragma omp parallel
*/
    ierr=VecGetArray(*x_c,&p_xc); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(this->m_FFTCoarsePlan,this->m_XHatCoarse,p_xc,timer);
    ierr=VecRestoreArray(*x_c,&p_xc); CHKERRQ(ierr);

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(timer);

    // increment counter
    this->m_Opt->IncrementCounter(FFT,2);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief do setup for applying restriction operator
 * @param nx_c grid size on coarse grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeGridChangeIndices"
PetscErrorCode PreProcReg::ComputeGridChangeIndices(IntType* nx_f,IntType* nx_c)
{
    PetscErrorCode ierr=0;
    int rank,nprocs,nowned,ncommunicate,nprocessed,xrank,cart_grid[2],p1,p2;
    IntType oend_c[3],osc_x2,osc_x3,i_f[3],k_f[3],k_c[3],nxhalf_c[3];
    ScalarType nc[2];
    bool locallyowned,oncoarsegrid;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr=Assert(nx_c[0] <= nx_f[0],"grid size in restriction wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_c[1] <= nx_f[1],"grid size in restriction wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_c[2] <= nx_f[2],"grid size in restriction wrong"); CHKERRQ(ierr);

    // allocate if necessary
    if (this->m_IndicesC.empty()){
        this->m_IndicesC.resize(nprocs);
    }
    if (this->m_IndicesF.empty()){
        this->m_IndicesF.resize(nprocs);
    }

    for (int i = 0; i < nprocs; ++i){
        if (!this->m_IndicesC[i].empty()){
            this->m_IndicesC[i].clear();
        }
        if (!this->m_IndicesF[i].empty()){
            this->m_IndicesF[i].clear();
        }
    }

    for(int i = 0; i < 3; ++i){
        nxhalf_c[i] = static_cast<IntType>(std::ceil(static_cast<ScalarType>(nx_c[i])/2.0));
        oend_c[i] = this->m_ostartC[i] + this->m_osizeC[i];
    }

    // get cartesian grid (MPI)
    cart_grid[0] = this->m_Opt->GetNetworkDims(0);
    cart_grid[1] = this->m_Opt->GetNetworkDims(1);

    nc[0] = static_cast<ScalarType>(nx_c[1]);
    nc[1] = static_cast<ScalarType>(nx_c[2])/2.0 + 1.0;
    osc_x2 = static_cast<IntType>(std::ceil(nc[0]/static_cast<ScalarType>(cart_grid[0])));
    osc_x3 = static_cast<IntType>(std::ceil(nc[1]/static_cast<ScalarType>(cart_grid[1])));

    // for all points on fine grid
    nowned=0; ncommunicate=0; nprocessed=0;
    for (i_f[0] = 0; i_f[0] < this->m_osizeF[0]; ++i_f[0]){ // x1
        for (i_f[1] = 0; i_f[1] < this->m_osizeF[1]; ++i_f[1]){ // x2
            for (i_f[2] = 0; i_f[2] < this->m_osizeF[2]; ++i_f[2]){ // x3

                oncoarsegrid=true;

                for (int i = 0; i < 3; ++i){
                    // compute wave number index on fine grid
                    k_f[i] = i_f[i] + this->m_ostartF[i];

                    // only if current fourier entry is represented in spectral
                    // domain of coarse grid; we ignore the nyquist frequency nx_i/2
                    // because it's not informative
                    if ( k_f[i] >= nxhalf_c[i] && k_f[i] <= (nx_f[i]-nxhalf_c[i]) ) {
                        oncoarsegrid=false;
                    }
                }

                if (oncoarsegrid){

                    ++nprocessed;

                    locallyowned=true;

                    for (int i=0; i < 3; ++i){

                        // get wave number index on coarse grid from index on fine grid
                        k_c[i] = k_f[i] <= nxhalf_c[i] ? k_f[i] : nx_c[i] - nx_f[i] + k_f[i];

                        // sanity checks
                        if ( (k_c[i] < 0.0) || (k_c[i] > nx_c[i]) ){
                            std::cout<<"index out of bounds"<<std::endl;
                        }

                        // check if represented on current grid
                        if ( (k_c[i] < this->m_ostartC[i]) || (k_c[i] >= oend_c[i]) ){
                            locallyowned = false;
                        }

                    }

                    // compute processor id (we do this outside, to check if
                    // we indeed land on the current processor if the points
                    // are owned; sanity check)
                    p1=static_cast<int>(k_c[1]/osc_x2);
                    p2=static_cast<int>(k_c[2]/osc_x3);

                    // compute rank
                    xrank = p1*cart_grid[1] + p2;

                    if ( locallyowned ){ // if owned by local processor

                        // assign computed indices to array (for given rank)
                        for (int i=0; i < 3; ++i){
                            this->m_IndicesC[rank].push_back(k_c[i]);
                            this->m_IndicesF[rank].push_back(k_f[i]);
                        }

                        // check if woned is really owned
                        if (rank!=xrank){std::cout<<"rank not owned: "<<rank<<" "<<xrank<<std::endl;}
                        ++nowned;
                    }
                    else{

                        // assign computed indices to array (for given rank)
                        for (int i=0; i < 3; ++i){
                            this->m_IndicesC[xrank].push_back(k_c[i]);
                            this->m_IndicesF[xrank].push_back(k_f[i]);
                        }

                        if (rank==xrank){std::cout<<"rank owned: "<<rank<<" "<<xrank<<std::endl;}
                        ++ncommunicate;

                    }
                }
            } // i1
        } // i2
    } // i3

    ierr=this->GridChangeCommIndices(); CHKERRQ(ierr);

    this->m_GridChangeIndicesComputed = true;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief communicate indices
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GridChangeCommIndices"
PetscErrorCode PreProcReg::GridChangeCommIndices()
{
    PetscErrorCode ierr=0;
    int merr,nprocs,rank,i_recv,i_send;
    IntType n,k_c[3],k_f[3],os_send,os_recv,nr,ns,n_c,n_f;
    MPI_Request *send_request=NULL,*recv_request=NULL;
    MPI_Status status;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    if (this->m_OffsetSend==NULL){
        try{this->m_OffsetSend = new IntType[nprocs];}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_OffsetRecv==NULL){
        try{this->m_OffsetRecv = new IntType[nprocs];}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_NumSend==NULL){
        try{this->m_NumSend = new IntType[nprocs];}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if(this->m_NumRecv==NULL){
        try{this->m_NumRecv = new IntType[nprocs];}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    try{send_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{recv_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }


    // compute size to be allocated
    for (int p = 0; p < nprocs; ++p){

        this->m_NumSend[p] = 0;
        this->m_NumRecv[p] = 0;

        if (!this->m_IndicesF[p].empty()){
            n_f = this->m_IndicesF[p].size();
            n_c = this->m_IndicesC[p].size();
            ierr=Assert(n_f==n_c, "error in setup"); CHKERRQ(ierr);
            this->m_NumSend[p] = n_f/3;
        }

    }

    // communicate the amount of data we will send from one
    // processor to another
    MPI_Alltoall(this->m_NumSend,1,MPI_INT,this->m_NumRecv,1,MPI_INT,PETSC_COMM_WORLD);

    this->m_nAllocSend=this->m_NumSend[0];
    this->m_nAllocRecv=this->m_NumRecv[0];
    this->m_OffsetSend[0] = 0;
    this->m_OffsetRecv[0] = 0;
    for (int p = 1; p < nprocs; ++p){
        this->m_OffsetSend[p] = this->m_OffsetSend[p-1] + this->m_NumSend[p-1];
        this->m_OffsetRecv[p] = this->m_OffsetRecv[p-1] + this->m_NumRecv[p-1];
        this->m_nAllocSend += this->m_NumSend[p];
        this->m_nAllocRecv += this->m_NumRecv[p];
    }

    // if we actually need to allocate something
    if (this->m_nAllocSend > 0){

        if (this->m_FourierIndicesSendF==NULL){
            try{this->m_FourierIndicesSendF = new IntType[this->m_nAllocSend*3];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        if (this->m_FourierIndicesSendC==NULL){
            try{this->m_FourierIndicesSendC = new IntType[this->m_nAllocSend*3];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        // assign indices accross all procs
        for (int p = 0; p < nprocs; ++p){

            if (this->m_NumSend[p] != 0){

                // sanity check
                n = this->m_IndicesF[p].size()/3;

                if (n != this->m_NumSend[p]){
                    ss << "size mismatch " << n << "!=" << this->m_NumSend[p];
                    ierr=ThrowError(ss.str()); CHKERRQ(ierr);
                }

                // do the assignment
                for (int j = 0; j < n; ++j){
                    // get index
                    for (int i = 0; i < 3; ++i){
                        k_f[i] = this->m_IndicesF[p][j*3 + i];
                        k_c[i] = this->m_IndicesC[p][j*3 + i];
                    }

                    // get offset
                    os_send = this->m_OffsetSend[p];

                    // assign to flat array
                    this->m_FourierIndicesSendF[3*j + 0 + 3*os_send] = k_f[0];
                    this->m_FourierIndicesSendF[3*j + 1 + 3*os_send] = k_f[1];
                    this->m_FourierIndicesSendF[3*j + 2 + 3*os_send] = k_f[2];

                    // assign to flat array
                    this->m_FourierIndicesSendC[3*j + 0 + 3*os_send] = k_c[0];
                    this->m_FourierIndicesSendC[3*j + 1 + 3*os_send] = k_c[1];
                    this->m_FourierIndicesSendC[3*j + 2 + 3*os_send] = k_c[2];

                } // for all points
            } // if indices are not empty
        } // for all procs

    } // alloc

    // allocate receiving array
    if (this->m_nAllocRecv > 0){

        if (this->m_FourierIndicesRecvF==NULL){
            try{this->m_FourierIndicesRecvF = new IntType[this->m_nAllocRecv*3];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        if (this->m_FourierIndicesRecvC==NULL){
            try{this->m_FourierIndicesRecvC = new IntType[this->m_nAllocRecv*3];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

    } // alloc

    // for all procs, send indices
    for (int i = 0; i < nprocs; ++i){

        i_send=i; i_recv=i;
        send_request[i_send] = MPI_REQUEST_NULL;
        recv_request[i_recv] = MPI_REQUEST_NULL;

        ns = this->m_NumSend[i];
        os_send = this->m_OffsetSend[i];
        if (ns > 0){
            merr=MPI_Isend(&this->m_FourierIndicesSendF[3*os_send],3*ns,MPI_INT,i_send,0,PETSC_COMM_WORLD,&send_request[i_send]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }

        nr = this->m_NumRecv[i];
        os_recv = this->m_OffsetRecv[i];
        if (nr > 0){
            merr=MPI_Irecv(&this->m_FourierIndicesRecvF[3*os_recv],3*nr,MPI_INT,i_recv,0,PETSC_COMM_WORLD,&recv_request[i_recv]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    // we have to wait until all communication is
    // finished before we proceed
    for (int i=0; i < nprocs; ++i){
        if(send_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&send_request[i], &status); }
        if(recv_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&recv_request[i], &status); }
    }


    // for all procs, send indices
    for (int i = 0; i < nprocs; ++i){

        i_send=i; i_recv=i;
        send_request[i_send] = MPI_REQUEST_NULL;
        recv_request[i_recv] = MPI_REQUEST_NULL;

        ns = this->m_NumSend[i];
        os_send = this->m_OffsetSend[i];
        if (ns > 0){
            merr=MPI_Isend(&this->m_FourierIndicesSendC[3*os_send],3*ns,MPI_INT,i_send,0,PETSC_COMM_WORLD,&send_request[i_send]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }

        nr = this->m_NumRecv[i];
        os_recv = this->m_OffsetRecv[i];
        if (nr > 0){
            merr=MPI_Irecv(&this->m_FourierIndicesRecvC[3*os_recv],3*nr,MPI_INT,i_recv,0,PETSC_COMM_WORLD,&recv_request[i_recv]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    // we have to wait until all communication is
    // finished before we proceed
    for (int i=0; i < nprocs; ++i){
        if(send_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&send_request[i], &status); }
        if(recv_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&recv_request[i], &status); }
    }


    // we only have to communicate these indices once
    this->m_IndicesCommunicated=true;

    if (send_request!=NULL) { delete [] send_request; send_request=NULL; }
    if (recv_request!=NULL) { delete [] recv_request; recv_request=NULL; }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup for applying restriction operator
 * @param nx_c grid size on coarse grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GridChangeCommDataRestrict"
PetscErrorCode PreProcReg::GridChangeCommDataRestrict()
{
    PetscErrorCode ierr=0;
    int merr,nprocs,rank,i_recv,i_send;
    IntType n,l,i_f[3],os_send,os_recv,nr,ns;
    MPI_Request *send_request=NULL,*recv_request=NULL;
    MPI_Status status;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr=Assert(this->m_OffsetSend!=NULL,"error in setup"); CHKERRQ(ierr);
    ierr=Assert(this->m_OffsetRecv!=NULL,"error in setup"); CHKERRQ(ierr);
    ierr=Assert(this->m_NumSend!=NULL,"error in setup"); CHKERRQ(ierr);
    ierr=Assert(this->m_NumRecv!=NULL,"error in setup"); CHKERRQ(ierr);

    try{send_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{recv_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }


    // if we actually need to allocate something
    if (this->m_nAllocSend > 0){

        if (this->m_FourierCoeffSendF==NULL){
            try{this->m_FourierCoeffSendF = new ScalarType[this->m_nAllocSend*2];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        for (int p = 0; p < nprocs; ++p){
            n = this->m_NumSend[p];

            if (n != 0){

                os_send = this->m_OffsetSend[p];

                for (int j = 0; j < n; ++j){
                    for (int i = 0; i < 3; ++i){

                        // get index (in local space)
                        i_f[i] = this->m_FourierIndicesSendF[3*j + i + 3*os_send] - this->m_ostartF[i];

                        // check if we're inside expected range
                        if ( (i_f[i] >= this->m_osizeF[i]) || (i_f[i] < 0) ){
                            std::cout<<" r "<<rank<<" "<<i_f[i]<<">="<<this->m_osizeF[i]<<"   "<<i_f[i]<<"<0"<< std::endl;
                        }
                    }

                    // compute flat index
                    l = GetLinearIndex(i_f[0],i_f[1],i_f[2],this->m_osizeF);

                    // assign values to coarse grid
                    this->m_FourierCoeffSendF[2*j + 0 + 2*os_send] = this->m_XHatFine[l][0];
                    this->m_FourierCoeffSendF[2*j + 1 + 2*os_send] = this->m_XHatFine[l][1];

                } // for all points
            } // if indices are not empty
        } // for all procs
    }

    if (this->m_nAllocRecv > 0){

        if (this->m_FourierCoeffRecvF==NULL){
            try{this->m_FourierCoeffRecvF = new ScalarType[this->m_nAllocRecv*2];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

    }

    // send and recv fourier coefficients on fine grid
    for (int i = 0; i < nprocs; ++i){

        i_send=i; i_recv=i;
        send_request[i_send] = MPI_REQUEST_NULL;
        recv_request[i_recv] = MPI_REQUEST_NULL;

        os_send = this->m_OffsetSend[i];
        ns = this->m_NumSend[i];
        if (ns > 0){
            merr=MPI_Isend(&this->m_FourierCoeffSendF[2*os_send],2*ns,MPI_DOUBLE,i_send,0,PETSC_COMM_WORLD,&send_request[i_send]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }

        os_recv = this->m_OffsetRecv[i];
        nr = this->m_NumRecv[i];
        if (nr > 0){
            merr=MPI_Irecv(&this->m_FourierCoeffRecvF[2*os_recv],2*nr,MPI_DOUBLE,i_recv,0,PETSC_COMM_WORLD,&recv_request[i_recv]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    for (int i=0; i < nprocs; ++i){
        if(send_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&send_request[i], &status); }
        if(recv_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&recv_request[i], &status); }
    }

    if (send_request!=NULL) { delete [] send_request; send_request=NULL; }
    if (recv_request!=NULL) { delete [] recv_request; recv_request=NULL; }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief communicate data for prolongation (this is the transpose
 * of the restriction operator); we define the restriction as the
 * forward operation and this as the adjoint operation; therefore,
 * we send here, what has been received on the coarse grid
 * @param nx_c grid size on coarse grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GridChangeCommDataProlong"
PetscErrorCode PreProcReg::GridChangeCommDataProlong()
{
    PetscErrorCode ierr=0;
    int merr,nprocs,rank,i_recv,i_send;
    IntType n,l,i_c[3],os_send,os_recv,nr,ns;
    MPI_Request *send_request=NULL,*recv_request=NULL;
    MPI_Status status;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr=Assert(this->m_OffsetSend!=NULL,"error in setup"); CHKERRQ(ierr);
    ierr=Assert(this->m_OffsetRecv!=NULL,"error in setup"); CHKERRQ(ierr);
    ierr=Assert(this->m_NumSend!=NULL,"error in setup"); CHKERRQ(ierr);
    ierr=Assert(this->m_NumRecv!=NULL,"error in setup"); CHKERRQ(ierr);

    try{send_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{recv_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }


    // if we actually need to allocate something
    if (this->m_nAllocRecv > 0){

        if (this->m_FourierCoeffSendC==NULL){
            try{this->m_FourierCoeffSendC = new ScalarType[this->m_nAllocRecv*2];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        for (int p = 0; p < nprocs; ++p){
            n = this->m_NumRecv[p];

            if (n != 0){

                os_recv = this->m_OffsetRecv[p];

                for (int j = 0; j < n; ++j){
                    for (int i = 0; i < 3; ++i){

                        // get index (in local space)
                        i_c[i] = this->m_FourierIndicesRecvC[3*j + i + 3*os_recv] - this->m_ostartC[i];

                        // check if we're inside expected range
                        if ( (i_c[i] >= this->m_osizeC[i]) || (i_c[i] < 0) ){
                            std::cout<<" r "<<rank<<" "<<i_c[i]<<">="<<this->m_osizeC[i]<<"   "<<i_c[i]<<"<0"<< std::endl;
                        }
                    }

                    // compute flat index
                    l = GetLinearIndex(i_c[0],i_c[1],i_c[2],this->m_osizeC);

                    // assign values to coarse grid
                    this->m_FourierCoeffSendC[2*j + 0 + 2*os_recv] = this->m_XHatCoarse[l][0];
                    this->m_FourierCoeffSendC[2*j + 1 + 2*os_recv] = this->m_XHatCoarse[l][1];

                } // for all points
            } // if indices are not empty
        } // for all procs
    }

    if (this->m_nAllocSend > 0){

        if (this->m_FourierCoeffRecvC==NULL){
            try{this->m_FourierCoeffRecvC = new ScalarType[this->m_nAllocSend*2];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

    }

    // send and recv fourier coefficients on fine grid
    for (int i = 0; i < nprocs; ++i){

        i_send=i; i_recv=i;
        send_request[i_send] = MPI_REQUEST_NULL;
        recv_request[i_recv] = MPI_REQUEST_NULL;

        os_send = this->m_OffsetRecv[i];
        ns = this->m_NumRecv[i];
        if (ns > 0){
            merr=MPI_Isend(&this->m_FourierCoeffSendC[2*os_send],2*ns,MPI_DOUBLE,i_send,0,PETSC_COMM_WORLD,&send_request[i_send]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }

        os_recv = this->m_OffsetSend[i];
        nr = this->m_NumSend[i];
        if (nr > 0){
            merr=MPI_Irecv(&this->m_FourierCoeffRecvC[2*os_recv],2*nr,MPI_DOUBLE,i_recv,0,PETSC_COMM_WORLD,&recv_request[i_recv]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    for (int i=0; i < nprocs; ++i){
        if(send_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&send_request[i], &status); }
        if(recv_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&recv_request[i], &status); }
    }

    if (send_request!=NULL) { delete [] send_request; send_request=NULL; }
    if (recv_request!=NULL) { delete [] recv_request; recv_request=NULL; }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief prolong vector field
 * @param vcoarse input vector field
 * @param vfine output vector field vfine = P[vcoarse]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Prolong"
PetscErrorCode PreProcReg::Prolong(VecField* v_f, VecField* v_c, IntType* nx_f, IntType* nx_c)
{
    PetscErrorCode ierr;

    ierr=Assert(v_f!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(v_c!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->Prolong(&v_f->m_X1,v_c->m_X1,nx_f,nx_c); CHKERRQ(ierr);
    ierr=this->Prolong(&v_f->m_X2,v_c->m_X2,nx_f,nx_c); CHKERRQ(ierr);
    ierr=this->Prolong(&v_f->m_X3,v_c->m_X3,nx_f,nx_c); CHKERRQ(ierr);

    PetscFunctionReturn(0);

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
#undef __FUNCT__
#define __FUNCT__ "Prolong"
PetscErrorCode PreProcReg::Prolong(Vec* x_f, Vec x_c, IntType* nx_f, IntType* nx_c)
{
    PetscErrorCode ierr=0;
    int rank,nprocs;
    IntType l,n,ns,os_send,k_f[3],i_f[3];
    ScalarType *p_xf=NULL,*p_xc=NULL,scale,coeff[2];
    double timer[5]={0,0,0,0,0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("applying prolongation operator"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr=Assert(x_c!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(x_f!=NULL,"null pointer"); CHKERRQ(ierr);

    if ( (nx_c[0] == nx_f[0]) && (nx_c[1] == nx_f[1]) && (nx_c[2] == nx_f[2]) ){
        ierr=VecCopy(x_c,*x_f); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    for (int i = 0; i < 3; ++i){
        this->m_nxC[i] = nx_c[i];
        this->m_nxF[i] = nx_f[i];
    }

    if (this->m_ResetGridChangeOps){
        this->m_GridChangeOpsSet=false;
    }

    // set up fft operators
    if (!this->m_GridChangeOpsSet){
        ierr=this->SetupGridChangeOps(nx_f,nx_c); CHKERRQ(ierr);
    }

    n  = this->m_osizeF[0];
    n *= this->m_osizeF[1];
    n *= this->m_osizeF[2];

#pragma omp parallel
{
    IntType l;
#pragma omp for
    // set freqencies to zero
    for (l = 0; l < n; ++l){
        this->m_XHatFine[l][0] = 0.0;
        this->m_XHatFine[l][1] = 0.0;
    }
} // pragma omp parallel

    // compute fft of data on fine grid
    ierr=VecGetArray(x_c,&p_xc); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_FFTCoarsePlan,p_xc,this->m_XHatCoarse,timer);
    ierr=VecRestoreArray(x_c,&p_xc); CHKERRQ(ierr);

    // compute indices for mapping from coarse grid to fine grid
    if (!this->m_GridChangeIndicesComputed){
        ierr=this->ComputeGridChangeIndices(nx_f,nx_c); CHKERRQ(ierr);
    }
    ierr=this->GridChangeCommDataProlong(); CHKERRQ(ierr);

    // get grid sizes/fft scales
    scale = this->m_FFTCoarseScale;

    // get number of entries we are going to assign
    for (int p = 0; p < nprocs; ++p){

        ns = this->m_NumSend[p];
        os_send = this->m_OffsetSend[p];

        for (int j = 0; j < ns; ++j){

            bool outofbounds=false;

            for ( int i = 0; i < 3; ++i ){

                k_f[i] = this->m_FourierIndicesSendF[3*j + i + 3*os_send] ;

                // get wave number index on coarse grid from index on fine grid
//                k_c[i] = k_f[i] <= nxhalf_c[i] ? k_f[i] : nx_c[i] - nx_f[i] + k_f[i];
                i_f[i] = k_f[i] - this->m_ostartF[i];

                if ( (k_f[i] < this->m_ostartF[i]) || (k_f[i] > this->m_ostartF[i] + this->m_osizeF[i]) ){
                    outofbounds=true;
                }
                if ( (i_f[i] < 0) || (i_f[i] > this->m_osizeF[i]) ){
                    outofbounds=true;
                }
            }
            if (outofbounds) std::cout << i_f[0] << " " << i_f[1] << " " << i_f[2] << std::endl;
            if (!outofbounds){

                // compute flat index
                l = GetLinearIndex(i_f[0],i_f[1],i_f[2],this->m_osizeF);

                // get fourier coefficients
                coeff[0] = this->m_FourierCoeffRecvC[2*j + 0 + 2*os_send];
                coeff[1] = this->m_FourierCoeffRecvC[2*j + 1 + 2*os_send];

                // assign values to coarse grid
                this->m_XHatFine[l][0] = scale*coeff[0];
                this->m_XHatFine[l][1] = scale*coeff[1];
            }
        }
    }





/*
    // get number of entries we are going to assign
    n = this->m_IndicesC[rank].size()/3;

#pragma omp parallel
{
    IntType l_f,l_c,k_c,k_f,i_c[3],i_f[3];
#pragma omp for
    for (IntType j = 0; j < n; ++j){

        // compute local index
        for (int i = 0; i < 3; ++i){

            k_f = this->m_IndicesF[rank][j*3 + i];
            k_c = this->m_IndicesC[rank][j*3 + i];

            i_c[i] = k_c - this->m_ostartC[i];
            i_f[i] = k_f - this->m_ostartF[i];

            if ( (i_c[i] >= this->m_osizeC[i]) || (i_c[i] < 0.0) ){
                std::cout<<" r "<<rank<<" "<<i_c[i]<<">="<<this->m_osizeC[i]<<"   "<<i_c[i]<<"<0"<<std::endl;
            }
            if ( (i_f[i] >= this->m_osizeF[i]) || (i_f[i] < 0.0) ){
                std::cout<<" r "<<rank<<" "<<i_f[i]<<">="<<this->m_osizeF[i]<<"   "<<i_f[i]<<"<0"<< std::endl;
            }
        }
        // compute flat index
        l_c = GetLinearIndex(i_c[0],i_c[1],i_c[2],this->m_osizeC);
        l_f = GetLinearIndex(i_f[0],i_f[1],i_f[2],this->m_osizeF);

        // assign values to coarse grid
        this->m_XHatFine[l_f][0] = scale*this->m_XHatCoarse[l_c][0];
        this->m_XHatFine[l_f][1] = scale*this->m_XHatCoarse[l_c][1];

    }

} // pragma omp parallel
*/


    ierr=VecGetArray(*x_f,&p_xf); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(this->m_FFTFinePlan,this->m_XHatFine,p_xf,timer);
    ierr=VecRestoreArray(*x_f,&p_xf); CHKERRQ(ierr);

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(timer);

    // increment counter
    this->m_Opt->IncrementCounter(FFT,2);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief apply cutoff frequency filter
 * @param xflt output/filtered x
 * @param x input
 * @param pct cut off precentage (provide 0.5 for 50%)
 * @param lowpass flag to switch on low pass filter; default is true
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyRectFreqFilter"
PetscErrorCode PreProcReg::ApplyRectFreqFilter(VecField* vflt,VecField* v,ScalarType pct,bool lowpass)
{
    PetscErrorCode ierr=0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=this->ApplyRectFreqFilter(vflt->m_X1,v->m_X1,pct,lowpass); CHKERRQ(ierr);
    ierr=this->ApplyRectFreqFilter(vflt->m_X2,v->m_X2,pct,lowpass); CHKERRQ(ierr);
    ierr=this->ApplyRectFreqFilter(vflt->m_X3,v->m_X3,pct,lowpass); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply cutoff frequency filter
 * @param xflt output/filtered x
 * @param x input
 * @param pct cut off precentage (provide 0.5 for 50%)
 * @param lowpass flag to switch on low pass filter; default is true
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyRectFreqFilter"
PetscErrorCode PreProcReg::ApplyRectFreqFilter(Vec xflt,Vec x,ScalarType pct,bool lowpass)
{
    PetscErrorCode ierr=0;
    IntType nalloc;
    ScalarType *p_x=NULL,*p_xflt=NULL;
    ScalarType nxhalf[3],scale,cfreq[3][2],indicator[2],indic;
    int nx[3];
    double timer[5]={0,0,0,0,0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(xflt!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(pct >= 0.0 && pct <= 1.0,"parameter error"); CHKERRQ(ierr);

    if (pct == 1.0 && lowpass){
        ierr=VecCopy(x,xflt); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }
    if (pct == 1.0 && !lowpass){
        ierr=VecSet(xflt,0.0); CHKERRQ(ierr);
        PetscFunctionReturn(ierr);
    }

    indicator[0] = 1;
    indicator[1] = 0;

    if(!lowpass){
        indicator[0] = 0;
        indicator[1] = 1;
    }

    // get local pencil size and allocation size
    nalloc=this->m_Opt->GetFFT().nalloc;

    // allocate
    if(this->m_xhat==NULL){ this->m_xhat=(ScalarTypeFD*)accfft_alloc(nalloc); }

    // get parameters
    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<int>(this->m_Opt->GetDomainPara().nx[i]);
        nxhalf[i] = static_cast<ScalarType>(nx[i]/2.0);
    }
    scale = this->m_Opt->ComputeFFTScale();

    // compute fft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_Opt->GetFFT().plan,p_x,this->m_xhat,timer);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

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
    for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1){ // x1
        for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2){ // x2
            for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                w[0] = static_cast<ScalarType>(i1 + this->m_Opt->GetFFT().ostart[0]);
                w[1] = static_cast<ScalarType>(i2 + this->m_Opt->GetFFT().ostart[1]);
                w[2] = static_cast<ScalarType>(i3 + this->m_Opt->GetFFT().ostart[2]);

                bool inside = true;
                inside = ( ( (w[0] < cfreq[0][0]) || (w[0] > cfreq[0][1]) ) && inside ) ? true : false;
                inside = ( ( (w[1] < cfreq[1][0]) || (w[1] > cfreq[1][1]) ) && inside ) ? true : false;
                inside = ( ( (w[2] < cfreq[2][0]) || (w[2] > cfreq[2][1]) ) && inside ) ? true : false;

                indic = inside ? indicator[0] : indicator[1];

                // compute linear / flat index
                li = GetLinearIndex(i1,i2,i3,this->m_Opt->GetFFT().osize);

                if ( indic == 0 ){
                    this->m_xhat[li][0] = 0.0;
                    this->m_xhat[li][1] = 0.0;
                }
                else{
                    this->m_xhat[li][0] *= scale;
                    this->m_xhat[li][1] *= scale;
                }

            } // i1
        } // i2
    } // i3

} // pragma omp parallel


    // compute inverse fft
    ierr=VecGetArray(xflt,&p_xflt); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(this->m_Opt->GetFFT().plan,this->m_xhat,p_xflt,timer);
    ierr=VecRestoreArray(xflt,&p_xflt); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(timer);

    // increment counter
    this->m_Opt->IncrementCounter(FFT,2);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);

}




/********************************************************************
 * @brief apply gaussian smoothing operator to input data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplySmoothing"
PetscErrorCode PreProcReg::ApplySmoothing(Vec xsmooth, Vec x)
{
    PetscErrorCode ierr;
    IntType nalloc;
    ScalarType *p_x=NULL,*p_xsmooth=NULL,nx[3],c[3],scale;
    double timer[5]={0,0,0,0,0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(xsmooth!=NULL,"null pointer"); CHKERRQ(ierr);

    // get local pencil size and allocation size
    nalloc=this->m_Opt->GetFFT().nalloc;

    if(this->m_xhat == NULL){
        this->m_xhat=(ScalarTypeFD*)accfft_alloc(nalloc);
    }

    // get parameters
    for (int i = 0; i < 3; ++i){

        nx[i] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[i]);

        // sigma is provided by user in # of grid points
        c[i] = this->m_Opt->GetSigma(i)*this->m_Opt->GetDomainPara().hx[i];
        c[i] *= c[i];

    }

    scale = this->m_Opt->ComputeFFTScale();

    // compute fft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_Opt->GetFFT().plan,p_x,this->m_xhat,timer);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);


#pragma omp parallel
{
    IntType i1,i2,i3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1){ // x1
        for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2){ // x2
            for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                ScalarType k1 = static_cast<ScalarType>(i1 + this->m_Opt->GetFFT().ostart[0]);
                ScalarType k2 = static_cast<ScalarType>(i2 + this->m_Opt->GetFFT().ostart[1]);
                ScalarType k3 = static_cast<ScalarType>(i3 + this->m_Opt->GetFFT().ostart[2]);

                // check if grid index is larger or smaller then
                // half of the total grid size
                bool flagx1 = (k1 <= nx[0]*0.5);
                bool flagx2 = (k2 <= nx[1]*0.5);
                bool flagx3 = (k3 <= nx[2]*0.5);

                k1 = flagx1 ? k1 : -nx[0] + k1;
                k2 = flagx2 ? k2 : -nx[1] + k2;
                k3 = flagx3 ? k3 : -nx[2] + k3;

                ScalarType sik = 0.5*( (k1*k1*c[0]) + (k2*k2*c[1]) + (k3*k3*c[2]) );
                sik = exp(-sik);

                // compute linear / flat index
                IntType li = GetLinearIndex(i1,i2,i3,this->m_Opt->GetFFT().osize);

                this->m_xhat[li][0] *= scale*sik;
                this->m_xhat[li][1] *= scale*sik;

            } // i1
        } // i2
    } // i3

} // pragma omp parallel


    // compute inverse fft
    ierr=VecGetArray(xsmooth,&p_xsmooth); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(this->m_Opt->GetFFT().plan,this->m_xhat,p_xsmooth,timer);
    ierr=VecRestoreArray(xsmooth,&p_xsmooth); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




} // end of namespace

#endif // _PREPROCREG_CPP_
