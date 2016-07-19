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

    this->m_IndicesComputed=false;
    this->m_ResetGridChangeOperators=false;
    this->m_GridChangeOperatorsSet=false;

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

    this->m_NumSendF=NULL;
    this->m_NumRecvF=NULL;

    this->m_OffsetSendF=NULL;
    this->m_OffsetRecvF=NULL;

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

    if (this->m_NumSendF==NULL){
        delete [] this->m_NumSendF;
        this->m_NumSendF=NULL;
    }
    if (this->m_NumRecvF==NULL){
        delete [] this->m_NumRecvF;
        this->m_NumRecvF=NULL;
    }

    if (this->m_OffsetSendF==NULL){
        delete [] this->m_OffsetSendF;
        this->m_OffsetSendF=NULL;
    }
    if (this->m_OffsetRecvF==NULL){
        delete [] this->m_OffsetRecvF;
        this->m_OffsetRecvF=NULL;
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
#define __FUNCT__ "SetupGridChangeOperators"
PetscErrorCode PreProcReg::SetupGridChangeOperators(IntType* nx_f,IntType* nx_c)
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
    this->m_GridChangeOperatorsSet=true;

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
    ScalarType *p_xf=NULL,*p_xc=NULL,scale;
    IntType n;
    int rank;
    double timer[5]={0,0,0,0,0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("applying restriction operator"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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

    if (this->m_ResetGridChangeOperators){
        this->m_GridChangeOperatorsSet=false;
    }

    // set up fft operators
    if (!this->m_GridChangeOperatorsSet){
        ierr=this->SetupGridChangeOperators(nx_f,nx_c); CHKERRQ(ierr);
    }


    n  = this->m_osizeC[0];
    n *= this->m_osizeC[1];
    n *= this->m_osizeC[2];

#pragma omp parallel
{
    IntType l;
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
    if (this->m_IndicesComputed==false){
        ierr=this->ComputeIndices(nx_f,nx_c); CHKERRQ(ierr);
    }

    // get grid sizes/fft scales
    scale = this->m_FFTFineScale;

    // get number of entries we are going to assign
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
#define __FUNCT__ "ComputeIndices"
PetscErrorCode PreProcReg::ComputeIndices(IntType* nx_f,IntType* nx_c)
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


    ierr=this->CommunicateData(); CHKERRQ(ierr);
    ierr=this->AssignFourierCoeff(); CHKERRQ(ierr);

    this->m_IndicesComputed = true;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup for applying restriction operator
 * @param nx_c grid size on coarse grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "AssignFourierCoeff"
PetscErrorCode PreProcReg::AssignFourierCoeff()
{
    PetscErrorCode ierr=0;
    int merr,nprocs,rank,i_recv,i_send;
    IntType n,k_f[3],l,i_f[3],nalloc_send,nalloc_recv,os_send,os_recv,nr,ns,nxhalf_c[3];
    MPI_Request *send_request=NULL,*recv_request=NULL;
    MPI_Status status;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    if (this->m_OffsetSendF==NULL){
        try{this->m_OffsetSendF = new IntType[nprocs];}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_OffsetRecvF==NULL){
        try{this->m_OffsetRecvF = new IntType[nprocs];}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_NumSendF==NULL){
        try{this->m_NumSendF = new IntType[nprocs];}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if(this->m_NumRecvF==NULL){
        try{this->m_NumRecvF = new IntType[nprocs];}
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
        this->m_NumSendF[p] = 0;
        this->m_NumRecvF[p] = 0;
        if (!this->m_IndicesF[p].empty()){
            this->m_NumSendF[p] = this->m_IndicesF[p].size()/3;
        }

    }

    for(int i = 0; i < 3; ++i){
        nxhalf_c[i] = static_cast<IntType>(std::ceil(static_cast<ScalarType>(this->m_nxC[i])/2.0));
    }

    // communicate the amount of data we will send from one
    // processor to another
    MPI_Alltoall(this->m_NumSendF,1,MPI_INT,this->m_NumRecvF,1,MPI_INT,PETSC_COMM_WORLD);

    nalloc_send=this->m_NumSendF[0];
    nalloc_recv=this->m_NumRecvF[0];
    this->m_OffsetSendF[0] = 0;
    this->m_OffsetRecvF[0] = 0;
    for (int p = 1; p < nprocs; ++p){
        this->m_OffsetSendF[p] = this->m_OffsetSendF[p-1] + this->m_NumSendF[p-1];
        this->m_OffsetRecvF[p] = this->m_OffsetRecvF[p-1] + this->m_NumRecvF[p-1];

        nalloc_send += this->m_NumSendF[p];
        nalloc_recv += this->m_NumRecvF[p];
    }
    std::cout<< "assigning fourier coeffiecients" <<std::endl;

    // if we actually need to allocate something
    if (nalloc_send > 0){

        if (this->m_FourierCoeffSendF==NULL){
            try{this->m_FourierCoeffSendF = new ScalarType[nalloc_send*2];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        if (this->m_FourierIndicesSendF==NULL){
            try{this->m_FourierIndicesSendF = new IntType[nalloc_send*3];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        for (int p = 0; p < nprocs; ++p){

            if (this->m_NumSendF[p] != 0){

                n = this->m_IndicesF[p].size()/3;
                if (n != this->m_NumSendF[p]){
                    ss << "size mismatch " << n << "!=" << this->m_NumSendF[p];
                    ierr=ThrowError(ss.str()); CHKERRQ(ierr);
                }

                for (int j = 0; j < n; ++j){
                    for (int i = 0; i < 3; ++i){

                        k_f[i] = this->m_IndicesF[p][j*3 + i];

                        // convert global index into local index
                        i_f[i] = k_f[i] - this->m_ostartF[i];

                        // check if we're inside expected range
                        if ( (i_f[i] >= this->m_osizeF[i]) || (i_f[i] < 0) ){
                            std::cout<<" r "<<rank<<" "<<i_f[i]<<">="<<this->m_osizeF[i]<<"   "<<i_f[i]<<"<0"<< std::endl;
                        }

                    }

                    // compute flat index
                    l = GetLinearIndex(i_f[0],i_f[1],i_f[2],this->m_osizeF);

                    os_send = this->m_OffsetSendF[p];
                    os_recv = this->m_OffsetRecvF[p];

                    // assign values to coarse grid
                    this->m_FourierCoeffSendF[2*j + 0 + 2*os_send] = this->m_XHatFine[l][0];
                    this->m_FourierCoeffSendF[2*j + 1 + 2*os_send] = this->m_XHatFine[l][1];

                    this->m_FourierIndicesSendF[3*j + 0 + 3*os_send] = k_f[0];
                    this->m_FourierIndicesSendF[3*j + 1 + 3*os_send] = k_f[1];
                    this->m_FourierIndicesSendF[3*j + 2 + 3*os_send] = k_f[2];


                } // for all points
            } // if indices are not empty
        } // for all procs
    }

    if (nalloc_recv > 0){
        if (this->m_FourierCoeffRecvF==NULL){
            try{this->m_FourierCoeffRecvF = new ScalarType[nalloc_recv*2];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        if (this->m_FourierIndicesRecvF==NULL){
            try{this->m_FourierIndicesRecvF = new IntType[nalloc_recv*3];}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        for (int i = 0; i < 2*nalloc_recv; ++i){
            this->m_FourierCoeffRecvF[i] = -1;
        }

        for (int i = 0; i < 3*nalloc_recv; ++i){
            this->m_FourierIndicesRecvF[i] = -1;
        }

    }

    // send and recv fourier coefficients on fine grid
    for (int i = 0; i < nprocs; ++i){

        i_send=i; i_recv=i;
        send_request[i_send] = MPI_REQUEST_NULL;
        recv_request[i_recv] = MPI_REQUEST_NULL;

        os_send = this->m_OffsetSendF[i];
        ns = this->m_NumSendF[i];
        if (ns > 0){
            merr=MPI_Isend(&this->m_FourierCoeffSendF[2*os_send],2*ns,MPI_DOUBLE,i_send,0,PETSC_COMM_WORLD,&send_request[i_send]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }

        os_recv = this->m_OffsetRecvF[i];
        nr = this->m_NumRecvF[i];
        if (nr > 0){
            merr=MPI_Irecv(&this->m_FourierCoeffRecvF[2*os_recv],2*nr,MPI_DOUBLE,i_recv,0,PETSC_COMM_WORLD,&recv_request[i_recv]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    for (int i=0; i < nprocs; ++i){
        if(send_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&send_request[i], &status); }
        if(recv_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&recv_request[i], &status); }
    }


    // send and recv indices on fine grid
    for (int i = 0; i < nprocs; ++i){

        i_send=i; i_recv=i;
        send_request[i_send] = MPI_REQUEST_NULL;
        recv_request[i_recv] = MPI_REQUEST_NULL;

        os_send = this->m_OffsetSendF[i];
        ns = this->m_NumSendF[i];
        if (ns > 0){
            merr=MPI_Isend(&this->m_FourierIndicesSendF[3*os_send],3*ns,MPI_INT,i_send,0,PETSC_COMM_WORLD,&send_request[i_send]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }

        os_recv = this->m_OffsetRecvF[i];
        nr = this->m_NumRecvF[i];
        if (nr > 0){
            merr=MPI_Irecv(&this->m_FourierIndicesRecvF[3*os_recv],3*nr,MPI_INT,i_recv,0,PETSC_COMM_WORLD,&recv_request[i_recv]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }
    }

    for (int i=0; i < nprocs; ++i){
        if(send_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&send_request[i], &status); }
        if(recv_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&recv_request[i], &status); }
    }

    IntType k_c[3];
    for (int p = 0; p < nprocs; ++p){
        nr = this->m_NumRecvF[p];
        os_recv = this->m_OffsetRecvF[p];
        for (int j = 0; j < nr; ++j){
            bool outofbounds=false;
            for (int i = 0; i < 3; ++i){
                k_f[i] = this->m_FourierIndicesRecvF[3*j + i + 3*os_recv];

                // get wave number index on coarse grid from index on fine grid
                k_c[i] = k_f[i] <= nxhalf_c[i] ? k_f[i] : this->m_nxC[i] - this->m_nxF[i] + k_f[i];

                if ( (k_c[i] < this->m_ostartC[i]) || (k_c[i] > this->m_ostartC[i] + this->m_osizeC[i]) ){
                    outofbounds=true;
                }
            }
            if (outofbounds) std::cout << k_f[0] << " " << k_f[1] << " " << k_f[2] << std::endl;
        }
    }


    if (send_request!=NULL) { delete [] send_request; send_request=NULL; }
    if (recv_request!=NULL) { delete [] recv_request; recv_request=NULL; }
//    if (n_offset_send!=NULL) { delete [] n_offset_send; n_offset_send=NULL; }
//    if (n_offset_recv!=NULL) { delete [] n_offset_recv; n_offset_recv=NULL; }
//    if (n_send!=NULL) { delete [] n_send; n_send=NULL; }
//    if (n_recv!=NULL) { delete [] n_recv; n_recv=NULL; }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief do setup for applying restriction operator
 * @param nx_c grid size on coarse grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CommunicateData"
PetscErrorCode PreProcReg::CommunicateData()
{
    PetscErrorCode ierr=0;
    int rank,nprocs,merr,i_send,i_recv,ns,nr,os_send,os_recv;
    int *n_send=NULL,*n_recv=NULL,*offset_send=NULL,*offset_recv=NULL;
    IntType nalloc_recv,nalloc_send;
    double *send_val=NULL,*recv_val=NULL;
    MPI_Status status;
    MPI_Request *send_request=NULL,*recv_request=NULL;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);


    try{n_send = new IntType[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{n_recv = new IntType[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{offset_send = new IntType[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{offset_recv = new IntType[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{send_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{recv_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    for (int i = 0; i < nprocs; ++i){
        n_send[i] = this->m_IndicesF[i].size()/3;
    }

    // communicate the amount of data we will send from one
    // processor to another
    MPI_Alltoall(n_send,1,MPI_INT,n_recv,1,MPI_INT,PETSC_COMM_WORLD);

    // compute total number of values we'll send
    // and receive as well as offset within arrays
    offset_send[0] = 0;
    offset_recv[0] = 0;
    nalloc_recv = n_recv[0];
    nalloc_send = n_send[0];
    for (int i = 1; i < nprocs; ++i){
        offset_send[i] = offset_send[i-1] + n_send[i-1];
        offset_recv[i] = offset_recv[i-1] + n_recv[i-1];
        nalloc_recv += n_recv[i];
        nalloc_send += n_send[i];
    }
//    std::cout<<" alloc recv "<< nalloc_recv<<std::endl;
//    std::cout<<" alloc send "<< nalloc_send<<std::endl;

    // allocate container to hold values (send and receive)
    try{send_val = new double[nalloc_send];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{recv_val = new double[nalloc_recv];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }


    //if (rank == 0) std::cout<< "starting communication" << std::endl;
    for (int i = 0; i < nalloc_send; ++i){
        send_val[i] = rank;
    }

    for (int i = 0; i < nalloc_recv; ++i){
        recv_val[i] = -1;
    }


    for (int i = 0; i < nprocs; ++i){

        i_send=i; i_recv=i;
        send_request[i_send] = MPI_REQUEST_NULL;
        recv_request[i_recv] = MPI_REQUEST_NULL;

        os_send = offset_send[i];
        ns = n_send[i];
        if (ns > 0){
            merr=MPI_Isend(&send_val[os_send],ns,MPI_DOUBLE,i_send,0,PETSC_COMM_WORLD,&send_request[i_send]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }

        os_recv = offset_recv[i];
        nr = n_recv[i];
        if (nr > 0){
            merr=MPI_Irecv(&recv_val[os_recv],nr,MPI_DOUBLE,i_recv,0,PETSC_COMM_WORLD,&recv_request[i_recv]);
            ierr=MPIERRQ(merr); CHKERRQ(ierr);
        }
    }


    for (int i=0; i < nprocs; ++i){
        if(recv_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&recv_request[i], &status); }
        if(send_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&send_request[i], &status); }
    }

/*  if (rank == 0){
        std::cout << " rank " << rank << std::endl;
        for (int i=0; i<nalloc_recv; ++i){ std::cout << std::setprecision(1) << (int)recv_val[i] << std::endl; }
    }
    if (rank == 1){
        std::cout << " rank " << rank << std::endl;
        for (int i=0; i<nalloc_recv; ++i){ std::cout << std::setprecision(1) << (int)recv_val[i] << std::endl; }
    } */

    if(n_send!=NULL) delete [] n_send;
    if(n_recv!=NULL) delete [] n_recv;

    if(recv_val!=NULL) delete [] recv_val;
    if(send_val!=NULL) delete [] send_val;

    if(offset_send!=NULL) delete [] offset_send;
    if(offset_recv!=NULL) delete [] offset_recv;

    if(send_request!=NULL) delete [] send_request;
    if(recv_request!=NULL) delete [] recv_request;

//    if (rank == 0) std::cout<< "communication done" << std::endl;

    PetscFunctionReturn(ierr);

}





/********************************************************************
 * @brief do setup for applying restriction operator
 * @param nx_c grid size on coarse grid
 *******************************************************************/
/*
#undef __FUNCT__
#define __FUNCT__ "CommunicateData"
PetscErrorCode PreProcReg::CommunicateData()
{
    PetscErrorCode ierr=0;
    int nprocs,rank,nrecvalloc,nsendalloc,cr,cs,ns,nr,i_recv,i_send,merr;
    //IntType *n_send=NULL,*n_recv=NULL,*send_offset=NULL,*recv_offset=NULL;
    IntType n_send[2],n_recv[2],send_offset[2],recv_offset[2];
    //ScalarType *recv_val=NULL,*send_val=NULL;
    ScalarType recv_val[20],send_val[20];
    MPI_Status status;
    MPI_Request *send_request=NULL,*recv_request=NULL;
    bool isempty;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    // figure out how many points each process owns that need to be sent

    try{n_send = new IntType[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{n_recv = new IntType[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    // for all processors
    for (int i=0; i < nprocs; ++i){

        // check if list is empty
//        isempty = this->m_IndicesC[i].empty();

//        if ( isempty ) n_send[i] = 0;
//        else           n_send[i] = this->m_IndicesC[i].size();
        n_send[i] = 10;
    }

    // communicate
    MPI_Alltoall(n_send,1,MPI_INT,n_recv,1,MPI_INT,PETSC_COMM_WORLD);
    if (nprocs > 1){
        if (rank==0){
            std::cout << " rank: " << rank << std::endl;
            std::cout << " [0] " << n_send[0] << " [1] " << n_send[1]  << std::endl;
            std::cout << " [0] " << n_recv[0] << " [1] " << n_recv[1] << std::endl;
        }
        if (rank==1){
            std::cout << " rank: " << rank << std::endl;
            std::cout << " [0] " << n_send[0] << " [1] " << n_send[1] << std::endl;
            std::cout << " [0] " << n_recv[0] << " [1] " << n_recv[1] << std::endl;
        }
    }
    else{
            std::cout << " rank: " << rank << std::endl;
            std::cout << " [0] " << n_send[0] << std::endl;
            std::cout << " [0] " << n_recv[0] << std::endl;
    }
    // compute how many points we receive (this includes the points
    // we also have on this node)
    nrecvalloc=0; nsendalloc=0;
    for (int i=0; i < nprocs; ++i){
        nrecvalloc += n_recv[i];
        nsendalloc += n_send[i];
    }

//    std::cout << " allocation size recv: " << nrecvalloc << std::endl;
//    std::cout << " allocation size send: " << nsendalloc << std::endl;
    try{recv_val = new ScalarType[nrecvalloc];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    try{send_val = new ScalarType[nsendalloc];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{send_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    try{recv_request = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{send_offset = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{recv_offset = new MPI_Request[nprocs];}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    recv_offset[0] = 0;
    send_offset[0] = 0;
    for (int i = 1; i < nprocs; ++i){
       recv_offset[i] = recv_offset[i-1] + n_recv[i-1];
       send_offset[i] = send_offset[i-1] + n_send[i-1];
    }


    for (int i = 0; i < nrecvalloc; ++i){
        recv_val[i] = 10 + rank;
    }

    for (int i = 0; i < nsendalloc; ++i){
        send_val[i] = rank;
    }

    // send all indices
    for (int i=0; i < nprocs; ++i){
        i_recv = i;
        i_send = i;

        recv_request[i_recv] = MPI_REQUEST_NULL;
        send_request[i_send] = MPI_REQUEST_NULL;

        ns = n_send[i_send]; // number to send
        cs = send_offset[i];
        //std::cout<<ns<<std::endl;
        if (ns != 0){
//            merr=MPI_Isend(&send_val[cs],ns,MPI_DOUBLE,i_send,MPI_ANY_TAG,PETSC_COMM_WORLD,&send_request[i_send]);
//            ierr=MPIERRQ(merr); CHKERRQ(ierr);
            //merr=MPI_Send(&send_val,ns,MPI_DOUBLE,i_send,MPI_ANY_TAG,PETSC_COMM_WORLD);
        }

        nr = n_recv[i_recv]; // number to receive
        cr = recv_offset[i];
        //std::cout<<nr<<std::endl;

        if (nr != 0){
//            merr=MPI_Irecv(&recv_val[cr],nr,MPI_DOUBLE,i_recv,MPI_ANY_TAG,PETSC_COMM_WORLD,&recv_request[i_recv]);
//            ierr=MPIERRQ(merr); CHKERRQ(ierr);
            //merr=MPI_Recv(&recv_val,nr,MPI_DOUBLE,i_recv,MPI_ANY_TAG,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
        }

    }

    merr=MPI_Isend(&send_val,10,MPI_DOUBLE,rank,MPI_ANY_TAG,PETSC_COMM_WORLD,&send_request[rank]);
    ierr=MPIERRQ(merr); CHKERRQ(ierr);
    merr=MPI_Irecv(&recv_val,10,MPI_DOUBLE,rank,MPI_ANY_TAG,PETSC_COMM_WORLD,&recv_request[rank]);
    ierr=MPIERRQ(merr); CHKERRQ(ierr);

    for (int i=0; i < nprocs; ++i){
//        if(recv_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&recv_request[i], &status); }
//        if(send_request[i]!=MPI_REQUEST_NULL) { MPI_Wait(&send_request[i], &status); }
    }

    //std::cout<< "communication done " << rank << std::endl;



//    if(n_send!=NULL) delete [] n_send;
//    if(n_recv!=NULL) delete [] n_recv;
    if(recv_request!=NULL) delete [] recv_request;
    if(send_request!=NULL) delete [] send_request;
//    if(recv_val!=NULL) delete [] recv_val;
//    if(send_val!=NULL) delete [] send_val;

    PetscFunctionReturn(ierr);
}
*/


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
    int rank;
    IntType n;
    ScalarType *p_xf=NULL,*p_xc=NULL,scale;
    double timer[5]={0,0,0,0,0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("applying prolongation operator"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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

    if (this->m_ResetGridChangeOperators){
        this->m_GridChangeOperatorsSet=false;
    }

    // set up fft operators
    if (!this->m_GridChangeOperatorsSet){
        ierr=this->SetupGridChangeOperators(nx_f,nx_c); CHKERRQ(ierr);
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

    // get size
    if (this->m_IndicesComputed==false){
        ierr=this->ComputeIndices(nx_f,nx_c); CHKERRQ(ierr);
    }

    // get grid sizes/fft scales
    scale = this->m_FFTCoarseScale;

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
