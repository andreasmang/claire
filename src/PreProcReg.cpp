#ifndef _PREPROCREG_CPP_
#define _PREPROCREG_CPP_


#include "PreProcReg.hpp"


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
    this->m_ReadWrite = NULL;
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

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set io interface for data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetIO"
PetscErrorCode PreProcReg::SetReadWrite(PreProcReg::ReadWriteType* io)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr=Assert(io != NULL,"null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite=io;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief prolong data
 * @param x input vector
 * @param y output vector y = P[x]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Prolong"
PetscErrorCode PreProcReg::Prolong(Vec y, Vec x, IntType* nxpro)
{
    PetscErrorCode ierr;

    ierr=Assert(y!=NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(x!=NULL, "null pointer"); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief restrict data
 * @param x input vector
 * @param y output vector xl = R[x]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Restrict"
PetscErrorCode PreProcReg::Restrict(Vec xcoarse, Vec x, IntType* nx_c)
{
    PetscErrorCode ierr;
    accfft_plan* plan=NULL;
    ScalarType *p_xdummy=NULL,*p_x=NULL,*p_xcoarse=NULL,scale;
    ScalarTypeFD* p_xcoarsehat=NULL;
    Complex *p_xhatdummy=NULL;
    IntType nalloc,n;
    int _nx_c[3],_ostart_c[3],_osize_c[3],_isize_c[3],_istart_c[3],rank;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("applying restriction operator"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(xcoarse!=NULL,"null pointer"); CHKERRQ(ierr);

    // everything we do not set is zero
    ierr=VecSet(xcoarse,0.0); CHKERRQ(ierr);

    // restriction operator
    ierr=this->RestrictionGetPoints(nx_c); CHKERRQ(ierr);

    if(this->m_xhat == NULL){
        this->m_xhat=(ScalarTypeFD*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
    }

    // get grid sizes/fft scales
    scale = this->m_Opt->ComputeFFTScale();

    // parse input sizes
    _nx_c[0] = static_cast<int>(nx_c[0]);
    _nx_c[1] = static_cast<int>(nx_c[1]);
    _nx_c[2] = static_cast<int>(nx_c[2]);

    // allocate container for inverse FFT
    nalloc=accfft_local_size_dft_r2c(_nx_c,_isize_c,_istart_c,_osize_c,_ostart_c,this->m_Opt->GetFFT().mpicomm);
    p_xcoarsehat=(ScalarTypeFD*)accfft_alloc(nalloc);

    // compute fft of input data
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_Opt->GetFFT().plan,p_x,this->m_xhat,ffttimers);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

    n = this->m_IndicesC[rank].size();
    ierr=Assert(n == this->m_IndicesF[rank].size(),"size error"); CHKERRQ(ierr);

#pragma omp parallel
{
    IntType li,li_c;
#pragma omp for
    for (IntType i = 0; i < n; ++i){

        li   = this->m_IndicesF[rank][i];
        li_c = this->m_IndicesC[rank][i];

        p_xcoarsehat[li_c][0] = scale*this->m_xhat[li][0];
        p_xcoarsehat[li_c][1] = scale*this->m_xhat[li][1];

    }

} // pragma omp parallel

    // allocate fft
    p_xdummy = (ScalarType*)accfft_alloc(nalloc);
    p_xhatdummy = (Complex*)accfft_alloc(nalloc);
    plan=accfft_plan_dft_3d_r2c(_nx_c,p_xdummy,(double*)p_xhatdummy,this->m_Opt->GetFFT().mpicomm,ACCFFT_MEASURE);

    ierr=VecGetArray(xcoarse,&p_xcoarse); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(plan,p_xcoarsehat,p_xcoarse,ffttimers);
    ierr=VecRestoreArray(xcoarse,&p_xcoarse); CHKERRQ(ierr);

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    // clean up
    if(plan != NULL){ accfft_destroy_plan(plan); plan=NULL; }
    if(p_xcoarsehat != NULL){ accfft_free(p_xcoarsehat); p_xcoarsehat=NULL; }
    if(p_xdummy != NULL){ accfft_free(p_xdummy); p_xdummy=NULL; }
    if(p_xhatdummy != NULL){ accfft_free(p_xhatdummy); p_xhatdummy=NULL; }

    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief restrict data
 * @param x input vector
 * @param y output vector y = R[x]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RestrictionGetPoints"
PetscErrorCode PreProcReg::RestrictionGetPoints(IntType* nx_c)
{
    PetscErrorCode ierr;
    int _nx_c[3],_ostart_c[3],_osize_c[3],_isize_c[3],_istart_c[3];
    int rank,nprocs,nowned,nsend,nprocessed,xrank,c_grid[2],p1,p2;
    IntType nx[3],ostart_c[3],osize_c[3],oend_c[3],osize[3],ostart[3];
    IntType li,li_c, i1,i2,i3, i1_c,i2_c,i3_c;
    ScalarType k1,k2,k3, k1_c,k2_c,k3_c, nxhalf_c[3];
    bool owned;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr=Assert(nx_c[0] < this->m_Opt->GetDomainPara().nx[0],"grid size in restriction wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_c[1] < this->m_Opt->GetDomainPara().nx[1],"grid size in restriction wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_c[2] < this->m_Opt->GetDomainPara().nx[2],"grid size in restriction wrong"); CHKERRQ(ierr);

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

    // parse input sizes (we have to use int, to parse
    // the arguments to accfft)
    _nx_c[0] = static_cast<int>(nx_c[0]); // get coarse grid size
    _nx_c[1] = static_cast<int>(nx_c[1]); // get coarse grid size
    _nx_c[2] = static_cast<int>(nx_c[2]); // get coarse grid size

    // get sizes for coarse grid size
    accfft_local_size_dft_r2c(_nx_c,_isize_c,_istart_c,_osize_c,_ostart_c,this->m_Opt->GetFFT().mpicomm);

    for(int i = 0; i < 3; ++i){

        // assign sizes for coarse grid
        osize_c[i]  = static_cast<IntType>(_osize_c[i]);
        ostart_c[i] = static_cast<IntType>(_ostart_c[i]);
        nxhalf_c[i] = std::floor(static_cast<ScalarType>(nx_c[i])/2.0);
        oend_c[i]   = ostart_c[i] + osize_c[i];

        // assign sizes for fine grid
        nx[i]     = this->m_Opt->GetDomainPara().nx[i];
        ostart[i] = this->m_Opt->GetFFT().ostart[i];
        osize[i]  = this->m_Opt->GetFFT().osize[i];
    }

    // get cartesian grid (MPI)
    c_grid[0] = this->m_Opt->GetNetworkDims(0);
    c_grid[1] = this->m_Opt->GetNetworkDims(1);

    IntType sizex2 = std::ceil( static_cast<ScalarType>(nx_c[1])/static_cast<ScalarType>(c_grid[0]));
    IntType sizex3 = std::ceil( (static_cast<ScalarType>(nx_c[2])/2.0 + 1.0)/static_cast<ScalarType>(c_grid[1]));


    nowned=0; nsend=0; nprocessed=0;
    for (i1 = 0; i1 < osize[0]; ++i1){ // x1
        for (i2 = 0; i2 < osize[1]; ++i2){ // x2
            for (i3 = 0; i3 < osize[2]; ++i3){ // x3

                // compute wave number index on fine grid
                k1 = static_cast<IntType>(i1 + ostart[0]);
                k2 = static_cast<IntType>(i2 + ostart[1]);
                k3 = static_cast<IntType>(i3 + ostart[2]);

                // only if current fourier entry is represented in
                // spectral domain of coarse grid
                if (  ( k1 <= nxhalf_c[0] || k1 > (nx[0] - nxhalf_c[0]) )
                   && ( k2 <= nxhalf_c[1] || k2 > (nx[1] - nxhalf_c[1]) )
                   && ( k3 <= nxhalf_c[2] || k3 > (nx[2] - nxhalf_c[2]) ) ){

                    ++nprocessed;
                    //std::cout << "(" << i1 << "," << i2 << "," << i3 << ")   ";

                    // get wave number index on coarse grid
                    k1_c = k1 <= nxhalf_c[0] ? k1 : nx_c[0] - nx[0] + k1;
                    k2_c = k2 <= nxhalf_c[1] ? k2 : nx_c[1] - nx[1] + k2;
                    k3_c = k3 <= nxhalf_c[2] ? k3 : nx_c[2] - nx[2] + k3;

                    if ( (k1_c < 0.0) || (k2_c < 0.0) || (k3_c < 0.0) ){
                        std::cout<<" index out of bounds (smaller than zero)"<<std::endl;
                    }
                    if ( (k1_c > nx_c[0]) || (k2_c > nx_c[1]) || (k3_c > nx_c[2]) ){
                        std::cout<<" index out of bounds (larger than nx)"<<std::endl;
                    }

                    owned=true;
                    if ( (k1_c < ostart_c[0]) || (k1_c >= oend_c[0]) ) owned = false;
                    if ( (k2_c < ostart_c[1]) || (k2_c >= oend_c[1]) ) owned = false;
                    if ( (k3_c < ostart_c[2]) || (k3_c >= oend_c[2]) ) owned = false;

                    // compute processor id
                    p1=static_cast<int>(k2_c/sizex2);
                    p2=static_cast<int>(k3_c/sizex3);

                    xrank = p1*c_grid[1] + p2;

                    if ( owned ){

                        // compute local index
                        i1_c = static_cast<IntType>(k1_c) - ostart_c[0];
                        i2_c = static_cast<IntType>(k2_c) - ostart_c[1];
                        i3_c = static_cast<IntType>(k3_c) - ostart_c[2];

                        if ( (i1_c >= osize_c[0]) || (i2_c >= osize_c[1]) || (i3_c >= osize_c[2]) ){
                            std::cout<<" index out of bounds (larger than osize)"<<std::endl;
                        }
                        if ( (i1_c < 0.0) || (i2_c < 0.0) || (i3_c < 0.0) ){
                            std::cout<<" index out of bounds (larger than osize)"<<std::endl;
                        }
                        // compute flat index
                        li_c = GetLinearIndex(i1_c,i2_c,i3_c,osize_c);
                        li   = GetLinearIndex(i1,i2,i3,osize);

                        // assign computed indices to array (for given rank)
                        this->m_IndicesC[rank].push_back(li_c);
                        this->m_IndicesF[rank].push_back(li);

                        // check if woned is really owned
                        if ( rank != xrank ){
                            std::cout<<" rank not owned: " << rank << " " << xrank <<std::endl;
                        }
                        ++nowned;
                    }
                    else{

                        if ( rank == xrank ){
                            std::cout<<" rank owned: " << rank << " " << xrank <<std::endl;
                        }

                        this->m_IndicesC[xrank].push_back(0);
                        this->m_IndicesF[xrank].push_back(0);

                        ++nsend;

                    }
                }

            } // i1
        } // i2
    } // i3

    // TODO: send data that does not belong to current proc to
    // other procs

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief prolong data
 * @param x input vector field
 * @param y output vector field y = P[x]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Prolong"
PetscErrorCode PreProcReg::Prolong(VecField* y, VecField* x, IntType* nxpro)
{
    PetscErrorCode ierr;

    ierr=Assert(y!=NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(x!=NULL, "null pointer"); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief prolong data
 * @param x input vector field
 * @param y output vector field y = P[x]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Prolong"
PetscErrorCode PreProcReg::Restrict(VecField* y, VecField* x, IntType* nxpro)
{
    PetscErrorCode ierr;

    ierr=Assert(y!=NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(x!=NULL, "null pointer"); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}



/********************************************************************
 * @brief apply gaussian smoothing operator to input data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplySmoothing"
PetscErrorCode PreProcReg::ApplySmoothing(Vec y, Vec x)
{
    PetscErrorCode ierr;
    int isize[3],osize[3],istart[3],ostart[3],n[3];
    IntType iosize[3];
    IntType nalloc;
    ScalarType hx[3],nx[3],c[3],scale;
    ScalarType *p_x=NULL, *p_y=NULL;
    double ffttimers[5]={0,0,0,0,0};
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    n[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
    n[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
    n[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

    // get local pencil size and allocation size
    nalloc=accfft_local_size_dft_r2c_t<ScalarType>(n,isize,istart,osize,ostart,
                                                   this->m_Opt->GetFFT().mpicomm);
    if(this->m_xhat == NULL){
        this->m_xhat=(ScalarTypeFD*)accfft_alloc(nalloc);
    }

    for (int i = 0; i < 3; ++i){

        hx[i] = this->m_Opt->GetDomainPara().hx[i];
        nx[i] = static_cast<ScalarType>(n[i]);

        // sigma is provided by user in # of grid points
        c[i] = this->m_Opt->GetSigma(i)*hx[i];
        c[i] *= c[i];

        iosize[i] = osize[i];
    }

    scale = this->m_Opt->ComputeFFTScale();

    // compute fft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_Opt->GetFFT().plan,p_x,this->m_xhat,ffttimers);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);


#pragma omp parallel
{
#pragma omp for
    for (IntType i1 = 0; i1 < iosize[0]; ++i1){ // x1
        for (IntType i2 = 0; i2 < iosize[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < iosize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                ScalarType k1 = static_cast<ScalarType>(i1 + ostart[0]);
                ScalarType k2 = static_cast<ScalarType>(i2 + ostart[1]);
                ScalarType k3 = static_cast<ScalarType>(i3 + ostart[2]);

                // check if grid index is larger or smaller then
                // half of the total grid size
                bool flagx1 = (k1 <= nx[0]*0.5);
                bool flagx2 = (k2 <= nx[1]*0.5);
                bool flagx3 = (k3 <= nx[2]*0.5);

                k1 = flagx1 ? k1 : -nx[0] + k1;
                k2 = flagx2 ? k2 : -nx[1] + k2;
                k3 = flagx3 ? k3 : -nx[2] + k3;

                ScalarType sik = 0.5*((k1*k1*c[0]) + (k2*k2*c[1]) + (k3*k3*c[2]));
                sik = exp(-sik);

                // compute linear / flat index
                IntType li = GetLinearIndex(i1,i2,i3,iosize);

                this->m_xhat[li][0] *= scale*sik;
                this->m_xhat[li][1] *= scale*sik;

            } // i1
        } // i2
    } // i3

} // pragma omp parallel

    // compute inverse fft
    ierr=VecGetArray(y,&p_y); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(this->m_Opt->GetFFT().plan,this->m_xhat,p_y,ffttimers);
    ierr=VecRestoreArray(y,&p_y); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);

}




} // end of namespace

#endif // _PREPROCREG_CPP_
