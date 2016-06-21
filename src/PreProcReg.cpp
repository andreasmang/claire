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
PetscErrorCode PreProcReg::Restrict(Vec xl, Vec x, IntType* nx_l)
{
    PetscErrorCode ierr;
    accfft_plan* plan=NULL;
    ScalarType *p_xdummy=NULL,*p_x=NULL,*p_xl=NULL,scale;
    ScalarTypeFD* p_xlhat=NULL;
    Complex *p_xhatdummy=NULL;
    IntType nalloc;
    int _nx_l[3],_ostart_l[3],_osize_l[3],_isize_l[3],_istart_l[3];
    IntType nx[3],ostart_l[3],osize_l[3];//oend_l[3];
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(xl!=NULL,"null pointer"); CHKERRQ(ierr);

    if(this->m_xhat == NULL){
        this->m_xhat=(ScalarTypeFD*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
    }

    // get grid sizes/fft scales
    scale = this->m_Opt->ComputeFFTScale();

    // parse input sizes
    _nx_l[0] = static_cast<int>(nx_l[0]);
    _nx_l[1] = static_cast<int>(nx_l[1]);
    _nx_l[2] = static_cast<int>(nx_l[2]);

    ierr=Assert(nx_l[0] < this->m_Opt->GetDomainPara().nx[0],"grid sizes in restriction wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_l[1] < this->m_Opt->GetDomainPara().nx[1],"grid sizes in restriction wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_l[2] < this->m_Opt->GetDomainPara().nx[2],"grid sizes in restriction wrong"); CHKERRQ(ierr);

    // allocate container for inverse FFT
    nalloc=accfft_local_size_dft_r2c(_nx_l,_isize_l,_istart_l,_osize_l,_ostart_l,this->m_Opt->GetFFT().mpicomm);
    p_xlhat=(ScalarTypeFD*)accfft_alloc(nalloc);

    for(int i = 0; i < 3; ++i){
        nx[i] = this->m_Opt->GetDomainPara().nx[i];
        osize_l[i] = static_cast<IntType>(_osize_l[i]);
        ostart_l[i] = static_cast<IntType>(_ostart_l[i]);
        //oend_l[i] = ostart_l[i] + osize_l[i];
    }

    // compute fft of input data
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_Opt->GetFFT().plan,p_x,this->m_xhat,ffttimers);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

#pragma omp parallel
{
    IntType li=0,lj=0;
    ScalarType k1,k2,k3,j1,j2,j3;
#pragma omp for
    for (IntType i1 = 0; i1 < osize_l[0]; ++i1){ // x1
        for (IntType i2 = 0; i2 < osize_l[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < osize_l[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                k1 = static_cast<ScalarType>(i1 + ostart_l[0]);
                k2 = static_cast<ScalarType>(i2 + ostart_l[1]);
                k3 = static_cast<ScalarType>(i3 + ostart_l[2]);

                j1 = i1;
                j2 = i2;
                j3 = i3;

                if (k1 > nx_l[0]/2) j1 = (nx[0]-1) + k1 - nx_l[0];
                if (k2 > nx_l[1]/2) j2 = (nx[1]-1) + k2 - nx_l[1];
                if (k3 > nx_l[2]/2) j3 = (nx[2]-1) + k3 - nx_l[2];

                li = GetLinearIndex(i1,i2,i3,osize_l);
                lj = GetLinearIndex(j1,j2,j3,this->m_Opt->GetFFT().osize);

                p_xlhat[li][0] = scale*this->m_xhat[lj][0];
                p_xlhat[li][1] = scale*this->m_xhat[lj][1];

            } // i1
        } // i2
    } // i3

} // pragma omp parallel


    // allocate fft
    p_xdummy = (ScalarType*)accfft_alloc(nalloc);
    p_xhatdummy = (Complex*)accfft_alloc(nalloc);
    plan=accfft_plan_dft_3d_r2c(_nx_l,p_xdummy,(double*)p_xhatdummy,this->m_Opt->GetFFT().mpicomm,ACCFFT_MEASURE);

    ierr=VecGetArray(xl,&p_xl); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(plan,p_xlhat,p_xl,ffttimers);
    ierr=VecRestoreArray(xl,&p_xl); CHKERRQ(ierr);

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    // clean up
    if(plan != NULL){ accfft_destroy_plan(plan); plan=NULL; }
    if(p_xlhat != NULL){ accfft_free(p_xlhat); p_xlhat=NULL; }
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
#define __FUNCT__ "Restrict"
PetscErrorCode PreProcReg::RestrictionGetPoints(IntType* nx_)
{
    PetscErrorCode ierr;
    int _nx_[3],_ostart_[3],_osize_[3],_isize_[3],_istart_[3];
    int rank,nprocs,ni,no,orank,c_grid[2],p1,p2;
    IntType nx[3],ostart_[3],osize_[3],osize[3],ostart[3],oend[3];
    IntType li,li_,i1_,i2_,i3_;
    ScalarType k1,k2,k3,k1_,k2_,k3_,nxhalf_[3];
    bool owned;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

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
    _nx_[0] = static_cast<int>(nx_[0]);
    _nx_[1] = static_cast<int>(nx_[1]);
    _nx_[2] = static_cast<int>(nx_[2]);

    // allocate container for inverse FFT
    accfft_local_size_dft_r2c(_nx_,_isize_,_istart_,_osize_,_ostart_,this->m_Opt->GetFFT().mpicomm);

    for(int i = 0; i < 3; ++i){

        nx[i] = this->m_Opt->GetDomainPara().nx[i];
        osize_[i] = static_cast<IntType>(_osize_[i]);
        ostart_[i] = static_cast<IntType>(_ostart_[i]);

        nxhalf_[i] = std::ceil(static_cast<ScalarType>(nx_[i])/2.0);

        ostart[i] = this->m_Opt->GetFFT().ostart[i];
        osize[i] = this->m_Opt->GetFFT().osize[i];
        oend[i] = ostart[i] + osize[i];
//        std::cout<<oend[i]<<std::endl;
    }

    c_grid[0] = this->m_Opt->GetNetworkDims(0);
    c_grid[1] = this->m_Opt->GetNetworkDims(1);

    ni=0;
    no=0;
    for (i1_ = 0; i1_ < osize_[0]; ++i1_){ // x1
        for (i2_ = 0; i2_ < osize_[1]; ++i2_){ // x2
            for (i3_ = 0; i3_ < osize_[2]; ++i3_){ // x3

                // compute coordinates (nodal grid)
                k1_ = static_cast<ScalarType>(i1_ + ostart_[0]);
                k2_ = static_cast<ScalarType>(i2_ + ostart_[1]);
                k3_ = static_cast<ScalarType>(i3_ + ostart_[2]);

                // compute index of fine grid according to quadrant
                k1 = k1_ < nxhalf_[0] ? k1_ : (nx[0]-1) + k1_ - nx_[0];
                k2 = k2_ < nxhalf_[1] ? k2_ : (nx[1]-1) + k2_ - nx_[1];
                k3 = k3_ < nxhalf_[2] ? k3_ : (nx[2]-1) + k3_ - nx_[2];


                owned=true;
                if (k1 < ostart[0] || k1 >= oend[0] ) owned = false;
                if (k2 < ostart[1] || k2 >= oend[1] ) owned = false;
                if (k3 < ostart[2] || k3 >= 2*oend[2] ) owned = false;

                li_ = GetLinearIndex(k1_,k2_,k3_,nx_);
                li  = GetLinearIndex(k1,k2,k3,nx);

                if (owned){
                    this->m_IndicesC[rank].push_back(li_);
                    this->m_IndicesF[rank].push_back(li);
                    ++ni;
                }
                else{

                    // compute processor id
                    p1=static_cast<int>(std::ceil(k1/static_cast<ScalarType>(osize[0])));
                    p2=static_cast<int>(std::ceil(k2/static_cast<ScalarType>(osize[1])));

                    orank = p1*c_grid[0]+p2;
//                    std::cout<<k1<<" "<<k2<<" "<<k3<<" "<<std::endl;
//                    ierr=Assert(orank < nprocs,"rank larger than number of procs"); CHKERRQ(ierr);
//                    this->m_IndicesC[orank].push_back(li_);
//                    this->m_IndicesF[orank].push_back(li);
                    ++no;

                }

            } // i1
        } // i2
    } // i3


    std::cout<< rank << " " << no << " " << ni <<std::endl;

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
#define __FUNCT__ "ApplyGaussianSmoothing"
PetscErrorCode PreProcReg::ApplyGaussianSmoothing(Vec y, Vec x)
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
