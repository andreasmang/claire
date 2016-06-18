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
    this->m_Kxhat = NULL;

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
    if(this->m_Kxhat!=NULL){
        accfft_free(this->m_Kxhat);
        this->m_Kxhat = NULL;
    }

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set io interface for data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetIO"
PetscErrorCode PreProcReg::SetIO(PreProcReg::ReadWriteType* io)
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
 * @param y output vector y = R[x]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Restrict"
PetscErrorCode PreProcReg::Restrict(Vec y, Vec x, IntType* nxres)
{
    PetscErrorCode ierr;
    accfft_plan* fftplan=NULL;
    ScalarType *u=NULL,*p_x=NULL,*p_y=NULL;
    FFTScaType* p_yhat=NULL;
    Complex *uk=NULL;
    IntType alloc_max;
    int nx[3],ostart[3],osize[3],isize[3],istart[3],cgrid[2];
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(y!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    if(this->m_xhat == NULL){
        this->m_xhat=(FFTScaType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);
    }

    cgrid[0] = this->m_Opt->GetNetworkDims(0);
    cgrid[1] = this->m_Opt->GetNetworkDims(1);

    // parse input sizes
    nx[0] = nxres[0];
    nx[1] = nxres[1];
    nx[2] = nxres[2];

    // allocate container for inverse FFT
    alloc_max=accfft_local_size_dft_r2c(nx,isize,istart,osize,ostart,this->m_Opt->GetFFT().mpicomm);
    p_yhat=(FFTScaType*)accfft_alloc(alloc_max);

    // compute fft of input data
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x,this->m_xhat,ffttimers);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

#pragma omp parallel
{
    IntType li=0,lj=0;
    ScalarType k1,k2,k3;
#pragma omp for
    for (IntType i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1){ // x1
        for (IntType i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                k1 = static_cast<ScalarType>(i1 + this->m_Opt->GetFFT().ostart[0]);
                k2 = static_cast<ScalarType>(i2 + this->m_Opt->GetFFT().ostart[1]);
                k3 = static_cast<ScalarType>(i3 + this->m_Opt->GetFFT().ostart[2]);

                // compute linear / flat index
                li = GetLinearIndex(i1,i2,i3,this->m_Opt->GetFFT().osize);

//                p_yhat[lj][0] = this->m_xhat[li][0];
//                p_yhat[lj][1] = this->m_xhat[li][1];

            } // i1
        } // i2
    } // i3

} // pragma omp parallel

    // allocate fft
    u = (ScalarType*)accfft_alloc(alloc_max);
    uk = (Complex*)accfft_alloc(alloc_max);
    fftplan=accfft_plan_dft_3d_r2c(nx,u,(double*)uk,this->m_Opt->GetFFT().mpicomm,ACCFFT_MEASURE);

    // compute fft of input data
    ierr=VecSet(y,0.0); CHKERRQ(ierr);
    ierr=VecGetArray(y,&p_y); CHKERRQ(ierr);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(fftplan,p_yhat,p_y,ffttimers);
    ierr=VecRestoreArray(y,&p_y); CHKERRQ(ierr);

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    // clean up
    if(u != NULL){ accfft_free(u); u=NULL; }
    if(uk != NULL){ accfft_free(uk); uk=NULL; }
    if(p_yhat != NULL){ accfft_free(p_yhat); p_yhat=NULL; }
    if(fftplan != NULL){ accfft_destroy_plan(fftplan); fftplan = NULL; }

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
    size_t alloc_max;
    ScalarType hx[3],nx[3],c[3],scale;
    ScalarType *p_x=NULL, *p_y=NULL;
    double ffttimers[5]={0,0,0,0,0};
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    n[0] = static_cast<int>(this->m_Opt->GetNumGridPoints(0));
    n[1] = static_cast<int>(this->m_Opt->GetNumGridPoints(1));
    n[2] = static_cast<int>(this->m_Opt->GetNumGridPoints(2));

    // get local pencil size and allocation size
    alloc_max=accfft_local_size_dft_r2c_t<ScalarType>(n,isize,istart,osize,ostart,
                                                        this->m_Opt->GetFFT().mpicomm);
    if(this->m_xhat == NULL){
        this->m_xhat=(FFTScaType*)accfft_alloc(alloc_max);
    }

    for (int i = 0; i < 3; ++i){

        hx[i] = this->m_Opt->GetSpatialStepSize(i);
        nx[i] = static_cast<ScalarType>(n[i]);

        // sigma is provided by user in # of grid points
        c[i] = this->m_Opt->GetSigma(i)*hx[i];
        c[i] *= c[i];

        iosize[i] = osize[i];
    }

    scale = this->m_Opt->ComputeFFTScale();

    // compute fft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x,this->m_xhat,ffttimers);
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
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_xhat,p_y,ffttimers);
    ierr=VecRestoreArray(y,&p_y); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);

}




} // end of namespace

#endif // _PREPROCREG_CPP_
