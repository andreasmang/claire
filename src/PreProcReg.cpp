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
