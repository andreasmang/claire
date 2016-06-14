
#ifndef _MULTILEVELPYRAMID_CPP_
#define _MULTILEVELPYRAMID_CPP_



#include "MultiLevelPyramid.hpp"



namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MultiLevelPyramid"
MultiLevelPyramid::MultiLevelPyramid()
{
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MultiLevelPyramid"
MultiLevelPyramid::MultiLevelPyramid(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;

}




/********************************************************************
 * @brief destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~MultiLevelPyramid"
MultiLevelPyramid::~MultiLevelPyramid()
{
    this->ClearMemory();
}


/********************************************************************
 * @brief initialize class
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode MultiLevelPyramid::Initialize()
{

    this->m_DataL01=NULL;
    this->m_DataL02=NULL;
    this->m_DataL03=NULL;
    this->m_DataL04=NULL;
    this->m_DataL05=NULL;
    this->m_DataL06=NULL;
    this->m_DataL07=NULL;
    this->m_DataL08=NULL;
    this->m_DataL09=NULL;
    this->m_DataL10=NULL;
    this->m_DataL11=NULL;
    this->m_DataL12=NULL;
    this->m_DataL13=NULL;
    this->m_DataL14=NULL;
    this->m_DataL15=NULL;

    this->m_NumLevels=0;

}




/********************************************************************
 * @brief clean up class
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode MultiLevelPyramid::ClearMemory()
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    if(this->m_DataL01!=NULL){
        ierr=VecDestroy(&this->m_DataL01); CHKERRQ(ierr);
        this->m_DataL01=NULL;
    }
    if(this->m_DataL02!=NULL){
        ierr=VecDestroy(&this->m_DataL02); CHKERRQ(ierr);
        this->m_DataL02=NULL;
    }
    if(this->m_DataL03!=NULL){
        ierr=VecDestroy(&this->m_DataL03); CHKERRQ(ierr);
        this->m_DataL03=NULL;
    }
    if(this->m_DataL04!=NULL){
        ierr=VecDestroy(&this->m_DataL04); CHKERRQ(ierr);
        this->m_DataL04=NULL;
    }
    if(this->m_DataL05!=NULL){
        ierr=VecDestroy(&this->m_DataL05); CHKERRQ(ierr);
        this->m_DataL05=NULL;
    }
    if(this->m_DataL06!=NULL){
        ierr=VecDestroy(&this->m_DataL06); CHKERRQ(ierr);
        this->m_DataL06=NULL;
    }
    if(this->m_DataL07!=NULL){
        ierr=VecDestroy(&this->m_DataL07); CHKERRQ(ierr);
        this->m_DataL07=NULL;
    }
    if(this->m_DataL08!=NULL){
        ierr=VecDestroy(&this->m_DataL08); CHKERRQ(ierr);
        this->m_DataL08=NULL;
    }
    if(this->m_DataL09!=NULL){
        ierr=VecDestroy(&this->m_DataL09); CHKERRQ(ierr);
        this->m_DataL09=NULL;
    }
    if(this->m_DataL10!=NULL){
        ierr=VecDestroy(&this->m_DataL10); CHKERRQ(ierr);
        this->m_DataL10=NULL;
    }
    if(this->m_DataL11!=NULL){
        ierr=VecDestroy(&this->m_DataL11); CHKERRQ(ierr);
        this->m_DataL11=NULL;
    }
    if(this->m_DataL12!=NULL){
        ierr=VecDestroy(&this->m_DataL12); CHKERRQ(ierr);
        this->m_DataL12=NULL;
    }
    if(this->m_DataL13!=NULL){
        ierr=VecDestroy(&this->m_DataL13); CHKERRQ(ierr);
        this->m_DataL13=NULL;
    }
    if(this->m_DataL14!=NULL){
        ierr=VecDestroy(&this->m_DataL14); CHKERRQ(ierr);
        this->m_DataL14=NULL;
    }
    if(this->m_DataL15!=NULL){
        ierr=VecDestroy(&this->m_DataL15); CHKERRQ(ierr);
        this->m_DataL15=NULL;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief display level message
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DisplayLevelMsg"
PetscErrorCode MultiLevelPyramid::DisplayLevelMsg(int level)
{
    PetscErrorCode ierr;
    int rank;
    std::stringstream ss;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // display message to user
    ss << std::scientific << "level " << std::setw(3) << level + 1
       <<" of " << this->m_NumLevels
       <<"    nx=("<< this->m_nx[level][0]
       <<","    << this->m_nx[level][1]
       <<","    << this->m_nx[level][2]
       << "); (nl,ng)=("<< this->m_nlocal[level]
       << "," << this->m_nglobal[level] << ")";

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ierr=Msg(ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss.str( std::string() ); ss.clear();

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief compute grid size and number of levels
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeGridSize"
PetscErrorCode MultiLevelPyramid::ComputeGridSize()
{
    PetscErrorCode ierr;
    IntType nxi,nxmin,nl,ng,nalloc;
    int level,j,nx[3],isize[3],istart[3],ostart[3],osize[3],maxlevel,minlevel;
    ScalarType n;

    PetscFunctionBegin;

    // compute number of levels
    nxmin = this->m_Opt->GetDomainPara().nx[0];
    for (int i = 1; i < 3; ++i){
        nxi = this->m_Opt->GetDomainPara().nx[i];
        nxmin = nxmin < nxi ? nxmin : nxi;
    }
    maxlevel = static_cast<int>(std::ceil(std::log2(static_cast<ScalarType>(nxmin))));
    minlevel = static_cast<int>(this->m_Opt->GetGridContPara().minlevel);
    this->m_NumLevels = maxlevel-minlevel;

    ierr=Assert(this->m_NumLevels > 0,"error in size"); CHKERRQ(ierr);

    // allocate memory
    this->m_nx.resize(this->m_NumLevels); // grid size per level
    this->m_isize.resize(this->m_NumLevels); // grid size per level (spatial domain)
    this->m_istart.resize(this->m_NumLevels); // start index per level (spatial domain)
    this->m_osize.resize(this->m_NumLevels); // grid size per level (frequency domain)
    this->m_ostart.resize(this->m_NumLevels); // start index per level (frequency domain)

   for (int i = 0; i < this->m_NumLevels; ++i){
        this->m_nx[i].resize(3);

        this->m_istart[i].resize(3);
        this->m_isize[i].resize(3);

        this->m_ostart[i].resize(3);
        this->m_osize[i].resize(3);
    }

    this->m_nlocal.resize(this->m_NumLevels); // local points (MPI task) per level
    this->m_nglobal.resize(this->m_NumLevels); // global points per level
    this->m_nallocfd.resize(this->m_NumLevels); // alloc size in fourier domain

    level=0;
    while (level < this->m_NumLevels){

        j = this->m_NumLevels - (level+1);

        nl = 1; // reset local size
        ng = 1; // reset global size

        // compute number of grid points for current level
        for (int i = 0; i < 3; ++i){
            if (level == 0) this->m_nx[j][i] = this->m_Opt->GetDomainPara().nx[i];
            else{
                n = static_cast<ScalarType>(this->m_nx[j+1][i]);
                this->m_nx[j][i] = static_cast<IntType>( std::ceil(n/2.0) );
            }

            // compute global size
            ng *= this->m_nx[j][i];
            nx[i] = static_cast<int>(this->m_nx[j][i]);
        }
        this->m_nglobal[j] = ng;

        // get the local sizes
        this->m_nallocfd[j]=accfft_local_size_dft_r2c(nx,isize,istart,osize,ostart,this->m_Opt->GetFFT().mpicomm);

        // compute local sizes
        for (int i = 0; i < 3; ++i){

            nl *= static_cast<IntType>(isize[i]);

            this->m_isize[j][i] = static_cast<IntType>(isize[i]);
            this->m_istart[j][i] = static_cast<IntType>(istart[i]);

            this->m_osize[j][i] = static_cast<IntType>(osize[i]);
            this->m_ostart[j][i] = static_cast<IntType>(ostart[i]);

        }
        this->m_nlocal[j] = nl;

        ++level; // increment

    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief allocate entire pyramid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "AllocatePyramid"
PetscErrorCode MultiLevelPyramid::AllocatePyramid()
{
    PetscErrorCode ierr;
    IntType nl,ng;
    std::stringstream ss;
    int level;

    PetscFunctionBegin;

    ierr=this->ComputeGridSize(); CHKERRQ(ierr);

    level=0;
    while (level < this->m_NumLevels){

        if (this->m_Opt->GetVerbosity() > 2){
            ss << std::scientific << "allocating ML data: level " << std::setw(3) << level + 1
               <<" of " << this->m_NumLevels
               <<" nx=("<< this->m_nx[level][0]
               <<","    << this->m_nx[level][1]
               <<","    << this->m_nx[level][2]
               << "); (nl,ng)=("<< this->m_nlocal[level]
               << "," << this->m_nglobal[level] << ")";

            ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();
        }

        nl = this->m_nlocal[level];
        ng = this->m_nglobal[level];

        if (level == 0){
            ierr=this->Allocate(&this->m_DataL01,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL01,0.0); CHKERRQ(ierr);
        }
        else if (level == 1){
            ierr=this->Allocate(&this->m_DataL02,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL02,0.0); CHKERRQ(ierr);
        }
        else if (level == 2){
            ierr=this->Allocate(&this->m_DataL03,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL03,0.0); CHKERRQ(ierr);
        }
        else if (level == 3){
            ierr=this->Allocate(&this->m_DataL04,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL04,0.0); CHKERRQ(ierr);
        }
        else if (level == 4){
            ierr=this->Allocate(&this->m_DataL05,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL05,0.0); CHKERRQ(ierr);
        }
        else if (level == 5){
            ierr=this->Allocate(&this->m_DataL06,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL06,0.0); CHKERRQ(ierr);
        }
        else if (level == 6){
            ierr=this->Allocate(&this->m_DataL07,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL07,0.0); CHKERRQ(ierr);
        }
        else if (level == 7){
            ierr=this->Allocate(&this->m_DataL08,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL08,0.0); CHKERRQ(ierr);
        }
        else if (level == 8){
            ierr=this->Allocate(&this->m_DataL09,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL09,0.0); CHKERRQ(ierr);
        }
        else if (level == 9){
            ierr=this->Allocate(&this->m_DataL10,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL10,0.0); CHKERRQ(ierr);
        }
        else if (level == 10){
            ierr=this->Allocate(&this->m_DataL11,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL11,0.0); CHKERRQ(ierr);
        }
        else if (level == 11){
            ierr=this->Allocate(&this->m_DataL12,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL12,0.0); CHKERRQ(ierr);
        }
        else if (level == 12){
            ierr=this->Allocate(&this->m_DataL13,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL13,0.0); CHKERRQ(ierr);
        }
        else if (level == 13){
            ierr=this->Allocate(&this->m_DataL14,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL14,0.0); CHKERRQ(ierr);
        }
        else if (level == 14){
            ierr=this->Allocate(&this->m_DataL15,nl,ng); CHKERRQ(ierr);
            ierr=VecSet(this->m_DataL15,0.0); CHKERRQ(ierr);
        }
        else{ ierr=ThrowError("level not accessible"); CHKERRQ(ierr); }

        ++level;

    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief allocate entire pyramid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "AllocatePyramid"
PetscErrorCode MultiLevelPyramid::SetData(Vec x, int level)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (level == 0){
        ierr=VecCopy(x,this->m_DataL01); CHKERRQ(ierr);
    }
    else if (level == 1){
        ierr=VecCopy(x,this->m_DataL02); CHKERRQ(ierr);
    }
    else if (level == 2){
        ierr=VecCopy(x,this->m_DataL03); CHKERRQ(ierr);
    }
    else if (level == 3){
        ierr=VecCopy(x,this->m_DataL04); CHKERRQ(ierr);
    }
    else if (level == 4){
        ierr=VecCopy(x,this->m_DataL05); CHKERRQ(ierr);
    }
    else if (level == 5){
        ierr=VecCopy(x,this->m_DataL06); CHKERRQ(ierr);
    }
    else if (level == 6){
        ierr=VecCopy(x,this->m_DataL07); CHKERRQ(ierr);
    }
    else if (level == 7){
        ierr=VecCopy(x,this->m_DataL08); CHKERRQ(ierr);
    }
    else if (level == 8){
        ierr=VecCopy(x,this->m_DataL09); CHKERRQ(ierr);
    }
    else if (level == 9){
        ierr=VecCopy(x,this->m_DataL10); CHKERRQ(ierr);
    }
    else if (level == 10){
        ierr=VecCopy(x,this->m_DataL11); CHKERRQ(ierr);
    }
    else if (level == 11){
        ierr=VecCopy(x,this->m_DataL12); CHKERRQ(ierr);
    }
    else if (level == 12){
        ierr=VecCopy(x,this->m_DataL13); CHKERRQ(ierr);
    }
    else if (level == 13){
        ierr=VecCopy(x,this->m_DataL14); CHKERRQ(ierr);
    }
    else if (level == 14){
        ierr=VecCopy(x,this->m_DataL15); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("level not accessible"); CHKERRQ(ierr); }


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief allocate
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode MultiLevelPyramid::Allocate(Vec* x, IntType nl, IntType ng)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (*x!=NULL){ ierr=VecDestroy(x); CHKERRQ(ierr); *x=NULL; }

    ierr=VecCreate(PETSC_COMM_WORLD,x); CHKERRQ(ierr);
    ierr=VecSetSizes(*x,nl,ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(*x); CHKERRQ(ierr);
    ierr=VecSet(*x,0.0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetUp"
PetscErrorCode MultiLevelPyramid::SetUp(Vec x)
{
    PetscErrorCode ierr;
    accfft_plan *plan=NULL;
    int level,nxresfft[3];
    Vec xres;
    ScalarType *u=NULL,*p_x=NULL,*p_xres=NULL;
    typedef ScalarType FFTScalarType[2];
    FFTScalarType *p_xhat=NULL, *p_xreshat=NULL;
    double ffttimers[5]={0,0,0,0,0};
    IntType nxres[3],nx[3];
    PetscFunctionBegin;

    // allocate the data pyramid
    ierr=this->AllocatePyramid(); CHKERRQ(ierr);

    // set data on finest grid
    ierr=this->SetData(x,this->m_NumLevels-1); CHKERRQ(ierr);

    // allocate data for fourier domain
    p_xhat=(FFTScalarType*)accfft_alloc(this->m_Opt->GetFFT().nalloc);

    // compute fft of input data
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,FFTScalarType>(this->m_Opt->GetFFT().plan,p_x,p_xhat,ffttimers);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

    // get grid size
    nx[0] = this->m_Opt->GetDomainPara().nx[0];
    nx[1] = this->m_Opt->GetDomainPara().nx[1];
    nx[2] = this->m_Opt->GetDomainPara().nx[2];

    level=this->m_NumLevels-1;
    while (level > 0){

        // get number of grid points for current level
        for (int i=0; i<3; ++i){
            nxres[i] = this->m_nx[level][i];
            nxresfft[i] = static_cast<int>(this->m_nx[level][i]);
        }

        // allocate array for restricted data (in spectral domain)
        p_xreshat=(FFTScalarType*)accfft_alloc(this->m_nallocfd[level]);

#pragma omp parallel
{
        IntType li=0,lj=0,k[3];
        bool inbounds;
#pragma omp for

        for (IntType i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1){ // x1
            for (IntType i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2){ // x2
                for (IntType i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3){ // x3

                    // compute coordinates (nodal grid)
                    k[0] = i1 + this->m_Opt->GetFFT().ostart[0];
                    k[1] = i2 + this->m_Opt->GetFFT().ostart[1];
                    k[2] = i3 + this->m_Opt->GetFFT().ostart[2];

                    if (k[0] > nx[0]/2) k[0]-=nx[0];
                    if (k[1] > nx[1]/2) k[1]-=nx[1];
                    if (k[2] > nx[2]/2) k[2]-=nx[2];

                    inbounds=true;
                    inbounds =            (abs(k[0]) < nxres[0]/2);
                    inbounds = inbounds ? (abs(k[1]) < nxres[1]/2) : inbounds;
                    inbounds = inbounds ? (abs(k[2]) < nxres[2]/2) : inbounds;

                    if (inbounds){
                        // compute linear / flat index
                        li = GetLinearIndex(i1,i2,i3,this->m_Opt->GetFFT().osize);
                        ++lj;
                        p_xreshat[lj][0] = p_xhat[li][0];
                        p_xreshat[lj][1] = p_xhat[li][1];
                    }
                } // i1
            } // i2
        } // i3
} // pragma omp parallel

        // get pointer to level
        xres=NULL;
        ierr=this->GetLevel(&xres,level); CHKERRQ(ierr);
        ierr=Assert(xres!=NULL, "null pointer"); CHKERRQ(ierr);

        // compute inverse fft of restricted data
        ierr=VecGetArray(xres,&p_xres); CHKERRQ(ierr);
        plan=accfft_plan_dft_3d_r2c(nxresfft,p_xres,(double*)p_xreshat,this->m_Opt->GetFFT().mpicomm,ACCFFT_MEASURE);
        ierr=VecRestoreArray(xres,&p_xres); CHKERRQ(ierr);

        // delete data
        accfft_free(p_xreshat); p_xreshat=NULL;
        accfft_destroy_plan(plan); plan = NULL;

        --level; // decrement
    }

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(ffttimers);


    // clean up
    if(p_xhat != NULL){ accfft_free(p_xhat); p_xhat=NULL; }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief get data at specific level
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetLevel"
PetscErrorCode MultiLevelPyramid::GetLevel(Vec* x, int level)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (level == 0){
        ierr=VecDuplicate(this->m_DataL01,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL01,*x); CHKERRQ(ierr);
//        x = &this->m_DataL01;
    }
    else if (level == 1){
        ierr=VecDuplicate(this->m_DataL02,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL02,*x); CHKERRQ(ierr);
//        x = &this->m_DataL02;
    }
    else if (level == 2){
        ierr=VecDuplicate(this->m_DataL03,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL03,*x); CHKERRQ(ierr);
//        x = &this->m_DataL03;
    }
    else if (level == 3){
        ierr=VecDuplicate(this->m_DataL04,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL04,*x); CHKERRQ(ierr);
//        x = &this->m_DataL04;
    }
    else if (level == 4){
        ierr=VecDuplicate(this->m_DataL05,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL05,*x); CHKERRQ(ierr);
//        x = &this->m_DataL05;
    }
    else if (level == 5){
        ierr=VecDuplicate(this->m_DataL06,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL06,*x); CHKERRQ(ierr);
//        x = &this->m_DataL06;
    }
    else if (level == 6){
        ierr=VecDuplicate(this->m_DataL07,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL07,*x); CHKERRQ(ierr);
//        x = &this->m_DataL07;
    }
    else if (level == 7){
        ierr=VecDuplicate(this->m_DataL08,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL08,*x); CHKERRQ(ierr);
//        x = &this->m_DataL08;
    }
    else if (level == 8){
        ierr=VecDuplicate(this->m_DataL09,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL09,*x); CHKERRQ(ierr);
//        x = &this->m_DataL09;
    }
    else if (level == 9){
        ierr=VecDuplicate(this->m_DataL10,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL10,*x); CHKERRQ(ierr);
//        x = &this->m_DataL10;
    }
    else if (level == 10){
        ierr=VecDuplicate(this->m_DataL11,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL11,*x); CHKERRQ(ierr);
//        x = &this->m_DataL11;
    }
    else if (level == 11){
        ierr=VecDuplicate(this->m_DataL12,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL12,*x); CHKERRQ(ierr);
//        x = &this->m_DataL12;
    }
    else if (level == 12){
        ierr=VecDuplicate(this->m_DataL13,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL13,*x); CHKERRQ(ierr);
//        x = &this->m_DataL13;
    }
    else if (level == 13){
        ierr=VecDuplicate(this->m_DataL14,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL14,*x); CHKERRQ(ierr);
//        x = &this->m_DataL14;
    }
    else if (level == 14){
        ierr=VecDuplicate(this->m_DataL15,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL15,*x); CHKERRQ(ierr);
//        x = &this->m_DataL15;
    }
    else{ ierr=ThrowError("level not accessible"); CHKERRQ(ierr); }




    PetscFunctionReturn(0);
}




} // end of namespace

#endif // _MULTILEVELPYRAMID_CPP_
