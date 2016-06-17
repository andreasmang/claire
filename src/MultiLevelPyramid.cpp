
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
    this->m_Domain.nx.resize(this->m_NumLevels); // grid size per level
    this->m_Domain.isize.resize(this->m_NumLevels); // grid size per level (spatial domain)
    this->m_Domain.istart.resize(this->m_NumLevels); // start index per level (spatial domain)
    this->m_Domain.osize.resize(this->m_NumLevels); // grid size per level (frequency domain)
    this->m_Domain.ostart.resize(this->m_NumLevels); // start index per level (frequency domain)

   for (int i = 0; i < this->m_NumLevels; ++i){
        this->m_Domain.nx[i].resize(3);

        this->m_Domain.istart[i].resize(3);
        this->m_Domain.isize[i].resize(3);

        this->m_Domain.ostart[i].resize(3);
        this->m_Domain.osize[i].resize(3);
    }

    this->m_Domain.nlocal.resize(this->m_NumLevels); // local points (MPI task) per level
    this->m_Domain.nglobal.resize(this->m_NumLevels); // global points per level
    this->m_Domain.nallocfd.resize(this->m_NumLevels); // alloc size in fourier domain

    level=0;
    while (level < this->m_NumLevels){

        j = this->m_NumLevels - (level+1);

        nl = 1; // reset local size
        ng = 1; // reset global size

        // compute number of grid points for current level
        for (int i = 0; i < 3; ++i){

            if (level==0){
                this->m_Domain.nx[j][i] = this->m_Opt->GetDomainPara().nx[i];
            }
            else{
                n = static_cast<ScalarType>(this->m_Domain.nx[j+1][i]);
                this->m_Domain.nx[j][i] = static_cast<IntType>( std::ceil(n/2.0) );
            }

            // compute global size
            ng *= this->m_Domain.nx[j][i];
            nx[i] = static_cast<int>(this->m_Domain.nx[j][i]);

        }
        this->m_Domain.nglobal[j] = ng;

        // get the local sizes
        this->m_Domain.nallocfd[j]=accfft_local_size_dft_r2c(nx,isize,istart,osize,ostart,this->m_Opt->GetFFT().mpicomm);

        // compute local sizes
        for (int i = 0; i < 3; ++i){

            nl *= static_cast<IntType>(isize[i]);

            this->m_Domain.isize[j][i] = static_cast<IntType>(isize[i]);
            this->m_Domain.istart[j][i] = static_cast<IntType>(istart[i]);

            this->m_Domain.osize[j][i] = static_cast<IntType>(osize[i]);
            this->m_Domain.ostart[j][i] = static_cast<IntType>(ostart[i]);

        }
        this->m_Domain.nlocal[j] = nl;

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
               <<" nx=("<< this->m_Domain.nx[level][0]
               <<","    << this->m_Domain.nx[level][1]
               <<","    << this->m_Domain.nx[level][2]
               << "); (nl,ng)=("<< this->m_Domain.nlocal[level]
               << "," << this->m_Domain.nglobal[level] << ")";

            ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();
        }

        nl = this->m_Domain.nlocal[level];
        ng = this->m_Domain.nglobal[level];

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
    int level,inxlevel[3];
    Vec *xlevel;
    ScalarType *p_x=NULL,*p_xl=NULL,*p_xdummy=NULL,scale;
    typedef ScalarType FFTScalarType[2];
    FFTScalarType *p_xhat=NULL,*p_xlhat=NULL,*p_xhatdummy=NULL;
    double ffttimers[5]={0,0,0,0,0};
    IntType osizel[3];
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

    // get grid sizes/fft scales
    scale = this->m_Opt->ComputeFFTScale();

    level=0;
    while (level < this->m_NumLevels-1){

        for (int i=0; i<3; ++i){
            inxlevel[i] = static_cast<int>(this->m_Domain.nx[level][i]);
            osizel[i] = static_cast<IntType>(this->m_Domain.osize[level][i]);
        }

        // allocate array for restricted data (in spectral domain)
        p_xlhat=(FFTScalarType*)accfft_alloc(this->m_Domain.nallocfd[level]);

        p_xdummy=(ScalarType*)accfft_alloc(this->m_Domain.nallocfd[level]);
        p_xhatdummy=(FFTScalarType*)accfft_alloc(this->m_Domain.nallocfd[level]);

#pragma omp parallel
{
        IntType il,i,k1l,k2l,k3l,i1,i2,i3;
#pragma omp for

        for (IntType i1l = 0; i1l < this->m_Domain.osize[level][0]; ++i1l){ // x1
            for (IntType i2l = 0; i2l < this->m_Domain.osize[level][1]; ++i2l){ // x2
                for (IntType i3l = 0; i3l < this->m_Domain.osize[level][2]; ++i3l){ // x3

                    // compute grid index
                    k1l = i1l + this->m_Domain.ostart[level][0];
                    k2l = i2l + this->m_Domain.ostart[level][1];
                    k3l = i3l + this->m_Domain.ostart[level][2];

                    i1 = i1l;
                    i2 = i2l;
                    i3 = i3l;

                    if (k1l > this->m_Domain.nx[level][0]/2){
                        i1 = (this->m_Opt->GetDomainPara().nx[0]-1) + k1l - this->m_Domain.nx[level][0];
                    }
                    if (k2l > this->m_Domain.nx[level][1]/2){
                        i2 = (this->m_Opt->GetDomainPara().nx[1]-1) + k2l - this->m_Domain.nx[level][1];
                    }
                    if (k3l > this->m_Domain.nx[level][2]/2){
                        i3 = (this->m_Opt->GetDomainPara().nx[2]-1) + k3l - this->m_Domain.nx[level][2];
                    }

                    il = GetLinearIndex(i1l,i2l,i3l,osizel);
                    i  = GetLinearIndex(i1,i2,i3,this->m_Opt->GetFFT().osize);

                    p_xlhat[il][0] = scale*p_xhat[i][0];
                    p_xlhat[il][1] = scale*p_xhat[i][1];

                } // i1l
            } // i2l
        } // i3l

} // pragma omp parallel

        // get pointer to level
        xlevel=NULL;
        ierr=this->GetDataPointer(&xlevel,level); CHKERRQ(ierr);
        ierr=Assert(*xlevel!=NULL, "pointer is null pointer"); CHKERRQ(ierr);

        // compute inverse fft of restricted data

        ierr=VecGetArray(*xlevel,&p_xl); CHKERRQ(ierr);
        //plan=accfft_plan_dft_3d_r2c(inxlevel,p_xl,(double*)p_xlhat,this->m_Opt->GetFFT().mpicomm,ACCFFT_MEASURE);
        plan=accfft_plan_dft_3d_r2c(inxlevel,p_xdummy,(double*)p_xhatdummy,this->m_Opt->GetFFT().mpicomm,ACCFFT_MEASURE);
        accfft_execute_c2r_t<FFTScalarType,ScalarType>(plan,p_xlhat,p_xl,ffttimers);
        ierr=VecRestoreArray(*xlevel,&p_xl); CHKERRQ(ierr);

        // delete data
        if (p_xlhat!=NULL) { accfft_free(p_xlhat); p_xlhat=NULL; }
        if (p_xdummy!=NULL) { accfft_free(p_xdummy); p_xdummy=NULL; }
        if (p_xhatdummy!=NULL) { accfft_free(p_xhatdummy); p_xhatdummy=NULL; }
        if (plan!=NULL){ accfft_destroy_plan(plan); plan=NULL; }

        ++level;
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


/********************************************************************
 * @brief get data at specific level
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetLevel"
PetscErrorCode MultiLevelPyramid::GetDataPointer(Vec** x, int level)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (level == 0){
        *x = &this->m_DataL01;
    }
    else if (level == 1){
        *x = &this->m_DataL02;
    }
    else if (level == 2){
        *x = &this->m_DataL03;
    }
    else if (level == 3){
        *x = &this->m_DataL04;
    }
    else if (level == 4){
        *x = &this->m_DataL05;
    }
    else if (level == 5){
        *x = &this->m_DataL06;
    }
    else if (level == 6){
        *x = &this->m_DataL07;
    }
    else if (level == 7){
        *x = &this->m_DataL08;
    }
    else if (level == 8){
        *x = &this->m_DataL09;
    }
    else if (level == 9){
        *x = &this->m_DataL10;
    }
    else if (level == 10){
        *x = &this->m_DataL11;
    }
    else if (level == 11){
        *x = &this->m_DataL12;
    }
    else if (level == 12){
        *x = &this->m_DataL13;
    }
    else if (level == 13){
        *x = &this->m_DataL14;
    }
    else if (level == 14){
        *x = &this->m_DataL15;
    }
    else{ ierr=ThrowError("level not accessible"); CHKERRQ(ierr); }


    PetscFunctionReturn(0);
}


} // end of namespace

#endif // _MULTILEVELPYRAMID_CPP_
