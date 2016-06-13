
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

    this->m_MaxLevel=0;
    this->m_MinLevel=0;
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
    ss << std::scientific << "level " << std::setw(3) << level
       <<" nx=("<< std::setw(4) << this->m_nx[level][0]
       <<","    << std::setw(4) << this->m_nx[level][1]
       <<","    << std::setw(4) << this->m_nx[level][2]
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
    IntType nxi,nxmin,nl,ng;
    int level,j,nx[3],isize[3],istart[3],ostart[3],osize[3];
    ScalarType n;

    PetscFunctionBegin;

    // compute number of levels
    nxmin = this->m_Opt->GetDomainPara().nx[0];
    for (int i = 1; i < 3; ++i){
        nxi = this->m_Opt->GetDomainPara().nx[i];
        nxmin = nxmin < nxi ? nxmin : nxi;
    }
    this->m_MaxLevel  = static_cast<int>(std::ceil(std::log2(static_cast<ScalarType>(nxmin))));
    this->m_MinLevel  = static_cast<int>(this->m_Opt->GetGridContPara().minlevel);
    this->m_NumLevels = this->m_MaxLevel-this->m_MinLevel;
    ierr=Assert(this->m_NumLevels > 0,"error in size"); CHKERRQ(ierr);

    // allocate memory
    this->m_nx.resize(this->m_NumLevels); // grid size per level
    for (int i = 0; i < this->m_NumLevels; ++i) this->m_nx[i].resize(3);
    this->m_nlocal.resize(this->m_NumLevels); // local points (MPI task) per level
    this->m_nglobal.resize(this->m_NumLevels); // global points per level

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
        accfft_local_size_dft_r2c(nx,isize,istart,osize,ostart,this->m_Opt->GetFFT().mpicomm);

        // compute local sizes
        for (int i = 0; i < 3; ++i){
            nl *= static_cast<IntType>(isize[i]);
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
    int level;

    PetscFunctionBegin;

    ierr=this->ComputeGridSize(); CHKERRQ(ierr);

    level=0;
    while (level < this->m_NumLevels){

        // display level message
        ierr=this->DisplayLevelMsg(level); CHKERRQ(ierr);

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


    }

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

    if (x!=NULL){
        ierr=VecDestroy(x); CHKERRQ(ierr);
        x=NULL;
    }

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
PetscErrorCode MultiLevelPyramid::SetUp()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=this->AllocatePyramid(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of namespace

#endif // _MULTILEVELPYRAMID_CPP_
