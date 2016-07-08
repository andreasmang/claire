
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

    this->m_PreProc=NULL;

    PetscFunctionReturn(0);
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
 * @brief set read write operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetPreProc"
PetscErrorCode MultiLevelPyramid::SetPreProc(PreProcReg* ppr)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(ppr != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_PreProc = ppr;

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
    int level,nlevels;

    PetscFunctionBegin;

    // set up parametes for grid continuation
    ierr=this->m_Opt->SetupGridCont(); CHKERRQ(ierr);

    // get number of levels
    nlevels = this->m_Opt->GetGridContPara().nlevels;

    level=0;
    while (level < nlevels){

        nl = this->m_Opt->GetGridContPara().nlocal[level];
        ng = this->m_Opt->GetGridContPara().nglobal[level];

        if (this->m_Opt->GetVerbosity() > 2){
            ss << std::scientific << "allocating ML data: level " << std::setw(3) << level + 1
               <<" of " << nlevels
               <<" nx=("<< this->m_Opt->GetGridContPara().nx[level][0]
               <<","    << this->m_Opt->GetGridContPara().nx[level][1]
               <<","    << this->m_Opt->GetGridContPara().nx[level][2]
               << "); (nl,ng)=("<< nl << "," << ng << ")";

            ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();
        }


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
#define __FUNCT__ "SetData"
PetscErrorCode MultiLevelPyramid::SetData(Vec x, int level)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if      (level ==  0){ ierr=VecCopy(x,this->m_DataL01); CHKERRQ(ierr); }
    else if (level ==  1){ ierr=VecCopy(x,this->m_DataL02); CHKERRQ(ierr); }
    else if (level ==  2){ ierr=VecCopy(x,this->m_DataL03); CHKERRQ(ierr); }
    else if (level ==  3){ ierr=VecCopy(x,this->m_DataL04); CHKERRQ(ierr); }
    else if (level ==  4){ ierr=VecCopy(x,this->m_DataL05); CHKERRQ(ierr); }
    else if (level ==  5){ ierr=VecCopy(x,this->m_DataL06); CHKERRQ(ierr); }
    else if (level ==  6){ ierr=VecCopy(x,this->m_DataL07); CHKERRQ(ierr); }
    else if (level ==  7){ ierr=VecCopy(x,this->m_DataL08); CHKERRQ(ierr); }
    else if (level ==  8){ ierr=VecCopy(x,this->m_DataL09); CHKERRQ(ierr); }
    else if (level ==  9){ ierr=VecCopy(x,this->m_DataL10); CHKERRQ(ierr); }
    else if (level == 10){ ierr=VecCopy(x,this->m_DataL11); CHKERRQ(ierr); }
    else if (level == 11){ ierr=VecCopy(x,this->m_DataL12); CHKERRQ(ierr); }
    else if (level == 12){ ierr=VecCopy(x,this->m_DataL13); CHKERRQ(ierr); }
    else if (level == 13){ ierr=VecCopy(x,this->m_DataL14); CHKERRQ(ierr); }
    else if (level == 14){ ierr=VecCopy(x,this->m_DataL15); CHKERRQ(ierr); }
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
#define __FUNCT__ "DoSetup"
PetscErrorCode MultiLevelPyramid::DoSetup(Vec x)
{
    PetscErrorCode ierr;
    IntType nxlevel[3],nx[3];
    int nlevels;
    Vec *xlevel;
    PetscFunctionBegin;

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_PreProc!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_PreProc->ResetGridChangeOperators(true);

    // allocate the data pyramid
    ierr=this->AllocatePyramid(); CHKERRQ(ierr);

    // get number of levels
    nlevels=this->m_Opt->GetGridContPara().nlevels;

    // set data on finest grid
    ierr=this->SetData(x,nlevels-1); CHKERRQ(ierr);

    nx[0] = this->m_Opt->GetDomainPara().nx[0];
    nx[1] = this->m_Opt->GetDomainPara().nx[1];
    nx[2] = this->m_Opt->GetDomainPara().nx[2];

    // for all levels
    for (int l = 0; l < nlevels-1; ++l){

        // get grid size
        for (int i=0; i<3; ++i){
            nxlevel[i] = this->m_Opt->GetGridContPara().nx[l][i];
        }

        // get pointer to level
        ierr=this->GetData(&xlevel,l); CHKERRQ(ierr);
        ierr=Assert(*xlevel!=NULL,"null pointer"); CHKERRQ(ierr);

        // restrict data
        ierr=this->m_PreProc->Restrict(xlevel,x,nxlevel,nx); CHKERRQ(ierr);

    }

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
    }
    else if (level == 1){
        ierr=VecDuplicate(this->m_DataL02,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL02,*x); CHKERRQ(ierr);
    }
    else if (level == 2){
        ierr=VecDuplicate(this->m_DataL03,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL03,*x); CHKERRQ(ierr);
    }
    else if (level == 3){
        ierr=VecDuplicate(this->m_DataL04,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL04,*x); CHKERRQ(ierr);
    }
    else if (level == 4){
        ierr=VecDuplicate(this->m_DataL05,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL05,*x); CHKERRQ(ierr);
    }
    else if (level == 5){
        ierr=VecDuplicate(this->m_DataL06,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL06,*x); CHKERRQ(ierr);
    }
    else if (level == 6){
        ierr=VecDuplicate(this->m_DataL07,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL07,*x); CHKERRQ(ierr);
    }
    else if (level == 7){
        ierr=VecDuplicate(this->m_DataL08,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL08,*x); CHKERRQ(ierr);
    }
    else if (level == 8){
        ierr=VecDuplicate(this->m_DataL09,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL09,*x); CHKERRQ(ierr);
    }
    else if (level == 9){
        ierr=VecDuplicate(this->m_DataL10,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL10,*x); CHKERRQ(ierr);
    }
    else if (level == 10){
        ierr=VecDuplicate(this->m_DataL11,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL11,*x); CHKERRQ(ierr);
    }
    else if (level == 11){
        ierr=VecDuplicate(this->m_DataL12,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL12,*x); CHKERRQ(ierr);
    }
    else if (level == 12){
        ierr=VecDuplicate(this->m_DataL13,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL13,*x); CHKERRQ(ierr);
    }
    else if (level == 13){
        ierr=VecDuplicate(this->m_DataL14,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL14,*x); CHKERRQ(ierr);
    }
    else if (level == 14){
        ierr=VecDuplicate(this->m_DataL15,x); CHKERRQ(ierr);
        ierr=VecCopy(this->m_DataL15,*x); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("level not accessible"); CHKERRQ(ierr); }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief get data at specific level
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetData"
PetscErrorCode MultiLevelPyramid::GetData(Vec** x, int level)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if      (level ==  0) *x = &this->m_DataL01;
    else if (level ==  1) *x = &this->m_DataL02;
    else if (level ==  2) *x = &this->m_DataL03;
    else if (level ==  3) *x = &this->m_DataL04;
    else if (level ==  4) *x = &this->m_DataL05;
    else if (level ==  5) *x = &this->m_DataL06;
    else if (level ==  6) *x = &this->m_DataL07;
    else if (level ==  7) *x = &this->m_DataL08;
    else if (level ==  8) *x = &this->m_DataL09;
    else if (level ==  9) *x = &this->m_DataL10;
    else if (level == 10) *x = &this->m_DataL11;
    else if (level == 11) *x = &this->m_DataL12;
    else if (level == 12) *x = &this->m_DataL13;
    else if (level == 13) *x = &this->m_DataL14;
    else if (level == 14) *x = &this->m_DataL15;
    else{ ierr=ThrowError("level not accessible"); CHKERRQ(ierr); }

    PetscFunctionReturn(0);
}




} // end of namespace

#endif // _MULTILEVELPYRAMID_CPP_
