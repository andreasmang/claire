#ifndef _VECFIELD_CPP_
#define _VECFIELD_CPP_



#include "VecField.hpp"



namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecField"
VecField::VecField()
{
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~VecField"
VecField::~VecField(void)
{
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecField"
VecField::VecField(RegOpt* opt)
{
    this->Initialize();
    this->SetOpt(opt);
    this->Allocate();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecField"
VecField::VecField(RegOpt* opt,int level)
{
    this->Initialize();
    this->SetOpt(opt);
    this->Allocate(level);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "VecField"
VecField::VecField(IntType nl, IntType ng)
{
    this->Initialize();
    this->Allocate(nl,ng);

}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode VecField::Initialize(void)
{
    PetscFunctionBegin;

    this->m_Opt=NULL;

    this->m_X1=NULL;
    this->m_X2=NULL;
    this->m_X3=NULL;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode VecField::ClearMemory(void)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if(this->m_X1 != NULL) { ierr=VecDestroy(&this->m_X1); CHKERRQ(ierr); this->m_X1 = NULL; }
    if(this->m_X2 != NULL) { ierr=VecDestroy(&this->m_X2); CHKERRQ(ierr); this->m_X2 = NULL; }
    if(this->m_X3 != NULL) { ierr=VecDestroy(&this->m_X3); CHKERRQ(ierr); this->m_X3 = NULL; }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set the options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetOpt"
PetscErrorCode VecField::SetOpt(RegOpt* opt)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(opt!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_Opt = opt;


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode VecField::Allocate()
{
    PetscErrorCode ierr;
    IntType nl,ng;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr=this->ClearMemory(); CHKERRQ(ierr);

    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    // allocate vector field
    ierr=VecCreate(PETSC_COMM_WORLD,&this->m_X1); CHKERRQ(ierr);
    ierr=VecSetSizes(this->m_X1,nl,ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(this->m_X1); CHKERRQ(ierr);

    // pass options
    ierr=VecDuplicate(this->m_X1,&this->m_X2); CHKERRQ(ierr);
    ierr=VecDuplicate(this->m_X1,&this->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode VecField::Allocate(int level)
{
    PetscErrorCode ierr;
    IntType nl,ng;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr=this->ClearMemory(); CHKERRQ(ierr);

    nl = this->m_Opt->GetGridContPara().nlocal[level];
    ng = this->m_Opt->GetGridContPara().nglobal[level];

    // allocate vector field
    ierr=VecCreate(PETSC_COMM_WORLD,&this->m_X1); CHKERRQ(ierr);
    ierr=VecSetSizes(this->m_X1,nl,ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(this->m_X1); CHKERRQ(ierr);

    // pass options
    ierr=VecDuplicate(this->m_X1,&this->m_X2); CHKERRQ(ierr);
    ierr=VecDuplicate(this->m_X1,&this->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief function to allocate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Allocate"
PetscErrorCode VecField::Allocate(IntType nl, IntType ng)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // make sure, that all pointers are deallocated
    ierr=this->ClearMemory(); CHKERRQ(ierr);

    // allocate vector field
    ierr=VecCreate(PETSC_COMM_WORLD,&this->m_X1); CHKERRQ(ierr);
    ierr=VecSetSizes(this->m_X1,nl,ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(this->m_X1); CHKERRQ(ierr);

    // pass options
    ierr=VecDuplicate(this->m_X1,&this->m_X2); CHKERRQ(ierr);
    ierr=VecDuplicate(this->m_X1,&this->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}






/********************************************************************
 * @brief Copy
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Copy"
PetscErrorCode VecField::Copy(VecField* v)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=VecCopy(v->m_X1,this->m_X1); CHKERRQ(ierr);
    ierr=VecCopy(v->m_X2,this->m_X2); CHKERRQ(ierr);
    ierr=VecCopy(v->m_X3,this->m_X3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set value
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetValue"
PetscErrorCode VecField::SetValue(ScalarType value)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=VecSet(this->m_X1,value); CHKERRQ(ierr);
    ierr=VecSet(this->m_X2,value); CHKERRQ(ierr);
    ierr=VecSet(this->m_X3,value); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief sets the individual components of a vector field;
 * the input is a flat petsc array
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetComponents"
PetscErrorCode VecField::SetComponents(Vec w)
{
    PetscErrorCode ierr;
    IntType nl,n;
    ScalarType *p_x1,*p_x2,*p_x3;
    const ScalarType *p_w;

    PetscFunctionBegin;

    // get local size of vector field
    ierr=VecGetLocalSize(w,&n); CHKERRQ(ierr);

    ierr=Assert(this->m_X1!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_X2!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_X3!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=VecGetArrayRead(w,&p_w); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X3,&p_x3); CHKERRQ(ierr);

    //compute size of each individual component
    nl = n / 3;
    //ierr=Assert(nl==this->m_Opt->GetDomainPara().nlocal,"dimension mismatch"); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i){
        p_x1[i] = p_w[i     ];
        p_x2[i] = p_w[i+  nl];
        p_x3[i] = p_w[i+2*nl];
    }
} // pragma omp parallel


    ierr=VecRestoreArrayRead(w,&p_w); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X3,&p_x3); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief get components of vector field and store them
 * in a flat vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetComponents"
PetscErrorCode VecField::GetComponents(Vec w)
{
    PetscErrorCode ierr;
    IntType nl,n;
    ScalarType *p_x1,*p_x2,*p_x3,*p_w;

    PetscFunctionBegin;

    // get local size of vector field
    ierr=VecGetLocalSize(w,&n); CHKERRQ(ierr);

    ierr=VecGetArray(w,&p_w); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X3,&p_x3); CHKERRQ(ierr);

    //compute size of each individual component
    nl = n / 3;

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i){
        p_w[i     ] = p_x1[i];
        p_w[i+  nl] = p_x2[i];
        p_w[i+2*nl] = p_x3[i];

    }
} // pragma omp parallel

    ierr=VecRestoreArray(w,&p_w); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X3,&p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief scale vector by scalar value
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Scale"
PetscErrorCode VecField::Scale(ScalarType value)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=VecScale(this->m_X1,value); CHKERRQ(ierr);
    ierr=VecScale(this->m_X2,value); CHKERRQ(ierr);
    ierr=VecScale(this->m_X3,value); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief pointwise scale of vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Scale"
PetscErrorCode VecField::Scale(Vec s)
{
    PetscErrorCode ierr;
    IntType nl;
    ScalarType *p_vx1,*p_vx2,*p_vx3,*p_s;

    PetscFunctionBegin;

    // get pointers
    ierr=VecGetArray(s,&p_s); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X3,&p_vx3); CHKERRQ(ierr);

    // get local size of vector field
    ierr=VecGetLocalSize(s,&nl); CHKERRQ(ierr);
    //nl = this->m_Opt->GetDomainPara().nlocal;

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i){
        ScalarType scale = p_s[i];
        p_vx1[i] = scale*p_vx1[i];
        p_vx2[i] = scale*p_vx2[i];
        p_vx3[i] = scale*p_vx3[i];
    }
} // pragma omp parallel

    // get pointers
    ierr=VecRestoreArray(s,&p_s); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X3,&p_vx3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief pointwise scale of a vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Scale"
PetscErrorCode VecField::Scale(VecField* v,Vec s)
{
    PetscErrorCode ierr;
    IntType nl;
    ScalarType *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_s=NULL,*p_svx1=NULL,*p_svx2=NULL,
                *p_svx3=NULL;

    PetscFunctionBegin;

    // get pointers
    ierr=VecGetArray(s,&p_s); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecGetArray(v->m_X1,&p_svx1); CHKERRQ(ierr);
    ierr=VecGetArray(v->m_X2,&p_svx2); CHKERRQ(ierr);
    ierr=VecGetArray(v->m_X3,&p_svx3); CHKERRQ(ierr);

    // get local size of vector field
    ierr=VecGetLocalSize(s,&nl); CHKERRQ(ierr);
    //nl = this->m_Opt->GetDomainPara().nlocal;

#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i){
        ScalarType scale = p_s[i];
        p_svx1[i] = scale*p_vx1[i];
        p_svx2[i] = scale*p_vx2[i];
        p_svx3[i] = scale*p_vx3[i];
    }
} // pragma omp parallel

    // get pointers
    ierr=VecRestoreArray(s,&p_s); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(v->m_X1,&p_svx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(v->m_X2,&p_svx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(v->m_X3,&p_svx3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



} // end of name space



#endif // _VECFIELD_CPP_
