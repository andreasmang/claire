#ifndef _SEMILAGRANGIAN_CPP_
#define _SEMILAGRANGIAN_CPP_

#include "SemiLagrangian.hpp"


namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SemiLagrangian"
SemiLagrangian::SemiLagrangian()
{
    this->Initialize();
}




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SemiLagrangian"
SemiLagrangian::SemiLagrangian(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~SemiLagrangian"
SemiLagrangian::~SemiLagrangian()
{
    this->ClearMemory();
}




/********************************************************************
 * @brief init class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode SemiLagrangian::Initialize()
{
    PetscFunctionBegin;

    this->m_ReadWrite=NULL;

    this->m_WorkVecField = NULL;

    this->m_XA=NULL;
    this->m_XS=NULL;
    this->m_X=NULL;

    this->m_AdjointPlan = NULL;
    this->m_AdjointPlanVec = NULL;

    this->m_StatePlan = NULL;
    this->m_StatePlanVec = NULL;

    this->m_VecFieldPlan = NULL;

    this->m_ScaFieldGhost=NULL;
    this->m_VecFieldGhost=NULL;

    this->m_Opt=NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clears memory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode SemiLagrangian::ClearMemory()
{

    PetscFunctionBegin;

    if(this->m_WorkVecField!=NULL){
        delete this->m_WorkVecField;
        this->m_WorkVecField = NULL;
    }

    if(this->m_X!=NULL){
        delete this->m_X; this->m_X = NULL;
    }

    if(this->m_XS!=NULL){
        delete this->m_XS; this->m_XS = NULL;
    }
    if(this->m_XA!=NULL){
        delete this->m_XA; this->m_XA = NULL;
    }

    if(this->m_AdjointPlan!=NULL){
        delete this->m_AdjointPlan;
        this->m_AdjointPlan = NULL;
    }
    if(this->m_StatePlan!=NULL){
        delete this->m_StatePlan;
        this->m_StatePlan = NULL;
    }
    if(this->m_StatePlanVec!=NULL){
        delete this->m_StatePlanVec;
        this->m_StatePlanVec = NULL;
    }
    if(this->m_AdjointPlanVec!=NULL){
        delete this->m_AdjointPlanVec;
        this->m_AdjointPlanVec = NULL;
    }
    if(this->m_VecFieldPlan!=NULL){
        delete this->m_VecFieldPlan;
        this->m_VecFieldPlan = NULL;
    }

    if (this->m_ScaFieldGhost!=NULL){
        accfft_free(this->m_ScaFieldGhost);
        this->m_ScaFieldGhost=NULL;
    }

    if (this->m_VecFieldGhost!=NULL){
        accfft_free(this->m_VecFieldGhost);
        this->m_VecFieldGhost=NULL;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReadWrite"
PetscErrorCode SemiLagrangian::SetReadWrite(ReadWriteReg* rw)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(rw!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite = rw;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief compute the trajectory from the velocity field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeTrajectory"
PetscErrorCode SemiLagrangian::ComputeTrajectory(VecField* v, std::string flag)
{

    PetscErrorCode ierr=0;
    ScalarType ht,hthalf,hx[3],x1,x2,x3;
    ScalarType *p_v1=NULL,*p_v2=NULL,*p_v3=NULL,
                *p_vX1=NULL,*p_vX2=NULL,*p_vX3=NULL;
    IntType isize[3],istart[3],l,i1,i2,i3,nl;
    ScalarType* X=NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    if (this->m_WorkVecField==NULL){
        try{this->m_WorkVecField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    nl = this->m_Opt->GetDomainPara().nlocal;

    if (strcmp(flag.c_str(),"state")==0){

        if (this->m_XS == NULL){
            try{ this->m_XS = new double [3*nl]; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        X = this->m_XS;

    }
    else if (strcmp(flag.c_str(),"adjoint")==0){

        if (this->m_XA == NULL){
            try{ this->m_XA = new double [3*nl]; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        X = this->m_XA;

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    for (int i = 0; i < 3; ++i){
        hx[i] = this->m_Opt->GetDomainPara().hx[i];
        isize[i] = this->m_Opt->GetDomainPara().isize[i];
        istart[i] = this->m_Opt->GetDomainPara().istart[i];
    }


    // \tilde{X} = x - ht v
    ierr=v->GetArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
    for (i1 = 0; i1 < isize[0]; ++i1){  // x1
        for (i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (i3 = 0; i3 < isize[2]; ++i3){ // x3

                // compute linear / flat index
                l = GetLinearIndex(i1,i2,i3,isize);

                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                X[l*3+0] = x1 - ht*p_v1[l];
                X[l*3+1] = x2 - ht*p_v2[l];
                X[l*3+2] = x3 - ht*p_v3[l];

            } // i1
        } // i2
    } // i3
}// pragma omp for
    ierr=v->RestoreArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);

    // normalize to [0,1]
    for(IntType i = 0; i < 3*nl; ++i) X[i] /= (2.0*PETSC_PI);

    // interpolate velocity field v(X)
    ierr=this->MapCoordinateVector(flag); CHKERRQ(ierr);
    ierr=this->Interpolate(this->m_WorkVecField,v,flag); CHKERRQ(ierr);

    // normalize to [0,2*pi]
    for(IntType i = 0; i < 3*nl; ++i) X[i] *= (2.0*PETSC_PI);


    // X = x - 0.5*ht*(v + v(x - ht v))
    ierr=v->GetArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);
    ierr=this->m_WorkVecField->GetArrays(p_vX1,p_vX2,p_vX3); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
    for (i1 = 0; i1 < isize[0]; ++i1){  // x1
        for (i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (i3 = 0; i3 < isize[2]; ++i3){ // x3

                // compute linear / flat index
                l = GetLinearIndex(i1,i2,i3,isize);

                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                X[l*3+0] = x1 - hthalf*(p_vX1[l]+p_v1[l]);
                X[l*3+1] = x2 - hthalf*(p_vX2[l]+p_v2[l]);
                X[l*3+2] = x3 - hthalf*(p_vX3[l]+p_v3[l]);

            } // i1
        } // i2
    } // i3
}// pragma omp for
    ierr=this->m_WorkVecField->RestoreArrays(p_vX1,p_vX2,p_vX3); CHKERRQ(ierr);
    ierr=v->RestoreArrays(p_v1,p_v2,p_v3); CHKERRQ(ierr);

    // normalize to [0,1]
    for(IntType i = 0; i < 3*nl; ++i) X[i] /= (2.0*PETSC_PI);

    ierr=this->MapCoordinateVector(flag); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief interpolate scalar field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangian::Interpolate(Vec* w,Vec v,std::string flag)
{
    PetscErrorCode ierr;
    ScalarType *p_w=NULL,*p_v=NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(*w!=NULL,"output is null pointer"); CHKERRQ(ierr);
    ierr=Assert(v!=NULL,"input is null pointer"); CHKERRQ(ierr);

    ierr=VecGetArray(v,&p_v); CHKERRQ(ierr);
    ierr=VecGetArray(*w,&p_w); CHKERRQ(ierr);

    ierr=this->Interpolate(p_w,p_v,flag); CHKERRQ(ierr);

    ierr=VecRestoreArray(*w,&p_w); CHKERRQ(ierr);
    ierr=VecRestoreArray(v,&p_v); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief interpolate scalar field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangian::Interpolate(ScalarType* w,ScalarType* v,std::string flag)
{
    PetscErrorCode ierr;
    int nx[3],isize_g[3],isize[3],istart_g[3],istart[3],c_dims[2],nl;
    accfft_plan* plan=NULL;
    IntType g_alloc_max;
    double timers[4]={0,0,0,0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(w!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(v!=NULL,"null pointer"); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<int>(this->m_Opt->GetDomainPara().nx[i]);
        isize[i] = static_cast<int>(this->m_Opt->GetDomainPara().isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->GetDomainPara().istart[i]);
    }

    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    nl =static_cast<int>( this->m_Opt->GetDomainPara().nlocal );
    ierr=Assert(nl!=0,"size problem"); CHKERRQ(ierr);

    // deal with ghost points
    plan = this->m_Opt->GetFFT().plan;
    g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(plan,this->m_GhostSize,isize_g,istart_g);

    if(this->m_ScaFieldGhost==NULL){
        this->m_ScaFieldGhost = (ScalarType*)accfft_alloc(g_alloc_max);
    }

    accfft_get_ghost_xyz(plan,this->m_GhostSize,isize_g,v,this->m_ScaFieldGhost);

    if (strcmp(flag.c_str(),"state")==0){

        ierr=Assert(this->m_XS!=NULL,"state X is null pointer"); CHKERRQ(ierr);

        this->m_StatePlan->interpolate(this->m_ScaFieldGhost,1,nx,isize,istart,
                                        nl,this->m_GhostSize,w,c_dims,
                                        this->m_Opt->GetFFT().mpicomm,timers);
    }
    else if (strcmp(flag.c_str(),"adjoint")==0){

        ierr=Assert(this->m_XA!=NULL,"adjoint X is null pointer"); CHKERRQ(ierr);

        this->m_AdjointPlan->interpolate(this->m_ScaFieldGhost,1,nx,isize,istart,
                                        nl,this->m_GhostSize,w,c_dims,
                                        this->m_Opt->GetFFT().mpicomm,timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IP);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief interpolate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangian::Interpolate(VecField* w, VecField* v, std::string flag)
{
    PetscErrorCode ierr;
    ScalarType *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_wx1=NULL,*p_wx2=NULL,*p_wx3=NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(v!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(w!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=w->GetArrays(p_wx1,p_wx2,p_wx3); CHKERRQ(ierr);
    ierr=v->GetArrays(p_vx1,p_vx2,p_vx3); CHKERRQ(ierr);

    ierr=this->Interpolate(p_wx1,p_wx2,p_wx3,p_vx1,p_vx2,p_vx3,flag); CHKERRQ(ierr);

    ierr=w->RestoreArrays(p_wx1,p_wx2,p_wx3); CHKERRQ(ierr);
    ierr=v->RestoreArrays(p_vx1,p_vx2,p_vx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief interpolate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangian::Interpolate( ScalarType* wx1,
                                            ScalarType* wx2,
                                            ScalarType* wx3,
                                            ScalarType* vx1,
                                            ScalarType* vx2,
                                            ScalarType* vx3,
                                            std::string flag)
{
    PetscErrorCode ierr;
    int nx[3],isize_g[3],isize[3],istart_g[3],istart[3],c_dims[2];
    double timers[4] = {0,0,0,0};
    accfft_plan* plan=NULL;
    IntType g_alloc_max;
    IntType nl,nlghost;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(vx1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vx2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vx3!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx3!=NULL,"input is null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->GetDomainPara().nlocal;

    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<int>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<int>(this->m_Opt->GetDomainPara().isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->GetDomainPara().istart[i]);
    }

    // get network dimensions
    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    if (this->m_X==NULL){
        try{ this->m_X = new double [3*nl]; }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    for (IntType i = 0; i < nl; ++i){
        this->m_X[0*nl+i] = vx1[i];
        this->m_X[1*nl+i] = vx2[i];
        this->m_X[2*nl+i] = vx3[i];
    }

    // get ghost sizes
    plan = this->m_Opt->GetFFT().plan;
    g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(plan,this->m_GhostSize,isize_g,istart_g);
    ierr=Assert(g_alloc_max!=0,"alloc problem"); CHKERRQ(ierr);

    // get nlocal for ghosts
    nlghost = 1;
    for (IntType i = 0; i < 3; ++i){
        nlghost *= static_cast<IntType>(isize_g[i]);
    }

    // deal with ghost points
    if(this->m_VecFieldGhost==NULL){
        this->m_VecFieldGhost = (ScalarType*)accfft_alloc(3*g_alloc_max);
    }


    // do the communication for the ghost points
    for (int i = 0; i < 3; i++){
        accfft_get_ghost_xyz(plan,this->m_GhostSize,isize_g,
                                 &this->m_X[i*nl],
                                 &this->m_VecFieldGhost[i*nlghost]);
    }

    if (strcmp(flag.c_str(),"state")==0){

        ierr=Assert(this->m_XS!=NULL,"state X null pointer"); CHKERRQ(ierr);
        ierr=Assert(this->m_StatePlanVec!=NULL,"state X null pointer"); CHKERRQ(ierr);

        this->m_StatePlanVec->interpolate(this->m_VecFieldGhost,3,nx,isize,istart,
                                            nl,this->m_GhostSize,this->m_X,c_dims,
                                            this->m_Opt->GetFFT().mpicomm,timers);

    }
    else if (strcmp(flag.c_str(),"adjoint")==0){

        ierr=Assert(this->m_XA!=NULL,"adjoint X null pointer"); CHKERRQ(ierr);
        ierr=Assert(this->m_AdjointPlanVec!=NULL,"state X null pointer"); CHKERRQ(ierr);

        this->m_AdjointPlanVec->interpolate(this->m_VecFieldGhost,3,nx,isize,istart,
                                            nl,this->m_GhostSize,this->m_X,c_dims,
                                            this->m_Opt->GetFFT().mpicomm,timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }


    for (IntType i = 0; i < nl; ++i){
        wx1[i] = this->m_X[0*nl+i];
        wx2[i] = this->m_X[1*nl+i];
        wx3[i] = this->m_X[2*nl+i];
    }

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IPVEC);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}





/********************************************************************
 * @brief interpolate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangian::Interpolate( ScalarType* wx1,
                                            ScalarType* wx2,
                                            ScalarType* wx3,
                                            ScalarType* vx1,
                                            ScalarType* vx2,
                                            ScalarType* vx3,
                                            ScalarType* yx1,
                                            ScalarType* yx2,
                                            ScalarType* yx3)
{
    PetscErrorCode ierr;
    int nx[3],isize_g[3],isize[3],istart_g[3],istart[3],c_dims[2];
    double timers[4] = {0,0,0,0};
    IntType nl,nlghost,nalloc;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(vx1!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(vx2!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(vx3!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=Assert(wx1!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx2!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx3!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=Assert(yx1!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(yx2!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(yx3!=NULL,"null pointer"); CHKERRQ(ierr);


    nl = this->m_Opt->GetDomainPara().nlocal;

    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<int>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<int>(this->m_Opt->GetDomainPara().isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->GetDomainPara().istart[i]);
    }

    // get network dimensions
    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    if (this->m_X==NULL){
        try{ this->m_X = new double [3*nl]; }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // create planer
    if (this->m_VecFieldPlan==NULL){
        try{ this->m_VecFieldPlan = new Interp3_Plan; }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        this->m_VecFieldPlan->allocate(nl,3);
    }

    for (IntType i = 0; i < nl; ++i){
        this->m_X[i*3+0] = yx1[i]/(2.0*PETSC_PI);
        this->m_X[i*3+1] = yx2[i]/(2.0*PETSC_PI);
        this->m_X[i*3+2] = yx3[i]/(2.0*PETSC_PI);
    }

    // scatter
    this->m_VecFieldPlan->scatter(3,nx,isize,istart,nl,this->m_GhostSize,this->m_X,
                                  c_dims,this->m_Opt->GetFFT().mpicomm,timers);

    for (IntType i = 0; i < nl; ++i){
        this->m_X[0*nl+i] = vx1[i];
        this->m_X[1*nl+i] = vx2[i];
        this->m_X[2*nl+i] = vx3[i];
    }

    // get ghost sizes
    nalloc=accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->GetFFT().plan,
                                            this->m_GhostSize,isize_g,istart_g);

    // get nlocal for ghosts
    nlghost = 1;
    for (IntType i = 0; i < 3; ++i){
        nlghost *= static_cast<IntType>(isize_g[i]);
    }

    // deal with ghost points
    if(this->m_VecFieldGhost==NULL){
        this->m_VecFieldGhost = (ScalarType*)accfft_alloc(3*nalloc);
    }


    // do the communication for the ghost points
    for (int i = 0; i < 3; i++){
        accfft_get_ghost_xyz(this->m_Opt->GetFFT().plan,this->m_GhostSize,isize_g,
                            &this->m_X[i*nl],&this->m_VecFieldGhost[i*nlghost]);
    }

    this->m_VecFieldPlan->interpolate(this->m_VecFieldGhost,3,nx,isize,istart,
                                      nl,this->m_GhostSize,this->m_X,c_dims,
                                      this->m_Opt->GetFFT().mpicomm,timers);

    for (IntType i = 0; i < nl; ++i){
        wx1[i] = this->m_X[0*nl+i];
        wx2[i] = this->m_X[1*nl+i];
        wx3[i] = this->m_X[2*nl+i];
    }

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IPVEC);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief change from lexicographical ordering to xyz
 * @param flag flag to switch between forward and adjoint solves
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MapCoordinateVector"
PetscErrorCode SemiLagrangian::MapCoordinateVector(std::string flag)
{
    PetscErrorCode ierr;
    int nx[3],nlocal,isize[3],istart[3];
    int c_dims[2];
    double timers[4] = {0,0,0,0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // get sizes
    nlocal = static_cast<int>(this->m_Opt->GetDomainPara().nlocal);

    for (int i = 0; i < 3; ++i){
        nx[i]     = static_cast<int>(this->m_Opt->GetDomainPara().nx[i]);
        isize[i]  = static_cast<int>(this->m_Opt->GetDomainPara().isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->GetDomainPara().istart[i]);
    }

    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    if (strcmp(flag.c_str(),"state")==0){

        // characteristic for state equation should
        // have been computed already
        ierr=Assert(this->m_XS != NULL,"null pointer"); CHKERRQ(ierr);

        // create planer
        if (this->m_StatePlan == NULL){
            try{ this->m_StatePlan = new Interp3_Plan; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_StatePlan->allocate(nlocal,1);
        }
        // scatter
        this->m_StatePlan->scatter(1,nx,isize,istart,nlocal,this->m_GhostSize,this->m_XS,
                                    c_dims,this->m_Opt->GetFFT().mpicomm,timers);


        // create planer
        if (this->m_StatePlanVec == NULL){
            try{ this->m_StatePlanVec = new Interp3_Plan; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_StatePlanVec->allocate(nlocal,3);
        }
        // scatter
        this->m_StatePlanVec->scatter(3,nx,isize,istart,nlocal,this->m_GhostSize,this->m_XS,
                                        c_dims,this->m_Opt->GetFFT().mpicomm,timers);

    }
    else if (strcmp(flag.c_str(),"adjoint")==0){

        // characteristic for adjoint equation should
        // have been computed already
        ierr=Assert(this->m_XA != NULL,"null pointer"); CHKERRQ(ierr);

        // create planer
        if (this->m_AdjointPlan == NULL){
            try{ this->m_AdjointPlan = new Interp3_Plan; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_AdjointPlan->allocate(nlocal,1);
        }

        // scatter
        this->m_AdjointPlan->scatter(1,nx,isize,istart,nlocal,this->m_GhostSize,this->m_XA,
                                    c_dims,this->m_Opt->GetFFT().mpicomm,timers);

        // create planer
        if (this->m_AdjointPlanVec == NULL){
            try{ this->m_AdjointPlanVec = new Interp3_Plan; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_AdjointPlanVec->allocate(nlocal,3);
        }

        // scatter
        this->m_AdjointPlanVec->scatter(3,nx,isize,istart,nlocal,this->m_GhostSize,this->m_XA,
                                        c_dims,this->m_Opt->GetFFT().mpicomm,timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

    this->m_Opt->IncreaseInterpTimers(timers);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




} // namespace




#endif // _SEMILAGRANGIAN_CPP_
