#ifndef _SEMILAGRANGIANGPU_CPP_
#define _SEMILAGRANGIANGPU_CPP_

#include "SemiLagrangianGPU.h"


namespace reg
{




/********************************************************************
 * Name: SemiLagrangianGPU
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SemiLagrangianGPU"
SemiLagrangianGPU::SemiLagrangianGPU()
{
    this->Initialize();
}



/********************************************************************
 * Name: SemiLagrangianGPU
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SemiLagrangianGPU"
SemiLagrangianGPU::SemiLagrangianGPU(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}



/********************************************************************
 * Name: SemiLagrangianGPU
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~SemiLagrangianGPU"
SemiLagrangianGPU::~SemiLagrangianGPU()
{
    this->ClearMemory();
}




/********************************************************************
 * Name: Initialize
 * Description: init class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode SemiLagrangianGPU::Initialize()
{
    PetscFunctionBegin;

    this->m_TrajectoryA = NULL;
    this->m_TrajectoryS = NULL;
    this->m_InitialTrajectory = NULL;
    this->m_WorkVecField = NULL;

    this->m_XA=NULL;
    this->m_XS=NULL;
    this->m_iVecField=NULL;
    this->m_xVecField=NULL;

    this->m_AdjointPlan = NULL;
    this->m_AdjointPlanVec = NULL;

    this->m_StatePlan = NULL;
    this->m_StatePlanVec = NULL;

    this->m_ScaFieldGhost=NULL;
    this->m_VecFieldGhost=NULL;

    this->m_Opt=NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clears memory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode SemiLagrangianGPU::ClearMemory()
{

    PetscFunctionBegin;

    if(this->m_TrajectoryS != NULL){
        delete this->m_TrajectoryS;
        this->m_TrajectoryS = NULL;
    }

    if(this->m_TrajectoryA != NULL){
        delete this->m_TrajectoryA;
        this->m_TrajectoryA = NULL;
    }

    if(this->m_InitialTrajectory!=NULL){
        delete this->m_InitialTrajectory;
        this->m_InitialTrajectory = NULL;
    }

    if(this->m_WorkVecField!=NULL){
        delete this->m_WorkVecField;
        this->m_WorkVecField = NULL;
    }

    if(this->m_xVecField!=NULL){
        delete this->m_xVecField;
        this->m_xVecField = NULL;
    }

    if(this->m_iVecField!=NULL){
        delete this->m_iVecField;
        this->m_iVecField = NULL;
    }

    if(this->m_XS!=NULL){
        delete this->m_XS;
        this->m_XS = NULL;
    }
    if(this->m_XA!=NULL){
        delete this->m_XA;
        this->m_XA = NULL;
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

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeTrajectory
 * Description: compute the trajectory from the velocity field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeTrajectory"
PetscErrorCode SemiLagrangianGPU::ComputeTrajectory(VecField* v, std::string flag)
{

    PetscErrorCode ierr;
    ScalarType ht,hthalf;
    VecField* X;

    PetscFunctionBegin;

    if (this->m_InitialTrajectory == NULL){
        ierr=this->ComputeInitialCondition(); CHKERRQ(ierr);
    }

    if (this->m_WorkVecField==NULL){
        try{this->m_WorkVecField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (strcmp(flag.c_str(),"state")!=0){

        if (this->m_TrajectoryS==NULL){
            try{this->m_TrajectoryS = new VecField(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        X = this->m_TrajectoryS;

    }
    else if (strcmp(flag.c_str(),"adjoint")!=0){

        if (this->m_TrajectoryA==NULL){
            try{this->m_TrajectoryA = new VecField(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        X = this->m_TrajectoryA;

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    // \tilde{X} = x - ht v
    ierr=VecWAXPY(X->m_X1,-ht,v->m_X1,this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr=VecWAXPY(X->m_X2,-ht,v->m_X2,this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr=VecWAXPY(X->m_X3,-ht,v->m_X3,this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);

    ierr=this->MapCoordinateVector(flag); CHKERRQ(ierr);
    ierr=this->Interpolate(this->m_WorkVecField,v,flag); CHKERRQ(ierr);

    // compute F0 + F1 = v + v(\tilde{X})
    ierr=VecAXPY(this->m_WorkVecField->m_X1,1.0,v->m_X1); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField->m_X2,1.0,v->m_X2); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField->m_X3,1.0,v->m_X3); CHKERRQ(ierr);

    // compyte X = x - 0.5*ht*(v + v(x - ht v))
    ierr=VecWAXPY(X->m_X1,-hthalf,this->m_WorkVecField->m_X1,this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr=VecWAXPY(X->m_X2,-hthalf,this->m_WorkVecField->m_X2,this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr=VecWAXPY(X->m_X3,-hthalf,this->m_WorkVecField->m_X3,this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);

    ierr=this->MapCoordinateVector(flag); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Interpolate
 * Description: interpolate scalar field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangianGPU::Interpolate(Vec* w,Vec v,std::string flag)
{
    PetscErrorCode ierr;
    ScalarType *p_w=NULL,*p_v=NULL;

    PetscFunctionBegin;
    ierr=Assert(*w!=NULL,"output is null pointer"); CHKERRQ(ierr);
    ierr=Assert(v!=NULL,"input is null pointer"); CHKERRQ(ierr);

    ierr=VecGetArray(v,&p_v); CHKERRQ(ierr);
    ierr=VecGetArray(*w,&p_w); CHKERRQ(ierr);

    ierr=this->Interpolate(p_w,p_v,flag); CHKERRQ(ierr);

    ierr=VecRestoreArray(*w,&p_w); CHKERRQ(ierr);
    ierr=VecRestoreArray(v,&p_v); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: Interpolate
 * Description: interpolate scalar field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangianGPU::Interpolate(ScalarType* w,ScalarType* v,std::string flag)
{
    PetscErrorCode ierr;
    int nx[3],isize_g[3],istart_g[3],c_dims[2];
    accfft_plan* plan=NULL;
    IntType g_alloc_max;
    double timers[4]={0,0,0,0};
    unsigned long long nl=0;

    PetscFunctionBegin;

    ierr=Assert(w!=NULL,"output is null pointer"); CHKERRQ(ierr);
    ierr=Assert(v!=NULL,"input is null pointer"); CHKERRQ(ierr);

    nx[0] = this->m_Opt->m_MiscOpt->N[0];
    nx[1] = this->m_Opt->m_MiscOpt->N[1];
    nx[2] = this->m_Opt->m_MiscOpt->N[2];

    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    nl = this->m_Opt->GetNLocal();

    // deal with ghost points
    plan = this->m_Opt->m_MiscOpt->plan;
    g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(plan,this->m_GhostSize,isize_g,istart_g);

    if(this->m_ScaFieldGhost==NULL){
        this->m_ScaFieldGhost = (ScalarType*)accfft_alloc(g_alloc_max);
    }

    accfft_get_ghost_xyz(plan,this->m_GhostSize,isize_g,v,this->m_ScaFieldGhost);

    if (strcmp(flag.c_str(),"state")!=0){

        ierr=Assert(this->m_XS!=NULL,"state X is null pointer"); CHKERRQ(ierr);

        this->m_StatePlan->interpolate(this->m_ScaFieldGhost,1,nx,
                            this->m_Opt->m_MiscOpt->isize,
                            this->m_Opt->m_MiscOpt->istart,
                            nl,this->m_GhostSize,w,c_dims,
                            this->m_Opt->m_MiscOpt->c_comm,timers);
    }
    else if (strcmp(flag.c_str(),"adjoint")!=0){

        ierr=Assert(this->m_XA!=NULL,"adjoint X is null pointer"); CHKERRQ(ierr);
        this->m_AdjointPlan->interpolate(this->m_ScaFieldGhost,1,nx,
                            this->m_Opt->m_MiscOpt->isize,
                            this->m_Opt->m_MiscOpt->istart,
                            nl,this->m_GhostSize,w,c_dims,
                            this->m_Opt->m_MiscOpt->c_comm,timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

//    ierr=this->m_Opt->StopTimer(IPEXEC); CHKERRQ(ierr);

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IP);

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: Interpolate
 * Description: interpolate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangianGPU::Interpolate(VecField* w, VecField* v, std::string flag)
{
    PetscErrorCode ierr;
    ScalarType *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_wx1=NULL,*p_wx2=NULL,*p_wx3=NULL;
    PetscFunctionBegin;

    ierr=Assert(v!=NULL,"output is null pointer"); CHKERRQ(ierr);
    ierr=Assert(w!=NULL,"input is null pointer"); CHKERRQ(ierr);

    ierr=VecGetArray(w->m_X1,&p_wx1); CHKERRQ(ierr);
    ierr=VecGetArray(w->m_X2,&p_wx2); CHKERRQ(ierr);
    ierr=VecGetArray(w->m_X3,&p_wx3); CHKERRQ(ierr);

    ierr=VecGetArray(v->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(v->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(v->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=this->Interpolate(p_wx1,p_wx2,p_wx3,
                           p_vx1,p_vx2,p_vx3,flag); CHKERRQ(ierr);

    ierr=VecRestoreArray(v->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(v->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(v->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(w->m_X1,&p_wx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(w->m_X2,&p_wx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(w->m_X3,&p_wx3); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: Interpolate
 * Description: interpolate vector field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangianGPU::Interpolate( ScalarType* vox1,
                                                      ScalarType* vox2,
                                                      ScalarType* vox3,
                                                      ScalarType* vix1,
                                                      ScalarType* vix2,
                                                      ScalarType* vix3,
                                                      std::string flag)
{
    PetscErrorCode ierr;
    int nx[3],isize_g[3],istart_g[3],c_dims[2];
    double timers[4] = {0,0,0,0};
    accfft_plan* plan=NULL;
    IntType g_alloc_max;
    IntType nl,nlghost;
    PetscFunctionBegin;

    ierr=Assert(vix1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vix2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vix3!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vox1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vox2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vox3!=NULL,"input is null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->GetNLocal();
    nx[0] = this->m_Opt->m_MiscOpt->N[0];
    nx[1] = this->m_Opt->m_MiscOpt->N[1];
    nx[2] = this->m_Opt->m_MiscOpt->N[2];

    // get network dimensions
    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    if (this->m_iVecField==NULL){
        try{ this->m_iVecField = new double [3*nl]; }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_xVecField == NULL){
        try{ this->m_xVecField = new double [3*nl]; }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    for (unsigned long i = 0; i < nl; ++i){
        this->m_iVecField[0*nl + i] = vix1[i];
        this->m_iVecField[1*nl + i] = vix2[i];
        this->m_iVecField[2*nl + i] = vix3[i];
    }

    // get ghost sizes
    plan = this->m_Opt->m_MiscOpt->plan;
    g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(plan,this->m_GhostSize,isize_g,istart_g);

    // get nlocal for ghosts
    nlghost = 1;
    for (unsigned int i = 0; i < 3; ++i){
        nlghost *= static_cast<unsigned long>(isize_g[i]);
    }

    // deal with ghost points
    if(this->m_VecFieldGhost==NULL){
        this->m_VecFieldGhost = (ScalarType*)accfft_alloc(3*g_alloc_max);
    }


    // do the communication for the ghost points
    for (unsigned int i = 0; i < 3; i++){
        accfft_get_ghost_xyz(plan,this->m_GhostSize,isize_g,
                                 &this->m_iVecField[i*nl],
                                 &this->m_VecFieldGhost[i*nlghost]);
    }

    if (strcmp(flag.c_str(),"state")!=0){

        ierr=Assert(this->m_XS!=NULL,"state X null pointer"); CHKERRQ(ierr);
        ierr=Assert(this->m_StatePlanVec!=NULL,"state X null pointer"); CHKERRQ(ierr);


        this->m_StatePlanVec->interpolate(this->m_VecFieldGhost,3,nx,
                                        this->m_Opt->m_MiscOpt->isize,
                                        this->m_Opt->m_MiscOpt->istart,
                                        nl,this->m_GhostSize,this->m_xVecField,c_dims,
                                        this->m_Opt->m_MiscOpt->c_comm,timers);


    }
    else if (strcmp(flag.c_str(),"adjoint")!=0){

        ierr=Assert(this->m_XA!=NULL,"adjoint X null pointer"); CHKERRQ(ierr);
        ierr=Assert(this->m_AdjointPlanVec!=NULL,"state X null pointer"); CHKERRQ(ierr);

        this->m_AdjointPlanVec->interpolate(this->m_VecFieldGhost,3,nx,
                                        this->m_Opt->m_MiscOpt->isize,
                                        this->m_Opt->m_MiscOpt->istart,
                                        nl,this->m_GhostSize,this->m_xVecField,c_dims,
                                        this->m_Opt->m_MiscOpt->c_comm,timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

//    ierr=this->m_Opt->StopTimer(IPVECEXEC); CHKERRQ(ierr);

    for (unsigned long i = 0; i < nl; ++i){
        vox1[i] = this->m_xVecField[0*nl + i];
        vox2[i] = this->m_xVecField[1*nl + i];
        vox3[i] = this->m_xVecField[2*nl + i];
    }

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IPVEC);

    PetscFunctionReturn(0);

}


/********************************************************************
 * Name: ComputeInitialCondition
 * Description: compute initial condition for trajectory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeInitialCondition"
PetscErrorCode SemiLagrangianGPU::ComputeInitialCondition()
{
    PetscErrorCode ierr;
    ScalarType *p_x1=NULL,*p_x2=NULL,*p_x3=NULL;
    unsigned int isize[3],istart[3];
    ScalarType hx[3];
    PetscFunctionBegin;

    if(this->m_InitialTrajectory==NULL){
        try{this->m_InitialTrajectory = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    for (unsigned int i = 0; i < 3; ++i){
        hx[i]     = this->m_Opt->GetSpatialStepSize(i);
        isize[i]  = this->m_Opt->GetISize(i);
        istart[i] = this->m_Opt->GetIStart(i);
    }

    ierr=VecGetArray(this->m_InitialTrajectory->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_InitialTrajectory->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_InitialTrajectory->m_X3,&p_x3); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (unsigned int i1 = 0; i1 < isize[0]; ++i1){  // x1
        for (unsigned int i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (unsigned int i3 = 0; i3 < isize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                ScalarType x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                ScalarType x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                ScalarType x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                IntType linidx = GetLinearIndex(i1,i2,i3,isize);

                // assign values
                p_x1[linidx] = x1;
                p_x2[linidx] = x2;
                p_x3[linidx] = x3;

            } // i1
        } // i2
    } // i3
}// pragma omp for

    ierr=VecRestoreArray(this->m_InitialTrajectory->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_InitialTrajectory->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_InitialTrajectory->m_X3,&p_x3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: map coordinates
 * Description: change from lexicographical ordering to xyz
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MapCoordinateVector"
PetscErrorCode SemiLagrangianGPU::MapCoordinateVector(std::string flag)
{
    PetscErrorCode ierr;
    const ScalarType *p_x1=NULL,*p_x2=NULL,*p_x3=NULL;
    int nx[3],isize_g[3],istart_g[3],c_dims[2],isize[3],istart[3];
    accfft_plan* plan=NULL;
    IntType g_alloc_max;
    IntType nl;
    VecField *X=NULL;
    double timers[4] = {0,0,0,0};

    PetscFunctionBegin;

    for (int i = 0; i < 3; ++i){
        nx[i] = this->m_Opt->m_MiscOpt->N[i];
        isize[i] = this->m_Opt->m_MiscOpt->isize[i];
        istart[i] = this->m_Opt->m_MiscOpt->istart[i];
    }
    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    nl = this->m_Opt->GetNLocal();

    if (strcmp(flag.c_str(),"state")!=0){

        if (this->m_XS == NULL){
            try{ this->m_XS = new double [3*nl]; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        ierr=Assert(this->m_TrajectoryS != NULL,"trajectory S null pointer"); CHKERRQ(ierr);
        ierr=VecGetArrayRead(this->m_TrajectoryS->m_X1,&p_x1); CHKERRQ(ierr);
        ierr=VecGetArrayRead(this->m_TrajectoryS->m_X2,&p_x2); CHKERRQ(ierr);
        ierr=VecGetArrayRead(this->m_TrajectoryS->m_X3,&p_x3); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i){
            this->m_XS[i*3+0] = p_x1[i]/(2.0*PETSC_PI);
            this->m_XS[i*3+1] = p_x2[i]/(2.0*PETSC_PI);
            this->m_XS[i*3+2] = p_x3[i]/(2.0*PETSC_PI);
        }
} // pragma omp parallel

        ierr=VecRestoreArrayRead(this->m_TrajectoryS->m_X1,&p_x1); CHKERRQ(ierr);
        ierr=VecRestoreArrayRead(this->m_TrajectoryS->m_X2,&p_x2); CHKERRQ(ierr);
        ierr=VecRestoreArrayRead(this->m_TrajectoryS->m_X3,&p_x3); CHKERRQ(ierr);

        // create planer
        if (this->m_StatePlan == NULL){
            g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->m_MiscOpt->plan,
                                                            this->m_GhostSize,isize_g,istart_g);
            try{ this->m_StatePlan = new Interp3_Plan_GPU(g_alloc_max); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_StatePlan->allocate(nl,1);
        }
        // scatter
        this->m_StatePlan->scatter(1,nx,this->m_Opt->m_MiscOpt->isize,
                                    this->m_Opt->m_MiscOpt->istart,nl,
                                    this->m_GhostSize,this->m_XS,
                                    c_dims,this->m_Opt->m_MiscOpt->c_comm,timers);


        // create planer
        if (this->m_StatePlanVec == NULL){
            g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->m_MiscOpt->plan,
                                                            this->m_GhostSize,isize_g,istart_g);
            try{ this->m_StatePlanVec = new Interp3_Plan_GPU(g_alloc_max); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_StatePlanVec->allocate(nl,3);
        }
        // scatter
        this->m_StatePlanVec->scatter(3,nx,this->m_Opt->m_MiscOpt->isize,
                                        this->m_Opt->m_MiscOpt->istart,nl,
                                        this->m_GhostSize,this->m_XS,
                                        c_dims,this->m_Opt->m_MiscOpt->c_comm,timers);



    }
    else if (strcmp(flag.c_str(),"adjoint")!=0){

        if (this->m_XA == NULL){
            try{ this->m_XA = new double [3*nl]; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        ierr=Assert(this->m_TrajectoryA != NULL,"trajectory A is null pointer"); CHKERRQ(ierr);
        ierr=VecGetArrayRead(this->m_TrajectoryA->m_X1,&p_x1); CHKERRQ(ierr);
        ierr=VecGetArrayRead(this->m_TrajectoryA->m_X2,&p_x2); CHKERRQ(ierr);
        ierr=VecGetArrayRead(this->m_TrajectoryA->m_X3,&p_x3); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i){
            this->m_XA[i*3+0] = p_x1[i]/(2.0*PETSC_PI);
            this->m_XA[i*3+1] = p_x2[i]/(2.0*PETSC_PI);
            this->m_XA[i*3+2] = p_x3[i]/(2.0*PETSC_PI);
        }
} // pragma omp parallel

        ierr=VecRestoreArrayRead(this->m_TrajectoryA->m_X1,&p_x1); CHKERRQ(ierr);
        ierr=VecRestoreArrayRead(this->m_TrajectoryA->m_X2,&p_x2); CHKERRQ(ierr);
        ierr=VecRestoreArrayRead(this->m_TrajectoryA->m_X3,&p_x3); CHKERRQ(ierr);

        // create planer
        if (this->m_AdjointPlan == NULL){
            g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->m_MiscOpt->plan,
                                                            this->m_GhostSize,isize_g,istart_g);
            try{ this->m_AdjointPlan = new Interp3_Plan_GPU(g_alloc_max); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_AdjointPlan->allocate(nl,1);
        }

        // scatter
        this->m_AdjointPlan->scatter(1,nx,this->m_Opt->m_MiscOpt->isize,
                                    this->m_Opt->m_MiscOpt->istart,nl,
                                    this->m_GhostSize,this->m_XA,
                                    c_dims,this->m_Opt->m_MiscOpt->c_comm,timers);

        // create planer
        if (this->m_AdjointPlanVec == NULL){
            g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->m_MiscOpt->plan,
                                                            this->m_GhostSize,isize_g,istart_g);
            try{ this->m_AdjointPlanVec = new Interp3_Plan_GPU(g_alloc_max); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_AdjointPlanVec->allocate(nl,3);
        }

        // scatter
        this->m_AdjointPlanVec->scatter(3,nx,this->m_Opt->m_MiscOpt->isize,
                                    this->m_Opt->m_MiscOpt->istart,nl,
                                    this->m_GhostSize,this->m_XA,
                                    c_dims,this->m_Opt->m_MiscOpt->c_comm,timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

    this->m_Opt->IncreaseInterpTimers(timers);

    PetscFunctionReturn(0);
}



} // namespace

#endif // _SEMILAGRANGIAN_CPP_
