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

    this->m_TrajectoryA = NULL;
    this->m_TrajectoryS = NULL;
    this->m_InitialTrajectory = NULL;

    this->m_WorkVecField1 = NULL;
    this->m_WorkVecField2 = NULL;
    this->m_WorkVecField3 = NULL;
    this->m_WorkVecField4 = NULL;

    this->m_XA=NULL;
    this->m_XS=NULL;
    this->m_iVecField=NULL;
    this->m_xVecField=NULL;

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

    if(this->m_WorkVecField1!=NULL){
        delete this->m_WorkVecField1;
        this->m_WorkVecField1 = NULL;
    }

    if(this->m_WorkVecField2!=NULL){
        delete this->m_WorkVecField2;
        this->m_WorkVecField2 = NULL;
    }

    if(this->m_WorkVecField3!=NULL){
        delete this->m_WorkVecField3;
        this->m_WorkVecField3 = NULL;
    }

    if(this->m_WorkVecField4!=NULL){
        delete this->m_WorkVecField4;
        this->m_WorkVecField4 = NULL;
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
 * @brief compute the trajectory from the velocity field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeTrajectory"
PetscErrorCode SemiLagrangian::ComputeTrajectory(VecField* v, std::string flag)
{

    PetscErrorCode ierr;
    ScalarType ht,hthalf;
    VecField* X;

    PetscFunctionBegin;

    if (this->m_InitialTrajectory == NULL){
        ierr=this->ComputeInitialCondition(); CHKERRQ(ierr);
    }

    if (this->m_WorkVecField1==NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
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
    ierr=this->Interpolate(this->m_WorkVecField1,v,flag); CHKERRQ(ierr);

    // compute F0 + F1 = v + v(\tilde{X})
    ierr=VecAXPY(this->m_WorkVecField1->m_X1,1.0,v->m_X1); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField1->m_X2,1.0,v->m_X2); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField1->m_X3,1.0,v->m_X3); CHKERRQ(ierr);

    // compyte X = x - 0.5*ht*(v + v(x - ht v))
    ierr=VecWAXPY(X->m_X1,-hthalf,this->m_WorkVecField1->m_X1,this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr=VecWAXPY(X->m_X2,-hthalf,this->m_WorkVecField1->m_X2,this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr=VecWAXPY(X->m_X3,-hthalf,this->m_WorkVecField1->m_X3,this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);

    ierr=this->MapCoordinateVector(flag); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute the deformation map (we integrate the
 * characteristic backward in time)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMap"
PetscErrorCode SemiLagrangian::ComputeDeformationMap(VecField *y,VecField* v)
{

    PetscErrorCode ierr;
    IntType nl,nt;
    ScalarType ht,hthalf;
    ScalarType *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_vyx1=NULL,*p_vyx2=NULL,*p_vyx3=NULL,
                *p_vytildex1=NULL,*p_vytildex2=NULL,*p_vytildex3=NULL,
                *p_ytildex1=NULL,*p_ytildex2=NULL,*p_ytildex3=NULL,
                *p_yx1=NULL,*p_yx2=NULL,*p_yx3=NULL;

    PetscFunctionBegin;

    ierr=Assert(y != NULL, "input map is null"); CHKERRQ(ierr);
    ierr=Assert(v != NULL, "input velocity field is null"); CHKERRQ(ierr);

    if (this->m_InitialTrajectory == NULL){
        ierr=this->ComputeInitialCondition(); CHKERRQ(ierr);
    }


    if (this->m_WorkVecField1==NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2==NULL){
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3==NULL){
        try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    ierr=y->Copy(this->m_InitialTrajectory); CHKERRQ(ierr);

    ierr=VecGetArray(v->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(v->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(v->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecGetArray(y->m_X1,&p_yx1); CHKERRQ(ierr);
    ierr=VecGetArray(y->m_X2,&p_yx2); CHKERRQ(ierr);
    ierr=VecGetArray(y->m_X3,&p_yx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_ytildex1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_ytildex2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_ytildex3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField2->m_X1,&p_vyx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField2->m_X2,&p_vyx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField2->m_X3,&p_vyx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField3->m_X1,&p_vytildex1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField3->m_X2,&p_vytildex2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField3->m_X3,&p_vytildex3); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j){

        // evaluate v(y)
        ierr=this->Interpolate( p_vyx1,p_vyx2,p_vyx3,
                                p_vx1,p_vx2,p_vx3,
                                p_yx1,p_yx2,p_yx3 ); CHKERRQ(ierr);

        // compute intermediate variable (fist stage of RK2)
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i){

            p_ytildex1[i] = p_yx1[i] - ht*p_vyx1[i];
            p_ytildex2[i] = p_yx2[i] - ht*p_vyx2[i];
            p_ytildex3[i] = p_yx3[i] - ht*p_vyx3[i];

        }
}// end of pragma omp parallel

        // evaluate v(ytilde)
        ierr=this->Interpolate( p_vytildex1,p_vytildex2,p_vytildex3,
                                p_vx1,p_vx2,p_vx3,
                                p_ytildex1,p_ytildex2,p_ytildex3 ); CHKERRQ(ierr);

        // update deformation map (second stage of RK2)
#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i){
            p_yx1[i] = p_yx1[i] - hthalf*(p_vytildex1[i] + p_vyx1[i]);
            p_yx2[i] = p_yx2[i] - hthalf*(p_vytildex2[i] + p_vyx2[i]);
            p_yx3[i] = p_yx3[i] - hthalf*(p_vytildex3[i] + p_vyx3[i]);
        }
}// end of pragma omp parallel


    } // for all time points

    ierr=VecRestoreArray(v->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(v->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(v->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(y->m_X1,&p_yx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(y->m_X2,&p_yx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(y->m_X3,&p_yx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_ytildex1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_ytildex2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_ytildex3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField2->m_X1,&p_vyx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField2->m_X2,&p_vyx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField2->m_X3,&p_vyx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField3->m_X1,&p_vytildex1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField3->m_X2,&p_vytildex2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField3->m_X3,&p_vytildex3); CHKERRQ(ierr);


    PetscFunctionReturn(0);
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
 * @brief interpolate scalar field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscErrorCode SemiLagrangian::Interpolate(ScalarType* w,ScalarType* v,std::string flag)
{
    PetscErrorCode ierr;
    int nx[3],isize_g[3],isize[3],istart_g[3],istart[3],c_dims[2];
    accfft_plan* plan=NULL;
    IntType g_alloc_max,nl;
    double timers[4]={0,0,0,0};

    PetscFunctionBegin;

    ierr=Assert(w!=NULL,"output is null pointer"); CHKERRQ(ierr);
    ierr=Assert(v!=NULL,"input is null pointer"); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<int>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<int>(this->m_Opt->GetISize(i));
        istart[i] = static_cast<int>(this->m_Opt->GetIStart(i));
    }

    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    nl = this->m_Opt->GetNLocal();

    // deal with ghost points
    plan = this->m_Opt->GetFFTPlan();
    g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(plan,this->m_GhostSize,isize_g,istart_g);

    if(this->m_ScaFieldGhost==NULL){
        this->m_ScaFieldGhost = (ScalarType*)accfft_alloc(g_alloc_max);
    }

    accfft_get_ghost_xyz(plan,this->m_GhostSize,isize_g,v,this->m_ScaFieldGhost);

    if (strcmp(flag.c_str(),"state")!=0){

        ierr=Assert(this->m_XS!=NULL,"state X is null pointer"); CHKERRQ(ierr);

        this->m_StatePlan->interpolate(this->m_ScaFieldGhost,1,nx,isize,istart,
                            nl,this->m_GhostSize,w,c_dims,
                            this->m_Opt->GetComm(),timers);
    }
    else if (strcmp(flag.c_str(),"adjoint")!=0){

        ierr=Assert(this->m_XA!=NULL,"adjoint X is null pointer"); CHKERRQ(ierr);
        this->m_AdjointPlan->interpolate(this->m_ScaFieldGhost,1,nx,isize,istart,
                            nl,this->m_GhostSize,w,c_dims,
                            this->m_Opt->GetComm(),timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

//    ierr=this->m_Opt->StopTimer(IPEXEC); CHKERRQ(ierr);

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IP);

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

    ierr=Assert(vx1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vx2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vx3!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx3!=NULL,"input is null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->GetNLocal();

    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<int>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<int>(this->m_Opt->GetISize(i));
        istart[i] = static_cast<int>(this->m_Opt->GetIStart(i));
    }

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

    for (IntType i = 0; i < nl; ++i){
        this->m_iVecField[0*nl + i] = vx1[i];
        this->m_iVecField[1*nl + i] = vx2[i];
        this->m_iVecField[2*nl + i] = vx3[i];
    }

    // get ghost sizes
    plan = this->m_Opt->GetFFTPlan();
    g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(plan,this->m_GhostSize,isize_g,istart_g);

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
                                 &this->m_iVecField[i*nl],
                                 &this->m_VecFieldGhost[i*nlghost]);
    }

    if (strcmp(flag.c_str(),"state")!=0){

        ierr=Assert(this->m_XS!=NULL,"state X null pointer"); CHKERRQ(ierr);
        ierr=Assert(this->m_StatePlanVec!=NULL,"state X null pointer"); CHKERRQ(ierr);


        this->m_StatePlanVec->interpolate(this->m_VecFieldGhost,3,nx,isize,istart,
                                        nl,this->m_GhostSize,this->m_xVecField,c_dims,
                                        this->m_Opt->GetComm(),timers);


    }
    else if (strcmp(flag.c_str(),"adjoint")!=0){

        ierr=Assert(this->m_XA!=NULL,"adjoint X null pointer"); CHKERRQ(ierr);
        ierr=Assert(this->m_AdjointPlanVec!=NULL,"state X null pointer"); CHKERRQ(ierr);

        this->m_AdjointPlanVec->interpolate(this->m_VecFieldGhost,3,nx,isize,istart,
                                        nl,this->m_GhostSize,this->m_xVecField,c_dims,
                                        this->m_Opt->GetComm(),timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }


    for (IntType i = 0; i < nl; ++i){
        wx1[i] = this->m_xVecField[0*nl + i];
        wx2[i] = this->m_xVecField[1*nl + i];
        wx3[i] = this->m_xVecField[2*nl + i];
    }

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IPVEC);

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
    accfft_plan* plan=NULL;
    IntType nl,nlghost,g_alloc_max;
    PetscFunctionBegin;

    ierr=Assert(vx1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vx2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(vx3!=NULL,"input is null pointer"); CHKERRQ(ierr);

    ierr=Assert(wx1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(wx3!=NULL,"input is null pointer"); CHKERRQ(ierr);

    ierr=Assert(yx1!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(yx2!=NULL,"input is null pointer"); CHKERRQ(ierr);
    ierr=Assert(yx3!=NULL,"input is null pointer"); CHKERRQ(ierr);


    nl = this->m_Opt->GetNLocal();

    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<int>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<int>(this->m_Opt->GetISize(i));
        istart[i] = static_cast<int>(this->m_Opt->GetIStart(i));
    }

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

    // create planer
    if (this->m_VecFieldPlan == NULL){
        try{ this->m_VecFieldPlan = new Interp3_Plan; }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        this->m_VecFieldPlan->allocate(nl,3);
    }

    for (IntType i = 0; i < nl; ++i){
        this->m_iVecField[i*3+0] = yx1[i]/(2.0*PETSC_PI);
        this->m_iVecField[i*3+1] = yx2[i]/(2.0*PETSC_PI);
        this->m_iVecField[i*3+2] = yx3[i]/(2.0*PETSC_PI);
    }

    // scatter
    this->m_VecFieldPlan->scatter(3,nx,isize,istart,nl,
                                  this->m_GhostSize,this->m_iVecField,
                                  c_dims,this->m_Opt->GetComm(),timers);

    for (IntType i = 0; i < nl; ++i){
        this->m_iVecField[0*nl + i] = vx1[i];
        this->m_iVecField[1*nl + i] = vx2[i];
        this->m_iVecField[2*nl + i] = vx3[i];
    }

    // get ghost sizes
    plan = this->m_Opt->GetFFTPlan();
    g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c(plan,this->m_GhostSize,isize_g,istart_g);

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
                                 &this->m_iVecField[i*nl],
                                 &this->m_VecFieldGhost[i*nlghost]);
    }

    this->m_VecFieldPlan->interpolate(this->m_VecFieldGhost,3,nx,isize,istart,
                                      nl,this->m_GhostSize,this->m_xVecField,c_dims,
                                      this->m_Opt->GetComm(),timers);

    for (IntType i = 0; i < nl; ++i){
        wx1[i] = this->m_xVecField[0*nl + i];
        wx2[i] = this->m_xVecField[1*nl + i];
        wx3[i] = this->m_xVecField[2*nl + i];
    }

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IPVEC);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief compute initial condition for trajectory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeInitialCondition"
PetscErrorCode SemiLagrangian::ComputeInitialCondition()
{
    PetscErrorCode ierr;
    ScalarType *p_x1=NULL,*p_x2=NULL,*p_x3=NULL;
    IntType isize[3],istart[3];
    ScalarType hx[3];

    PetscFunctionBegin;

    // allocate initial trajectory
    if(this->m_InitialTrajectory==NULL){
        try{this->m_InitialTrajectory = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("slm: computing initial condition"); CHKERRQ(ierr);
    }

    for (int i = 0; i < 3; ++i){
        hx[i] = this->m_Opt->GetSpatialStepSize(i);
        isize[i] = static_cast<int>(this->m_Opt->GetISize(i));
        istart[i] = static_cast<int>(this->m_Opt->GetIStart(i));
    }

    ierr=VecGetArray(this->m_InitialTrajectory->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_InitialTrajectory->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_InitialTrajectory->m_X3,&p_x3); CHKERRQ(ierr);

#pragma omp parallel
{
    IntType li;
    ScalarType x1,x2,x3;
#pragma omp for
    for (IntType i1 = 0; i1 < isize[0]; ++i1){  // x1
        for (IntType i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < isize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                li = GetLinearIndex(i1,i2,i3,isize);
                // assign values
                p_x1[li] = x1;
                p_x2[li] = x2;
                p_x3[li] = x3;

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
 * @brief change from lexicographical ordering to xyz
 * @param flag flag to switch between forward and adjoint solves
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MapCoordinateVector"
PetscErrorCode SemiLagrangian::MapCoordinateVector(std::string flag)
{
    PetscErrorCode ierr;
    const ScalarType *p_x1=NULL,*p_x2=NULL,*p_x3=NULL;
    int nx[3],nl,isize[3],istart[3];
    int c_dims[2];
    double timers[4] = {0,0,0,0};

    PetscFunctionBegin;

    for (int i = 0; i < 3; ++i){
        nx[i] = static_cast<int>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<int>(this->m_Opt->GetISize(i));
        istart[i] = static_cast<int>(this->m_Opt->GetIStart(i));
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
            try{ this->m_StatePlan = new Interp3_Plan; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_StatePlan->allocate(nl,1);
        }
        // scatter
        this->m_StatePlan->scatter(1,nx,isize,istart,nl,
                                    this->m_GhostSize,this->m_XS,
                                    c_dims,this->m_Opt->GetComm(),timers);


        // create planer
        if (this->m_StatePlanVec == NULL){
            try{ this->m_StatePlanVec = new Interp3_Plan; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_StatePlanVec->allocate(nl,3);
        }
        // scatter
        this->m_StatePlanVec->scatter(3,nx,isize,istart,nl,
                                        this->m_GhostSize,this->m_XS,
                                        c_dims,this->m_Opt->GetComm(),timers);



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
            try{ this->m_AdjointPlan = new Interp3_Plan; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_AdjointPlan->allocate(nl,1);
        }

        // scatter
        this->m_AdjointPlan->scatter(1,nx,isize,istart,nl,
                                    this->m_GhostSize,this->m_XA,
                                    c_dims,this->m_Opt->GetComm(),timers);

        // create planer
        if (this->m_AdjointPlanVec == NULL){
            try{ this->m_AdjointPlanVec = new Interp3_Plan; }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            this->m_AdjointPlanVec->allocate(nl,3);
        }

        // scatter
        this->m_AdjointPlanVec->scatter(3,nx,isize,istart,nl,
                                    this->m_GhostSize,this->m_XA,
                                    c_dims,this->m_Opt->GetComm(),timers);

    }
    else { ierr=ThrowError("flag wrong"); CHKERRQ(ierr); }

    this->m_Opt->IncreaseInterpTimers(timers);

    PetscFunctionReturn(0);
}



} // namespace

#endif // _SEMILAGRANGIAN_CPP_
