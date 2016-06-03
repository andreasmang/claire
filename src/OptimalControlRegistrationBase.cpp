#ifndef _OPTIMALCONTROLREGISTRATIONBASE_CPP_
#define _OPTIMALCONTROLREGISTRATIONBASE_CPP_

#include "OptimalControlRegistrationBase.hpp"




namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistrationBase"
OptimalControlRegistrationBase::OptimalControlRegistrationBase() : SuperClass()
{
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~OptimalControlRegistrationBase"
OptimalControlRegistrationBase::~OptimalControlRegistrationBase(void)
{
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistrationBase"
OptimalControlRegistrationBase::OptimalControlRegistrationBase(RegOpt* opt) : SuperClass(opt)
{
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode OptimalControlRegistrationBase::Initialize(void)
{
    PetscFunctionBegin;

    // pointer for container of velocity field
    this->m_VelocityField = NULL;

    // pointer for container of incremental velocity field
    this->m_IncVelocityField = NULL;

    // pointers to images
    this->m_TemplateImage = NULL;
    this->m_ReferenceImage = NULL;

    this->m_IO=NULL; ///< read / write object
    this->m_SL=NULL; ///< semi lagranigan

    this->m_WorkScaField1 = NULL;
    this->m_WorkScaField2 = NULL;
    this->m_WorkScaField3 = NULL;
    this->m_WorkScaField4 = NULL;

    this->m_WorkVecField1 = NULL;
    this->m_WorkVecField2 = NULL;
    this->m_WorkVecField3 = NULL;
    this->m_WorkVecField4 = NULL;

    this->m_VelocityIsZero = false; // flag for velocity field

    this->m_Regularization = NULL; // pointer for regularization class

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode OptimalControlRegistrationBase::ClearMemory(void)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_VelocityField != NULL){
        delete this->m_VelocityField;
        this->m_VelocityField = NULL;
    }

    if (this->m_IncVelocityField != NULL){
        delete this->m_IncVelocityField;
        this->m_IncVelocityField = NULL;
    }

    if (this->m_Regularization != NULL){
        delete this->m_Regularization;
        this->m_Regularization = NULL;
    }

    if (this->m_WorkScaField1 != NULL){
        ierr=VecDestroy(&this->m_WorkScaField1); CHKERRQ(ierr);
        this->m_WorkScaField1 = NULL;
    }

    if (this->m_WorkScaField2 != NULL){
        ierr=VecDestroy(&this->m_WorkScaField2); CHKERRQ(ierr);
        this->m_WorkScaField2 = NULL;
    }

    if (this->m_WorkScaField3 != NULL){
        ierr=VecDestroy(&this->m_WorkScaField3); CHKERRQ(ierr);
        this->m_WorkScaField3 = NULL;
    }

    if (this->m_WorkScaField4 != NULL){
        ierr=VecDestroy(&this->m_WorkScaField4); CHKERRQ(ierr);
        this->m_WorkScaField4 = NULL;
    }

    if (this->m_WorkVecField1 != NULL){
        delete this->m_WorkVecField1;
        this->m_WorkVecField1 = NULL;
    }

    if (this->m_WorkVecField2 != NULL){
        delete this->m_WorkVecField2;
        this->m_WorkVecField2 = NULL;
    }

    if (this->m_WorkVecField3 != NULL){
        delete this->m_WorkVecField3;
        this->m_WorkVecField3 = NULL;
    }

    if (this->m_WorkVecField4 != NULL){
        delete this->m_WorkVecField4;
        this->m_WorkVecField4 = NULL;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetIO"
PetscErrorCode OptimalControlRegistrationBase::SetIO(ReadWriteReg* io)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(io != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_IO = io;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReferenceImage"
PetscErrorCode OptimalControlRegistrationBase::SetReferenceImage(Vec mR)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(mR != NULL, "input reference image is null pointer"); CHKERRQ(ierr);
    this->m_ReferenceImage = mR;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set template image (i.e., the image to be registered)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetTemplateImage"
PetscErrorCode OptimalControlRegistrationBase::SetTemplateImage(Vec mT)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(mT != NULL, "input template image is null pointer"); CHKERRQ(ierr);
    this->m_TemplateImage = mT;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set velocity field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetVelocityField"
PetscErrorCode OptimalControlRegistrationBase::SetVelocityField(Vec v)
{
    PetscErrorCode ierr;
    IntType nl;
    PetscInt nlv;
    PetscFunctionBegin;

    // check for null pointer
    ierr=Assert(v != NULL, "input velocity field is null pointer"); CHKERRQ(ierr);

    // check input vector size
    nl = this->m_Opt->GetNLocal();
    ierr=VecGetLocalSize(v,&nlv); CHKERRQ(ierr);
    ierr=Assert(3*nl == nlv, "size mismatch"); CHKERRQ(ierr);

    // allocate velocity field
    if(this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // copy buffer
    ierr=this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set the registration options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetVelocity2Zero"
PetscErrorCode OptimalControlRegistrationBase::SetVelocity2Zero()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if(this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr=this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief check if velocity field is zero
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "IsVelocityZero"
PetscErrorCode OptimalControlRegistrationBase::IsVelocityZero()
{
    PetscErrorCode ierr;
    ScalarType normv1,normv2,normv3;
    PetscFunctionBegin;

    this->m_VelocityIsZero = false;
    ierr=Assert(this->m_VelocityField!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=VecNorm(this->m_VelocityField->m_X1,NORM_INFINITY,&normv1); CHKERRQ(ierr);
    ierr=VecNorm(this->m_VelocityField->m_X2,NORM_INFINITY,&normv2); CHKERRQ(ierr);
    ierr=VecNorm(this->m_VelocityField->m_X3,NORM_INFINITY,&normv3); CHKERRQ(ierr);

    this->m_VelocityIsZero = (normv1 == 0.0) && (normv2 == 0.0) && (normv3 == 0.0);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute a synthetic test problem
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupSyntheticProb"
PetscErrorCode OptimalControlRegistrationBase::SetupSyntheticProb()
{
    PetscErrorCode ierr;
    IntType nl,ng;
    IntType isize[3],istart[3];
    ScalarType *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,*p_mt=NULL,hx[3];
    int problem=3;

    PetscFunctionBegin;

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("setting up synthetic test problem"); CHKERRQ(ierr);
    }
    nl = this->m_Opt->GetNLocal();
    ng = this->m_Opt->GetNGlobal();

    for (int i = 0; i < 3; ++i){
        hx[i] = this->m_Opt->GetSpatialStepSize(i);
        isize[i]  = this->m_Opt->GetISize(i);
        istart[i] = this->m_Opt->GetIStart(i);
    }

    // allocate vector fields
    if(this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    if(this->m_Opt->GetRegModel() == STOKES){
        problem=4;
    }

    // allocate reference image
    if(this->m_ReferenceImage == NULL){
        // create an extra array for initial guess (has to be flat for optimizer)
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_ReferenceImage); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_ReferenceImage,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_ReferenceImage); CHKERRQ(ierr);
        ierr=VecSet(this->m_ReferenceImage,0.0); CHKERRQ(ierr);
    }

    // allocate template image
    if(this->m_TemplateImage == NULL){
        // create an extra array for initial guess (has to be flat for optimizer)
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_TemplateImage); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_TemplateImage,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_TemplateImage); CHKERRQ(ierr);
        ierr=VecSet(this->m_TemplateImage,0.0); CHKERRQ(ierr);
    }


    ierr=VecGetArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_TemplateImage,&p_mt); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (IntType i1 = 0; i1 < isize[0]; ++i1){  // x1
        for (IntType i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < isize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                ScalarType x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                ScalarType x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                ScalarType x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                IntType i = GetLinearIndex(i1,i2,i3,isize);
                p_mt[i] = (sin(x1)*sin(x1) + sin(x2)*sin(x2) + sin(x3)*sin(x3))/3.0;

                if (problem == 0){
                    ScalarType v0 = 0.5;
                    // compute the velocity field
                    p_vx1[i] = v0*sin(x3)*cos(x2)*sin(x2);
                    p_vx2[i] = v0*sin(x1)*cos(x3)*sin(x3);
                    p_vx3[i] = v0*sin(x2)*cos(x1)*sin(x1);
                }
                else if (problem == 1){
                    // compute the velocity field
                    p_vx1[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                    p_vx2[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                    p_vx3[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                }
                else if (problem == 3){
                    p_vx1[i] = cos(x1)*sin(x2);
                    p_vx2[i] = cos(x2)*sin(x1);
                    p_vx3[i] = cos(x1)*sin(x3);
                }
                else if (problem == 4){
                    p_vx1[i] = cos(x2)*cos(x3);
                    p_vx2[i] = sin(x3)*sin(x1);
                    p_vx3[i] = cos(x1)*cos(x2);
                }
            } // i1
        } // i2
    } // i3
} // pragma omp parallel

    ierr=VecRestoreArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_TemplateImage,&p_mt); CHKERRQ(ierr);
    ierr=Rescale(this->m_TemplateImage,0.0,1.0); CHKERRQ(ierr);

    // solve the forward problem using the computed
    // template image and the computed velocity field as input
    ierr=this->SolveForwardProblem(this->m_ReferenceImage); CHKERRQ(ierr);
    ierr=Rescale(this->m_ReferenceImage,0.0,1.0); CHKERRQ(ierr);

    ierr=VecSet(this->m_VelocityField->m_X1,0.0); CHKERRQ(ierr);
    ierr=VecSet(this->m_VelocityField->m_X2,0.0); CHKERRQ(ierr);
    ierr=VecSet(this->m_VelocityField->m_X3,0.0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief copies some input data field to all time points
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CopyToAllTimePoints"
PetscErrorCode OptimalControlRegistrationBase::CopyToAllTimePoints(Vec u, Vec uj)
{
    PetscErrorCode ierr;
    ScalarType *p_u=NULL,*p_uj=NULL;
    IntType nl,nt;
    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();

    // get pointers
    ierr=VecGetArray(u,&p_u); CHKERRQ(ierr);  ///< vec for entire time horizon
    ierr=VecGetArray(uj,&p_uj); CHKERRQ(ierr); ///< vec at single point in time

    // for all time points
    for (IntType j = 0; j <= nt; ++j){

        // copy data to all time points
        try{ std::copy(p_uj,p_uj+nl,p_u+j*nl); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }

    }

    // restore pointers
    ierr=VecRestoreArray(u,&p_u); CHKERRQ(ierr);
    ierr=VecRestoreArray(uj,&p_uj); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief copies some input data field to all time points
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeCFLCondition"
PetscErrorCode OptimalControlRegistrationBase::ComputeCFLCondition()
{
    PetscErrorCode ierr;
    ScalarType vmax,vmaxscaled;
    std::stringstream ss;
    ScalarType hx[3],c;
    IntType ntcfl;
    PetscFunctionBegin;

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }

    if (this->m_WorkVecField1 == NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr=this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i){
        hx[i] = this->m_Opt->GetSpatialStepSize(i);
    }

    ierr=VecAbs(this->m_WorkVecField1->m_X1); CHKERRQ(ierr);
    ierr=VecAbs(this->m_WorkVecField1->m_X2); CHKERRQ(ierr);
    ierr=VecAbs(this->m_WorkVecField1->m_X3); CHKERRQ(ierr);

    // compute max( |v_1| + |v_2| + |v_3| )
    ierr=VecSet(this->m_WorkScaField1,0.0); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkScaField1,1.0,this->m_WorkVecField1->m_X1);
    ierr=VecAXPY(this->m_WorkScaField1,1.0,this->m_WorkVecField1->m_X2);
    ierr=VecAXPY(this->m_WorkScaField1,1.0,this->m_WorkVecField1->m_X3);
    ierr=VecMax(this->m_WorkScaField1,NULL,&vmax); CHKERRQ(ierr);

    // compute max( |v_1|/hx1 + |v_2|/hx2 + |v_3|/hx3 )
    ierr=VecSet(this->m_WorkScaField1,0.0); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkScaField1,1.0/hx[0],this->m_WorkVecField1->m_X1);
    ierr=VecAXPY(this->m_WorkScaField1,1.0/hx[1],this->m_WorkVecField1->m_X2);
    ierr=VecAXPY(this->m_WorkScaField1,1.0/hx[2],this->m_WorkVecField1->m_X3);
    ierr=VecMax(this->m_WorkScaField1,NULL,&vmaxscaled); CHKERRQ(ierr);

    // if we have a zero velocity field, we do not have to worry
    ntcfl = this->m_Opt->GetNumTimePoints();
    if ( vmaxscaled != 0.0 ){

        c  = 1.0;//this->m_opt->cflnum;
        // compute min number of time steps
        ntcfl = static_cast<IntType>(ceil(vmaxscaled/c));

    }

    ss  << "||v||_infty = " << std::scientific << std::fixed << vmax
        << " nt_CFL = " << ntcfl;
    ierr=DbgMsg(ss.str()); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CheckBounds"
PetscErrorCode OptimalControlRegistrationBase::CheckBounds(Vec v, bool& boundreached)
{
    PetscErrorCode ierr;
    ScalarType jmin,jmax,jbound;
    bool minboundreached,maxboundreached;
    std::stringstream ss;
    PetscFunctionBegin;

    if (this->m_VelocityField == NULL){
       try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    minboundreached = false;
    maxboundreached = false;

    // parse input velocity field
    ierr=this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // compute determinant of deformation gradient
    ierr=this->ComputeDetDefGrad(); CHKERRQ(ierr);

    jmin   = this->m_Opt->GetJacMin();
    jmax   = this->m_Opt->GetJacMax();
    jbound = this->m_Opt->GetJacBound();

    // check if jmin < bound and 1/jmax < bound
    minboundreached = jmin <= jbound;
    maxboundreached = 1.0/jmax <= jbound;

    boundreached = (minboundreached || maxboundreached) ? true : false;
    if (boundreached) ierr=WrngMsg("jacobian bound reached"); CHKERRQ(ierr);

    // display what's going on
    if (this->m_Opt->GetVerbosity() > 1){

        if(minboundreached){
            ss << std::scientific
            << "min(det(grad(y))) = "<< jmin << " <= " << jbound;
            ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();
        }

        if(maxboundreached){
            ss << std::scientific
            << "max(det(grad(y))) = "<< jmax << " >= " << 1.0/jbound
            << " ( = 1/bound )";
            ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();
        }

    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDetDefGrad"
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGrad()
{
    PetscErrorCode ierr;
    ScalarType minddg,maxddg,meanddg;
    std::stringstream ss, ssnum;

    PetscFunctionBegin;

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("computing deformation gradient"); CHKERRQ(ierr);
    }

    if (this->m_VelocityField == NULL){
       try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }

    // set initial condition
    ierr=VecSet(this->m_WorkScaField1,1.0); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr=this->IsVelocityZero(); CHKERRQ(ierr);

    if(this->m_VelocityIsZero == false){

        // call the solver
        switch (this->m_Opt->GetPDESolver()){
            case RK2:
            {
                ierr=this->ComputeDetDefGradRK2(); CHKERRQ(ierr);
                break;
            }
            case SL:
            {
                ierr=this->ComputeDetDefGradSL(); CHKERRQ(ierr);
                break;
            }
            default:
            {
                ierr=ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                break;
            }
        }
    }

    ierr=VecMin(this->m_WorkScaField1,NULL,&minddg); CHKERRQ(ierr);
    ierr=VecMax(this->m_WorkScaField1,NULL,&maxddg); CHKERRQ(ierr);
    ierr=VecSum(this->m_WorkScaField1,&meanddg); CHKERRQ(ierr);
    meanddg /= static_cast<ScalarType>(this->m_Opt->GetNGlobal());

    // remember
    this->m_Opt->SetJacMin(minddg);
    this->m_Opt->SetJacMax(maxddg);
    this->m_Opt->SetJacMean(meanddg);


    if (this->m_Opt->GetVerbosity() > 1 || this->m_Opt->MonitorJacobian()){
        ss  << std::scientific << "det(grad(y)) : (min, mean, max)="
            << "(" << minddg << ", " << meanddg << ", " << maxddg<<")";
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str( std::string() ); ss.clear();
    }


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDetDefGradRK2"
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradRK2()
{
    PetscErrorCode ierr;
    ScalarType *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,*p_jbar=NULL,
                *p_gx1=NULL,*p_gx2=NULL,*p_gx3=NULL,*p_divv=NULL,
                *p_jac=NULL, *p_rhs0=NULL;
    ScalarType ht,hthalf;
    IntType nl,nt;
    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;
    double timings[5]={0,0,0,0,0};
    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkVecField1 == NULL){
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField4 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField4); CHKERRQ(ierr);
    }

    // set initial condition
    ierr=VecSet(this->m_WorkScaField1,1.0); CHKERRQ(ierr);

    // get pointers
    ierr=VecGetArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_gx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_gx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_gx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField1,&p_jac); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField2,&p_divv); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField3,&p_jbar); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField4,&p_rhs0); CHKERRQ(ierr);

    // compute div(v)
    accfft_divergence(p_divv,p_vx1,p_vx2,p_vx3,this->m_Opt->GetFFTPlan(),timings);

    // for all time points
    for (IntType j = 0; j <= nt; ++j){

        accfft_grad(p_gx1,p_gx2,p_gx3,p_jac,this->m_Opt->GetFFTPlan(),&XYZ,timings);

#pragma omp parallel
{
#pragma omp  for
        for (IntType i=0; i < nl; ++i){ // for all grid points

            // \bar{j} = -(\vect{v} \cdot \igrad) j + j (\idiv \vect{v})
            p_rhs0[i] = -( p_vx1[i]*p_gx1[i]
                         + p_vx2[i]*p_gx2[i]
                         + p_vx3[i]*p_gx3[i] )
                         + p_jac[i]*p_divv[i];

            p_jbar[i] = p_jac[i] + ht*p_rhs0[i];

        }
} // pragma omp

        accfft_grad(p_gx1,p_gx2,p_gx3,p_jbar,this->m_Opt->GetFFTPlan(),&XYZ,timings);

#pragma omp parallel
{
#pragma omp  for
        for (IntType i=0; i < nl; ++i){ // for all grid points

            // \bar{j} = -(\vect{v} \cdot \igrad) j + j (\idiv \vect{v})
            ScalarType rhs1 = -( p_vx1[i]*p_gx1[i]
                               + p_vx2[i]*p_gx2[i]
                               + p_vx3[i]*p_gx3[i] )
                               + p_jbar[i]*p_divv[i];

            p_jac[i] = p_jac[i] + hthalf*(p_rhs0[i] + rhs1);

        }
} // pragma omp


    }

    ierr=VecRestoreArray(this->m_WorkScaField1,&p_jac); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_divv); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField3,&p_jbar); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField4,&p_rhs0); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_gx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_gx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_gx3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute determinant of deformation gradient
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDetDefGradSL"
PetscErrorCode OptimalControlRegistrationBase::ComputeDetDefGradSL()
{
    PetscErrorCode ierr;
    ScalarType *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_gjx1=NULL,*p_gjx2=NULL,*p_gjx3=NULL,
                *p_divjacv=NULL,*p_divjacvX=NULL,
                *p_cgradvj=NULL,*p_cgradvjX=NULL,
                *p_jvx1=NULL,*p_jvx2=NULL,*p_jvx3=NULL,
                *p_jac=NULL,*p_jacX=NULL,*p_rhs0=NULL;
    ScalarType ht,hthalf;
    IntType nl,nt;
    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;
    double timings[5]={0,0,0,0,0};
    PetscFunctionBegin;

    if (this->m_SL == NULL){
        try{this->m_SL = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL){
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL){
       try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3 == NULL){
       try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField4 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField4); CHKERRQ(ierr);
    }

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    // compute trajectory
    ierr=this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr=this->m_SL->ComputeTrajectory(this->m_WorkVecField1,"state"); CHKERRQ(ierr);


    // get pointers
    ierr=VecGetArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_gjx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_gjx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_gjx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField2->m_X1,&p_jvx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField2->m_X2,&p_jvx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField2->m_X3,&p_jvx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField1,&p_jac); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField2,&p_jacX); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField3,&p_cgradvj); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField4,&p_divjacv); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField3->m_X1,&p_rhs0); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField3->m_X2,&p_cgradvjX); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField3->m_X3,&p_divjacvX); CHKERRQ(ierr);


    for( IntType j = 0; j < nt; ++j ){ // for all time points

        // compute J(X,t^j)
        ierr=this->m_SL->Interpolate(p_jacX,p_jac,"state"); CHKERRQ(ierr);

        // compute grad(jac) for convective derivative
        accfft_grad(p_gjx1,p_gjx2,p_gjx3,p_jac,this->m_Opt->GetFFTPlan(),&XYZ,timings);

#pragma omp parallel
{
        ScalarType jac;
#pragma omp  for
        for (IntType i=0; i < nl; ++i){ // for all grid points

            //(\vect{v} \cdot \igrad) jac
            p_cgradvj[i] =  p_vx1[i]*p_gjx1[i]
                          + p_vx2[i]*p_gjx2[i]
                          + p_vx3[i]*p_gjx3[i];

            jac = p_jac[i];

            p_jvx1[i] = jac*p_vx1[i];
            p_jvx2[i] = jac*p_vx2[i];
            p_jvx3[i] = jac*p_vx3[i];

        }
} // pragma omp

        // compute div(jac v)
        accfft_divergence(p_divjacv,p_jvx1,p_jvx2,p_jvx3,this->m_Opt->GetFFTPlan(),timings);

        // compute J(X,t^j)
        ierr=this->m_SL->Interpolate(p_cgradvjX,p_cgradvj,"state"); CHKERRQ(ierr);

        // compute J(X,t^j)
        ierr=this->m_SL->Interpolate(p_divjacvX,p_divjacv,"state"); CHKERRQ(ierr);


#pragma omp parallel
{
        ScalarType jactilde;
#pragma omp  for
        for (IntType i=0; i < nl; ++i){ // for all grid points

            p_rhs0[i] = (p_divjacvX[i] - p_cgradvjX[i]);

            jactilde = p_jacX[i] + ht*p_rhs0[i];

            p_jac[i] = jactilde;

            p_jvx1[i] = jactilde*p_vx1[i];
            p_jvx2[i] = jactilde*p_vx2[i];
            p_jvx3[i] = jactilde*p_vx3[i];


        }
} // pragma omp

        // compute div(jactilde v)
        accfft_divergence(p_divjacv,p_jvx1,p_jvx2,p_jvx3,this->m_Opt->GetFFTPlan(),timings);

        // compute grad(jactilde) for convective derivative
        accfft_grad(p_gjx1,p_gjx2,p_gjx3,p_jac,this->m_Opt->GetFFTPlan(),&XYZ,timings);

#pragma omp parallel
{
        ScalarType rhs1,cgradvjtilde;
#pragma omp  for
        for (IntType i=0; i < nl; ++i){ // for all grid points

            //(\vect{v} \cdot \igrad) jac
            cgradvjtilde =  p_vx1[i]*p_gjx1[i]
                          + p_vx2[i]*p_gjx2[i]
                          + p_vx3[i]*p_gjx3[i];

            rhs1 = (p_divjacv[i] - cgradvjtilde);

            p_jac[i] = p_jacX[i] + hthalf*(p_rhs0[i] + rhs1);

        }
} // pragma omp

    }

    ierr=VecRestoreArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_gjx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_gjx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_gjx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField2->m_X1,&p_jvx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField2->m_X2,&p_jvx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField2->m_X3,&p_jvx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkScaField1,&p_jac); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_jacX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField3,&p_cgradvj); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField4,&p_divjacv); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField3->m_X1,&p_rhs0); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField3->m_X2,&p_cgradvjX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField3->m_X3,&p_divjacvX); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMap"
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMap()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_VelocityField == NULL){
       try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=this->m_VelocityField->SetValue(0.0); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("computing deformation map"); CHKERRQ(ierr);
    }

    // call the solver
    switch (this->m_Opt->GetPDESolver()){
        case RK2:
        {
            ierr=this->ComputeDeformationMapRK2(); CHKERRQ(ierr);
            break;
        }
        case SL:
        {
            ierr=this->ComputeDeformationMapSL(); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
            break;
        }
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMapRK2"
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapRK2()
{
    //PetscErrorCode ierr;
    PetscFunctionBegin;


    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief compute deformation map
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDeformationMapSL"
PetscErrorCode OptimalControlRegistrationBase::ComputeDeformationMapSL()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(this->m_VelocityField != NULL, "velocity field is null pointer"); CHKERRQ(ierr);

    // allocate vector fields
    if (this->m_WorkVecField1 == NULL){
       try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate semi-lagrangian solver
    if(this->m_SL == NULL){
        try{this->m_SL = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // compute deformation map y using an SL time integrator
    ierr=this->m_WorkVecField1->SetValue(0.0); CHKERRQ(ierr);
    ierr=this->m_SL->ComputeDeformationMap(this->m_WorkVecField1,this->m_VelocityField); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of namespace



#endif// _OPTIMALCONTROLREGISTRATIONBASE_CPP_
