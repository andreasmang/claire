#ifndef _OPTIMALCONTROLREGISTRATION_CPP_
#define _OPTIMALCONTROLREGISTRATION_CPP_

#include <fstream>

// local includes
#include "OptimalControlRegistration.hpp"
#include "RegularizationRegistrationH1.hpp"
#include "RegularizationRegistrationH2.hpp"
#include "RegularizationRegistrationH1SN.hpp"
#include "RegularizationRegistrationH2SN.hpp"




namespace reg
{




/********************************************************************
 * Name: OptimalControlRegistration
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistration"
OptimalControlRegistration::OptimalControlRegistration() : SuperClass()
{
    this->Initialize();
}




/********************************************************************
 * Name: OptimalControlRegistration
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~OptimalControlRegistration"
OptimalControlRegistration::~OptimalControlRegistration(void)
{
    this->ClearMemory();
}




/********************************************************************
 * Name: OptimalControlRegistration
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistration"
OptimalControlRegistration::OptimalControlRegistration(RegOpt* opt) : SuperClass(opt)
{
    this->Initialize();
}




/********************************************************************
 * Name: Initialize
 * Description: init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode OptimalControlRegistration::Initialize(void)
{
    PetscFunctionBegin;

    this->m_TemplateImage = NULL;  ///< reference iamge
    this->m_ReferenceImage = NULL;  ///< template image

    this->m_VelocityField = NULL;  ///< control variable
    this->m_IncVelocityField = NULL;  ///< incremental control variable

    this->m_StateVariable = NULL; ///< state variable
    this->m_AdjointVariable = NULL; ///< adjoin variable
    this->m_IncStateVariable = NULL; ///< incremental state variable
    this->m_IncAdjointVariable = NULL; ///< incremental adjoint variable

    this->m_WorkScaField1 = NULL;
    this->m_WorkScaField2 = NULL;
    this->m_WorkScaField3 = NULL;
    this->m_WorkScaField4 = NULL;

    this->m_WorkVecField1 = NULL;
    this->m_WorkVecField2 = NULL;
    this->m_WorkVecField3 = NULL;
    this->m_WorkVecField4 = NULL;

    this->m_SL = NULL;
    this->m_PCKSP = NULL;

    // we have not done anything
    this->m_InFirstIteration = true;

    this->m_PCMatVec = NULL;
    this->m_Regularization = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode OptimalControlRegistration::ClearMemory(void)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // delete all variables
    if (this->m_StateVariable != NULL){
        ierr=VecDestroy(&this->m_StateVariable); CHKERRQ(ierr);
        this->m_StateVariable=NULL;
    }

    if (this->m_AdjointVariable != NULL){
        ierr=VecDestroy(&this->m_AdjointVariable); CHKERRQ(ierr);
        this->m_AdjointVariable=NULL;
    }

    if (this->m_IncStateVariable != NULL){
        ierr=VecDestroy(&this->m_IncStateVariable); CHKERRQ(ierr);
        this->m_IncStateVariable=NULL;
    }

    if (this->m_IncAdjointVariable != NULL){
        ierr=VecDestroy(&this->m_IncAdjointVariable); CHKERRQ(ierr);
        this->m_IncAdjointVariable=NULL;
    }

    if (this->m_WorkScaField1 != NULL){
        ierr=VecDestroy(&this->m_WorkScaField1); CHKERRQ(ierr);
        this->m_WorkScaField1=NULL;
    }
    if (this->m_WorkScaField2 != NULL){
        ierr=VecDestroy(&this->m_WorkScaField2); CHKERRQ(ierr);
        this->m_WorkScaField2=NULL;
    }
    if (this->m_WorkScaField3 != NULL){
        ierr=VecDestroy(&this->m_WorkScaField3); CHKERRQ(ierr);
        this->m_WorkScaField3=NULL;
    }
    if (this->m_WorkScaField4 != NULL){
        ierr=VecDestroy(&this->m_WorkScaField4); CHKERRQ(ierr);
        this->m_WorkScaField4=NULL;
    }

    if (this->m_VelocityField != NULL){
        delete this->m_VelocityField;
        this->m_VelocityField = NULL;
    }

    if (this->m_IncVelocityField != NULL){
        delete this->m_IncVelocityField;
        this->m_IncVelocityField = NULL;
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

    // delete class for regularization model
    if (this->m_Regularization != NULL){
        delete this->m_Regularization;
        this->m_Regularization = NULL;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: AllocateRegularization
 * Description: allocate regularization model
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "AllocateRegularization"
PetscErrorCode OptimalControlRegistration::AllocateRegularization()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    switch(this->m_Opt->GetRegNorm()){
        case H1:
        {
            try{ this->m_Regularization = new RegularizationRegistrationH1(this->m_Opt); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H2:
        {
            try{ this->m_Regularization = new RegularizationRegistrationH2(this->m_Opt); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H1SN:
        {
            try{ this->m_Regularization = new RegularizationRegistrationH1SN(this->m_Opt); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        case H2SN:
        {
            try{ this->m_Regularization = new RegularizationRegistrationH2SN(this->m_Opt); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
            break;
        }
        default: { ierr=reg::ThrowError("regularization model not defined"); CHKERRQ(ierr); }
    }

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: SolveForwardProblem
 * Description: solve the forward problem (we assume the user has
 * set the template image and the velocity field using the assocated
 * set functions)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveForwardProblem"
PetscErrorCode OptimalControlRegistration::SolveForwardProblem(Vec m)
{
    PetscErrorCode ierr;
    ScalarType *p_m1=NULL,*p_m=NULL;
    IntType nt,nl;
    PetscFunctionBegin;

    ierr=Assert(m!=NULL,"null pointer"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("computing forward solve"); CHKERRQ(ierr);
    }

    // compute solution of state equation
    ierr=this->SolveStateEquation(); CHKERRQ(ierr);

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();

    // copy memory for m_1
    ierr=VecGetArray(m,&p_m1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    try{ std::copy(p_m+nt*nl,p_m+(nt+1)*nl,p_m1); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }
    ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    ierr=VecRestoreArray(m,&p_m1); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateL2Distance
 * Description: evaluate the l2 distance between m_R and m_1
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateL2Distance"
PetscErrorCode OptimalControlRegistration::EvaluateL2Distance(ScalarType* D)
{
    PetscErrorCode ierr;
    ScalarType *p_m1=NULL,*p_m=NULL;
    IntType nt,nl;
    ScalarType hd=0.0,dr;

    PetscFunctionBegin;

    ierr=Assert(this->m_ReferenceImage!=NULL,"reference image is null pointer"); CHKERRQ(ierr);

    if(this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }

    // compute hd
    hd = this->m_Opt->GetLebesqueMeasure();
    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();

    // compute solution of state equation
    ierr=this->SolveStateEquation(); CHKERRQ(ierr);

    // copy memory for m_1
    ierr=VecGetArray(this->m_WorkScaField2,&p_m1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    try{ std::copy(p_m+nt*nl,p_m+(nt+1)*nl,p_m1); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_m1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);

    // compute ||m_R - m_1||_2^2
    ierr=VecWAXPY(this->m_WorkScaField1,-1.0,this->m_WorkScaField2,this->m_ReferenceImage); CHKERRQ(ierr);
    ierr=VecTDot(this->m_WorkScaField1,this->m_WorkScaField1,&dr); CHKERRQ(ierr);

    // objective value
    *D = 0.5*hd*dr;

    PetscFunctionReturn(0);
}







/********************************************************************
 * Name: EvaluateObjective
 * Description: evaluates the objective value
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateObjective"
PetscErrorCode OptimalControlRegistration::EvaluateObjective(ScalarType* J, Vec v)
{
    PetscErrorCode ierr;
    ScalarType D=0.0,R=0.0;
    PetscFunctionBegin;

    // allocate velocity field
    if (this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL){
        ierr=this->AllocateRegularization(); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("evaluating objective functional"); CHKERRQ(ierr);
    }

    ierr=this->m_Opt->StartTimer(OBJEXEC); CHKERRQ(ierr);

    // set components of velocity field
    ierr=this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // evaluate the regularization model
    ierr=this->EvaluateL2Distance(&D); CHKERRQ(ierr);

    // evaluate the regularization model
    ierr=this->m_Regularization->EvaluateFunctional(&R,this->m_VelocityField); CHKERRQ(ierr);

    // add up the contributions
    *J = D + R;

    ierr=this->m_Opt->StopTimer(OBJEXEC); CHKERRQ(ierr);

    this->m_Opt->IncrementCounter(OBJEVAL);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: EvaluateGradient
 * Description: evaluates the reduced gradient of the lagrangian
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EvaluateGradient"
PetscErrorCode OptimalControlRegistration::EvaluateGradient(Vec dvJ, Vec v)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // allocate velocity field
    if (this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate container for reduced gradient
    if (this->m_WorkVecField1 == NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    // allocate container for reduced gradient
    if (this->m_WorkVecField2 == NULL){
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL){
        ierr=this->AllocateRegularization(); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("evaluating reduced gradient"); CHKERRQ(ierr);
    }

    ierr=this->m_Opt->StartTimer(GRADEXEC); CHKERRQ(ierr);

    // parse input arguments
    ierr=this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // compute solution of adjoint equation (i.e., \lambda(x,t))
    ierr=this->SolveAdjointEquation(); CHKERRQ(ierr);

    // compute body force \int_0^1 grad(m)\lambda dt
    // assigned to work vecfield 2
    ierr=this->ComputeBodyForce(); CHKERRQ(ierr);

    // evaluate gradient of regularization model
    ierr=this->m_Regularization->EvaluateGradient(this->m_WorkVecField1,this->m_VelocityField); CHKERRQ(ierr);

    // \vect{g}_v = \beta_v \D{A}[\vect{v}] + \D{K}[\vect{b}]
    ierr=VecAXPY(this->m_WorkVecField1->m_X1,1.0,this->m_WorkVecField2->m_X1); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField1->m_X2,1.0,this->m_WorkVecField2->m_X2); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField1->m_X3,1.0,this->m_WorkVecField2->m_X3); CHKERRQ(ierr);

    // parse to output
    ierr=this->m_WorkVecField1->GetComponents(dvJ); CHKERRQ(ierr);

    ierr=this->m_Opt->StopTimer(GRADEXEC); CHKERRQ(ierr);

    this->m_Opt->IncrementCounter(GRADEVAL);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeBodyForce
 * Description: compute the body force
 * b = \int_0^1 grad(m) \lambda dt
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeBodyForce"
PetscErrorCode OptimalControlRegistration::ComputeBodyForce()
{
    PetscErrorCode ierr;
    IntType nt,nl;
    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;
    ScalarType *p_mj=NULL,*p_m=NULL,*p_l=NULL,*p_l0=NULL,
               *p_gradm1=NULL,*p_gradm2=NULL,*p_gradm3=NULL,
               *p_b1=NULL, *p_b2=NULL, *p_b3=NULL;
    ScalarType ht,scale,hd;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    // check for null pointers
    ierr=Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

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
    if(this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("computing body force"); CHKERRQ(ierr);
    }

    // get problem dimensions and weights
    nl = this->m_Opt->GetNLocal();
    nt = this->m_Opt->GetNumTimePoints();
    ht = this->m_Opt->GetTimeStepSize();
    hd = this->m_Opt->GetLebesqueMeasure();

    ierr=Assert(nt > 0,"number of time points < 0"); CHKERRQ(ierr);
    ierr=Assert(ht > 0,"time step size <= 0"); CHKERRQ(ierr);
    ierr=Assert(hd > 0,"spatial step size <= 0"); CHKERRQ(ierr);

    // init body force for numerical integration
    ierr=VecSet(this->m_WorkVecField2->m_X1,0.0); CHKERRQ(ierr);
    ierr=VecSet(this->m_WorkVecField2->m_X2,0.0); CHKERRQ(ierr);
    ierr=VecSet(this->m_WorkVecField2->m_X3,0.0); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr=this->IsVelocityZero(); CHKERRQ(ierr);
    this->m_VelocityIsZero = false;
    if(this->m_VelocityIsZero){

        // m and \lambda are constant in time
        ierr=VecGetArray(this->m_TemplateImage,&p_mj); CHKERRQ(ierr);

        ierr=VecGetArray(this->m_WorkVecField2->m_X1,&p_gradm1); CHKERRQ(ierr); // for grad(m)_1
        ierr=VecGetArray(this->m_WorkVecField2->m_X2,&p_gradm2); CHKERRQ(ierr); // for grad(m)_2
        ierr=VecGetArray(this->m_WorkVecField2->m_X3,&p_gradm3); CHKERRQ(ierr); // for grad(m)_3

        // computing gradient of m
        accfft_grad(p_gradm1,p_gradm2,p_gradm3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

        ierr=VecRestoreArray(this->m_WorkVecField2->m_X1,&p_gradm1); CHKERRQ(ierr); // for grad(m)_1
        ierr=VecRestoreArray(this->m_WorkVecField2->m_X2,&p_gradm2); CHKERRQ(ierr); // for grad(m)_2
        ierr=VecRestoreArray(this->m_WorkVecField2->m_X3,&p_gradm3); CHKERRQ(ierr); // for grad(m)_3

        ierr=VecRestoreArray(this->m_TemplateImage,&p_mj); CHKERRQ(ierr);

        // compute \igrad(m_0)\lambda_0
        ierr=VecGetArray(this->m_WorkScaField1,&p_l0); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);
        try{ std::copy(p_l,p_l+nl,p_l0); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr=VecRestoreArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkScaField1,&p_l0); CHKERRQ(ierr);

        // scale by hd (ht cancels, as we sum up the same function)
        ierr=VecScale(this->m_WorkScaField1,hd); CHKERRQ(ierr);
        ierr=this->m_WorkVecField2->Scale(this->m_WorkScaField1); CHKERRQ(ierr);

    }
    else{ /// non zero velocity field

        scale = hd*ht;

        // get arrays
        ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr); // state variable for all t^j
        ierr=VecGetArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr); // adjoint variable for all t^j

        ierr=VecGetArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr); // state variable m at time point t^j

        ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_gradm1); CHKERRQ(ierr); // for grad(m)_1
        ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_gradm2); CHKERRQ(ierr); // for grad(m)_2
        ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_gradm3); CHKERRQ(ierr); // for grad(m)_3

        ierr=VecGetArray(this->m_WorkVecField2->m_X1,&p_b1); CHKERRQ(ierr); // for \vect{b}_1
        ierr=VecGetArray(this->m_WorkVecField2->m_X2,&p_b2); CHKERRQ(ierr); // for \vect{b}_2
        ierr=VecGetArray(this->m_WorkVecField2->m_X3,&p_b3); CHKERRQ(ierr); // for \vect{b}_3

        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j){

            IntType k     = j*nl;
            IntType knext = (j+1)*nl;

            // copy memory for m_j
            try{ std::copy(p_m+k,p_m+knext,p_mj); }
            catch(std::exception&){
                ierr=ThrowError("copy failed"); CHKERRQ(ierr);
            }

            // grad(m^j)
            accfft_grad(p_gradm1,p_gradm2,p_gradm3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

            // trapezoidal rule
            if ((j == 0) || (j == nt)){scale*=0.5;};

#pragma omp parallel
{
#pragma omp for

            // for all grid points at tj
            for (IntType i = 0; i < nl; ++i){

                // get \lambda(x_i,t^j)
                ScalarType l = p_l[k+i];

                // \vect{b}_i += h_d*ht*\lambda^j (\grad m^j)_i
                p_b1[i] = p_b1[i] + scale*p_gradm1[i]*l;
                p_b2[i] = p_b2[i] + scale*p_gradm2[i]*l;
                p_b3[i] = p_b3[i] + scale*p_gradm3[i]*l;

            }

} // parallel

            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)){scale*=2.0;};

        }

        // restore arrays
        ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr); // state variable for all t^j
        ierr=VecRestoreArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr); // adjoint variable for all t^j

        ierr=VecRestoreArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr); // state variable m at time point t^j

        ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_gradm1); CHKERRQ(ierr); // for grad(m)_1
        ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_gradm2); CHKERRQ(ierr); // for grad(m)_2
        ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_gradm3); CHKERRQ(ierr); // for grad(m)_3

        ierr=VecRestoreArray(this->m_WorkVecField2->m_X1,&p_b1); CHKERRQ(ierr); // for \vect{b}_1
        ierr=VecRestoreArray(this->m_WorkVecField2->m_X2,&p_b2); CHKERRQ(ierr); // for \vect{b}_2
        ierr=VecRestoreArray(this->m_WorkVecField2->m_X3,&p_b3); CHKERRQ(ierr); // for \vect{b}_3

    } // else zero velocity field

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: HessianMatvec
 * Description: applies the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode OptimalControlRegistration::HessianMatVec(Vec Hvtilde, Vec vtilde)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // check for null pointers
    ierr=Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate container for incremental velocity field
    if (this->m_IncVelocityField == NULL){
        try{this->m_IncVelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate container for hessian matvec
    if (this->m_WorkVecField1 == NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL){
        ierr=this->AllocateRegularization(); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("computing hessian matvec"); CHKERRQ(ierr);
    }

    ierr=this->m_Opt->StartTimer(HMVEXEC); CHKERRQ(ierr);

    // parse input
    ierr=this->m_IncVelocityField->SetComponents(vtilde); CHKERRQ(ierr);

    // compute \tilde{m}(x,t)
    ierr=this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t)
    ierr=this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // compute incremental body force
    ierr=this->ComputeIncBodyForce(); CHKERRQ(ierr);

    // apply 2nd variation of regularization model to
    // incremental control variable: \beta \D{A}[\vect{\tilde{v}}]
    ierr=this->m_Regularization->HessianMatVec(this->m_WorkVecField1,this->m_IncVelocityField); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \beta \D{A}[\vect{\tilde{v}}] + \D{K}[\vect{\tilde{b}}]
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    ierr=VecAXPY(this->m_WorkVecField1->m_X1,1.0,this->m_WorkVecField2->m_X1); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField1->m_X2,1.0,this->m_WorkVecField2->m_X2); CHKERRQ(ierr);
    ierr=VecAXPY(this->m_WorkVecField1->m_X3,1.0,this->m_WorkVecField2->m_X3); CHKERRQ(ierr);

    // pass to output
    ierr=this->m_WorkVecField1->GetComponents(Hvtilde); CHKERRQ(ierr);

    ierr=this->m_Opt->StopTimer(HMVEXEC); CHKERRQ(ierr);

    this->m_Opt->IncrementCounter(HESSMATVEC);

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: ComputeInitialCondition
 * Description: piccard iteration
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeInitialCondition"
PetscErrorCode OptimalControlRegistration::ComputeInitialCondition(Vec m, Vec lambda)
{
    PetscErrorCode ierr;
    IntType nt, nl,ng;
    ScalarType hd;
    PetscFunctionBegin;

    ierr=Assert(this->m_IO!=NULL,"null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetNumTimePoints();
    ng = this->m_Opt->GetNGlobal();
    nl = this->m_Opt->GetNLocal();
    hd = this->m_Opt->GetLebesqueMeasure();

    // allocate container for incremental velocity field
    if (this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL){
        ierr=this->AllocateRegularization(); CHKERRQ(ierr);
    }

    // allocate state and adjoint variables
    if (this->m_StateVariable == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_StateVariable); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_StateVariable,(nt+1)*nl,(nt+1)*ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_StateVariable); CHKERRQ(ierr);
    }

    // allocate state and adjoint variables
    if (this->m_AdjointVariable == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_AdjointVariable); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_AdjointVariable,(nt+1)*nl,(nt+1)*ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_AdjointVariable); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("piccard iteration"); CHKERRQ(ierr);
    }

    // copy input state and adjoint variable to class variables
    ierr=VecCopy(m,this->m_StateVariable); CHKERRQ(ierr);
    ierr=VecCopy(lambda,this->m_AdjointVariable); CHKERRQ(ierr);

    // compute body force (assigned to work vec field 2)
    ierr=this->ComputeBodyForce(); CHKERRQ(ierr);

    // piccard step: solve A[v] = - ht \sum_j \lambda^j grad(m^j)
    ierr=this->m_WorkVecField2->Scale(-1.0/hd); CHKERRQ(ierr);

    ierr=this->m_Regularization->ApplyInverseOperator(this->m_VelocityField,
                                                      this->m_WorkVecField2); CHKERRQ(ierr);

    // reset the adjoint variables
    ierr=VecSet(this->m_StateVariable,0.0); CHKERRQ(ierr);
    ierr=VecSet(this->m_AdjointVariable,0.0); CHKERRQ(ierr);

    ierr=this->m_IO->Write(this->m_VelocityField,
                            "initial-condition-x1.nii.gz",
                            "initial-condition-x2.nii.gz",
                            "initial-condition-x3.nii.gz"); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: PiccardIteration
 * Description: piccard iteration
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PiccardIteration"
PetscErrorCode OptimalControlRegistration::PiccardIteration(Vec v)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // allocate container for incremental velocity field
    if (this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL){
        ierr=this->AllocateRegularization(); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("piccard iteration"); CHKERRQ(ierr);
    }

    // parse input arguments
    ierr=this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // compute solution of adjoint equation (i.e., \lambda(x,t))
    ierr=this->SolveStateEquation(); CHKERRQ(ierr);
    ierr=this->SolveAdjointEquation(); CHKERRQ(ierr);

    // compute body force (assigned to work vec field 2)
    ierr=this->ComputeBodyForce(); CHKERRQ(ierr);

    // piccard step: solve A[v] = - ht \sum_j \lambda^j grad(m^j)
    ierr=this->m_WorkVecField2->Scale(-1.0); CHKERRQ(ierr);
    ierr=this->m_Regularization->ApplyInverseOperator(this->m_VelocityField,this->m_WorkVecField2); CHKERRQ(ierr);

    this->m_Opt->IncrementCounter(GRADEVAL);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: PrecondMatVec
 * Description: applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PrecondMatVec"
PetscErrorCode OptimalControlRegistration::PrecondMatVec(Vec Px, Vec x)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);

    // switch case for choice of preconditioner
    switch(this->m_Opt->GetPrecondMeth()){
        case NOPC:
        {
            ierr=WrngMsg("no preconditioner used"); CHKERRQ(ierr);
            ierr=VecCopy(x,Px); CHKERRQ(ierr);
            break;
        }
        case INVREG:
        {
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

            ierr=this->m_WorkVecField1->SetComponents(x); CHKERRQ(ierr);
            ierr=this->m_Regularization->ApplyInverseOperator(this->m_WorkVecField2,this->m_WorkVecField1); CHKERRQ(ierr);
            ierr=this->m_WorkVecField2->GetComponents(Px); CHKERRQ(ierr);
            break;

        }
        case TWOLEVEL:
        {
            ierr=this->Apply2LevelPrecond(Px,x); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner not defined"); CHKERRQ(ierr);
            break;
        }
    }

    ierr=this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    this->m_Opt->IncrementCounter(PCMATVEC);


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Apply2LevelPrecond
 * Description: applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Apply2LevelPrecond"
PetscErrorCode OptimalControlRegistration::Apply2LevelPrecond(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr=this->Setup2LevelPrecond(); CHKERRQ(ierr);

    ierr=KSPSolve(this->m_PCKSP,x,Px); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: EstimateHessianEigVals
 * Description: estimate hessian eigenvalues
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EstimateHessianEigVals"
PetscErrorCode OptimalControlRegistration::EstimateHessianEigVals()
{
    PetscErrorCode ierr;
    ScalarType emax,emin;
    PetscFunctionBegin;

    emax = 1.0;
    emin = 0.1;

    ierr=KSPChebyshevSetEigenvalues(this->m_PCKSP,emax,emin); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Setup2LevelPrecond
 * Description: applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Setup2LevelPrecond"
PetscErrorCode OptimalControlRegistration::Setup2LevelPrecond()
{
    PetscErrorCode ierr;
    PC pc=NULL;
    ScalarType reltol,abstol,divtol;
    IntType maxit;
    IntType nl,ng;

    PetscFunctionBegin;

    divtol = 1E+06;
    abstol = 1E-16;
    reltol = 1E-16;
    maxit  = 1000;

    nl = this->m_Opt->GetNLocal();
    ng = this->m_Opt->GetNGlobal();

    if (this->m_PCKSP == NULL){
        ierr=KSPCreate(PETSC_COMM_WORLD,&this->m_PCKSP); CHKERRQ(ierr);
    }


    switch (this->m_Opt->GetPCSolverType()){
        case PCCHEB:
        {
            // chebyshev iteration
            ierr=KSPSetType(this->m_PCKSP,KSPCHEBYSHEV); CHKERRQ(ierr);
            maxit  = 10;
            break;
        }
        case PCPCG:
        {
            // preconditioned conjugate gradient
            ierr=KSPSetType(this->m_PCKSP,KSPCG); CHKERRQ(ierr);
            reltol = this->GetKSPTolerance();
            reltol *= 1E-1;
            break;
        }
        case PCFCG:
        {
            // flexible conjugate gradient
            ierr=KSPSetType(this->m_PCKSP,KSPFCG); CHKERRQ(ierr);
            maxit  = 10;
            break;
        }
        case PCGMRES:
        {
            // GMRES
            ierr=KSPSetType(this->m_PCKSP,KSPGMRES); CHKERRQ(ierr);
            reltol = this->GetKSPTolerance();
            reltol *= 1E-1;
            break;
        }
        case PCFGMRES:
        {
            // flexible GMRES
            ierr=KSPSetType(this->m_PCKSP,KSPFGMRES); CHKERRQ(ierr);
            maxit  = 10;
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }
    reltol = std::max(reltol,1E-16); // make sure tolerance is non-zero
    reltol = std::min(reltol,0.25); // make sure tolerance smaller than 0.25


    ierr=KSPSetTolerances(this->m_PCKSP,reltol,abstol,divtol,maxit); CHKERRQ(ierr);

    //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
    //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
    //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
    ierr=KSPSetNormType(this->m_PCKSP,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    ierr=KSPSetInitialGuessNonzero(this->m_PCKSP,PETSC_TRUE); CHKERRQ(ierr);

    // set up matvec for preconditioner
    if (this->m_PCMatVec != NULL){
        ierr=MatDestroy(&this->m_PCMatVec); CHKERRQ(ierr);
        this->m_PCMatVec = NULL;
     }

    ierr=MatCreateShell(PETSC_COMM_WORLD,3*nl,3*nl,3*ng,3*ng,this,&this->m_PCMatVec); CHKERRQ(ierr);
    ierr=MatShellSetOperation(this->m_PCMatVec,MATOP_MULT,(void(*)(void))TwoLevelPCMatVec); CHKERRQ(ierr);

    // set operator
    ierr=KSPSetOperators(this->m_PCKSP,this->m_PCMatVec,this->m_PCMatVec);CHKERRQ(ierr);
    ierr=KSPMonitorSet(this->m_PCKSP,PrecondMonitor,this,PETSC_NULL); CHKERRQ(ierr);
    // remove preconditioner
    ierr=KSPGetPC(this->m_PCKSP,&pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCNONE); CHKERRQ(ierr); ///< set no preconditioner

    // finish
    ierr=KSPSetUp(this->m_PCKSP); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: TwoLevelPrecondMatVec
 * Description: applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "TwoLevelPrecondMatVec"
PetscErrorCode OptimalControlRegistration::TwoLevelPrecondMatVec(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr=VecCopy(x,Px); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}





/********************************************************************
 * Name: ComputeIncBodyForce
 * Description: compute the incremental body force
 * \tilde{\vect{b}} = \int_0^1 \igrad\tilde{m}\lambda
 *                           + \igrad m\tilde{\lambda} dt
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeIncBodyForce"
PetscErrorCode OptimalControlRegistration::ComputeIncBodyForce()
{
    PetscErrorCode ierr;
    IntType nt=0;
    IntType nl=0;

    ScalarType *p_m=NULL,*p_mj=NULL,*p_mt=NULL,*p_mtj=NULL,
                *p_l=NULL,*p_lt=NULL,*p_bt1=NULL,*p_bt2=NULL,*p_bt3=NULL,
                *p_gradm1=NULL,*p_gradm2=NULL,*p_gradm3=NULL,
                *p_gradmt1=NULL,*p_gradmt2=NULL,*p_gradmt3=NULL;

    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;
    ScalarType ht=0.0,scale=0.0,hd=0.0;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("computing incremental body force"); CHKERRQ(ierr);
    }

    ierr=Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
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
    if (this->m_WorkVecField3 == NULL){
        try{this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hd = this->m_Opt->GetLebesqueMeasure();

    ierr=Assert(nt > 0,"number of time points <= 0"); CHKERRQ(ierr);
    ierr=Assert(ht > 0,"time step size <= 0"); CHKERRQ(ierr);
    ierr=Assert(hd > 0,"spatial step size <= 0"); CHKERRQ(ierr);
    scale = ht*hd;

    ierr=this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField2->m_X1,&p_bt1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField2->m_X2,&p_bt2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField2->m_X3,&p_bt3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_gradm1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_gradm2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_gradm3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr); // state variable for all t^j
    ierr=VecGetArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr); // temporary vector
    ierr=VecGetArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr); // incremental adjoint variable for all t^j

    // TODO: add case for zero velocity field
    if (this->m_Opt->GetOptMeth() == FULLNEWTON){

        ierr=Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
        ierr=Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);

        if (this->m_WorkScaField2 == NULL){
            ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
        }

        ierr=VecGetArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr); // adjoint variable for all t^j
        ierr=VecGetArray(this->m_IncStateVariable,&p_mt); CHKERRQ(ierr); // incremental state variable for all t^j
        ierr=VecGetArray(this->m_WorkScaField2,&p_mtj); CHKERRQ(ierr); // incremental adjoint variable for all t^j

        ierr=VecGetArray(this->m_WorkVecField3->m_X1,&p_gradmt1); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField3->m_X2,&p_gradmt2); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField3->m_X3,&p_gradmt3); CHKERRQ(ierr);

        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j){

            // copy memory for m_j
            try{ std::copy(p_m+j*nl,p_m+(j+1)*nl,p_mj); }
            catch(std::exception&){
                ierr=ThrowError("copy failed"); CHKERRQ(ierr);
            }

            // copy memory for \tilde{m}_j
            try{ std::copy(p_mt+j*nl,p_mt+(j+1)*nl,p_mtj); }
            catch(std::exception&){
                ierr=ThrowError("copy failed"); CHKERRQ(ierr);
            }

            // computing gradient of m
            accfft_grad(p_gradm1,p_gradm2,p_gradm3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

            // computing gradient of \tilde{m}
            accfft_grad(p_gradmt1,p_gradmt2,p_gradmt3,p_mtj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)){scale*=0.5;};

#pragma omp parallel
{
#pragma omp for
            // compute \vect{\tilde{b}}^k_i
            // += h_d*ht*(\tilde{\lambda}^j (\grad m^j)^k
            //    + \lambda^j (\grad \tilde{m}^j)^k)_i
            for (IntType i=0; i < nl; ++i){ // for all grid points

                // get \lambda(x_i,t^j) and \tilde{\lambda}(x_i,t^j)
                ScalarType lj  = p_l[j*nl+i];
                ScalarType ltj = p_lt[j*nl+i];

                p_bt1[i] = p_bt1[i] + scale*(p_gradm1[i]*ltj + p_gradmt1[i]*lj);
                p_bt2[i] = p_bt2[i] + scale*(p_gradm2[i]*ltj + p_gradmt2[i]*lj);
                p_bt3[i] = p_bt3[i] + scale*(p_gradm3[i]*ltj + p_gradmt3[i]*lj);

            } // for all grid points

} // pragma omp parallel

            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)){scale*=2.0;};

        } // for all time points

        ierr=VecRestoreArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr); // adjoint variable for all t^j
        ierr=VecRestoreArray(this->m_IncStateVariable,&p_mt); CHKERRQ(ierr); // incremental state variable for all t^j
        ierr=VecRestoreArray(this->m_WorkScaField2,&p_mtj); CHKERRQ(ierr); // incremental adjoint variable for all t^j

        ierr=VecRestoreArray(this->m_WorkVecField3->m_X1,&p_gradmt1); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField3->m_X2,&p_gradmt2); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField3->m_X3,&p_gradmt3); CHKERRQ(ierr);

    }// full newton
    else if (this->m_Opt->GetOptMeth() == GAUSSNEWTON){

        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j){ // for all time points

            IntType k     = j*nl;
            IntType knext = (j+1)*nl;

            // copy memory for m^j
            try{ std::copy(p_m+k,p_m+knext,p_mj); }
            catch(std::exception&){
                ierr=ThrowError("copy failed"); CHKERRQ(ierr);
            }

            // compute gradient of m^j
            accfft_grad(p_gradm1,p_gradm2,p_gradm3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)){scale*=0.5;};

#pragma omp parallel
{
#pragma omp for
            // compute \vect{\tilde{b}}^k_i += h_d*ht*(\tilde{\lambda}^j (\grad m^j)^k
            for (IntType i=0; i < nl; ++i){ // for all grid points

                // get \tilde{\lambda}(x_i,t^j)
                ScalarType ltj = p_lt[k+i];

                p_bt1[i] = p_bt1[i] + scale*p_gradm1[i]*ltj;
                p_bt2[i] = p_bt2[i] + scale*p_gradm2[i]*ltj;
                p_bt3[i] = p_bt3[i] + scale*p_gradm3[i]*ltj;

            } // for all grid points
} // pragma omp parallel

            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)){scale*=2.0;};

        } // for all time points

    } // gauss newton
    else{ ierr=ThrowError("hessian approximation not implemented"); CHKERRQ(ierr); }

    // restore all arrays
    ierr=VecRestoreArray(this->m_WorkVecField2->m_X1,&p_bt1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField2->m_X2,&p_bt2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField2->m_X3,&p_bt3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_gradm1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_gradm2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_gradm3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr); // state variable for all t^j
    ierr=VecRestoreArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr); // temporary vector
    ierr=VecRestoreArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr); // incremental adjoint variable for all t^j

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveStateEquation
 * Description: solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveStateEquation"
PetscErrorCode OptimalControlRegistration::SolveStateEquation(void)
{
    PetscErrorCode ierr=0;
    IntType nl,ng,nt;
    ScalarType *p_m=NULL,*p_m0=NULL;
    std::stringstream ss;
    PetscFunctionBegin;

    ierr=Assert(this->m_VelocityField!=NULL,"velocity is null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_TemplateImage!=NULL,"template image is null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetNumTimePoints();
    ng = this->m_Opt->GetNGlobal();
    nl = this->m_Opt->GetNLocal();

    ierr=Assert(nt > 0, "number of time points <= 0"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2){
        ss << "solving state equation (nt="<<nt<<")";
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
    }

    if (this->m_Opt->MonitorCFLCondition()){
        ierr=this->ComputeCFLCondition(); CHKERRQ(ierr);
    }

    ierr=this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // allocate state and adjoint variables
    if (this->m_StateVariable == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_StateVariable); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_StateVariable,(nt+1)*nl,(nt+1)*ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_StateVariable); CHKERRQ(ierr);
        ierr=VecSet(this->m_StateVariable,0.0); CHKERRQ(ierr);
    }

    // check if velocity field is zero
    ierr=this->IsVelocityZero(); CHKERRQ(ierr);
    if(this->m_VelocityIsZero){

        // we copy m_0 to all t for v=0
        ierr=this->CopyToAllTimePoints(this->m_StateVariable,this->m_TemplateImage); CHKERRQ(ierr);

    }
    else{

        // copy initial condition m_0 = m_T
        ierr=VecGetArray(this->m_TemplateImage,&p_m0); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
        try{ std::copy(p_m0,p_m0+nl,p_m); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr=VecRestoreArray(this->m_TemplateImage,&p_m0); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);

        // call the solver
        switch (this->m_Opt->GetPDESolver()){
            case RK2:
            {
                ierr=this->SolveStateEquationRK2(); CHKERRQ(ierr);
                break;
            }
            case SL:
            {
                ierr=this->SolveStateEquationSL(); CHKERRQ(ierr);
                break;
            }
            default:
            {
                ierr=ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                break;
            }
        }

    } // velocity field is zero

    ierr=this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveStateEquationRK2
 * Description: solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveStateEquationRK2"
PetscErrorCode OptimalControlRegistration::SolveStateEquationRK2(void)
{
    PetscErrorCode ierr;
    IntType nl,nt;
    ScalarType *p_mj=NULL,*p_m=NULL,*p_mbar=NULL,*p_rhs0=NULL,
                *p_gmx1=NULL,*p_gmx2=NULL,*p_gmx3=NULL,
                *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL;
    ScalarType ht=0.0,hthalf=0.0;
    double ffttimers[5]={0,0,0,0,0};
    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;

    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if(this->m_WorkVecField1 == NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
    }

    // copy initial condition to buffer
    ierr=VecGetArray(this->m_WorkScaField1,&p_mbar); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_gmx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_gmx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_gmx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField2,&p_mj); CHKERRQ(ierr);

    // copy memory (m_0 to m_j)
    try{ std::copy(p_m,p_m+nl,p_mj); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j){

        // compute gradient of m_j
        accfft_grad(p_gmx1,p_gmx2,p_gmx3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
        for (IntType i=0; i < nl; ++i){

             p_rhs0[i] = -p_gmx1[i]*p_vx1[i]
                         -p_gmx2[i]*p_vx2[i]
                         -p_gmx3[i]*p_vx3[i];

             // compute intermediate result
             p_mbar[i] = p_mj[i] + ht*p_rhs0[i];

        }
} // pragma omp parallel

        // compute gradient of \bar{m}
        accfft_grad(p_gmx1,p_gmx2,p_gmx3,p_mbar,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
        for (IntType i=0; i < nl; ++i){

            ScalarType rhs1 = -p_gmx1[i]*p_vx1[i]
                              -p_gmx2[i]*p_vx2[i]
                              -p_gmx3[i]*p_vx3[i];

            // we have overwritten m_j with intermediate result
            // m_j = m_{j-1} + 0.5*ht*(RHS0 + RHS1)
            p_mj[i] = p_mj[i] + hthalf*(p_rhs0[i] + rhs1);
        }
} // parallel

        // copy to buffer
        try{ std::copy(p_mj,p_mj+nl,p_m+(j+1)*nl); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }

    } // for all time points

    // copy initial condition to buffer
    ierr=VecGetArray(this->m_WorkScaField2,&p_mj); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField1,&p_mbar); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_gmx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_gmx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_gmx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveStateEquationSL
 * Description: solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveStateEquationSL"
PetscErrorCode OptimalControlRegistration::SolveStateEquationSL(void)
{
    PetscErrorCode ierr;
    IntType nl,nt;
    ScalarType *p_m=NULL,*p_mj=NULL,*p_mjX=NULL;
    std::stringstream ss;
    std::string filename;
    PetscFunctionBegin;

    if(this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }
    if (this->m_WorkVecField1 == NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_SL == NULL){
        try{this->m_SL = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();

    // compute trajectory
    ierr=this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr=this->m_SL->ComputeTrajectory(this->m_WorkVecField1,"state"); CHKERRQ(ierr);


    // copy memory (m_0 to m_j)
    ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);
    try{ std::copy(p_m,p_m+nl,p_mj); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }

    ierr=VecGetArray(this->m_WorkScaField2,&p_mjX); CHKERRQ(ierr);

    for( IntType j = 0; j < nt; ++j ){ // for all time points

        // compute m(X,t^{j+1}) (interpolate state variable)
        ierr=this->m_SL->Interpolate(p_mjX,p_mj,"state"); CHKERRQ(ierr);

        // store m(X,t^{j+1})
        try{ std::copy(p_mjX,p_mjX+nl,p_mj); }
        catch(std::exception&){
            ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
        }

        // store time history
        try{ std::copy(p_mj,p_mj+nl,p_m+(j+1)*nl); }
        catch(std::exception&){
            ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
        }

    }// for all time points

    ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_mjX); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveAdjointEquation
 * Description: solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveAdjointEquation"
PetscErrorCode OptimalControlRegistration::SolveAdjointEquation(void)
{
    PetscErrorCode ierr=0;
    IntType nl,ng,nt;
    ScalarType *p_m=NULL,*p_m1=NULL,*p_l=NULL,*p_l0=NULL;
    std::stringstream ss;

    PetscFunctionBegin;

    ierr=Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetNumTimePoints();
    ng = this->m_Opt->GetNGlobal();
    nl = this->m_Opt->GetNLocal();
    ierr=Assert(nt > 0, "number of time points < 0"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2){
        ss << "solving adjoint equation (nt="<<nt<<")";
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
    }

    ierr=this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // allocate variables
    if (this->m_AdjointVariable == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_AdjointVariable); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_AdjointVariable,(nt+1)*nl,(nt+1)*ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_AdjointVariable); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField1==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField2==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }


    // copy memory for m_1
    ierr=VecGetArray(this->m_WorkScaField1,&p_m1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    try{ std::copy(p_m+nt*nl,p_m+(nt+1)*nl,p_m1); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }
    ierr=VecRestoreArray(this->m_WorkScaField1,&p_m1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);


    // compute terminal condition \lambda_1 = -(m_1 - m_R) = m_R - m_1
    ierr=VecWAXPY(this->m_WorkScaField2,-1.0,this->m_WorkScaField1,this->m_ReferenceImage); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr=this->IsVelocityZero(); CHKERRQ(ierr);
    if(this->m_VelocityIsZero){

        // we copy \lambda_1 to all t for v=0
        ierr=this->CopyToAllTimePoints(this->m_AdjointVariable,this->m_WorkScaField2); CHKERRQ(ierr);

    }
    else{

        // copy final condition \lambda = (m_R - m) at t=1
        ierr=VecGetArray(this->m_WorkScaField2,&p_l0); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);

        try{ std::copy(p_l0,p_l0+nl,p_l+nt*nl); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }

        ierr=VecRestoreArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkScaField2,&p_l0); CHKERRQ(ierr);

        // call the solver
        switch (this->m_Opt->GetPDESolver()){
            case RK2:
            {
                ierr=this->SolveAdjointEquationRK2(); CHKERRQ(ierr);
                break;
            }
            case SL:
            {
                ierr=this->SolveAdjointEquationSL(); CHKERRQ(ierr);
                break;
            }
            default:
            {
                ierr=ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                break;
            }
        }

    }

    ierr=this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveAdjointEquationRK2
 * Description: solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveAdjointEquationRK2"
PetscErrorCode OptimalControlRegistration::SolveAdjointEquationRK2(void)
{
    PetscErrorCode ierr;
    IntType nl;
    IntType nt;
    ScalarType *p_l=NULL,*p_rhs0=NULL,*p_rhs1=NULL,*p_lj=NULL,
                *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_ljvx1=NULL,*p_ljvx2=NULL,*p_ljvx3=NULL;
    ScalarType hthalf,ht;
    double ffttimers[5]={0,0,0,0,0};
    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    ierr=Assert(nt > 0, "number of time points < 0"); CHKERRQ(ierr);
    ierr=Assert(ht > 0, "number of time points < 0"); CHKERRQ(ierr);

    if(this->m_WorkScaField1==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField3==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField4==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField4); CHKERRQ(ierr);
    }
    if(this->m_WorkVecField1==NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr=VecGetArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField4,&p_rhs1); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_ljvx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_ljvx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_ljvx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField1,&p_lj); CHKERRQ(ierr);
    try{ std::copy(p_l+nt*nl,p_l+(nt+1)*nl,p_lj); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j){ // for all time points

#pragma omp parallel
{
#pragma omp for
        // scale \vect{v} by \lambda
        for(IntType i=0; i < nl; ++i){ // for all grid points

            ScalarType lambda = p_lj[i];

            p_ljvx1[i] = lambda*p_vx1[i];
            p_ljvx2[i] = lambda*p_vx2[i];
            p_ljvx3[i] = lambda*p_vx3[i];

        }// for all grid points
} // pragma omp parallel

        // compute \idiv(\lambda\vect{v})
        accfft_divergence(p_rhs0,p_ljvx1,p_ljvx2,p_ljvx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
        for(IntType i=0; i < nl; ++i){ // for all grid points

            // compute \bar{\lambda} = \lambda_j + ht*\idiv(\lambda\vect{v})
            ScalarType lbar = p_lj[i] + ht*p_rhs0[i];

            // scale \vect{v} by \bar{\lambda}
            p_ljvx1[i] = p_vx1[i]*lbar;
            p_ljvx2[i] = p_vx2[i]*lbar;
            p_ljvx3[i] = p_vx3[i]*lbar;

        }
} // pragma omp parallel

        // compute \idiv(\bar{\lambda}\vect{v})
        accfft_divergence(p_rhs1,p_ljvx1,p_ljvx2,p_ljvx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
        for(IntType i=0; i < nl; ++i){ // for all grid points
            p_lj[i] = p_lj[i] + hthalf*(p_rhs0[i]+p_rhs1[i]);
        }
} // pragma omp parallel

        try{ std::copy(p_lj,p_lj+nl,p_l+(nt-(j+1))*nl); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }

    } // for all time points

    ierr=VecRestoreArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField4,&p_rhs1); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_ljvx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_ljvx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_ljvx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveAdjointEquationSL
 * Description: solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveAdjointEquationSL"
PetscErrorCode OptimalControlRegistration::SolveAdjointEquationSL()
{
    PetscErrorCode ierr;
    double ffttimers[5]={0,0,0,0,0};
    ScalarType *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_divv=NULL,*p_divvX=NULL,
                *p_l=NULL,*p_lj=NULL,*p_ljX=NULL;
    ScalarType ht;
    IntType nl;
    IntType nt;

    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();

    ierr=Assert(this->m_VelocityField!=NULL,"null pointer"); CHKERRQ(ierr);

    if(this->m_WorkScaField1==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField2==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField3==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
    }
    if(this->m_WorkScaField4==NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField4); CHKERRQ(ierr);
    }
    if(this->m_WorkVecField1==NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_SL == NULL){
        try{this->m_SL = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }


    // scale v by -1
    ierr=this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr=this->m_WorkVecField1->Scale(-1.0); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_vx3); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField3,&p_divv); CHKERRQ(ierr);

    // compute \idiv(\tilde{\lambda}\vect{v})
    accfft_divergence(p_divv,p_vx1,p_vx2,p_vx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
    this->m_Opt->IncrementCounter(FFT,4);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_vx3); CHKERRQ(ierr);

    // compute trajectory
    ierr=this->m_SL->ComputeTrajectory(this->m_WorkVecField1,"adjoint"); CHKERRQ(ierr);

    // evaluate div(v) at X
    ierr=VecGetArray(this->m_WorkScaField4,&p_divvX); CHKERRQ(ierr);
    ierr=this->m_SL->Interpolate(p_divvX,p_divv,"adjoint"); CHKERRQ(ierr);

    // copy final condition
    ierr=VecGetArray(this->m_WorkScaField1,&p_lj); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);
    try{ std::copy(p_l+nt*nl,p_l+(nt+1)*nl,p_lj); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }

    ierr=VecGetArray(this->m_WorkScaField2,&p_ljX); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j){

        // compute lambda(t^j,X)
        ierr=this->m_SL->Interpolate(p_ljX,p_lj,"adjoint"); CHKERRQ(ierr);

#pragma omp parallel
{
        IntType knext = (nt-(j+1))*nl;
#pragma omp for
        for (IntType i = 0; i < nl; ++i){

            ScalarType ljX  = p_ljX[i];

            ScalarType rhs0 = -ljX*p_divvX[i];
            ScalarType rhs1 = -(ljX + ht*rhs0)*p_divv[i];

            // compute \lambda(x,t^{j+1})
            p_lj[i] = ljX + 0.5*ht*(rhs0 + rhs1);

            // TODO: check if it gets faster if I do a memcopy outside
            p_l[knext + i] = p_lj[i];

        }
} // pragma omp parallel

    }

    ierr=VecRestoreArray(this->m_WorkScaField4,&p_divvX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField1,&p_lj); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_ljX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField3,&p_divv); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);

}





/********************************************************************
 * Name: SolveIncStateEquation
 * Description: solve the incremental state equation
 * \p_t \tilde{m} + \igrad m \cdot \vect{\tilde{v}}
 *                + \igrad \tilde{m} \cdot \vect{v} = 0
 * subject to \tilde{m}_0 = 0
 * solved forward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveIncStateEquation"
PetscErrorCode OptimalControlRegistration::SolveIncStateEquation(void)
{
    PetscErrorCode ierr=0;
    IntType nl,ng,nt;
    std::stringstream ss;
    PetscFunctionBegin;

    ierr=Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetNumTimePoints();
    ng = this->m_Opt->GetNGlobal();
    nl = this->m_Opt->GetNLocal();
    ierr=Assert(nt > 0, "number of time points < 0"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2){
        ss << "solving incremental state equation (nt="<<nt<<")";
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
    }

    ierr=this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // allocate variables
    if (this->m_IncStateVariable == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_IncStateVariable); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_IncStateVariable,(nt+1)*nl,(nt+1)*ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_IncStateVariable); CHKERRQ(ierr);
    }

    // set initial value
    ierr=VecSet(this->m_IncStateVariable,0.0); CHKERRQ(ierr);

    // call the solver
    switch (this->m_Opt->GetPDESolver()){
        case RK2:
        {
            ierr=this->SolveIncStateEquationRK2(); CHKERRQ(ierr);
            break;
        }
        case SL:
        {
            ierr=this->SolveIncStateEquationSL(); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
            break;
        }
    }

    ierr=this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveIncStateEquationRK2
 * Description: solve the incremental state equation
 * \p_t \tilde{m} + \igrad m \cdot \vect{\tilde{v}}
 *                + \igrad \tilde{m} \cdot \vect{v} = 0
 * subject to \tilde{m}_0 = 0
 * solved forward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveIncStateEquationRK2"
PetscErrorCode OptimalControlRegistration::SolveIncStateEquationRK2(void)
{
    PetscErrorCode ierr=0;
    IntType nl=0;
    IntType nt=0;
    ScalarType *p_m=NULL,*p_mj=NULL,*p_mt=NULL,*p_mtj=NULL,*p_mtbar=NULL,
                *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_gmx1=NULL,*p_gmx2=NULL,*p_gmx3=NULL,
                *p_gmtx1=NULL,*p_gmtx2=NULL,*p_gmtx3=NULL,
                *p_vtx1=NULL,*p_vtx2=NULL,*p_vtx3=NULL,*p_rhs0=NULL;
    ScalarType ht=0.0,hthalf=0.0;
    double ffttimers[5]={0,0,0,0,0};
    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;
    PetscFunctionBegin;


    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    ierr=Assert(nt > 0, "number of time points < 0"); CHKERRQ(ierr);
    ierr=Assert(ht > 0, "time step size < 0"); CHKERRQ(ierr);

    // allocate variables
    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
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

    ierr=VecGetArray(this->m_IncStateVariable,&p_mt); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_gmx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_gmx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_gmx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_IncVelocityField->m_X1,&p_vtx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_IncVelocityField->m_X2,&p_vtx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_IncVelocityField->m_X3,&p_vtx3); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr=this->IsVelocityZero(); CHKERRQ(ierr);
    if(this->m_VelocityIsZero){

        // compute gradient of m_1 (m is constant)
        ierr=VecGetArray(this->m_TemplateImage,&p_mj); CHKERRQ(ierr);
        accfft_grad(p_gmx1,p_gmx2,p_gmx3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);
        ierr=VecRestoreArray(this->m_TemplateImage,&p_mj); CHKERRQ(ierr);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j){

            IntType k     = j*nl;
            IntType knext = (j+1)*nl;

            // the right hand side remains constant;
            // we can reduce the 2 RK2 steps to a single one
#pragma omp parallel
{
#pragma omp for
            for(IntType i=0; i < nl; ++i){
                 p_mt[knext+i] = p_mt[k+i] - ht*(p_gmx1[i]*p_vtx1[i]
                                                +p_gmx2[i]*p_vtx2[i]
                                                +p_gmx3[i]*p_vtx3[i]);
            }
} // pragma omp parallel

        } // for all time points

    }
    else{ // velocity field is non-zero

        ierr=VecGetArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkScaField2,&p_mtbar); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkScaField4,&p_mtj); CHKERRQ(ierr);

        ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);

        ierr=VecGetArray(this->m_WorkVecField2->m_X1,&p_gmtx1); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField2->m_X2,&p_gmtx2); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_WorkVecField2->m_X3,&p_gmtx3); CHKERRQ(ierr);

        ierr=VecGetArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

        // copy initial condition
        try{ std::copy(p_mt,p_mt+nl,p_mtj); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j){

            IntType k     = j*nl;
            IntType knext = (j+1)*nl;

            // get m_j
            try{ std::copy(p_m+k,p_m+knext,p_mj); }
            catch(std::exception&){
                ierr=ThrowError("copy failed"); CHKERRQ(ierr);
            }

            // compute gradient of m_j
            accfft_grad(p_gmx1,p_gmx2,p_gmx3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

            // compute gradient of \tilde{m}_j
            accfft_grad(p_gmtx1,p_gmtx2,p_gmtx3,p_mtj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
            for(IntType i=0; i < nl; ++i){

                 p_rhs0[i] = -p_gmtx1[i]*p_vx1[i] - p_gmtx2[i]*p_vx2[i] - p_gmtx3[i]*p_vx3[i]
                             -p_gmx1[i]*p_vtx1[i] - p_gmx2[i]*p_vtx2[i] - p_gmx3[i]*p_vtx3[i];

                 // compute intermediate result
                 p_mtbar[i] = p_mtj[i] + ht*p_rhs0[i];

            }
} // pragma omp parallel

            // get m_{j+1}
            try{ std::copy(p_m+knext,p_m+(j+2)*nl,p_mj); }
            catch(std::exception&){
                ierr=ThrowError("copy failed"); CHKERRQ(ierr);
            }

            // compute gradient of m_{j+1}
            accfft_grad(p_gmx1,p_gmx2,p_gmx3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

            // compute gradient of \tilde{m}_j
            accfft_grad(p_gmtx1,p_gmtx2,p_gmtx3,p_mtbar,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
            for(IntType i=0; i < nl; ++i){

                // evaluate right hand side
                ScalarType rhs1 = -p_gmtx1[i]*p_vx1[i] - p_gmtx2[i]*p_vx2[i] - p_gmtx3[i]*p_vx3[i]
                                  -p_gmx1[i]*p_vtx1[i] - p_gmx2[i]*p_vtx2[i] - p_gmx3[i]*p_vtx3[i];

                // compute intermediate result
                p_mtj[i] = p_mtj[i] + hthalf*(rhs1 + p_rhs0[i]);

            }
} // pragma omp parallel

            // store time history
            try{ std::copy(p_mtj,p_mtj+nl,p_mt+knext); }
            catch(std::exception&){
                ierr=ThrowError("copy failed"); CHKERRQ(ierr);
            }

        } // for all time points

        // copy initial condition to buffer
        ierr=VecRestoreArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkScaField2,&p_mtbar); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkScaField4,&p_mtj); CHKERRQ(ierr);

        ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);

        ierr=VecRestoreArray(this->m_WorkVecField2->m_X1,&p_gmtx1); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField2->m_X2,&p_gmtx2); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_WorkVecField2->m_X3,&p_gmtx3); CHKERRQ(ierr);

        ierr=VecRestoreArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);


    }// velzero


    ierr=VecRestoreArray(this->m_IncVelocityField->m_X1,&p_vtx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_IncVelocityField->m_X2,&p_vtx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_IncVelocityField->m_X3,&p_vtx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_gmx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_gmx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_gmx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_IncStateVariable,&p_mt); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveIncStateEquationSL
 * Description: solve the incremental state equation
 * \p_t \tilde{m} + \igrad m \cdot \vect{\tilde{v}}
 *                + \igrad \tilde{m} \cdot \vect{v} = 0
 * subject to \tilde{m}_0 = 0
 * solved forward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveIncStateEquationSL"
PetscErrorCode OptimalControlRegistration::SolveIncStateEquationSL(void)
{
    PetscErrorCode ierr=0;
    IntType nl=0;
    IntType nt=0;
    std::bitset<3> XYZ; XYZ[0]=1;XYZ[1]=1;XYZ[2]=1;
    ScalarType ht,hthalf;
    double ffttimers[5]={0,0,0,0,0};
    ScalarType *p_gmjx1=NULL,*p_gmjx2=NULL,*p_gmjx3=NULL,
                *p_gmjnextx1=NULL,*p_gmjnextx2=NULL,*p_gmjnextx3=NULL,
                *p_gmjXx1=NULL,*p_gmjXx2=NULL,*p_gmjXx3=NULL,
                *p_mtilde=NULL,*p_mtj=NULL,*p_mtjX=NULL,*p_m=NULL,*p_mj=NULL;
    const ScalarType *p_vtildeXx1=NULL,*p_vtildeXx2=NULL,*p_vtildeXx3=NULL,
                *p_vtildex1=NULL,*p_vtildex2=NULL,*p_vtildex3=NULL;

    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    ierr=Assert(this->m_StateVariable!=NULL, "state variable is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_IncStateVariable!=NULL, "incremental state variable is null"); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
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
    if (this->m_WorkVecField4 == NULL){
        try{this->m_WorkVecField4 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_SL == NULL){
        try{this->m_SL = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=this->m_SL->ComputeTrajectory(this->m_VelocityField,"state"); CHKERRQ(ierr);
    }

    ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_IncStateVariable,&p_mtilde); CHKERRQ(ierr);

    // compute \tilde{\vect{v}}(X)
    ierr=this->m_SL->Interpolate(this->m_WorkVecField2,this->m_IncVelocityField,"state"); CHKERRQ(ierr);

    // get \tilde{\vect{v}}(X)
    ierr=VecGetArrayRead(this->m_WorkVecField2->m_X1,&p_vtildeXx1); CHKERRQ(ierr);
    ierr=VecGetArrayRead(this->m_WorkVecField2->m_X2,&p_vtildeXx2); CHKERRQ(ierr);
    ierr=VecGetArrayRead(this->m_WorkVecField2->m_X3,&p_vtildeXx3); CHKERRQ(ierr);

    // get \tilde{\vect{v}}
    ierr=VecGetArrayRead(this->m_IncVelocityField->m_X1,&p_vtildex1); CHKERRQ(ierr);
    ierr=VecGetArrayRead(this->m_IncVelocityField->m_X2,&p_vtildex2); CHKERRQ(ierr);
    ierr=VecGetArrayRead(this->m_IncVelocityField->m_X3,&p_vtildex3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField3->m_X1,&p_gmjx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField3->m_X2,&p_gmjx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField3->m_X3,&p_gmjx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField4->m_X1,&p_gmjXx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField4->m_X2,&p_gmjXx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField4->m_X3,&p_gmjXx3); CHKERRQ(ierr);

    // copy gradient for m_0 = m_T
    ierr=VecGetArray(this->m_TemplateImage,&p_mj); CHKERRQ(ierr);
    accfft_grad(p_gmjx1,p_gmjx2,p_gmjx3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
    this->m_Opt->IncrementCounter(FFT,4);
    ierr=VecRestoreArray(this->m_TemplateImage,&p_mj); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_gmjnextx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_gmjnextx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_gmjnextx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField2,&p_mtjX); CHKERRQ(ierr);

    // copy initial condition
    ierr=VecGetArray(this->m_WorkScaField3,&p_mtj); CHKERRQ(ierr);
    try{ std::copy(p_mtilde,p_mtilde+nl,p_mtj); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }

    for (IntType j = 0; j < nt; ++j){

        // get m_{j+1}
        try{ std::copy(p_m+(j+1)*nl,p_m+(j+2)*nl,p_mj); }
        catch(std::exception&){
            ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
        }

        // compute gradient
        accfft_grad(p_gmjnextx1,p_gmjnextx2,p_gmjnextx3,p_mj,this->m_Opt->m_MiscOpt->plan,&XYZ,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

        // interpolate gradient
        ierr=this->m_SL->Interpolate(p_gmjXx1,p_gmjXx2,p_gmjXx3,p_gmjx1,p_gmjx2,p_gmjx3,"state");

        // interpolate gradient
        ierr=this->m_SL->Interpolate(p_mtjX,p_mtj,"state"); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i){

            p_mtj[i] = p_mtjX[i] - hthalf*( p_gmjXx1[i]*p_vtildeXx1[i]
                                          + p_gmjXx2[i]*p_vtildeXx2[i]
                                          + p_gmjXx3[i]*p_vtildeXx3[i]
                                          + p_gmjnextx1[i]*p_vtildex1[i]
                                          + p_gmjnextx2[i]*p_vtildex2[i]
                                          + p_gmjnextx3[i]*p_vtildex3[i] );
        }

} // parallel

        try{ std::copy(p_mtj,p_mtj+nl,p_mtilde+(j+1)*nl); }
        catch(std::exception&){
            ierr=ThrowError("copying failed"); CHKERRQ(ierr);
        }

        try{ std::copy(p_gmjnextx1,p_gmjnextx1+nl,p_gmjx1); }
        catch(std::exception&){
            ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
        }
        try{ std::copy(p_gmjnextx2,p_gmjnextx2+nl,p_gmjx2); }
        catch(std::exception&){
            ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
        }
        try{ std::copy(p_gmjnextx3,p_gmjnextx3+nl,p_gmjx3); }
        catch(std::exception&){
            ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
        }
    }

    ierr=VecRestoreArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_mtjX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField3,&p_mtj); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_gmjnextx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_gmjnextx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_gmjnextx3); CHKERRQ(ierr);

    ierr=VecRestoreArrayRead(this->m_WorkVecField2->m_X1,&p_vtildeXx1); CHKERRQ(ierr);
    ierr=VecRestoreArrayRead(this->m_WorkVecField2->m_X2,&p_vtildeXx2); CHKERRQ(ierr);
    ierr=VecRestoreArrayRead(this->m_WorkVecField2->m_X3,&p_vtildeXx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField3->m_X1,&p_gmjx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField3->m_X2,&p_gmjx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField3->m_X3,&p_gmjx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField4->m_X1,&p_gmjXx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField4->m_X2,&p_gmjXx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField4->m_X3,&p_gmjXx3); CHKERRQ(ierr);

    ierr=VecRestoreArrayRead(this->m_IncVelocityField->m_X1,&p_vtildex1); CHKERRQ(ierr);
    ierr=VecRestoreArrayRead(this->m_IncVelocityField->m_X2,&p_vtildex2); CHKERRQ(ierr);
    ierr=VecRestoreArrayRead(this->m_IncVelocityField->m_X3,&p_vtildex3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_IncStateVariable,&p_mtilde); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: SolveIncAdjointEquation
 * Description: solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveIncAdjointEquation"
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquation(void)
{
    PetscErrorCode ierr=0;
    ScalarType *p_lt=NULL,*p_mt=NULL,*p_ltj=NULL;
    IntType nl,ng,k,nt;
    std::stringstream ss;

    PetscFunctionBegin;

    ierr=Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }

    nt = this->m_Opt->GetNumTimePoints();
    ng = this->m_Opt->GetNGlobal();
    nl = this->m_Opt->GetNLocal();
    ierr=Assert(nt > 0, "number of time points not set correctly"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2){
        ss << "solving incremental adjoint equation (nt="<<nt<<")";
        ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
    }

    ierr=this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // allocate state and adjoint variables
    if (this->m_IncAdjointVariable == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_IncAdjointVariable); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_IncAdjointVariable,(nt+1)*nl,(nt+1)*ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_IncAdjointVariable); CHKERRQ(ierr);
    }

    // set terminal condition \tilde{\lambda}_1 = -\tilde{m}_1
    ierr=VecGetArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_IncStateVariable,&p_mt); CHKERRQ(ierr);
    k = nt*nl;
    for (IntType i = 0; i < nl; ++i){
        p_lt[k+i] = -p_mt[k+i];
    }
    ierr=VecRestoreArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_IncStateVariable,&p_mt); CHKERRQ(ierr);


    // check if velocity field is zero
    if(this->m_Opt->GetOptMeth() == GAUSSNEWTON){

        ierr=this->IsVelocityZero(); CHKERRQ(ierr);
        if(this->m_VelocityIsZero){

            // set terminal condition \tilde{\lambda}_1 = -\tilde{m}_1
            ierr=VecGetArray(this->m_WorkScaField1,&p_ltj); CHKERRQ(ierr);
            ierr=VecGetArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr);
            try{ std::copy(p_lt+nt*nl,p_lt+(nt+1)*nl,p_ltj); }
            catch(std::exception&){
                ierr=ThrowError("copy failed"); CHKERRQ(ierr);
            }
            ierr=VecRestoreArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr);
            ierr=VecRestoreArray(this->m_WorkScaField1,&p_ltj); CHKERRQ(ierr);

            // we copy \tilde{\lambda}_1 to all t for v=0
            ierr=this->CopyToAllTimePoints(this->m_IncAdjointVariable,this->m_WorkScaField1); CHKERRQ(ierr);

        }
        else{

            // call the solver
            switch (this->m_Opt->GetPDESolver()){
                case RK2:
                {
                    ierr=this->SolveIncAdjointEquationGNRK2(); CHKERRQ(ierr);
                    break;
                }
                case SL:
                {
                    ierr=this->SolveIncAdjointEquationGNSL(); CHKERRQ(ierr);
                    break;
                }
                default:
                {
                    ierr=ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                    break;
                }
            }
        } // zero velocity field

    } // gauss newton
    else if(this->m_Opt->GetOptMeth() == FULLNEWTON){

        // call the solver
        switch (this->m_Opt->GetPDESolver()){
            case RK2:
            {
                ierr=ThrowError("not tested"); CHKERRQ(ierr);
                ierr=this->SolveIncAdjointEquationFNRK2(); CHKERRQ(ierr);
                break;
            }
            case SL:
            {
                ierr=ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                //ierr=this->SolveIncAdjointEquationRNRK2(); CHKERRQ(ierr);
                break;
            }
            default:
            {
                ierr=ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                break;
            }
        }

    }
    else{ierr=ThrowError("update method not defined"); CHKERRQ(ierr);}

    ierr=this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: SolveIncAdjointEquationGNRK2
 * Description: solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveIncAdjointEquationGNRK2"
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquationGNRK2(void)
{
    PetscErrorCode ierr;
    IntType nl=0;
    IntType nt=0;
    ScalarType *p_lt=NULL,*p_ltj=NULL,*p_rhs0=NULL,*p_rhs1=NULL,
                *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_ltjvx1=NULL,*p_ltjvx2=NULL,*p_ltjvx3=NULL;
    double ffttimers[5]={0,0,0,0,0};
    ScalarType ht,hthalf;
    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
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

    ierr=VecGetArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField4,&p_rhs1); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_ltjvx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_ltjvx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_ltjvx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField1,&p_ltj); CHKERRQ(ierr);

    // get terminal condition \tilde{\lambda}_1 = -\tilde{m}_1
    try{ std::copy(p_lt+nt*nl,p_lt+(nt+1)*nl,p_ltj); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }


    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j){ // for all time points

#pragma omp parallel
{
#pragma omp for
        // scale \vect{v} by \lambda
        for(IntType i=0; i < nl; ++i){ // for all grid points

            ScalarType lt = p_ltj[i];

            p_ltjvx1[i] = p_vx1[i]*lt;
            p_ltjvx2[i] = p_vx2[i]*lt;
            p_ltjvx3[i] = p_vx3[i]*lt;

        }// for all grid points
} // pragma omp parallel

        // compute \idiv(\tilde{\lambda}\vect{v})
        accfft_divergence(p_rhs0,p_ltjvx1,p_ltjvx2,p_ltjvx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
        for(IntType i=0; i < nl; ++i){ // for all grid points

            // compute \bar{\tilde{\lambda}} = \tilde{\lambda}^j + ht*\idiv(\tilde{\lambda}^j\vect{v})
            ScalarType ltbar = p_ltj[i] + ht*p_rhs0[i];

            // scale \vect{v} by \bar{\lambda}
            p_ltjvx1[i] = p_vx1[i]*ltbar;
            p_ltjvx2[i] = p_vx2[i]*ltbar;
            p_ltjvx3[i] = p_vx3[i]*ltbar;

        }
} // pragma omp parallel

        // compute \idiv(\bar{\lambda}\vect{v})
        accfft_divergence(p_rhs1,p_ltjvx1,p_ltjvx2,p_ltjvx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
        for(IntType i=0; i < nl; ++i){ // for all grid points
            p_ltj[i] = p_ltj[i] + hthalf*(p_rhs0[i]+p_rhs1[i]);
        }
} // pragma omp parallel

        try{ std::copy(p_ltj,p_ltj+nl,p_lt+(nt-(j+1))*nl); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }

    } // for all time points

    ierr=VecRestoreArray(this->m_WorkScaField1,&p_ltj); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField4,&p_rhs1); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_ltjvx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_ltjvx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_ltjvx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);

}



/********************************************************************
 * Name: SolveIncAdjointEquationFNRK2
 * Description: solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveIncAdjointEquationFNRK2"
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquationFNRK2(void)
{
    PetscErrorCode ierr;
    IntType nl=0;
    IntType nt=0;
    ScalarType *p_l=NULL,*p_lj=NULL,*p_lt=NULL,
                *p_rhs0=NULL,*p_rhs1=NULL,
                *p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL,
                *p_vtx1=NULL,*p_vtx2=NULL,*p_vtx3=NULL,
                *p_ltjvx1=NULL,*p_ltjvx2=NULL,*p_ltjvx3=NULL;
    double ffttimers[5]={0,0,0,0,0};
    ScalarType ht,hthalf;

    PetscFunctionBegin;

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkScaField3 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
    }

    ierr=VecGetArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr);


    ierr=VecGetArray(this->m_IncVelocityField->m_X1,&p_vtx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_IncVelocityField->m_X2,&p_vtx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_IncVelocityField->m_X3,&p_vtx3); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_ltjvx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_ltjvx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_ltjvx3); CHKERRQ(ierr);

    ierr=this->IsVelocityZero(); CHKERRQ(ierr);
    if(this->m_VelocityIsZero){

        // lambda is constant
        ierr=VecGetArray(this->m_AdjointVariable,&p_lj); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
        // lambda  and v are constant in time
        // scale \vect{v} by \lambda
        for(IntType i=0; i < nl; ++i){ // for all grid points

            ScalarType lambda = p_lj[i];

            p_ltjvx1[i] = p_vtx1[i]*lambda;
            p_ltjvx2[i] = p_vtx2[i]*lambda;
            p_ltjvx3[i] = p_vtx3[i]*lambda;

        }// for all grid points
} // pragma omp parallel

        ierr=VecRestoreArray(this->m_AdjointVariable,&p_lj); CHKERRQ(ierr);

        // compute \idiv(\tilde{\lambda}\vect{v})
        accfft_divergence(p_rhs0,p_ltjvx1,p_ltjvx2,p_ltjvx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
        this->m_Opt->IncrementCounter(FFT,4);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j){ // for all time points

            IntType k = (nt-j)*nl;
            IntType knext = (nt-(j+1))*nl;
#pragma omp parallel
{
#pragma omp for
            for(IntType i = 0; i < nl; ++i){ // for all grid points
                p_lt[knext+i] = p_lt[k+i] + ht*p_rhs0[i];
            }
} // pragma omp parallel

        } // for all time points

    } // velocity is zero
    else{

        ierr=VecGetArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);

        ierr=VecGetArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

        ierr=VecGetArray(this->m_WorkScaField4,&p_rhs1); CHKERRQ(ierr);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j){ // for all time points

            IntType k     = (nt-j)*nl;
            IntType knext = (nt-(j+1))*nl;
#pragma omp parallel
{
#pragma omp for
            for(IntType i=0; i < nl; ++i){ // for all grid points

                ScalarType l  = p_l[k+i];
                ScalarType lt = p_lt[k+i];

                p_ltjvx1[i] = p_vx1[i]*lt + p_vtx1[i]*l;
                p_ltjvx2[i] = p_vx2[i]*lt + p_vtx2[i]*l;
                p_ltjvx3[i] = p_vx3[i]*lt + p_vtx3[i]*l;

            }// for all grid points
} // pragma omp parallel

            // compute \idiv(\tilde{\lambda}\vect{v})
            accfft_divergence(p_rhs0,p_ltjvx1,p_ltjvx2,p_ltjvx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
            for(IntType i=0; i < nl; ++i){ // for all grid points

                // \bar{\lambda} = \tilde{\lambda}^j + ht*\idiv(\lambda^j\vect{v})
                ScalarType ltbar = p_lt[k+i] + ht*p_rhs0[i];
                ScalarType l     = p_l[knext+i];

                // v \bar{\lambda} + \vect{\tilde{v}}\lambda^{j+1}
                p_ltjvx1[i] = p_vx1[i]*ltbar + p_vtx1[i]*l;
                p_ltjvx2[i] = p_vx2[i]*ltbar + p_vtx2[i]*l;
                p_ltjvx3[i] = p_vx3[i]*ltbar + p_vtx3[i]*l;

            }
} // pragma omp parallel

            // compute \idiv(\bar{\lambda}\vect{v})
            accfft_divergence(p_rhs1,p_ltjvx1,p_ltjvx2,p_ltjvx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
            this->m_Opt->IncrementCounter(FFT,4);

#pragma omp parallel
{
#pragma omp for
            for(IntType i=0; i < nl; ++i){ // for all grid points
                p_lt[knext+i] = p_lt[k+i] + hthalf*(p_rhs0[i]+p_rhs1[i]);
            }
} // pragma omp parallel


        } // for all time points

        ierr=VecRestoreArray(this->m_WorkScaField4,&p_rhs1); CHKERRQ(ierr);

        ierr=VecRestoreArray(this->m_VelocityField->m_X1,&p_vx1); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_VelocityField->m_X2,&p_vx2); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_VelocityField->m_X3,&p_vx3); CHKERRQ(ierr);

        ierr=VecRestoreArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);

    } // velzero

    ierr=VecRestoreArray(this->m_IncAdjointVariable,&p_lt); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_IncVelocityField->m_X1,&p_vtx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_IncVelocityField->m_X2,&p_vtx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_IncVelocityField->m_X3,&p_vtx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_ltjvx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_ltjvx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_ltjvx3); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_WorkScaField3,&p_rhs0); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: SolveIncAdjointEquationGNSL
 * Description: solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveIncAdjointEquationGNSL"
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquationGNSL(void)
{
    PetscErrorCode ierr;
    IntType nl=0;
    IntType nt=0;
    double ffttimers[5]={0,0,0,0,0};
    ScalarType *p_ltilde=NULL,*p_ltildej=NULL,*p_ltildejX=NULL,
                *p_divv=NULL,*p_divvX=NULL,*p_vx1=NULL,*p_vx2=NULL,*p_vx3=NULL;
    ScalarType ht=0.0,hthalf=0.0;

    PetscFunctionBegin;

    if (this->m_WorkVecField1 == NULL){
        try{this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField4 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField4); CHKERRQ(ierr);
    }
    if (this->m_SL == NULL){
        try{this->m_SL = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=this->m_SL->ComputeTrajectory(this->m_VelocityField,"adjoint"); CHKERRQ(ierr);
    }

    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    // remember time history (i.e. copy final condition
    // $\tilde{\lambda}_1 = -\tilde{m}_1$ into buffer for $\tilde{\lambda}
    ierr=VecGetArray(this->m_IncAdjointVariable,&p_ltilde); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField1,&p_ltildej); CHKERRQ(ierr);
    try{ std::copy(p_ltilde+nt*nl,p_ltilde+(nt+1)*nl,p_ltildej); }
    catch(std::exception&){
        ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
    }

    // set v to -v
    ierr=this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr=this->m_WorkVecField1->Scale(-1.0); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField3,&p_divv); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkVecField1->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkVecField1->m_X3,&p_vx3); CHKERRQ(ierr);

    // compute div(v)
    accfft_divergence(p_divv,p_vx1,p_vx2,p_vx3,this->m_Opt->m_MiscOpt->plan,ffttimers);
    this->m_Opt->IncrementCounter(FFT,4);

    ierr=VecRestoreArray(this->m_WorkVecField1->m_X1,&p_vx1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X2,&p_vx2); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkVecField1->m_X3,&p_vx3); CHKERRQ(ierr);

    // evaluate div(v) on characteristic
    ierr=VecGetArray(this->m_WorkScaField4,&p_divvX); CHKERRQ(ierr);

    ierr=this->m_SL->Interpolate(p_divvX,p_divv,"adjoint"); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField2,&p_ltildejX); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j){

        ierr=this->m_SL->Interpolate(p_ltildejX,p_ltildej,"adjoint"); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i){

            // get \tilde{\lambda}(X(t^j))
            ScalarType ltildejX = p_ltildejX[i];

            // scale v(X) by \tilde{\lambda}(X(t^j))
            ScalarType rhs0 = -ltildejX*p_divvX[i];

            // scale v by \lamba{\lambda}
            ScalarType rhs1 = -(ltildejX + ht*rhs0)*p_divv[i];

            // final rk2 step
            p_ltildej[i] = ltildejX + hthalf*(rhs0 + rhs1);

        }
} // pragma omp parallel

        // store time history (necessary for optimization)
        try{ std::copy(p_ltildej,p_ltildej+nl,p_ltilde+(nt-(j+1))*nl); }
        catch(std::exception&){
            ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
        }

    } // for all time points

    ierr=VecRestoreArray(this->m_WorkScaField1,&p_ltildej); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_ltildejX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField4,&p_divvX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField3,&p_divv); CHKERRQ(ierr);

    ierr=VecRestoreArray(this->m_IncAdjointVariable,&p_ltilde); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: SolveIncAdjointEquationFNSL
 * Description: solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveIncAdjointEquationFNSL"
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquationFNSL(void)
{
    PetscErrorCode ierr;
    //IntType nl=0;
    //IntType nt=0;
    PetscFunctionBegin;

/*
    if (this->m_SL == NULL){
        try{this->m_SL = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=this->m_SL->ComputeTrajectory(this->m_VelocityField,"state"); CHKERRQ(ierr);
    }
*/
    ierr=ThrowError("not implemented"); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: FinalizeIteration
 * Description: finalize the current iteration
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "FinalizeIteration"
PetscErrorCode OptimalControlRegistration::FinalizeIteration(Vec v)
{

    PetscErrorCode ierr;
    int rank;
    IntType nl,nt;
    std::string filename,fnx1,fnx2,fnx3;
    std::stringstream ss;
    std::ofstream logwriter;
    ScalarType *p_m1=NULL,*p_m=NULL;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }

    // compute hd
    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();

    // if not yet allocted, do so
    if (this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    ierr=this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // store iterates
    if ( this->m_Opt->StoreIterates() ){

        // copy memory for m_1
        ierr=VecGetArray(this->m_WorkScaField1,&p_m1); CHKERRQ(ierr);
        ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
        try{ std::copy(p_m+nt*nl,p_m+(nt+1)*nl,p_m1); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr=VecRestoreArray(this->m_WorkScaField1,&p_m1); CHKERRQ(ierr);
        ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);

        ss  << "deformed-template-image-i="
            << std::setw(3) << std::setfill('0')
            << this->m_NumOuterIter  << ".nii.gz";
        ierr=this->m_IO->Write(this->m_WorkScaField1,ss.str()); CHKERRQ(ierr);
        ss.str("");

        // construct file names for velocity field components
        ss  << "velocity-field-i="
            << std::setw(3) << std::setfill('0')
            << this->m_NumOuterIter  << "-x1.nii.gz";
        fnx1 = ss.str();
        ss.str("");

        ss  << "velocity-field-i="
            << std::setw(3) << std::setfill('0')
            << this->m_NumOuterIter  << "-x2.nii.gz";
        fnx2 = ss.str();
        ss.str("");

        ss  << "velocity-field-i="
            << std::setw(3) << std::setfill('0')
            << this->m_NumOuterIter  << "-x3.nii.gz";
        fnx3 = ss.str();
        ss.str("");

        // velocity field out
        ierr=this->m_IO->Write(this->m_VelocityField,fnx1,fnx2,fnx3); CHKERRQ(ierr);

    } // store iterates

    // compute determinant of deformation gradient and write it to file
    if ( this->m_Opt->MonitorJacobian() ){

        ierr=this->ComputeDetDefGrad(); CHKERRQ(ierr);

        if (rank == 0){

            filename = this->m_Opt->GetXFolder()
                        + "registration-performance-jacobians.log";
            // create output file or append to output file
            logwriter.open(filename.c_str(), std::ofstream::out | std::ofstream::app );
            ierr=Assert(logwriter.is_open(),"could not open file for writing"); CHKERRQ(ierr);
            ss  << std::scientific
                <<  "iter " << std::setw(3)
                << std::right << this->m_NumOuterIter << "    " << std::left
                << std::setw(20) << this->m_Opt->GetJacMin() << " "
                << std::setw(20) << this->m_Opt->GetJacMean() <<" "
                << std::setw(20) << this->m_Opt->GetJacMax();
            logwriter << ss.str() << std::endl;

        }

    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Finalize
 * Description: finalize the registration
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Finalize"
PetscErrorCode OptimalControlRegistration::Finalize(Vec v)
{
    PetscErrorCode ierr;
    std::string filename,fn,ext;
    IntType nt;
    int rank,nproc,nstr,nnum;
    IntType nl;
    std::ofstream logwriter;
    std::stringstream ss, ssnum;
    ScalarType mRmT_2,mRmT_infty,mRm1_2,
                mRm1_infty,mR_2,mR_infty,
                drrel_infty,drrel_2;
    ScalarType *p_m1=NULL, *p_m=NULL;

    PetscFunctionBegin;

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nproc);

    if (this->m_Opt->GetVerbosity() >= 2){
        ierr=DbgMsg("finalizing registration"); CHKERRQ(ierr);
    }

    // if not yet allocted, do so
    if (this->m_VelocityField == NULL){
        try{this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField3); CHKERRQ(ierr);
    }

    // process timers
    ierr=this->m_Opt->ProcessTimers(); CHKERRQ(ierr);

    // write log file
    if (this->m_Opt->LoggingEnabled()){
        ierr=this->m_Opt->WriteLogFile(); CHKERRQ(ierr);
    }

    // set components of velocity field
    ierr=this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // parse extension
    ext = ".nii.gz";//this->m_Opt->GetXExtension();

    ierr=VecWAXPY(this->m_WorkScaField1,-1.0,this->m_TemplateImage,this->m_ReferenceImage); CHKERRQ(ierr);
    ierr=VecNorm(this->m_WorkScaField1,NORM_2,&mRmT_2); CHKERRQ(ierr);
    ierr=VecNorm(this->m_WorkScaField1,NORM_INFINITY,&mRmT_infty); CHKERRQ(ierr);

    if( this->m_Opt->WriteImages() == true ){

        ierr=Assert(this->m_IO != NULL,"null pointer"); CHKERRQ(ierr);

        //  write reference and template image
        ierr=this->m_IO->Write(this->m_ReferenceImage,"reference-image"+ext); CHKERRQ(ierr);
        ierr=this->m_IO->Write(this->m_TemplateImage,"template-image"+ext); CHKERRQ(ierr);

        //ierr=VecAbs(this->m_WorkScaField1); CHKERRQ(ierr);
        ierr=this->m_IO->Write(this->m_WorkScaField1,"residual-before"+ext); CHKERRQ(ierr);

    }

    // deformed template out
    // compute solution of state equation
    ierr=this->SolveStateEquation(); CHKERRQ(ierr);

    // compute hd
    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();

    // copy memory for m_1
    ierr=VecGetArray(this->m_WorkScaField1,&p_m1); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);
    try{ std::copy(p_m+nt*nl,p_m+(nt+1)*nl,p_m1); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }
    ierr=VecRestoreArray(this->m_WorkScaField1,&p_m1); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_StateVariable,&p_m); CHKERRQ(ierr);

    // ||m_R - m_1 ||
    ierr=VecWAXPY(this->m_WorkScaField2,-1.0,this->m_WorkScaField1,this->m_ReferenceImage); CHKERRQ(ierr);
    ierr=VecNorm(this->m_WorkScaField2,NORM_2,&mRm1_2); CHKERRQ(ierr);
    ierr=VecNorm(this->m_WorkScaField2,NORM_INFINITY,&mRm1_infty); CHKERRQ(ierr);

    if( this->m_Opt->WriteImages() == true ){

        ierr=Rescale(this->m_WorkScaField1,0,1); CHKERRQ(ierr);

        ierr=this->m_IO->Write(this->m_WorkScaField1,"deformed-template-image"+ext); CHKERRQ(ierr);
        ierr=this->m_IO->Write(this->m_WorkScaField2,"residual-after"+ext); CHKERRQ(ierr);

        // velocity field out
        ierr=this->m_IO->Write(this->m_VelocityField,"velocity-field-x1"+ext,
                                                     "velocity-field-x2"+ext,
                                                     "velocity-field-x3"+ext); CHKERRQ(ierr);

        ierr=VecPointwiseMult(this->m_WorkScaField1,this->m_VelocityField->m_X1,this->m_VelocityField->m_X1); CHKERRQ(ierr);
        ierr=VecPointwiseMult(this->m_WorkScaField2,this->m_VelocityField->m_X2,this->m_VelocityField->m_X2); CHKERRQ(ierr);
        ierr=VecPointwiseMult(this->m_WorkScaField3,this->m_VelocityField->m_X3,this->m_VelocityField->m_X3); CHKERRQ(ierr);

        ierr=VecAXPY(this->m_WorkScaField1,1.0,this->m_WorkScaField2); CHKERRQ(ierr);
        ierr=VecAXPY(this->m_WorkScaField1,1.0,this->m_WorkScaField3); CHKERRQ(ierr);
        ierr=VecSqrtAbs(this->m_WorkScaField1); CHKERRQ(ierr);
        ierr=this->m_IO->Write(this->m_WorkScaField1,"velocity-field-2norm"+ext); CHKERRQ(ierr);

        // determinant of deformation gradient out
        ierr=this->ComputeDetDefGrad(); CHKERRQ(ierr);
        ierr=Assert( this->m_WorkScaField1 != NULL, "null pointer"); CHKERRQ(ierr);
        ierr=this->m_IO->Write(this->m_WorkScaField1,"det-deformation-grad"+ext); CHKERRQ(ierr);

        ierr=this->ComputeDeformationMap(); CHKERRQ(ierr);
        ierr=Assert( this->m_WorkVecField1 != NULL, "null pointer"); CHKERRQ(ierr);
        ierr=this->m_IO->Write(this->m_WorkVecField1,"deformation-map-x1"+ext,
                                                     "deformation-map-x2"+ext,
                                                     "deformation-map-x3"+ext); CHKERRQ(ierr);

    }

    ierr=VecNorm(this->m_ReferenceImage,NORM_2,&mR_2); CHKERRQ(ierr);
    ierr=VecNorm(this->m_ReferenceImage,NORM_INFINITY,&mR_infty); CHKERRQ(ierr);

    mRmT_infty = (mRmT_infty > 0.0) ? mRmT_infty : 1.0;
    mRmT_2     = (mRmT_2     > 0.0) ? mRmT_2     : 1.0;

    drrel_infty=mRm1_infty/mRmT_infty;
    drrel_2=mRm1_2/mRmT_2;

    if (this->m_Opt->LoggingEnabled()){

        if (rank == 0){

            nnum = 20; nstr = 20;
            filename = this->m_Opt->GetXFolder() + "registration-performance-residuals";
            fn = filename + ".log";

            // create output file
            logwriter.open(fn.c_str());
            ierr=Assert(logwriter.is_open(),"could not open file for writing"); CHKERRQ(ierr);

            ss  << std::scientific << std::left
                << std::setw(nstr) << "||mR-mT||_2" << std::right
                << std::setw(nnum) << mRmT_2;
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(nstr) << "||mR-mT||_infty" << std::right
                << std::setw(nnum) << mRmT_infty;
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(nstr) << "||mR-m1||_2" << std::right
                << std::setw(nnum) << mRm1_2;
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(nstr) << "||mR-m1||_infty" << std::right
                << std::setw(nnum) << mRm1_infty;
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(nstr) << "||mR-m1||_2,rel" << std::right
                << std::setw(nnum) << drrel_2;
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(nstr) << "||mR-m1||_infty,rel" << std::right
                << std::setw(nnum) << drrel_infty;
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

        }
    }

    PetscFunctionReturn(0);

}



} // end of name space

#endif // _OPTIMALCONTROLREGISTRATION_CPP_

