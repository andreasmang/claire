#ifndef _OPTIMALCONTROLREGISTRATIONIC_CPP_
#define _OPTIMALCONTROLREGISTRATIONIC_CPP_

#include <math.h>
#include "OptimalControlRegistrationIC.hpp"




namespace reg
{




/********************************************************************
 * Name: OptimalControlRegistrationIC
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistrationIC"
OptimalControlRegistrationIC::OptimalControlRegistrationIC() : SuperClass()
{
    this->Initialize();
}




/********************************************************************
 * Name: OptimalControlRegistrationIC
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~OptimalControlRegistrationIC"
OptimalControlRegistrationIC::~OptimalControlRegistrationIC(void)
{
    this->ClearMemory();
}




/********************************************************************
 * Name: OptimalControlRegistrationIC
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistrationIC"
OptimalControlRegistrationIC::OptimalControlRegistrationIC(RegOpt* opt) : SuperClass(opt)
{
    this->Initialize();
}




/********************************************************************
 * Name: Initialize
 * Description: init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode OptimalControlRegistrationIC::Initialize(void)
{
    PetscFunctionBegin;

    this->m_x1hat=NULL;
    this->m_x2hat=NULL;
    this->m_x3hat=NULL;

    this->m_Kx1hat=NULL;
    this->m_Kx2hat=NULL;
    this->m_Kx3hat=NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode OptimalControlRegistrationIC::ClearMemory(void)
{
//    PetscErrorCode ierr;
    PetscFunctionBegin;

    if(this->m_x1hat!=NULL){
        accfft_free(this->m_x1hat);
        this->m_x1hat = NULL;
    }
    if(this->m_x2hat!=NULL){
        accfft_free(this->m_x2hat);
        this->m_x2hat = NULL;
    }
    if(this->m_x3hat!=NULL){
        accfft_free(this->m_x3hat);
        this->m_x3hat = NULL;
    }
    if(this->m_Kx1hat!=NULL){
        accfft_free(this->m_Kx1hat);
        this->m_Kx1hat = NULL;
    }
    if(this->m_Kx2hat!=NULL){
        accfft_free(this->m_Kx2hat);
        this->m_Kx2hat = NULL;
    }
    if(this->m_Kx3hat!=NULL){
        accfft_free(this->m_Kx3hat);
        this->m_Kx3hat = NULL;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeBodyForce
 * Description: compute the body force
 * b = K[\int_0^1 \igrad m \lambda d t],
 * where K is an operator that projects v onto the manifold of
 * divergence free velocity fields
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeBodyForce"
PetscErrorCode OptimalControlRegistrationIC::ComputeBodyForce()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_WorkVecField1 == NULL){
        this->m_WorkVecField1 = new VecField(this->m_Opt);
    }
    if (this->m_WorkVecField2 == NULL){
        this->m_WorkVecField2 = new VecField(this->m_Opt);
    }

    // assigned to work vec field 2
    ierr=SuperClass::ComputeBodyForce(); CHKERRQ(ierr);

    ierr=this->m_WorkVecField1->Copy(this->m_WorkVecField2); CHKERRQ(ierr);
    ierr=this->ApplyProjection(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr=this->m_WorkVecField2->AXPY(1.0,this->m_WorkVecField1); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeIncBodyForce
 * Description: compute the body force
 * b = K[\int_0^1\igrad\tilde{m}\lambda+\igrad m \tilde{\lambda} dt]
 * where K is an operator that projects \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeIncBodyForce"
PetscErrorCode OptimalControlRegistrationIC::ComputeIncBodyForce()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_WorkVecField1 == NULL){
        this->m_WorkVecField1 = new VecField(this->m_Opt);
    }
    if (this->m_WorkVecField2 == NULL){
        this->m_WorkVecField2 = new VecField(this->m_Opt);
    }

    // assigned to work vec field 2
    ierr=SuperClass::ComputeIncBodyForce(); CHKERRQ(ierr);

    ierr=this->m_WorkVecField1->Copy(this->m_WorkVecField2); CHKERRQ(ierr);
    ierr=this->ApplyProjection(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr=this->m_WorkVecField2->AXPY(1.0,this->m_WorkVecField1); CHKERRQ(ierr);

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
PetscErrorCode OptimalControlRegistrationIC::SolveAdjointEquationSL()
{
    PetscErrorCode ierr;
    ScalarType *p_l=NULL,*p_lj=NULL,*p_ljX=NULL;
    IntType nl;
    IntType nt;

    PetscFunctionBegin;

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;

    ierr=Assert(this->m_VelocityField!=NULL,"null pointer"); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }
    if(this->m_WorkVecField1==NULL){
        this->m_WorkVecField1 = new VecField(this->m_Opt);
    }

    // scale v by -1
    ierr=this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr=this->m_WorkVecField1->Scale(-1.0); CHKERRQ(ierr);

    // compute trajectory
    ierr=this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1,"adjoint"); CHKERRQ(ierr);

    ierr=VecGetArray(this->m_WorkScaField1,&p_lj); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField2,&p_ljX); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);

    // copy initial condition \lambda = (m_R - m) at t=1
    try{ std::copy(p_l+nt*nl,p_l+(nt+1)*nl,p_lj); }
    catch(std::exception&){
        ierr=ThrowError("copy failed"); CHKERRQ(ierr);
    }


    for (IntType j = 0; j < nt; ++j){

        // compute lambda(t^j,X)
        ierr=this->m_SemiLagrangianMethod->Interpolate(p_ljX,p_lj,"adjoint"); CHKERRQ(ierr);

        // store \lambda(X,t^{j+1})
        try{ std::copy(p_ljX,p_ljX+nl,p_lj); }
        catch(std::exception&){
            ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
        }
        // store time history
        try{ std::copy(p_lj,p_lj+nl,p_l+(nt-(j+1))*nl); }
        catch(std::exception&){
            ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
        }

    }

    ierr=VecRestoreArray(this->m_AdjointVariable,&p_l); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_ljX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField1,&p_lj); CHKERRQ(ierr);

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
PetscErrorCode OptimalControlRegistrationIC::SolveIncAdjointEquationGNSL(void)
{
    PetscErrorCode ierr;
    IntType nl, nt;
    ScalarType *p_ltilde=NULL,*p_ltildej=NULL,*p_ltildejX=NULL;

    PetscFunctionBegin;

    if (this->m_WorkScaField1 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField1); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL){
        ierr=VecDuplicate(this->m_ReferenceImage,&this->m_WorkScaField2); CHKERRQ(ierr);
    }

    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nlocal;

    ierr=VecGetArray(this->m_WorkScaField1,&p_ltildej); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_WorkScaField2,&p_ltildejX); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_IncAdjointVariable,&p_ltilde); CHKERRQ(ierr);

    // remember time history (i.e. copy final condition
    // $\tilde{\lambda}_1 = -\tilde{m}_1$ into buffer for $\tilde{\lambda}
    try{ std::copy(p_ltilde+nt*nl,p_ltilde+(nt+1)*nl,p_ltildej); }
    catch(std::exception&){
        ierr=ThrowError("copying of data failed"); CHKERRQ(ierr);
    }

    // for all time points
    for (IntType j = 0; j < nt; ++j){

        ierr=this->m_SemiLagrangianMethod->Interpolate(p_ltildejX,p_ltildej,"adjoint"); CHKERRQ(ierr);

        // store time history (necessary for optimization)
        try{ std::copy(p_ltildejX,p_ltildejX+nl,p_ltildej); }
        catch(std::exception&){
            ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        // store time history (necessary for optimization)
        try{ std::copy(p_ltildej,p_ltildej+nl,p_ltilde+(nt-(j+1))*nl); }
        catch(std::exception&){
            ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
        }

    }

    ierr=VecRestoreArray(this->m_IncAdjointVariable,&p_ltilde); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField2,&p_ltildejX); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_WorkScaField1,&p_ltildej); CHKERRQ(ierr);


    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: ApplyProjection
 * Description: apply projection to map \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyProjection"
PetscErrorCode OptimalControlRegistrationIC::ApplyProjection(VecField* x)
{
    PetscErrorCode ierr;
    ScalarType *p_x1=NULL, *p_x2=NULL, *p_x3=NULL, scale;
    long int nx[3];
    IntType nalloc;
    double ffttimers[5]={0,0,0,0,0};

    PetscFunctionBegin;

    nx[0] = static_cast<long int>(this->m_Opt->GetNumGridPoints(0));
    nx[1] = static_cast<long int>(this->m_Opt->GetNumGridPoints(1));
    nx[2] = static_cast<long int>(this->m_Opt->GetNumGridPoints(2));


    nalloc = this->m_Opt->GetFFT().nalloc;
    scale = this->m_Opt->ComputeFFTScale();

    if(this->m_x1hat == NULL){
        this->m_x1hat=(FFTScaType*)accfft_alloc(nalloc);
    }
    if(this->m_x2hat == NULL){
        this->m_x2hat=(FFTScaType*)accfft_alloc(nalloc);
    }
    if(this->m_x3hat == NULL){
        this->m_x3hat=(FFTScaType*)accfft_alloc(nalloc);
    }

    if(this->m_Kx1hat == NULL){
        this->m_Kx1hat=(FFTScaType*)accfft_alloc(nalloc);
    }
    if(this->m_Kx2hat == NULL){
        this->m_Kx2hat=(FFTScaType*)accfft_alloc(nalloc);
    }
    if(this->m_Kx3hat == NULL){
        this->m_Kx3hat=(FFTScaType*)accfft_alloc(nalloc);
    }

    ierr=VecGetArray(x->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecGetArray(x->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecGetArray(x->m_X3,&p_x3); CHKERRQ(ierr);

    // compute forward fft
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x1,this->m_x1hat,ffttimers);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x2,this->m_x2hat,ffttimers);
    accfft_execute_r2c_t<ScalarType,FFTScaType>(this->m_Opt->GetFFT().plan,p_x3,this->m_x3hat,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);


#pragma omp parallel
{
    long int x1,x2,x3,wx1,wx2,wx3;
    ScalarType lapinvik,gradik1,gradik2,gradik3;
    IntType i,i1,i2,i3;
#pragma omp for
    for (i1 = 0; i1 < this->m_Opt->GetFFT().osize[0]; ++i1){
        for (i2 = 0; i2 < this->m_Opt->GetFFT().osize[1]; ++i2){
            for (i3 = 0; i3 < this->m_Opt->GetFFT().osize[2]; ++i3){

                x1 = static_cast<long int>(i1 + this->m_Opt->GetFFT().ostart[0]);
                x2 = static_cast<long int>(i2 + this->m_Opt->GetFFT().ostart[1]);
                x3 = static_cast<long int>(i3 + this->m_Opt->GetFFT().ostart[2]);

                // set wavenumber
                wx1 = x1;
                wx2 = x2;
                wx3 = x3;

                if(x1 > nx[0]/2) wx1-=nx[0];
                if(x2 > nx[1]/2) wx2-=nx[1];
                if(x3 > nx[2]/2) wx3-=nx[2];

                // compute inverse laplacian operator
                lapinvik = static_cast<ScalarType>(wx1*wx1 + wx2*wx2 + wx3*wx3);
                //lapinvik = round(lapinvik) == 0.0 ? -1.0 : 1.0/lapinvik;
                lapinvik = lapinvik == 0.0 ? -1.0 : -1.0/lapinvik;

                if(x1 == nx[0]/2) wx1 = 0;
                if(x2 == nx[1]/2) wx2 = 0;
                if(x3 == nx[2]/2) wx3 = 0;

                // compute gradient operator
                gradik1 = static_cast<ScalarType>(wx1);
                gradik2 = static_cast<ScalarType>(wx2);
                gradik3 = static_cast<ScalarType>(wx3);

                i=GetLinearIndex(i1,i2,i3,this->m_Opt->GetFFT().osize);

                // compute div(b)
                this->m_Kx1hat[i][0] = -scale*(gradik1*this->m_x1hat[i][0]
                                             + gradik2*this->m_x2hat[i][0]
                                             + gradik3*this->m_x3hat[i][0]);

                this->m_Kx1hat[i][1] =  scale*(gradik1*this->m_x1hat[i][1]
                                             + gradik2*this->m_x2hat[i][1]
                                             + gradik3*this->m_x3hat[i][1]);

                // compute lap^{-1} div(b)
                this->m_Kx1hat[i][0] *= lapinvik;
                this->m_Kx1hat[i][1] *= lapinvik;

                // compute x2 gradient of lab^{-1} div(b)
                this->m_Kx2hat[i][0] = -gradik2*this->m_Kx1hat[i][0];
                this->m_Kx2hat[i][1] =  gradik2*this->m_Kx1hat[i][1];

                // compute x3 gradient of lab^{-1} div(b)
                this->m_Kx3hat[i][0] = -gradik3*this->m_Kx1hat[i][0];
                this->m_Kx3hat[i][1] =  gradik3*this->m_Kx1hat[i][1];

                // compute x1 gradient of lab^{-1} div(b)
                this->m_Kx1hat[i][0] *= -gradik1;
                this->m_Kx1hat[i][1] *=  gradik1;

            }
        }
    }
}// pragma omp parallel


    // compute inverse fft
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Kx1hat,p_x1,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Kx2hat,p_x2,ffttimers);
    accfft_execute_c2r_t<FFTScaType,ScalarType>(this->m_Opt->GetFFT().plan,this->m_Kx3hat,p_x3,ffttimers);
    this->m_Opt->IncrementCounter(FFT,3);

    ierr=VecRestoreArray(x->m_X1,&p_x1); CHKERRQ(ierr);
    ierr=VecRestoreArray(x->m_X2,&p_x2); CHKERRQ(ierr);
    ierr=VecRestoreArray(x->m_X3,&p_x3); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(ffttimers);

    PetscFunctionReturn(0);
}


} // end of namespace


#endif // _OPTIMALCONTROLREGISTRATIONIC_CPP_
