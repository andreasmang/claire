#ifndef _SYNPROBREGISTRATION_CPP_
#define _SYNPROBREGISTRATION_CPP_

#include "SynProbRegistration.hpp"


namespace reg
{




/********************************************************************
 * Name: SynProbRegistration
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SynProbRegistration"
SynProbRegistration::SynProbRegistration()
{
    this->Initialize();
}


/********************************************************************
 * Name: SynProbRegistration
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SynProbRegistration"
SynProbRegistration::SynProbRegistration(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * Name: SynProbRegistration
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~SynProbRegistration"
SynProbRegistration::~SynProbRegistration()
{
    this->ClearMemory();
}


/********************************************************************
 * Name: Initialize
 * Description: initialize class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode SynProbRegistration::Initialize()
{

    this->m_Opt = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clear memory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode SynProbRegistration::ClearMemory()
{

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: ComputeSmoothScalarField
 * Description: clear memory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeSmoothScalarField"
PetscErrorCode SynProbRegistration::ComputeSmoothScalarField(Vec m,const unsigned int id)
{
    PetscErrorCode ierr;
    unsigned int isize[3],istart[3];
    ScalarType *p_m=NULL;
    ScalarType value,hx[3];
    const ScalarType twopi = 2.0*PETSC_PI;
    const ScalarType sigma = 160.0;

    PetscFunctionBegin;

    ierr=Assert(m!= NULL,"null pointer"); CHKERRQ(ierr);

    for (unsigned int i = 0; i < 3; ++i){
        hx[i]     = this->m_Opt->m_MiscOpt->h[i];
        isize[i]  = this->m_Opt->m_MiscOpt->isize[i];
        istart[i] = this->m_Opt->m_MiscOpt->istart[i];
    }

    ierr=VecGetArray(m,&p_m); CHKERRQ(ierr);

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
                IntType i = GetLinearIndex(i1,i2,i3,isize);

                if (id == 0){
                    ScalarType s1 = sin(x1);
                    ScalarType s2 = sin(x2);
                    ScalarType s3 = sin(x3);
                    value = (s1*s1 + s2*s2 + s3*s3)/3.0;
                }
                else if (id == 1){
                   // first derivative of id 0 with respect to x1
                    value = 2.0*sin(x1)*cos(x1)/3.0;
                }
                else if (id == 2){
                    // first derivative of id 0 with respect to x2
                    value = 2.0*sin(x2)*cos(x2)/3.0;
                }
                else if (id == 3){
                    // first derivative of id 0 with respect to x3
                    value = 2.0*sin(x3)*cos(x3)/3.0;
                }
                else if (id == 4){
                    // laplacian of id 0
                    ScalarType s1 = 2.0*cos(x1)*cos(x1) - 2.0*sin(x1)*sin(x1);
                    ScalarType s2 = 2.0*cos(x2)*cos(x2) - 2.0*sin(x2)*sin(x2);
                    ScalarType s3 = 2.0*cos(x3)*cos(x3) - 2.0*sin(x3)*sin(x3);
                    value = (s1 + s2 + s3)/3.0;
                }
                else if (id == 5){
                    // biharmonic of id 0
                    ScalarType s1 = 8.0*sin(x1)*sin(x1) - 8.0*cos(x1)*cos(x1);
                    ScalarType s2 = 8.0*sin(x2)*sin(x2) - 8.0*cos(x2)*cos(x2);
                    ScalarType s3 = 8.0*sin(x3)*sin(x3) - 8.0*cos(x3)*cos(x3);
                    value = (s1 + s2 + s3)/3.0;
                }
                else if (id == 6){
                    // rhs for poisson problem
                    ScalarType c = sigma/(twopi*twopi);
                    ScalarType x1s = (x1-PETSC_PI)*(x1-PETSC_PI);
                    ScalarType x2s = (x2-PETSC_PI)*(x2-PETSC_PI);
                    ScalarType x3s = (x3-PETSC_PI)*(x3-PETSC_PI);
                    ScalarType y = - exp( -c * (x1s + x2s + x3s) );
                    value = - ( -6.0*c + 4.0*c*c*(x1s + x2s + x3s) )*y;
                }
                else if (id == 7){

                    // analytic solution for rhs of poisson problem (id 6)
                    // rhs for poisson problem
                    ScalarType c = sigma/(twopi*twopi);
                    ScalarType x1s = (x1-PETSC_PI)*(x1-PETSC_PI);
                    ScalarType x2s = (x2-PETSC_PI)*(x2-PETSC_PI);
                    ScalarType x3s = (x3-PETSC_PI)*(x3-PETSC_PI);

                    ScalarType y = -exp(-c * (x1s + x2s + x3s) );
                    value = y + 1.0;
                }
                else if (id == 8){
                    value  = 0.5 + 0.5*sin(x1)*sin(x2)*sin(x3);
                }
                else if (id == 9){
                    value = sin(2.0*x1)*sin(2.0*x2)*sin(2.0*x3);
                    value *= value;
                }
                else { value = 0.0; }

                p_m[i] = value;

            } // i1
        } // i2
    } // i3
} // pragma omp parallel

    ierr=VecRestoreArray(m,&p_m); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeSquare
 * Description: compute square
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeSquare"
PetscErrorCode SynProbRegistration::ComputeSquare(Vec m)
{
    PetscErrorCode ierr;
    unsigned int isize[3],istart[3];
    ScalarType *p_m=NULL,hx[3];

    PetscFunctionBegin;

    for (unsigned int i = 0; i < 3; ++i){
        hx[i] = this->m_Opt->m_MiscOpt->h[i];
        isize[i] = this->m_Opt->m_MiscOpt->isize[i];
        istart[i] = this->m_Opt->m_MiscOpt->istart[i];
    }

    ierr=VecGetArray(m,&p_m); CHKERRQ(ierr);

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

                x1 -= PETSC_PI; x1 = fabs(x1);
                x2 -= PETSC_PI; x2 = fabs(x2);
                x3 -= PETSC_PI; x3 = fabs(x3);

                // compute linear / flat index
                IntType i = GetLinearIndex(i1,i2,i3,isize);

                p_m[i] = (std::max(std::max(x1,x2),x3) < 2.0*PETSC_PI/4.0) ? 1.0 : 0.0;

            } // i1
        } // i2
    } // i3
} // pragma omp parallel

    ierr=VecRestoreArray(m,&p_m); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeExpSin
 * Description: compute exp sin field
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeExpSin"
PetscErrorCode SynProbRegistration::ComputeExpSin(Vec m)
{
    PetscErrorCode ierr;
    unsigned int isize[3],istart[3];
    const ScalarType twopi = 2.0*PETSC_PI;
    ScalarType *p_m=NULL,hx[3],sigma[3];

    PetscFunctionBegin;

    for (unsigned int i = 0; i < 3; ++i){
        hx[i] = this->m_Opt->m_MiscOpt->h[i];
        isize[i] = this->m_Opt->m_MiscOpt->isize[i];
        istart[i] = this->m_Opt->m_MiscOpt->istart[i];
    }

    ierr=VecGetArray(m,&p_m); CHKERRQ(ierr);

    sigma[0] = twopi/8.0;
    sigma[1] = twopi/8.0;
    sigma[2] = twopi/8.0;

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
                IntType i = GetLinearIndex(i1,i2,i3,isize);

                ScalarType s1 = sin(2.0*x1);
                ScalarType s2 = sin(2.0*x2);
                ScalarType s3 = sin(2.0*x3);

                ScalarType x1s = x1 - twopi/2.0;
                ScalarType x2s = x2 - twopi/2.0;
                ScalarType x3s = x3 - twopi/2.0;
                x1s = x1s*x1s / (2.0*sigma[0]*sigma[0]);
                x2s = x2s*x2s / (2.0*sigma[1]*sigma[1]);
                x3s = x3s*x3s / (2.0*sigma[2]*sigma[2]);

                p_m[i] = exp( - (x1s + x2s + x3s) ) * (s1*s1 + s2*s2 + s3*s3);


            } // i1
        } // i2
    } // i3
} // pragma omp parallel

    ierr=VecRestoreArray(m,&p_m); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ComputeDiamond
 * Description: compute diamond
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDiamond"
PetscErrorCode SynProbRegistration::ComputeDiamond(Vec m,const unsigned int id)
{
    PetscErrorCode ierr;
    unsigned int isize[3],istart[3];
    ScalarType *p_m=NULL,hx[3];

    PetscFunctionBegin;

    for (unsigned int i = 0; i < 3; ++i){
        hx[i] = this->m_Opt->m_MiscOpt->h[i];
        isize[i] = this->m_Opt->m_MiscOpt->isize[i];
        istart[i] = this->m_Opt->m_MiscOpt->istart[i];
    }

    ierr=VecGetArray(m,&p_m); CHKERRQ(ierr);


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
                IntType i = GetLinearIndex(i1,i2,i3,isize);

                x1 = fabs(x1-PETSC_PI);
                x2 = fabs(x2-PETSC_PI);
                x3 = fabs(x3-PETSC_PI);

                ScalarType x=0;

                if (id == 1){

                    x = x1 + x2 + x3;

                    if(x < sqrt(3.0)*PETSC_PI/2.0) p_m[i] =  0.4;
                    if(x < sqrt(3.0)*PETSC_PI/3.0) p_m[i] =  0.6;
                    if(x < sqrt(3.0)*PETSC_PI/4.0) p_m[i] =  0.8;
                    if(x < sqrt(3.0)*PETSC_PI/8.0) p_m[i] =  1.0;

                }
                else{

                    x = std::max(std::max(x1,x2),x3);

                    if(x < PETSC_PI/2.0) p_m[i] =  0.4;
                    if(x < PETSC_PI/3.0) p_m[i] =  0.6;
                    if(x < PETSC_PI/4.0) p_m[i] =  0.8;
                    if(x < PETSC_PI/8.0) p_m[i] =  1.0;

                }

            } // i1
        } // i2
    } // i3
} // pragma omp parallel

    ierr=VecRestoreArray(m,&p_m); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of name space


#endif // _SYNPROBREGISTRATION_CPP_
