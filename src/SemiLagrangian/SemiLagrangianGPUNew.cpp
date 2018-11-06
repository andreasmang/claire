/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CLAIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _SEMILAGRANGIANGPUNEW_CPP_
#define _SEMILAGRANGIANGPUNEW_CPP_

#include "SemiLagrangianGPUNew.hpp"
#include <petsc/private/vecimpl.h>



namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
SemiLagrangianGPUNew::SemiLagrangianGPUNew() {
    this->Initialize();
}




/********************************************************************
 * @brief default constructor
 *******************************************************************/
SemiLagrangianGPUNew::SemiLagrangianGPUNew(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
    this->InitializeInterpolationTexture();
    this->ComputeInitialTrajectory();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
SemiLagrangianGPUNew::~SemiLagrangianGPUNew() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init class variables
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_X = NULL;
    this->m_WorkVecField1 = NULL;
    this->m_InitialTrajectory = NULL;
    this->m_texture = 0;

    this->m_Opt = NULL;
    this->m_Dofs[0] = 1;
    this->m_Dofs[1] = 3;


    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief init empty texture for interpolation on GPU 
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::InitializeInterpolationTexture() {
    PetscErrorCode ierr = 0;
    int nx[3];
    PetscFunctionBegin;
    
    for (unsigned int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
    }

    this->m_texture = gpuInitEmptyTexture(nx);
    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief clears memory //TODO this function will become obsolete in the case we
 * perform everything on the GPU
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_X != NULL) {
        delete this->m_X;
        this->m_X = NULL;
    }

    if (this->m_InitialTrajectory != NULL) {
        delete this->m_InitialTrajectory;
        this->m_InitialTrajectory = NULL;
    }

    if (this->m_texture != 0) {
        cudaDestroyTextureObject(this->m_texture);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set work vector field to not have to allocate it locally
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::SetWorkVecField(VecField* x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_WorkVecField1 = x;

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief compute the trajectory from the velocity field based
 * on an rk2 scheme (todo: make the velocity field a const vector)
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::ComputeTrajectory(VecField* v, std::string flag) {
    PetscErrorCode ierr = 0;
    IntType nl;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->m_Domain.nl;

    // compute trajectory by calling a CUDA kernel
    if (this->m_Opt->m_PDESolver.rkorder == 2) {
        ierr = this->ComputeTrajectoryRK2(v, flag); CHKERRQ(ierr);
    } else if (this->m_Opt->m_PDESolver.rkorder == 4) {
        ierr = this->ComputeTrajectoryRK4(v, flag); CHKERRQ(ierr);
    } else {
        ierr = ThrowError("rk order not implemented"); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief compute the initial trajectory
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::ComputeInitialTrajectory() {
    PetscErrorCode ierr;
    ScalarType *p_x1=NULL,*p_x2=NULL,*p_x3=NULL;
    IntType isize[3],istart[3];
    int nx[3];
    
    ScalarType hx[3];
    IntType l,i1,i2,i3;
    PetscFunctionBegin;

    if(this->m_InitialTrajectory==NULL){
        try{this->m_InitialTrajectory = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    for (unsigned int i = 0; i < 3; ++i){
        hx[i]     = this->m_Opt->m_Domain.hx[i];
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i]  = this->m_Opt->m_Domain.isize[i];
        istart[i] = this->m_Opt->m_Domain.istart[i];
    }
    
    // allocate 1 time 3D Cuda Array for texture interpolation
    
    ierr = VecGetArray(this->m_InitialTrajectory->m_X1, &p_x1); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_InitialTrajectory->m_X2, &p_x2); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_InitialTrajectory->m_X3, &p_x3); CHKERRQ(ierr);


#pragma omp parallel
{
#pragma omp for
    for (unsigned int i1 = 0; i1 < isize[0]; ++i1){  // x1
        for (unsigned int i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (unsigned int i3 = 0; i3 < isize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                ScalarType x1 = static_cast<ScalarType>(i1 + istart[0]);
                ScalarType x2 = static_cast<ScalarType>(i2 + istart[1]);
                ScalarType x3 = static_cast<ScalarType>(i3 + istart[2]);

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


//   ierr = PrintVectorMemoryLocation(this->m_InitialTrajectory->m_X1, "Initial trajectory before"); CHKERRQ(ierr);
/*    ierr = GetRawPointer(this->m_InitialTrajectory->m_X1, &p_x1); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_InitialTrajectory->m_X2, &p_x2); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_InitialTrajectory->m_X3, &p_x3); CHKERRQ(ierr);
    
    getSemiLagrangianInitialCondition(p_x1, p_x2, p_x3, nx, &(this->m_Opt->m_GPUtime));
    
    ierr = RestoreRawPointer(this->m_InitialTrajectory->m_X1, &p_x1); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_InitialTrajectory->m_X2, &p_x2); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_InitialTrajectory->m_X3, &p_x3); CHKERRQ(ierr);
*/
//    ierr = PrintVectorMemoryLocation(this->m_InitialTrajectory->m_X1, "Initial trajectory after"); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief compute the trajectory from the velocity field based
 * on an rk2 scheme (todo: make the velocity field a const vector)
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::ComputeTrajectoryRK2(VecField* v, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType ht, hthalf, hx[3], invhx[3], x1, x2, x3, scale = 0.0;
    const ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL;
    ScalarType *w_X1 = NULL, *w_X2 = NULL, *w_X3 = NULL;
    ScalarType *p_vX1 = NULL, *p_vX2 = NULL, *p_vX3 = NULL;
    IntType isize[3], istart[3], l, i1, i2, i3;
    std::stringstream ss;
    double runtime=0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if(this->m_X == NULL){
        try{this->m_X = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr = Assert(this->m_WorkVecField1 != NULL, "null pointer"); CHKERRQ(ierr);
    
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_InitialTrajectory == NULL){
        ierr=this->ComputeInitialTrajectory(); CHKERRQ(ierr);
    }
    
    // switch between state and adjoint variable
    if (strcmp(flag.c_str(), "state") == 0) {
        scale =  1.0;
    } else if (strcmp(flag.c_str(), "adjoint") == 0) {
        scale = -1.0;
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }


    for (int i = 0; i < 3; ++i) {
        hx[i]     = this->m_Opt->m_Domain.hx[i];
        invhx[i]  = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[i])/(PETSC_PI*2.0);
        isize[i]  = this->m_Opt->m_Domain.isize[i];
        istart[i] = this->m_Opt->m_Domain.istart[i];
    }
    
    // X = x - ht v
    ierr = VecWAXPY(this->m_X->m_X1, -scale*ht*invhx[0], v->m_X1, this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X->m_X2, -scale*ht*invhx[0], v->m_X2, this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X->m_X3, -scale*ht*invhx[0], v->m_X3, this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);
    
    // interpolate velocity field v(X)
    ierr = this->Interpolate(this->m_WorkVecField1, v, flag); CHKERRQ(ierr);
    
    // X = x - 0.5*ht*(v + v(x - ht v))
    // F = F0 + F1 = v + v(x-ht*v)
    ierr = VecAXPY(this->m_WorkVecField1->m_X1, 1.0, v->m_X1); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_WorkVecField1->m_X2, 1.0, v->m_X2); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_WorkVecField1->m_X3, 1.0, v->m_X3); CHKERRQ(ierr);

    // X = x - 0.5*ht*F
    ierr = VecWAXPY(this->m_X->m_X1, -scale*hthalf*invhx[0], this->m_WorkVecField1->m_X1, this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X->m_X2, -scale*hthalf*invhx[0], this->m_WorkVecField1->m_X2, this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(this->m_X->m_X3, -scale*hthalf*invhx[0], this->m_WorkVecField1->m_X3, this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);
  
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief compute the trajectory from the velocity field based
 * on an rk2 scheme (todo: make the velocity field a const vector)
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::ComputeTrajectoryRK4(VecField* v, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType ht, hthalf, hx[3], x1, x2, x3, scale = 0.0;
    const ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL;
    ScalarType *p_vX1 = NULL, *p_vX2 = NULL, *p_vX3 = NULL,
               *p_f1 = NULL, *p_f2 = NULL, *p_f3 = NULL;
    IntType isize[3], istart[3], l, i1, i2, i3;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    /*
    ierr = Assert(this->m_WorkVecField1 != NULL, "null pointer"); CHKERRQ(ierr);
    if (this->m_WorkVecField2 == NULL) {
        try{this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    // switch between state and adjoint variable
    if (strcmp(flag.c_str(), "state") == 0) {
        scale =  1.0;
    } else if (strcmp(flag.c_str(), "adjoint") == 0) {
        scale = -1.0;
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }


    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->m_Domain.hx[i];
        isize[i] = this->m_Opt->m_Domain.isize[i];
        istart[i] = this->m_Opt->m_Domain.istart[i];
    }

    ierr = this->m_WorkVecField2->GetArrays(p_f1, p_f2, p_f3); CHKERRQ(ierr);

    // first stage of rk4
    ierr = v->GetArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    for (i1 = 0; i1 < isize[0]; ++i1) {   // x1
        for (i2 = 0; i2 < isize[1]; ++i2) {   // x2
            for (i3 = 0; i3 < isize[2]; ++i3) {   // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                l = GetLinearIndex(i1, i2, i3, isize);

                p_f1[l] = p_v1[l];
                p_f2[l] = p_v2[l];
                p_f3[l] = p_v3[l];

                this->m_X[l*3+0] = (x1 - scale*hthalf*p_v1[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+1] = (x2 - scale*hthalf*p_v2[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+2] = (x3 - scale*hthalf*p_v3[l])/(2.0*PETSC_PI); // normalized to [0,1]
            }  // i1
        }  // i2
    }  // i3
    ierr = v->RestoreArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    // evaluate right hand side
    ierr = this->CommunicateCoord(flag); CHKERRQ(ierr);
    ierr = this->Interpolate(this->m_WorkVecField1, v, flag); CHKERRQ(ierr);

    // second stage of rk4
    ierr = this->m_WorkVecField1->GetArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
    for (i1 = 0; i1 < isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                l = GetLinearIndex(i1, i2, i3, isize);

                p_f1[l] += 2.0*p_vX1[l];
                p_f2[l] += 2.0*p_vX2[l];
                p_f3[l] += 2.0*p_vX3[l];

                this->m_X[l*3+0] = (x1 - scale*hthalf*p_vX1[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+1] = (x2 - scale*hthalf*p_vX2[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+2] = (x3 - scale*hthalf*p_vX3[l])/(2.0*PETSC_PI); // normalized to [0,1]
            }  // i1
        }  // i2
    }  // i3
    ierr = this->m_WorkVecField1->RestoreArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);

    // evaluate right hand side
    ierr = this->CommunicateCoord(flag); CHKERRQ(ierr);
    ierr = this->Interpolate(this->m_WorkVecField1, v, flag); CHKERRQ(ierr);

    // third stage of rk4
    ierr = this->m_WorkVecField1->GetArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
    for (i1 = 0; i1 < isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                l = GetLinearIndex(i1, i2, i3, isize);

                p_f1[l] += 2.0*p_vX1[l];
                p_f2[l] += 2.0*p_vX2[l];
                p_f3[l] += 2.0*p_vX3[l];

                this->m_X[l*3+0] = (x1 - scale*ht*p_vX1[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+1] = (x2 - scale*ht*p_vX2[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+2] = (x3 - scale*ht*p_vX3[l])/(2.0*PETSC_PI); // normalized to [0,1]
            }  // i1
        }  // i2
    }  // i3
    ierr = this->m_WorkVecField1->RestoreArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);

    // evaluate right hand side
    ierr = this->CommunicateCoord(flag); CHKERRQ(ierr);
    ierr = this->Interpolate(this->m_WorkVecField1, v, flag); CHKERRQ(ierr);

    // fourth stage of rk4
    ierr = this->m_WorkVecField1->GetArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
    for (i1 = 0; i1 < isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                l = GetLinearIndex(i1, i2, i3, isize);

                p_f1[l] += p_vX1[l];
                p_f2[l] += p_vX2[l];
                p_f3[l] += p_vX3[l];

                this->m_X[l*3+0] = (x1 - scale*(ht/6.0)*p_f1[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+1] = (x2 - scale*(ht/6.0)*p_f2[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+2] = (x3 - scale*(ht/6.0)*p_f3[l])/(2.0*PETSC_PI); // normalized to [0,1]
            }  // i1
        }  // i2
    }  // i3
    ierr = this->m_WorkVecField1->RestoreArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);

    ierr = this->m_WorkVecField2->RestoreArrays(p_f1, p_f2, p_f3); CHKERRQ(ierr);

    // communicate the final characteristic
    ierr = this->CommunicateCoord(flag); CHKERRQ(ierr);
    */
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief interpolate scalar field
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Interpolate(Vec* xo, Vec xi, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType *p_xo = NULL, *p_xi = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(*xo != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert( xi != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = GetRawPointerReadWrite( xi, &p_xi); CHKERRQ(ierr);
    ierr = GetRawPointerReadWrite(*xo, &p_xo); CHKERRQ(ierr);

    ierr = this->Interpolate(p_xo, p_xi, flag); CHKERRQ(ierr);

    ierr = RestoreRawPointerReadWrite(*xo, &p_xo); CHKERRQ(ierr);
    ierr = RestoreRawPointerReadWrite( xi, &p_xi); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief interpolate scalar field
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Interpolate(ScalarType* xo, ScalarType* xi, std::string flag) {
    PetscErrorCode ierr = 0;
    int nx[3], isize_g[3], isize[3], istart_g[3], istart[3], c_dims[2], neval, order, nghost;
    IntType nl, nalloc;
    std::stringstream ss;
    double timers[4] = {0, 0, 0, 0};
    const ScalarType *xq1, *xq2, *xq3;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(xi != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xo != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    nl     = this->m_Opt->m_Domain.nl;
    order  = this->m_Opt->m_PDESolver.iporder;
    nghost = order;
    neval  = static_cast<int>(nl);

    for (int i = 0; i < 3; ++i) {
        nx[i]     = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i]  = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }

    ierr = this->m_X->GetArraysRead(xq1, xq2, xq3);
    
     
    // compute interpolation for all components of the input scalar field
    if (strcmp(flag.c_str(), "state") == 0) {
        gpuInterp3D(xi, xq1, xq2, xq3, xo, nx, this->m_texture, &(this->m_Opt->m_GPUtime));
    } else if (strcmp(flag.c_str(), "adjoint") == 0) {
        gpuInterp3D(xi, xq1, xq2, xq3, xo, nx, this->m_texture, &(this->m_Opt->m_GPUtime));
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }
    
    ierr = this->m_X->RestoreArraysRead(xq1, xq2, xq3);
    
    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);
    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IP);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief interpolate vector field
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Interpolate(VecField* vo, VecField* vi, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType *p_vix1 = NULL, *p_vix2 = NULL, *p_vix3 = NULL;
    ScalarType *p_vox1 = NULL, *p_vox2 = NULL, *p_vox3 = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vi != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vo != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = vi->GetArraysReadWrite(p_vix1, p_vix2, p_vix3); CHKERRQ(ierr);
    ierr = vo->GetArraysReadWrite(p_vox1, p_vox2, p_vox3); CHKERRQ(ierr);
    
    ierr = this->Interpolate(p_vox1, p_vox2, p_vox3, p_vix1, p_vix2, p_vix3, flag); CHKERRQ(ierr);
    
    ierr = vi->RestoreArraysReadWrite(p_vix1, p_vix2, p_vix3); CHKERRQ(ierr);
    ierr = vo->RestoreArraysReadWrite(p_vox1, p_vox2, p_vox3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief interpolate vector field
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Interpolate(ScalarType* wx1, ScalarType* wx2, ScalarType* wx3,
                                                 ScalarType* vx1, ScalarType* vx2, ScalarType* vx3, std::string flag) {
    PetscErrorCode ierr = 0;
    int nx[3], isize_g[3], isize[3], istart_g[3], istart[3], c_dims[2], nghost, order;
    double timers[4] = {0, 0, 0, 0};
    std::stringstream ss;
    IntType nl, nlghost, nalloc;
    const ScalarType *xq1, *xq2, *xq3;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vx1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vx2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vx3 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx3 != NULL, "null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->m_Domain.nl;
    order = this->m_Opt->m_PDESolver.iporder;
    nghost = order;

    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i] = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }

    
    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    // Get the query points
    ierr = this->m_X->GetArraysRead(xq1, xq2, xq3); CHKERRQ(ierr);

    if (strcmp(flag.c_str(),"state") == 0) {
        gpuInterp3D(vx1, xq1, xq2, xq3, wx1, nx, this->m_texture, &(this->m_Opt->m_GPUtime));
        gpuInterp3D(vx2, xq1, xq2, xq3, wx2, nx, this->m_texture, &(this->m_Opt->m_GPUtime));
        gpuInterp3D(vx3, xq1, xq2, xq3, wx3, nx, this->m_texture, &(this->m_Opt->m_GPUtime));
    } else if (strcmp(flag.c_str(),"adjoint") == 0) {
        gpuInterp3D(vx1, xq1, xq2, xq3, wx1, nx, this->m_texture, &(this->m_Opt->m_GPUtime));
        gpuInterp3D(vx2, xq1, xq2, xq3, wx2, nx, this->m_texture, &(this->m_Opt->m_GPUtime));
        gpuInterp3D(vx3, xq1, xq2, xq3, wx3, nx, this->m_texture, &(this->m_Opt->m_GPUtime));
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

    ierr = this->m_X->RestoreArraysRead(xq1, xq2, xq3); CHKERRQ(ierr);

    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IPVEC);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}


/********************************************************************
 * @brief set coordinate vector and communicate to interpolation plan
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::SetQueryPoints(ScalarType* y1, ScalarType* y2, ScalarType* y3, std::string flag) {
    PetscErrorCode ierr = 0;
    IntType nl;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->m_Domain.nl;

    // if query points have not yet been allocated
    /*
    if (this->m_X == NULL) {
        try {this->m_X = new ScalarType[3*nl];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // copy data to a flat vector
    for (IntType i = 0; i < nl; ++i) {
        this->m_X[0*nl+i] = y1[i];
        this->m_X[1*nl+i] = y2[i];
        this->m_X[2*nl+i] = y3[i];
    }

    // evaluate right hand side
    ierr = this->CommunicateCoord(flag); CHKERRQ(ierr);
    */
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief communicate the coordinate vector (query points)
 * @param flag to switch between forward and adjoint solves
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::CommunicateCoord(std::string flag) {
    PetscErrorCode ierr;
    int nx[3], nl, isize[3], istart[3], nghost;
    int c_dims[2];
    double timers[4] = {0, 0, 0, 0};
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    
    /*
    ierr = Assert(this->m_Opt->m_FFT.mpicomm != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nl     = static_cast<int>(this->m_Opt->m_Domain.nl);
    nghost = this->m_Opt->m_PDESolver.iporder;
    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i] = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }

    // get network dimensions
    c_dims[0] = this->m_Opt->m_CartGridDims[0];
    c_dims[1] = this->m_Opt->m_CartGridDims[1];

    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    if (strcmp(flag.c_str(), "state") == 0) {
        // characteristic for state equation should have been computed already
        ierr = Assert(this->m_X != NULL, "null pointer"); CHKERRQ(ierr);
        // create planer
        if (this->m_StatePlan == NULL) {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("allocating state plan"); CHKERRQ(ierr);
            }
            try {this->m_StatePlan = new Interp3_Plan();}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
            this->m_StatePlan->allocate(nl, this->m_Dofs, 2);
        }

        // scatter
        this->m_StatePlan->scatter(nx, isize, istart, nl, nghost, this->m_X,
                                   c_dims, this->m_Opt->m_FFT.mpicomm, timers);
    } else if (strcmp(flag.c_str(), "adjoint") == 0) {
        // characteristic for adjoint equation should have been computed already
        ierr = Assert(this->m_X != NULL, "null pointer"); CHKERRQ(ierr);
        // create planer
        if (this->m_AdjointPlan == NULL) {
            if (this->m_Opt->m_Verbosity > 2) {
                ierr = DbgMsg("allocating adjoint plan"); CHKERRQ(ierr);
            }
            try {this->m_AdjointPlan = new Interp3_Plan();}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
            this->m_AdjointPlan->allocate(nl, this->m_Dofs, 2);
        }

        // communicate coordinates
        this->m_AdjointPlan->scatter(nx, isize, istart, nl, nghost, this->m_X,
                                     c_dims, this->m_Opt->m_FFT.mpicomm, timers);
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);

    this->m_Opt->IncreaseInterpTimers(timers);
    */
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}

}  // namespace reg




#endif  // _SEMILAGRANGIAN_CPP_
