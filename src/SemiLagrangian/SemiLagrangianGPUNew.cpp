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
#include "SemiLagrangianKernel.hpp"
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
    this->m_Opt = opt;
    this->Initialize();
    
    if (opt->m_Verbosity > 2) {
      DbgMsg("SemiLagrangianGPUNew created");
    }
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
    IntType isize[3];
    PetscFunctionBegin;
    
    this->m_Xstate = nullptr;
    this->m_Xadjoint = nullptr;
    this->m_X = nullptr;
    this->m_WorkScaField1 = nullptr;
    this->m_WorkScaField2 = nullptr;
    this->m_ScaFieldGhost = nullptr;
    this->m_VecFieldGhost = nullptr;
    this->m_WorkVecField1 = nullptr;
    this->m_InitialTrajectory = nullptr;
    this->m_texture = 0;

    this->m_Dofs[0] = 1;
    this->m_Dofs[1] = 3;
    
    this->m_tmpInterpol1 = nullptr;
    this->m_tmpInterpol2 = nullptr;
    
    this->m_StatePlan = nullptr;
    this->m_AdjointPlan = nullptr;
    this->m_GhostPlan = nullptr;

    if (this->m_Opt->rank_cnt > 1) {
      IntType nl = this->m_Opt->m_Domain.nl;
      IntType ng = this->m_Opt->m_Domain.ng;
      for (int i=0; i<3; i++) {
        isize[i] = this->m_Opt->m_Domain.isize[i];
      }
      
      this->nghost = this->m_Opt->m_PDESolver.iporder;
      ierr = AllocateOnce(this->m_GhostPlan, this->m_Opt, this->nghost); CHKERRQ(ierr);
      this->g_alloc_max = this->m_GhostPlan->get_ghost_local_size_x(this->isize_g, this->istart_g);
      this->nlghost  = this->isize_g[0];
      this->nlghost *= this->isize_g[1];
      this->nlghost *= this->isize_g[2];
      
      cudaMalloc((void**)&this->m_VecFieldGhost, this->g_alloc_max); 
      //cudaMalloc((void**)&this->m_VecFieldGhost, 3*this->g_alloc_max); 
      //cudaMalloc((void**)&this->m_ScaFieldGhost, this->g_alloc_max); 
      //cudaMalloc((void**)&this->m_WorkScaField1, 3*nl*sizeof(ScalarType));
      //cudaMalloc((void**)&this->m_WorkScaField2, nl*sizeof(ScalarType));

      ierr = AllocateOnce(this->m_StatePlan, this->g_alloc_max, this->cuda_aware);
      this->m_StatePlan->allocate(nl, this->m_Dofs, 2, this->nghost, this->isize_g);

      //ierr = AllocateOnce(this->m_AdjointPlan, this->g_alloc_max, this->cuda_aware);
      //this->m_AdjointPlan->allocate(nl, this->m_Dofs, 2);
    } else {
      ierr = AllocateOnce(this->m_Xstate, this->m_Opt); CHKERRQ(ierr);
      ierr = AllocateOnce(this->m_Xadjoint, this->m_Opt); CHKERRQ(ierr);
    }
    
    //ierr = AllocateOnce(this->m_InitialTrajectory, m_Opt); CHKERRQ(ierr);
    //ierr = this->ComputeInitialTrajectory(); CHKERRQ(ierr);
    ierr = this->InitializeInterpolationTexture(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief clears memory
 * perform everything on the GPU
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Xstate != nullptr) {
        delete this->m_Xstate;
        this->m_Xstate = nullptr;
    }
    
    if (this->m_Xadjoint != nullptr) {
        delete this->m_Xadjoint;
        this->m_Xadjoint = nullptr;
    }

    if (this->m_InitialTrajectory != nullptr) {
        delete this->m_InitialTrajectory;
        this->m_InitialTrajectory = nullptr;
    }

    if (this->m_X != nullptr) {
        ierr = VecDestroy(&this->m_X); CHKERRQ(ierr);
        this->m_X = nullptr;
    }
    
    if (this->m_WorkScaField1 != nullptr) {
      cudaFree(this->m_WorkScaField1);
      this->m_WorkScaField1 = nullptr;
    }
    
    if (this->m_WorkScaField2 != nullptr) {
      cudaFree(this->m_WorkScaField2);
      this->m_WorkScaField2 = nullptr;
    }
    
    if (this->m_ScaFieldGhost != nullptr) {
      cudaFree(this->m_ScaFieldGhost);
      this->m_ScaFieldGhost = nullptr;
    }

    if (this->m_VecFieldGhost != nullptr) {
      cudaFree(this->m_VecFieldGhost); 
      this->m_VecFieldGhost = nullptr;
    }

    if (this->m_texture != 0) {
      cudaDestroyTextureObject(this->m_texture);
    }
  
    if (this->m_tmpInterpol1 != nullptr) {
      cudaFree(this->m_tmpInterpol1);
      this->m_tmpInterpol1 = nullptr;
    }
    
    if (this->m_tmpInterpol2 != nullptr) {
      cudaFree(this->m_tmpInterpol2);
    }

    Free(this->m_StatePlan);
    Free(this->m_AdjointPlan);
    Free(this->m_GhostPlan);
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief init empty texture for interpolation on GPU 
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::InitializeInterpolationTexture() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Opt->rank_cnt == 1) {
      int isize[3];
      isize[0] = static_cast<int>(this->m_Opt->m_Domain.isize[0]);
      isize[1] = static_cast<int>(this->m_Opt->m_Domain.isize[1]);
      isize[2] = static_cast<int>(this->m_Opt->m_Domain.isize[2]);
      this->m_texture = gpuInitEmptyTexture(isize);
      if (this->m_Opt->m_PDESolver.iporder == 3) {
        cudaMalloc((void**) &this->m_tmpInterpol1, sizeof(float)*this->m_Opt->m_Domain.nl);
        cudaMalloc((void**) &this->m_tmpInterpol2, sizeof(float)*this->m_Opt->m_Domain.nl);
      }
    } else {
      this->m_texture = gpuInitEmptyTexture(this->isize_g);
    }
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set work vector field to not have to allocate it locally
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::SetWorkVecField(VecField* x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(x != nullptr, "null pointer"); CHKERRQ(ierr);
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
    ScalarType *p_x1=nullptr,*p_x2=nullptr,*p_x3=nullptr;
    IntType isize[3],istart[3];
    IntType nx[3];
    
    ScalarType hx[3];
    IntType l,i1,i2,i3;
    PetscFunctionBegin;


    for (unsigned int i = 0; i < 3; ++i){
      if (this->m_Opt->rank_cnt == 1) {
        hx[i]     = 1.;
      } else {
        hx[i]     = 1./static_cast<ScalarType>(this->m_Opt->m_Domain.nx[i]);//this->m_Opt->m_Domain.hx[i];
      }
        nx[i]     = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i]  = this->m_Opt->m_Domain.isize[i];
        istart[i] = this->m_Opt->m_Domain.istart[i];
    }

    ierr = VecGetArray(this->m_InitialTrajectory->m_X1, &p_x1); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_InitialTrajectory->m_X2, &p_x2); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_InitialTrajectory->m_X3, &p_x3); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
    for (int i1 = 0; i1 < isize[0]; ++i1){  // x1
        for (int i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (int i3 = 0; i3 < isize[2]; ++i3){ // x3

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
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief set the initial trajectory from outside the class
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::SetInitialTrajectory(const ScalarType* pX) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    
    ierr = this->m_InitialTrajectory->SetComponents(pX, "stride"); CHKERRQ(ierr);
    if (this->m_Opt->rank_cnt == 1) {
      ierr = VecScale(this->m_InitialTrajectory->m_X1, 1./this->m_Opt->m_Domain.hx[0]); CHKERRQ(ierr);
      ierr = VecScale(this->m_InitialTrajectory->m_X2, 1./this->m_Opt->m_Domain.hx[1]); CHKERRQ(ierr);
      ierr = VecScale(this->m_InitialTrajectory->m_X3, 1./this->m_Opt->m_Domain.hx[2]); CHKERRQ(ierr);
    } else {
      ierr = VecScale(this->m_InitialTrajectory->m_X1, this->m_Opt->m_Domain.nx[0]); CHKERRQ(ierr);
      ierr = VecScale(this->m_InitialTrajectory->m_X2, this->m_Opt->m_Domain.nx[1]); CHKERRQ(ierr);
      ierr = VecScale(this->m_InitialTrajectory->m_X3, this->m_Opt->m_Domain.nx[2]); CHKERRQ(ierr);
    }
    
    if (m_Opt->m_Verbosity > 1) {
      DbgMsgCall("SetInitialTrajectory completed"); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute the trajectory from the velocity field based
 * on an rk2 scheme (todo: make the velocity field a const vector)
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::ComputeTrajectoryRK2(VecField* v, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType ht, hthalf, hx[3], invhx[3], scale = 0.0;
    IntType isize[3], istart[3], nx[3];
    std::stringstream ss;
    ScalarType* xq = nullptr;
    IntType nl=1, ng=1;
    const ScalarType *vx, *vy, *vz;
    double runtime=0;
    const ScalarType* p;
    VecField *X = nullptr;
    ScalarType norm;
    
    TrajectoryKernel kernel;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_WorkVecField1 != nullptr, "null pointer"); CHKERRQ(ierr);
    
    kernel.isize[0] = this->m_Opt->m_Domain.isize[0];
    kernel.isize[1] = this->m_Opt->m_Domain.isize[1];
    kernel.isize[2] = this->m_Opt->m_Domain.isize[2];
    kernel.istart[0] = this->m_Opt->m_Domain.istart[0];
    kernel.istart[1] = this->m_Opt->m_Domain.istart[1];
    kernel.istart[2] = this->m_Opt->m_Domain.istart[2];
    
    if (this->m_Opt->rank_cnt == 1) { // Coordinates in [0, Ni)
      kernel.ix[0] = 1.; kernel.ix[1] = 1.; kernel.ix[2] = 1.;
      kernel.hx[0] = this->m_Opt->GetTimeStepSize()/this->m_Opt->m_Domain.hx[0];
      kernel.hx[1] = this->m_Opt->GetTimeStepSize()/this->m_Opt->m_Domain.hx[1];
      kernel.hx[2] = this->m_Opt->GetTimeStepSize()/this->m_Opt->m_Domain.hx[2];
    } else { // Coordinates in [0, 1)
      kernel.ix[0] = 1./static_cast<ScalarType>(this->m_Opt->m_Domain.nx[0]);
      kernel.ix[1] = 1./static_cast<ScalarType>(this->m_Opt->m_Domain.nx[1]);
      kernel.ix[2] = 1./static_cast<ScalarType>(this->m_Opt->m_Domain.nx[2]);
      kernel.hx[0] = this->m_Opt->GetTimeStepSize()/(2.*PETSC_PI);
      kernel.hx[1] = this->m_Opt->GetTimeStepSize()/(2.*PETSC_PI);
      kernel.hx[2] = this->m_Opt->GetTimeStepSize()/(2.*PETSC_PI);
    }
    
    if (flag.compare("state") == 0) {
        X = this->m_Xstate;
    } else if (flag.compare("adjoint") == 0) {
        X = this->m_Xadjoint;
        kernel.hx[0] *= -1.;
        kernel.hx[1] *= -1.;
        kernel.hx[2] *= -1.;
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }
    
    // RK2 stage 1
    if (this->m_Opt->rank_cnt == 1) {
      ierr = X->GetArraysWrite(kernel.pX); CHKERRQ(ierr);
    } else {
      //kernel.pX[0] = &this->m_VecFieldGhost[0*this->m_Opt->m_Domain.nl];
      //kernel.pX[1] = &this->m_VecFieldGhost[1*this->m_Opt->m_Domain.nl];
      //kernel.pX[2] = &this->m_VecFieldGhost[2*this->m_Opt->m_Domain.nl];
      ierr = this->m_WorkVecField1->GetArrays(kernel.pX); CHKERRQ(ierr);
    }
    ierr = v->GetArraysRead(kernel.pV); CHKERRQ(ierr);
    
    ierr = kernel.RK2_Step1(); CHKERRQ(ierr);
    
    if (this->m_Opt->rank_cnt == 1) {
      ierr = X->RestoreArrays(); CHKERRQ(ierr);
    } else {
      ierr = this->m_WorkVecField1->RestoreArrays(); CHKERRQ(ierr);
      ierr = this->MapCoordinateVector(flag);
    }
    ierr = v->RestoreArrays(); CHKERRQ(ierr);
    
    // Interpolate on Euler coordinates
    ierr = this->Interpolate(this->m_WorkVecField1, v, flag); CHKERRQ(ierr);
    
    // RK2 stage 2 
    ierr = v->GetArraysRead(kernel.pV); CHKERRQ(ierr);
    //ierr = this->m_WorkVecField1->GetArraysRead(kernel.pVx); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(kernel.pVx); CHKERRQ(ierr);

    if (this->m_Opt->rank_cnt == 1) {
      ierr = X->GetArraysWrite(kernel.pX); CHKERRQ(ierr);
    } 
    
    //ierr = kernel.RK2_Step2(); CHKERRQ(ierr);
    ierr = kernel.RK2_Step2_inplace(); CHKERRQ(ierr);
    
    ierr = this->m_WorkVecField1->RestoreArrays(); CHKERRQ(ierr);
    ierr = v->RestoreArrays(); CHKERRQ(ierr);
    
    if (this->m_Opt->rank_cnt > 1) {
      ierr = this->MapCoordinateVector(flag);
    } else {
      ierr = X->RestoreArrays(); CHKERRQ(ierr);
    }
    
    
    /*ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;
    
    if (this->m_Opt->m_Verbosity > 2) {
        std::string str = "update trajectory: ";
        str += flag;
        ierr = DbgMsg2(str); CHKERRQ(ierr);
        ierr = v->DebugInfo("SL v", __LINE__, __FILE__); CHKERRQ(ierr);
    }
    
    // switch between state and adjoint variable
    if (flag.compare("state") == 0) {
        X = this->m_Xstate;
        scale =  1.0;
    } else if (flag.compare("adjoint") == 0) {
        X = this->m_Xadjoint;
        scale = -1.0;
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }


    for (int i = 0; i < 3; ++i) {
        hx[i]     = this->m_Opt->m_Domain.hx[i];
        nx[i]     = this->m_Opt->m_Domain.nx[i];
        isize[i]  = this->m_Opt->m_Domain.isize[i];
        istart[i] = this->m_Opt->m_Domain.istart[i];
        if (this->m_Opt->rank_cnt == 1) {
          invhx[i]  = 1./hx[i];
        } else {
          invhx[i]  = 1./(2.0*PETSC_PI);
        }
        nl       *= isize[i];
        ng       *= nx[i];
    }
  
    ierr = Assert(X != nullptr, "null ptr"); CHKERRQ(ierr);
    ierr = Assert(this->m_Xstate != nullptr, "null ptr"); CHKERRQ(ierr);
    ierr = Assert(this->m_Xadjoint != nullptr, "null ptr"); CHKERRQ(ierr);
    ierr = Assert(v != nullptr, "nullptr"); CHKERRQ(ierr);
    ierr = Assert(this->m_InitialTrajectory != nullptr, "nullptr"); CHKERRQ(ierr);
    
    // X = x - ht*v ; single GPU: X \in [0, N], multi GPU: X \in [0, 1]
    ierr = VecWAXPY(X->m_X1, -scale*ht*invhx[0], v->m_X1, this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X2, -scale*ht*invhx[1], v->m_X2, this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X3, -scale*ht*invhx[2], v->m_X3, this->m_InitialTrajectory->m_X3); CHKERRQ(ierr); 

    if (this->m_Opt->rank_cnt > 1) {
      ierr = this->MapCoordinateVector(flag);
    }
    
    ierr = this->Interpolate(this->m_WorkVecField1, v, flag); CHKERRQ(ierr);

    // X = x - 0.5*ht*(v + v(x - ht v))
    // F = F0 + F1 = v + v(x-ht*v)
    ierr = VecAXPY(this->m_WorkVecField1->m_X1, 1.0, v->m_X1); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_WorkVecField1->m_X2, 1.0, v->m_X2); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_WorkVecField1->m_X3, 1.0, v->m_X3); CHKERRQ(ierr);

    // X = x - 0.5*ht*F
    ierr = VecWAXPY(X->m_X1, -scale*hthalf*invhx[0], this->m_WorkVecField1->m_X1, this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X2, -scale*hthalf*invhx[1], this->m_WorkVecField1->m_X2, this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X3, -scale*hthalf*invhx[2], this->m_WorkVecField1->m_X3, this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);

    if (this->m_Opt->rank_cnt > 1) {
      ierr = this->MapCoordinateVector(flag);
    }*/
    
    if (this->m_Opt->m_Verbosity > 2) {
      DbgMsgCall("Trajectory computed");
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


PetscErrorCode SemiLagrangianGPUNew::ComputeTrajectoryRK4(VecField* v, std::string flag) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief interpolate scalar field
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Interpolate(Vec* xo, Vec xi, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType *p_xo = nullptr, *p_xi = nullptr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(*xo != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert( xi != nullptr, "null pointer"); CHKERRQ(ierr);
    
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
    int nx[3], isize_g[3], isize[3], istart_g[3], istart[3], c_dims[2], neval;
    IntType nl, nalloc;
    std::stringstream ss;
    double timers[4] = {0, 0, 0, 0};
    const ScalarType *xq1, *xq2, *xq3;
    Interp3_Plan_GPU* interp_plan = nullptr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(xi != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xo != nullptr, "null pointer"); CHKERRQ(ierr);
    
    nl     = this->m_Opt->m_Domain.nl;
    neval  = static_cast<int>(nl);

    for (int i = 0; i < 3; ++i) {
        nx[i]     = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i]  = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }
    
    c_dims[0] = this->m_Opt->m_CartGridDims[0];
    c_dims[1] = this->m_Opt->m_CartGridDims[1];

    ZeitGeist_define(SL_INTERPOL);
    ZeitGeist_tick(SL_INTERPOL);
    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    if (this->m_Opt->rank_cnt > 1) {
      ZeitGeist_define(INTERPOL_COMM);
      ZeitGeist_tick(INTERPOL_COMM);
      this->m_GhostPlan->share_ghost_x(xi, this->m_VecFieldGhost);
      ZeitGeist_tock(INTERPOL_COMM);
     /* 
      if (flag.compare("state") == 0) {
          interp_plan = this->m_StatePlan;
      } else if (flag.compare("adjoint") == 0) {
          interp_plan = this->m_AdjointPlan;
      } else {
          ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
      }
      */
      ScalarType *wout[1] = {xo};
      
      interp_plan = this->m_StatePlan;

      interp_plan->interpolate( this->m_VecFieldGhost, 
                                this->isize_g, 
                                this->nlghost,
                                nl, 
                                wout,
                                this->m_Opt->m_Domain.mpicomm, 
                                this->m_tmpInterpol1, 
                                this->m_tmpInterpol2, 
                                this->m_texture, 
                                this->m_Opt->m_PDESolver.iporder, 
                                &(this->m_Opt->m_GPUtime), 0, flag);
    } else {
      // compute interpolation for all components of the input scalar field
      if (flag.compare("state") == 0) {
          ierr = Assert(this->m_Xstate != nullptr, "null pointer"); CHKERRQ(ierr);
          ierr = this->m_Xstate->GetArraysRead(xq1, xq2, xq3);
          const ScalarType* xq[3] = {xq1, xq2, xq3};
          gpuInterp3D(xi, xq, xo, this->m_tmpInterpol1, this->m_tmpInterpol2, nx, static_cast<long int>(nl), this->m_texture, this->m_Opt->m_PDESolver.iporder, &(this->m_Opt->m_GPUtime));
          ierr = this->m_Xstate->RestoreArrays(); CHKERRQ(ierr);
      } else if (flag.compare("adjoint") == 0) {
          ierr = Assert(this->m_Xadjoint != nullptr, "null pointer"); CHKERRQ(ierr);
          ierr = this->m_Xadjoint->GetArraysRead(xq1, xq2, xq3);
          const ScalarType* xq[3] = {xq1, xq2, xq3};
          gpuInterp3D(xi, xq, xo, this->m_tmpInterpol1, this->m_tmpInterpol2, nx, static_cast<long int>(nl), this->m_texture, this->m_Opt->m_PDESolver.iporder, &(this->m_Opt->m_GPUtime));
          ierr = this->m_Xadjoint->RestoreArrays(); CHKERRQ(ierr);
      } else {
          ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
      }
    }

    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);
    ZeitGeist_tock(SL_INTERPOL);
    
    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IP);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief interpolate vector field
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Interpolate(VecField* vo, VecField* vi, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType *p_vix1 = nullptr, *p_vix2 = nullptr, *p_vix3 = nullptr;
    ScalarType *p_vox1 = nullptr, *p_vox2 = nullptr, *p_vox3 = nullptr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vi != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vo != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = vi->GetArraysReadWrite(p_vix1, p_vix2, p_vix3); CHKERRQ(ierr);
    ierr = vo->GetArraysReadWrite(p_vox1, p_vox2, p_vox3); CHKERRQ(ierr);
    ierr = this->Interpolate(p_vox1, p_vox2, p_vox3, p_vix1, p_vix2, p_vix3, flag); CHKERRQ(ierr);
    ierr = vi->RestoreArrays(); CHKERRQ(ierr);
    ierr = vo->RestoreArrays(); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief interpolate vector field - single GPU optimised version
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Interpolate(ScalarType* wx1, ScalarType* wx2, ScalarType* wx3,
                                                 ScalarType* vx1, ScalarType* vx2, ScalarType* vx3, std::string flag) {
    PetscErrorCode ierr = 0;
    int nx[3], isize_g[3], isize[3], istart_g[3], istart[3], c_dims[2], nghost, order;
    double timers[4] = {0, 0, 0, 0};
    std::stringstream ss;
    IntType nl, nlghost, g_alloc_max;
    Interp3_Plan_GPU* interp_plan = nullptr;
    const ScalarType *xq1, *xq2, *xq3;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vx1 != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vx2 != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vx3 != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx1 != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx2 != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx3 != nullptr, "null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->m_Domain.nl;
    order = this->m_Opt->m_PDESolver.iporder;
    nghost = order;

    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i] = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }
    
    c_dims[0] = this->m_Opt->m_CartGridDims[0];
    c_dims[1] = this->m_Opt->m_CartGridDims[1];

    ZeitGeist_define(SL_INTERPOL);
    ZeitGeist_tick(SL_INTERPOL);
    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    interp_plan = this->m_StatePlan;
    
    ScalarType* vin[3] = {vx1, vx2, vx3};
    ScalarType* wout[3] = {wx1, wx2, wx3};

    if (this->m_Opt->rank_cnt > 1) {
      ZeitGeist_define(INTERPOL_COMM);
      for (int i=0; i<3; i++) { 
        ZeitGeist_tick(INTERPOL_COMM);
        this->m_GhostPlan->share_ghost_x(vin[i], &this->m_VecFieldGhost[0*this->nlghost]);
        //this->m_GhostPlan->share_ghost_x(vx1, &this->m_VecFieldGhost[0*this->nlghost]);
        //this->m_GhostPlan->share_ghost_x(vx2, &this->m_VecFieldGhost[1*this->nlghost]);
        //this->m_GhostPlan->share_ghost_x(vx3, &this->m_VecFieldGhost[2*this->nlghost]);
        ZeitGeist_tock(INTERPOL_COMM);

        //ScalarType *wout[3] = {wx1, wx2, wx3};
        
        // do interpolation
        interp_plan->interpolate( this->m_VecFieldGhost, 
                                  this->isize_g, 
                                  this->nlghost,
                                  nl, 
                                  &wout[i],
                                  this->m_Opt->m_Domain.mpicomm, 
                                  this->m_tmpInterpol1, 
                                  this->m_tmpInterpol2, 
                                  this->m_texture, 
                                  this->m_Opt->m_PDESolver.iporder, 
                                  &(this->m_Opt->m_GPUtime), 0, flag);

        //ierr = cudaMemcpy((void*)wx1, (const void*)&this->m_WorkScaField1[0*nl], nl*sizeof(ScalarType), cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
        //ierr = cudaMemcpy((void*)wx2, (const void*)&this->m_WorkScaField1[1*nl], nl*sizeof(ScalarType), cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
        //ierr = cudaMemcpy((void*)wx3, (const void*)&this->m_WorkScaField1[2*nl], nl*sizeof(ScalarType), cudaMemcpyDeviceToDevice); CHKERRCUDA(ierr);
      }

    } else {
      

      if (flag.compare("state") == 0) {

          ierr = this->m_Xstate->GetArraysRead(xq1, xq2, xq3);
          const ScalarType* xq[3] = {xq1, xq2, xq3};
          gpuInterpVec3D(vx1, vx2, vx3, xq, wx1, wx2, wx3, this->m_tmpInterpol1, this->m_tmpInterpol2, nx, static_cast<long int>(nl), this->m_texture, this->m_Opt->m_PDESolver.iporder, &(this->m_Opt->m_GPUtime));
          ierr = this->m_Xstate->RestoreArrays(); CHKERRQ(ierr);

      } else if (flag.compare("adjoint") == 0) {
          
          ierr = this->m_Xadjoint->GetArraysRead(xq1, xq2, xq3);
          const ScalarType* xq[3] = {xq1, xq2, xq3};
          gpuInterpVec3D(vx1, vx2, vx3, xq, wx1, wx2, wx3, this->m_tmpInterpol1, this->m_tmpInterpol2, nx, static_cast<long int>(nl), this->m_texture, this->m_Opt->m_PDESolver.iporder, &(this->m_Opt->m_GPUtime));
          ierr = this->m_Xadjoint->RestoreArrays(); CHKERRQ(ierr);
         
      } else {
          ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
      }
    }
    
    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);
    ZeitGeist_tock(SL_INTERPOL);
    ZeitGeist_inc(SL_INTERPOL);
    ZeitGeist_inc(SL_INTERPOL);

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IPVEC);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief set coordinate vector and communicate to interpolation plan
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::SetQueryPoints(ScalarType* y1, ScalarType* y2, ScalarType* y3, std::string flag) {
    PetscErrorCode ierr = 0;
    IntType nl;
    VecField* X = nullptr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->m_Domain.nl;

    if (flag.compare("state") == 0) {
        X = this->m_Xstate;
    } else if (flag.compare("adjoint") == 0) {
        X = this->m_Xadjoint;
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }
    
    ierr = X->SetComponents(y1, y2, y3); 
    
    ScalarType invhx[3];
    for (int i=0; i<3; i++) invhx[i] = 1./this->m_Opt->m_Domain.hx[i];
    ierr = X->Scale(invhx);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

PetscErrorCode SemiLagrangianGPUNew::GetQueryPoints(ScalarType* y1, ScalarType* y2, ScalarType* y3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_Xstate != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = this->m_Xstate->GetComponents(y1, y2, y3); CHKERRQ(ierr);
#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
#pragma omp parallel for
    for (IntType i = 0; i < this->m_Opt->m_Domain.nl; ++i) {
      y1[i] *= this->m_Opt->m_Domain.hx[0];
      y2[i] *= this->m_Opt->m_Domain.hx[1];
      y3[i] *= this->m_Opt->m_Domain.hx[2];
    }
#endif
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * Name: GetQueryPoints
 * Description: get the query points in the pointer y
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::GetQueryPoints(ScalarType* y) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_Xstate != nullptr, "null pointer"); CHKERRQ(ierr);

    ierr = this->m_Xstate->GetComponents(y, "stride"); CHKERRQ(ierr);
#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
#pragma omp parallel for
    for (IntType i = 0; i < this->m_Opt->m_Domain.nl; ++i) {
      y[3*i+0] *= this->m_Opt->m_Domain.hx[0];
      y[3*i+1] *= this->m_Opt->m_Domain.hx[1];
      y[3*i+2] *= this->m_Opt->m_Domain.hx[2];
    }
#endif
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


PetscErrorCode SemiLagrangianGPUNew::CommunicateCoord(std::string flag) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  PetscFunctionReturn(ierr);
}

/********************************************************************
 * Name: map coordinates
 * Description: change from lexicographical ordering to xyz
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::MapCoordinateVector(std::string flag) {
    PetscErrorCode ierr;
    int nx[3], isize_g[3], istart_g[3], c_dims[2], isize[3], istart[3];
    IntType nl;
    double timers[4] = {0,0,0,0};
    ScalarType* p_X[3] = {nullptr, nullptr, nullptr};

    PetscFunctionBegin;
    
    if (this->m_Opt->m_Verbosity > 2) {
      ierr = DbgMsgCall("Mapping query points"); CHKERRQ(ierr);
    }

    for (int i = 0; i < 3; ++i){
        nx[i] = this->m_Opt->m_Domain.nx[i];
        isize[i] = this->m_Opt->m_Domain.isize[i];
        istart[i] = this->m_Opt->m_Domain.istart[i];
    }
    
    c_dims[0] = this->m_Opt->m_CartGridDims[0];
    c_dims[1] = this->m_Opt->m_CartGridDims[1];
  
    nl = this->m_Opt->m_Domain.nl;
    
    ZeitGeist_define(INTERPOL_COMM);
    ZeitGeist_tick(INTERPOL_COMM);
    
    //p_X[0] = &this->m_VecFieldGhost[0*this->m_Opt->m_Domain.nl];
    //p_X[1] = &this->m_VecFieldGhost[1*this->m_Opt->m_Domain.nl];
    //p_X[2] = &this->m_VecFieldGhost[2*this->m_Opt->m_Domain.nl];
    ierr = this->m_WorkVecField1->GetArrays(p_X); CHKERRQ(ierr);
    
    this->m_StatePlan->scatter(nx, isize, istart, nl, this->nghost, p_X[0], p_X[1], p_X[2], c_dims, this->m_Opt->m_Domain.mpicomm, timers, flag);

    ierr = this->m_WorkVecField1->RestoreArrays(); CHKERRQ(ierr);

/*
    if (flag.compare("state")==0) {
        
        //ierr = this->m_Xstate->GetArrays(p_X);
        this->m_StatePlan->scatter(nx, isize, istart, nl, this->nghost, p_X[0], p_X[1], p_X[2], c_dims, this->m_Opt->m_Domain.mpicomm, timers);
        //ierr = this->m_Xstate->RestoreArrays(p_X);
      
    } else if (flag.compare("adjoint")==0) {

        //ierr = this->m_Xadjoint->GetArrays(p_X);
        this->m_AdjointPlan->scatter(nx, isize, istart, nl, this->nghost, p_X[0], p_X[1], p_X[2], c_dims, this->m_Opt->m_Domain.mpicomm, timers);
        //ierr = this->m_Xadjoint->RestoreArrays(p_X);

    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }
*/    
    ZeitGeist_tock(INTERPOL_COMM);

    this->m_Opt->IncreaseInterpTimers(timers);

    PetscFunctionReturn(ierr);
}

}  // namespace reg




#endif  // _SEMILAGRANGIAN_CPP_
