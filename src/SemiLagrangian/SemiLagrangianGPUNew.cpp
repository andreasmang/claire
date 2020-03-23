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
    this->m_Opt = opt;
    this->Initialize();
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
    int isize[3];
    PetscFunctionBegin;
    
    this->m_Xstate = nullptr;
    this->m_Xadjoint = nullptr;
    this->m_X = nullptr;
    this->m_WorkVec = nullptr;
    this->m_ScaFieldGhost = nullptr;
    this->m_VecFieldGhost = nullptr;
    this->m_WorkVecField1 = nullptr;
    this->m_InitialTrajectory = nullptr;
    this->m_texture = 0;

    this->m_Dofs[0] = 1;
    this->m_Dofs[1] = 3;
    
    this->m_tmpInterpol1 = nullptr;
    this->m_tmpInterpol2 = nullptr;
    

    this->m_GhostWork1 = nullptr;
    this->m_GhostWork2 = nullptr;
  
    this->m_StatePlan = nullptr;
    this->m_StatePlanVec = nullptr;
    this->m_AdjointPlan = nullptr;
    this->m_AdjointPlanVec = nullptr;
  
    ierr = AllocateOnce(this->m_Xstate, this->m_Opt); CHKERRQ(ierr);
    ierr = AllocateOnce(this->m_Xadjoint, this->m_Opt); CHKERRQ(ierr);

#if defined(REG_HAS_MPICUDA)
    int nl = this->m_Opt->m_Domain.nl;
    int ng = this->m_Opt->m_Domain.ng;
    for (int i=0; i<3; i++) {
      isize[i] = this->m_Opt->m_Domain.isize[i];
    }
    
    this->nghost = this->m_Opt->m_PDESolver.iporder;
    this->g_alloc_max = accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->m_FFT.fft->m_plan, this->nghost, this->isize_g, this->istart_g);
    this->nlghost = this->isize_g[0]*this->isize_g[1]*this->isize_g[2];
    this->m_GhostWork1 = pvfmm::aligned_new<ScalarType> (this->m_Opt->m_FFT.fft->m_plan->alloc_max + 2 * this->nghost * isize[2] * isize[0]);
    this->m_GhostWork2 = pvfmm::aligned_new<ScalarType> (this->m_Opt->m_FFT.fft->m_plan->alloc_max + 2 * this->nghost * isize[2] * isize[0] + 2 * this->nghost * isize[2] * this->isize_g[1]);
    this->m_VecFieldGhost = reinterpret_cast<ScalarType*> (accfft_alloc(3*this->g_alloc_max));
    
    ierr = AllocateOnce(this->m_StatePlan, this->m_Dofs[0]*this->g_alloc_max);
    this->m_StatePlan->allocate(nl, this->m_Dofs[0]);

    ierr = AllocateOnce(this->m_StatePlanVec, this->m_Dofs[1]*this->g_alloc_max);
    this->m_StatePlanVec->allocate(nl, this->m_Dofs[1]); 

    ierr = AllocateOnce(this->m_AdjointPlan, this->m_Dofs[0]*this->g_alloc_max);
    this->m_AdjointPlan->allocate(nl, this->m_Dofs[0]);

    ierr = AllocateOnce(this->m_AdjointPlanVec, this->m_Dofs[1]*this->g_alloc_max);
    this->m_AdjointPlanVec->allocate(nl, this->m_Dofs[1]);
    
    ierr = VecCreate(this->m_X, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecCreate(this->m_WorkVec, 3*nl, 3*ng); CHKERRQ(ierr);
#endif
    ierr = AllocateOnce(this->m_InitialTrajectory, m_Opt); CHKERRQ(ierr);
    ierr = this->ComputeInitialTrajectory(); CHKERRQ(ierr);
    ierr = this->InitializeInterpolationTexture(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief init empty texture for interpolation on GPU 
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::InitializeInterpolationTexture() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
    this->m_texture = gpuInitEmptyTexture(this->m_Opt->m_Domain.isize);
    if (this->m_Opt->m_PDESolver.iporder == 3) {
      cudaMalloc((void**) &this->m_tmpInterpol1, sizeof(float)*this->m_Opt->m_Domain.nl);
      cudaMalloc((void**) &this->m_tmpInterpol2, sizeof(float)*this->m_Opt->m_Domain.nl);
    }
#elif defined(REG_HAS_MPICUDA) 
    this->m_texture = gpuInitEmptyTexture(this->isize_g);
#else
    ierr = DebugNotImplemented(); CHKERRQ(ierr);
#endif
    
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
    
    if (this->m_WorkVec != nullptr) {
        ierr = VecDestroy(&this->m_WorkVec); CHKERRQ(ierr); 
        this->m_WorkVec = nullptr;
    }
  
    if (this->m_ScaFieldGhost != nullptr) {
        accfft_free(this->m_ScaFieldGhost);
        this->m_ScaFieldGhost = nullptr;
    }

    if (this->m_VecFieldGhost != nullptr) {
        accfft_free(this->m_VecFieldGhost);
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
      this->m_tmpInterpol2 = nullptr;
    }

    if (this->m_GhostWork1 != nullptr) {
      pvfmm::aligned_delete<ScalarType>(this->m_GhostWork1);
      this->m_GhostWork1 = nullptr;
    }

    if (this->m_GhostWork2 != nullptr) {
      pvfmm::aligned_delete<ScalarType>(this->m_GhostWork2);
      this->m_GhostWork2 = nullptr;
    }

    if (this->m_StatePlan != nullptr) {
      Free(this->m_StatePlan);
      this->m_StatePlan = nullptr;
    }

    if (this->m_StatePlanVec != nullptr) {
      Free(this->m_StatePlanVec);
      this->m_StatePlanVec = nullptr;
    }

    if (this->m_AdjointPlan != nullptr) {
      Free(this->m_AdjointPlan);
      this->m_AdjointPlan = nullptr;
    }

    if (this->m_AdjointPlanVec != nullptr) {
      Free(this->m_AdjointPlanVec);
      this->m_AdjointPlanVec = nullptr;
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
    int nx[3];
    
    ScalarType hx[3];
    IntType l,i1,i2,i3;
    PetscFunctionBegin;


    for (unsigned int i = 0; i < 3; ++i){
        hx[i]     = this->m_Opt->m_Domain.hx[i];
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
    
    PetscFunctionReturn(0);
}

/********************************************************************
 * @brief set the initial trajectory from outside the class
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::SetInitialTrajectory(const ScalarType* pX) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    
    ierr = this->m_InitialTrajectory->SetComponents(pX, "stride"); CHKERRQ(ierr);
#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA) 
    ierr = VecScale(this->m_InitialTrajectory->m_X1, 1./this->m_Opt->m_Domain.hx[0]); CHKERRQ(ierr);
    ierr = VecScale(this->m_InitialTrajectory->m_X2, 1./this->m_Opt->m_Domain.hx[1]); CHKERRQ(ierr);
    ierr = VecScale(this->m_InitialTrajectory->m_X3, 1./this->m_Opt->m_Domain.hx[2]); CHKERRQ(ierr);
#endif
    PetscFunctionReturn(0);
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
    VecField *X;
    ScalarType norm;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_WorkVecField1 != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;
    
    if (this->m_Opt->m_Verbosity > 2) {
        std::string str = "update trajectory: ";
        str += flag;
        ierr = DbgMsg2(str); CHKERRQ(ierr);
        ierr = v->DebugInfo("SL v", __LINE__, __FILE__); CHKERRQ(ierr);
    }

    // switch between state and adjoint variable
    if (strcmp(flag.c_str(), "state") == 0) {
        X = this->m_Xstate;
        scale =  1.0;
    } else if (strcmp(flag.c_str(), "adjoint") == 0) {
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
        invhx[i]  = static_cast<ScalarType>(this->m_Opt->m_Domain.nx[i])/(PETSC_PI*2.0);
        nl       *= isize[i];
        ng       *= nx[i];
    }
    

#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
    // X = x - (ht/hx) v
    ierr = VecWAXPY(X->m_X1, -scale*ht*invhx[0], v->m_X1, this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X2, -scale*ht*invhx[1], v->m_X2, this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X3, -scale*ht*invhx[2], v->m_X3, this->m_InitialTrajectory->m_X3); CHKERRQ(ierr); 
#else
    // X = x - ht*v
    // do not use invhx when scattering (will mess up with Amir's scattering function)
    ierr = VecWAXPY(X->m_X1, -scale*ht, v->m_X1, this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X2, -scale*ht, v->m_X2, this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X3, -scale*ht, v->m_X3, this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);
  
   // ierr = VecGetArrayRead(this->m_InitialTrajectory->m_X1, &p);
   // ierr = VecGetArrayRead(v->m_X1, &vx); 
   // for (int i = 0; i < nl; i++ ) { 
   //     printf("%d, x = %f \t v = %f \t X = %f\n", i, p[i], vx[i], p[i]-ht*vx[i]);
   // }
   // ierr = VecRestoreArrayRead(this->m_InitialTrajectory->m_X1, &p);
   // ierr = VecRestoreArrayRead(v->m_X1, &vx); 

    ierr = X->Scale(1.0f/(2.0f*PETSC_PI)); CHKERRQ(ierr);
    
    // flatten out the coordinate array
    ierr = X->GetComponents(this->m_X, "stride"); CHKERRQ(ierr); 
    
    //const ScalarType* p_mX;
    //ierr = VecGetArrayRead(this->m_X, &p_mX); CHKERRQ(ierr);
    //for (int i=0; i<3*nl; i+=3) { 
    //    PetscPrintf(PETSC_COMM_WORLD, "[%d]  X*[0] = %f, X*[1] = %f, X*[2] = %f\n", i/3, p_mX[i], p_mX[i+1], p_mX[i+2]);
    //}
    //ierr = VecRestoreArrayRead(this->m_X, &p_mX); CHKERRQ(ierr);

    // need to communicate the coordinates here before interpolation
    ierr = this->MapCoordinateVector(flag);
#endif
    // interpolate velocity field v(X)
    ierr = this->Interpolate(this->m_WorkVecField1, v, flag); CHKERRQ(ierr);
    
    //const ScalarType *p_mw, *p_v;
    //ierr = VecGetArrayRead(v->m_X1, &p_v); CHKERRQ(ierr);
    //ierr = VecGetArrayRead(this->m_WorkVecField1->m_X1, &p_mw); CHKERRQ(ierr);
    //for (int i=0; i<nl; i++) { 
    //    PetscPrintf(PETSC_COMM_WORLD, "[%d]  v = %f, mw = %f\n", i, p_v[i], p_mw[i]);
    //}
    //ierr = VecRestoreArrayRead(v->m_X1, &p_v); CHKERRQ(ierr);
    //ierr = VecRestoreArrayRead(this->m_WorkVecField1->m_X1, &p_mw); CHKERRQ(ierr);

    //ierr = VecNorm(this->m_WorkVecField1->m_X1, NORM_2, &norm); CHKERRQ(ierr);
    //PetscPrintf(PETSC_COMM_WORLD, "norm of interpolated work field = %f\n", norm);

    //ierr = VecNorm(this->m_WorkVecField1->m_X1, NORM_1, &norm); CHKERRQ(ierr);
    //std::cout << "v norm at euler point = " << norm << std::endl;

    // X = x - 0.5*ht*(v + v(x - ht v))
    // F = F0 + F1 = v + v(x-ht*v)
    ierr = VecAXPY(this->m_WorkVecField1->m_X1, 1.0, v->m_X1); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_WorkVecField1->m_X2, 1.0, v->m_X2); CHKERRQ(ierr);
    ierr = VecAXPY(this->m_WorkVecField1->m_X3, 1.0, v->m_X3); CHKERRQ(ierr);

    // X = x - 0.5*ht*F
#if defined(REG_HAS_CUDA) && !defined(REG_HAS_MPICUDA)
    ierr = VecWAXPY(X->m_X1, -scale*hthalf*invhx[0], this->m_WorkVecField1->m_X1, this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X2, -scale*hthalf*invhx[1], this->m_WorkVecField1->m_X2, this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X3, -scale*hthalf*invhx[2], this->m_WorkVecField1->m_X3, this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);
#else
    // no invhx factor when using MPI CUDA or MPI in general
    ierr = VecWAXPY(X->m_X1, -scale*hthalf, this->m_WorkVecField1->m_X1, this->m_InitialTrajectory->m_X1); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X2, -scale*hthalf, this->m_WorkVecField1->m_X2, this->m_InitialTrajectory->m_X2); CHKERRQ(ierr);
    ierr = VecWAXPY(X->m_X3, -scale*hthalf, this->m_WorkVecField1->m_X3, this->m_InitialTrajectory->m_X3); CHKERRQ(ierr);
    
    ierr = X->Scale(1.0f/(2.0f*PETSC_PI)); CHKERRQ(ierr);
    
    // flatten out the coordinate array
    ierr = X->GetComponents(this->m_X, "stride"); CHKERRQ(ierr); 
    // need to communicate the coordinates here before interpolation for the state or adjoint PDE
    ierr = this->MapCoordinateVector(flag);
#endif
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
    int nx[3], isize_g[3], isize[3], istart_g[3], istart[3], c_dims[2], neval, order, nghost;
    IntType nl, nalloc;
    std::stringstream ss;
    double timers[4] = {0, 0, 0, 0};
    const ScalarType *xq1, *xq2, *xq3;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(xi != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xo != nullptr, "null pointer"); CHKERRQ(ierr);

    ZeitGeist_define(SL_INTERPOL);
    ZeitGeist_tick(SL_INTERPOL);
    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    nl     = this->m_Opt->m_Domain.nl;
    nghost = this->m_Opt->m_PDESolver.iporder;
    neval  = static_cast<int>(nl);

    for (int i = 0; i < 3; ++i) {
        nx[i]     = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i]  = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }
    
     
    // compute interpolation for all components of the input scalar field
    if (strcmp(flag.c_str(), "state") == 0) {
        ierr = Assert(this->m_Xstate != nullptr, "null pointer"); CHKERRQ(ierr);
        ierr = this->m_Xstate->GetArraysRead(xq1, xq2, xq3);
        //gpuInterp3D(xi, xq1, xq2, xq3, xo, this->m_tmpInterpol1, this->m_tmpInterpol2, nx, this->m_texture, this->m_Opt->m_PDESolver.iporder, &(this->m_Opt->m_GPUtime));
        ierr = this->m_Xstate->RestoreArrays(); CHKERRQ(ierr);
    } else if (strcmp(flag.c_str(), "adjoint") == 0) {
        ierr = Assert(this->m_Xadjoint != nullptr, "null pointer"); CHKERRQ(ierr);
        ierr = this->m_Xadjoint->GetArraysRead(xq1, xq2, xq3);
        //gpuInterp3D(xi, xq1, xq2, xq3, xo, this->m_tmpInterpol1, this->m_tmpInterpol2, nx, this->m_texture, this->m_Opt->m_PDESolver.iporder, &(this->m_Opt->m_GPUtime));
        ierr = this->m_Xadjoint->RestoreArrays(); CHKERRQ(ierr);
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }
        
    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);
    ZeitGeist_tock(SL_INTERPOL);
    
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
    ScalarType *p_vix1 = nullptr, *p_vix2 = nullptr, *p_vix3 = nullptr;
    ScalarType *p_vox1 = nullptr, *p_vox2 = nullptr, *p_vox3 = nullptr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vi != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vo != nullptr, "null pointer"); CHKERRQ(ierr);
    
#ifdef REG_HAS_MPICUDA
    ierr = vi->GetComponents(this->m_WorkVec, "block"); CHKERRQ(ierr);
    ierr = this->Interpolate(flag); CHKERRQ(ierr);
    ierr = vo->SetComponents(this->m_WorkVec, "block"); CHKERRQ(ierr);
#else
    ierr = vi->GetArraysReadWrite(p_vix1, p_vix2, p_vix3); CHKERRQ(ierr);
    ierr = vo->GetArraysReadWrite(p_vox1, p_vox2, p_vox3); CHKERRQ(ierr);
    ierr = this->Interpolate(p_vox1, p_vox2, p_vox3, p_vix1, p_vix2, p_vix3, flag); CHKERRQ(ierr);
    ierr = vi->RestoreArraysReadWrite(p_vix1, p_vix2, p_vix3); CHKERRQ(ierr);
    ierr = vo->RestoreArraysReadWrite(p_vox1, p_vox2, p_vox3); CHKERRQ(ierr);
#endif
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief interpolate vector field - multi-GPU version
 *******************************************************************/
PetscErrorCode SemiLagrangianGPUNew::Interpolate(std::string flag) {
    PetscFunctionBegin;
    PetscErrorCode ierr = 0;
    
    int nx[3], isize_g[3], isize[3], istart_g[3], istart[3], c_dims[2], nghost, order;
    double timers[4] = {0, 0, 0, 0};
    std::stringstream ss;
    IntType nl, nlghost, g_alloc_max;
    Interp3_Plan_GPU* interp_plan = nullptr;
    ScalarType* p_WorkVec = nullptr;
    IntType i;

    this->m_Opt->Enter(__func__);
    
    nl = this->m_Opt->m_Domain.nl;
    order = this->m_Opt->m_PDESolver.iporder;
    

    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i] = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }

    c_dims[0] = this->m_Opt->m_CartGridDims[0];
    c_dims[1] = this->m_Opt->m_CartGridDims[1];

    ierr = VecGetArray(this->m_WorkVec, &p_WorkVec); CHKERRQ(ierr);
    for (unsigned int i = 0; i < 3; i++) {
      share_ghost_layer(this->m_Opt->m_FFT.fft->m_plan, this->nghost, this->isize_g, &p_WorkVec[i*nl], &this->m_VecFieldGhost[i*this->nlghost], this->m_GhostWork1, this->m_GhostWork2); 
    }
    ierr = VecRestoreArray(this->m_WorkVec, &p_WorkVec); CHKERRQ(ierr);
    
    //ierr = VecGetArray(this->m_WorkVec, &p_WorkVec);
    //// Small check for correctness of ghost vecfield
    //for (IntType k=0; k<3; k++) {
    //for (IntType i1 = 0; i1 < isize[0]; i1++) {
    //    //printf("\n\n");
    //    for (IntType i2 = 0; i2 < isize[1]; i2++) {
    //        //printf("\n");
    //        for (IntType i3 = 0; i3 < isize[2]; i3++) {
    //            i = reg::GetLinearIndex(i1, i2, i3, isize);
    //            printf("%f,", k*nl+i,p_WorkVec[k*nl+i]);
    //        }
    //        printf("\n");
    //    }
    //    printf("\n\n");
    //}
    //printf("\n\n\n");
    //}
    //ierr = VecRestoreArray(this->m_WorkVec, &p_WorkVec);

    
    // Small check for correctness of ghost vecfield - THIS WORKS
    //std::cout << "=================================================================================================" << std::endl;
    //for (IntType k=2; k<3; k++) {
    //    for (IntType i1 = 0; i1 < this->isize_g[0]; i1++) {
    //        //printf("\n\n");
    //        for (IntType i2 = 0; i2 < this->isize_g[1]; i2++) {
    //            //printf("\n");
    //            for (IntType i3 = 0; i3 < this->isize_g[2]; i3++) {
    //                int idx = reg::GetLinearIndex(i1, i2, i3, this->isize_g);
    //                printf("%f,", this->m_VecFieldGhost[k*this->nlghost+idx]);
    //            }
    //            printf("\n");
    //        }
    //        printf("\n\n");
    //    }
    //    printf("\n\n\n\n");
    //}
    
    ZeitGeist_define(SL_INTERPOL);
    ZeitGeist_tick(SL_INTERPOL);
    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    if (strcmp(flag.c_str(),"state") == 0) {
        interp_plan = this->m_StatePlanVec;
    } else if (strcmp(flag.c_str(),"adjoint") == 0) {
        interp_plan = this->m_AdjointPlanVec;
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }
    
    ierr = VecGetArray(this->m_WorkVec, &p_WorkVec); CHKERRQ(ierr);
    interp_plan->interpolate( this->m_VecFieldGhost, 
                              3, 
                              nx,
                              isize,
                              istart,
                              this->isize_g, 
                              this->nlghost,
                              nl, 
                              this->nghost, 
                              p_WorkVec,
                              c_dims,
                              this->m_Opt->m_FFT.mpicomm, 
                              timers, 
                              this->m_tmpInterpol1, 
                              this->m_tmpInterpol2, 
                              this->m_texture, 
                              this->m_Opt->m_PDESolver.iporder, 
                              &(this->m_Opt->m_GPUtime));
    ierr = VecRestoreArray(this->m_WorkVec, &p_WorkVec);
    
    //std::cout << "================================ interpolated velocity field ====================================" << std::endl;
    //ierr = VecGetArray(this->m_WorkVec, &p_WorkVec);
    //// Small check for correctness of ghost vecfield
    //for (IntType i1 = 0; i1 < isize[0]; i1++) {
    //    for (IntType i2 = 0; i2 < isize[1]; i2++) {
    //        for (IntType i3 = 0; i3 < isize[2]; i3++) {
    //            i = reg::GetLinearIndex(i1, i2, i3, isize);
    //            printf("[%d], %f, %f, %f\n", i, p_WorkVec[0*nl+i], p_WorkVec[1*nl+i], p_WorkVec[2*nl+i]);
    //        }
    //    }
    //}
    //ierr = VecRestoreArray(this->m_WorkVec, &p_WorkVec);
    
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

    ZeitGeist_define(SL_INTERPOL);
    ZeitGeist_tick(SL_INTERPOL);
    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);


    if (strcmp(flag.c_str(),"state") == 0) {

        ierr = this->m_Xstate->GetArraysRead(xq1, xq2, xq3);
        gpuInterpVec3D(vx1, vx2, vx3, xq1, xq2, xq3, wx1, wx2, wx3, this->m_tmpInterpol1, this->m_tmpInterpol2, nx, static_cast<long int>(nl), this->m_texture, this->m_Opt->m_PDESolver.iporder, &(this->m_Opt->m_GPUtime));
        ierr = this->m_Xstate->RestoreArrays(); CHKERRQ(ierr);

    } else if (strcmp(flag.c_str(),"adjoint") == 0) {
        
        ierr = this->m_Xadjoint->GetArraysRead(xq1, xq2, xq3);
        gpuInterpVec3D(vx1, vx2, vx3, xq1, xq2, xq3, wx1, wx2, wx3, this->m_tmpInterpol1, this->m_tmpInterpol2, nx, static_cast<long int>(nl), this->m_texture, this->m_Opt->m_PDESolver.iporder, &(this->m_Opt->m_GPUtime));
        ierr = this->m_Xadjoint->RestoreArrays(); CHKERRQ(ierr);
       
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }
    
    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);
    ZeitGeist_tock(SL_INTERPOL);
    ZeitGeist_inc(SL_INTERPOL);
    ZeitGeist_inc(SL_INTERPOL);

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
    if (this->m_X == nullptr) {
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
    const ScalarType *p_x1=nullptr,*p_x2=nullptr,*p_x3=nullptr;
    int nx[3], isize_g[3], istart_g[3], c_dims[2], isize[3], istart[3], order;
    IntType nl;
    double timers[4] = {0,0,0,0};
    ScalarType* p_X;

    PetscFunctionBegin;
    
    if (this->m_Opt->m_Verbosity > 1) {
      std::string str = "Mapping query points\n";
      ierr = DbgMsgCall(str); CHKERRQ(ierr);
    }

    
    order  = this->m_Opt->m_PDESolver.iporder;

    for (int i = 0; i < 3; ++i){
        nx[i] = this->m_Opt->m_Domain.nx[i];
        isize[i] = this->m_Opt->m_Domain.isize[i];
        istart[i] = this->m_Opt->m_Domain.istart[i];
    }
    
    c_dims[0] = this->m_Opt->m_CartGridDims[0];
    c_dims[1] = this->m_Opt->m_CartGridDims[1];

    nl = this->m_Opt->m_Domain.nl;

    ierr = Assert(this->m_X != nullptr, "null pointer"); CHKERRQ(ierr);
     
    ierr = VecGetArray(this->m_X, &p_X); CHKERRQ(ierr);
    
    if (strcmp(flag.c_str(),"state")==0) {

        //this->m_StatePlan->scatter(this->m_Dofs[0], nx, isize, istart, nl, this->nghost, p_X, c_dims, this->m_Opt->m_FFT.mpicomm, timers);

        this->m_StatePlanVec->scatter(this->m_Dofs[1], nx, isize, istart, nl, this->nghost, p_X, c_dims, this->m_Opt->m_FFT.mpicomm, timers);

    } else if (strcmp(flag.c_str(),"adjoint")==0) {

        //this->m_AdjointPlanscatter(this->m_Dofs[0], nx, isize, istart, nl, this->nghost, p_X, c_dims, this->m_Opt->m_FFT.mpicomm, timers);

        this->m_AdjointPlanVec->scatter(this->m_Dofs[1], nx, isize, istart, nl, this->nghost, p_X, c_dims, this->m_Opt->m_FFT.mpicomm, timers);

    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(this->m_X, &p_X); CHKERRQ(ierr);

    this->m_Opt->IncreaseInterpTimers(timers);

    PetscFunctionReturn(ierr);
}

}  // namespace reg




#endif  // _SEMILAGRANGIAN_CPP_
