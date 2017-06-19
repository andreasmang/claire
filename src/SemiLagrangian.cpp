/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _SEMILAGRANGIAN_CPP_
#define _SEMILAGRANGIAN_CPP_

#include "SemiLagrangian.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
SemiLagrangian::SemiLagrangian() {
    this->Initialize();
}




/********************************************************************
 * @brief default constructor
 *******************************************************************/
SemiLagrangian::SemiLagrangian(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
SemiLagrangian::~SemiLagrangian() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init class variables
 *******************************************************************/
PetscErrorCode SemiLagrangian::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_X = NULL;
    this->m_rkorder = 2;
    this->m_WorkVecField1 = NULL;
    this->m_WorkVecField2 = NULL;

    this->m_StatePlan = NULL;
    this->m_AdjointPlan = NULL;

    this->m_ScaFieldGhost = NULL;
    this->m_VecFieldGhost = NULL;

    this->m_Opt = NULL;
    this->m_Dofs[0] = 1;
    this->m_Dofs[1] = 3;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clears memory
 *******************************************************************/
PetscErrorCode SemiLagrangian::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_X != NULL) {
        delete [] this->m_X;
        this->m_X = NULL;
    }

    if (this->m_AdjointPlan != NULL) {
        delete this->m_AdjointPlan;
        this->m_AdjointPlan = NULL;
    }
    if (this->m_StatePlan != NULL) {
        delete this->m_StatePlan;
        this->m_StatePlan = NULL;
    }

    if (this->m_ScaFieldGhost != NULL) {
        accfft_free(this->m_ScaFieldGhost);
        this->m_ScaFieldGhost = NULL;
    }

    if (this->m_VecFieldGhost != NULL) {
        accfft_free(this->m_VecFieldGhost);
        this->m_VecFieldGhost = NULL;
    }

    if (this->m_WorkVecField2 != NULL) {
        delete this->m_WorkVecField2;
        this->m_WorkVecField2 = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set work vector field to not have to allocate it locally
 *******************************************************************/
PetscErrorCode SemiLagrangian::SetWorkVecField(VecField* x) {
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
PetscErrorCode SemiLagrangian::ComputeTrajectory(VecField* v, std::string flag) {
    PetscErrorCode ierr = 0;
    IntType nl;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->m_Domain.nl;

    // if trajectory has not yet been allocated, allocate
    if (this->m_X == NULL) {
        try {this->m_X = new ScalarType[3*nl];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // compute trajectory
    if (this->m_rkorder == 2) {
        ierr = this->ComputeTrajectoryRK2(v, flag); CHKERRQ(ierr);
    } else if (this->m_rkorder == 4) {
        ierr = this->ComputeTrajectoryRK4(v, flag); CHKERRQ(ierr);
    } else {
        ierr = ThrowError("rk order not implemented"); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute the trajectory from the velocity field based
 * on an rk2 scheme (todo: make the velocity field a const vector)
 *******************************************************************/
PetscErrorCode SemiLagrangian::ComputeTrajectoryRK2(VecField* v, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType ht, hthalf, hx[3], x1, x2, x3, scale;
    const ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL;
    ScalarType *p_vX1 = NULL, *p_vX2 = NULL, *p_vX3 = NULL;
    IntType isize[3], istart[3], l, i1, i2, i3;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_WorkVecField1 != NULL, "null pointer"); CHKERRQ(ierr);

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
        hx[i]     = this->m_Opt->m_Domain.hx[i];
        isize[i]  = this->m_Opt->m_Domain.isize[i];
        istart[i] = this->m_Opt->m_Domain.istart[i];
    }


    // \tilde{X} = x - ht v
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
                this->m_X[l*3+0] = (x1 - scale*ht*p_v1[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+1] = (x2 - scale*ht*p_v2[l])/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+2] = (x3 - scale*ht*p_v3[l])/(2.0*PETSC_PI); // normalized to [0,1]
            }  // i1
        }  // i2
    }  // i3
    ierr = v->RestoreArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);


    // communicate the characteristic
    ierr = this->CommunicateCoord(flag); CHKERRQ(ierr);

    // interpolate velocity field v(X)
    ierr = this->Interpolate(this->m_WorkVecField1, v, flag); CHKERRQ(ierr);

    // X = x - 0.5*ht*(v + v(x - ht v))
    ierr = v->GetArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);
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
                this->m_X[l*3+0] = (x1 - scale*hthalf*(p_vX1[l] + p_v1[l]))/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+1] = (x2 - scale*hthalf*(p_vX2[l] + p_v2[l]))/(2.0*PETSC_PI); // normalized to [0,1]
                this->m_X[l*3+2] = (x3 - scale*hthalf*(p_vX3[l] + p_v3[l]))/(2.0*PETSC_PI); // normalized to [0,1]
            }  // i1
        }  // i2
    }  // i3
    ierr = this->m_WorkVecField1->RestoreArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
    ierr = v->RestoreArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    // communicate the characteristic
    ierr = this->CommunicateCoord(flag); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute the trajectory from the velocity field based
 * on an rk2 scheme (todo: make the velocity field a const vector)
 *******************************************************************/
PetscErrorCode SemiLagrangian::ComputeTrajectoryRK4(VecField* v, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType ht, hthalf, hx[3], x1, x2, x3, scale;
    const ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL;
    ScalarType *p_vX1 = NULL, *p_vX2 = NULL, *p_vX3 = NULL,
               *p_f1 = NULL, *p_f2 = NULL, *p_f3 = NULL;
    IntType isize[3], istart[3], l, i1, i2, i3;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

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
        hx[i]     = this->m_Opt->m_Domain.hx[i];
        isize[i]  = this->m_Opt->m_Domain.isize[i];
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

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief interpolate scalar field
 *******************************************************************/
PetscErrorCode SemiLagrangian::Interpolate(Vec* xo, Vec xi, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType *p_xo = NULL, *p_xi = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(*xo != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert( xi != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = VecGetArray( xi, &p_xi); CHKERRQ(ierr);
    ierr = VecGetArray(*xo, &p_xo); CHKERRQ(ierr);

    ierr = this->Interpolate(p_xo, p_xi, flag); CHKERRQ(ierr);

    ierr = VecRestoreArray(*xo, &p_xo); CHKERRQ(ierr);
    ierr = VecRestoreArray( xi, &p_xi); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief interpolate scalar field
 *******************************************************************/
PetscErrorCode SemiLagrangian::Interpolate(ScalarType* xo, ScalarType* xi, std::string flag) {
    PetscErrorCode ierr = 0;
    int nx[3], isize_g[3], isize[3], istart_g[3], istart[3], c_dims[2], neval, order, nghost;
    IntType nl, nalloc;
    std::stringstream ss;
    double timers[4] = {0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(xi != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xo != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    nl     = this->m_Opt->m_Domain.nl;
    order  = this->m_Opt->m_PDESolver.interpolationorder;
    nghost = order;
    neval  = static_cast<int>(nl);

    for (int i = 0; i < 3; ++i) {
        nx[i]     = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i]  = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }

    c_dims[0] = this->m_Opt->m_CartGridDims[0];
    c_dims[1] = this->m_Opt->m_CartGridDims[1];

    // deal with ghost points
    nalloc = accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->m_FFT.plan, nghost, isize_g, istart_g);

    // if scalar field with ghost points has not been allocated
    if (this->m_ScaFieldGhost == NULL) {
        this->m_ScaFieldGhost = reinterpret_cast<ScalarType*>(accfft_alloc(nalloc));
    }

    // assign ghost points based on input scalar field
    accfft_get_ghost_xyz(this->m_Opt->m_FFT.plan, nghost, isize_g, xi, this->m_ScaFieldGhost);

    // compute interpolation for all components of the input scalar field
    if (strcmp(flag.c_str(), "state") == 0) {
        this->m_StatePlan->interpolate(this->m_ScaFieldGhost, nx, isize, istart,
                                       neval, nghost, xo, c_dims, this->m_Opt->m_FFT.mpicomm, timers, 0);
    } else if (strcmp(flag.c_str(), "adjoint") == 0) {
        this->m_AdjointPlan->interpolate(this->m_ScaFieldGhost, nx, isize, istart,
                                       neval, nghost, xo, c_dims, this->m_Opt->m_FFT.mpicomm, timers, 0);
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }
    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);
    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IP);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief interpolate vector field
 *******************************************************************/
PetscErrorCode SemiLagrangian::Interpolate(VecField* vo, VecField* vi, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType *p_vix1 = NULL, *p_vix2 = NULL, *p_vix3 = NULL,
               *p_vox1 = NULL, *p_vox2 = NULL, *p_vox3 = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vi != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vo != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = vi->GetArrays(p_vix1, p_vix2, p_vix3); CHKERRQ(ierr);
    ierr = vo->GetArrays(p_vox1, p_vox2, p_vox3); CHKERRQ(ierr);

    ierr = this->Interpolate(p_vox1, p_vox2, p_vox3, p_vix1, p_vix2, p_vix3, flag); CHKERRQ(ierr);

    ierr = vo->RestoreArrays(p_vox1, p_vox2, p_vox3); CHKERRQ(ierr);
    ierr = vi->RestoreArrays(p_vix1, p_vix2, p_vix3); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief interpolate vector field
 *******************************************************************/
PetscErrorCode SemiLagrangian::Interpolate(ScalarType* wx1, ScalarType* wx2, ScalarType* wx3,
                                           ScalarType* vx1, ScalarType* vx2, ScalarType* vx3, std::string flag) {
    PetscErrorCode ierr = 0;
    int nx[3], isize_g[3], isize[3], istart_g[3], istart[3], c_dims[2], nghost, order;
    double timers[4] = {0, 0, 0, 0};
    std::stringstream ss;
    IntType nl, nlghost, nalloc;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vx1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vx2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vx3 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx3 != NULL, "null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->m_Domain.nl;
    order = this->m_Opt->m_PDESolver.interpolationorder;
    nghost = order;

    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->m_Domain.nx[i]);
        isize[i] = static_cast<int>(this->m_Opt->m_Domain.isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->m_Domain.istart[i]);
    }

    // get network dimensions
    c_dims[0] = this->m_Opt->m_CartGridDims[0];
    c_dims[1] = this->m_Opt->m_CartGridDims[1];

    if (this->m_X == NULL) {
        try {this->m_X = new ScalarType [3*nl];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // copy data to a flat vector
    for (IntType i = 0; i < nl; ++i) {
        this->m_X[0*nl+i] = vx1[i];
        this->m_X[1*nl+i] = vx2[i];
        this->m_X[2*nl+i] = vx3[i];
    }

    ierr = this->m_Opt->StartTimer(IPSELFEXEC); CHKERRQ(ierr);

    // get ghost sizes
    nalloc = accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->m_FFT.plan, nghost, isize_g, istart_g);

    // get nl for ghosts
    nlghost = 1;
    for (int i = 0; i < 3; ++i) {
        nlghost *= static_cast<IntType>(isize_g[i]);
    }

    // deal with ghost points
    if (this->m_VecFieldGhost == NULL) {
        this->m_VecFieldGhost = reinterpret_cast<ScalarType*>(accfft_alloc(3*nalloc));
    }


    // do the communication for the ghost points
    for (int i = 0; i < 3; i++) {
        accfft_get_ghost_xyz(this->m_Opt->m_FFT.plan, nghost, isize_g, &this->m_X[i*nl],
                             &this->m_VecFieldGhost[i*nlghost]);
    }

    if (strcmp(flag.c_str(),"state") == 0) {
        ierr = Assert(this->m_StatePlan != NULL, "null pointer"); CHKERRQ(ierr);
        this->m_StatePlan->interpolate(this->m_VecFieldGhost, nx, isize, istart,
                                       nl, nghost, this->m_X, c_dims, this->m_Opt->m_FFT.mpicomm, timers, 1);
    } else if (strcmp(flag.c_str(),"adjoint") == 0) {
        ierr = Assert(this->m_AdjointPlan != NULL, "null pointer"); CHKERRQ(ierr);
        this->m_AdjointPlan->interpolate(this->m_VecFieldGhost, nx, isize, istart,
                                         nl, nghost, this->m_X, c_dims, this->m_Opt->m_FFT.mpicomm, timers, 1);
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

    ierr = this->m_Opt->StopTimer(IPSELFEXEC); CHKERRQ(ierr);

    for (IntType i = 0; i < nl; ++i) {
        wx1[i] = this->m_X[0*nl+i];
        wx2[i] = this->m_X[1*nl+i];
        wx3[i] = this->m_X[2*nl+i];
    }

    this->m_Opt->IncreaseInterpTimers(timers);
    this->m_Opt->IncrementCounter(IPVEC);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief communicate the coordinate vector (query points)
 * @param flag to switch between forward and adjoint solves
 *******************************************************************/
PetscErrorCode SemiLagrangian::CommunicateCoord(std::string flag) {
    PetscErrorCode ierr;
    int nx[3], nl, isize[3], istart[3], nghost;
    int c_dims[2];
    double timers[4] = {0, 0, 0, 0};
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_Opt->m_FFT.mpicomm != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nl     = static_cast<int>(this->m_Opt->m_Domain.nl);
    nghost = this->m_Opt->m_PDESolver.interpolationorder;
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

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




}  // namespace reg




#endif  // _SEMILAGRANGIAN_CPP_
