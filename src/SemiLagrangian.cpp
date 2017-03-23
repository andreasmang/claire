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

    this->m_ReadWrite = NULL;
    this->m_WorkVecField = NULL;

    this->m_X = NULL;
    this->m_XA = NULL;
    this->m_XS = NULL;

    this->m_StatePlan = NULL;
    this->m_AdjointPlan = NULL;
    this->m_StatePlanVec = NULL;
    this->m_AdjointPlanVec = NULL;

    this->m_VecFieldPlan = NULL;

    this->m_ScaFieldGhost = NULL;
    this->m_VecFieldGhost = NULL;

    this->m_Opt = NULL;

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
    if (this->m_XS != NULL) {
        delete [] this->m_XS;
        this->m_XS = NULL;
    }
    if (this->m_XA != NULL) {
        delete [] this->m_XA;
        this->m_XA = NULL;
    }

    if (this->m_AdjointPlan != NULL) {
        delete this->m_AdjointPlan;
        this->m_AdjointPlan = NULL;
    }
    if (this->m_StatePlan != NULL) {
        delete this->m_StatePlan;
        this->m_StatePlan = NULL;
    }
    if (this->m_StatePlanVec != NULL) {
        delete this->m_StatePlanVec;
        this->m_StatePlanVec = NULL;
    }
    if (this->m_AdjointPlanVec != NULL) {
        delete this->m_AdjointPlanVec;
        this->m_AdjointPlanVec = NULL;
    }
    if (this->m_VecFieldPlan!=NULL) {
        delete this->m_VecFieldPlan;
        this->m_VecFieldPlan = NULL;
    }

    if (this->m_ScaFieldGhost != NULL) {
        accfft_free(this->m_ScaFieldGhost);
        this->m_ScaFieldGhost = NULL;
    }

    if (this->m_VecFieldGhost != NULL) {
        accfft_free(this->m_VecFieldGhost);
        this->m_VecFieldGhost = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
PetscErrorCode SemiLagrangian::SetReadWrite(ReadWriteReg* readwrite) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(readwrite != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite = readwrite;

    PetscFunctionReturn(ierr);

}




/********************************************************************
 * @brief set work vector field to not have to allocate it locally
 *******************************************************************/
PetscErrorCode SemiLagrangian::SetWorkVecField(VecField* v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_WorkVecField = v;

    PetscFunctionReturn(ierr);

}





/********************************************************************
 * @brief compute the trajectory from the velocity field based
 * on an rk2 scheme (todo: make the velocity field a const vector)
 *******************************************************************/
PetscErrorCode SemiLagrangian::ComputeTrajectory(VecField* v, std::string flag) {
    PetscErrorCode ierr = 0;
    ScalarType ht, hthalf, hx[3], x1, x2, x3, scale;
    const ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL;
    ScalarType *p_vX1 = NULL, *p_vX2 = NULL, *p_vX3 = NULL;
    IntType isize[3], istart[3], l, i1, i2, i3, nl;
    std::stringstream ss;
    ScalarType* X = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_WorkVecField != NULL, "null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->GetDomainPara().nl;

    // switch between state and adjoint variable
    if (strcmp(flag.c_str(),"state") == 0) {
        if (this->m_XS == NULL) {
            if (this->m_Opt->GetVerbosity() > 2) {
                ss << " >> " << __func__ << ": allocation (size = " << 3*nl << ")";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.clear(); ss.str(std::string());
            }
            try {this->m_XS = new ScalarType [3*nl];}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        X = this->m_XS;
        scale = 1.0;
    } else if (strcmp(flag.c_str(),"adjoint") == 0) {
        if (this->m_XA == NULL) {
            if (this->m_Opt->GetVerbosity() > 2) {
                ss << " >> " << __func__ << ": allocation (size = " << 3*nl << ")";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.clear(); ss.str(std::string());
            }
            try {this->m_XA = new ScalarType [3*nl];}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        X = this->m_XA;
        scale = -1.0;
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    for (int i = 0; i < 3; ++i) {
        hx[i] = this->m_Opt->GetDomainPara().hx[i];
        isize[i] = this->m_Opt->GetDomainPara().isize[i];
        istart[i] = this->m_Opt->GetDomainPara().istart[i];
    }

    // \tilde{X} = x - ht v
    ierr = v->GetArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
    for (i1 = 0; i1 < isize[0]; ++i1) {   // x1
        for (i2 = 0; i2 < isize[1]; ++i2) {   // x2
            for (i3 = 0; i3 < isize[2]; ++i3) {   // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                l = GetLinearIndex(i1, i2, i3, isize);

                X[l*3+0] = x1 - scale*ht*p_v1[l];
                X[l*3+1] = x2 - scale*ht*p_v2[l];
                X[l*3+2] = x3 - scale*ht*p_v3[l];
            }  // i1
        }  // i2
    }  // i3
}  // pragma omp for
    ierr = v->RestoreArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    // normalize to [0,1]
    for (IntType i = 0; i < 3*nl; ++i) {
        X[i] /= (2*PETSC_PI);
    }

    // interpolate velocity field v(X)
    ierr = this->MapCoordinateVector(flag); CHKERRQ(ierr);
    ierr = this->Interpolate(this->m_WorkVecField, v, flag); CHKERRQ(ierr);

    // normalize to [0,2*pi]
    for (IntType i = 0; i < 3*nl; ++i) {
        X[i] *= (2*PETSC_PI);
    }

    // X = x - 0.5*ht*(v + v(x - ht v))
    ierr = v->GetArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField->GetArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
    for (i1 = 0; i1 < isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + istart[2]);

                // compute linear / flat index
                l = GetLinearIndex(i1, i2, i3, isize);

                X[l*3+0] = x1 - scale*hthalf*(p_vX1[l] + p_v1[l]);
                X[l*3+1] = x2 - scale*hthalf*(p_vX2[l] + p_v2[l]);
                X[l*3+2] = x3 - scale*hthalf*(p_vX3[l] + p_v3[l]);
            }  // i1
        }  // i2
    }  // i3
}  // pragma omp for
    ierr = this->m_WorkVecField->RestoreArrays(p_vX1, p_vX2, p_vX3); CHKERRQ(ierr);
    ierr = v->RestoreArraysRead(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    // normalize to [0,1]
    for (IntType i = 0; i < 3*nl; ++i) {
        X[i] /= (2*PETSC_PI);
    }

    ierr =this->MapCoordinateVector(flag); CHKERRQ(ierr);

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
    IntType nl, g_alloc_max;
    std::stringstream ss;
    accfft_plan_t<ScalarType, ComplexType, FFTWPlanType>* plan = NULL;  ///< accfft plan
//    accfft_plan* plan = NULL;
    double timers[4] = {0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(xi != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(xo != NULL, "null pointer"); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->GetDomainPara().nx[i]);
        isize[i] = static_cast<int>(this->m_Opt->GetDomainPara().isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->GetDomainPara().istart[i]);
    }

    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    nl = this->m_Opt->GetDomainPara().nl;
    order = this->m_Opt->GetPDESolverPara().interpolationorder;
    nghost = order;

    neval = static_cast<int>(nl);
    ierr = Assert(neval != 0, "size problem"); CHKERRQ(ierr);

    // deal with ghost points
    plan = this->m_Opt->GetFFT().plan;
    g_alloc_max = accfft_ghost_xyz_local_size_dft_r2c(plan, nghost, isize_g, istart_g);
    if (this->m_ScaFieldGhost == NULL) {
        ierr = Assert(g_alloc_max > 0 && g_alloc_max < std::numeric_limits<int>::max(), "allocation error"); CHKERRQ(ierr);
        if (this->m_Opt->GetVerbosity() > 2) {
            ss << " >> " << __func__ << ": allocation (size = " << g_alloc_max << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
        this->m_ScaFieldGhost = reinterpret_cast<ScalarType*>(accfft_alloc(g_alloc_max));
    }

    // compute interpolation for all components of the input scalar field
    accfft_get_ghost_xyz(plan, nghost, isize_g, xi, this->m_ScaFieldGhost);

    if (strcmp(flag.c_str(), "state") == 0) {
        ierr = Assert(this->m_XS != NULL, "state X is null pointer"); CHKERRQ(ierr);
        this->m_StatePlan->interpolate(this->m_ScaFieldGhost, 1, nx, isize, istart,
                                       neval, nghost, xo, c_dims,
                                       this->m_Opt->GetFFT().mpicomm, timers);
 //       this->m_StatePlan->high_order_interpolate(this->m_ScaFieldGhost, 1, nx, isize, istart,
 //                                      neval, nghost, xo, c_dims,
 //                                      this->m_Opt->GetFFT().mpicomm, timers, order);
    } else if (strcmp(flag.c_str(), "adjoint") == 0) {
        ierr = Assert(this->m_XA != NULL, "adjoint X is null pointer"); CHKERRQ(ierr);
        this->m_AdjointPlan->interpolate(this->m_ScaFieldGhost, 1, nx, isize, istart,
                                         neval, nghost, xo, c_dims,
                                         this->m_Opt->GetFFT().mpicomm, timers);
 //       this->m_AdjointPlan->high_order_interpolate(this->m_ScaFieldGhost, 1, nx, isize, istart,
 //                                      neval, nghost, xo, c_dims,
 //                                      this->m_Opt->GetFFT().mpicomm, timers, order);
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

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
    accfft_plan_t<ScalarType, ComplexType, FFTWPlanType>* plan = NULL;  ///< accfft plan
//    accfft_plan* plan;
    std::stringstream ss;
    IntType nl, nlghost, g_alloc_max;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(vx1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vx2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(vx3 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(wx3 != NULL, "null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->GetDomainPara().nl;
    order = this->m_Opt->GetPDESolverPara().interpolationorder;
    nghost = order;

    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<int>(this->m_Opt->GetDomainPara().isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->GetDomainPara().istart[i]);
    }

    // get network dimensions
    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    if (this->m_X == NULL) {
        ierr = Assert(3*nl > 0 && 3*nl < std::numeric_limits<IntType>::max(), "allocation error"); CHKERRQ(ierr);
        if (this->m_Opt->GetVerbosity() > 2) {
            ss << " >> " << __func__ << ": allocation (size = " << 3*nl << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
        try {this->m_X = new ScalarType [3*nl];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    for (IntType i = 0; i < nl; ++i) {
        this->m_X[0*nl+i] = vx1[i];
        this->m_X[1*nl+i] = vx2[i];
        this->m_X[2*nl+i] = vx3[i];
    }

    // get ghost sizes
    plan = this->m_Opt->GetFFT().plan;
    g_alloc_max = accfft_ghost_xyz_local_size_dft_r2c(plan, nghost, isize_g, istart_g);
    ierr = Assert(g_alloc_max != 0, "alloc problem"); CHKERRQ(ierr);

    // get nl for ghosts
    nlghost = 1;
    for (int i = 0; i < 3; ++i) {
        nlghost *= static_cast<IntType>(isize_g[i]);
    }

    // deal with ghost points
    if (this->m_VecFieldGhost == NULL) {
        ierr = Assert(g_alloc_max > 0 && g_alloc_max < std::numeric_limits<int>::max(), "allocation error"); CHKERRQ(ierr);
        if (this->m_Opt->GetVerbosity() > 2) {
            ss << " >> " << __func__ << ": allocation (size = " << g_alloc_max << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
        this->m_VecFieldGhost = reinterpret_cast<ScalarType*>(accfft_alloc(3*g_alloc_max));
    }


    // do the communication for the ghost points
    for (int i = 0; i < 3; i++) {
        accfft_get_ghost_xyz(plan, nghost, isize_g, &this->m_X[i*nl],
                             &this->m_VecFieldGhost[i*nlghost]);
    }

    if (strcmp(flag.c_str(),"state") == 0) {
        ierr = Assert(this->m_XS != NULL, "state X null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_StatePlanVec != NULL, "state X null pointer"); CHKERRQ(ierr);
        this->m_StatePlanVec->interpolate(this->m_VecFieldGhost, 3, nx, isize, istart,
                                          nl, nghost, this->m_X, c_dims,
                                          this->m_Opt->GetFFT().mpicomm, timers);
 //       this->m_StatePlanVec->high_order_interpolate(this->m_VecFieldGhost, 3, nx, isize, istart,
 //                                         nl, nghost, this->m_X, c_dims,
 //                                         this->m_Opt->GetFFT().mpicomm, timers, order);
    } else if (strcmp(flag.c_str(),"adjoint") == 0) {
        ierr = Assert(this->m_XA != NULL, "adjoint X null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_AdjointPlanVec != NULL, "state X null pointer"); CHKERRQ(ierr);

        this->m_AdjointPlanVec->interpolate(this->m_VecFieldGhost, 3, nx, isize, istart,
                                            nl, nghost, this->m_X, c_dims,
                                            this->m_Opt->GetFFT().mpicomm, timers);
 //       this->m_AdjointPlanVec->high_order_interpolate(this->m_VecFieldGhost, 3, nx, isize, istart,
 //                                         nl, nghost, this->m_X, c_dims,
 //                                         this->m_Opt->GetFFT().mpicomm, timers, order);
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

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
 * @brief interpolate vector field
 *******************************************************************/
PetscErrorCode SemiLagrangian::Interpolate( ScalarType* wx1, ScalarType* wx2, ScalarType* wx3,
                                            ScalarType* vx1, ScalarType* vx2, ScalarType* vx3,
                                            ScalarType* yx1, ScalarType* yx2, ScalarType* yx3) {
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

    ierr = Assert(yx1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(yx2 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(yx3 != NULL, "null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->GetDomainPara().nl;
    order = this->m_Opt->GetPDESolverPara().interpolationorder;
    nghost = order;

    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->GetNumGridPoints(i));
        isize[i] = static_cast<int>(this->m_Opt->GetDomainPara().isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->GetDomainPara().istart[i]);
    }

    // get network dimensions
    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    if (this->m_X == NULL) {
        ierr = Assert(3*nl > 0 && 3*nl < std::numeric_limits<IntType>::max(), "allocation error"); CHKERRQ(ierr);
        if (this->m_Opt->GetVerbosity() > 2) {
            ss << " >> " << __func__ << ": allocation (size = " << 3*nl << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
        try {this->m_X = new ScalarType [3*nl];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // create planer
    if (this->m_VecFieldPlan == NULL) {
        if (this->m_Opt->GetVerbosity() > 2) {
            ss << " >> " << __func__ << ": allocation of plan";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
        try {this->m_VecFieldPlan = new Interp3_Plan;}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        this->m_VecFieldPlan->allocate(nl, 3);
    }

    for (IntType i = 0; i < nl; ++i) {
        this->m_X[i*3+0] = yx1[i]/(2*PETSC_PI);
        this->m_X[i*3+1] = yx2[i]/(2*PETSC_PI);
        this->m_X[i*3+2] = yx3[i]/(2*PETSC_PI);
    }

    // scatter
    this->m_VecFieldPlan->scatter(3, nx, isize, istart, nl, nghost,
                                  this->m_X, c_dims, this->m_Opt->GetFFT().mpicomm, timers);

    for (IntType i = 0; i < nl; ++i) {
        this->m_X[0*nl+i] = vx1[i];
        this->m_X[1*nl+i] = vx2[i];
        this->m_X[2*nl+i] = vx3[i];
    }

    // get ghost sizes
    nalloc = accfft_ghost_xyz_local_size_dft_r2c(this->m_Opt->GetFFT().plan,
                                                 nghost, isize_g, istart_g);

    // get nl for ghosts
    nlghost = 1;
    for (int i = 0; i < 3; ++i) {
        nlghost *= static_cast<IntType>(isize_g[i]);
    }

    // deal with ghost points
    if (this->m_VecFieldGhost == NULL) {
        ierr = Assert(3*nalloc > 0 && 3*nalloc < std::numeric_limits<int>::max(), "allocation error"); CHKERRQ(ierr);
        if (this->m_Opt->GetVerbosity() > 2) {
            ss << " >> " << __func__ << ": allocation (size = " << 3*nalloc << ")";
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
        this->m_VecFieldGhost = reinterpret_cast<ScalarType*>(accfft_alloc(3*nalloc));
    }


    // do the communication for the ghost points
    for (int i = 0; i < 3; i++) {
        accfft_get_ghost_xyz(this->m_Opt->GetFFT().plan, nghost, isize_g,
                             &this->m_X[i*nl], &this->m_VecFieldGhost[i*nlghost]);
    }

    this->m_VecFieldPlan->interpolate(this->m_VecFieldGhost, 3, nx, isize, istart,
                                      nl, nghost, this->m_X, c_dims,
                                      this->m_Opt->GetFFT().mpicomm, timers);
 //   this->m_VecFieldPlan->high_order_interpolate(this->m_VecFieldGhost, 3, nx, isize, istart,
 //                                     nl, nghost, this->m_X, c_dims,
 //                                     this->m_Opt->GetFFT().mpicomm, timers, order);

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
 * @brief change from lexicographical ordering to xyz
 * @param flag flag to switch between forward and adjoint solves
 *******************************************************************/
PetscErrorCode SemiLagrangian::MapCoordinateVector(std::string flag) {
    PetscErrorCode ierr;
    int nx[3], nl, isize[3], istart[3], nghost;
    int c_dims[2];
    double timers[4] = {0, 0, 0, 0};
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // get sizes
    nl = static_cast<int>(this->m_Opt->GetDomainPara().nl);
    nghost = this->m_Opt->GetPDESolverPara().interpolationorder;

    ierr = Assert(this->m_Opt->GetFFT().mpicomm != NULL, "null pointer"); CHKERRQ(ierr);

    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Opt->GetDomainPara().nx[i]);
        isize[i] = static_cast<int>(this->m_Opt->GetDomainPara().isize[i]);
        istart[i] = static_cast<int>(this->m_Opt->GetDomainPara().istart[i]);
    }

    c_dims[0] = this->m_Opt->GetNetworkDims(0);
    c_dims[1] = this->m_Opt->GetNetworkDims(1);

    if (strcmp(flag.c_str(),"state") == 0) {
        // characteristic for state equation should have been computed already
        ierr = Assert(this->m_XS != NULL, "null pointer"); CHKERRQ(ierr);

        // create planer
        if (this->m_StatePlan == NULL) {
            if (this->m_Opt->GetVerbosity() > 2) {
                ss << " >> " << __func__ << ": allocation of plan (state, nl=" << nl << ")";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.clear(); ss.str(std::string());
            }
            try {this->m_StatePlan = new Interp3_Plan();}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
            this->m_StatePlan->allocate(nl, 1);
        }

        // scatter
        this->m_StatePlan->scatter(1, nx, isize, istart, nl, nghost, this->m_XS,
                                   c_dims, this->m_Opt->GetFFT().mpicomm, timers);

        // create planer
        if (this->m_StatePlanVec == NULL) {
            if (this->m_Opt->GetVerbosity() > 2) {
                ss << " >> " << __func__ << ": allocation of plan (state; vector field, nl=" << nl << ")";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.clear(); ss.str(std::string());
            }
            try {this->m_StatePlanVec = new Interp3_Plan();}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
            this->m_StatePlanVec->allocate(nl, 3);
        }

        // scatter
        this->m_StatePlanVec->scatter(3, nx, isize, istart, nl, nghost, this->m_XS,
                                      c_dims, this->m_Opt->GetFFT().mpicomm, timers);
    } else if (strcmp(flag.c_str(),"adjoint") == 0) {
        // characteristic for adjoint equation should
        // have been computed already
        ierr = Assert(this->m_XA != NULL, "null pointer"); CHKERRQ(ierr);

        // create planer
        if (this->m_AdjointPlan == NULL) {
            if (this->m_Opt->GetVerbosity() > 2) {
                ss << " >> " << __func__ << ": allocation of plan (adjoint, nl=" << nl << ")";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.clear(); ss.str(std::string());
            }
            try {this->m_AdjointPlan = new Interp3_Plan();}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
            this->m_AdjointPlan->allocate(nl, 1);
        }

        // scatter
        this->m_AdjointPlan->scatter(1, nx, isize, istart, nl, nghost, this->m_XA,
                                     c_dims, this->m_Opt->GetFFT().mpicomm, timers);

        // create planer
        if (this->m_AdjointPlanVec == NULL) {
            if (this->m_Opt->GetVerbosity() > 2) {
                ss << " >> " << __func__ << ": allocation of plan (adjoint; vector field, nl=" << nl << ")";
                ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.clear(); ss.str(std::string());
            }
            try {this->m_AdjointPlanVec = new Interp3_Plan();}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
            this->m_AdjointPlanVec->allocate(nl, 3);
        }

        // scatter
        this->m_AdjointPlanVec->scatter(3, nx, isize, istart, nl, nghost, this->m_XA,
                                        c_dims, this->m_Opt->GetFFT().mpicomm, timers);
    } else {
        ierr = ThrowError("flag wrong"); CHKERRQ(ierr);
    }

    this->m_Opt->IncreaseInterpTimers(timers);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




}  // namespace reg




#endif  // _SEMILAGRANGIAN_CPP_
