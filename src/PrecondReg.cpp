/*************************************************************************
 *  Copyright (c) 2015-2016.
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

#ifndef _PRECONDREG_CPP_
#define _PRECONDREG_CPP_

#include "PrecondReg.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
PrecondReg::PrecondReg() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
PrecondReg::~PrecondReg() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
PrecondReg::PrecondReg(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode PrecondReg::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt = NULL;     ///< options (default; on grid we solve)
    this->m_CoarseGrid.m_Opt = NULL;   ///< options for coarse grid

    this->m_MatVec = NULL;              ///< pointer to matvec in krylov method
    this->m_MatVecEigEst = NULL;        ///< pointer to matvec in krylov method
    this->m_KrylovMethod = NULL;        ///< pointer to krylov method
    this->m_KrylovMethodEigEst = NULL;  ///< pointer to krylov method
    this->m_RandomNumGen = NULL;        ///< pointer to krylov method

    this->m_CoarseGrid.x = NULL;     ///< container for input to hessian matvec on coarse grid
    this->m_CoarseGrid.y = NULL;    ///< container for hessian matvec on coarse grid

    this->m_PreProc = NULL;         ///< pointer to preprocessing operator
    this->m_CoarseGrid.m_OptimizationProblem = NULL;   ///< options for coarse grid

    this->m_ControlVariable = NULL;     ///< control variable on fine grid
    this->m_IncControlVariable = NULL;  ///< incremental control variable on fine grid

    this->m_CoarseGrid.m_StateVariable = NULL;         ///< state variable on coarse grid
    this->m_CoarseGrid.m_AdjointVariable = NULL;       ///< adjoint variable on coarse grid
    this->m_CoarseGrid.m_ControlVariable = NULL;       ///< control variable on coarse grid
    this->m_CoarseGrid.m_IncControlVariable = NULL;    ///< incremental control variable on coarse grid

    this->m_WorkVecField = NULL;            ///< temporary vector field
    this->m_WorkScaField1 = NULL;           ///< temporary scalar field
    this->m_WorkScaField2 = NULL;           ///< temporary scalar field
    this->m_CoarseGrid.m_WorkScaField1 = NULL;     ///< temporary scalar field (coarse level)
    this->m_CoarseGrid.m_WorkScaField2 = NULL;     ///< temporary scalar field (coarse level)

    this->m_CoarseGrid.setupdone = false;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode PrecondReg::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_KrylovMethod != NULL) {
        ierr = KSPDestroy(&this->m_KrylovMethod); CHKERRQ(ierr);
        this->m_KrylovMethod = NULL;
    }
    if (this->m_MatVec != NULL) {
        ierr = MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }
    if (this->m_MatVecEigEst != NULL) {
        ierr = MatDestroy(&this->m_MatVecEigEst); CHKERRQ(ierr);
        this->m_MatVecEigEst = NULL;
    }


    if (this->m_CoarseGrid.x != NULL) {
        ierr = VecDestroy(&this->m_CoarseGrid.x); CHKERRQ(ierr);
        this->m_CoarseGrid.x = NULL;
    }
    if (this->m_CoarseGrid.y != NULL) {
        ierr = VecDestroy(&this->m_CoarseGrid.y); CHKERRQ(ierr);
        this->m_CoarseGrid.y = NULL;
    }
    if (this->m_CoarseGrid.m_StateVariable != NULL) {
        ierr = VecDestroy(&this->m_CoarseGrid.m_StateVariable); CHKERRQ(ierr);
        this->m_CoarseGrid.m_StateVariable = NULL;
    }
    if (this->m_CoarseGrid.m_AdjointVariable != NULL) {
        ierr = VecDestroy(&this->m_CoarseGrid.m_AdjointVariable); CHKERRQ(ierr);
        this->m_CoarseGrid.m_AdjointVariable = NULL;
    }


    if (this->m_WorkVecField != NULL) {
        delete this->m_WorkVecField;
        this->m_WorkVecField=NULL;
    }
    if (this->m_WorkScaField1 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField1); CHKERRQ(ierr);
        this->m_WorkScaField1 = NULL;
    }
    if (this->m_WorkScaField2 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField2); CHKERRQ(ierr);
        this->m_WorkScaField2 = NULL;
    }
    if (this->m_CoarseGrid.m_WorkScaField1 != NULL) {
        ierr = VecDestroy(&this->m_CoarseGrid.m_WorkScaField1); CHKERRQ(ierr);
        this->m_CoarseGrid.m_WorkScaField1 = NULL;
    }
    if (this->m_CoarseGrid.m_WorkScaField2 != NULL) {
        ierr = VecDestroy(&this->m_CoarseGrid.m_WorkScaField2); CHKERRQ(ierr);
        this->m_CoarseGrid.m_WorkScaField2 = NULL;
    }


    if (this->m_CoarseGrid.m_OptimizationProblem != NULL) {
        delete this->m_CoarseGrid.m_OptimizationProblem;
        this->m_CoarseGrid.m_OptimizationProblem = NULL;
    }

    if (this->m_ControlVariable != NULL) {
        delete this->m_ControlVariable;
        this->m_ControlVariable=NULL;
    }
    if (this->m_IncControlVariable != NULL) {
        delete this->m_IncControlVariable;
        this->m_IncControlVariable = NULL;
    }


    if (this->m_CoarseGrid.m_ControlVariable != NULL) {
        delete this->m_CoarseGrid.m_ControlVariable;
        this->m_CoarseGrid.m_ControlVariable = NULL;
    }
    if (this->m_CoarseGrid.m_IncControlVariable != NULL) {
        delete this->m_CoarseGrid.m_IncControlVariable;
        this->m_CoarseGrid.m_IncControlVariable = NULL;
    }

    if (this->m_CoarseGrid.m_Opt != NULL) {
        delete this->m_CoarseGrid.m_Opt;
        this->m_CoarseGrid.m_Opt = NULL;
    }

    if (this->m_RandomNumGen != NULL) {
        ierr = PetscRandomDestroy(&this->m_RandomNumGen); CHKERRQ(ierr);
        this->m_RandomNumGen=NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them)
 *******************************************************************/
PetscErrorCode PrecondReg::SetProblem(PrecondReg::OptProbType* optprob) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(optprob != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_OptimizationProblem = optprob;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them)
 *******************************************************************/
PetscErrorCode PrecondReg::SetPreProc(PreProcReg* preproc) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(preproc != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_PreProc = preproc;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief make sure that we re-initiate/recompute important
 * quantities; implemented to allow multiple calls of the solver
 * without destroying it;
 *******************************************************************/
PetscErrorCode PrecondReg::Reset() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // switch case for choice of preconditioner
    switch(this->m_Opt->GetKrylovSolverPara().pctype) {
        case NOPC:
        {
            // no need to do anything
            break;
        }
        case INVREG:
        {
            // no need to do anything
            break;
        }
        case TWOLEVEL:
        {
            // in case we call the solver multiple times (for
            // instance when we do parameter continuation) without
            // destroying the preconditioner, it is necessary to
            // recompute eigenvalues
            this->m_Opt->KrylovMethodEigValsEstimated(false);
            break;
        }
        default:
        {
            ierr = ThrowError("preconditioner not defined"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief setup phase of preconditioner
 *******************************************************************/
PetscErrorCode PrecondReg::DoSetup() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // start timer
    ierr = this->m_Opt->StartTimer(PMVSETUP); CHKERRQ(ierr);

    // switch case for choice of preconditioner
    if (this->m_Opt->GetKrylovSolverPara().pctype == TWOLEVEL) {
        // apply restriction to adjoint, state and control variable
        ierr = this->ApplyRestriction(); CHKERRQ(ierr);
    }
    this->m_Opt->PrecondSetupDone(true);

    // stop timer
    ierr = this->m_Opt->StopTimer(PMVSETUP); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief setup phase of preconditioner
 *******************************************************************/
PetscErrorCode PrecondReg::SetupCoarseGrid() {
    PetscErrorCode ierr = 0;
    IntType nt, nc;
    ScalarType scale, value;
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // check if optimization problem is set up
    ierr = Assert(this->m_Opt != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_OptimizationProblem != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_PreProc != NULL, "null pointer"); CHKERRQ(ierr);

    nt  = this->m_Opt->GetDomainPara().nt;
    nc  = this->m_Opt->GetDomainPara().nc;

    // set up options for coarse grid (copy all parameters
    // but the grid resolution and do setup of all plans)
    if (this->m_CoarseGrid.m_Opt != NULL) {
        delete this->m_CoarseGrid.m_Opt;
        this->m_CoarseGrid.m_Opt = NULL;
    }
    try {this->m_CoarseGrid.m_Opt = new RegOpt(*this->m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // get grid scale and compute number of grid points
    scale = this->m_Opt->GetKrylovSolverPara().pcgridscale;
    for (int i = 0; i < 3; ++i) {
        value = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[i])/scale;
        this->m_CoarseGrid.m_Opt->SetNumGridPoints(i,static_cast<IntType>(std::ceil(value)));
    }
    ierr = this->m_CoarseGrid.m_Opt->DoSetup(false); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2) {
        ss  << "setup of preconditioner (data allocation) "
            << "nx (f): (" << this->m_Opt->GetDomainPara().nx[0]
            << "," << this->m_Opt->GetDomainPara().nx[1]
            << "," << this->m_Opt->GetDomainPara().nx[2] << "); "
            << "nx (c): (" << this->m_CoarseGrid.m_Opt->GetDomainPara().nx[0]
            << "," << this->m_CoarseGrid.m_Opt->GetDomainPara().nx[1]
            << "," << this->m_CoarseGrid.m_Opt->GetDomainPara().nx[2] << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // allocate optimization problem
    if (this->m_CoarseGrid.m_OptimizationProblem != NULL) {
        delete this->m_CoarseGrid.m_OptimizationProblem;
        this->m_CoarseGrid.m_OptimizationProblem = NULL;
    }

    // allocate class for registration
    if (this->m_Opt->GetRegModel() == COMPRESSIBLE) {
        try {this->m_CoarseGrid.m_OptimizationProblem = new OptimalControlRegistration(this->m_CoarseGrid.m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    } else if (this->m_Opt->GetRegModel() == STOKES) {
        try {this->m_CoarseGrid.m_OptimizationProblem = new OptimalControlRegistrationIC(this->m_CoarseGrid.m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    } else if (this->m_Opt->GetRegModel() == RELAXEDSTOKES) {
        try {this->m_CoarseGrid.m_OptimizationProblem  = new OptimalControlRegistrationRelaxedIC(this->m_CoarseGrid.m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    } else {
        ierr = ThrowError("registration model not defined"); CHKERRQ(ierr);
    }

    // create vector fields
    ierr = VecCreate(this->m_WorkScaField1, this->m_Opt->GetDomainPara().nl, this->m_Opt->GetDomainPara().ng); CHKERRQ(ierr);
    ierr = VecCreate(this->m_WorkScaField2, this->m_Opt->GetDomainPara().nl, this->m_Opt->GetDomainPara().ng); CHKERRQ(ierr);

    try {this->m_ControlVariable = new VecField(this->m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {this->m_IncControlVariable = new VecField(this->m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ierr = VecCreate(this->m_CoarseGrid.m_StateVariable, (nt+1)*nc*this->m_CoarseGrid.nl(), (nt+1)*nc*this->m_CoarseGrid.ng()); CHKERRQ(ierr);
    ierr = VecCreate(this->m_CoarseGrid.m_AdjointVariable, (nt+1)*nc*this->m_CoarseGrid.nl(), (nt+1)*nc*this->m_CoarseGrid.ng()); CHKERRQ(ierr);

    ierr = VecCreate(this->m_CoarseGrid.m_WorkScaField1, this->m_CoarseGrid.nl(), this->m_CoarseGrid.ng()); CHKERRQ(ierr);
    ierr = VecCreate(this->m_CoarseGrid.m_WorkScaField2, this->m_CoarseGrid.nl(), this->m_CoarseGrid.ng()); CHKERRQ(ierr);

    try {this->m_CoarseGrid.m_ControlVariable = new VecField(this->m_CoarseGrid.m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {this->m_CoarseGrid.m_IncControlVariable = new VecField(this->m_CoarseGrid.m_Opt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ierr = VecCreate(this->m_CoarseGrid.x, 3*this->m_CoarseGrid.nl(), 3*this->m_CoarseGrid.ng()); CHKERRQ(ierr);
    ierr = VecCreate(this->m_CoarseGrid.y, 3*this->m_CoarseGrid.nl(), 3*this->m_CoarseGrid.ng()); CHKERRQ(ierr);

    this->m_CoarseGrid.setupdone = true;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
PetscErrorCode PrecondReg::MatVec(Vec Px, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // switch case for choice of preconditioner
    switch(this->m_Opt->GetKrylovSolverPara().pctype) {
        case NOPC:
        {
            ierr = WrngMsg("no preconditioner used"); CHKERRQ(ierr);
            ierr = VecCopy(x,Px); CHKERRQ(ierr);
            break;
        }
        case INVREG:
        {
            ierr = this->ApplyInvRegPrecond(Px, x); CHKERRQ(ierr);
            break;
        }
        case TWOLEVEL:
        {
            ierr = this->Apply2LevelPrecond(Px, x); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("preconditioner not defined"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply inverse of regularization operator as preconditioner
 *******************************************************************/
PetscErrorCode PrecondReg::ApplyInvRegPrecond(Vec Px, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // check if optimization problem is set up
    ierr = Assert(this->m_OptimizationProblem != NULL, "null pointer"); CHKERRQ(ierr);

    // start timer
    ierr = this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);

    // apply inverse regularization operator
    ierr = this->m_OptimizationProblem->ApplyInvRegOp(Px, x); CHKERRQ(ierr);

    // stop timer
    ierr = this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PCMATVEC);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
PetscErrorCode PrecondReg::Apply2LevelPrecond(Vec Px, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    ScalarType pct, value;
    this->m_Opt->Enter(__func__);

    // do allocation of coarse grid
    if (!this->m_CoarseGrid.setupdone) {
        ierr = this->SetupCoarseGrid(); CHKERRQ(ierr);
    }

    // do setup
    if (this->m_KrylovMethod == NULL) {
        ierr = this->SetupKrylovMethod(this->m_CoarseGrid.nl(), this->m_CoarseGrid.ng()); CHKERRQ(ierr);
    }

    // check if all the necessary pointers have been initialized
    ierr = Assert(this->m_PreProc != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_CoarseGrid.x != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_CoarseGrid.y  != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_CoarseGrid.m_Opt != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncControlVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_CoarseGrid.m_IncControlVariable != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate vector field
    if (this->m_WorkVecField == NULL) {
        try{ this->m_WorkVecField = new VecField(this->m_Opt); }
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    pct = 0; // set to zero, cause we search for a max
    for (int i = 0; i < 3; ++i) {
        value = static_cast<ScalarType>(this->m_CoarseGrid.m_Opt->GetDomainPara().nx[i])
                            /static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[i]);
        pct = value > pct ? value : pct;
    }

    // set components
    ierr = this->m_WorkVecField->SetComponents(x); CHKERRQ(ierr);

    // apply low pass filter before we restrict
    // incremental control variable to coarse grid
    ierr = this->m_PreProc->ApplyRectFreqFilter(this->m_IncControlVariable,
                                                this->m_WorkVecField, pct); CHKERRQ(ierr);

    // apply restriction operator to incremental control variable
    ierr = this->m_PreProc->Restrict(this->m_CoarseGrid.m_IncControlVariable,
                                     this->m_IncControlVariable,
                                     this->m_CoarseGrid.m_Opt->GetDomainPara().nx,
                                     this->m_Opt->GetDomainPara().nx); CHKERRQ(ierr);

    // get the components to interface hessian mat vec
    ierr = this->m_CoarseGrid.m_IncControlVariable->GetComponents(this->m_CoarseGrid.x); CHKERRQ(ierr);


    // invert preconditioner
    ierr = this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);
    ierr = KSPSolve(this->m_KrylovMethod, this->m_CoarseGrid.x, this->m_CoarseGrid.y); CHKERRQ(ierr);
    ierr = this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);


    // get components (for interface of hessian matvec)
    ierr = this->m_CoarseGrid.m_IncControlVariable->SetComponents(this->m_CoarseGrid.y); CHKERRQ(ierr);

    // apply prolongation operator
    ierr = this->m_PreProc->Prolong(this->m_IncControlVariable,
                                    this->m_CoarseGrid.m_IncControlVariable,
                                    this->m_Opt->GetDomainPara().nx,
                                    this->m_CoarseGrid.m_Opt->GetDomainPara().nx); CHKERRQ(ierr);

    // apply low pass filter to output of hessian matvec
    ierr = this->m_PreProc->ApplyRectFreqFilter(this->m_IncControlVariable,
                                                this->m_IncControlVariable, pct); CHKERRQ(ierr);

    // apply high-pass filter to input
    ierr = this->m_PreProc->ApplyRectFreqFilter(this->m_WorkVecField,
                                                this->m_WorkVecField, pct, false); CHKERRQ(ierr);

    // add up high and low frequency components
    this->m_IncControlVariable->AXPY(1.0, this->m_WorkVecField); CHKERRQ(ierr);

    // parse to output
    ierr = this->m_IncControlVariable->GetComponents(Px); CHKERRQ(ierr);


    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the restriction operator to the state, adjoint,
 * and control variable (setup phase of 2level preconditioner)
 *******************************************************************/
PetscErrorCode PrecondReg::ApplyRestriction() {
    PetscErrorCode ierr = 0;
    IntType nl_f, nl_c, nt, nc, l_f, l_c, lnext_f;
    std::stringstream ss;
    Vec m = NULL, lambda = NULL;
    ScalarType *p_mj = NULL, *p_m = NULL, *p_mjcoarse = NULL, *p_mcoarse = NULL,
                *p_lj = NULL, *p_l = NULL, *p_ljcoarse = NULL, *p_lcoarse = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // check if optimization problem is set up
    ierr = Assert(this->m_OptimizationProblem != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_PreProc != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_CoarseGrid.m_Opt != NULL, "null pointer"); CHKERRQ(ierr);

    nt  = this->m_Opt->GetDomainPara().nt;
    nc  = this->m_Opt->GetDomainPara().nc;

    if (this->m_Opt->GetVerbosity() > 1) {
        ss  << "applying restriction to variables "
            << " (" << this->m_Opt->GetDomainPara().nx[0]
            << ","  << this->m_Opt->GetDomainPara().nx[1]
            << ","  << this->m_Opt->GetDomainPara().nx[2] << ") ->"
            << " (" << this->m_CoarseGrid.m_Opt->GetDomainPara().nx[0]
            << ","  << this->m_CoarseGrid.m_Opt->GetDomainPara().nx[1]
            << ","  << this->m_CoarseGrid.m_Opt->GetDomainPara().nx[2] << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // if parameter continuation is enabled, parse regularization weight
    if (this->m_Opt->GetParaCont().enabled) {
        this->m_CoarseGrid.m_Opt->SetRegularizationWeight(0, this->m_Opt->GetRegNorm().beta[0]);
        this->m_CoarseGrid.m_Opt->SetRegularizationWeight(1, this->m_Opt->GetRegNorm().beta[1]);
    }

    // get variables from optimization problem on fine level
    ierr = this->m_OptimizationProblem->GetControlVariable(this->m_ControlVariable); CHKERRQ(ierr);
    ierr = this->m_OptimizationProblem->GetStateVariable(m); CHKERRQ(ierr);
    ierr = this->m_OptimizationProblem->GetAdjointVariable(lambda); CHKERRQ(ierr);

    // restrict control variable
    ierr = this->m_PreProc->Restrict(this->m_CoarseGrid.m_ControlVariable,
                                     this->m_ControlVariable,
                                     this->m_CoarseGrid.m_Opt->GetDomainPara().nx,
                                     this->m_Opt->GetDomainPara().nx); CHKERRQ(ierr);

    ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(lambda, &p_l); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_CoarseGrid.m_StateVariable, &p_mcoarse); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_CoarseGrid.m_AdjointVariable, &p_lcoarse); CHKERRQ(ierr);

    nl_c = this->m_CoarseGrid.m_Opt->GetDomainPara().nl;
    nl_f = this->m_Opt->GetDomainPara().nl;

    // apply restriction operator to time series of images
    for (IntType j = 0; j <= nt; ++j) {  // for all time points
        for (IntType k = 0; k < nc; ++k) {  // for all components
            l_f = j*nl_f*nc + k*nl_f;
            lnext_f = j*nl_f*nc + (k+1)*nl_f;

            // get time point of state variable on fine grid
            ierr = VecGetArray(this->m_WorkScaField1, &p_mj); CHKERRQ(ierr);
            try {std::copy(p_m+l_f, p_m+lnext_f, p_mj); }
            catch (std::exception&) {
                ierr = ThrowError("copy failed"); CHKERRQ(ierr);
            }
            ierr = VecRestoreArray(this->m_WorkScaField1, &p_mj); CHKERRQ(ierr);

            // apply restriction operator to m_j
            ierr = this->m_PreProc->Restrict(&this->m_CoarseGrid.m_WorkScaField1,
                                              this->m_WorkScaField1,
                                              this->m_CoarseGrid.m_Opt->GetDomainPara().nx,
                                              this->m_Opt->GetDomainPara().nx); CHKERRQ(ierr);

            // store restricted state variable
            l_c = j*nl_c*nc + k*nl_c;
            ierr = VecGetArray(this->m_CoarseGrid.m_WorkScaField1, &p_mjcoarse); CHKERRQ(ierr);
            try {std::copy(p_mjcoarse, p_mjcoarse+nl_c, p_mcoarse+l_c);}
            catch (std::exception&) {
                ierr = ThrowError("copy failed"); CHKERRQ(ierr);
            }
            ierr = VecRestoreArray(this->m_CoarseGrid.m_WorkScaField1, &p_mjcoarse); CHKERRQ(ierr);

            // get time point of adjoint variable on fine grid
            ierr = VecGetArray(this->m_WorkScaField2, &p_lj); CHKERRQ(ierr);
            try {std::copy(p_l+l_f, p_l+lnext_f, p_lj);}
            catch(std::exception&) {
                ierr = ThrowError("copy failed"); CHKERRQ(ierr);
            }
            ierr = VecRestoreArray(this->m_WorkScaField2, &p_lj); CHKERRQ(ierr);

            // apply restriction operator
            ierr = this->m_PreProc->Restrict(&this->m_CoarseGrid.m_WorkScaField2,
                                              this->m_WorkScaField2,
                                              this->m_CoarseGrid.m_Opt->GetDomainPara().nx,
                                              this->m_Opt->GetDomainPara().nx); CHKERRQ(ierr);

            // store restricted adjoint variable
            ierr = VecGetArray(this->m_CoarseGrid.m_WorkScaField2, &p_ljcoarse); CHKERRQ(ierr);
            try {std::copy(p_ljcoarse, p_ljcoarse+nl_c, p_lcoarse+l_c);}
            catch (std::exception&) {
                ierr = ThrowError("copy failed"); CHKERRQ(ierr);
            }
            ierr = VecRestoreArray(this->m_CoarseGrid.m_WorkScaField2, &p_ljcoarse); CHKERRQ(ierr);
        }  // for all components
    }  // for all time points

    ierr = VecRestoreArray(this->m_CoarseGrid.m_AdjointVariable, &p_lcoarse); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_CoarseGrid.m_StateVariable, &p_mcoarse); CHKERRQ(ierr);
    ierr = VecRestoreArray(lambda, &p_l); CHKERRQ(ierr);
    ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);

    // parse variables to optimization problem on coarse level (we have to set the control
    // variable first)
    ierr = this->m_CoarseGrid.m_OptimizationProblem->SetControlVariable(this->m_CoarseGrid.m_ControlVariable); CHKERRQ(ierr);
    ierr = this->m_CoarseGrid.m_OptimizationProblem->SetStateVariable(this->m_CoarseGrid.m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_CoarseGrid.m_OptimizationProblem->SetAdjointVariable(this->m_CoarseGrid.m_AdjointVariable); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup for krylov method
 *******************************************************************/
PetscErrorCode PrecondReg::SetupKrylovMethod(IntType nl, IntType ng) {
    PetscErrorCode ierr = 0;
    PC pc = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_KrylovMethod == NULL, "expecting null pointer"); CHKERRQ(ierr);
    ierr = KSPCreate(PETSC_COMM_WORLD, &this->m_KrylovMethod); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("preconditioner: setup krylovmethod"); CHKERRQ(ierr);
    }

    switch (this->m_Opt->GetKrylovSolverPara().pcsolver) {
        case CHEB:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("preconditioner: semi-iterative chebyshev method selected"); CHKERRQ(ierr);
            }
            // chebyshev iteration
            ierr = KSPSetType(this->m_KrylovMethod, KSPCHEBYSHEV); CHKERRQ(ierr);
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
            ierr = KSPChebyshevEstEigSetUseRandom(this->m_KrylovMethod, PETSC_TRUE); CHKERRQ(ierr);
#endif
            break;
        }
        case PCG:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("preconditioner: cg selected"); CHKERRQ(ierr);
            }
            // preconditioned conjugate gradient
            ierr = KSPSetType(this->m_KrylovMethod, KSPCG); CHKERRQ(ierr);
            break;
        }
        case FCG:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("preconditioner: flexible cg selected"); CHKERRQ(ierr);
            }
            // flexible conjugate gradient
            ierr = KSPSetType(this->m_KrylovMethod, KSPFCG); CHKERRQ(ierr);
            break;
        }
        case GMRES:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("preconditioner: gmres selected"); CHKERRQ(ierr);
            }
            // generalized minimal residual method
            ierr = KSPSetType(this->m_KrylovMethod, KSPGMRES); CHKERRQ(ierr);
            break;
        }
        case FGMRES:
        {
            if (this->m_Opt->GetVerbosity() > 2) {
                ierr = DbgMsg("preconditioner: flexible gmres selected"); CHKERRQ(ierr);
            }
            // flexible generalized minimal residual method
            ierr = KSPSetType(this->m_KrylovMethod, KSPFGMRES); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }

    //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
    //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
    //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
    ierr = KSPSetNormType(this->m_KrylovMethod, KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    //ierr = KSPSetNormType(this->m_KrylovMethod,KSP_NORM_PRECONDITIONED); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethod, PETSC_FALSE); CHKERRQ(ierr);

    //ierr = KSPSetPostSolve(this->m_KrylovMethod,PostKrylovSolve,this);
    ierr = KSPSetPreSolve(this->m_KrylovMethod, InvertPrecondPreKrylovSolve, this); CHKERRQ(ierr);

    // set up matvec for preconditioner
    if (this->m_MatVec != NULL) {
        ierr = MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    ierr = MatCreateShell(PETSC_COMM_WORLD, 3*nl, 3*nl, 3*ng, 3*ng, this, &this->m_MatVec); CHKERRQ(ierr);
    ierr = MatShellSetOperation(this->m_MatVec, MATOP_MULT, (void(*)(void))InvertPrecondMatVec); CHKERRQ(ierr);
    ierr = KSPSetOperators(this->m_KrylovMethod, this->m_MatVec, this->m_MatVec);CHKERRQ(ierr);

    // TODO: make sure this make sense; we will have to switch the hessian
    // operator, which currently is not done
//   if (     (this->m_Opt->GetKrylovSolverPara().pcsolver == GMRES)
//        ||  (this->m_Opt->GetKrylovSolverPara().pcsolver == FGMRES) ) {
//        ierr = MatSetOption(this->m_MatVec, MAT_SYMMETRIC, PETSC_FALSE); CHKERRQ(ierr);
//    }

    ierr = MatSetOption(this->m_MatVec, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(ierr);
    //ierr = MatSetOption(this->m_MatVec,MAT_SYMMETRIC,PETSC_FALSE); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 1) {
        ierr = KSPMonitorSet(this->m_KrylovMethod, InvertPrecondKrylovMonitor, this, NULL); CHKERRQ(ierr);
    }

    // remove preconditioner (the hessian is a preconditioned
    // representation
    ierr = KSPGetPC(this->m_KrylovMethod, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCNONE); CHKERRQ(ierr);  ///< set no preconditioner

    // finish
    ierr = KSPSetFromOptions(this->m_KrylovMethod); CHKERRQ(ierr);
    ierr = KSPSetUp(this->m_KrylovMethod); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup for two level preconditioner
 *******************************************************************/
PetscErrorCode PrecondReg::HessianMatVec(Vec Hx, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // apply hessian (hessian matvec)
    if (this->m_Opt->GetKrylovSolverPara().pctype == TWOLEVEL) {
        if (this->m_Opt->GetVerbosity() > 2) {
            ierr = DbgMsg("preconditioner: (H^coarse + Q^coarse)[x^coarse]"); CHKERRQ(ierr);
        }
        ierr = this->m_CoarseGrid.m_OptimizationProblem->HessianMatVec(Hx, x, false); CHKERRQ(ierr);
    } else {
        ierr = this->m_OptimizationProblem->HessianMatVec(Hx, x); CHKERRQ(ierr);
    }

    // increment counter
    this->m_Opt->IncrementCounter(PCMATVEC);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief this is an interface to compute the eigenvalues needed
 * when considering a chebyshev method to invert the preconditioner;
 * the eigenvalues are estimated using the Lanczo (KSPCG) or
 * Arnoldi (KSPGMRES) process using a random right hand side vector
 *******************************************************************/
PetscErrorCode PrecondReg::EstimateEigenValues() {
    PetscErrorCode ierr = 0;
    IntType i, n, neig, nl, ng;
    std::stringstream ss;
    Vec b = NULL, x = NULL;
    ScalarType *re = NULL, *im = NULL, eigmin, eigmax, emin, emax;

    PetscFunctionBegin;

/*
    itermax = 20;

    // get iteration number; if we need to many iterations, we
    // might want to reestimate the eigenvalues
    n = this->m_Opt->GetKrylovSolverPara().iter;
    if ( n > itermax ) {

        ss << "iter="<< n << " > itermax=" << itermax <<" re-estimating eigenvalues";
        ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
        this->m_EigenValuesEstimated=false;

        // reset krylov iterations
        this->m_Opt->SetKrylovIterations(0);

    }
*/
    if (!this->m_Opt->GetKrylovSolverPara().eigvalsestimated) {
        if (this->m_Opt->GetKrylovSolverPara().usepetsceigest) {
            if (this->m_Opt->GetVerbosity() > 1) {
                ierr = DbgMsg("estimating eigenvalues (petsc)"); CHKERRQ(ierr);
            }
            // default interface for chebyshev method to estimate eigenvalues
            ierr = KSPChebyshevEstEigSet(this->m_KrylovMethod, PETSC_DECIDE,
                                                               PETSC_DECIDE,
                                                               PETSC_DECIDE,
                                                               PETSC_DECIDE); CHKERRQ(ierr);
        } else {
            if (this->m_Opt->GetVerbosity() > 1) {
                ierr = DbgMsg("estimating eigenvalues"); CHKERRQ(ierr);
            }
            // get sizes
            nl = this->m_Opt->GetDomainPara().nl;
            ng = this->m_Opt->GetDomainPara().ng;

            ierr = VecCreate(x, 3*nl, 3*ng); CHKERRQ(ierr);
            ierr = VecCreate(b, 3*nl, 3*ng); CHKERRQ(ierr);

            // use random right hand side
            if (this->m_RandomNumGen == NULL) {
                ierr = PetscRandomCreate(PetscObjectComm((PetscObject)b), &this->m_RandomNumGen); CHKERRQ(ierr);
            }
            ierr = VecSetRandom(b, this->m_RandomNumGen); CHKERRQ(ierr);

            // do setup
            if (this->m_KrylovMethodEigEst == NULL) {
                ierr = this->SetupKrylovMethodEigEst(); CHKERRQ(ierr);
            }
            ierr = Assert(this->m_KrylovMethodEigEst != NULL, "null pointer"); CHKERRQ(ierr);

            ierr = KSPSolve(this->m_KrylovMethodEigEst, b, x); CHKERRQ(ierr);
            ierr = KSPGetIterationNumber(this->m_KrylovMethodEigEst, &n); CHKERRQ(ierr);

            ierr = PetscMalloc2(n, &re, n, &im); CHKERRQ(ierr);
            ierr = KSPComputeEigenvalues(this->m_KrylovMethodEigEst, n, re, im, &neig); CHKERRQ(ierr);

            eigmin = PETSC_MAX_REAL;
            eigmax = PETSC_MIN_REAL;

            for (i = 0; i < neig; ++i) {
                eigmin = PetscMin(eigmin, re[i]);
                eigmax = PetscMax(eigmax, re[i]);
            }

            // clear memory
            ierr = PetscFree2(re, im); CHKERRQ(ierr);
            ierr = this->m_CoarseGrid.m_OptimizationProblem->EstimateExtremalHessEigVals(emin, emax); CHKERRQ(ierr);

            ierr = KSPChebyshevSetEigenvalues(this->m_KrylovMethod, eigmax, eigmin); CHKERRQ(ierr);
        }   // switch between eigenvalue estimators

        // set flag
        this->m_Opt->KrylovMethodEigValsEstimated(true);
    }

    if (x != NULL) {ierr = VecDestroy(&x); CHKERRQ(ierr);}
    if (b != NULL) {ierr = VecDestroy(&b); CHKERRQ(ierr);}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief do setup for krylov method to estimate eigenvalues
 *******************************************************************/
PetscErrorCode PrecondReg::SetupKrylovMethodEigEst() {
    PetscErrorCode ierr = 0;
    PC pc = NULL;
    IntType nl, ng;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // get sizes
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    // create krylov method
    if (this->m_KrylovMethodEigEst != NULL) {
        ierr = KSPDestroy(&this->m_KrylovMethodEigEst); CHKERRQ(ierr);
        this->m_KrylovMethodEigEst = NULL;
    }
    ierr = KSPCreate(PETSC_COMM_WORLD,&this->m_KrylovMethodEigEst); CHKERRQ(ierr);

    // preconditioned conjugate gradient
    ierr = KSPSetType(this->m_KrylovMethodEigEst, KSPCG); CHKERRQ(ierr);

    //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
    //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
    //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
    ierr = KSPSetNormType(this->m_KrylovMethodEigEst, KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    //ierr = KSPSetNormType(this->m_KrylovMethod,KSP_NORM_PRECONDITIONED); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethodEigEst, PETSC_FALSE); CHKERRQ(ierr);

    // set up matvec for preconditioner
    if (this->m_MatVecEigEst != NULL) {
        ierr = MatDestroy(&this->m_MatVecEigEst); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    ierr = MatCreateShell(PETSC_COMM_WORLD, 3*nl, 3*nl, 3*ng, 3*ng, this, &this->m_MatVecEigEst); CHKERRQ(ierr);
    ierr = MatShellSetOperation(this->m_MatVecEigEst, MATOP_MULT, (void(*)(void))InvertPrecondMatVec); CHKERRQ(ierr);
    ierr = KSPSetOperators(this->m_KrylovMethodEigEst, this->m_MatVecEigEst, this->m_MatVecEigEst); CHKERRQ(ierr);
    ierr = MatSetOption(this->m_MatVecEigEst, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(ierr);
    //ierr = MatSetOption(this->m_MatVec,MAT_SYMMETRIC,PETSC_FALSE); CHKERRQ(ierr);

    // remove preconditioner
    ierr = KSPGetPC(this->m_KrylovMethodEigEst, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCNONE); CHKERRQ(ierr); ///< set no preconditioner

    ierr = KSPSetTolerances(this->m_KrylovMethodEigEst, 1E-12, 1E-12, 1E+6, 10); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethodEigEst, PETSC_FALSE); CHKERRQ(ierr);
    ierr = KSPAppendOptionsPrefix(this->m_KrylovMethodEigEst, "esteig_"); CHKERRQ(ierr);

    // we are going to estimate eigenvalues with this
    ierr = KSPSetComputeEigenvalues(this->m_KrylovMethodEigEst, PETSC_TRUE); CHKERRQ(ierr);

    // finish
    ierr = KSPSetUp(this->m_KrylovMethodEigEst); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}   // namespace reg




#endif   // _PRECONDREG_CPP_
