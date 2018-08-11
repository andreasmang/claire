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

#ifndef _CLAIRE_CPP_
#define _CLAIRE_CPP_

// global includes
#include <string>
#include <algorithm>

// local includes
#include "CLAIRE.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
CLAIRE::CLAIRE() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
CLAIRE::CLAIRE(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
CLAIRE::~CLAIRE() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode CLAIRE::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_StateVariable = NULL;       ///< state variable
    this->m_AdjointVariable = NULL;     ///< adjoint variable
    this->m_IncStateVariable = NULL;    ///< incremental state variable
    this->m_IncAdjointVariable = NULL;  ///< incremental adjoint variable

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode CLAIRE::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // delete all variables
    ierr = this->ClearVariables(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode CLAIRE::ClearVariables(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // delete all variables
    if (this->m_StateVariable != NULL) {
        ierr = VecDestroy(&this->m_StateVariable); CHKERRQ(ierr);
        this->m_StateVariable = NULL;
    }
    if (this->m_AdjointVariable != NULL) {
        ierr = VecDestroy(&this->m_AdjointVariable); CHKERRQ(ierr);
        this->m_AdjointVariable = NULL;
    }
    if (this->m_IncStateVariable != NULL) {
        ierr = VecDestroy(&this->m_IncStateVariable); CHKERRQ(ierr);
        this->m_IncStateVariable = NULL;
    }
    if (this->m_IncAdjointVariable != NULL) {
        ierr = VecDestroy(&this->m_IncAdjointVariable); CHKERRQ(ierr);
        this->m_IncAdjointVariable = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief setup solver (to not have to allocate things on the
 * fly; this allows us to essentially do a warm start)
 *******************************************************************/
PetscErrorCode CLAIRE::InitializeSolver(void) {
    PetscErrorCode ierr = 0;
    IntType nt, nl, nc, ng;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    if (this->m_StateVariable == NULL) {
        ierr = VecCreate(this->m_StateVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        if (this->m_AdjointVariable == NULL) {
            ierr = VecCreate(this->m_AdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
        }
        if (this->m_IncAdjointVariable == NULL) {
            ierr = VecCreate(this->m_IncAdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
        }
        if (this->m_IncStateVariable == NULL) {
            ierr = VecCreate(this->m_IncStateVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
        }
    } else {
        if (this->m_AdjointVariable == NULL) {
            ierr = VecCreate(this->m_AdjointVariable, nc*nl, nc*ng); CHKERRQ(ierr);
        }
        if (this->m_IncAdjointVariable == NULL) {
            ierr = VecCreate(this->m_IncAdjointVariable, nc*nl, nc*ng); CHKERRQ(ierr);
        }
        if (this->m_IncStateVariable == NULL) {
            ierr = VecCreate(this->m_IncStateVariable, nc*nl, nc*ng); CHKERRQ(ierr);
        }
    }

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL) {
        ierr = VecCreate(this->m_WorkScaField3, nl, ng); CHKERRQ(ierr);
    }

    if (this->m_Opt->m_PDESolver.type == SL) {
        if (this->m_SemiLagrangianMethod == NULL) {
            try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }

    if (this->m_Differentiation == NULL) {
        try {this->m_Differentiation = new DifferentiationSM(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }


    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief initialize the optimization (we essentially evaluate
 * the objective functional and the gradient for a given initial
 * guess)
 *******************************************************************/
PetscErrorCode CLAIRE::InitializeOptimization() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    std::stringstream ss;
    ScalarType value, alpha, jvt, jv, lsred, descent;
    Vec g = NULL, dv = NULL, v = NULL, vtilde = NULL;
    bool lssuccess, restoreinitialguess = false;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // if velocity field is null pointer, we did not set
    // any initial guess
    if (this->m_VelocityField == NULL) {
        if (this->m_Opt->m_Verbosity > 2) {
            ierr = DbgMsg("allocating velocity field"); CHKERRQ(ierr);
        }
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    } else { // user might have provided initial guess
        ierr = this->IsVelocityZero(); CHKERRQ(ierr);
        if (!this->m_VelocityIsZero) {
            restoreinitialguess = true;
            if (this->m_WorkVecField1 == NULL) {
                try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
                catch (std::bad_alloc& err) {
                    ierr = reg::ThrowError(err); CHKERRQ(ierr);
                }
            }
            ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
        }
    }
    ierr = VecCreate(g, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecCreate(v, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecSet(v, 0.0); CHKERRQ(ierr);

    // if we use a non-zero initial guess, we compute
    // the first velocity using a steepest descent approach
    if (!this->m_Opt->m_OptPara.usezeroinitialguess) {
        ierr = VecCreate(dv, 3*nl, 3*ng); CHKERRQ(ierr);
        ierr = VecCreate(vtilde, 3*nl, 3*ng); CHKERRQ(ierr);

        lsred = 1E-4;  // reduction rate for line search
        for (int l = 0; l < 1; ++l) {
            // evaluate objective function
            ierr = this->EvaluateObjective(&jv, v); CHKERRQ(ierr);

            // compute gradient
            ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);

            // compute search direction (gradient in sobolev space)
            ierr = this->EvaluateGradient(dv, v); CHKERRQ(ierr);

            // inner product between gradient and search direction
            ierr = VecTDot(g, dv, &descent); CHKERRQ(ierr);

            alpha = 1.0; lssuccess = false;
            for (int i = 0; i < 20; ++i) {
                // compute trial velocity
                ierr = VecWAXPY(vtilde, -alpha, dv, v); CHKERRQ(ierr);

                // evaluate objective function
                ierr = this->EvaluateObjective(&jvt, vtilde); CHKERRQ(ierr);

                // armijo rule
                if (jvt < jv + alpha*lsred*descent) {
                    lssuccess = true;
                    break;
                }
                alpha /= 2.0;
            }
            if (lssuccess) {
                if (this->m_Opt->m_Verbosity > 1) {
                    ss << "line search successful (initialization; alpha=" << std::scientific << alpha << ")";
                    ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                    ss.clear(); ss.str(std::string());
                }
                ierr = VecCopy(vtilde, v); CHKERRQ(ierr);
            } else {
                if (this->m_Opt->m_Verbosity > 1) {
                    ss << "line search failed (initialization; alpha=" << std::scientific << alpha << ")";
                    ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                    ss.clear(); ss.str(std::string());
                }
                break;
            }
        }
    }
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // evaluate distance measure
    ierr = this->EvaluateDistanceMeasure(&value); CHKERRQ(ierr);
    this->m_Opt->m_Monitor.dval0 = value;

    // evaluate objective functional
    ierr = this->EvaluateObjective(&value, v); CHKERRQ(ierr);
    this->m_Opt->m_Monitor.jval0 = value;

    // compute gradient
    ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);

    // compute gradient norm
    ierr = VecNorm(g, NORM_2, &value); CHKERRQ(ierr);
    this->m_Opt->m_Monitor.gradnorm0 = value;

    if (this->m_Opt->m_Verbosity > 0) {
        ss << "initial gradient norm: "<< std::scientific << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    // if we had a non-zero initial velocity, we'll restore it
    if (restoreinitialguess) {
        ierr = this->m_VelocityField->Copy(this->m_WorkVecField1); CHKERRQ(ierr);
    }

    // clean up
    if (v != NULL) {ierr = VecDestroy(&v); CHKERRQ(ierr); v = NULL;}
    if (g != NULL) {ierr = VecDestroy(&g); CHKERRQ(ierr); g = NULL;}
    if (dv != NULL) {ierr = VecDestroy(&dv); CHKERRQ(ierr); dv = NULL;}
    if (vtilde != NULL) {ierr = VecDestroy(&vtilde); CHKERRQ(ierr); vtilde = NULL;}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set m at t=0
 * @param[in] m0 density/image at t=0 (initial condition of forward
 * problem)
 *******************************************************************/
PetscErrorCode CLAIRE::SetInitialState(Vec m0) {
    PetscErrorCode ierr = 0;
    ScalarType *p_m0 = NULL, *p_m = NULL;
    IntType nt, nl, nc, ng;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m0 != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // allocate state variable
    if (this->m_StateVariable == NULL) {
        if (this->m_Opt->m_RegFlags.runinversion) {
            ierr = VecCreate(this->m_StateVariable, (nt+1)*nl*nc, (nt+1)*ng*nc); CHKERRQ(ierr);
        } else {
            ierr = VecCreate(this->m_StateVariable, nl*nc, ng*nc); CHKERRQ(ierr);
        }
    }

    // copy m_0 to m(t=0)
    ierr = GetRawPointer(m0, &p_m0); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    try {std::copy(p_m0, p_m0+nl*nc, p_m);}
    catch (std::exception& err) {
        ierr = ThrowError(err); CHKERRQ(ierr);
    }
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointer(m0, &p_m0); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get m at t=1
 * @param[out] m1 density/image at t=1 (solution of forward problem)
 *******************************************************************/
PetscErrorCode CLAIRE::GetFinalState(Vec m1) {
    PetscErrorCode ierr = 0;
    ScalarType *p_m1 = NULL, *p_m = NULL;
    IntType nt, nl, nc;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    if (!this->m_Opt->m_RegFlags.runinversion) {
        nt = 0; // we did not store the time history
    }

    // copy m(t=1) to m_1
    ierr = GetRawPointer(m1, &p_m1); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
    catch (std::exception& err) {
        ierr = ThrowError(err); CHKERRQ(ierr);
    }
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointer(m1, &p_m1); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set lambda at t=1
 * @param[out] l1 adjoint variable at t=1 (final condition of adjoint
 * problem)
 *******************************************************************/
PetscErrorCode CLAIRE::SetFinalAdjoint(Vec l1) {
    PetscErrorCode ierr = 0;
    ScalarType *p_l1 = NULL, *p_l = NULL;
    IntType nt, nl, ng, nc;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(l1 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // we do not store the time history for a gauss-newton approximation
    if (this->m_Opt->m_OptPara.method == GAUSSNEWTON) {
        nt = 0;
    }

    // allocate pointer if not done so already
    if (this->m_AdjointVariable == NULL) {
        ierr = VecCreate(this->m_AdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }

    // copy l1 to lambda(t=1)
    ierr = GetRawPointer(l1, &p_l1); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    try {std::copy(p_l+nt*nl*nc, p_l+(nt+1)*nl*nc, p_l1);}
    catch (std::exception& err) {
        ierr = ThrowError(err); CHKERRQ(ierr);
    }
    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = RestoreRawPointer(l1, &p_l1); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem (we assume the user has
 * set the velocity field)
 * @param[in] m0 density/image at t=0
 * @param[out] m1 density/image at t=1 (solution of transport
 * equation)
 *******************************************************************/
PetscErrorCode CLAIRE::SolveForwardProblem(Vec m1, Vec m0) {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m0 != NULL, "null pointer"); CHKERRQ(ierr);

    // set initial condition
    this->m_TemplateImage = m0;

    // compute solution of state equation
    ierr = this->SolveStateEquation(); CHKERRQ(ierr);

    // only copy if output is necessary
    if (m1 != NULL) {
        ierr = this->GetFinalState(m1); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem (we assume the user has
 * set the velocity field)
 * @param[in] m1 density/image at t=1
 * @param[out] l0 adjoint variable at t=0 (solution of transport
 * equation)
 *******************************************************************/
PetscErrorCode CLAIRE::SolveAdjointProblem(Vec l0, Vec m1) {
    PetscErrorCode ierr = 0;
    ScalarType *p_m = NULL, *p_m1 = NULL, *p_l = NULL, *p_l0 = NULL;
    IntType nt, nl, nc, ng;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m1 != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // allocate state variable
    if (this->m_StateVariable == NULL) {
        ierr = VecCreate(this->m_StateVariable, (nt+1)*nl*nc, (nt+1)*ng*nc); CHKERRQ(ierr);
        ierr = VecSet(this->m_StateVariable, 0); CHKERRQ(ierr);
    }

    // copy memory for m_1
    ierr = GetRawPointer(m1, &p_m1); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    try {std::copy(p_m1, p_m1+nl*nc, p_m+nt*nl*nc);}
    catch (std::exception& err) {
        ierr = ThrowError(err); CHKERRQ(ierr);
    }
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointer(m1, &p_m1); CHKERRQ(ierr);

    // compute solution of state equation
    ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);

    // copy memory for lambda0
    ierr = GetRawPointer(l0, &p_l0); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    try {std::copy(p_l, p_l+nl*nc, p_l0);}
    catch (std::exception& err) {
        ierr = ThrowError(err); CHKERRQ(ierr);
    }
    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = RestoreRawPointer(l0, &p_l0); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set state variable from externally
 *******************************************************************/
PetscErrorCode CLAIRE::SetStateVariable(Vec m) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // we have to allocate the variable, because we delete it
    // at the end once we're done; since it comes from external
    // we need to make sure that we don't delete the external
    // pointer
    if (this->m_StateVariable == NULL) {
        ierr = VecCreate(this->m_StateVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }
    ierr = VecCopy(m, this->m_StateVariable); CHKERRQ(ierr);

    // if semi lagrangian pde solver is used,
    // we have to initialize it here
    if (this->m_Opt->m_PDESolver.type == SL) {
        ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
        if (this->m_SemiLagrangianMethod == NULL) {
            try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        // compute trajectory
        if (this->m_WorkVecField1 == NULL) {
            try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set state variable from externally
 *******************************************************************/
PetscErrorCode CLAIRE::GetStateVariable(Vec& m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m == NULL, "null pointer expected"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    m = this->m_StateVariable;

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set adjoint variable from externally
 *******************************************************************/
PetscErrorCode CLAIRE::SetAdjointVariable(Vec lambda) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nc = this->m_Opt->m_Domain.nc;
    nt = this->m_Opt->m_Domain.nt;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // we have to allocate the variable, because we delete it
    // at the end once we're done; since it comes from external
    // we need to make sure that we don't delete the external pointer
    if (this->m_AdjointVariable == NULL) {
        if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
            ierr = VecCreate(this->m_AdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
        } else {
            ierr = VecCreate(this->m_AdjointVariable, nc*nl, nc*ng); CHKERRQ(ierr);
        }
    }
    ierr = VecCopy(lambda, this->m_AdjointVariable); CHKERRQ(ierr);

    if (this->m_Opt->m_PDESolver.type == SL) {
        ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
        if (this->m_SemiLagrangianMethod == NULL) {
            try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        // compute trajectory
        if (this->m_WorkVecField1 == NULL) {
            try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief get adjoint variable
 *******************************************************************/
PetscErrorCode CLAIRE::GetAdjointVariable(Vec& lambda) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda == NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    lambda = this->m_AdjointVariable;
//    ierr = VecCopy(this->m_AdjointVariable, lambda); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluate the l2 distance between m_R and m_1
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateDistanceMeasure(ScalarType* D) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // compute solution of state equation
    ierr = this->SolveStateEquation(); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate distance measure
    if (this->m_Opt->m_Distance.reset) {
        ierr = this->SetupDistanceMeasure(); CHKERRQ(ierr);
        this->m_Opt->m_Distance.reset = false;
    }
    if (this->m_DistanceMeasure == NULL) {
        ierr = this->SetupDistanceMeasure(); CHKERRQ(ierr);
    }

    // set state variable
    ierr = this->m_DistanceMeasure->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);

    // evaluate distance measure
    ierr = this->m_DistanceMeasure->EvaluateFunctional(D); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief evaluates the objective value
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateObjective(ScalarType* J, Vec v) {
    PetscErrorCode ierr = 0;
    ScalarType D = 0.0, R = 0.0;
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // start timer
    ierr = this->m_Opt->StartTimer(OBJEXEC); CHKERRQ(ierr);

    // set components of velocity field
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // evaluate the regularization model
    ierr = this->EvaluateDistanceMeasure(&D); CHKERRQ(ierr);

    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (!this->m_VelocityIsZero) {
        // evaluate the regularization model
        if (this->m_WorkVecField1 == NULL) {
            try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        ierr = this->m_Regularization->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_Regularization->EvaluateFunctional(&R, this->m_VelocityField); CHKERRQ(ierr);
    }

    // add up the contributions
    *J = D + R;

    // store for access (e.g., used in coupling)
    this->m_Opt->m_Monitor.jval = *J;
    this->m_Opt->m_Monitor.dval = D;
    this->m_Opt->m_Monitor.rval = R;

    if (this->m_Opt->m_Verbosity > 1) {
        ss << "J(v) = D(v) + R(v) = " << std::scientific
           << this->m_Opt->m_Monitor.dval << " + "
           << this->m_Opt->m_Monitor.rval;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
    }

    // stop timer
    ierr = this->m_Opt->StopTimer(OBJEXEC); CHKERRQ(ierr);

    // increment counter for objective evaluations
    this->m_Opt->IncrementCounter(OBJEVAL);


    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates the reduced gradient of the lagrangian
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateGradient(Vec g, Vec v) {
    PetscErrorCode ierr = 0;
    ScalarType  value, nvx1, nvx2, nvx3;
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg("evaluating gradient"); CHKERRQ(ierr);
    }
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // start timer
    ierr = this->m_Opt->StartTimer(GRADEXEC); CHKERRQ(ierr);

    // parse input arguments
    if (v != NULL) {
        ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);
    }

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = this->m_VelocityField->Norm(nvx1, nvx2, nvx3); CHKERRQ(ierr);
        ss  << "||v||_2 = (" << std::scientific
            << nvx1 << "," << nvx2 << "," << nvx3 << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    // compute solution of adjoint equation (i.e., \lambda(x,t))
    // and compute body force \int_0^1 grad(m)\lambda dt
    // which is assigned to work vecfield 2
    ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);

    // evaluate gradient of regularization model
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        // \vect{g}_v = \D{K}[\vect{b}]
        if (g != NULL) {
            ierr = this->m_WorkVecField2->GetComponents(g); CHKERRQ(ierr);
        }
    } else {
        switch (this->m_Opt->m_OptPara.gradtype) {
            case L2GRAD:
            {
                // evaluate l2 gradient
                ierr = this->EvaluateL2Gradient(g); CHKERRQ(ierr);
                break;
            }
/*
            case SGRAD:
            {
                // evaluate sobolev gradient
                ierr = this->EvaluateSobolevGradient(g, false); CHKERRQ(ierr);
                break;
            }
            case SYMSGRAD:
            {
                // evaluate sobolev gradient
                ierr = this->EvaluateSobolevGradient(g, true); CHKERRQ(ierr);
                break;
            }
*/
            default:
            {
                ierr = ThrowError("operator not implemented"); CHKERRQ(ierr);
                break;
            }
        }
    }


    // parse to output
    if (g != NULL) {
        // get and scale by lebesque measure
//        hd = this->m_Opt->GetLebesgueMeasure();
//        ierr = VecScale(g, hd); CHKERRQ(ierr);

        if (this->m_Opt->m_Verbosity > 2) {
            ierr = VecNorm(g, NORM_2, &value); CHKERRQ(ierr);
            ss << "||g||_2 = " << std::scientific << value;
            ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.clear(); ss.str(std::string());
        }
    }

    // stop timer
    ierr = this->m_Opt->StopTimer(GRADEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(GRADEVAL);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates the reduced gradient of the lagrangian (l2)
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateL2Gradient(Vec g) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // evaluate / apply gradient operator for regularization
    ierr = this->m_Regularization->EvaluateGradient(this->m_WorkVecField1, this->m_VelocityField); CHKERRQ(ierr);

    // \vect{g}_v = \beta_v \D{A}[\vect{v}] + \D{K}[\vect{b}]
    ierr = this->m_WorkVecField1->AXPY(1.0, this->m_WorkVecField2); CHKERRQ(ierr);

    // copy
    ierr = this->m_WorkVecField1->GetComponents(g); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates the reduced gradient of the lagrangian (in
 * sobolev space incuded by regularization operator)
 *******************************************************************/
PetscErrorCode CLAIRE::EvaluateSobolevGradient(Vec g, bool flag) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // evaluate / apply gradient operator for regularization
    ierr = this->m_Regularization->ApplyInverse(this->m_WorkVecField1, this->m_WorkVecField2, flag); CHKERRQ(ierr);

    // \vect{g}_v = \vect{v} + (\beta_v \D{A})^{-1}\D{K}[\vect{b}]
    ierr = this->m_WorkVecField1->AXPY(1.0, this->m_VelocityField); CHKERRQ(ierr);

    // copy to output
    ierr = this->m_WorkVecField1->GetComponents(g); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute the body force
 * b = \int_0^1 grad(m) \lambda dt
 *******************************************************************/
PetscErrorCode CLAIRE::ComputeBodyForce() {
    PetscErrorCode ierr = 0;
    IntType nt, nl, nc, l;
    std::stringstream ss;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    ScalarType *p_mt = NULL, *p_m = NULL, *p_l = NULL,
               *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL,
               *p_b1 = NULL, *p_b2 = NULL, *p_b3 = NULL;
    ScalarType ht, scale, lambda, value;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // check for null pointers
    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    // get problem dimensions and weights
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ht = this->m_Opt->GetTimeStepSize();
    scale = ht;

    ierr = Assert(nt > 0, "nt<=0"); CHKERRQ(ierr);
    ierr = Assert(ht > 0, "ht<=0"); CHKERRQ(ierr);

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // init body force for numerical integration
    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);

    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        // m and \lambda are constant in time
        ierr = GetRawPointer(this->m_TemplateImage, &p_mt); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {  // for all components
            // compute gradient of m
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_mt+k*nl, this->m_Opt->m_FFT.plan, &xyz, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

            // b = \sum_k\int_{\Omega} \lambda_k \grad m_k dt
#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                lambda = p_l[k*nl+i];
                p_b1[i] += lambda*p_gradm1[i]/static_cast<ScalarType>(nc);
                p_b2[i] += lambda*p_gradm2[i]/static_cast<ScalarType>(nc);
                p_b3[i] += lambda*p_gradm3[i]/static_cast<ScalarType>(nc);
            }
        }
}  // omp
        ierr = RestoreRawPointer(this->m_TemplateImage, &p_mt); CHKERRQ(ierr);
    } else {  // non zero velocity field
        ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j) {
            // scaling for trapezoidal rule
            if ((j == 0) || (j == nt)) scale *= 0.5;
            for (IntType k = 0; k < nc; ++k) {  // for all components
                l = j*nl*nc + k*nl;

                // grad(m^j)
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m+l, this->m_Opt->m_FFT.plan, &xyz, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
                // \vect{b}_i += h_d*ht*\lambda^j (\grad m^j)_i
                for (IntType i = 0; i < nl; ++i) {
                    lambda = p_l[l+i];  // get \lambda(x_i,t^j)
                    p_b1[i] += scale*p_gradm1[i]*lambda/static_cast<ScalarType>(nc);
                    p_b2[i] += scale*p_gradm2[i]*lambda/static_cast<ScalarType>(nc);
                    p_b3[i] += scale*p_gradm3[i]*lambda/static_cast<ScalarType>(nc);
                }
            }
}  // omp
            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)) scale *= 2.0;
        }
        // restore arrays
        ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    }  // else zero velocity field

    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);   // adjoint variable for all t^j
    ierr = this->m_WorkVecField1->RestoreArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = this->m_WorkVecField2->Norm(value); CHKERRQ(ierr);
        ss << "||b||_2 = " << std::scientific << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the hessian to a vector
 * @param[in] vtilde incremental velocity field
 * @param[in] scale flag to switch on scaling by lebesgue measure
 * @param[out] Hvtilde hessian applied to vector
 *******************************************************************/
PetscErrorCode CLAIRE::HessianMatVec(Vec Hvtilde, Vec vtilde, bool scale) {
    PetscErrorCode ierr = 0;
   ScalarType hd; //, gamma;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg("computing hessian matvec"); CHKERRQ(ierr);
    }

    ierr = this->m_Opt->StartTimer(HMVEXEC); CHKERRQ(ierr);

    // switch between hessian operators
    switch (this->m_Opt->m_KrylovMethod.matvectype) {
        case DEFAULTMATVEC:
        {
            // apply hessian H to \tilde{v}
            ierr = this->HessMatVec(Hvtilde, vtilde); CHKERRQ(ierr);
            //ierr = VecCopy(vtilde, Hvtilde); CHKERRQ(ierr);
            break;
        }
        case PRECONDMATVEC:
        {
            // apply analytically preconditioned hessian H to \tilde{v}
            ierr = this->PrecondHessMatVec(Hvtilde, vtilde); CHKERRQ(ierr);
            break;
        }
        case PRECONDMATVECSYM:
        {
            // apply analytically preconditioned hessian H to \tilde{v}
            // hessian operator is symmetrized
            ierr = this->PrecondHessMatVecSym(Hvtilde, vtilde); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("operator not implemented"); CHKERRQ(ierr);
            break;
        }
    }


    if (Hvtilde != NULL) {
        // TODO @ Andreas: fix for two-level precond
        // scale by lebesgue measure
        if (scale == false) {
            hd = this->m_Opt->GetLebesgueMeasure();
            ierr = VecScale(Hvtilde, 1/hd); CHKERRQ(ierr);
        }

//        gamma = this->m_Opt->m_KrylovMethod.hessshift;
//        if (gamma > 0.0) {
//            ierr = VecAXPY(Hvtilde, gamma, vtilde); CHKERRQ(ierr);
//        }
    }

    // stop hessian matvec timer
    ierr = this->m_Opt->StopTimer(HMVEXEC); CHKERRQ(ierr);

    // increment matvecs
    this->m_Opt->IncrementCounter(HESSMATVEC);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the hessian to a vector (default way of doing this)
 *******************************************************************/
PetscErrorCode CLAIRE::HessMatVec(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate container for incremental velocity field
    if (this->m_IncVelocityField == NULL) {
        try {this->m_IncVelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // parse input
    if (vtilde != NULL) {
        ierr = this->m_IncVelocityField->SetComponents(vtilde); CHKERRQ(ierr);
    }

    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t)
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // compute incremental body force
//    ierr = this->ComputeIncBodyForce(); CHKERRQ(ierr);

    // apply 2nd variation of regularization model to
    // incremental control variable: \beta*\D{A}[\vect{\tilde{v}}]
    ierr = this->m_Regularization->HessianMatVec(this->m_WorkVecField1, this->m_IncVelocityField); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \beta*\D{A}[\vect{\tilde{v}}] + \D{K}[\vect{\tilde{b}}]
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    ierr = this->m_WorkVecField1->AXPY(1.0, this->m_WorkVecField2); CHKERRQ(ierr);
    // pass to output
    if (Hvtilde != NULL) {
        ierr = this->m_WorkVecField1->GetComponents(Hvtilde); CHKERRQ(ierr);
    }


    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the analytically (spectrally) preconditioned
 * hessian, i.e.
 * P(H[\tilde{v}]) = (\beta A)^{-1}(\beta A[\tilde{v}] + b[\tilde{v}])
 *                 = \tilde{v} + (\beta A)^{-1}(b[\tilde{v}])
 * it is important to note, that this matrix is no longer symmetric;
 * we therefore can't use pcg
 *******************************************************************/
PetscErrorCode CLAIRE::PrecondHessMatVec(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType hd;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate container for incremental velocity field
    if (this->m_IncVelocityField == NULL) {
        try {this->m_IncVelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // parse input
    if (vtilde != NULL) {
        ierr = this->m_IncVelocityField->SetComponents(vtilde); CHKERRQ(ierr);
    }

    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t) and incremental body force
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // apply inverse of 2nd variation of regularization model to
    // incremental body force: (\beta \D{A})^{-1}\D{K}[\vect{\tilde{b}}]
    ierr = this->m_Regularization->ApplyInverse(this->m_WorkVecField1, this->m_WorkVecField2, false); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \vect{\tilde{v}} + (\beta \D{A})^{-1} \D{K}[\vect{\tilde{b}}]
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    hd = this->m_Opt->GetLebesgueMeasure();
    ierr = this->m_WorkVecField2->WAXPY(hd, this->m_IncVelocityField, this->m_WorkVecField1); CHKERRQ(ierr);

    // pass to output
    if (Hvtilde != NULL) {
        ierr = this->m_WorkVecField2->GetComponents(Hvtilde); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the analytically (spectrally) preconditioned
 * hessian, i.e.
 * P(H[\tilde{v}]) = \tilde{v} + P^{1/2}(b[\tilde{v}])P^{1/2}
 * it is important to note, that this matrix is no longer symmetric;
 * we therefore can't use pcg
 *******************************************************************/
PetscErrorCode CLAIRE::PrecondHessMatVecSym(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    ScalarType hd;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate container for incremental velocity field
    if (this->m_IncVelocityField == NULL) {
        try {this->m_IncVelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // allocate work vec field 5 (1,2,3, and 4 are overwritten
    // during the computation of the incremental forward and adjoint
    // solve and the computation of the incremental body force)
    if (this->m_WorkVecField5 == NULL) {
        try {this->m_WorkVecField5 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // parse input (store incremental velocity field \tilde{v})
    if (vtilde != NULL) {
        ierr = this->m_WorkVecField5->SetComponents(vtilde); CHKERRQ(ierr);
    } else {
        ierr = this->m_WorkVecField5->Copy(this->m_IncVelocityField); CHKERRQ(ierr);
    }
    // apply inverse of 2nd variation of regularization model to
    // incremental body force: (\beta\D{A})^{-1/2}\D{K}[\vect{\tilde{b}}](\beta\D{A})^{-1/2}

    // apply (\beta\D{A})^{-1/2} to incremental velocity field
    ierr = this->m_Regularization->ApplyInverse(this->m_IncVelocityField, this->m_WorkVecField5, true); CHKERRQ(ierr);

    // now solve the PDEs given the preconditioned incremental velocity field

    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t) and compute incremental body force
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // apply (\beta\D{A})^{-1/2} to incremental body force
    ierr = this->m_Regularization->ApplyInverse(this->m_WorkVecField1, this->m_WorkVecField2, true); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \vect{\tilde{v}} + (\beta \D{A})^{-1/2}\D{K}[\vect{\tilde{b}}](\beta \D{A})^{-1/2}
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    hd = this->m_Opt->GetLebesgueMeasure();
    ierr = this->m_WorkVecField5->Scale(hd); CHKERRQ(ierr);
    ierr = this->m_WorkVecField5->AXPY(1.0, this->m_WorkVecField1); CHKERRQ(ierr);

    // pass to output
    if (Hvtilde != NULL) {
        ierr = this->m_WorkVecField5->GetComponents(Hvtilde); CHKERRQ(ierr);
    }
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute initial condition given some initial guess for
 * the state variable $m$ and the adjoint variable $\lambda$
 * @param[in] m state variable
 * @param[in] lambda adjoint variable
 *******************************************************************/
PetscErrorCode CLAIRE::ComputeInitialCondition(Vec m, Vec lambda) {
    PetscErrorCode ierr = 0;
    IntType nt, nl, nc, ng;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    ext = this->m_Opt->m_FileNames.extension;

    // allocate container for incremental velocity field
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL) {
        ierr = this->SetupRegularization(); CHKERRQ(ierr);
    }

    // allocate state and adjoint variables
    if (this->m_StateVariable == NULL) {
        ierr = VecCreate(this->m_StateVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }

    // allocate state and adjoint variables
    if (this->m_AdjointVariable == NULL) {
        if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
            ierr = VecCreate(this->m_AdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
        } else {
            ierr = VecCreate(this->m_AdjointVariable, nc*nl, nc*ng); CHKERRQ(ierr);
        }
    }

    if (this->m_Opt->m_Verbosity > 2) {
        ierr = DbgMsg("piccard iteration"); CHKERRQ(ierr);
    }

    // copy input state and adjoint variable to class variables
    ierr = VecCopy(m, this->m_StateVariable); CHKERRQ(ierr);
    ierr = VecCopy(lambda, this->m_AdjointVariable); CHKERRQ(ierr);

    // compute body force (assigned to work vec field 2)
    // TODO: this will crash for GAUSS NEWTON
    ierr = this->ComputeBodyForce(); CHKERRQ(ierr);

    // piccard step: solve A[v] = - ht \sum_j \lambda^j grad(m^j)
    ierr = this->m_WorkVecField2->Scale(-1.0); CHKERRQ(ierr);

    // apply inverse regularization operator / spectral preconditioning
    ierr = this->m_Regularization->ApplyInverse(this->m_VelocityField, this->m_WorkVecField2); CHKERRQ(ierr);

    // reset the adjoint variables
    ierr = VecSet(this->m_StateVariable, 0.0); CHKERRQ(ierr);
    ierr = VecSet(this->m_AdjointVariable, 0.0); CHKERRQ(ierr);

    ierr = this->m_ReadWrite->Write(this->m_VelocityField, "initial-condition"+ext); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute the incremental body force
 * \tilde{\vect{b}} = \int_0^1 \igrad\tilde{m}\lambda
 *                           + \igrad m\tilde{\lambda} dt
 *******************************************************************/
PetscErrorCode CLAIRE::ComputeIncBodyForce() {
    PetscErrorCode ierr = 0;
    IntType nt, nl, nc, l;
    ScalarType *p_m = NULL, *p_mt = NULL, *p_l = NULL, *p_lt = NULL,
                *p_bt1 = NULL, *p_bt2 = NULL, *p_bt3 = NULL,
                *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL,
                *p_gradmt1 = NULL, *p_gradmt2 = NULL, *p_gradmt3 = NULL;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    ScalarType ht, scale, lj, ltj;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ht = this->m_Opt->GetTimeStepSize();
    ierr = Assert(nt > 0, "nt <= 0"); CHKERRQ(ierr);
    ierr = Assert(ht > 0, "ht <= 0"); CHKERRQ(ierr);
    ierr = Assert(nc > 0, "nc <= 0"); CHKERRQ(ierr);
    scale = ht;

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // init array
    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_bt1, p_bt2, p_bt3); CHKERRQ(ierr);

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);  // state variable for all t^j
    ierr = GetRawPointer(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);  // incremental adjoint variable for all t^j

    // TODO(andreas): add case for zero velocity field
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);

        if (this->m_WorkVecField3 == NULL) {
            try {this->m_WorkVecField3 = new VecField(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }

        ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);  // adjoint variable for all t^j
        ierr = GetRawPointer(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);  // incremental state variable for all t^j

        ierr = this->m_WorkVecField3->GetArrays(p_gradmt1, p_gradmt2, p_gradmt3); CHKERRQ(ierr);

        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j) {  // for all time points
            // trapezoidal rule (apply scaling)
            if ((j == 0) || (j == nt)) scale *= 0.5;
            for (IntType k = 0; k < nc; ++k) {  // for all image components
                l = j*nl*nc + k*nl;

                // computing gradient of m
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m+l, this->m_Opt->m_FFT.plan, &XYZ, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

                // computing gradient of \tilde{m}
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gradmt1, p_gradmt2, p_gradmt3, p_mt+l, this->m_Opt->m_FFT.plan, &XYZ, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
                // compute \vect{\tilde{b}}^k_i
                // += h_d*ht*(\tilde{\lambda}^j (\grad m^j)^k
                //    + \lambda^j (\grad \tilde{m}^j)^k)_i
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    // get \lambda(x_i,t^j) and \tilde{\lambda}(x_i,t^j)
                    lj  = p_l[l+i];
                    ltj = p_lt[l+i];

                    p_bt1[i] += scale*(p_gradm1[i]*ltj + p_gradmt1[i]*lj)/static_cast<ScalarType>(nc);
                    p_bt2[i] += scale*(p_gradm2[i]*ltj + p_gradmt2[i]*lj)/static_cast<ScalarType>(nc);
                    p_bt3[i] += scale*(p_gradm3[i]*ltj + p_gradmt3[i]*lj)/static_cast<ScalarType>(nc);
                }  // for all grid points
            }  // for all image components
}  // omp

            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)) scale *= 2.0;
        }  // for all time points

        ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);  // adjoint variable for all t^j
        ierr = RestoreRawPointer(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);  // incremental state variable for all t^j
        ierr = this->m_WorkVecField3->RestoreArrays(p_gradmt1, p_gradmt2, p_gradmt3); CHKERRQ(ierr);
    } else if (this->m_Opt->m_OptPara.method == GAUSSNEWTON) {  // gauss newton approximation
        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j) {  // for all time points
            // trapezoidal rule (apply scaling)
            if ((j == 0) || (j == nt)) scale *= 0.5;
            for (IntType k = 0; k < nc; ++k) {  // for all image components
                l = j*nl*nc + k*nl;

                // compute gradient of m^j
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m+l, this->m_Opt->m_FFT.plan, &XYZ, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
                // compute \vect{\tilde{b}}^k_i += h_d*ht*(\tilde{\lambda}^j (\grad m^j)^k
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    // get \tilde{\lambda}(x_i,t^j)
                    ltj = p_lt[l+i];
                    p_bt1[i] += scale*p_gradm1[i]*ltj/static_cast<ScalarType>(nc);
                    p_bt2[i] += scale*p_gradm2[i]*ltj/static_cast<ScalarType>(nc);
                    p_bt3[i] += scale*p_gradm3[i]*ltj/static_cast<ScalarType>(nc);
                }  // for all grid points
}  // omp

            }  // for all image components
            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)) scale *= 2.0;
        }   // for all time points
    } else {
        ierr = ThrowError("hessian approximation not implemented"); CHKERRQ(ierr);
    }

    // restore all arrays
    ierr = this->m_WorkVecField1->RestoreArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_bt1, p_bt2, p_bt3); CHKERRQ(ierr);

    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);  // state variable for all t^j
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);  // incremental adjoint variable for all t^j

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::StoreStateVariable() {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt;
    ScalarType *p_m = NULL, *p_mj = NULL;
    std::stringstream ss;
    std::string ext;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    ierr = Assert(nt > 0, "nt <= 0"); CHKERRQ(ierr);
    ext = this->m_Opt->m_FileNames.extension;

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

    // store time history
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    // store individual time points
    for (IntType j = 0; j <= nt; ++j) {
        for (IntType k = 0; k < nc; ++k) {
            ierr = GetRawPointer(this->m_WorkScaField1, &p_mj); CHKERRQ(ierr);
            try {std::copy(p_m+j*nl*nc + k*nl, p_m+j*nl*nc + (k+1)*nl, p_mj);}
            catch (std::exception& err) {
                ierr = ThrowError(err); CHKERRQ(ierr);
            }
            ierr = RestoreRawPointer(this->m_WorkScaField1, &p_mj); CHKERRQ(ierr);
            // write out
            ss.str(std::string()); ss.clear();

            if (nc > 1) {
                ss << "state-variable-k=" << std::setw(3) << std::setfill('0')  << k
                   << "-j=" << std::setw(3) << std::setfill('0') << j << ext;
            } else {
                ss << "state-variable-j=" << std::setw(3) << std::setfill('0') << j << ext;
            }
            ierr = this->m_ReadWrite->Write(this->m_WorkScaField1, ss.str()); CHKERRQ(ierr);
        }  // for number of vector components
    }  // for number of time points
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);


    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveStateEquation() {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt;
    ScalarType *p_m = NULL;
    std::stringstream ss;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);

    // check cfl condition / update time step
    if (this->m_Opt->m_PDESolver.monitorcflnumber ||
        this->m_Opt->m_PDESolver.adapttimestep) {
        ierr = this->ComputeCFLCondition(); CHKERRQ(ierr);
    }

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ierr = Assert(nt > 0, "nt <= 0"); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "solving state equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->m_Domain.nx[0]
           << "," << this->m_Opt->m_Domain.nx[1]
           << "," << this->m_Opt->m_Domain.nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }


    // set initial condition m_0 = m_T
    ierr = this->SetInitialState(this->m_TemplateImage); CHKERRQ(ierr);

    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        // we copy m_0 to all t for v=0
        if (this->m_Opt->m_RegFlags.runinversion) {
            ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
            for (IntType j = 1; j <= nt; ++j) {
                try {std::copy(p_m, p_m+nc*nl, p_m+j*nl*nc);}
                catch (std::exception& err) {
                    ierr = ThrowError(err); CHKERRQ(ierr);
                }
            }
            ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        }
    } else {
        // call the solver
        if (this->m_Opt->m_PDESolver.pdetype == CONTINUITYEQ) {
            ierr = this->SolveContinuityEquationSL(); CHKERRQ(ierr);
        } else if (this->m_Opt->m_PDESolver.pdetype == TRANSPORTEQ) {
            // call the solver
            switch (this->m_Opt->m_PDESolver.type) {
                case RK2:
                {
                    ierr = this->SolveStateEquationRK2(); CHKERRQ(ierr);
                    break;
                }
                case SL:
                {
                    ierr = this->SolveStateEquationSL(); CHKERRQ(ierr);
                    break;
                }
                default:
                {
                    ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                    break;
                }
            }
        } else {
            ierr = ThrowError("PDE type does not exist"); CHKERRQ(ierr);
        }

    }  // velocity field is zero

    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ScalarType maxval, minval, nvx1, nvx2, nvx3;
        ierr = VecMax(this->m_StateVariable, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(this->m_StateVariable, NULL, &minval); CHKERRQ(ierr);
        ss  << "state variable: [" << std::scientific
            << minval << "," << maxval << "]";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = this->m_VelocityField->Norm(nvx1, nvx2, nvx3); CHKERRQ(ierr);
        ss  << "velocity norm: (" << std::scientific
            << nvx1 << "," << nvx2 << "," << nvx3 <<")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }


    // store time series
    if (this->m_Opt->m_ReadWriteFlags.timeseries) {
        ierr = this->StoreStateVariable(); CHKERRQ(ierr);
    }

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveStateEquationRK2(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, l, lnext;
    ScalarType *p_m = NULL, *p_mbar = NULL, *p_rhs0 = NULL,
                *p_gmx1 = NULL, *p_gmx2 = NULL, *p_gmx3 = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL;
    ScalarType ht = 0.0, hthalf = 0.0, rhs1;
    bool store = true;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // flag to identify if we store the time history
    store = this->m_Opt->m_RegFlags.runinversion;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }

    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gmx1, p_gmx2, p_gmx3); CHKERRQ(ierr);

    // copy initial condition to buffer
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField1, &p_mbar); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField2, &p_rhs0); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {
        if (store) {
            l = j*nl*nc; lnext = (j+1)*nl*nc;
        } else {
            l = 0; lnext = 0;
        }

        for (IntType k = 0; k < nc; ++k) {
            // compute gradient of k-th component of m_j
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gmx1, p_gmx2, p_gmx3, p_m + l + k*nl, this->m_Opt->m_FFT.plan, &XYZ, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

            // evaluate right hand side and compute intermediate rk2 step
#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                 p_rhs0[i] = -p_gmx1[i]*p_vx1[i]
                             -p_gmx2[i]*p_vx2[i]
                             -p_gmx3[i]*p_vx3[i];

                 // compute intermediate result
                 p_mbar[i] = p_m[l + k*nl + i] + ht*p_rhs0[i];
            }
}  // omp
            // compute gradient of \bar{m}
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gmx1, p_gmx2, p_gmx3, p_mbar, this->m_Opt->m_FFT.plan, &XYZ, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

            // evaluate right hand side and wrap up integration
#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                rhs1 = -p_gmx1[i]*p_vx1[i]
                       -p_gmx2[i]*p_vx2[i]
                       -p_gmx3[i]*p_vx3[i];

                // we have overwritten m_j with intermediate result
                // m_{j+1} = m_j + 0.5*ht*(RHS0 + RHS1)
                p_m[lnext + k*nl + i] = p_m[l + k*nl + i] + hthalf*(p_rhs0[i] + rhs1);
            }
}  // omp
        }  // for all components
    }  // for all time points

    // copy initial condition to buffer
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_mbar); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField2, &p_rhs0); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(p_gmx1, p_gmx2, p_gmx3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveStateEquationSL(void) {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt, l, lnext;
    ScalarType *p_m = NULL;
    bool store = true;
    std::stringstream ss;
    std::string filename;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    // flag to identify if we store the time history
    store = this->m_Opt->m_RegFlags.runinversion;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);

    // get state variable m
    ierr = GetRawPointerReadWrite(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        if (store) {
            l = j*nl*nc; lnext = (j+1)*nl*nc;
        } else {
            l = 0; lnext = 0;
        }
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            // compute m(X,t^{j+1}) (interpolate state variable)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_m + lnext + k*nl, p_m + l + k*nl, "state"); CHKERRQ(ierr);
        }
    }
    
    ierr = RestoreRawPointerReadWrite(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveAdjointEquation() {
    PetscErrorCode ierr = 0;
    IntType nl, nc, ng, nt;
    ScalarType *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL,
               *p_b1 = NULL, *p_b2 = NULL, *p_b3 = NULL, *p_m = NULL, *p_l = NULL;
    ScalarType hd;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    double timer[NFFTTIMERS] = {0};
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    hd  = this->m_Opt->GetLebesgueMeasure();

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "solving adjoint equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->m_Domain.nx[0]
           << "," << this->m_Opt->m_Domain.nx[1]
           << "," << this->m_Opt->m_Domain.nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        if (this->m_AdjointVariable == NULL) {
            ierr = VecCreate(this->m_AdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
        }
    } else {
        if (this->m_AdjointVariable == NULL) {
            ierr = VecCreate(this->m_AdjointVariable, nc*nl, nc*ng); CHKERRQ(ierr);
        }
    }

    if (this->m_DistanceMeasure == NULL) {
        ierr = this->SetupDistanceMeasure(); CHKERRQ(ierr);
    }
    ierr = this->m_DistanceMeasure->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetAdjointVariable(this->m_AdjointVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetFinalConditionAE(); CHKERRQ(ierr);

    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
        // adjoint variable is constant in time
        if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
            for (IntType j = 1; j <= nt; ++j) {
                try {std::copy(p_l + nt*nc*nl, p_l + (nt+1)*nc*nl, p_l + (nt-j)*nl*nc);}
                catch (std::exception& err) {
                    ierr = ThrowError(err); CHKERRQ(ierr);
                }
            }
        }

        // allocate the memory for the computation of the body force
        if (this->m_WorkVecField1 == NULL) {
            try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }
        if (this->m_WorkVecField2 == NULL) {
            try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
            catch (std::bad_alloc& err) {
                ierr = reg::ThrowError(err); CHKERRQ(ierr);
            }
        }

        // init body force for numerical integration
        ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);

        // m and \lambda are constant in time
        ierr = GetRawPointer(this->m_TemplateImage, &p_m); CHKERRQ(ierr);
        ierr = this->m_WorkVecField1->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
        ierr = this->m_WorkVecField2->GetArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);

        for (IntType k = 0; k < nc; ++k) {  // for all components
            // compute gradient of m
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m + k*nl, this->m_Opt->m_FFT.plan, &xyz, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

            // b = \sum_k\int_{\Omega} \lambda_k \grad m_k dt
#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                ScalarType lambda = p_l[k*nl+i];
                p_b1[i] += lambda*p_gradm1[i]/static_cast<ScalarType>(nc);
                p_b2[i] += lambda*p_gradm2[i]/static_cast<ScalarType>(nc);
                p_b3[i] += lambda*p_gradm3[i]/static_cast<ScalarType>(nc);
            }
        }
}  // pragma
        ierr = RestoreRawPointer(this->m_TemplateImage, &p_m); CHKERRQ(ierr);
        ierr = this->m_WorkVecField1->RestoreArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
        ierr = this->m_WorkVecField2->RestoreArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);

        // for full newton method we have to store the adjoint variable
        ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);

        // increment timers
        this->m_Opt->IncreaseFFTTimers(timer);
    } else {
        // call the solver
        switch (this->m_Opt->m_PDESolver.type) {
            case RK2:
            {
                ierr = this->SolveAdjointEquationRK2(); CHKERRQ(ierr);
                break;
            }
            case SL:
            {
                ierr = this->SolveAdjointEquationSL(); CHKERRQ(ierr);
                break;
            }
            default:
            {
                ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                break;
            }
        }
    }

    // apply projection
    ierr = this->ApplyProjection(); CHKERRQ(ierr);

    // scale result by hd
    ierr = this->m_WorkVecField2->Scale(hd); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ScalarType maxval, minval;
        ierr = VecMin(this->m_AdjointVariable, NULL, &minval); CHKERRQ(ierr);
        ierr = VecMax(this->m_AdjointVariable, NULL, &maxval); CHKERRQ(ierr);
        ss << "adjoint variable: [" << std::scientific << minval << "," << maxval << "]";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0 (solved backward in time)
 *******************************************************************/
PetscErrorCode CLAIRE::SolveAdjointEquationRK2(void) {
    PetscErrorCode ierr;
    IntType nl, ng, nc, nt, ll, lm, llnext;
    ScalarType *p_l = NULL, *p_m = NULL, *p_rhs0 = NULL, *p_rhs1 = NULL,
               *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
               *p_vec1 = NULL, *p_vec2 = NULL, *p_vec3 = NULL,
               *p_b1 = NULL, *p_b2 = NULL, *p_b3 = NULL;
    ScalarType hthalf, ht, lambdabar, lambda, scale;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    double timer[NFFTTIMERS] = {0};
    bool fullnewton = false;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;
    scale = ht;

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);
    ierr = Assert(ht > 0, "ht < 0"); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }


    // for full newton we store $\lambda$
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        fullnewton = true;
    }
    ierr = GetRawPointer(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);

    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_vec1, p_vec2, p_vec3); CHKERRQ(ierr);

    // init body force for numerical integration
    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        lm = (nt-j)*nc*nl;
        if (fullnewton) {
            ll = (nt-j)*nc*nl; llnext = (nt-(j+1))*nc*nl;
        } else {
            ll = 0; llnext = 0;
        }

        // scaling for trapezoidal rule
        if (j == 0) scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            // scale \vect{v} by \lambda

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                lambda = p_l[ll + k*nl + i];
                p_vec1[i] = lambda*p_v1[i];
                p_vec2[i] = lambda*p_v2[i];
                p_vec3[i] = lambda*p_v3[i];
            }  // for all grid points
}  // omp
            // compute \idiv(\lambda\vect{v})
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_divergence_t(p_rhs0, p_vec1, p_vec2, p_vec3, this->m_Opt->m_FFT.plan, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTDIV);

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                // compute \bar{\lambda} = \lambda_j + ht*\idiv(\lambda\vect{v})
                lambdabar = p_l[ll + k*nl +i] + ht*p_rhs0[i];

                // scale \vect{v} by \bar{\lambda}
                p_vec1[i] = p_v1[i]*lambdabar;
                p_vec2[i] = p_v2[i]*lambdabar;
                p_vec3[i] = p_v3[i]*lambdabar;
            }
}  // omp
            // compute \idiv(\bar{\lambda}\vect{v})
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_divergence_t(p_rhs1, p_vec1, p_vec2, p_vec3, this->m_Opt->m_FFT.plan, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTDIV);

            // grad(m^j)
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_vec1, p_vec2, p_vec3, p_m + lm + k*nl, this->m_Opt->m_FFT.plan, &xyz, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                lambda = p_l[ll + k*nl + i];
                // second step of rk2 time integration
                p_l[llnext + k*nl + i] = lambda + hthalf*(p_rhs0[i] + p_rhs1[i]);

                // compute bodyforce
                p_b1[i] += scale*p_vec1[i]*lambda/static_cast<ScalarType>(nc);
                p_b2[i] += scale*p_vec2[i]*lambda/static_cast<ScalarType>(nc);
                p_b3[i] += scale*p_vec3[i]*lambda/static_cast<ScalarType>(nc);
            }
}  // omp
        }  // for all image components
        // trapezoidal rule (revert scaling)
        if (j == 0) scale *= 2.0;
    }  // for all time points

    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ll = k*nl; lm = k*nl;

        // compute gradient of m (for incremental body force)
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_vec1, p_vec2, p_vec3, p_m + lm, this->m_Opt->m_FFT.plan, &xyz, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            lambda = p_l[ll + i];
            // compute bodyforce
            p_b1[i] += 0.5*scale*p_vec1[i]*lambda/static_cast<ScalarType>(nc);
            p_b2[i] += 0.5*scale*p_vec2[i]*lambda/static_cast<ScalarType>(nc);
            p_b3[i] += 0.5*scale*p_vec3[i]*lambda/static_cast<ScalarType>(nc);
        }
    }

    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_vec1, p_vec2, p_vec3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);

    ierr = RestoreRawPointer(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveAdjointEquationSL() {
    PetscErrorCode ierr = 0;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_divv = NULL, *p_divvx = NULL,
                *p_l = NULL, *p_lx = NULL, *p_m = NULL,
                *p_vec1 = NULL, *p_vec2 = NULL, *p_vec3 = NULL,
                *p_b1 = NULL, *p_b2 = NULL, *p_b3 = NULL;
    ScalarType ht, lambdax, lambda, rhs0, rhs1, scale;
    IntType nl, ng, nc, nt, ll, lm, llnext;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    bool fullnewton = false;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    scale = ht;

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL) {
        ierr = VecCreate(this->m_WorkScaField3, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // compute trajectory for adjoint equations
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    // for full newton we store the adjoint variable
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        fullnewton = true;
    }
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField3, &p_lx); CHKERRQ(ierr);

    // compute divergence of velocity field
    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(p_divv, p_v1, p_v2, p_v3, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    // interpolate velocity field v(X)
//    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField1, this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    // compute divergence of velocity field at X
    ierr = this->m_WorkVecField1->GetArrays(p_vec1, p_vec2, p_vec3); CHKERRQ(ierr);
//    this->m_Opt->StartTimer(FFTSELFEXEC);
//    accfft_divergence_t(p_divvx, p_vec1, p_vec2, p_vec3, this->m_Opt->m_FFT.plan, timer);
//    this->m_Opt->StopTimer(FFTSELFEXEC);
//    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    // evaluate div(v) along characteristic X
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_divvx, p_divv, "adjoint"); CHKERRQ(ierr);

    // init body force for numerical integration
    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);

    // perform numerical time integration for adjoint variable and
    // add up body force
    for (IntType j = 0; j < nt; ++j) {
        lm = (nt-j)*nc*nl;
        if (fullnewton) {
            ll = (nt-j)*nc*nl; llnext = (nt-(j+1))*nc*nl;
        } else {
            ll = 0; llnext = 0;
        }

        // scaling for trapezoidal rule (for body force)
        if (j == 0) scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {
            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_lx, p_l + ll + k*nl, "adjoint"); CHKERRQ(ierr);

            // compute gradient of m (for incremental body force)
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_vec1, p_vec2, p_vec3, p_m + lm + k*nl, this->m_Opt->m_FFT.plan, &xyz, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);
#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                lambda  = p_l[ll + k*nl + i];
                lambdax = p_lx[i];

                rhs0 = lambdax*p_divvx[i];
                rhs1 = (lambdax + ht*rhs0)*p_divv[i];

                // compute \lambda(x,t^{j+1})
                p_l[llnext + k*nl + i] = lambdax + 0.5*ht*(rhs0 + rhs1);

                // compute bodyforce
                p_b1[i] += scale*p_vec1[i]*lambda/static_cast<ScalarType>(nc);
                p_b2[i] += scale*p_vec2[i]*lambda/static_cast<ScalarType>(nc);
                p_b3[i] += scale*p_vec3[i]*lambda/static_cast<ScalarType>(nc);
            }
        }
}  // omp
        // trapezoidal rule (revert scaling; for body force)
        if (j == 0) scale *= 2.0;
    }

    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ll = k*nl; lm = k*nl;

        // compute gradient of m (for incremental body force)
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_vec1, p_vec2, p_vec3, p_m + lm, this->m_Opt->m_FFT.plan, &xyz, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            lambda = p_l[ll + i];
            // compute bodyforce
            p_b1[i] += 0.5*scale*p_vec1[i]*lambda/static_cast<ScalarType>(nc);
            p_b2[i] += 0.5*scale*p_vec2[i]*lambda/static_cast<ScalarType>(nc);
            p_b3[i] += 0.5*scale*p_vec3[i]*lambda/static_cast<ScalarType>(nc);
        }
    }
}  // omp
    ierr = this->m_WorkVecField2->RestoreArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_vec1, p_vec2, p_vec3); CHKERRQ(ierr);

    ierr = RestoreRawPointer(this->m_WorkScaField3, &p_lx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \idiv m\vect{v} = 0  with initial condition m_0 = m_T
 * (solved forward in time)
 *******************************************************************/
PetscErrorCode CLAIRE::SolveContinuityEquationSL() {
    PetscErrorCode ierr = 0;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_divv = NULL, *p_divvx = NULL,
                *p_m = NULL, *p_mx = NULL;
    ScalarType mx, rhs0, rhs1, ht, hthalf;
    IntType nl, ng, nc, nt, l, lnext;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    bool store;

    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // flag to identify if we store the time history
    store = this->m_Opt->m_RegFlags.runinversion;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL) {
        ierr = VecCreate(this->m_WorkScaField3, nl, ng); CHKERRQ(ierr);
    }


    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);


    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField3, &p_mx); CHKERRQ(ierr);

    // compute divergence of velocity field
    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(p_divv, p_v1, p_v2, p_v3, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    // interpolate velocity field v(X)
//    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField1, this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    // compute divergence of velocity field at X
//    ierr = this->m_WorkVecField1->GetArrays(p_vec1, p_vec2, p_vec3); CHKERRQ(ierr);
//    this->m_Opt->StartTimer(FFTSELFEXEC);
//    accfft_divergence_t(p_divvx, p_vec1, p_vec2, p_vec3, this->m_Opt->m_FFT.plan, timer);
//    this->m_Opt->StopTimer(FFTSELFEXEC);
//    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    // evaluate div(v) along characteristic X
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_divvx, p_divv, "state"); CHKERRQ(ierr);

    // perform numerical time integration for state variable
    for (IntType j = 0; j < nt; ++j) {
        if (store) {
            l = j*nl*nc; lnext = (j+1)*nl*nc;
        } else {
            l = 0; lnext = 0;
        }

        // scaling for trapezoidal rule (for body force)
        for (IntType k = 0; k < nc; ++k) {
            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_mx, p_m + l + k*nl, "state"); CHKERRQ(ierr);

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                mx = p_mx[i];

                rhs0 = -mx*p_divvx[i];
                rhs1 = -(mx + ht*rhs0)*p_divv[i];
                //if (std::abs(p_divv[i]) > 0.1) { std::cout << p_divv[i] << " ";}
                // compute \lambda(x,t^{j+1})
                p_m[lnext + k*nl + i] = mx + hthalf*(rhs0 + rhs1);
            }
        }
}  // omp
    }

    ierr = RestoreRawPointer(this->m_WorkScaField3, &p_mx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the incremental state equation
 * \p_t \tilde{m} + \igrad m \cdot \vect{\tilde{v}}
 *                + \igrad \tilde{m} \cdot \vect{v} = 0
 * subject to \tilde{m}_0 = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncStateEquation(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "solving incremental state equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->m_Domain.nx[0]
           << "," << this->m_Opt->m_Domain.nx[1]
           << "," << this->m_Opt->m_Domain.nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // allocate variables
    if (this->m_IncStateVariable == NULL) {
        if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
            ierr = VecCreate(this->m_IncStateVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
        } else {
            ierr = VecCreate(this->m_IncStateVariable, nc*nl, nc*ng); CHKERRQ(ierr);
        }
    }

    // set initial value
    ierr = VecSet(this->m_IncStateVariable, 0.0); CHKERRQ(ierr);


    // start timer
    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);


    // call the solver
    switch (this->m_Opt->m_PDESolver.type) {
        case RK2:
        {
            ierr = this->SolveIncStateEquationRK2(); CHKERRQ(ierr);
            break;
        }
        case SL:
        {
            ierr = this->SolveIncStateEquationSL(); CHKERRQ(ierr);
            //ierr = this->SolveIncStateEquationRK2(); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
            break;
        }
    }

    if (this->m_Opt->m_Verbosity > 2) {
        ScalarType maxval, minval;
        ierr = VecMax(this->m_IncStateVariable, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(this->m_IncStateVariable, NULL, &minval); CHKERRQ(ierr);
        ss << "incremental state variable: [" << std::scientific << minval << "," << maxval << "]";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // stop timer
    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief solve the incremental state equation
 * \p_t \tilde{m} + \igrad m \cdot \vect{\tilde{v}}
 *                + \igrad \tilde{m} \cdot \vect{v} = 0
 * subject to \tilde{m}_0 = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncStateEquationRK2(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, lm, lmnext, lmt, lmtnext;
    ScalarType *p_m = NULL, *p_mt = NULL, *p_mtbar = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_gmx1 = NULL, *p_gmx2 = NULL, *p_gmx3 = NULL,
                *p_gmtx1 = NULL, *p_gmtx2 = NULL, *p_gmtx3 = NULL,
                *p_vtx1 = NULL, *p_vtx2 = NULL, *p_vtx3 = NULL, *p_rhs0 = NULL;
    ScalarType ht, hthalf;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    double timer[NFFTTIMERS] = {0};
    bool fullnewton = false;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);
    ierr = Assert(ht > 0, "time step size < 0"); CHKERRQ(ierr);

    // allocate variables
    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {   // gauss newton
        fullnewton = true;
    }

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->GetArrays(p_gmx1, p_gmx2, p_gmx3); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArrays(p_vtx1, p_vtx2, p_vtx3); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        // compute gradient of first time point of image component
        for (IntType k = 0; k < nc; ++k) {
            // template image is constant in time
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gmx1, p_gmx2, p_gmx3, p_m + k*nl, this->m_Opt->m_FFT.plan, &XYZ, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

            // compute incremental state variable for all time points
            // note: we do not need to store the time history for
            // \tilde{m} if we consider a Gauss--Newton approximation
            for (IntType j = 0; j < nt; ++j) {
                if (fullnewton) {
                    lmt = j*nl*nc; lmtnext = (j+1)*nl*nc;
                } else {
                    lmt = 0; lmtnext = 0;
                }
#pragma omp parallel
{
#pragma omp for
                // the right hand side remains constant;
                // we can reduce the 2 RK2 steps to a single one
                for (IntType i = 0; i < nl; ++i) {
                     p_mt[lmtnext + k*nl + i] = p_mt[lmt + k*nl + i] - ht*(p_gmx1[i]*p_vtx1[i]
                                                                         + p_gmx2[i]*p_vtx2[i]
                                                                         + p_gmx3[i]*p_vtx3[i]);
                }
}  // omp
            }  // for all time points
        }  // for all image components
    } else {  // velocity field is non-zero
        ierr = GetRawPointer(this->m_WorkScaField1, &p_mtbar); CHKERRQ(ierr);
        ierr = GetRawPointer(this->m_WorkScaField2, &p_rhs0); CHKERRQ(ierr);

        ierr = this->m_WorkVecField2->GetArrays(p_gmtx1, p_gmtx2, p_gmtx3); CHKERRQ(ierr);
        ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j) {
            for (IntType k = 0; k < nc; ++k) {
                lm = j*nl*nc; lmnext = (j+1)*nl*nc;
                if (fullnewton) {
                    lmt = j*nl*nc; lmtnext = (j+1)*nl*nc;
                } else {
                    lmt = 0; lmtnext = 0;
                }

                // compute gradient of m_j
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gmx1, p_gmx2, p_gmx3, p_m + lm + k*nl, this->m_Opt->m_FFT.plan, &XYZ, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

                // compute gradient of \tilde{m}_j
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gmtx1, p_gmtx2, p_gmtx3, p_mt + lmt + k*nl, this->m_Opt->m_FFT.plan, &XYZ, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

                for (IntType i = 0; i < nl; ++i) {
                     p_rhs0[i] = -p_gmtx1[i]*p_vx1[i] - p_gmtx2[i]*p_vx2[i] - p_gmtx3[i]*p_vx3[i]
                                 -p_gmx1[i]*p_vtx1[i] - p_gmx2[i]*p_vtx2[i] - p_gmx3[i]*p_vtx3[i];
                     // compute intermediate result
                     p_mtbar[i] = p_mt[lmt + k*nl + i] + ht*p_rhs0[i];
                }

                // compute gradient of m_{j+1}
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gmx1, p_gmx2, p_gmx3, p_m + lmnext, this->m_Opt->m_FFT.plan, &XYZ, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

                // compute gradient of \tilde{m}_j
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gmtx1, p_gmtx2, p_gmtx3, p_mtbar, this->m_Opt->m_FFT.plan, &XYZ, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
                for (IntType i = 0; i < nl; ++i) {
                    // evaluate right hand side
                    ScalarType rhs1 = -p_gmtx1[i]*p_vx1[i] - p_gmtx2[i]*p_vx2[i] - p_gmtx3[i]*p_vx3[i]
                                      -p_gmx1[i]*p_vtx1[i] - p_gmx2[i]*p_vtx2[i] - p_gmx3[i]*p_vtx3[i];

                    // compute intermediate result
                    p_mt[lmtnext + k*nl + i] = p_mt[lmt + k*nl + i] + hthalf*(rhs1 + p_rhs0[i]);
                }
            }  // for all image components
}  // omp
        }  // for all time points

        // copy initial condition to buffer
        ierr = RestoreRawPointer(this->m_WorkScaField1, &p_mtbar); CHKERRQ(ierr);
        ierr = RestoreRawPointer(this->m_WorkScaField2, &p_rhs0); CHKERRQ(ierr);

        ierr = this->m_WorkVecField2->RestoreArrays(p_gmtx1, p_gmtx2, p_gmtx3); CHKERRQ(ierr);
        ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    }  // velzero

    ierr = this->m_IncVelocityField->RestoreArrays(p_vtx1, p_vtx2, p_vtx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_gmx1, p_gmx2, p_gmx3); CHKERRQ(ierr);

    ierr = RestoreRawPointer(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief solve the incremental state equation
 * \p_t \tilde{m} + \igrad m \cdot \vect{\tilde{v}}
 *                + \igrad \tilde{m} \cdot \vect{v} = 0
 * subject to \tilde{m}_0 = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncStateEquationSL(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nt, nc, lm, lmnext, lmt, lmtnext;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    ScalarType ht, hthalf;
    ScalarType *p_gm1 = NULL, *p_gm2 = NULL, *p_gm3 = NULL,
                *p_mtilde = NULL, *p_m = NULL, *p_mx = NULL;
    const ScalarType *p_vtilde1 = NULL, *p_vtilde2 = NULL, *p_vtilde3 = NULL,
                     *p_vtildex1 = NULL, *p_vtildex2 = NULL, *p_vtildex3 = NULL;
    double timer[NFFTTIMERS] = {0};
    bool fullnewton = false;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    }

    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {   // gauss newton
        fullnewton = true;
    }

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField1, &p_mx); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gm1, p_gm2, p_gm3); CHKERRQ(ierr);

    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField2, this->m_IncVelocityField, "state"); CHKERRQ(ierr);

    ierr = this->m_WorkVecField2->GetArraysRead(p_vtildex1, p_vtildex2, p_vtildex3); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArraysRead(p_vtilde1, p_vtilde2, p_vtilde3); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j) {  // for all time points
        lm = j*nl*nc; lmnext = (j+1)*nl*nc;
        if (fullnewton) {   // full newton
            lmt = j*nl*nc; lmtnext = (j+1)*nl*nc;
        } else {
            lmt = 0; lmtnext = 0;
        }

        for (IntType k = 0; k < nc; ++k) {  // for all image components
            // interpolate incremental adjoint variable \tilde{m}^j(X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_mtilde + lmtnext + k*nl, p_mtilde + lmt + k*nl, "state"); CHKERRQ(ierr);
            // interpolate m
//            ierr = this->m_SemiLagrangianMethod->Interpolate(p_mx, p_m + lm + k*nl, "state"); CHKERRQ(ierr);

            // compute gradient for state variable
//            this->m_Opt->StartTimer(FFTSELFEXEC);
//            accfft_grad_t(p_gm1, p_gm2, p_gm3, p_mx, this->m_Opt->m_FFT.plan, &XYZ, timer);
//            this->m_Opt->StopTimer(FFTSELFEXEC);
//            this->m_Opt->IncrementCounter(FFT, FFTGRAD);


            // compute gradient for state variable
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gm1, p_gm2, p_gm3, p_m + lm + k*nl, this->m_Opt->m_FFT.plan, &XYZ, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

            ierr = this->m_SemiLagrangianMethod->Interpolate(p_gm1, p_gm2, p_gm3, p_gm1, p_gm2, p_gm3, "state"); CHKERRQ(ierr);

            // first part of time integration
#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                p_mtilde[lmtnext + k*nl + i] -= hthalf*(p_gm1[i]*p_vtildex1[i]
                                                      + p_gm2[i]*p_vtildex2[i]
                                                      + p_gm3[i]*p_vtildex3[i]);
            }
}  // omp
            // compute gradient for state variable at next time time point
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gm1, p_gm2, p_gm3, p_m + lmnext + k*nl, this->m_Opt->m_FFT.plan, &XYZ, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
            // second part of time integration
            for (IntType i = 0; i < nl; ++i) {
                p_mtilde[lmtnext + k*nl + i] -= hthalf*(p_gm1[i]*p_vtilde1[i]
                                                      + p_gm2[i]*p_vtilde2[i]
                                                      + p_gm3[i]*p_vtilde3[i]);
            }
        }  // for all image components
}  // omp
    }  // for all time points

    ierr = this->m_IncVelocityField->RestoreArraysRead(p_vtilde1, p_vtilde2, p_vtilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArraysRead(p_vtildex1, p_vtildex2, p_vtildex3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_gm1, p_gm2, p_gm3); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_mx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncAdjointEquation(void) {
    PetscErrorCode ierr = 0;
    ScalarType *p_ltilde = NULL, *p_m = NULL,
               *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL,
               *p_btilde1 = NULL, *p_btilde2 = NULL, *p_btilde3 = NULL;
    IntType nl, ng, nc, nt;
    ScalarType hd;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    double timer[NFFTTIMERS] = {0};
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    hd  = this->m_Opt->GetLebesgueMeasure();
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);

    if (this->m_Opt->m_Verbosity > 2) {
        ss << "solving incremental adjoint equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->m_Domain.nx[0]
           << "," << this->m_Opt->m_Domain.nx[1]
           << "," << this->m_Opt->m_Domain.nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // allocate state and adjoint variables
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        if (this->m_IncAdjointVariable == NULL) {
            ierr = VecCreate(this->m_IncAdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
        }
    } else {
        if (this->m_IncAdjointVariable == NULL) {
            ierr = VecCreate(this->m_IncAdjointVariable, nc*nl, nc*ng); CHKERRQ(ierr);
        }
    }

    if (this->m_DistanceMeasure == NULL) {
        ierr = this->SetupDistanceMeasure(); CHKERRQ(ierr);
    }
    ierr = this->m_DistanceMeasure->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetStateVariable(this->m_StateVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetIncStateVariable(this->m_IncStateVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetIncAdjointVariable(this->m_IncAdjointVariable); CHKERRQ(ierr);
    ierr = this->m_DistanceMeasure->SetFinalConditionIAE(); CHKERRQ(ierr);

    // check if velocity field is zero
    if (this->m_Opt->m_OptPara.method == GAUSSNEWTON) {   // gauss newton
        ierr = this->IsVelocityZero(); CHKERRQ(ierr);
        if (this->m_VelocityIsZero) {
            // since we're already in gauss newton mode, we do not
            // need to copy anything
            // copy terminal condition \tilde{\lambda}_1 = -\tilde{m}_1 to all time points
            ierr = GetRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);

            if (this->m_WorkVecField1 == NULL) {
                try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
                catch (std::bad_alloc& err) {
                    ierr = reg::ThrowError(err); CHKERRQ(ierr);
                }
            }
            if (this->m_WorkVecField2 == NULL) {
                try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
                catch (std::bad_alloc& err) {
                    ierr = reg::ThrowError(err); CHKERRQ(ierr);
                }
            }

            // m and \lambda are constant in time
            ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
            ierr = this->m_WorkVecField1->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);

            // init body force for numerical integration
            ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);
            ierr = this->m_WorkVecField2->GetArrays(p_btilde1, p_btilde2, p_btilde3); CHKERRQ(ierr);

            // $m$ and $\tilde{\lambda}$ are constant
            for (IntType k = 0; k < nc; ++k) {  // for all components
                // compute gradient of m
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m + k*nl, this->m_Opt->m_FFT.plan, &xyz, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
                // b = \sum_k\int_{\Omega} \lambda_k \grad m_k dt
                for (IntType i = 0; i < nl; ++i) {
                    ScalarType ltilde = p_ltilde[k*nl + i];
                    p_btilde1[i] += ltilde*p_gradm1[i]/static_cast<ScalarType>(nc);
                    p_btilde2[i] += ltilde*p_gradm2[i]/static_cast<ScalarType>(nc);
                    p_btilde3[i] += ltilde*p_gradm3[i]/static_cast<ScalarType>(nc);
                }
            }
}  // omp
            ierr = this->m_WorkVecField2->RestoreArrays(p_btilde1, p_btilde2, p_btilde3); CHKERRQ(ierr);
            ierr = this->m_WorkVecField1->RestoreArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
            ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

            ierr = RestoreRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);

            // increment timers
            this->m_Opt->IncreaseFFTTimers(timer);
        } else {
            // call the solver
            switch (this->m_Opt->m_PDESolver.type) {
                case RK2:
                {
                    ierr = this->SolveIncAdjointEquationGNRK2(); CHKERRQ(ierr);
                    break;
                }
                case SL:
                {
                    ierr = this->SolveIncAdjointEquationGNSL(); CHKERRQ(ierr);
                    break;
                }
                default:
                {
                    ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                    break;
                }
            }
        }  // zero velocity field
    } else if (this->m_Opt->m_OptPara.method == FULLNEWTON) {   // full newton
        // call the solver
        switch (this->m_Opt->m_PDESolver.type) {
            case RK2:
            {
                ierr = ThrowError("not tested"); CHKERRQ(ierr);
                ierr = this->SolveIncAdjointEquationFNRK2(); CHKERRQ(ierr);
                break;
            }
            case SL:
            {
                ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
//                ierr = this->SolveIncAdjointEquationRNRK2(); CHKERRQ(ierr);
                break;
            }
            default:
            {
                ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
                break;
            }
        }
    } else {
        ierr = ThrowError("update method not defined"); CHKERRQ(ierr);
    }

    // apply K[\tilde{b}]
    ierr = this->ApplyProjection(); CHKERRQ(ierr);

    // scale result by hd
    ierr = this->m_WorkVecField2->Scale(hd); CHKERRQ(ierr);


    if (this->m_Opt->m_Verbosity > 2) {
        ScalarType maxval, minval;
        ierr = VecMax(this->m_IncAdjointVariable, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(this->m_IncAdjointVariable, NULL, &minval); CHKERRQ(ierr);
        ss << "incremental adjoint variable: [" << std::scientific << minval << "," << maxval << "]";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PDESOLVE);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncAdjointEquationGNRK2(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, ll, lm;
    ScalarType *p_ltilde = NULL, *p_rhs0 = NULL, *p_rhs1 = NULL, *p_m = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_bt1 = NULL, *p_bt2 = NULL, *p_bt3 = NULL,
                *p_ltjvx1 = NULL, *p_ltjvx2 = NULL, *p_ltjvx3 = NULL,
                *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    ScalarType ht, hthalf, scale, ltilde;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    scale = ht;
    hthalf = 0.5*ht;

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3 == NULL) {
        try {this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }


    ierr = GetRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_ltjvx1, p_ltjvx2, p_ltjvx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    // init body force for numerical integration
    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_bt1, p_bt2, p_bt3); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        lm = (nt-j)*nc*nl;
        if (j == 0) scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            ll = k*nl;
#pragma omp parallel
{
#pragma omp for
            // scale \vect{v} by \lambda
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                ScalarType lt = p_ltilde[ll + i];

                p_ltjvx1[i] = p_vx1[i]*lt;
                p_ltjvx2[i] = p_vx2[i]*lt;
                p_ltjvx3[i] = p_vx3[i]*lt;
            }  // for all grid points
}  // omp
            // compute \idiv(\tilde{\lambda}\vect{v})
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_divergence_t(p_rhs0, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->m_FFT.plan, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTDIV);

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                // compute \bar{\tilde{\lambda}} = \tilde{\lambda}^j + ht*\idiv(\tilde{\lambda}^j\vect{v})
                ScalarType ltbar = p_ltilde[ll + i] + ht*p_rhs0[i];

                // scale \vect{v} by \bar{\lambda}
                p_ltjvx1[i] = p_vx1[i]*ltbar;
                p_ltjvx2[i] = p_vx2[i]*ltbar;
                p_ltjvx3[i] = p_vx3[i]*ltbar;
            }
}  // omp
            // compute \idiv(\bar{\lambda}\vect{v})
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_divergence_t(p_rhs1, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->m_FFT.plan, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTDIV);

            // compute gradient of m^j
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m + lm + k*nl, this->m_Opt->m_FFT.plan, &xyz, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                ltilde = p_ltilde[ll + i];    // get \tilde{\lambda}(x)

                // compute integration
                p_ltilde[ll + i] = ltilde + hthalf*(p_rhs0[i]+p_rhs1[i]);

                // compute incremental body force
                p_bt1[i] += scale*p_gradm1[i]*ltilde/static_cast<ScalarType>(nc);
                p_bt2[i] += scale*p_gradm2[i]*ltilde/static_cast<ScalarType>(nc);
                p_bt3[i] += scale*p_gradm3[i]*ltilde/static_cast<ScalarType>(nc);
            }
}  // omp
        }  // for all image components
        if (j == 0) scale *= 2.0;
    }  // for all time points

    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ll = k*nl; lm = k*nl;

        // compute gradient of m (for incremental body force)
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m + lm, this->m_Opt->m_FFT.plan, &xyz, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            ltilde = p_ltilde[ll + i];
            // compute bodyforce
            p_bt1[i] += 0.5*scale*p_gradm1[i]*ltilde/static_cast<ScalarType>(nc);
            p_bt2[i] += 0.5*scale*p_gradm2[i]*ltilde/static_cast<ScalarType>(nc);
            p_bt3[i] += 0.5*scale*p_gradm3[i]*ltilde/static_cast<ScalarType>(nc);
        }
}  // omp
    }

    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_ltjvx1, p_ltjvx2, p_ltjvx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_bt1, p_bt2, p_bt3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncAdjointEquationFNRK2(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, l, lnext;
    ScalarType *p_l = NULL, *p_lt = NULL, *p_rhs0 = NULL, *p_rhs1 = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_vtx1 = NULL, *p_vtx2 = NULL, *p_vtx3 = NULL,
                *p_ltjvx1 = NULL, *p_ltjvx2 = NULL, *p_ltjvx3 = NULL;
    ScalarType ht, hthalf, lambda, lambdatilde, ltbar;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);

    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_ltjvx1, p_ltjvx2, p_ltjvx3); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArrays(p_vtx1, p_vtx2, p_vtx3); CHKERRQ(ierr);

    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            l = k*nl;
            // lambda and v are constant in time

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                ScalarType lambda = p_l[l+i];

                // scale \vect{v} by \lambda
                p_ltjvx1[i] = p_vtx1[i]*lambda;
                p_ltjvx2[i] = p_vtx2[i]*lambda;
                p_ltjvx3[i] = p_vtx3[i]*lambda;
            }  // for all grid points
}  // omp
            // compute \idiv(\tilde{\lambda}\vect{v})
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_divergence_t(p_rhs0, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->m_FFT.plan, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTDIV);

            // compute numerical time integration
            for (IntType j = 0; j < nt; ++j) {  // for all time points
                l = (nt-j)*nc*nl + k*nl;
                lnext = (nt-(j+1))*nc*nl + k*nl;

#pragma omp parallel
{
#pragma omp for
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    p_lt[lnext+i] = p_lt[l+i] + ht*p_rhs0[i];
                }
            }  // for all time points
}  // omp
        }  // for all image components
    } else {  // velocity is zero
        if (this->m_WorkScaField2 == NULL) {
            ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
        }

        ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
        ierr = GetRawPointer(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j) {  // for all time points
            for (IntType k = 0; k < nc; ++k) {  // for all image components
                l = (nt-j)*nc*nl + k*nl;
                lnext = (nt-(j+1))*nc*nl + k*nl;

#pragma omp parallel
{
#pragma omp for
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    lambda  = p_l[l+i];
                    lambdatilde = p_lt[l+i];

                    p_ltjvx1[i] = p_vx1[i]*lambdatilde + p_vtx1[i]*lambda;
                    p_ltjvx2[i] = p_vx2[i]*lambdatilde + p_vtx2[i]*lambda;
                    p_ltjvx3[i] = p_vx3[i]*lambdatilde + p_vtx3[i]*lambda;
                }  // for all grid points
}  // omp
                // compute \idiv(\tilde{\lambda}\vect{v})
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_divergence_t(p_rhs0, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->m_FFT.plan, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTDIV);

#pragma omp parallel
{
#pragma omp for
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    // \bar{\lambda} = \tilde{\lambda}^j + ht*\idiv(\lambda^j\vect{v})
                    ltbar = p_lt[l+i] + ht*p_rhs0[i];
                    lambda = p_l[lnext+i];

                    // v \bar{\lambda} + \vect{\tilde{v}}\lambda^{j+1}
                    p_ltjvx1[i] = p_vx1[i]*ltbar + p_vtx1[i]*lambda;
                    p_ltjvx2[i] = p_vx2[i]*ltbar + p_vtx2[i]*lambda;
                    p_ltjvx3[i] = p_vx3[i]*ltbar + p_vtx3[i]*lambda;
                }
}  // omp
                // compute \idiv(\bar{\lambda}\vect{v})
                this->m_Opt->StartTimer(FFTSELFEXEC);
                accfft_divergence_t(p_rhs1, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->m_FFT.plan, timer);
                this->m_Opt->StopTimer(FFTSELFEXEC);
                this->m_Opt->IncrementCounter(FFT, FFTDIV);

#pragma omp parallel
{
#pragma omp for
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    p_lt[lnext+i] = p_lt[l+i] + hthalf*(p_rhs0[i]+p_rhs1[i]);
                }
}  // omp
            }  // for all image components
        }  // for all time points
        ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
        ierr = RestoreRawPointer(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);
    }  // velzero

    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);

    ierr = this->m_IncVelocityField->RestoreArrays(p_vtx1, p_vtx2, p_vtx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_ltjvx1, p_ltjvx2, p_ltjvx3); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timer);
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncAdjointEquationGNSL(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, ll, lm;
    ScalarType *p_ltilde = NULL, *p_ltildex = NULL, *p_m = NULL,
                *p_divv = NULL, *p_divvx = NULL,
                *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL,
                *p_bt1 = NULL, *p_bt2 = NULL, *p_bt3 = NULL,
                *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL;
    ScalarType ht, hthalf, ltilde, ltildex, rhs0, rhs1, scale;
    std::bitset<3> xyz; xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    double timer[NFFTTIMERS] = {0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    ht = this->m_Opt->GetTimeStepSize();
    scale = ht;
    hthalf = 0.5*ht;

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField3 == NULL) {
        ierr = VecCreate(this->m_WorkScaField3, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField1); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }

    // compute divergence of velocity field
    ierr = GetRawPointer(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = this->m_VelocityField->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(p_divv, p_v1, p_v2, p_v3, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    ierr = this->m_VelocityField->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    // compute v(X)
//    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField1, this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    // evaluate divergence at X
//    ierr = GetRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
//    ierr = this->m_WorkVecField1->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
//    this->m_Opt->StartTimer(FFTSELFEXEC);
//    accfft_divergence_t(p_divvx, p_v1, p_v2, p_v3, this->m_Opt->m_FFT.plan, timer);
//    this->m_Opt->StopTimer(FFTSELFEXEC);
//    ierr = this->m_WorkVecField1->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);


    ierr = GetRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_divvx, p_divv, "adjoint"); CHKERRQ(ierr);

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField3, &p_ltildex); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);

    // initialize work vec field
    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_bt1, p_bt2, p_bt3); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j) {
        lm = (nt-j)*nc*nl;
        if (j == 0) scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {
            ll = k*nl;
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_ltildex, p_ltilde + ll, "adjoint"); CHKERRQ(ierr);

            // compute gradient of m^j
            this->m_Opt->StartTimer(FFTSELFEXEC);
            accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m + lm + k*nl, this->m_Opt->m_FFT.plan, &xyz, timer);
            this->m_Opt->StopTimer(FFTSELFEXEC);
            this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
            for (IntType i = 0; i < nl; ++i) {
                ltilde  = p_ltilde[ll + i];   // get \tilde{\lambda}(x)
                ltildex = p_ltildex[i];       // get \tilde{\lambda}(X) (interpolated)

                // scale div(v)(X) by \tilde{\lambda}(X)
                rhs0 = ltildex*p_divvx[i];

                // scale div(v) by \tilde{\lambda}*
                rhs1 = (ltildex + ht*rhs0)*p_divv[i];

                // final rk2 step
                p_ltilde[ll + i] = ltildex + hthalf*(rhs0 + rhs1);

                p_bt1[i] += scale*p_gradm1[i]*ltilde/static_cast<ScalarType>(nc);
                p_bt2[i] += scale*p_gradm2[i]*ltilde/static_cast<ScalarType>(nc);
                p_bt3[i] += scale*p_gradm3[i]*ltilde/static_cast<ScalarType>(nc);
            }
}  // omp
        }  // for all image components
        if (j == 0) scale *= 2.0;
    }  // for all time points

    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        ll = k*nl; lm = k*nl;

        // compute gradient of m (for incremental body force)
        this->m_Opt->StartTimer(FFTSELFEXEC);
        accfft_grad_t(p_gradm1, p_gradm2, p_gradm3, p_m + lm, this->m_Opt->m_FFT.plan, &xyz, timer);
        this->m_Opt->StopTimer(FFTSELFEXEC);
        this->m_Opt->IncrementCounter(FFT, FFTGRAD);

#pragma omp parallel
{
#pragma omp for
        for (IntType i = 0; i < nl; ++i) {  // for all grid points
            ltilde = p_ltilde[ll + i];
            // compute bodyforce
            p_bt1[i] += 0.5*scale*p_gradm1[i]*ltilde/static_cast<ScalarType>(nc);
            p_bt2[i] += 0.5*scale*p_gradm2[i]*ltilde/static_cast<ScalarType>(nc);
            p_bt3[i] += 0.5*scale*p_gradm3[i]*ltilde/static_cast<ScalarType>(nc);
        }
}  // omp
    }

    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField3, &p_ltildex); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField2, &p_divvx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_bt1, p_bt2, p_bt3); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(timer);
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the incremental adjoint problem (incremental
 * adjoint equation)
 * -\p_t \tilde{\lambda} - \idiv \tilde{\lambda}\vect{v}
 *                       - \idiv \lambda\tilde{\vect{v}} = 0
 * subject to \tilde{\lambda}_1 + \tilde{m}_1 = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode CLAIRE::SolveIncAdjointEquationFNSL(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = ThrowError("not implemented"); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply projection to map \tilde{v} onto the manifold
 * of divergence free velocity fields
 *******************************************************************/
PetscErrorCode CLAIRE::ApplyProjection() {
    PetscErrorCode ierr = 0;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief finalize the current iteration
 *******************************************************************/
PetscErrorCode CLAIRE::FinalizeIteration(Vec v) {
    PetscErrorCode ierr = 0;
    int rank;
    IntType nl, ng, nc, nt, iter;
    std::string filename, ext;
    std::stringstream ss;
    std::ofstream logwriter;
    ScalarType *p_m1 = NULL, *p_m = NULL, rval, dval, jval, gval;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get number of time points and grid points
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // parse extension
    ext = this->m_Opt->m_FileNames.extension;

    // if not yet allocted, do so
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    // set velocity field
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // log objective values
    if (this->m_Opt->m_Log.enabled[LOGCONV]) {
        iter = this->m_Opt->GetCounter(ITERATIONS);
        jval = this->m_Opt->m_Monitor.jval;
        dval = this->m_Opt->m_Monitor.dval;
        rval = this->m_Opt->m_Monitor.rval;
        this->m_Opt->LogConvergence(iter, jval, dval, rval);
    }

    // log norm of gradient
    if (this->m_Opt->m_Log.enabled[LOGGRAD]) {
        gval = this->m_Opt->m_Monitor.gradnorm;
        this->m_Opt->m_Log.gradnorm.push_back(gval);
    }

    // store iterates
    if (this->m_Opt->m_ReadWriteFlags.iterates) {
        // allocate
        if (this->m_WorkScaFieldMC == NULL) {
            ierr = VecCreate(this->m_WorkScaFieldMC, nl*nc, ng*nc); CHKERRQ(ierr);
        }

        iter = this->m_Opt->GetCounter(ITERATIONS);
        ierr = Assert(iter >= 0, "problem in counter"); CHKERRQ(ierr);

        // copy memory for m_1
        ierr = GetRawPointer(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception& err) {
            ierr = ThrowError(err); CHKERRQ(ierr);
        }
        ierr = RestoreRawPointer(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

        ss  << "deformed-template-image-i=" << std::setw(3) << std::setfill('0') << iter << ext;
        ierr = this->m_ReadWrite->Write(this->m_WorkScaFieldMC, ss.str(), nc > 1); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        // construct file names for velocity field components
        ss  << "velocity-field-i=" << std::setw(3) << std::setfill('0') << iter << ext;
        filename = ss.str();
        ss.str(std::string()); ss.clear();

        // velocity field out
        ierr = this->m_ReadWrite->Write(this->m_VelocityField, filename); CHKERRQ(ierr);
    }  // store iterates


    if (this->m_Opt->m_StoreCheckPoints) {
        ierr = this->m_ReadWrite->Write(this->m_VelocityField, "velocity-field-checkpoint"+ext); CHKERRQ(ierr);
    }


    // compute determinant of deformation gradient and write it to file
    if (this->m_Opt->m_Monitor.detdgradenabled) {
        ierr = this->ComputeDetDefGrad(); CHKERRQ(ierr);

        // if user enabled the logger
        if (this->m_Opt->m_Log.enabled[LOGJAC]) {
            if (rank == 0) {
                filename  = this->m_Opt->m_FileNames.xfolder;
                filename += "registration-performance-detdefgrad.log";

                // create output file or append to output file
                logwriter.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
                ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);
                ss  << std::scientific
                    <<  "iter = "     << this->m_Opt->GetCounter(ITERATIONS)
                    <<  "   betav = " << this->m_Opt->m_RegNorm.beta[0] << "    "
                    << std::left << std::setw(20) << this->m_Opt->m_Monitor.detdgradmin << " "
                                 << std::setw(20) << this->m_Opt->m_Monitor.detdgradmean <<" "
                                 << std::setw(20) << this->m_Opt->m_Monitor.detdgradmax;
                logwriter << ss.str() << std::endl;
                ss.str(std::string()); ss.clear();
            }  // if on master rank
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief finalize the registration
 *******************************************************************/
PetscErrorCode CLAIRE::Finalize(VecField* v) {
    PetscErrorCode ierr = 0;
    std::string filename, fn, ext;
    IntType nl, ng, nc, nt;
    int rank, nproc;
    std::ofstream logwriter;
    std::stringstream ss, ssnum;
    ScalarType value;
    ScalarType *p_m1 = NULL, *p_mt = NULL, *p_mr = NULL, *p_m = NULL, *p_dr = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

    // get sizes
    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    if (this->m_Opt->m_Verbosity >= 2) {
        ierr = DbgMsg("finalizing registration"); CHKERRQ(ierr);
    }

    // if not yet allocted, do so and copy input
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }
    ierr = this->m_VelocityField->Copy(v); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaFieldMC == NULL) {
        ierr = VecCreate(this->m_WorkScaFieldMC, nl*nc, ng*nc); CHKERRQ(ierr);
    }

    // process timer
    ierr = this->m_Opt->ProcessTimers(); CHKERRQ(ierr);

    // parse extension
    ext = this->m_Opt->m_FileNames.extension;

    // compute residuals
    if (this->m_Opt->m_Log.enabled[LOGDIST]) {
        ierr = VecWAXPY(this->m_WorkScaFieldMC, -1.0, this->m_TemplateImage, this->m_ReferenceImage); CHKERRQ(ierr);

        ierr = VecNorm(this->m_WorkScaFieldMC, NORM_2, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(0, value);

        ierr = VecNorm(this->m_WorkScaFieldMC, NORM_INFINITY, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(1, value);

        // deformed template out (compute solution of state equation)
        ierr = this->SolveStateEquation(); CHKERRQ(ierr);

        // copy memory for m_1
        ierr = GetRawPointer(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception& err) {
            ierr = ThrowError(err); CHKERRQ(ierr);
        }
        ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        ierr = RestoreRawPointer(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);

        // ||m_R - m_1||
        ierr = VecAXPY(this->m_WorkScaFieldMC, -1.0, this->m_ReferenceImage); CHKERRQ(ierr);

        ierr = VecNorm(this->m_WorkScaFieldMC, NORM_2, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(2, value);

        ierr = VecNorm(this->m_WorkScaFieldMC, NORM_INFINITY, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(3, value);
    }

    // write deformed template image to file
    if (this->m_Opt->m_ReadWriteFlags.deftemplate) {
        // copy memory for m_1
        ierr = GetRawPointer(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception& err) {
            ierr = ThrowError(err); CHKERRQ(ierr);
        }
        ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        ierr = RestoreRawPointer(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);

        if (this->m_Opt->m_RegFlags.registerprobmaps) {
            ierr = EnsurePartitionOfUnity(this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = ShowValues(this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = ComputeBackGround(this->m_WorkScaField1, this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->WriteT(this->m_WorkScaField1, "background-image" + ext, false); CHKERRQ(ierr);
        }

        ierr = this->m_ReadWrite->WriteT(this->m_WorkScaFieldMC, "deformed-template-image" + ext, nc > 1); CHKERRQ(ierr);
    }

    // write residual images to file
    if (this->m_Opt->m_ReadWriteFlags.residual || this->m_Opt->m_ReadWriteFlags.invresidual) {
        ierr = GetRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);

        // write residual at t = 0 to file
        ierr = GetRawPointer(this->m_TemplateImage, &p_mt); CHKERRQ(ierr);
        if (this->m_Opt->m_ReadWriteFlags.residual) {
            ierr = GetRawPointer(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            for (IntType i = 0; i < nl*nc; ++i) {
                p_dr[i] = PetscAbs(p_mt[i] - p_mr[i]);
            }
            ierr = RestoreRawPointer(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->Write(this->m_WorkScaFieldMC, "residual-t=0" + ext, nc > 1); CHKERRQ(ierr);
        }

        if (this->m_Opt->m_ReadWriteFlags.invresidual) {
            ierr = GetRawPointer(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            for (IntType i = 0; i < nl*nc; ++i) {
                p_dr[i] = 1.0 - PetscAbs(p_mt[i] - p_mr[i]);
            }
            ierr = RestoreRawPointer(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->Write(this->m_WorkScaFieldMC, "inv-residual-t=0" + ext, nc > 1); CHKERRQ(ierr);
        }
        ierr = RestoreRawPointer(this->m_TemplateImage, &p_mt); CHKERRQ(ierr);


        // write residual at t = 1 to file
        ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        if (this->m_Opt->m_ReadWriteFlags.residual) {
            ierr = GetRawPointer(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            for (IntType i = 0; i < nl*nc; ++i) {
                p_dr[i] = PetscAbs(p_m[nt*nl*nc + i] - p_mr[i]);
            }
            ierr = RestoreRawPointer(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
//            ierr = ShowValues(this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->Write(this->m_WorkScaFieldMC, "residual-t=1" + ext, nc > 1); CHKERRQ(ierr);
        }
        if (this->m_Opt->m_ReadWriteFlags.invresidual) {
            ierr = GetRawPointer(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
            for (IntType i = 0; i < nl*nc; ++i) {
                p_dr[i] = 1.0 - PetscAbs(p_m[nt*nl*nc + i] - p_mr[i]);
            }
            ierr = RestoreRawPointer(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);

//            ierr = ShowValues(this->m_WorkScaFieldMC, nc); CHKERRQ(ierr);
            ierr = this->m_ReadWrite->Write(this->m_WorkScaFieldMC, "inv-residual-t=1" + ext, nc > 1); CHKERRQ(ierr);
        }
        ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

        // restore reference image
        ierr = RestoreRawPointer(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    }

    // write velocity field to file
    if (this->m_Opt->m_ReadWriteFlags.velocity) {
        ierr = this->m_ReadWrite->Write(this->m_VelocityField, "velocity-field" + ext); CHKERRQ(ierr);
    }

    // write norm of velocity field to file
    if (this->m_Opt->m_ReadWriteFlags.velnorm) {
        ierr = this->m_VelocityField->Norm(this->m_WorkScaField1); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkScaField1, "velocity-field-norm" + ext); CHKERRQ(ierr);
    }

    // write determinant of deformation gradient to file
    if (this->m_Opt->m_ReadWriteFlags.detdefgrad) {
        ierr = this->ComputeDetDefGrad(true); CHKERRQ(ierr);
    }

    // write determinant of deformation gradient to file
    if (this->m_Opt->m_ReadWriteFlags.defgrad) {
        if (this->m_DeformationFields == NULL) {
            ierr = this->SetupDeformationField(); CHKERRQ(ierr);
        }
        ierr = this->m_DeformationFields->ComputeDefGrad(true); CHKERRQ(ierr);
    }

    // write deformation map to file
    if (this->m_Opt->m_ReadWriteFlags.defmap) {
        if (this->m_DeformationFields == NULL) {
            ierr = this->SetupDeformationField(); CHKERRQ(ierr);
        }
        ierr = this->m_DeformationFields->ComputeDeformationMap(true); CHKERRQ(ierr);
    }

    // write deformation field to file
    if (this->m_Opt->m_ReadWriteFlags.deffield) {
        if (this->m_DeformationFields == NULL) {
            ierr = this->SetupDeformationField(); CHKERRQ(ierr);
        }
        ierr = this->m_DeformationFields->ComputeDisplacementField(true); CHKERRQ(ierr);
    }

    // write template and reference image
    if (this->m_Opt->m_ReadWriteFlags.templateim) {
//        ierr = ShowValues(this->m_TemplateImage, nc); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->WriteT(this->m_TemplateImage, "template-image"+ext, nc > 1); CHKERRQ(ierr);
    }
    if (this->m_Opt->m_ReadWriteFlags.referenceim) {
//        ierr = ShowValues(this->m_ReferenceImage, nc); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->WriteR(this->m_ReferenceImage, "reference-image"+ext, nc > 1); CHKERRQ(ierr);
    }

    // write log file
    ierr = this->m_Opt->WriteLogFile(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _CLAIRE_CPP_
