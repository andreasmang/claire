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

#ifndef _OPTIMALCONTROLREGISTRATION_CPP_
#define _OPTIMALCONTROLREGISTRATION_CPP_

// global includes
#include <string>
#include <algorithm>

// local includes
#include "OptimalControlRegistration.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
OptimalControlRegistration::OptimalControlRegistration() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
OptimalControlRegistration::OptimalControlRegistration(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
OptimalControlRegistration::~OptimalControlRegistration() {
    this->ClearMemory();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::Initialize() {
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
PetscErrorCode OptimalControlRegistration::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // delete all variables
    ierr = this->ClearVariables(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::ClearVariables(void) {
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
 * @brief initialize the optimization (we essentially evaluate
 * the objective functional and the gradient for a given initial
 * guess)
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::InitializeOptimization(VecField* v0) {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    std::stringstream ss;
    ScalarType value, hd, alpha, jvt, jv, lsred, descent;
    Vec g = NULL, dv = NULL, v = NULL, vtilde = NULL;
    bool lssuccess;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    // if velocity field is null pointer, we did not set
    // any initial guess
    if (this->m_VelocityField == NULL) {
        if (this->m_Opt->GetVerbosity() > 2) {
            ierr = DbgMsg("allocating velocity field"); CHKERRQ(ierr);
        }
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    ierr = VecCreate(g, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecCreate(v, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecSet(v, 0.0); CHKERRQ(ierr);

    // if we use a non-zero initial guess, we compute
    // the first velocity using a steepest descent approach
    if (!this->m_Opt->GetOptPara().usezeroinitialguess) {
        ierr = VecCreate(dv, 3*nl, 3*ng); CHKERRQ(ierr);
        ierr = VecCreate(vtilde, 3*nl, 3*ng); CHKERRQ(ierr);

        lsred = 1E-4;  // reduction rate for line search
        for (int l = 0; l < 1; ++l) {
            // evaluate objective function
            ierr = this->EvaluateObjective(&jv, v); CHKERRQ(ierr);

            // compute gradient
            ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);

            // compute search direction
            ierr = this->EvaluatePrecondGradient(dv, v); CHKERRQ(ierr);

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
                if (this->m_Opt->GetVerbosity() > 1) {
                    ss << "line search successful (initialization; alpha=" << std::scientific << alpha << ")";
                    ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                    ss.clear(); ss.str(std::string());
                }
                ierr = VecCopy(vtilde, v); CHKERRQ(ierr);
            } else {
                if (this->m_Opt->GetVerbosity() > 1) {
                    ss << "line search failed (initialization; alpha=" << std::scientific << alpha << ")";
                    ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
                    ss.clear(); ss.str(std::string());
                }
                break;
            }
        }
    }
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // get lebesque measure
    hd = this->m_Opt->GetLebesqueMeasure();

    // evaluate distance measure
    ierr = this->EvaluateDistanceMeasure(&value); CHKERRQ(ierr);
    this->m_InitDistanceValue = hd*value;

    // evaluate objective functional
    ierr = this->EvaluateObjective(&value, v); CHKERRQ(ierr);
//    this->m_InitObjectiveVal = hd*value;
    this->m_InitObjectiveValue = value;

    // compute gradient
    ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);

    // compute gradient norm
    ierr = VecNorm(g, NORM_2, &value); CHKERRQ(ierr);
    this->m_InitGradientNorm = value;

    if (this->m_Opt->GetVerbosity() > 0) {
        ss << "initial gradient norm: "<< std::scientific << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    if (v0 != NULL) {
//      TODO: what's going on here
//    ierr = VecSet(v, 0.0); CHKERRQ(ierr);
//    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);
//    ierr = v0->Copy(this->m_VelocityField); CHKERRQ(ierr);
//    ierr = v0->SetValue(0.0); CHKERRQ(ierr);
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
 * @brief solve the forward problem (we assume the user has
 * set the velocity field)
 * @param[in] m0 density/image at t=0
 * @param[out] m1 density/image at t=1 (solution of transport
 * equation)
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::SolveForwardProblem(Vec m1, Vec m0) {
    PetscErrorCode ierr = 0;
    ScalarType *p_m1 = NULL, *p_m = NULL;
    IntType nt, nl, nc;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m0 != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(m1 != NULL, "null pointer"); CHKERRQ(ierr);

    // set initial condition
    this->m_TemplateImage = m0;

    // compute solution of state equation
    ierr = this->SolveStateEquation(); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;

    // copy m(t=1) to m_1
    ierr = VecGetArray(m1, &p_m1); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
    catch (std::exception&) {
        ierr = ThrowError("copy failed"); CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecRestoreArray(m1, &p_m1); CHKERRQ(ierr);

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
PetscErrorCode OptimalControlRegistration::SolveAdjointProblem(Vec l0, Vec m1) {
    PetscErrorCode ierr = 0;
    ScalarType *p_m = NULL, *p_m1 = NULL, *p_l = NULL, *p_l0 = NULL;
    IntType nt, nl, nc, ng;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m1 != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    // allocate state variable
    if (this->m_StateVariable == NULL) {
        ierr = VecCreate(this->m_StateVariable, (nt+1)*nl*nc, (nt+1)*ng*nc); CHKERRQ(ierr);
    }

    // copy memory for m_1
    ierr = VecGetArray(m1, &p_m1); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    try {std::copy(p_m1, p_m1+nl*nc, p_m+nt*nl*nc);}
    catch (std::exception&) {
        ierr = ThrowError("copy failed"); CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecRestoreArray(m1, &p_m1); CHKERRQ(ierr);

    // compute solution of state equation
    ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);

    // copy memory for m_1
    ierr = VecGetArray(l0, &p_l0); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    try {std::copy(p_l, p_l+nl*nc, p_l0);}
    catch (std::exception&) {
        ierr = ThrowError("copy failed"); CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = VecRestoreArray(l0, &p_l0); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief set state variable from externally
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::SetStateVariable(Vec m) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(m != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

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
    if (this->m_Opt->GetPDESolverPara().type == SL) {
        ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
        if (this->m_SemiLagrangianMethod == NULL) {
            try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        if (this->m_WorkVecField1 == NULL) {
            try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        // compute trajectory
        ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1, "state"); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set state variable from externally
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::GetStateVariable(Vec& m) {
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
PetscErrorCode OptimalControlRegistration::SetAdjointVariable(Vec lambda) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nc = this->m_Opt->GetDomainPara().nc;
    nt = this->m_Opt->GetDomainPara().nt;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    // we have to allocate the variable, because we delete it
    // at the end once we're done; since it comes from external
    // we need to make sure that we don't delete the external pointer
    if (this->m_AdjointVariable == NULL) {
        ierr = VecCreate(this->m_AdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }

    ierr = VecCopy(lambda, this->m_AdjointVariable); CHKERRQ(ierr);

    if (this->m_Opt->GetPDESolverPara().type == SL) {
        ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
        if (this->m_SemiLagrangianMethod == NULL) {
            try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        if (this->m_WorkVecField1 == NULL) {
            try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        // compute trajectory
        ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
        ierr = this->m_WorkVecField1->Scale(-1.0); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1, "adjoint"); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set state variable from externally
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::GetAdjointVariable(Vec& lambda) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(lambda == NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    lambda = this->m_AdjointVariable;
//    ierr = VecCopy(this->m_AdjointVariable,lambda); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluate the l2 distance between m_R and m_1
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::EvaluateDistanceMeasure(ScalarType* D) {
    PetscErrorCode ierr = 0;
    ScalarType *p_mr = NULL, *p_m = NULL;
    IntType nt, nc, nl, l;
    int rval;
    ScalarType dr, value, l2distance;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;

    // compute solution of state equation
    ierr = this->SolveStateEquation(); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    l = nt*nl*nc;
    value = 0.0;
    for (IntType i = 0; i < nc*nl; ++i) {
        dr = (p_mr[i] - p_m[l+i]);
        value += dr*dr;
    }

    rval = MPI_Allreduce(&value, &l2distance, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

    // objective value
    *D = 0.5*l2distance/static_cast<ScalarType>(nc);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluates the objective value
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::EvaluateObjective(ScalarType* J, Vec v) {
    PetscErrorCode ierr = 0;
    ScalarType D = 0.0, R = 0.0, hd;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // get lebesque measure
    hd = this->m_Opt->GetLebesqueMeasure();

    // allocate
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
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
        ierr = this->m_Regularization->EvaluateFunctional(&R, this->m_VelocityField); CHKERRQ(ierr);
    }

    // add up the contributions
    *J = hd*(D + R);

    // store for access
    this->m_Opt->SetJVal(*J);
    this->m_Opt->SetDVal(hd*D);
    this->m_Opt->SetRVal(hd*R);

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
PetscErrorCode OptimalControlRegistration::EvaluateGradient(Vec g, Vec v) {
    PetscErrorCode ierr = 0;
    ScalarType hd, value;
    std::stringstream ss;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("evaluating gradient"); CHKERRQ(ierr);
    }

    // start timer
    ierr = this->m_Opt->StartTimer(GRADEXEC); CHKERRQ(ierr);

    // parse input arguments
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);
    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = VecNorm(v, NORM_2, &value); CHKERRQ(ierr);
        ss << "||v||_2 = " << std::scientific << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    // compute solution of adjoint equation (i.e., \lambda(x,t))
    ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);

    // compute body force \int_0^1 grad(m)\lambda dt
    // assigned to work vecfield 2
    ierr = this->ComputeBodyForce(); CHKERRQ(ierr);

    // evaluate gradient of regularization model
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        // \vect{g}_v = \D{K}[\vect{b}]
        if (this->m_Opt->GetVerbosity() > 2) {
            ierr = DbgMsg("zero velocity field"); CHKERRQ(ierr);
        }
        ierr = this->m_WorkVecField2->GetComponents(g); CHKERRQ(ierr);
    } else {
        // evaluate / apply gradient operator for regularization
        ierr = this->m_Regularization->EvaluateGradient(this->m_WorkVecField1, this->m_VelocityField); CHKERRQ(ierr);

        // \vect{g}_v = \beta_v \D{A}[\vect{v}] + \D{K}[\vect{b}]
        ierr = this->m_WorkVecField1->AXPY(1.0, this->m_WorkVecField2); CHKERRQ(ierr);

        // parse to output
        ierr = this->m_WorkVecField1->GetComponents(g); CHKERRQ(ierr);
    }

    // get and scale by lebesque measure
    hd = this->m_Opt->GetLebesqueMeasure();
    ierr = VecScale(g, hd); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = VecNorm(g, NORM_2, &value); CHKERRQ(ierr);
        ss << "||g||_2 = " << std::scientific << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }

    // stop timer
    ierr = this->m_Opt->StopTimer(GRADEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(GRADEVAL);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute the body force
 * b = \int_0^1 grad(m) \lambda dt
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::ComputeBodyForce() {
    PetscErrorCode ierr = 0;
    IntType nt, nl, nc, l;
    std::stringstream ss;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    ScalarType *p_mt = NULL, *p_m = NULL, *p_l = NULL,
               *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL,
               *p_b1 = NULL, *p_b2 = NULL, *p_b3 = NULL;
    ScalarType ht, scale, lambda, value;
    double timers[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // check for null pointers
    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    // get problem dimensions and weights
    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ht = this->m_Opt->GetTimeStepSize();
    scale = ht;

    ierr = Assert(nt > 0, "nt<=0"); CHKERRQ(ierr);
    ierr = Assert(ht > 0, "ht<=0"); CHKERRQ(ierr);

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // init body force for numerical integration
    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        if (this->m_Opt->GetVerbosity() > 1) {
            ierr = DbgMsg("zero velocity field"); CHKERRQ(ierr);
        }
        // m and \lambda are constant in time
        ierr = VecGetArray(this->m_TemplateImage, &p_mt); CHKERRQ(ierr);
        for (IntType k = 0; k < nc; ++k) {  // for all components
            // compute gradient of m
            accfft_grad(p_gradm1, p_gradm2, p_gradm3, p_mt+k*nl, this->m_Opt->GetFFT().plan, &XYZ, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            // b = \sum_k\int_{\Omega} \lambda_k \grad m_k dt
            for (IntType i = 0; i < nl; ++i) {
                lambda = p_l[k*nl+i];
                p_b1[i] += lambda*p_gradm1[i]/static_cast<ScalarType>(nc);
                p_b2[i] += lambda*p_gradm2[i]/static_cast<ScalarType>(nc);
                p_b3[i] += lambda*p_gradm3[i]/static_cast<ScalarType>(nc);
            }
        }
        ierr = VecRestoreArray(this->m_TemplateImage, &p_mt); CHKERRQ(ierr);
    } else {  // non zero velocity field
        ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j) {
            // scaling for trapezoidal rule
            if ((j == 0) || (j == nt)) scale *= 0.5;
            for (IntType k = 0; k < nc; ++k) {  // for all components
                l = j*nl*nc + k*nl;

                // grad(m^j)
                accfft_grad(p_gradm1, p_gradm2, p_gradm3, p_m+l, this->m_Opt->GetFFT().plan, &XYZ, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                // \vect{b}_i += h_d*ht*\lambda^j (\grad m^j)_i
                for (IntType i = 0; i < nl; ++i) {
                    lambda = p_l[l+i];  // get \lambda(x_i,t^j)
                    p_b1[i] += scale*p_gradm1[i]*lambda/static_cast<ScalarType>(nc);
                    p_b2[i] += scale*p_gradm2[i]*lambda/static_cast<ScalarType>(nc);
                    p_b3[i] += scale*p_gradm3[i]*lambda/static_cast<ScalarType>(nc);
                }
            }
            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)) scale *= 2.0;
        }
        // restore arrays
        ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    }  // else zero velocity field

    ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);   // adjoint variable for all t^j
    ierr = this->m_WorkVecField1->RestoreArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArrays(p_b1, p_b2, p_b3); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = this->m_WorkVecField2->Norm(value); CHKERRQ(ierr);
        ss << "||b||_2 = " << std::scientific << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
    }



    this->m_Opt->IncreaseFFTTimers(timers);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the hessian to a vector
 * @param[in] vtilde incremental velocity field
 * @param[in] scale flag to switch on scaling by lebesque measure
 * @param[out] Hvtilde hessian applied to vector
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::HessianMatVec(Vec Hvtilde, Vec vtilde, bool scale) {
    PetscErrorCode ierr = 0;
    ScalarType hd, gamma;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);


    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("computing hessian matvec"); CHKERRQ(ierr);
    }

    ierr = this->m_Opt->StartTimer(HMVEXEC); CHKERRQ(ierr);

    // switch between hessian operators
    switch (this->m_Opt->GetKrylovSolverPara().matvectype) {
        case DEFAULTMATVEC:
        {
            // apply hessian H to \tilde{v}
            ierr = this->HessMatVec(Hvtilde, vtilde); CHKERRQ(ierr);
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
            ierr = ThrowError("setup problem"); CHKERRQ(ierr);
            break;
        }
    }


    // scale by lebesque measure
    if (scale) {
        hd = this->m_Opt->GetLebesqueMeasure();
        ierr = VecScale(Hvtilde, hd); CHKERRQ(ierr);
    }

    gamma = this->m_Opt->GetKrylovSolverPara().hessshift;
    if (gamma > 0.0) {
        ierr = VecAXPY(Hvtilde, gamma, vtilde); CHKERRQ(ierr);
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
PetscErrorCode OptimalControlRegistration::HessMatVec(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate container for incremental velocity field
    if (this->m_IncVelocityField == NULL) {
        try {this->m_IncVelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    // parse input
    ierr = this->m_IncVelocityField->SetComponents(vtilde); CHKERRQ(ierr);

    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t)
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // compute incremental body force
    ierr = this->ComputeIncBodyForce(); CHKERRQ(ierr);

    // apply 2nd variation of regularization model to
    // incremental control variable: \beta*\D{A}[\vect{\tilde{v}}]
    ierr = this->m_Regularization->HessianMatVec(this->m_WorkVecField1, this->m_IncVelocityField); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \beta*\D{A}[\vect{\tilde{v}}] + \D{K}[\vect{\tilde{b}}]
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    ierr = this->m_WorkVecField1->AXPY(1.0, this->m_WorkVecField2); CHKERRQ(ierr);

    // pass to output
    ierr = this->m_WorkVecField1->GetComponents(Hvtilde); CHKERRQ(ierr);

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
PetscErrorCode OptimalControlRegistration::PrecondHessMatVec(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate container for incremental velocity field
    if (this->m_IncVelocityField == NULL) {
        try {this->m_IncVelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    // parse input
    ierr = this->m_IncVelocityField->SetComponents(vtilde); CHKERRQ(ierr);

    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t)
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // compute incremental body force
    ierr = this->ComputeIncBodyForce(); CHKERRQ(ierr);

    // apply inverse of 2nd variation of regularization model to
    // incremental body force: (\beta \D{A})^{-1}\D{K}[\vect{\tilde{b}}]
    ierr = this->m_Regularization->ApplyInvOp(this->m_WorkVecField1, this->m_WorkVecField2); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \vect{\tilde{v}} + (\beta \D{A})^{-1} \D{K}[\vect{\tilde{b}}]
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    ierr = this->m_WorkVecField2->WAXPY(1.0, this->m_IncVelocityField, this->m_WorkVecField1); CHKERRQ(ierr);

    // pass to output
    ierr = this->m_WorkVecField2->GetComponents(Hvtilde); CHKERRQ(ierr);

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
PetscErrorCode OptimalControlRegistration::PrecondHessMatVecSym(Vec Hvtilde, Vec vtilde) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    // allocate container for incremental velocity field
    if (this->m_IncVelocityField == NULL) {
        try {this->m_IncVelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    // allocate work vec field 5 (1,2,3, and 4 are overwritten
    // during the computation of the incremental forward and adjoint
    // solve and the computation of the incremental body force)
    if (this->m_WorkVecField5 == NULL) {
        try {this->m_WorkVecField5 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // parse input (store incremental velocity field \tilde{v})
    ierr = this->m_WorkVecField5->SetComponents(vtilde); CHKERRQ(ierr);

    // apply inverse of 2nd variation of regularization model to
    // incremental body force: (\beta\D{A})^{-1/2}\D{K}[\vect{\tilde{b}}](\beta\D{A})^{-1/2}

    // apply (\beta\D{A})^{-1/2} to incremental velocity field
    ierr = this->m_Regularization->ApplyInvOp(this->m_IncVelocityField, this->m_WorkVecField5, true); CHKERRQ(ierr);

    // now solve the PDEs given the preconditoined incremental velocity field

    // compute \tilde{m}(x,t)
    ierr = this->SolveIncStateEquation(); CHKERRQ(ierr);

    // compute \tilde{\lambda}(x,t)
    ierr = this->SolveIncAdjointEquation(); CHKERRQ(ierr);

    // compute incremental body force
    ierr = this->ComputeIncBodyForce(); CHKERRQ(ierr);

    // apply (\beta\D{A})^{-1/2} to incremental body force
    ierr = this->m_Regularization->ApplyInvOp(this->m_WorkVecField1, this->m_WorkVecField2, true); CHKERRQ(ierr);

    // \D{H}\vect{\tilde{v}} = \vect{\tilde{v}} + (\beta \D{A})^{-1/2}\D{K}[\vect{\tilde{b}}](\beta \D{A})^{-1/2}
    // we use the same container for the bodyforce and the incremental body force to
    // save some memory
    ierr = this->m_WorkVecField5->AXPY(1.0, this->m_WorkVecField1); CHKERRQ(ierr);

    // pass to output
    ierr = this->m_WorkVecField5->GetComponents(Hvtilde); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute initial condition given some initial guess for
 * the state variable $m$ and the adjoint variable $\lambda$
 * @param[in] m state variable
 * @param[in] lambda adjoint variable
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::ComputeInitialCondition(Vec m, Vec lambda) {
    PetscErrorCode ierr = 0;
    IntType nt, nl, nc, ng;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    ext = this->m_Opt->GetReadWriteFlags().extension;

    // allocate container for incremental velocity field
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    // allocate state and adjoint variables
    if (this->m_StateVariable == NULL) {
        ierr = VecCreate(this->m_StateVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }

    // allocate state and adjoint variables
    if (this->m_AdjointVariable == NULL) {
        ierr = VecCreate(this->m_AdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("piccard iteration"); CHKERRQ(ierr);
    }

    // copy input state and adjoint variable to class variables
    ierr = VecCopy(m, this->m_StateVariable); CHKERRQ(ierr);
    ierr = VecCopy(lambda, this->m_AdjointVariable); CHKERRQ(ierr);

    // compute body force (assigned to work vec field 2)
    ierr = this->ComputeBodyForce(); CHKERRQ(ierr);

    // piccard step: solve A[v] = - ht \sum_j \lambda^j grad(m^j)
    ierr = this->m_WorkVecField2->Scale(-1.0); CHKERRQ(ierr);

    ierr = this->m_Regularization->ApplyInvOp(this->m_VelocityField,
                                              this->m_WorkVecField2); CHKERRQ(ierr);

    // reset the adjoint variables
    ierr = VecSet(this->m_StateVariable, 0.0); CHKERRQ(ierr);
    ierr = VecSet(this->m_AdjointVariable, 0.0); CHKERRQ(ierr);

    ierr = this->m_ReadWrite->Write(this->m_VelocityField, "initial-condition"+ext); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief evaluate preconditioned gradient of lagrangian functional
 * @param[in] v velocity field
 * @param[out] g gradient
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::EvaluatePrecondGradient(Vec g, Vec v) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(g != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate container for incremental velocity field
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // allocate regularization model
    if (this->m_Regularization == NULL) {
        ierr = this->AllocateRegularization(); CHKERRQ(ierr);
    }

    if (this->m_Opt->GetVerbosity() > 2) {
        ierr = DbgMsg("evaluating preconditioned gradient"); CHKERRQ(ierr);
    }

    // parse input arguments
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    // compute solution of state equation (i.e., m(x,t))
    ierr = this->SolveStateEquation(); CHKERRQ(ierr);

    // compute solution of adjoint equation (i.e., \lambda(x,t))
    ierr = this->SolveAdjointEquation(); CHKERRQ(ierr);

    // compute body force (assigned to work vec field 2)
    ierr = this->ComputeBodyForce(); CHKERRQ(ierr);

    // piccard step: solve A^{-1}[ht \sum_j \lambda^j grad(m^j)]
    ierr = this->m_Regularization->ApplyInvOp(this->m_WorkVecField1,
                                              this->m_WorkVecField2); CHKERRQ(ierr);

    // compute v + A^{-1}[ht \sum_j \lambda^j grad(m^j)]
    ierr = this->m_WorkVecField1->AXPY(1.0, this->m_VelocityField); CHKERRQ(ierr);

    // parse output arguments
    ierr = this->m_WorkVecField1->GetComponents(g); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(GRADEVAL);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute the incremental body force
 * \tilde{\vect{b}} = \int_0^1 \igrad\tilde{m}\lambda
 *                           + \igrad m\tilde{\lambda} dt
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::ComputeIncBodyForce() {
    PetscErrorCode ierr = 0;
    IntType nt, nl, nc, l;
    ScalarType *p_m = NULL, *p_mt = NULL, *p_l = NULL, *p_lt = NULL,
                *p_bt1 = NULL, *p_bt2 = NULL, *p_bt3 = NULL,
                *p_gradm1 = NULL, *p_gradm2 = NULL, *p_gradm3 = NULL,
                *p_gradmt1 = NULL, *p_gradmt2 = NULL, *p_gradmt3 = NULL;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    ScalarType ht, scale, lj, ltj;
    double timers[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ht = this->m_Opt->GetTimeStepSize();
    ierr = Assert(nt > 0, "nt <= 0"); CHKERRQ(ierr);
    ierr = Assert(ht > 0, "ht <= 0"); CHKERRQ(ierr);
    ierr = Assert(nc > 0, "nc <= 0"); CHKERRQ(ierr);
    scale = ht;


    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // init array
    ierr = this->m_WorkVecField2->SetValue(0.0); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->GetArrays(p_gradm1, p_gradm2, p_gradm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArrays(p_bt1, p_bt2, p_bt3); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);  // state variable for all t^j
    ierr = VecGetArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);  // incremental adjoint variable for all t^j

    // TODO(andreas): add case for zero velocity field
    if (this->m_Opt->GetOptPara().method == FULLNEWTON) {
        ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);

        if (this->m_WorkVecField3 == NULL) {
            try {this->m_WorkVecField3 = new VecField(this->m_Opt);}
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);  // adjoint variable for all t^j
        ierr = VecGetArray(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);  // incremental state variable for all t^j

        ierr = this->m_WorkVecField3->GetArrays(p_gradmt1, p_gradmt2, p_gradmt3); CHKERRQ(ierr);

        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j) {  // for all time points
            // trapezoidal rule (apply scaling)
            if ((j == 0) || (j == nt)) scale *= 0.5;
            for (IntType k = 0; k < nc; ++k) {  // for all image components
                l = j*nl*nc + k*nl;

                // computing gradient of m
                accfft_grad(p_gradm1, p_gradm2, p_gradm3, p_m+l, this->m_Opt->GetFFT().plan, &XYZ, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                // computing gradient of \tilde{m}
                accfft_grad(p_gradmt1, p_gradmt2, p_gradmt3, p_mt+l, this->m_Opt->GetFFT().plan, &XYZ, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

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
            // trapezoidal rule (revert scaling)
            if ((j == 0) || (j == nt)) scale *= 2.0;
        }  // for all time points

        ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);  // adjoint variable for all t^j
        ierr = VecRestoreArray(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);  // incremental state variable for all t^j
        ierr = this->m_WorkVecField3->RestoreArrays(p_gradmt1, p_gradmt2, p_gradmt3); CHKERRQ(ierr);
    } else if (this->m_Opt->GetOptPara().method == GAUSSNEWTON) {  // gauss newton approximation
        // compute numerical integration (trapezoidal rule)
        for (IntType j = 0; j <= nt; ++j) {  // for all time points
            // trapezoidal rule (apply scaling)
            if ((j == 0) || (j == nt)) scale *= 0.5;
            for (IntType k = 0; k < nc; ++k) {  // for all image components
                l = j*nl*nc + k*nl;

                // compute gradient of m^j
                accfft_grad(p_gradm1, p_gradm2, p_gradm3, p_m+l, this->m_Opt->GetFFT().plan, &XYZ, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                // compute \vect{\tilde{b}}^k_i += h_d*ht*(\tilde{\lambda}^j (\grad m^j)^k
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    // get \tilde{\lambda}(x_i,t^j)
                    ltj = p_lt[l+i];
                    p_bt1[i] += scale*p_gradm1[i]*ltj/static_cast<ScalarType>(nc);
                    p_bt2[i] += scale*p_gradm2[i]*ltj/static_cast<ScalarType>(nc);
                    p_bt3[i] += scale*p_gradm3[i]*ltj/static_cast<ScalarType>(nc);
                }  // for all grid points
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

    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);  // state variable for all t^j
    ierr = VecRestoreArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);  // incremental adjoint variable for all t^j

    this->m_Opt->IncreaseFFTTimers(timers);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::SolveStateEquation(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt;
    ScalarType *p_m = NULL, *p_m0 = NULL, *p_mj = NULL;
    std::stringstream ss;
    std::string ext;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);

    // check cfl condition / update time step
    if (this->m_Opt->GetPDESolverPara().monitorcflnumber ||
        this->m_Opt->GetPDESolverPara().adapttimestep) {
        ierr = this->ComputeCFLCondition(); CHKERRQ(ierr);
    }

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ierr = Assert(nt > 0, "nt <= 0"); CHKERRQ(ierr);

    ext = this->m_Opt->GetReadWriteFlags().extension;

    if (this->m_Opt->GetVerbosity() > 2) {
        ss << "solving state equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->GetDomainPara().nx[0]
           << "," << this->m_Opt->GetDomainPara().nx[1]
           << "," << this->m_Opt->GetDomainPara().nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // allocate state and adjoint variables
    if (this->m_StateVariable == NULL) {
        ierr = VecCreate(this->m_StateVariable, (nt+1)*nl*nc, (nt+1)*ng*nc); CHKERRQ(ierr);
    }

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        // we copy m_0 to all t for v=0
        ierr = this->CopyToAllTimePoints(this->m_StateVariable, this->m_TemplateImage); CHKERRQ(ierr);
    } else {
        // copy initial condition m_0 = m_T
        ierr = VecGetArray(this->m_TemplateImage, &p_m0); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m0, p_m0+nl*nc, p_m);}
        catch (std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_TemplateImage, &p_m0); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

        // call the solver
        switch (this->m_Opt->GetPDESolverPara().type) {
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
    }  // velocity field is zero

    ierr = this->m_Opt->StopTimer(PDEEXEC); CHKERRQ(ierr);


    if (this->m_Opt->GetVerbosity() > 2) {
        ScalarType maxval, minval, value;
        ierr = VecMax(this->m_StateVariable, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(this->m_StateVariable, NULL, &minval); CHKERRQ(ierr);
        ss << "state variable: [" << std::scientific << minval << "," << maxval << "]";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        ierr = VecNorm(this->m_VelocityField->m_X1, NORM_2, &value); CHKERRQ(ierr);
        ss << "velocity norm: " << std::scientific << value;
        ierr = VecNorm(this->m_VelocityField->m_X2, NORM_2, &value); CHKERRQ(ierr);
        ss << " " << value;
        ierr = VecNorm(this->m_VelocityField->m_X3, NORM_2, &value); CHKERRQ(ierr);
        ss << " " << value;
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }


    // store time series
    if (this->m_Opt->GetReadWriteFlags().timeseries) {
        if (this->m_WorkScaField1 == NULL) {
            ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
        }
        ierr = Assert(this->m_ReadWrite != NULL, "null pointer"); CHKERRQ(ierr);

        // store time history
        ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        // store individual time points
        for (IntType j = 0; j <= nt; ++j) {
            for (IntType k = 0; k < nc; ++k) {
                ierr = VecGetArray(this->m_WorkScaField1, &p_mj); CHKERRQ(ierr);
                try {std::copy(p_m+j*nl*nc + k*nl, p_m+j*nl*nc + (k+1)*nl, p_mj);}
                catch(std::exception&) {
                    ierr = ThrowError("copying of data failed"); CHKERRQ(ierr);
                }
                ierr = VecRestoreArray(this->m_WorkScaField1, &p_mj); CHKERRQ(ierr);
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
        ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
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
PetscErrorCode OptimalControlRegistration::SolveStateEquationRK2(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, l, lnext;
    ScalarType *p_m = NULL, *p_mbar = NULL, *p_rhs0 = NULL,
                *p_gmx1 = NULL, *p_gmx2 = NULL, *p_gmx3 = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL;
    ScalarType ht = 0.0, hthalf = 0.0, rhs1;
    double timers[5] = {0, 0, 0, 0, 0};
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
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
    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField1, &p_mbar); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2, &p_rhs0); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {
        for (IntType k = 0; k < nc; ++k) {
            l = j*nl*nc + k*nl;
            lnext = (j+1)*nl*nc + k*nl;

            // compute gradient of k-th component of m_j
            accfft_grad(p_gmx1, p_gmx2, p_gmx3, p_m+l, this->m_Opt->GetFFT().plan, &XYZ, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            // evaluate right hand side and compute intermediate rk2 step
            for (IntType i = 0; i < nl; ++i) {
                 p_rhs0[i] = -p_gmx1[i]*p_vx1[i]
                             -p_gmx2[i]*p_vx2[i]
                             -p_gmx3[i]*p_vx3[i];

                 // compute intermediate result
                 p_mbar[i] = p_m[l+i] + ht*p_rhs0[i];
            }

            // compute gradient of \bar{m}
            accfft_grad(p_gmx1, p_gmx2, p_gmx3, p_mbar, this->m_Opt->GetFFT().plan, &XYZ, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            // evaluate right hand side and wrap up integration
            for (IntType i = 0; i < nl; ++i) {
                rhs1 = -p_gmx1[i]*p_vx1[i]
                       -p_gmx2[i]*p_vx2[i]
                       -p_gmx3[i]*p_vx3[i];

                // we have overwritten m_j with intermediate result
                // m_{j+1} = m_j + 0.5*ht*(RHS0 + RHS1)
                p_m[lnext+i] = p_m[l+i] + hthalf*(p_rhs0[i] + rhs1);
            }
        }  // for all components
    }  // for all time points

    // copy initial condition to buffer
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_mbar); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_rhs0); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->RestoreArrays(p_gmx1, p_gmx2, p_gmx3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timers);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the forward problem (state equation)
 * \p_t m + \igrad m\cdot\vect{v} = 0
 * subject to m_0 - m_T = 0
 * solved forward in time
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::SolveStateEquationSL(void) {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt, l, lnext;
    ScalarType *p_m = NULL;
    std::stringstream ss;
    std::string filename;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // compute trajectory
    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1, "state"); CHKERRQ(ierr);

    // get state variable m
    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            l = j*nl*nc + k*nl;
            lnext = (j+1)*nl*nc + k*nl;
            // compute m(X,t^{j+1}) (interpolate state variable)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_m+lnext, p_m+l, "state"); CHKERRQ(ierr);
        }
    }
    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::SolveAdjointEquation(void) {
    PetscErrorCode ierr = 0;
    IntType nl, nc, ng, nt, l;
    ScalarType *p_m = NULL, *p_l = NULL, *p_mr = NULL;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    if (this->m_Opt->GetVerbosity() > 2) {
        ss << "solving adjoint equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->GetDomainPara().nx[0]
           << "," << this->m_Opt->GetDomainPara().nx[1]
           << "," << this->m_Opt->GetDomainPara().nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }
    if (this->m_AdjointVariable == NULL) {
        ierr = VecCreate(this->m_AdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }

    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // compute terminal condition \lambda_1 = -(m_1 - m_R) = m_R - m_1
    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    l = nt*nc*nl;  // index for final condition
    for (IntType i = 0; i < nc*nl; ++i) {
        p_l[l+i] = (p_mr[i] - p_m[l+i]); // / static_cast<ScalarType>(nc);
    }
    ierr = VecRestoreArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
        // copy adjoint variable to all time points
        for (IntType j = 1; j <= nt; ++j) {
            try {std::copy(p_l+nt*nc*nl, p_l+(nt+1)*nc*nl, p_l+(nt-j)*nl*nc);}
            catch (std::exception&) {
                ierr = ThrowError("copy failed"); CHKERRQ(ierr);
            }
        }
        ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    } else {
        // call the solver
        switch (this->m_Opt->GetPDESolverPara().type) {
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

    if (this->m_Opt->GetVerbosity() > 2) {
        ScalarType maxval, minval;
        ierr = VecMax(this->m_AdjointVariable, NULL, &maxval); CHKERRQ(ierr);
        ierr = VecMin(this->m_AdjointVariable, NULL, &minval); CHKERRQ(ierr);
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
PetscErrorCode OptimalControlRegistration::SolveAdjointEquationRK2(void) {
    PetscErrorCode ierr;
    IntType nl, ng, nc, nt, l, lnext;
    ScalarType *p_l = NULL, *p_rhs0 = NULL, *p_rhs1 = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_ljvx1 = NULL, *p_ljvx2 = NULL, *p_ljvx3 = NULL;
    ScalarType hthalf, ht, lambdabar, lambda;
    double timers[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

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
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);

    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_ljvx1, p_ljvx2, p_ljvx3); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            l = (nt-j)*nc*nl + k*nl;
            lnext = (nt-(j+1))*nc*nl + k*nl;

            // scale \vect{v} by \lambda
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                lambda = p_l[l+i];
                p_ljvx1[i] = lambda*p_vx1[i];
                p_ljvx2[i] = lambda*p_vx2[i];
                p_ljvx3[i] = lambda*p_vx3[i];
            }  // for all grid points

            // compute \idiv(\lambda\vect{v})
            accfft_divergence(p_rhs0, p_ljvx1, p_ljvx2, p_ljvx3, this->m_Opt->GetFFT().plan, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                // compute \bar{\lambda} = \lambda_j + ht*\idiv(\lambda\vect{v})
                lambdabar = p_l[l+i] + ht*p_rhs0[i];
                // scale \vect{v} by \bar{\lambda}
                p_ljvx1[i] = p_vx1[i]*lambdabar;
                p_ljvx2[i] = p_vx2[i]*lambdabar;
                p_ljvx3[i] = p_vx3[i]*lambdabar;
            }

            // compute \idiv(\bar{\lambda}\vect{v})
            accfft_divergence(p_rhs1, p_ljvx1, p_ljvx2, p_ljvx3, this->m_Opt->GetFFT().plan, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                p_l[lnext+i] = p_l[l+i] + hthalf*(p_rhs0[i] + p_rhs1[i]);
            }
        }  // for all image components
    }  // for all time points

    ierr = this->m_WorkVecField1->RestoreArrays(p_ljvx1, p_ljvx2, p_ljvx3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timers);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief solve the adjoint problem (adjoint equation)
 * -\p_t \lambda - \idiv \lambda\vect{v} = 0
 * subject to \lambda_1 + (m_R - m_1) = 0
 * solved backward in time
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::SolveAdjointEquationSL() {
    PetscErrorCode ierr = 0;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_divv = NULL, *p_divvX = NULL, *p_l = NULL, *p_ljX = NULL;
    ScalarType ht, ljX, rhs0, rhs1;
    IntType nl, ng, nc, nt, l, lnext;
    double timers[5] = {0, 0, 0, 0, 0};

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();

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
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // scale v by -1
    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->Scale(-1.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);

    // compute \idiv(\tilde{\lambda}\vect{v})
    accfft_divergence(p_divv, p_vx1, p_vx2, p_vx3, this->m_Opt->GetFFT().plan, timers);
    this->m_Opt->IncrementCounter(FFT, 4);

    ierr = this->m_WorkVecField1->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_WorkVecField1, "adjoint"); CHKERRQ(ierr);

    // evaluate div(v) at X
    ierr = VecGetArray(this->m_WorkScaField2, &p_divvX); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_divvX, p_divv, "adjoint"); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField3, &p_ljX); CHKERRQ(ierr);
    for (IntType j = 0; j < nt; ++j) {
        for (IntType k = 0; k < nc; ++k) {
            l = (nt-j)*nc*nl + k*nl;
            lnext = (nt-(j+1))*nc*nl + k*nl;

            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_ljX, p_l+l, "adjoint"); CHKERRQ(ierr);

            for (IntType i = 0; i < nl; ++i) {
                ljX = p_ljX[i];

                rhs0 = -ljX*p_divvX[i];
                rhs1 = -(ljX + ht*rhs0)*p_divv[i];

                // compute \lambda(x,t^{j+1})
                p_l[lnext+i] = ljX + 0.5*ht*(rhs0 + rhs1);
            }
        }
    }
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_divvX); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField3, &p_ljX); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timers);

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
PetscErrorCode OptimalControlRegistration::SolveIncStateEquation(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2) {
        ss << "solving incremental state equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->GetDomainPara().nx[0]
           << "," << this->m_Opt->GetDomainPara().nx[1]
           << "," << this->m_Opt->GetDomainPara().nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    // allocate variables
    if (this->m_IncStateVariable == NULL) {
        ierr = VecCreate(this->m_IncStateVariable, (nt+1)*nl*nc, (nt+1)*ng*nc); CHKERRQ(ierr);
    }

    // start timer
    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // set initial value
    ierr = VecSet(this->m_IncStateVariable, 0.0); CHKERRQ(ierr);

    // call the solver
    switch (this->m_Opt->GetPDESolverPara().type) {
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

    if (this->m_Opt->GetVerbosity() > 2) {
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
PetscErrorCode OptimalControlRegistration::SolveIncStateEquationRK2(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, l, lnext;
    ScalarType *p_m = NULL, *p_mt = NULL, *p_mtbar = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_gmx1 = NULL, *p_gmx2 = NULL, *p_gmx3 = NULL,
                *p_gmtx1 = NULL, *p_gmtx2 = NULL, *p_gmtx3 = NULL,
                *p_vtx1 = NULL, *p_vtx2 = NULL, *p_vtx3 = NULL, *p_rhs0 = NULL;
    ScalarType ht, hthalf;
    double timers[5] = {0, 0, 0, 0, 0};
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
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
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);

    ierr = this->m_WorkVecField1->GetArrays(p_gmx1, p_gmx2, p_gmx3); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArrays(p_vtx1, p_vtx2, p_vtx3); CHKERRQ(ierr);

    // check if velocity field is zero
    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        // compute gradient of first time point of image component
        for (IntType k = 0; k < nc; ++k) {
            // template image is constant in time
            accfft_grad(p_gmx1, p_gmx2, p_gmx3, p_m+k*nl, this->m_Opt->GetFFT().plan, &XYZ, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            // compute incremental state variable for all time points
            for (IntType j = 0; j < nt; ++j) {
                l = j*nl*nc + k*nl;
                lnext = (j+1)*nl*nc + k*nl;

                // the right hand side remains constant;
                // we can reduce the 2 RK2 steps to a single one
                for (IntType i = 0; i < nl; ++i) {
                     p_mt[lnext+i] = p_mt[l+i] - ht*(p_gmx1[i]*p_vtx1[i]
                                                    +p_gmx2[i]*p_vtx2[i]
                                                    +p_gmx3[i]*p_vtx3[i]);
                }
            }
        }  // for all time points
    } else {  // velocity field is non-zero
        ierr = VecGetArray(this->m_WorkScaField1, &p_mtbar); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_WorkScaField2, &p_rhs0); CHKERRQ(ierr);

        ierr = this->m_WorkVecField2->GetArrays(p_gmtx1, p_gmtx2, p_gmtx3); CHKERRQ(ierr);
        ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j) {
            for (IntType k = 0; k < nc; ++k) {
                l = j*nl*nc + k*nl;
                lnext = (j+1)*nl*nc + k*nl;

                // compute gradient of m_j
                accfft_grad(p_gmx1, p_gmx2, p_gmx3, p_m+l, this->m_Opt->GetFFT().plan, &XYZ, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                // compute gradient of \tilde{m}_j
                accfft_grad(p_gmtx1, p_gmtx2, p_gmtx3, p_mt+l, this->m_Opt->GetFFT().plan, &XYZ, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                for (IntType i = 0; i < nl; ++i) {
                     p_rhs0[i] = -p_gmtx1[i]*p_vx1[i] - p_gmtx2[i]*p_vx2[i] - p_gmtx3[i]*p_vx3[i]
                                 -p_gmx1[i]*p_vtx1[i] - p_gmx2[i]*p_vtx2[i] - p_gmx3[i]*p_vtx3[i];
                     // compute intermediate result
                     p_mtbar[i] = p_mt[l+i] + ht*p_rhs0[i];
                }

                // compute gradient of m_{j+1}
                accfft_grad(p_gmx1, p_gmx2, p_gmx3, p_m+lnext, this->m_Opt->GetFFT().plan, &XYZ, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                // compute gradient of \tilde{m}_j
                accfft_grad(p_gmtx1, p_gmtx2, p_gmtx3, p_mtbar, this->m_Opt->GetFFT().plan, &XYZ, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                for (IntType i = 0; i < nl; ++i) {
                    // evaluate right hand side
                    ScalarType rhs1 = -p_gmtx1[i]*p_vx1[i] - p_gmtx2[i]*p_vx2[i] - p_gmtx3[i]*p_vx3[i]
                                      -p_gmx1[i]*p_vtx1[i] - p_gmx2[i]*p_vtx2[i] - p_gmx3[i]*p_vtx3[i];

                    // compute intermediate result
                    p_mt[lnext+i] = p_mt[l+i] + hthalf*(rhs1 + p_rhs0[i]);
                }
            }  // for all image components
        }  // for all time points

        // copy initial condition to buffer
        ierr = VecRestoreArray(this->m_WorkScaField1, &p_mtbar); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_WorkScaField2, &p_rhs0); CHKERRQ(ierr);

        ierr = this->m_WorkVecField2->RestoreArrays(p_gmtx1, p_gmtx2, p_gmtx3); CHKERRQ(ierr);
        ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    }  // velzero

    ierr = this->m_IncVelocityField->RestoreArrays(p_vtx1, p_vtx2, p_vtx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_gmx1, p_gmx2, p_gmx3); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timers);

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
PetscErrorCode OptimalControlRegistration::SolveIncStateEquationSL(void) {
    PetscErrorCode ierr = 0;
    IntType nl, nt, nc, l, lnext;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    ScalarType ht, hthalf;
    double timers[5] = {0, 0, 0, 0, 0};
    ScalarType *p_gm1 = NULL, *p_gm2 = NULL, *p_gm3 = NULL,
                *p_gmn1 = NULL, *p_gmn2 = NULL, *p_gmn3 = NULL,
                *p_gmX1 = NULL, *p_gmX2 = NULL, *p_gmX3 = NULL,
                *p_mtilde = NULL, *p_m = NULL;
    const ScalarType *p_vtildeX1 = NULL, *p_vtildeX2 = NULL, *p_vtildeX3 = NULL,
                     *p_vtilde1 = NULL, *p_vtilde2 = NULL, *p_vtilde3 = NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField2 == NULL) {
        try {this->m_WorkVecField2 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField3 == NULL) {
        try {this->m_WorkVecField3 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_WorkVecField4 == NULL) {
        try {this->m_WorkVecField4 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    }

    ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);

    // compute \tilde{\vect{v}}(X)
    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField2, this->m_IncVelocityField, "state"); CHKERRQ(ierr);

    ierr = this->m_IncVelocityField->GetArraysRead(p_vtilde1, p_vtilde2, p_vtilde3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_gmn1, p_gmn2, p_gmn3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->GetArraysRead(p_vtildeX1, p_vtildeX2, p_vtildeX3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->GetArrays(p_gm1, p_gm2, p_gm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->GetArrays(p_gmX1, p_gmX2, p_gmX3); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j) {  // for all time points
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            l = j*nl*nc + k*nl;
            lnext = (j+1)*nl*nc + k*nl;

            // compute gradient for m_j
            accfft_grad(p_gm1, p_gm2, p_gm3, p_m+l, this->m_Opt->GetFFT().plan, &XYZ, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            // compute gradient for m_j+1
            accfft_grad(p_gmn1, p_gmn2, p_gmn3, p_m+lnext, this->m_Opt->GetFFT().plan, &XYZ, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            // interpolate gradient of state variable
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_gmX1, p_gmX2, p_gmX3,
                                                             p_gm1, p_gm2, p_gm3, "state");

            // interpolate incremental adjoint variable
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_mtilde+lnext, p_mtilde+l, "state"); CHKERRQ(ierr);

            for (IntType i = 0; i < nl; ++i) {
                p_mtilde[lnext+i] -= hthalf*(p_gmX1[i]*p_vtildeX1[i]
                                           + p_gmX2[i]*p_vtildeX2[i]
                                           + p_gmX3[i]*p_vtildeX3[i]
                                           + p_gmn1[i]*p_vtilde1[i]
                                           + p_gmn2[i]*p_vtilde2[i]
                                           + p_gmn3[i]*p_vtilde3[i]);
            }
        }  // for all image components
    }  // for all time points

    ierr = this->m_WorkVecField1->RestoreArrays(p_gmn1, p_gmn2, p_gmn3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField2->RestoreArraysRead(p_vtildeX1, p_vtildeX2, p_vtildeX3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField3->RestoreArrays(p_gm1, p_gm2, p_gm3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField4->RestoreArrays(p_gmX1, p_gmX2, p_gmX3); CHKERRQ(ierr);

    ierr = this->m_IncVelocityField->RestoreArraysRead(p_vtilde1, p_vtilde2, p_vtilde3); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timers);

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
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquation(void) {
    PetscErrorCode ierr = 0;
    ScalarType *p_lt = NULL, *p_mt = NULL;
    IntType nl, ng, nc, nt, l;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_VelocityField != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ierr = Assert(nt > 0, "nt < 0"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2) {
        ss << "solving incremental adjoint equation (nx1,nx2,nx3,nc,nt) = ("
                  << this->m_Opt->GetDomainPara().nx[0]
           << "," << this->m_Opt->GetDomainPara().nx[1]
           << "," << this->m_Opt->GetDomainPara().nx[2]
           << "," << nc << "," << nt << ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    ierr = this->m_Opt->StartTimer(PDEEXEC); CHKERRQ(ierr);

    // allocate state and adjoint variables
    if (this->m_IncAdjointVariable == NULL) {
        ierr = VecCreate(this->m_IncAdjointVariable, (nt+1)*nc*nl, (nt+1)*nc*ng); CHKERRQ(ierr);
    }

    // set terminal condition \tilde{\lambda}_1 = -\tilde{m}_1
    ierr = VecGetArray(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);
    l = nt*nl*nc;
    for (IntType i = 0; i < nl*nc; ++i) {
        p_lt[l+i] = -p_mt[l+i]; // / static_cast<ScalarType>(nc);
    }
    ierr = VecRestoreArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_IncStateVariable, &p_mt); CHKERRQ(ierr);

    // check if velocity field is zero
    if (this->m_Opt->GetOptPara().method == GAUSSNEWTON) {   // gauss newton
        ierr = this->IsVelocityZero(); CHKERRQ(ierr);
        if (this->m_VelocityIsZero) {
            // copy terminal condition \tilde{\lambda}_1 = -\tilde{m}_1 to all time points
            ierr = VecGetArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);
            for (IntType j = 1; j <= nt; ++j) {
                try {std::copy(p_lt+nt*nc*nl, p_lt+(nt+1)*nc*nl, p_lt+(nt-j)*nl*nc);}
                catch (std::exception&) {
                    ierr = ThrowError("copy failed"); CHKERRQ(ierr);
                }
            }
            ierr = VecRestoreArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);
        } else {
            // call the solver
            switch (this->m_Opt->GetPDESolverPara().type) {
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
    } else if (this->m_Opt->GetOptPara().method == FULLNEWTON) {   // full newton
        // call the solver
        switch (this->m_Opt->GetPDESolverPara().type) {
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

    if (this->m_Opt->GetVerbosity() > 2) {
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
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquationGNRK2(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, l, lnext;
    ScalarType *p_lt = NULL, *p_rhs0 = NULL, *p_rhs1 = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_ltjvx1 = NULL, *p_ltjvx2 = NULL, *p_ltjvx3 = NULL;
    double timers[5] = {0, 0, 0, 0, 0};
    ScalarType ht, hthalf;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaField2 == NULL) {
        ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
    }

    ierr = VecGetArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_ltjvx1, p_ltjvx2, p_ltjvx3); CHKERRQ(ierr);
    ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    // compute numerical time integration
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            l = (nt-j)*nc*nl + k*nl;
            lnext = (nt-(j+1))*nc*nl + k*nl;

            // scale \vect{v} by \lambda
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                ScalarType lt = p_lt[l+i];

                p_ltjvx1[i] = p_vx1[i]*lt;
                p_ltjvx2[i] = p_vx2[i]*lt;
                p_ltjvx3[i] = p_vx3[i]*lt;
            }  // for all grid points

            // compute \idiv(\tilde{\lambda}\vect{v})
            accfft_divergence(p_rhs0, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->GetFFT().plan, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                // compute \bar{\tilde{\lambda}} = \tilde{\lambda}^j + ht*\idiv(\tilde{\lambda}^j\vect{v})
                ScalarType ltbar = p_lt[l+i] + ht*p_rhs0[i];

                // scale \vect{v} by \bar{\lambda}
                p_ltjvx1[i] = p_vx1[i]*ltbar;
                p_ltjvx2[i] = p_vx2[i]*ltbar;
                p_ltjvx3[i] = p_vx3[i]*ltbar;
            }

            // compute \idiv(\bar{\lambda}\vect{v})
            accfft_divergence(p_rhs1, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->GetFFT().plan, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                p_lt[lnext+i] = p_lt[l+i] + hthalf*(p_rhs0[i]+p_rhs1[i]);
            }
        }  // for all image components
    }  // for all time points

    ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_ltjvx1, p_ltjvx2, p_ltjvx3); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timers);

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
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquationFNRK2(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, l, lnext;
    ScalarType *p_l = NULL, *p_lt = NULL, *p_rhs0 = NULL, *p_rhs1 = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL,
                *p_vtx1 = NULL, *p_vtx2 = NULL, *p_vtx3 = NULL,
                *p_ltjvx1 = NULL, *p_ltjvx2 = NULL, *p_ltjvx3 = NULL;
    double timers[5] = {0, 0, 0, 0, 0};
    ScalarType ht, hthalf, lambda, lambdatilde, ltbar;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr = VecGetArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);

    ierr = VecRestoreArray(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_ltjvx1, p_ltjvx2, p_ltjvx3); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArrays(p_vtx1, p_vtx2, p_vtx3); CHKERRQ(ierr);

    ierr = this->IsVelocityZero(); CHKERRQ(ierr);
    if (this->m_VelocityIsZero) {
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            l = k*nl;
            // lambda and v are constant in time
            for (IntType i = 0; i < nl; ++i) {  // for all grid points
                ScalarType lambda = p_l[l+i];

                // scale \vect{v} by \lambda
                p_ltjvx1[i] = p_vtx1[i]*lambda;
                p_ltjvx2[i] = p_vtx2[i]*lambda;
                p_ltjvx3[i] = p_vtx3[i]*lambda;
            }  // for all grid points

            // compute \idiv(\tilde{\lambda}\vect{v})
            accfft_divergence(p_rhs0, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->GetFFT().plan, timers);
            this->m_Opt->IncrementCounter(FFT, 4);

            // compute numerical time integration
            for (IntType j = 0; j < nt; ++j) {  // for all time points
                l = (nt-j)*nc*nl + k*nl;
                lnext = (nt-(j+1))*nc*nl + k*nl;
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    p_lt[lnext+i] = p_lt[l+i] + ht*p_rhs0[i];
                }
            }  // for all time points
        }  // for all image components
    } else {  // velocity is zero
        if (this->m_WorkScaField2 == NULL) {
            ierr = VecCreate(this->m_WorkScaField2, nl, ng); CHKERRQ(ierr);
        }

        ierr = this->m_VelocityField->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);

        // compute numerical time integration
        for (IntType j = 0; j < nt; ++j) {  // for all time points
            for (IntType k = 0; k < nc; ++k) {  // for all image components
                l = (nt-j)*nc*nl + k*nl;
                lnext = (nt-(j+1))*nc*nl + k*nl;
                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    lambda  = p_l[l+i];
                    lambdatilde = p_lt[l+i];

                    p_ltjvx1[i] = p_vx1[i]*lambdatilde + p_vtx1[i]*lambda;
                    p_ltjvx2[i] = p_vx2[i]*lambdatilde + p_vtx2[i]*lambda;
                    p_ltjvx3[i] = p_vx3[i]*lambdatilde + p_vtx3[i]*lambda;
                }  // for all grid points

                // compute \idiv(\tilde{\lambda}\vect{v})
                accfft_divergence(p_rhs0, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->GetFFT().plan, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    // \bar{\lambda} = \tilde{\lambda}^j + ht*\idiv(\lambda^j\vect{v})
                    ltbar = p_lt[l+i] + ht*p_rhs0[i];
                    lambda = p_l[lnext+i];

                    // v \bar{\lambda} + \vect{\tilde{v}}\lambda^{j+1}
                    p_ltjvx1[i] = p_vx1[i]*ltbar + p_vtx1[i]*lambda;
                    p_ltjvx2[i] = p_vx2[i]*ltbar + p_vtx2[i]*lambda;
                    p_ltjvx3[i] = p_vx3[i]*ltbar + p_vtx3[i]*lambda;
                }

                // compute \idiv(\bar{\lambda}\vect{v})
                accfft_divergence(p_rhs1, p_ltjvx1, p_ltjvx2, p_ltjvx3, this->m_Opt->GetFFT().plan, timers);
                this->m_Opt->IncrementCounter(FFT, 4);

                for (IntType i = 0; i < nl; ++i) {  // for all grid points
                    p_lt[lnext+i] = p_lt[l+i] + hthalf*(p_rhs0[i]+p_rhs1[i]);
                }
            }  // for all image components
        }  // for all time points
        ierr = this->m_VelocityField->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_WorkScaField2, &p_rhs1); CHKERRQ(ierr);
    }  // velzero

    ierr = VecRestoreArray(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_IncAdjointVariable, &p_lt); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_rhs0); CHKERRQ(ierr);

    ierr = this->m_IncVelocityField->RestoreArrays(p_vtx1, p_vtx2, p_vtx3); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->RestoreArrays(p_ltjvx1, p_ltjvx2, p_ltjvx3); CHKERRQ(ierr);

    this->m_Opt->IncreaseFFTTimers(timers);
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
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquationGNSL(void) {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nc, nt, l, lnext;
    double timers[5] = {0, 0, 0, 0, 0};
    ScalarType *p_ltilde = NULL, *p_ltildejX = NULL,
                *p_divv = NULL, *p_divvX = NULL,
                *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL;
    ScalarType ht, hthalf, ltildejX, rhs0, rhs1;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;
    ht = this->m_Opt->GetTimeStepSize();
    hthalf = 0.5*ht;

    if (this->m_WorkVecField1 == NULL) {
        try {this->m_WorkVecField1 = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
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
    if (this->m_SemiLagrangianMethod == NULL) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangianType(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }

    // set v to -v
    ierr = this->m_WorkVecField1->Copy(this->m_VelocityField); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->Scale(-1.0); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);
    ierr = this->m_WorkVecField1->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    // compute div(v)
    accfft_divergence(p_divv, p_vx1, p_vx2, p_vx3, this->m_Opt->GetFFT().plan, timers);
    this->m_Opt->IncrementCounter(FFT, 4);

    ierr = this->m_WorkVecField1->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    // evaluate div(v) on characteristic
    ierr = VecGetArray(this->m_WorkScaField2, &p_divvX); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->Interpolate(p_divvX, p_divv, "adjoint"); CHKERRQ(ierr);

    ierr = VecGetArray(this->m_WorkScaField3, &p_ltildejX); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j) {
        for (IntType k = 0; k < nc; ++k) {
            l = (nt-j)*nc*nl + k*nl;
            lnext = (nt-(j+1))*nc*nl + k*nl;

            ierr = this->m_SemiLagrangianMethod->Interpolate(p_ltildejX, p_ltilde + l, "adjoint"); CHKERRQ(ierr);

            for (IntType i = 0; i < nl; ++i) {
                // get \tilde{\lambda}(X(t^j))
                ltildejX = p_ltildejX[i];

                // scale v(X) by \tilde{\lambda}(X(t^j))
                rhs0 = -ltildejX*p_divvX[i];

                // scale v by \lamba{\lambda}
                rhs1 = -(ltildejX + ht*rhs0)*p_divv[i];

                // final rk2 step
                p_ltilde[lnext+i] = ltildejX + hthalf*(rhs0 + rhs1);
            }
        }  // for all image components
    }  // for all time points

    ierr = VecRestoreArray(this->m_IncAdjointVariable, &p_ltilde); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField3, &p_ltildejX); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField2, &p_divvX); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_WorkScaField1, &p_divv); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(timers);
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
PetscErrorCode OptimalControlRegistration::SolveIncAdjointEquationFNSL(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = ThrowError("not implemented"); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief finalize the current iteration
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::FinalizeIteration(Vec v) {
    PetscErrorCode ierr = 0;
    int rank;
    IntType nl, ng, nc, nt, iter;
    std::string filename, ext;
    std::stringstream ss;
    std::ofstream logwriter;
    ScalarType *p_m1 = NULL, *p_m = NULL, rval, dval, jval, hd;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = Assert(v != NULL, "null pointer"); CHKERRQ(ierr);

    // get number of time points and grid points
    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    // parse extension
    ext = this->m_Opt->GetReadWriteFlags().extension;

    // if not yet allocted, do so
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    // set velocity field
    ierr = this->m_VelocityField->SetComponents(v); CHKERRQ(ierr);

    if (this->m_Opt->GetLogger().enabled[LOGCONV]) {
        iter = this->m_Opt->GetCounter(ITERATIONS);
        ierr = Assert(iter >= 0, "problem in counter"); CHKERRQ(ierr);

        ierr = this->EvaluateDistanceMeasure(&dval); CHKERRQ(ierr);
        ierr = this->m_Regularization->EvaluateFunctional(&rval, this->m_VelocityField); CHKERRQ(ierr);

        // get lebesque measure
        hd = this->m_Opt->GetLebesqueMeasure();

        // add up the contributions
        jval = hd*(dval + rval);
        this->m_Opt->LogConvergence(iter, jval, dval, rval);
    }


    // store iterates
    if (this->m_Opt->GetReadWriteFlags().iterates) {
        // allocate
        if (this->m_WorkScaFieldMC == NULL) {
            ierr = VecCreate(this->m_WorkScaFieldMC, nl*nc, ng*nc); CHKERRQ(ierr);
        }

        iter = this->m_Opt->GetCounter(ITERATIONS);
        ierr = Assert(iter >= 0, "problem in counter"); CHKERRQ(ierr);

        // copy memory for m_1
        ierr = VecGetArray(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

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


    if (this->m_Opt->StoreCheckPoints()) {
        ierr = this->m_ReadWrite->Write(this->m_VelocityField, "velocity-field-checkpoint"+ext); CHKERRQ(ierr);
    }


    // compute determinant of deformation gradient and write it to file
    if (this->m_Opt->GetRegMonitor().JAC) {
        ierr = this->ComputeDetDefGrad(); CHKERRQ(ierr);
        // if user enabled the logger
        if (this->m_Opt->GetLogger().enabled[LOGJAC]) {
            if (rank == 0) {
                filename  = this->m_Opt->GetReadWriteFlags().xfolder;
                filename += "registration-performance-detdefgrad.log";

                // create output file or append to output file
                logwriter.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
                ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);
                ss  << std::scientific
                    <<  "iter = "     << this->m_Opt->GetCounter(ITERATIONS)
                    <<  "   betav = " << this->m_Opt->GetRegNorm().beta[0] << "    "
                    << std::left << std::setw(20) << this->m_Opt->GetRegMonitor().jacmin << " "
                                 << std::setw(20) << this->m_Opt->GetRegMonitor().jacmean <<" "
                                 << std::setw(20) << this->m_Opt->GetRegMonitor().jacmax;
                logwriter << ss.str() << std::endl;
                ss.str(std::string()); ss.clear();
            }   // if on master rank
        }
    }

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief finalize the registration
 *******************************************************************/
PetscErrorCode OptimalControlRegistration::Finalize(VecField* v) {
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
    nt = this->m_Opt->GetDomainPara().nt;
    nc = this->m_Opt->GetDomainPara().nc;
    nl = this->m_Opt->GetDomainPara().nl;
    ng = this->m_Opt->GetDomainPara().ng;

    if (this->m_Opt->GetVerbosity() >= 2) {
        ierr = DbgMsg("finalizing registration"); CHKERRQ(ierr);
    }

    // if not yet allocted, do so and copy input
    if (this->m_VelocityField == NULL) {
        try {this->m_VelocityField = new VecField(this->m_Opt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    ierr = this->m_VelocityField->Copy(v); CHKERRQ(ierr);

    if (this->m_WorkScaField1 == NULL) {
        ierr = VecCreate(this->m_WorkScaField1, nl, ng); CHKERRQ(ierr);
    }
    if (this->m_WorkScaFieldMC == NULL) {
        ierr = VecCreate(this->m_WorkScaFieldMC, nl*nc, ng*nc); CHKERRQ(ierr);
    }

    // process timers
    ierr = this->m_Opt->ProcessTimers(); CHKERRQ(ierr);

    // parse extension
    ext = this->m_Opt->GetReadWriteFlags().extension;

    // compute residuals
    if (this->m_Opt->GetLogger().enabled[LOGRES]) {
        ierr = VecWAXPY(this->m_WorkScaFieldMC, -1.0, this->m_TemplateImage, this->m_ReferenceImage); CHKERRQ(ierr);

        ierr = VecNorm(this->m_WorkScaFieldMC, NORM_2, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(0, value);

        ierr = VecNorm(this->m_WorkScaFieldMC, NORM_INFINITY, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(1, value);

        // deformed template out (compute solution of state equation)
        ierr = this->SolveStateEquation(); CHKERRQ(ierr);

        // copy memory for m_1
        ierr = VecGetArray(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);

        // ||m_R - m_1||
        ierr = VecAXPY(this->m_WorkScaFieldMC, -1.0, this->m_ReferenceImage); CHKERRQ(ierr);

        ierr = VecNorm(this->m_WorkScaFieldMC, NORM_2, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(2, value);

        ierr = VecNorm(this->m_WorkScaFieldMC, NORM_INFINITY, &value); CHKERRQ(ierr);
        this->m_Opt->LogFinalResidual(3, value);
    }

    // write deformed template image to file
    if (this->m_Opt->GetReadWriteFlags().deftemplate) {
        // copy memory for m_1
        ierr = VecGetArray(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_m1);}
        catch (std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_WorkScaFieldMC, &p_m1); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->WriteT(this->m_WorkScaFieldMC, "deformed-template-image"+ext, nc > 1); CHKERRQ(ierr);
    }

    // write residual images to file
    if (this->m_Opt->GetReadWriteFlags().residual) {
        ierr = VecGetArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);

        ierr = VecGetArray(this->m_TemplateImage, &p_mt); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
        for (IntType i = 0; i < nl*nc; ++i) {
            p_dr[i] = 1.0 - PetscAbs(p_mr[i] - p_mt[i]);
        }
        ierr = VecRestoreArray(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
        ierr = VecRestoreArray(this->m_TemplateImage, &p_mt); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkScaFieldMC, "residual-t=0"+ext, nc > 1); CHKERRQ(ierr);

        // copy memory for m_1
        ierr = VecGetArray(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
        ierr = VecGetArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);
        try {std::copy(p_m+nt*nl*nc, p_m+(nt+1)*nl*nc, p_dr);}
        catch (std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_StateVariable, &p_m); CHKERRQ(ierr);

        for (IntType i = 0; i < nl*nc; ++i) {
            p_dr[i] = 1.0 - PetscAbs(p_mr[i] - p_dr[i]);
        }
        ierr = VecRestoreArray(this->m_WorkScaFieldMC, &p_dr); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkScaFieldMC, "residual-t=1"+ext, nc > 1); CHKERRQ(ierr);

        ierr = VecRestoreArray(this->m_ReferenceImage, &p_mr); CHKERRQ(ierr);
    }

    // write velocity field to file
    if (this->m_Opt->GetReadWriteFlags().results) {
        ierr = this->m_ReadWrite->Write(this->m_VelocityField, "velocity-field"+ext); CHKERRQ(ierr);
    }

    // write norm of velocity field to file
    if (this->m_Opt->GetReadWriteFlags().velnorm) {
        ierr = this->m_VelocityField->Norm(this->m_WorkScaField1); CHKERRQ(ierr);
        ierr = this->m_ReadWrite->Write(this->m_WorkScaField1, "velocity-field-norm"+ext); CHKERRQ(ierr);
    }

    // write determinant of deformation gradient to file
    if (this->m_Opt->GetReadWriteFlags().detdefgrad) {
        ierr = this->ComputeDetDefGrad(true); CHKERRQ(ierr);
    }

    // write determinant of deformation gradient to file
    if (this->m_Opt->GetReadWriteFlags().defgrad) {
        ierr = this->ComputeDefGrad(true); CHKERRQ(ierr);
    }

    // write deformation map to file
    if (this->m_Opt->GetReadWriteFlags().defmap) {
        ierr = this->ComputeDeformationMap(true); CHKERRQ(ierr);
    }

    // write deformation field to file
    if (this->m_Opt->GetReadWriteFlags().deffield) {
        ierr = this->ComputeDisplacementField(true); CHKERRQ(ierr);
    }

    // write template and reference image
    if (this->m_Opt->GetReadWriteFlags().templateim) {
        ierr = this->m_ReadWrite->WriteT(this->m_TemplateImage, "template-image"+ext, nc > 1); CHKERRQ(ierr);
    }
    if (this->m_Opt->GetReadWriteFlags().referenceim) {
        ierr = this->m_ReadWrite->WriteR(this->m_ReferenceImage, "reference-image"+ext, nc > 1); CHKERRQ(ierr);
    }

    // write log file
    ierr = this->m_Opt->WriteLogFile(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _OPTIMALCONTROLREGISTRATION_CPP_
