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

#ifndef _OPTIMALCONTROLREGISTRATION_H_
#define _OPTIMALCONTROLREGISTRATION_H_

#include <fstream>

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "OptimalControlRegistrationBase.hpp"




namespace reg {




class OptimalControlRegistration : public OptimalControlRegistrationBase {
 public:
    typedef OptimalControlRegistrationBase SuperClass;
    typedef OptimalControlRegistration Self;

    OptimalControlRegistration(void);
    OptimalControlRegistration(RegOpt*);
    virtual ~OptimalControlRegistration(void);

    /*! evaluate objective, gradient and distance measure for initial guess */
    PetscErrorCode InitializeOptimization(VecField* v0 = NULL);

    /*! evaluate distance between observed and predicted state */
    PetscErrorCode EvaluateDistanceMeasure(ScalarType*);

    /*! evaluate objective value */
    PetscErrorCode EvaluateObjective(ScalarType*, Vec);

    /*! evaluate reduced gradient (first variation
        of lagrangian with respect to control variable(s) */
    PetscErrorCode EvaluateGradient(Vec, Vec);

    /*! evaluate preconditined reduced gradient of Lagrangian L(v) */
    PetscErrorCode EvaluatePrecondGradient(Vec, Vec);

    /*! compute Hessian matvec (second variation
        of lagrangian with respect to control variable(s) */
    PetscErrorCode HessianMatVec(Vec, Vec, bool scale = true);

    /*! get state variable */
    PetscErrorCode GetStateVariable(Vec&);

    /*! set state variable */
    PetscErrorCode SetStateVariable(Vec);

    /*! get state variable */
    PetscErrorCode GetAdjointVariable(Vec&);

    /*! set state variable */
    PetscErrorCode SetAdjointVariable(Vec);

    /*! solve the state equation */
    PetscErrorCode SolveForwardProblem(Vec, Vec);

    /*! solve the state equation */
    PetscErrorCode SolveAdjointProblem(Vec, Vec);

    /*! finalize iteration */
    PetscErrorCode FinalizeIteration(Vec);

    /*! finalize registration */
    PetscErrorCode Finalize(VecField*);

    /*! compute initial condition via piccard iteration */
    PetscErrorCode ComputeInitialCondition(Vec, Vec);

    /*! allocate all the memory we need */
    PetscErrorCode InitializeSolver();

 protected:
    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! delete / reset variables */
    PetscErrorCode ClearVariables(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);

    /*! solve state equation */
    virtual PetscErrorCode SolveStateEquation(void);

    /*! solve incremental state equation */
    virtual PetscErrorCode SolveIncStateEquation(void);

    /*! solve adjoint equation */
    virtual PetscErrorCode SolveAdjointEquation(void);

    /*! solve incremental adjoint equation */
    virtual PetscErrorCode SolveIncAdjointEquation(void);

    /*! compute body force */
    virtual PetscErrorCode ComputeBodyForce(void);

    /*! compute incremental body force */
    virtual PetscErrorCode ComputeIncBodyForce(void);

    /*! rk2 solver for state equation */
    PetscErrorCode SolveStateEquationRK2();

    /*! rk2 solver for adjoint equation */
    PetscErrorCode SolveAdjointEquationRK2();

    /*! rk2 solver for incremental state equation */
    PetscErrorCode SolveIncStateEquationRK2();

    /*! rk2 solver for incremental adjoint equation (Newton step) */
    PetscErrorCode SolveIncAdjointEquationFNRK2();

    /*! rk2 solver for incremental adjoint equation (Gauss--Newton approximation) */
    PetscErrorCode SolveIncAdjointEquationGNRK2();

    /*! sl solver for state equation */
    PetscErrorCode SolveStateEquationSL();

    /*! sl solver for adjoint equation */
    virtual PetscErrorCode SolveAdjointEquationSL();

    /*! sl solver for inc state equation */
    PetscErrorCode SolveIncStateEquationSL();

    /*! sl solver for inc adjoint equation */
    virtual PetscErrorCode SolveIncAdjointEquationGNSL();

    /*! sl solver for inc adjoint equation */
    PetscErrorCode SolveIncAdjointEquationFNSL();

    /*! apply the projection operator to the
        body force and the incremental body force */
    virtual PetscErrorCode ApplyProjection();

    Vec m_StateVariable;        ///< time dependent state variable m(x,t)
    Vec m_AdjointVariable;      ///< time dependent adjoint variable \lambda(x,t)
    Vec m_IncStateVariable;     ///< time dependent incremental state variable \tilde{m}(x,t)
    Vec m_IncAdjointVariable;   ///< time dependent incremental adjoint variable \tilde{\lambda}(x,t)

 private:
    /*! compute the initial guess for the velocity field */
    PetscErrorCode ComputeInitialVelocity(void);

    PetscErrorCode HessMatVec(Vec, Vec);
    PetscErrorCode PrecondHessMatVec(Vec, Vec);
    PetscErrorCode PrecondHessMatVecSym(Vec, Vec);
};




}  // namespace reg




#endif  // _OPTIMALCONTROLREGISTRATION_H_
