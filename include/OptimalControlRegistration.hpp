/*
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
 *  along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 *
 */



#ifndef _OPTIMALCONTROLREGISTRATION_H_
#define _OPTIMALCONTROLREGISTRATION_H_

#include <fstream>

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "RegularizationRegistration.hpp"
#include "RegularizationRegistrationH1.hpp"
#include "RegularizationRegistrationH2.hpp"
#include "RegularizationRegistrationH1SN.hpp"
#include "RegularizationRegistrationH2SN.hpp"
#include "OptimalControlRegistrationBase.hpp"

namespace reg
{


class OptimalControlRegistration : public OptimalControlRegistrationBase
{
public:
    typedef OptimalControlRegistrationBase SuperClass;
    typedef OptimalControlRegistration Self;

    OptimalControlRegistration(void);
    OptimalControlRegistration(RegOpt*);
    ~OptimalControlRegistration(void);

    /*! evaluate objective, gradient and distance measure for initial guess */
    PetscErrorCode InitializeOptimization();

    /*! evaluate distance between observed and predicted state */
    PetscErrorCode EvaluateDistanceMeasure(ScalarType*);

    /*! evaluate objective value */
    PetscErrorCode EvaluateObjective(ScalarType*,Vec);

    /*! evaluate reduced gradient (first variation
        of lagrangian with respect to control variable(s) */
    PetscErrorCode EvaluateGradient(Vec,Vec);

    /*! compute Hessian matvec (second variation
        of lagrangian with respect to control variable(s) */
    PetscErrorCode HessianMatVec(Vec,Vec);

    /*! apply preconditioner for KKT system */
    PetscErrorCode PrecondMatVec(Vec, Vec);

    /*! solve the state equation */
    PetscErrorCode SolveForwardProblem(Vec);

    /*! finalize iteration */
    PetscErrorCode FinalizeIteration(Vec);

    /*! finalize registration */
    PetscErrorCode Finalize(VecField*);

    /*! piccard iteration */
    PetscErrorCode PiccardIteration(Vec);

    /*! compute initial condition via piccard iteration */
    PetscErrorCode ComputeInitialCondition(Vec,Vec);

protected:

    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

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

    // RK2 solvers
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

    /*! setup two level preconditioner */
    PetscErrorCode AllocateRegularization();

    Vec m_StateVariable; ///< time dependent state variable m(x,t)
    Vec m_AdjointVariable; ///< time dependent adjoint variable \lambda(x,t)
    Vec m_IncStateVariable; ///< time dependent incremental state variable \tilde{m}(x,t)
    Vec m_IncAdjointVariable; ///< time dependent incremental adjoint variable \tilde{\lambda}(x,t)

private:


};

} // end of namespace


#endif
