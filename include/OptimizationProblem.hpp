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

#ifndef _OPTIMIZATIONPROBLEM_H_
#define _OPTIMIZATIONPROBLEM_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "VecField.hpp"




namespace reg {




class OptimizationProblem {
 public:
    typedef OptimizationProblem Self;
    OptimizationProblem(void);
    OptimizationProblem(RegOpt*);
    ~OptimizationProblem(void);

    PetscErrorCode SetOptions(RegOpt*);
    inline RegOpt* GetOptions(void) {return this->m_Opt;};

    inline ScalarType GetInitialObjVal() {return this->m_InitObjectiveVal;};
    inline ScalarType GetInitialDistanceVal() {return this->m_InitDistanceVal;};
    inline ScalarType GetInitialGradNorm() {return this->m_InitGradNorm;};

    inline void IncrementIterations() {this->m_Opt->IncrementCounter(ITERATIONS);};

    /*! evaluate objective, gradient and distance measure for initial guess */
    virtual PetscErrorCode InitializeOptimization(VecField* v0 = NULL) = 0;

    /*! evaluate distance between observed and predicted state */
    virtual PetscErrorCode EvaluateDistanceMeasure(ScalarType*) = 0;

    /*! evaluate objective functional J(x) */
    virtual PetscErrorCode EvaluateObjective(ScalarType*, Vec) = 0;

    /*! evaluate gradient of Lagrangian L(x) */
    virtual PetscErrorCode EvaluateGradient(Vec,Vec) = 0;

    /*! apply Hessian matvec H\tilde{\vect{x}} */
    virtual PetscErrorCode HessianMatVec(Vec,Vec,bool scale=true) = 0;

    /*! compute estimate of extremal eigenvalues of hessian */
    virtual PetscErrorCode EstimateExtremalHessEigVals(ScalarType&,ScalarType&) = 0;

    /*! pre processing before krylov solve */
    virtual PetscErrorCode PreKrylovSolve(Vec,Vec) = 0;

    /*! post processing after krylov solve */
    virtual PetscErrorCode PostKrylovSolve(Vec,Vec) = 0;

    /*! apply inverse regularization operator */
    virtual PetscErrorCode ApplyInvRegOp(Vec,Vec) = 0;

    /*! solve forward problem */
    virtual PetscErrorCode SolveForwardProblem(Vec,Vec) = 0;

    /*! set control variable */
    virtual PetscErrorCode SetControlVariable(VecField*) = 0;

    /*! get control variable */
    virtual PetscErrorCode GetControlVariable(VecField*&) = 0;

    /*! get state variable */
    virtual PetscErrorCode GetStateVariable(Vec&) = 0;

    /*! set state variable */
    virtual PetscErrorCode SetStateVariable(Vec) = 0;

    /*! get state variable */
    virtual PetscErrorCode GetAdjointVariable(Vec&) = 0;

    /*! set state variable */
    virtual PetscErrorCode SetAdjointVariable(Vec) = 0;

    /*! finalize iteration */
    virtual PetscErrorCode FinalizeIteration(Vec) = 0;

    /*! finalize registration */
    virtual PetscErrorCode Finalize(VecField*) = 0;

    /*! apply two level preconditioner */
    virtual PetscErrorCode CheckBounds(Vec,bool&) = 0;

    /*! check gradient (derivative check via tayler expansion) */
    PetscErrorCode DerivativeCheck(void);

    /*! check if hessian is symmetric */
    PetscErrorCode HessianSymmetryCheck(void);

protected:

    PetscErrorCode Initialize(void);
    virtual PetscErrorCode ClearMemory(void) = 0;

    RegOpt* m_Opt;

    ScalarType m_InitGradNorm;
    ScalarType m_InitObjectiveVal;
    ScalarType m_InitDistanceVal;

private:

};

} // end of namespace


#endif // _OPTIMIZATIONPROBLEM_H_

