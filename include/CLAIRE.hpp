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

#ifndef _CLAIRE_HPP_
#define _CLAIRE_HPP_

#include <fstream>

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "CLAIREBase.hpp"




namespace reg {




class CLAIRE : public CLAIREBase {
 public:
    typedef CLAIREBase SuperClass;
    typedef CLAIRE Self;

    CLAIRE(void);
    CLAIRE(RegOpt*);
    virtual ~CLAIRE(void);

    /*! evaluate objective, gradient and distance measure for initial guess */
    PetscErrorCode InitializeOptimization();

    /*! evaluate distance between observed and predicted state */
    PetscErrorCode EvaluateDistanceMeasure(ScalarType*);

    /*! evaluate objective value */
    PetscErrorCode EvaluateObjective(ScalarType*, Vec);

    /*! evaluate reduced gradient (first variation
        of lagrangian with respect to control variable(s) */
    PetscErrorCode EvaluateGradient(Vec, Vec);

    /*! compute Hessian matvec (second variation
        of lagrangian with respect to control variable(s) */
    PetscErrorCode HessianMatVec(Vec, Vec, bool scale = true);
    
    /*! apply inverse H(v=0) */
    PetscErrorCode ApplyInvHessian(Vec, Vec, VecField*, bool first=false, bool twolevel=false, Preprocessing *preproc=nullptr);

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
    
    /*! get the final state variable */
    PetscErrorCode GetFinalState(Vec);

    /*! solve state equation */
    virtual PetscErrorCode SolveStateEquation(void);

    PetscErrorCode SymTwoLevelHessMatVec(Vec, Vec);
 protected:
    /*! compute Hessian matvec (second variation
        of lagrangian with respect to control variable(s) for zero velocity */
    PetscErrorCode H0HessMatVec(Vec, Vec);
    
    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! delete / reset variables */
    PetscErrorCode ClearVariables(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);

    /*! get the final state variable */
    PetscErrorCode SetInitialState(Vec);

    /*! get the final state variable */
    PetscErrorCode SetFinalAdjoint(Vec);

    /*! solve incremental state equation */
    virtual PetscErrorCode SolveIncStateEquation(void);

    /*! solve adjoint equation */
    virtual PetscErrorCode SolveAdjointEquation(void);

    /*! solve incremental adjoint equation */
    virtual PetscErrorCode SolveIncAdjointEquation(void);

    /*! evaluate l2-gradient */
    virtual PetscErrorCode EvaluateL2Gradient(Vec);

    /*! evaluate sobolev gradient */
    virtual PetscErrorCode EvaluateSobolevGradient(Vec, bool flag = false);

    /*! sl solver for continuity equation */
    PetscErrorCode SolveContinuityEquationSL();

    /*! apply the projection operator to the
        body force and the incremental body force */
    virtual PetscErrorCode ApplyProjection();

    //ScaField* m_StateVariable;        ///< time dependent state variable m(x,t)
    //ScaField* m_AdjointVariable;      ///< time dependent adjoint variable \lambda(x,t)
    //ScaField* m_IncStateVariable;     ///< time dependent incremental state variable \tilde{m}(x,t)
    //ScaField* m_IncAdjointVariable;   ///< time dependent incremental adjoint variable \tilde{\lambda}(x,t)
    
    virtual PetscErrorCode CreateCoarseReg();
    virtual PetscErrorCode InitializeCoarseReg();
    
    //CLAIRE* m_CoarseReg;
    //RegOpt* m_CoarseRegOpt;

 private:
    /*! compute the initial guess for the velocity field */
    PetscErrorCode ComputeInitialVelocity(void);

    PetscErrorCode HessMatVec(Vec, Vec);
    PetscErrorCode PrecondHessMatVec(Vec, Vec);
    PetscErrorCode PrecondHessMatVecSym(Vec, Vec);

    PetscErrorCode StoreStateVariable();
};




}  // namespace reg




#endif  // _CLAIRE_H_
