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

#ifndef _OPTIMIZATIONPROBLEM_HPP_
#define _OPTIMIZATIONPROBLEM_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"




namespace reg {
  
class Preprocessing;




class OptimizationProblem {
 public:
    typedef OptimizationProblem Self;
    OptimizationProblem(void);
    OptimizationProblem(RegOpt*);
    virtual ~OptimizationProblem(void);

    PetscErrorCode SetOptions(RegOpt*);
    inline RegOpt* GetOptions(void) {return this->m_Opt;};

    inline void SetConvergenceMessage(std::string msg){this->m_ConvergenceMessage = msg;}
    inline void Converged(bool flag){this->m_Converged = flag;}
    inline bool Converged(){return this->m_Converged;}
    inline std::string GetConvergenceMessage(){return this->m_ConvergenceMessage;}
    inline void IncrementIterations() {this->m_Opt->IncrementCounter(ITERATIONS);}

    /*! evaluate objective, gradient and distance measure for initial guess */
    virtual PetscErrorCode InitializeOptimization() = 0;
    
    /*! allocate all the memory we need */
    virtual PetscErrorCode InitializeSolver() = 0;

    /*! evaluate distance between observed and predicted state */
    virtual PetscErrorCode EvaluateDistanceMeasure(ScalarType*) = 0;

    /*! evaluate objective functional J(x) */
    virtual PetscErrorCode EvaluateObjective(ScalarType*, Vec) = 0;

    /*! evaluate gradient of Lagrangian L(x) */
    virtual PetscErrorCode EvaluateGradient(Vec, Vec) = 0;

    /*! apply Hessian matvec H\tilde{\vect{x}} */
    virtual PetscErrorCode HessianMatVec(Vec, Vec, bool scale = true) = 0;
    
    virtual PetscErrorCode PreHessian() = 0;

    /*! evaluate regularization functional for given control variable */
    virtual PetscErrorCode EvaluateRegularizationFunctional(ScalarType*, VecField*) = 0;

    /*! compute estimate of extremal eigenvalues of hessian */
    virtual PetscErrorCode EstimateExtremalHessEigVals(ScalarType&, ScalarType&) = 0;

    /*! pre processing before krylov solve */
    virtual PetscErrorCode PreKrylovSolve(Vec, Vec) = 0;

    /*! post processing after krylov solve */
    virtual PetscErrorCode PostKrylovSolve(Vec, Vec) = 0;

    /*! apply inverse regularization operator */
    virtual PetscErrorCode ApplyInvRegularizationOperator(Vec, Vec, bool flag = false) = 0;
    
    /*! apply inverse H(v=0) */
    virtual PetscErrorCode ApplyInvHessian(Vec, Vec, VecField*, bool first=false, bool twolevel=false, Preprocessing *preproc=nullptr) = 0;
    
    virtual PetscErrorCode SymTwoLevelHessMatVec(Vec, Vec) = 0;

    /*! solve forward problem */
    virtual PetscErrorCode SolveForwardProblem(Vec, Vec) = 0;

    /*! set mask (objective masking) */
    virtual PetscErrorCode SetMask(Vec) = 0;

    /*! get mask (objective masking) */
    virtual PetscErrorCode GetMask(Vec&) = 0;

    /*! set reference image */
    virtual PetscErrorCode SetReferenceImage(Vec) = 0;

    /*! get reference image */
    virtual PetscErrorCode GetReferenceImage(Vec&) = 0;

    /*! set control variable */
    virtual PetscErrorCode SetControlVariable(VecField*) = 0;

    /*! get control variable */
    virtual PetscErrorCode GetControlVariable(VecField*&) = 0;

    /*! get state variable */
    virtual PetscErrorCode GetStateVariable(Vec&) = 0;

    /*! set state variable */
    virtual PetscErrorCode SetStateVariable(Vec) = 0;

    /*! get adjoint variable */
    virtual PetscErrorCode GetAdjointVariable(Vec&) = 0;

    /*! set adjoint variable */
    virtual PetscErrorCode SetAdjointVariable(Vec) = 0;

    /*! finalize iteration */
    virtual PetscErrorCode FinalizeIteration(Vec) = 0;

    /*! finalize registration */
    virtual PetscErrorCode Finalize(VecField*) = 0;

    /*! apply two level preconditioner */
    virtual PetscErrorCode CheckBounds(Vec, bool&) = 0;

    /*! check gradient (derivative check via taylor expansion) */
    PetscErrorCode DerivativeCheckGradient(void);
    PetscErrorCode DerivativeCheckHessian(void);
    PetscErrorCode DerivativeCheckHessianFD(void);

    /*! check if hessian is symmetric */
    PetscErrorCode HessianSymmetryCheck(void);

    /*! check if hessian is symmetric */
    PetscErrorCode ComputeUpdateNorm(Vec, ScalarType&, ScalarType&);
    
    VecField* m_WorkVecField1;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField2;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField3;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField4;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField5;  ///< data container for vector field (temporary variable)

 protected:
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);

    RegOpt* m_Opt;

 private:
    Vec m_Iterate;
    std::string m_ConvergenceMessage;
    bool m_Converged;
};




}  // end of namespace




#endif  // _OPTIMIZATIONPROBLEM_HPP_
