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

#ifndef _CLAIREBASE_HPP_
#define _CLAIREBASE_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"
#include "TenField.hpp"
#include "ScaField.hpp"
#include "Preprocessing.hpp"
#include "ReadWriteReg.hpp"
#include "DeformationFields.hpp"
#include "DistanceMeasure.hpp"
#include "DistanceMeasureSL2.hpp"
#include "DistanceMeasureSL2aux.hpp"
#include "DistanceMeasureNCC.hpp"
#include "Regularization.hpp"
#include "Regularization.hpp"
#include "RegularizationL2.hpp"
#include "RegularizationH1.hpp"
#include "RegularizationH2.hpp"
#include "RegularizationH3.hpp"
#include "RegularizationH1SN.hpp"
#include "RegularizationH2SN.hpp"
#include "RegularizationH3SN.hpp"
#include "OptimizationProblem.hpp"
#include "SemiLagrangian.hpp"
#include "Differentiation.hpp"
#include "DifferentiationFD.hpp"
#include "DifferentiationSM.hpp"
#include "TransportProblem.hpp"
#include "TransportEquationSL.hpp"
#include "TransportEquationRK2.hpp"




namespace reg {




class CLAIREBase : public OptimizationProblem {
 public:
    typedef CLAIREBase Self;
    typedef OptimizationProblem SuperClass;
    typedef Regularization RegularizationType;
    typedef SemiLagrangian SemiLagrangianType;

    CLAIREBase(void);
    CLAIREBase(RegOpt*);
    virtual ~CLAIREBase(void);

    /*! set io object */
    PetscErrorCode SetReadWrite(ReadWriteReg*);

    /*! set template image */
    PetscErrorCode SetTemplateImage(Vec);

    /*! set reference image */
    PetscErrorCode SetReferenceImage(Vec);

    /*! set mask (mask objective) */
    PetscErrorCode SetMask(Vec);

    /*! set mask (mask objective) */
    PetscErrorCode GetMask(Vec&);

    /*! set cell density c (coupled formulation) */
    PetscErrorCode SetCellDensity(Vec);

    /*! set auxilary variable (q; coupled formulation)*/
    PetscErrorCode SetAuxVariable(Vec);

    /*! get template image */
    PetscErrorCode GetTemplateImage(Vec&);

    /*! get reference image */
    PetscErrorCode GetReferenceImage(Vec&);

    /*! set control variable */
    PetscErrorCode SetControlVariable(VecField*);

    /*! set control variable */
    PetscErrorCode SetIncControlVariable(VecField*);

    /*! get control variable */
    PetscErrorCode GetControlVariable(VecField*&);

    /*! get state variable */
    virtual PetscErrorCode GetStateVariable(Vec&) = 0;

    /*! set state variable */
    virtual PetscErrorCode SetStateVariable(Vec) = 0;

    /*! get adjoint variable */
    virtual PetscErrorCode GetAdjointVariable(Vec&) = 0;

    /*! set state variable */
    virtual PetscErrorCode SetAdjointVariable(Vec) = 0;
    
    virtual PetscErrorCode PreHessian();

    /*! compute determinant of deformation gradient, i.e.
        the jacobian of the deformation map */
    PetscErrorCode ComputeDetDefGrad(bool write2file = false, Vec detj = NULL);

    /*! compute deformation map */
    PetscErrorCode ComputeDeformationMaps(bool write2file = false);

    /*! compute synthetic test problem */
    PetscErrorCode SetupSyntheticProb(Vec&, Vec&, VecField(*)=nullptr);

    /*! evaluate objective, gradient and distance measure for initial guess */
    virtual PetscErrorCode InitializeOptimization() = 0;

    /*! evaluate l2-distance between observed and predicted state */
    virtual PetscErrorCode EvaluateDistanceMeasure(ScalarType*) = 0;

    /*! evaluate objective functional J(v) */
    virtual PetscErrorCode EvaluateObjective(ScalarType*, Vec) = 0;

    /*! evaluate regularization functional R(v) */
    PetscErrorCode EvaluateRegularizationFunctional(ScalarType*, VecField*);

    /*! evaluate reduced gradient */
    virtual PetscErrorCode EvaluateGradient(Vec, Vec) = 0;

    /*! apply Hessian matvec H\tilde{\vect{v}} */
    virtual PetscErrorCode HessianMatVec(Vec, Vec, bool scale = true) = 0;

    /*! compute estimate of extremal eigenvalues of hessian */
    PetscErrorCode EstimateExtremalHessEigVals(ScalarType&, ScalarType&);

    /*! pre processing before krylov solve */
    PetscErrorCode PreKrylovSolve(Vec, Vec);

    /*! post processing after krylov solve */
    PetscErrorCode PostKrylovSolve(Vec, Vec);

    /*! apply inverse regularization operator */
    PetscErrorCode ApplyInvRegularizationOperator(Vec, Vec, bool flag = false);
    
    /*! apply inverse H(v=0) */
    virtual PetscErrorCode ApplyInvHessian(Vec, Vec, VecField*, bool first=false, bool twolevel=false, Preprocessing *preproc=nullptr) = 0;
    
    virtual PetscErrorCode SymTwoLevelHessMatVec(Vec, Vec) = 0;


    /*! solve forward problem */
    virtual PetscErrorCode SolveForwardProblem(Vec, Vec) = 0;

    /*! solve forward problem */
    virtual PetscErrorCode SolveAdjointProblem(Vec, Vec) = 0;

    /*! solve the current iteration */
    virtual PetscErrorCode FinalizeIteration(Vec) = 0;

    /*! finalize the registration */
    virtual PetscErrorCode Finalize(VecField*) = 0;

    /*! function that checks bounds in parameter continuation */
    virtual PetscErrorCode CheckBounds(Vec, bool&);

    /*! allocate all the memory we need */
    virtual PetscErrorCode InitializeSolver() = 0;
    
    /*! get the final state variable */
    virtual PetscErrorCode GetFinalState(Vec) = 0;

    /*! solve state equation */
    virtual PetscErrorCode SolveStateEquation(void) = 0;

 protected:
    PetscErrorCode Initialize(void);
    PetscErrorCode ComputeInitialGuess(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode CopyToAllTimePoints(Vec, Vec);
    PetscErrorCode IsVelocityZero(void);

    virtual PetscErrorCode ClearVariables(void) = 0;

    /*! evaluate l2-gradient */
    virtual PetscErrorCode EvaluateL2Gradient(Vec) = 0;

    /*! evaluate sobolev gradient */
    virtual PetscErrorCode EvaluateSobolevGradient(Vec, bool flag = false) = 0;

    /*! allocate regularization operator */
    PetscErrorCode SetupRegularization();
    PetscErrorCode SetupDeformationField();
    PetscErrorCode SetupSpectralData();
    PetscErrorCode SetupDistanceMeasure();
    PetscErrorCode SetupTransportProblem();

    /*! compute cfl condition */
    PetscErrorCode ComputeCFLCondition();

    ScaField* m_TemplateImage;           ///< data container for reference image mR
    ScaField* m_ReferenceImage;          ///< data container for template image mT
    ScaField* m_AuxVariable;             ///< auxilary variable
    ScaField* m_CellDensity;             ///< cell density
    ScaField* m_Mask;                    ///< mask for objective functional masking

    VecField* m_VelocityField;      ///< data container for velocity field (control variable)
    VecField* m_IncVelocityField;   ///< data container for incremental velocity field (incremental control variable)
    
    VecField** m_GradientState;   ///< Gradients of the state variable
    VecField** m_GradientXState;   ///< Gradients of the interpolated state variable

    ScaField* m_WorkScaField1;  ///< work scalar field
    ScaField* m_WorkScaField2;  ///< work scalar field
    ScaField* m_WorkScaField3;  ///< work scalar field
    ScaField* m_WorkScaField4;  ///< work scalar field
    ScaField* m_WorkScaField5;  ///< work scalar field

    ScaField* m_WorkScaFieldMC;  ///< work scalar field for multi-component/vector fields

    ReadWriteReg* m_ReadWrite;                   ///< io; set from outside (not to be delted)
    RegularizationType* m_Regularization;        ///< regularization functional
    DistanceMeasure* m_DistanceMeasure;          ///< distance measure
    SemiLagrangianType* m_SemiLagrangianMethod;  ///< semi-lagrangian method
    DeformationFields* m_DeformationFields;      ///< interface to compute deformation fields from velocity
    Differentiation* m_Differentiation;          ///< interface to evaluate / apply differential operators
    DifferentiationFD* m_DifferentiationFD;          ///< interface to evaluate / apply differential operators
    TransportProblem* m_TransportProblem;        ///< interface to compute the transport equation

    bool m_VelocityIsZero;
    bool m_StoreTimeHistory;

    ComplexType *m_x1hat;
    ComplexType *m_x2hat;
    ComplexType *m_x3hat;
    
    ScaField* m_StateVariable;        ///< time dependent state variable m(x,t)
    ScaField* m_AdjointVariable;      ///< time dependent adjoint variable \lambda(x,t)
    ScaField* m_IncStateVariable;     ///< time dependent incremental state variable \tilde{m}(x,t)
    ScaField* m_IncAdjointVariable;   ///< time dependent incremental adjoint variable \tilde{\lambda}(x,t)
    
    virtual PetscErrorCode CreateCoarseReg() = 0;
    virtual PetscErrorCode InitializeCoarseReg() = 0;
    
    CLAIREBase* m_CoarseReg;
    RegOpt* m_CoarseRegOpt;
    
        /*! solve incremental state equation */
    virtual PetscErrorCode SolveIncStateEquation(void) = 0;

    /*! solve incremental adjoint equation */
    virtual PetscErrorCode SolveIncAdjointEquation(void) = 0;

 private:
    bool m_DeleteControlVariable;
    bool m_DeleteIncControlVariable;
    
  friend class CLAIRE;
};




}  // end of namespace




#endif  // _CLAIREBASE_HPP_
