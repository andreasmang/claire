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

#ifndef _OPTIMALCONTROLREGISTRATIONBASE_H_
#define _OPTIMALCONTROLREGISTRATIONBASE_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "VecField.hpp"
#include "TenField.hpp"
#include "Preprocessing.hpp"
#include "ReadWriteReg.hpp"
#include "SemiLagrangian.hpp"
#include "RegularizationRegistration.hpp"
#include "RegularizationRegistration.hpp"
#include "RegularizationRegistrationL2.hpp"
#include "RegularizationRegistrationH1.hpp"
#include "RegularizationRegistrationH2.hpp"
#include "RegularizationRegistrationH3.hpp"
#include "RegularizationRegistrationH1SN.hpp"
#include "RegularizationRegistrationH2SN.hpp"
#include "RegularizationRegistrationH3SN.hpp"
#include "OptimizationProblem.hpp"
//#include "SemiLagrangianGPU.hpp"




namespace reg {




class OptimalControlRegistrationBase : public OptimizationProblem {
 public:
    typedef ReadWriteReg ReadWriteType;
    typedef OptimalControlRegistrationBase Self;
    typedef OptimizationProblem SuperClass;
    typedef RegularizationRegistration RegularizationType;
    typedef SemiLagrangian SemiLagrangianType;

    OptimalControlRegistrationBase(void);
    OptimalControlRegistrationBase(RegOpt*);
    ~OptimalControlRegistrationBase(void);

    /*! set io object */
    PetscErrorCode SetReadWrite(ReadWriteType*);

    /*! set template image */
    PetscErrorCode SetTemplateImage(Vec);

    /*! set reference image */
    PetscErrorCode SetReferenceImage(Vec);

    /*! set template image */
    PetscErrorCode GetTemplateImage(Vec&);

    /*! set reference image */
    PetscErrorCode GetReferenceImage(Vec&);

    /*! set control variable */
    PetscErrorCode SetControlVariable(VecField*);

    /*! get control variable */
    PetscErrorCode GetControlVariable(VecField*&);

    /*! get state variable */
    virtual PetscErrorCode GetStateVariable(Vec&) = 0;

    /*! set state variable */
    virtual PetscErrorCode SetStateVariable(Vec) = 0;

    /*! get state variable */
    virtual PetscErrorCode GetAdjointVariable(Vec&) = 0;

    /*! set state variable */
    virtual PetscErrorCode SetAdjointVariable(Vec) = 0;

    /*! compute deformation gradient, i.e. jacobian of deformation map */
    PetscErrorCode ComputeDefGrad(bool write2file = false);

    /*! compute determinant of deformation gradient, i.e.
        the jacobian of the deformation map */
    PetscErrorCode ComputeDetDefGrad(bool write2file = false, Vec detj = NULL);

    /*! compute deformation map */
    PetscErrorCode ComputeDeformationMap(bool write2file = false, VecField* y = NULL);

    /*! check deformation map computation */
    PetscErrorCode CheckDefMapConsistency();

    /*! compute displacement field */
    PetscErrorCode ComputeDisplacementField(bool write2file = false);

    /*! compute synthetic test problem */
    PetscErrorCode SetupSyntheticProb(Vec&, Vec&);

    /*! evaluate objective, gradient and distance measure for initial guess */
    virtual PetscErrorCode InitializeOptimization(VecField* v0 = NULL) = 0;

    /*! evaluate l2-distance between observed and predicted state */
    virtual PetscErrorCode EvaluateDistanceMeasure(ScalarType*) = 0;

    /*! evaluate objective functional J(v) */
    virtual PetscErrorCode EvaluateObjective(ScalarType*, Vec) = 0;

    /*! evaluate regularization functional R(v) */
    PetscErrorCode EvaluateRegularizationFunctional(ScalarType*, VecField*);

    /*! evaluate gradient of Lagrangian L(v) */
    virtual PetscErrorCode EvaluateGradient(Vec, Vec) = 0;

    /*! evaluate preconditined reduced gradient of Lagrangian L(v) */
    virtual PetscErrorCode EvaluatePrecondGradient(Vec, Vec) = 0;

    /*! apply Hessian matvec H\tilde{\vect{v}} */
    virtual PetscErrorCode HessianMatVec(Vec, Vec, bool scale = true) = 0;

    /*! compute estimate of extremal eigenvalues of hessian */
    PetscErrorCode EstimateExtremalHessEigVals(ScalarType&, ScalarType&);

    /*! pre processing before krylov solve */
    PetscErrorCode PreKrylovSolve(Vec, Vec);

    /*! post processing after krylov solve */
    PetscErrorCode PostKrylovSolve(Vec, Vec);

    /*! apply inverse regularization operator */
    PetscErrorCode ApplyInvRegOp(Vec, Vec);

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


 protected:
    PetscErrorCode Initialize(void);
    PetscErrorCode ComputeInitialGuess(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode CopyToAllTimePoints(Vec, Vec);
    PetscErrorCode IsVelocityZero(void);

    virtual PetscErrorCode ClearVariables(void) = 0;

    /*! allocate regularization operator */
    PetscErrorCode AllocateRegularization();
    PetscErrorCode ComputeDefMapFromDisplacement();  ///< compute deformation map from displacement
    PetscErrorCode ApplyInvRegOpSqrt(Vec);
    PetscErrorCode ComputeRegularGrid(VecField*);  ///< compute coordinates for regular grid

    /*! compute cfl condition */
    PetscErrorCode ComputeCFLCondition();

    Vec m_TemplateImage;            ///< data container for reference image mR
    Vec m_ReferenceImage;           ///< data container for template image mT

    VecField* m_VelocityField;      ///< data container for velocity field (control variable)
    VecField* m_IncVelocityField;   ///< data container for incremental velocity field (incremental control variable)

    Vec m_WorkScaField1;  ///< work scalar field
    Vec m_WorkScaField2;  ///< work scalar field
    Vec m_WorkScaField3;  ///< work scalar field
    Vec m_WorkScaField4;  ///< work scalar field
    Vec m_WorkScaField5;  ///< work scalar field

    Vec m_WorkScaFieldMC;  ///< work scalar field for multi-component/vector fields

    VecField* m_WorkVecField1;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField2;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField3;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField4;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField5;  ///< data container for vector field (temporary variable)

    TenField* m_WorkTenField1;  ///< data container for tensor field (temporary variable)
    TenField* m_WorkTenField2;  ///< data container for tensor field (temporary variable)
    TenField* m_WorkTenField3;  ///< data container for tensor field (temporary variable)
    TenField* m_WorkTenField4;  ///< data container for tensor field (temporary variable)

    ReadWriteType* m_ReadWrite;  ///< io; set from outside (not to be delted)
    RegularizationType* m_Regularization;   ///< regularization functional
    SemiLagrangianType* m_SemiLagrangianMethod;   ///< semi-lagrangian method

    bool m_VelocityIsZero;
    bool m_ComputeInverseDefMap;

 private:
    PetscErrorCode ComputeDefGradSL();                  ///< implemented via SL time integrator
    PetscErrorCode ComputeDetDefGradSL();               ///< implemented via SL time integrator
    PetscErrorCode ComputeDetDefGradRK2();              ///< implemented via RK2 time integrator
    PetscErrorCode ComputeDetDefGradRK2A();             ///< implemented via RK2 time integrator (assymetric form)
    PetscErrorCode ComputeDetDefGradViaDispField();     ///< implemented via RK2 time integrator (asymetric form)
    PetscErrorCode ComputeDeformationMapSLRK2();        ///< implementation via SL time integrator using RK2
    PetscErrorCode ComputeDeformationMapSLRK4();        ///< implementation via SL time integrator using RK4
    PetscErrorCode ComputeDeformationMapRK2();          ///< implementation via RK2 time integrator
    PetscErrorCode ComputeDeformationMapRK2A();         ///< implementation via RK2A time integrator
    PetscErrorCode ComputeDisplacementFieldSL();        ///< implementation via SL time integrator
    PetscErrorCode ComputeDisplacementFieldRK2();       ///< implementation via RK2 time integrator
};




}  // end of namespace




#endif  // _OPTIMALCONTROLREGISTRATIONBASE_H_
