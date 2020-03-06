/*************************************************************************
 *  Copyright (c) 2018.
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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _DEFORMATIONFIELDS_HPP_
#define _DEFORMATIONFIELDS_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"
#include "TenField.hpp"
#include "Differentiation.hpp"
#include "ReadWriteReg.hpp"
//#ifdef REG_HAS_CUDA
//#include "SemiLagrangianGPUNew.hpp"
//#endif
#include "SemiLagrangian.hpp"




namespace reg {




class DeformationFields {
 public:
    typedef DeformationFields Self;
//#ifdef REG_HAS_CUDA
//    typedef SemiLagrangianGPUNew SemiLagrangianType;
//#else
    typedef SemiLagrangian SemiLagrangianType;
//#endif

    DeformationFields();
    DeformationFields(RegOpt*);
    virtual ~DeformationFields();

    /*! set velocity field */
    PetscErrorCode SetVelocityField(VecField*);

    /*! set work vector fields */
    PetscErrorCode SetWorkScaField(Vec, int);

    /*! set work vector fields */
    PetscErrorCode SetWorkVecField(VecField*, int);

    /*! set io object */
    PetscErrorCode SetReadWrite(ReadWriteReg*);
    
    /*! set differentiation object */
    PetscErrorCode SetDifferentiation(Differentiation*);

    /*! set semi-lagrangian method */
    PetscErrorCode SetSLM(SemiLagrangianType*);

    /*! compute determinant of deformation gradient, i.e.
        the jacobian of the deformation map */
    PetscErrorCode ComputeDetDefGrad();

    /*! compute deformation map */
    PetscErrorCode ComputeDeformationMap(bool write2file = false, VecField* y = NULL);

    /*! compute displacement field */
    PetscErrorCode ComputeDisplacementField(bool write2file = false);

    /*! compute deformation gradient, i.e. jacobian of deformation map */
    PetscErrorCode ComputeDefGrad(bool write2file = false);

    /*! check deformation map computation */
    PetscErrorCode CheckDefMapConsistency();

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

private:
    PetscErrorCode ComputeRegularGrid(VecField*);
    PetscErrorCode ComputeDefMapFromDisplacement();  ///< compute deformation map from displacement
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

    SemiLagrangianType* m_SemiLagrangianMethod;  ///< semi-lagrangian method

    VecField* m_VelocityField;  ///< velocity field

    Vec m_WorkScaField1;  ///< work scalar field
    Vec m_WorkScaField2;  ///< work scalar field
    Vec m_WorkScaField3;  ///< work scalar field
    Vec m_WorkScaField4;  ///< work scalar field
    Vec m_WorkScaField5;  ///< work scalar field

    VecField* m_WorkVecField1;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField2;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField3;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField4;  ///< data container for vector field (temporary variable)
    VecField* m_WorkVecField5;  ///< data container for vector field (temporary variable)

    TenField* m_WorkTenField1;  ///< data container for tensor field (temporary variable)
    TenField* m_WorkTenField2;  ///< data container for tensor field (temporary variable)
    TenField* m_WorkTenField3;  ///< data container for tensor field (temporary variable)
    TenField* m_WorkTenField4;  ///< data container for tensor field (temporary variable)

    ReadWriteReg* m_ReadWrite; ///< io; set from outside (not to be delted)
    Differentiation* m_Differentiation; ///< Differentiation interface; set from outside (not to be delted)

    bool m_ComputeInverseDefMap;

    RegOpt* m_Opt;
};




}  // namespace reg




#endif  // _DEFORMATIONFIELDS_HPP_
