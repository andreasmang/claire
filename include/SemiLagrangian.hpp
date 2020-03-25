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

#ifndef _SEMILAGRANGIAN_HPP_
#define _SEMILAGRANGIAN_HPP_

#ifdef REG_HAS_CUDA
#include "SemiLagrangianGPUNew.hpp"
#else

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"
#include "ReadWriteReg.hpp"
#include "interp3.hpp"




namespace reg {




class SemiLagrangian {
 public:
    SemiLagrangian();
    SemiLagrangian(RegOpt*);
    virtual ~SemiLagrangian();

    PetscErrorCode ComputeTrajectory(VecField*, std::string, ScalarType *mX=nullptr);
    PetscErrorCode SetInitialTrajectory(const ScalarType*);

    /*! interpolate vector field */
    PetscErrorCode Interpolate(VecField*, VecField*, std::string);

    /*! interpolate scalar field */
    PetscErrorCode Interpolate(Vec*, Vec, std::string);

    /*! interpolate scalar field */
    virtual PetscErrorCode Interpolate(ScalarType*, ScalarType*, std::string);

    /*! interpolate vector field */
    virtual PetscErrorCode Interpolate(ScalarType*, ScalarType*, ScalarType*,
                                       ScalarType*, ScalarType*, ScalarType*,
                                       std::string);

    /*! set coordinate vector */
    PetscErrorCode SetQueryPoints(ScalarType*, ScalarType*, ScalarType*, std::string);
    PetscErrorCode GetQueryPoints(ScalarType*, ScalarType*, ScalarType*);
    PetscErrorCode GetQueryPoints(ScalarType*);



    PetscErrorCode SetReadWrite(ReadWriteReg*);
    PetscErrorCode SetWorkVecField(VecField*);

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    virtual PetscErrorCode CommunicateCoord(std::string);
    PetscErrorCode ComputeTrajectoryRK2(VecField*, std::string, ScalarType *mX=nullptr);
    PetscErrorCode ComputeTrajectoryRK4(VecField*, std::string, ScalarType *mX=nullptr);

    RegOpt* m_Opt;

    VecField* m_WorkVecField1;
    VecField* m_WorkVecField2;

    Interp3_Plan* m_AdjointPlan;
    Interp3_Plan* m_StatePlan;

    ScalarType* m_X;
    ScalarType* m_ScaFieldGhost;
    ScalarType* m_VecFieldGhost;

    int m_Dofs[2];

    struct GhostPoints {
        int isize[3];
        int istart[3];
        int nghost;
    };
};




}  // namespace

#endif // REG_HAS_CUDA


#endif
