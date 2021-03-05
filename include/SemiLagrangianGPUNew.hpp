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

#ifndef _SEMILAGRANGIANGPUNEW_HPP_
#define _SEMILAGRANGIANGPUNEW_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"
#include "ReadWriteReg.hpp"
#include "interp3.hpp"
#include "interp3_gpu_mpi.hpp"
#include "GhostPlan.hpp"




namespace reg {




class SemiLagrangianGPUNew {
 public:
    SemiLagrangianGPUNew();
    SemiLagrangianGPUNew(RegOpt*);
    virtual ~SemiLagrangianGPUNew();

    PetscErrorCode ComputeTrajectory(VecField*, std::string);
    PetscErrorCode ComputeInitialTrajectory();
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

    PetscErrorCode SetReadWrite(ReadWriteReg*);
    PetscErrorCode SetWorkVecField(VecField*);
    PetscErrorCode SetQueryPoints(ScalarType*, ScalarType*, ScalarType*, std::string);
    PetscErrorCode GetQueryPoints(ScalarType*, ScalarType*, ScalarType*);
    PetscErrorCode GetQueryPoints(ScalarType*);

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();
    PetscErrorCode InitializeInterpolationTexture();

    virtual PetscErrorCode CommunicateCoord(std::string);
    PetscErrorCode ComputeTrajectoryRK2(VecField*, std::string);
    PetscErrorCode ComputeTrajectoryRK4(VecField*, std::string);
    PetscErrorCode MapCoordinateVector(std::string);

    RegOpt* m_Opt;

    VecField* m_WorkVecField1;
    VecField* m_Xstate;
    VecField* m_Xadjoint;
    VecField* m_InitialTrajectory;
    
    Vec m_X;
    
    cudaTextureObject_t m_texture;
    
    ScalarType* m_tmpInterpol1;
    ScalarType* m_tmpInterpol2;
    
    Interp3_Plan_GPU* m_StatePlan;
    Interp3_Plan_GPU* m_AdjointPlan;
    
    ScalarType* m_ScaFieldGhost;
    ScalarType* m_VecFieldGhost;
    ScalarType* m_WorkScaField1;
    ScalarType* m_WorkScaField2;
    
    GhostPlan* m_GhostPlan;
    
    int m_Dofs[2];
    IntType nghost, nlghost;
    IntType istart_g[3], isize_g[3];
    size_t g_alloc_max; 

    bool cuda_aware = true;
    
};

typedef SemiLagrangianGPUNew SemiLagrangian;




}  // namespace




#endif
