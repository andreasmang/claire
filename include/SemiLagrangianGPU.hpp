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



#ifndef _SEMILAGRANGIANGPU_H_
#define _SEMILAGRANGIANGPU_H_

#include "RegOpt.h"
#include "RegUtils.h"
#include "VecField.h"


size_t accfft_ghost_local_size_dft_r2c(accfft_plan* plan,int g_size, int * isize_g, int* istart_g);
void accfft_get_ghost(accfft_plan* plan,int g_size,int* isize_g, double* data,double* ghost_data);

size_t accfft_ghost_xyz_local_size_dft_r2c(accfft_plan* plan,int g_size, int * isize_g, int* istart_g);
void accfft_get_ghost_xyz(accfft_plan* plan,int g_size,int* isize_g, double* data,double* ghost_data);


#include "interp3_gpu_mpi.hpp"
#include "interp3.hpp"

namespace reg
{


class SemiLagrangianGPU
{
public:

    typedef SemiLagrangianGPU Self;

    SemiLagrangianGPU();
    SemiLagrangianGPU(RegOpt*);
    virtual ~SemiLagrangianGPU();

    PetscErrorCode ComputeTrajectory(VecField*,std::string);

    /*! interpolate vector field */
    PetscErrorCode Interpolate(VecField*,VecField*,std::string);

    /*! interpolate scalar field */
    PetscErrorCode Interpolate(Vec*,Vec,std::string);

    /*! interpolate scalar field */
    virtual PetscErrorCode Interpolate(ScalarType*,ScalarType*,std::string);

    /*! interpolate vector field */
    virtual PetscErrorCode Interpolate(ScalarType*,ScalarType*,ScalarType*,
                                       ScalarType*,ScalarType*,ScalarType*,
                                       std::string);

protected:

    PetscErrorCode Initialize();
    PetscErrorCode ComputeInitialCondition();
    virtual PetscErrorCode MapCoordinateVector(std::string);
    PetscErrorCode ClearMemory();

    VecField* m_TrajectoryS;
    VecField* m_TrajectoryA;
    VecField* m_InitialTrajectory;
    VecField* m_WorkVecField;

    ScalarType* m_XA;
    ScalarType* m_XS;
    ScalarType* m_iVecField;
    ScalarType* m_xVecField;

    Interp3_Plan_GPU* m_AdjointPlan;
    Interp3_Plan_GPU* m_StatePlan;
    Interp3_Plan_GPU* m_AdjointPlanVec;
    Interp3_Plan_GPU* m_StatePlanVec;

    ScalarType* m_ScaFieldGhost;
    ScalarType* m_VecFieldGhost;
    static const int m_GhostSize = 3;

    RegOpt* m_Opt;

};

} // namespace


#endif
