/*************************************************************************
 *  Copyright (c) 2017.
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

#ifndef _UNITTESTOPT_H_
#define _UNITTESTOPT_H_

#include "RegOpt.hpp"
#include "VecField.hpp"

namespace reg {

class UnitTestOpt : public RegOpt {
 public:
    typedef UnitTestOpt Self;
    typedef RegOpt SuperClass;

    UnitTestOpt();
    UnitTestOpt(int, char**);
    UnitTestOpt(const UnitTestOpt&);
    virtual ~UnitTestOpt();

    virtual PetscErrorCode Run();
    
    int rep;

 protected:
    virtual PetscErrorCode Initialize(void);
    virtual PetscErrorCode ClearMemory(void);
    virtual PetscErrorCode ParseArguments(int, char**);
    virtual PetscErrorCode Usage(bool advanced = false);
    virtual PetscErrorCode CheckArguments(void);
    
    virtual PetscErrorCode PrintTests();
    virtual PetscErrorCode CheckTests(char*, bool&);
    
    typedef enum {
      None=-2, 
      All=-1, 
      Interpolate,
      ForwardSolver,
      Trajectory,
      Gradient,
      Hessian,
      Reg,
      Diff
    } TestType;
    
    TestType m_TestType;
};

namespace UnitTest {
  PetscErrorCode ComputeTrajectoryError(Vec &, reg::RegOpt*, int);
  PetscErrorCode ComputeSyntheticVelocity(ScalarType*, ScalarType*, int);
  PetscErrorCode ComputeSyntheticData(Vec& m, RegOpt* opt);
  PetscErrorCode ComputeSyntheticData(VecField*&, RegOpt*, IntType=1);
  PetscErrorCode ComputeDiffFunction(VecField *, VecField *, int, RegOpt*);
  PetscErrorCode ComputeGradSpectral(ScalarType, VecField *, VecField *, RegOpt*);
  
  PetscErrorCode TestInterpolation(RegOpt *m_Opt);
  PetscErrorCode TestInterpolationMultiGPU(RegOpt *m_Opt);
  PetscErrorCode TestVectorFieldInterpolationMultiGPU(RegOpt *m_Opt);
  PetscErrorCode TestForwardSolver(RegOpt *m_Opt);
  PetscErrorCode TestTrajectory(RegOpt *m_Opt);
  PetscErrorCode TestTrajectoryMultiGPU(RegOpt *m_Opt);
  PetscErrorCode TestGradient(RegOpt *m_Opt);
  PetscErrorCode TestHessian(RegOpt *m_Opt);
  PetscErrorCode TestRegularization(RegOpt *m_Opt);
  PetscErrorCode TestDifferentiation(UnitTestOpt *m_Opt);
}


}  // namespace reg




#endif  // _UNITTESTOPT_H_
