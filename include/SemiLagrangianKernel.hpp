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
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _SEMILAGRANGIANKERNEL_HPP_
#define _SEMILAGRANGIANKERNEL_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"

namespace reg {
 
struct TrajectoryKernel {
  const ScalarType *pV[3];
  ScalarType *pVx[3];
  ScalarType *pX[3];
  
  IntType isize[3];
  IntType istart[3];
  
  ScalarType ix[3];
  ScalarType hx[3];
  
  PetscErrorCode RK2_Step1();
  PetscErrorCode RK2_Step2();
  PetscErrorCode RK2_Step2_inplace();
};

}

#endif
