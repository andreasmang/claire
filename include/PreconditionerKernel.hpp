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

#ifndef _PRECONDITIONERKERNEL_HPP_
#define _PRECONDITIONERKERNEL_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"

namespace reg {

struct H0PrecondKernel {
  const ScalarType *pRHS[3];
  ScalarType *pVhat[3];
  ScalarType *pM[3];
  const ScalarType *pGmt[3];
  const ScalarType *pReg[3];
  
  ScalarType omg;
  ScalarType beta;
  
  IntType nl;
  
  PetscErrorCode Iteration();
  PetscErrorCode IterationPart1();
  PetscErrorCode IterationPart2();
  PetscErrorCode Diagonal();
};

} // namepsace reg

#endif // _PRECONDITIONERKERNEL_HPP_
