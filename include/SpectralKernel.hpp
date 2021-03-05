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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _SPECTRALKERNEL_HPP_
#define _SPECTRALKERNEL_HPP_

#include "CLAIREUtils.hpp"

namespace reg {

struct SpectralKernel {
  IntType nx[3];
  IntType nl[3];
  IntType nstart[3];
  
  ScalarType scale;
  
  ScalarType *pWS;
  
  PetscErrorCode LowPassFilter(ComplexType *pXHat, ScalarType pct);
  PetscErrorCode HighPassFilter(ComplexType *pXHat, ScalarType pct);
  
  PetscErrorCode Restrict(ComplexType *pXc, const ComplexType *pXf, 
                          const IntType nx_c[3], const IntType osize_c[3], const IntType ostart_c[3]);
  PetscErrorCode Prolong(ComplexType *pXf, const ComplexType *pXc, 
                         const IntType nx_c[3],  const IntType osize_c[3], const IntType ostart_c[3]);
  
  PetscErrorCode Scale(ComplexType *pX, ScalarType val);
  
  PetscErrorCode Norm(ScalarType &norm, ComplexType *pXHat, const IntType w[3]);
};

} // namespace reg

#endif
