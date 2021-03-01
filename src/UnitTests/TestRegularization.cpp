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

#ifndef _TESTREGULARIZATION_CPP
#define _TESTREGULARIZATION_CPP

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "UnitTestOpt.hpp"
#include "Regularization.hpp"
#include "RegularizationL2.hpp"
#include "RegularizationH1.hpp"
#include "RegularizationH2.hpp"
#include "RegularizationH3.hpp"
#include "RegularizationH1SN.hpp"
#include "RegularizationH2SN.hpp"
#include "RegularizationH3SN.hpp"
#include "VecField.hpp"

namespace reg {
namespace UnitTest {
  
PetscErrorCode TestRegularization(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << "starting regularization unit test" << std::endl;
  
  Regularization *m_reg = nullptr;
  VecField *v1 = nullptr;
  VecField *t = nullptr;
  ScaField *ts = nullptr;
  ScalarType value = 0;
  ScalarType ref = 0;
  ScalarType beta0 = 0;
  ScalarType beta1 = 0;
  
  ierr = AllocateOnce(v1, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(t, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(ts, m_Opt); CHKERRQ(ierr);
  
  // v = (sin(z)cos(y)sin(y), sin(x)cos(z)sin(z), sin(y)cos(x)sin(x))
  ierr = ComputeSyntheticData(v1, m_Opt, 1); CHKERRQ(ierr);
  
  beta0 = m_Opt->m_RegNorm.beta[0];
  beta1 = m_Opt->m_RegNorm.beta[1];
  
  switch (m_Opt->m_RegNorm.type) {
  case L2:
      ierr = AllocateOnce<RegularizationL2>(m_reg, m_Opt); CHKERRQ(ierr);
      ref = beta0*M_PI*M_PI*M_PI*3./2.;
      break;
  case H1:
      ierr = AllocateOnce<RegularizationH1>(m_reg, m_Opt); CHKERRQ(ierr);
      ref = beta0*M_PI*M_PI*M_PI/2.*(15. + beta1*3.);
      break;
  case H2:
      ierr = AllocateOnce<RegularizationH2>(m_reg, m_Opt); CHKERRQ(ierr);
      ref = beta0*M_PI*M_PI*M_PI/2.*(75. + beta1*15. + beta1*beta1*3.);
      break;
  case H3:
      ierr = AllocateOnce<RegularizationH3>(m_reg, m_Opt); CHKERRQ(ierr);
      break;
  case H1SN:
      ierr = AllocateOnce<RegularizationH1SN>(m_reg, m_Opt); CHKERRQ(ierr);
      ref = beta0*M_PI*M_PI*M_PI*15./2.;
      break;
  case H2SN:
      ierr = AllocateOnce<RegularizationH2SN>(m_reg, m_Opt); CHKERRQ(ierr);
      ref = beta0*M_PI*M_PI*M_PI*75./2.;
      break;
  case H3SN:
      ierr = AllocateOnce<RegularizationH3SN>(m_reg, m_Opt); CHKERRQ(ierr);
      ref = beta0*M_PI*M_PI*M_PI*375./2.;
      break;
  default:
      ierr = reg::ThrowError("regularization model not defined"); CHKERRQ(ierr);
  };
  ierr = m_reg->SetDifferentiation(Differentiation::Type::Spectral); CHKERRQ(ierr);
  ierr = m_reg->SetSpectralData(); CHKERRQ(ierr);
  ierr = m_reg->SetWorkVecField(t); CHKERRQ(ierr);
  ierr = m_reg->SetWorkScaField(ts); CHKERRQ(ierr);
  
  ierr = m_reg->EvaluateFunctional(&value, v1); CHKERRQ(ierr);
  ref *= .5;
  std::cout << "reference is " << ref << std::endl;
  std::cout << "reg functional is " << value << std::endl;
  std::cout << "rel error is " << std::abs(value - ref)/ref << std::endl;
  
  ierr = Free(v1); CHKERRQ(ierr);
  ierr = Free(t); CHKERRQ(ierr);
  ierr = Free(m_reg); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

}} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

