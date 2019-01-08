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

#ifndef _TESTDIFFERENTIATION_CPP_
#define _TESTDIFFERENTIATION_CPP_

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "UnitTestOpt.hpp"
#include "Differentiation.hpp"
#include "DifferentiationSM.hpp"
#include "DifferentiationFD.hpp"
#include "VecField.hpp"

namespace reg {
namespace UnitTest {
  
PetscErrorCode TestDifferentiation(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << "starting differentiation unit test" << std::endl;
  
  DifferentiationSM *m_dif = nullptr;
  DifferentiationFD *m_fd = nullptr;
  VecField *v = nullptr;
  VecField *dv = nullptr;
  VecField *ref = nullptr;
  VecField *t = nullptr;
  ScalarType value = 0;
  ScalarType vnorm = 0;
  
  ierr = AllocateOnce(v, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(dv, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(ref, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(t, m_Opt); CHKERRQ(ierr);
  
  ierr = AllocateOnce(m_dif, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(m_fd, m_Opt); CHKERRQ(ierr);
  m_dif->SetupSpectralData();
  m_fd->SetupData();

  
  ierr = ComputeDiffFunction(v, ref, 0, m_Opt); CHKERRQ(ierr); // FD Grad
  ierr = m_fd->Gradient(dv, v->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "grad_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "grad_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "grad_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "grad_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "grad_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "grad_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
 
  ierr = ComputeDiffFunction(v, ref, 0, m_Opt); CHKERRQ(ierr); // FFT grad 
  ierr = m_dif->Gradient(dv, v->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "grad_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "grad_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "grad_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "grad_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "grad_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "grad_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;

  std::cout << std::endl;
  ierr = ComputeDiffFunction(v, ref, 1, m_Opt); CHKERRQ(ierr); // Div
  ierr = m_dif->Divergence(dv->m_X1, v); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "div error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "div error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
  std::cout << std::endl;
  ierr = ComputeDiffFunction(v, ref, 2, m_Opt); CHKERRQ(ierr); // Lap
  ierr = m_dif->Laplacian(dv->m_X1, v->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "lap scalar error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "lap scalar linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
  std::cout << std::endl;
  ierr = m_dif->Laplacian(dv, v); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "lap vector_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "lap vector_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "lap vector_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "lap vector_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "lap vector_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "lap vector_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
  std::cout << std::endl;
  ierr = ComputeDiffFunction(ref, v, 2, m_Opt); CHKERRQ(ierr); // Lap
  ierr = m_dif->InvRegLapOp(dv, v, false, -1.); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "inv lap vector_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "inv lap vector_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  std::cout << "inv lap vector_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "inv lap vector_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "inv lap vector_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  std::cout << "inv lap vector_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
*/  
  ierr = Free(v); CHKERRQ(ierr);
  ierr = Free(dv); CHKERRQ(ierr);
  ierr = Free(m_dif); CHKERRQ(ierr);
  ierr = Free(t); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

}} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

