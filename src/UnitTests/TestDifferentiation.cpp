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
#include <sstream>
#include "UnitTestOpt.hpp"
#include "Differentiation.hpp"
#include "DifferentiationSM.hpp"
#include "DifferentiationFD.hpp"
#include "VecField.hpp"

#ifdef REG_HAS_CUDA
//#define TEST_FD
#endif
//#define SPECTRAL_ANALYSIS

namespace reg {
namespace UnitTest {
  
PetscErrorCode TestDifferentiationMultiGPU(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::ostringstream ss;
  ss << "starting differentiation unit test" << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");
  
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
  
  ierr = AllocateOnce(m_fd, m_Opt); CHKERRQ(ierr);
  
  // Finite Difference Gradient
  ierr = ComputeDiffFunction(v, ref, 0, m_Opt); CHKERRQ(ierr); 
  ierr = m_fd->Gradient(dv, v->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  ss << "FD grad_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");
  
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  ss << "FD grad_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");

  
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  ss << "FD grad_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");
  
  
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  ss << "FD grad_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");
  
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  ss << "FD grad_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");
  
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  ss << "FD grad_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");
  
  // Finite Difference Divergence
  ierr = ComputeDiffFunction(v, ref, 1, m_Opt); CHKERRQ(ierr); 
  ierr = m_fd->Divergence(dv->m_X1, v); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  ss << "FD div error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");
  
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  ss << "FD div error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = Msg(ss.str()); CHKERRQ(ierr);
  ss.str("");
  
  ierr = Free(v); CHKERRQ(ierr);
  ierr = Free(dv); CHKERRQ(ierr);
  ierr = Free(m_fd); CHKERRQ(ierr);
  ierr = Free(t); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode TestDifferentiation(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  if (rank == 0) std::cout << "starting differentiation unit test" << std::endl;
  
#ifdef REG_HAS_CUDA
  printf("rank %2i uses GPU %2i\n", rank, m_Opt->m_gpu_id);
#endif
  
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
  //m_dif->SetupSpectralData();
  //m_fd->SetupData();
  
  ierr = ComputeDiffFunction(v, ref, 0, m_Opt); CHKERRQ(ierr); // Grad
#ifdef TEST_FD  
  ierr = m_fd->Gradient(dv, v->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD grad_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD grad_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD grad_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD grad_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD grad_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD grad_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
#endif
  ierr = m_dif->Gradient(dv, v->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "grad_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << " ref: " << vnorm << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "grad_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "grad_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "grad_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "grad_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "grad_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;

  if (rank == 0) std::cout << std::endl;
  ierr = ComputeDiffFunction(v, ref, 1, m_Opt); CHKERRQ(ierr); // Div
#ifdef TEST_FD
  ierr = m_fd->Divergence(dv->m_X1, v); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD div error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD div error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
#endif
  ierr = m_dif->Divergence(dv->m_X1, v); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "div error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "div error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
  if (rank == 0) std::cout << std::endl;
  ierr = ComputeDiffFunction(v, ref, 2, m_Opt); CHKERRQ(ierr); // Lap
#ifdef TEST_FD
  ierr = m_fd->Laplacian(dv->m_X1, v->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD lap error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD lap error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
#endif
  ierr = m_dif->Laplacian(dv->m_X1, v->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap scalar error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap scalar linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
  if (rank == 0)std::cout << std::endl;
  ierr = ComputeDiffFunction(v, ref, 2, m_Opt); CHKERRQ(ierr); // Lap
#ifdef TEST_FD
  ierr = m_fd->RegLapOp(dv, v, 1.); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, 1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, 1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, 1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD lap vector_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD lap vector_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD lap vector_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD lap vector_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD lap vector_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD lap vector_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
#endif
  ierr = m_dif->RegLapOp(dv, v, 1.); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, 1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, 1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, 1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap vector_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap vector_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap vector_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap vector_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap vector_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap vector_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
  if (rank == 0)std::cout << std::endl;
  ierr = ComputeDiffFunction(ref, v, 2, m_Opt); CHKERRQ(ierr); // Lap
#ifdef TEST_FD
  ierr = m_fd->InvRegLapOp(dv, v, false, -1.); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD inv lap vector_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD inv lap vector_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD inv lap vector_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD inv lap vector_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD inv lap vector_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "FD inv lap vector_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
#endif
  ierr = m_dif->InvRegLapOp(dv, v, false, -1.); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "inv lap vector_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "inv lap vector_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "inv lap vector_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "inv lap vector_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "inv lap vector_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(dv->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "inv lap vector_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
  if (rank == 0)std::cout << std::endl;
  ierr = ComputeDiffFunction(ref, v, 2, m_Opt); CHKERRQ(ierr); // Lap
  ref->Copy(v);
  ierr = m_dif->InvRegLapOp(dv, v, false, -1.); CHKERRQ(ierr);
  ierr = m_dif->RegLapOp(v, dv, -1.); CHKERRQ(ierr);
  ierr = VecAXPY(v->m_X1, -1., ref->m_X1); CHKERRQ(ierr);
  ierr = VecAXPY(v->m_X2, -1., ref->m_X2); CHKERRQ(ierr);
  ierr = VecAXPY(v->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
  ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(v->m_X1, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap inv lap vector_1 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(v->m_X2, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap inv lap vector_2 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(v->m_X3, NORM_2, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap inv lap vector_3 error l2: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X1, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(v->m_X1, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap inv lap vector_1 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X2, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(v->m_X2, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap inv lap vector_2 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  ierr = VecNorm(ref->m_X3, NORM_INFINITY, &vnorm); CHKERRQ(ierr);
  ierr = VecNorm(v->m_X3, NORM_INFINITY, &value); CHKERRQ(ierr);
  if (rank == 0) std::cout << "lap inv lap vector_3 error linf: " << value/(vnorm==0.0?1.:vnorm) << std::endl;
  
#ifdef SPECTRAL_ANALYSIS
  printf("\nSpectral Analysis 1st and 2nd derivative\n  w, FD8, FFT\n");
  for (int w=1; w<m_Opt->m_Domain.nx[2]/2; ++w) {
    printf("%3i",w);
    ierr = ComputeGradSpectral(w, v, ref, m_Opt); CHKERRQ(ierr);
#ifdef TEST_FD
    ierr = m_fd->Gradient(dv, v->m_X1); CHKERRQ(ierr);
    ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
    ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
    ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
    printf(", %.5e", value/(vnorm==0.0?1.:vnorm));
#endif
    ierr = m_dif->Gradient(dv, v->m_X1); CHKERRQ(ierr);
    ierr = VecAXPY(dv->m_X3, -1., ref->m_X3); CHKERRQ(ierr);
    ierr = VecNorm(ref->m_X3, NORM_2, &vnorm); CHKERRQ(ierr);
    ierr = VecNorm(dv->m_X3, NORM_2, &value); CHKERRQ(ierr);
    printf(", %.5e", value/(vnorm==0.0?1.:vnorm));
#ifdef TEST_FD
    ierr = m_fd->Laplacian(dv->m_X1, v->m_X1); CHKERRQ(ierr);
    ierr = VecAXPY(dv->m_X1, -1, ref->m_X1); CHKERRQ(ierr);
    ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
    ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
    printf(", %.5e", value/(vnorm==0.0?1.:vnorm));
#endif
    ierr = m_dif->Laplacian(dv->m_X1, v->m_X1); CHKERRQ(ierr);
    ierr = VecAXPY(dv->m_X1, -1, ref->m_X1); CHKERRQ(ierr);
    ierr = VecNorm(ref->m_X1, NORM_2, &vnorm); CHKERRQ(ierr);
    ierr = VecNorm(dv->m_X1, NORM_2, &value); CHKERRQ(ierr);
    printf(", %.5e\n", value/(vnorm==0.0?1.:vnorm));
  }
#endif

  ierr = Free(v); CHKERRQ(ierr);
  ierr = Free(dv); CHKERRQ(ierr);
  ierr = Free(m_dif); CHKERRQ(ierr);
  ierr = Free(t); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

}} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

