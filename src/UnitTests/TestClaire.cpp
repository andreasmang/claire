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

#ifndef _TESTINTERPOLATION_CPP_
#define _TESTINTERPOLATION_CPP_

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "UnitTestOpt.hpp"
#include "CLAIRE.hpp"
#include "VecField.hpp"

namespace reg {
namespace UnitTest {

PetscErrorCode ComputeTrajectoryError(Vec &X, VecField* vinterp, reg::RegOpt* m_Opt, int vcase) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  const ScalarType *p_cX;
  const ScalarType *p_v1 = nullptr, *p_v2=nullptr, *p_v3=nullptr;
  double e1=0,e2=0,e3=0;
  ScalarType x[3], xstar[3], Xtrue[3], v[3], vstar[3];
  ScalarType ht = 1.0/m_Opt->m_Domain.nt;
  ScalarType hthalf = ht/2.0;

  ierr = VecGetArrayRead(X, &p_cX); CHKERRQ(ierr);
  ierr = VecGetArrayRead(vinterp->m_X1, &p_v1); CHKERRQ(ierr);
  ierr = VecGetArrayRead(vinterp->m_X2, &p_v2); CHKERRQ(ierr);
  ierr = VecGetArrayRead(vinterp->m_X3, &p_v3); CHKERRQ(ierr);

  
  for (int i1=0; i1<m_Opt->m_Domain.isize[0]; i1++) {
    for (int i2=0; i2<m_Opt->m_Domain.isize[1]; i2++) {
      for (int i3=0; i3<m_Opt->m_Domain.isize[2]; i3++) {
        int i = GetLinearIndex(i1,i2,i3,m_Opt->m_Domain.isize);
        x[0] = m_Opt->m_Domain.hx[0]*static_cast<ScalarType>(i1 + m_Opt->m_Domain.istart[0]);
        x[1] = m_Opt->m_Domain.hx[1]*static_cast<ScalarType>(i2 + m_Opt->m_Domain.istart[1]);
        x[2] = m_Opt->m_Domain.hx[2]*static_cast<ScalarType>(i3 + m_Opt->m_Domain.istart[2]);
        ierr = ComputeSyntheticVelocity(v, x, vcase); CHKERRQ(ierr);
        xstar[0] = x[0] - ht*v[0];
        xstar[1] = x[1] - ht*v[1];
        xstar[2] = x[2] - ht*v[2];
        ierr = ComputeSyntheticVelocity(vstar, xstar, vcase); CHKERRQ(ierr);
        Xtrue[0] = x[0] - hthalf*(v[0] + vstar[0]);
        Xtrue[1] = x[1] - hthalf*(v[1] + vstar[1]);
        Xtrue[2] = x[2] - hthalf*(v[2] + vstar[2]);
        //if (i > 2e5 && i <(2e5+20)) {
        //  PetscPrintf(PETSC_COMM_WORLD, "Xtrue[0] = %f, x[0] = %f\n", vstar[0], p_v1[i]);
        //}
        e1 = std::max(e1, std::abs(static_cast<double>(Xtrue[0] - p_cX[3*i])));
        e2 = std::max(e2, std::abs(static_cast<double>(Xtrue[1] - p_cX[3*i+1])));
        e3 = std::max(e3, std::abs(static_cast<double>(Xtrue[2] - p_cX[3*i+2])));
        //e1 = std::max(e1, std::abs(static_cast<double>(xstar[0] - p_cX[3*i])));
        //e2 = std::max(e2, std::abs(static_cast<double>(xstar[1] - p_cX[3*i+1])));
        //e3 = std::max(e3, std::abs(static_cast<double>(xstar[2] - p_cX[3*i+2])));
        //e1 = std::max(e1, std::abs(static_cast<double>(vstar[0] - p_v1[i])));
        //e2 = std::max(e2, std::abs(static_cast<double>(vstar[1] - p_v2[i])));
        //e3 = std::max(e3, std::abs(static_cast<double>(vstar[2] - p_v3[i])));
      }
    }
  }
  
  ierr = VecRestoreArrayRead(vinterp->m_X1, &p_v1); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(vinterp->m_X2, &p_v2); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(vinterp->m_X3, &p_v3); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(X, &p_cX); CHKERRQ(ierr);
  
  e1 = std::max(e1, std::max(e2, e3));
  double error = 0;
  MPI_Reduce(&e1, &error, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
  std::ostringstream ss;
  ss << "Linf error for trajectory computation is " << e1 << std::endl;
  Msg(ss.str());
  ss.str("");
  PetscFunctionReturn(ierr);
}

PetscErrorCode TestTrajectoryMultiGPU(reg::RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  std::ostringstream ss;
  ss << "starting trajectory test " << std::endl;
  Msg(ss.str());
  ss.str("");

  reg::VecField* v = nullptr;
  reg::VecField* vwork = nullptr;
  reg::VecField* vinterp = nullptr;
  reg::SemiLagrangian* sl = nullptr;
  Vec X = nullptr;
  ScalarType norm, vnorm;
  ScalarType *p_X;
  const ScalarType *p_cX;
  int vcase = 0;
  
  ierr = VecCreate(X, 3*m_Opt->m_Domain.nl, 3*m_Opt->m_Domain.ng); CHKERRQ(ierr);
  ierr = AllocateOnce(v, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(vinterp, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(vwork, m_Opt); CHKERRQ(ierr);
  try { sl = new SemiLagrangian(m_Opt);}
  catch (std::bad_alloc& err) {
    ierr = reg::ThrowError(err); CHKERRQ(ierr);
  }

  ierr = ComputeSyntheticData(v, m_Opt, vcase); CHKERRQ(ierr);
  ierr = sl->SetWorkVecField(vwork); CHKERRQ(ierr);
  
  // first compute trajectory using given velocity field
  ierr = sl->ComputeTrajectory(v, "state"); CHKERRQ(ierr);
  //ierr = sl->Interpolate(vinterp, v, "state"); CHKERRQ(ierr);

  // extract the trajectory
  ierr = GetRawPointer(X, &p_X); CHKERRQ(ierr);
  ierr = sl->GetQueryPoints(p_X); CHKERRQ(ierr);
  ierr = RestoreRawPointer(X, &p_X); CHKERRQ(ierr);
  // rescale it back to 2*pi cordinates
  ierr = VecScale(X, 2.0*PETSC_PI); CHKERRQ(ierr);
  // everything until here is fine  
  
  ierr = ComputeTrajectoryError(X, vwork, m_Opt, vcase);

  ierr = Free(v);  CHKERRQ(ierr);
  ierr = Free(vinterp); CHKERRQ(ierr);
  ierr = Free(vwork); CHKERRQ(ierr);
  ierr = VecDestroy(&X); CHKERRQ(ierr);
  if (sl != nullptr) { delete sl; sl = nullptr; }

  PetscFunctionReturn(ierr);
}


PetscErrorCode TestTrajectory(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << "starting trajectory solver unit test" << std::endl;
  reg::VecField* v = nullptr;
  reg::VecField* t = nullptr;
//#ifdef REG_HAS_CUDA
//  reg::SemiLagrangianGPUNew *sl = nullptr;
//#else
  reg::SemiLagrangian *sl = nullptr;
//#endif
  //ScalarType *X1 = nullptr;
  //ScalarType *X2 = nullptr;
  //ScalarType *X3 = nullptr;
  //ScalarType *X = nullptr;
  Vec X = NULL;
  ScalarType* pX = nullptr;


  
  IntType nl,ng;
  nl = m_Opt->m_Domain.nl;
  ng = m_Opt->m_Domain.ng;

  ierr = AllocateOnce(v, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(t, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(sl, m_Opt); CHKERRQ(ierr);

  //ierr = AllocateArrayOnce(X1, m_Opt->m_Domain.nl); CHKERRQ(ierr);
  //ierr = AllocateArrayOnce(X2, m_Opt->m_Domain.nl); CHKERRQ(ierr);
  //ierr = AllocateArrayOnce(X3, m_Opt->m_Domain.nl); CHKERRQ(ierr);
  //ierr = AllocateArrayOnce(X, m_Opt->m_Domain.nl*3); CHKERRQ(ierr);
  ierr = VecCreate(X, 3*nl, 3*ng); CHKERRQ(ierr);
  
  ierr = ComputeSyntheticData(v, m_Opt, 3); CHKERRQ(ierr);
  
  ierr = sl->SetWorkVecField(t); CHKERRQ(ierr);
   
  // set X on the CPU
  ierr = VecGetArray(X, &pX); CHKERRQ(ierr);
  for (IntType i1 = 0; i1 < m_Opt->m_Domain.isize[0]; ++i1) {  // x1
    for (IntType i2 = 0; i2 < m_Opt->m_Domain.isize[1]; ++i2) {  // x2
      for (IntType i3 = 0; i3 < m_Opt->m_Domain.isize[2]; ++i3) {  // x3
        ScalarType x1, x2, x3;
        IntType l = GetLinearIndex(i1, i2, i3, m_Opt->m_Domain.isize);
        pX[l*3+0] = m_Opt->m_Domain.hx[0]*static_cast<ScalarType>(i1 + m_Opt->m_Domain.istart[0]);
        pX[l*3+1] = m_Opt->m_Domain.hx[1]*static_cast<ScalarType>(i2 + m_Opt->m_Domain.istart[1]);
        pX[l*3+2] = m_Opt->m_Domain.hx[2]*static_cast<ScalarType>(i3 + m_Opt->m_Domain.istart[2]);
      }  // i1
    }  // i2
  }  // i3
 
 ierr = VecRestoreArray(X, &pX);
 

#if defined(REG_HAS_CUDA) || defined(REG_HAS_MPICUDA)

  ierr = VecCUDAGetArray(X, &pX); CHKERRQ(ierr);
  ierr = sl->SetInitialTrajectory(pX); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(X, &pX); CHKERRQ(ierr);
  //std::cout << "Initial Trajectory Set" << std::endl;

  ierr = sl->ComputeTrajectory(v, "state"); CHKERRQ(ierr);
  //std::cout << "Trajectory Computed" << std::endl;
  
  ierr = sl->Interpolate(t, v, "state"); CHKERRQ(ierr);
  //std::cout << "Velocity interpolated" << std::endl;
  
  ierr = v->Copy(t); CHKERRQ(ierr);
  //std::cout << "copied t -> v" << std::endl;

  ierr = VecCUDAGetArray(X, &pX); CHKERRQ(ierr);
  ierr = sl->GetQueryPoints(pX);
  ierr = VecCUDARestoreArray(X, &pX); CHKERRQ(ierr);
  //std::cout << "Got query points" << std::endl;

  ierr = v->Scale(-1.); CHKERRQ(ierr);
  //std::cout << "v scaled by -1" << std::endl;

  ierr = VecCUDAGetArray(X, &pX); CHKERRQ(ierr);
  ierr = sl->SetInitialTrajectory(pX); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(X, &pX); CHKERRQ(ierr);
  //std::cout << "Initial trajectory reset" << std::endl;

  ierr = sl->ComputeTrajectory(v, "state"); CHKERRQ(ierr);
  //std::cout << "trajectory computed again with -v" << std::endl;

#else
  ierr = sl->ComputeTrajectory(v, "state", X); CHKERRQ(ierr);
  ierr = sl->Interpolate(t, v, "state"); CHKERRQ(ierr);
  ierr = v->Copy(t); CHKERRQ(ierr);
  ierr = sl->GetQueryPoints(X);
  ierr = v->Scale(-1.); CHKERRQ(ierr);
  ierr = sl->ComputeTrajectory(v, "state", X); CHKERRQ(ierr);
#endif
    
  ierr = VecCUDAGetArray(X, &pX); CHKERRQ(ierr);
  ierr = sl->GetQueryPoints(pX);
  ierr = VecCUDARestoreArray(X, &pX); CHKERRQ(ierr);

  ierr = VecGetArray(X, &pX); CHKERRQ(ierr);

  double ex = 0.0;
  double ey = 0.0;
  double ez = 0.0;
  for (IntType i1 = 0; i1 < m_Opt->m_Domain.isize[0]; ++i1) {  // x1
    for (IntType i2 = 0; i2 < m_Opt->m_Domain.isize[1]; ++i2) {  // x2
      for (IntType i3 = 0; i3 < m_Opt->m_Domain.isize[2]; ++i3) {  // x3
        // compute coordinates (nodal grid)
        ScalarType x1, x2, x3;
        x1 = m_Opt->m_Domain.hx[0]*static_cast<ScalarType>(i1 + m_Opt->m_Domain.istart[0]);
        x2 = m_Opt->m_Domain.hx[1]*static_cast<ScalarType>(i2 + m_Opt->m_Domain.istart[1]);
        x3 = m_Opt->m_Domain.hx[2]*static_cast<ScalarType>(i3 + m_Opt->m_Domain.istart[2]);

        // compute linear / flat index
        IntType l = GetLinearIndex(i1, i2, i3, m_Opt->m_Domain.isize);

        ex = std::max(ex,std::abs(static_cast<double>(pX[l*3+0] - x1)));
        ey = std::max(ey,std::abs(static_cast<double>(pX[l*3+1] - x2)));
        ez = std::max(ez,std::abs(static_cast<double>(pX[l*3+2] - x3)));

      }  // i1
    }  // i2
  }  // i3
    
  ierr = VecRestoreArray(X, &pX); CHKERRQ(ierr);
  std::cout << "Inf-norm on trajectory: " << std::max(ex,std::max(ey,ez)) << std::endl;
  
  
  ierr = Free(sl); CHKERRQ(ierr);
  ierr = Free(t); CHKERRQ(ierr);
  ierr = Free(v); CHKERRQ(ierr);
  //ierr = FreeArray(X1); CHKERRQ(ierr);
  //ierr = FreeArray(X2); CHKERRQ(ierr);
  //ierr = FreeArray(X3); CHKERRQ(ierr);
  //ierr = FreeArray(X); CHKERRQ(ierr);
  ierr = VecDestroy(&X); CHKERRQ(ierr);
  pX = nullptr;
    
  PetscFunctionReturn(ierr);
}

PetscErrorCode TestForwardSolver(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = reg::DbgMsg("starting forward solver unit test"); CHKERRQ(ierr);
  
  reg::ScaField *m0 = NULL, *m0true = NULL;
  reg::VecField* v = NULL;
  reg::VecField* t = NULL;
  reg::VecField* t2 = NULL;
  reg::TransportProblem* solver = NULL;
  ScalarType val, val0, relval;
  ScalarType global_val, global_val0, global_relval;

  // make sure we do not store time history
  m_Opt->m_RegFlags.runinversion = false;
  
  ierr = AllocateOnce(v, m_Opt);
  ierr = AllocateOnce(t, m_Opt);
  ierr = AllocateOnce(t2, m_Opt);
  ierr = AllocateOnce(m0, m_Opt);
  ierr = AllocateOnce(m0true, m_Opt);

  ierr = ComputeSyntheticData(*m0true, m_Opt); CHKERRQ(ierr);
  ierr = m0->Copy(m0true); CHKERRQ(ierr);
  ierr = ComputeSyntheticData(v, m_Opt); CHKERRQ(ierr);

  ierr = AllocateOnce<TransportEquationSL>(solver, m_Opt); CHKERRQ(ierr);
  
  ierr = solver->SetControlVariable(v); CHKERRQ(ierr);
  ierr = solver->SetStateVariable(m0); CHKERRQ(ierr);
  ierr = solver->SetWorkVecField(t, 1); CHKERRQ(ierr);
  ierr = solver->SetWorkVecField(t2, 2); CHKERRQ(ierr);

  ierr = reg::DbgMsg("computing error for forward solver"); CHKERRQ(ierr);

  ierr = solver->SolveForwardProblem(); CHKERRQ(ierr);
  ierr = v->Scale(-1.);                 CHKERRQ(ierr);
  ierr = solver->SetControlVariable(v); CHKERRQ(ierr);
  ierr = solver->SolveForwardProblem(); CHKERRQ(ierr);

  ierr = VecAXPY(*m0, -1.0, *m0true); CHKERRQ(ierr);
  ierr = VecNorm(*m0, NORM_2, &val); CHKERRQ(ierr);
  ierr = VecNorm(*m0true, NORM_2, &val0); CHKERRQ(ierr);
  
  relval = val;
  relval /= (val0 > 0.0 ? val0 : 1.0);
  
  MPI_Reduce(&relval, &global_relval, 1, MPI_FLOAT, MPI_MAX, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&val, &global_val, 1, MPI_FLOAT, MPI_MAX, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&val0, &global_val0, 1, MPI_FLOAT, MPI_MAX, 0, PETSC_COMM_WORLD);
  
  std::ostringstream ss;
  ss << "numerical error: "<< std::scientific << relval
         << " (absolute " << val << ",ref: " << val0 << ")" << std::endl;
  Msg(ss.str());
  
  double local_runtime, global_runtime;
  Msg("ZeitGeist:");
  for (auto zg : ZeitGeist::zgMap()) {
    char txt[120];
    local_runtime = zg.second.Total_s();
    MPI_Reduce(&local_runtime, &global_runtime, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
    sprintf(txt, "  %16s: %5lix, %0.10lf",zg.first.c_str(), zg.second.Count(), global_runtime);
    Msg(txt);
  }
  Msg("-----------------------------------------------------------------------------------------------------");

  ierr = Free(v); CHKERRQ(ierr);
  ierr = Free(t); CHKERRQ(ierr);
  ierr = Free(t2); CHKERRQ(ierr);
  ierr = Free(m0); CHKERRQ(ierr);
  ierr = Free(m0true); CHKERRQ(ierr);
  ierr = Free(solver); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode TestGradient(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << "starting gradient solver unit test" << std::endl;
  
  reg::ScaField* m = nullptr;
  reg::VecField* v = nullptr;
  reg::CLAIRE* registration = nullptr;

  ierr = AllocateOnce(registration, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(v, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(m, m_Opt); CHKERRQ(ierr);
  ierr = ComputeSyntheticData(m->m_X, m_Opt); CHKERRQ(ierr);
  ierr = ComputeSyntheticData(v, m_Opt); CHKERRQ(ierr);

  ierr = registration->SetControlVariable(v); CHKERRQ(ierr);
  ierr = registration->InitializeSolver(); CHKERRQ(ierr);
  ierr = registration->SetReferenceImage(m->m_X); CHKERRQ(ierr);
  ierr = registration->SolveForwardProblem(NULL, m->m_X); CHKERRQ(ierr);
  ierr = registration->EvaluateGradient(NULL, NULL); CHKERRQ(ierr);

  ierr = registration->DerivativeCheckGradient(); CHKERRQ(ierr);

  if (m != NULL) {ierr = Free(m); CHKERRQ(ierr); m = NULL;}
  if (v != NULL) {ierr = Free(v); CHKERRQ(ierr); v = NULL;}
  
  ierr = Free(registration); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}


PetscErrorCode TestHessian(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << "starting gradient solver unit test" << std::endl;
  
  reg::ScaField* m = nullptr;
  reg::VecField* v = nullptr;
  reg::CLAIRE* registration = nullptr;

  ierr = AllocateOnce(registration, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(m, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(v, m_Opt);
  ierr = ComputeSyntheticData(m->m_X, m_Opt); CHKERRQ(ierr);
  ierr = ComputeSyntheticData(v, m_Opt); CHKERRQ(ierr);

  ierr = registration->SetControlVariable(v); CHKERRQ(ierr);
  ierr = registration->InitializeSolver(); CHKERRQ(ierr);
  ierr = registration->SetReferenceImage(m->m_X); CHKERRQ(ierr);
  ierr = registration->SolveForwardProblem(NULL, m->m_X); CHKERRQ(ierr);
  ierr = registration->EvaluateGradient(NULL, NULL); CHKERRQ(ierr);
  ierr = registration->HessianMatVec(NULL, NULL); CHKERRQ(ierr);

  ierr = registration->DerivativeCheckHessian(); CHKERRQ(ierr);

  if (m != NULL) {ierr = Free(m); CHKERRQ(ierr); m = NULL;}
  if (v != NULL) {ierr = Free(v); CHKERRQ(ierr); v = NULL;}
  
  ierr = Free(registration); CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}

}} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

