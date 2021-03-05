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
#include <sstream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "UnitTestOpt.hpp"
#include "interp3.hpp"
#if defined(REG_HAS_MPICUDA) || defined(REG_HAS_CUDA)
#include "interp3_gpu_mpi.hpp"
#endif
#include "zeitgeist.hpp"
#include "GhostPlan.hpp"
#define CONSTANT

void TestFunction(ScalarType &val, const ScalarType x1, const ScalarType x2, const ScalarType x3) {
  val =  ( PetscSinReal(8*x1)*PetscSinReal(8*x1)
        + PetscSinReal(2*x2)*PetscSinReal(2*x2)
        + PetscSinReal(4*x3)*PetscSinReal(4*x3) )/3.0;
}

inline void printOnMaster(std::string msg) {
  PetscPrintf(PETSC_COMM_WORLD, msg.c_str());
}

static void printGPUMemory(int rank) {
    if (rank == 0) {
      size_t free, used;
      cudaMemGetInfo(&free, &used);
      used -= free;
      std::string msg = "Used mem = " + std::to_string(used/1E9) + " GB, Free mem = " + std::to_string(free/1E9) + " GB\n";
      PetscPrintf(PETSC_COMM_WORLD, msg.c_str());
    }
}

template<typename T>
void printVector(T* x, int n, std::string name) {
    std::string msg;
    msg = name +  ": ";
    for (int i=0; i<n; i++) {
        msg += (std::to_string(x[i]) + ",");
    }
    msg += "\n";
    PetscPrintf(PETSC_COMM_WORLD, msg.c_str());
}

template<typename T>
void printScalar(T x, std::string name) {
    std::string msg = name + ": " + std::to_string(x) + "\n";
    PetscPrintf(PETSC_COMM_WORLD, msg.c_str());
}

void TestError(ScalarType *ref, ScalarType *eval, IntType nl, double *err, double *maxval) {
  double error = 0;
  double max = 0;
  for (IntType i = 0; i < nl; ++i) {
    double local = static_cast<double>(ref[i]) - static_cast<double>(eval[i]);
    //PetscPrintf(PETSC_COMM_WORLD, "ref[i] = %f, eval[i] = %f\n", ref[i], eval[i]);
    double lmax = std::abs(eval[i]);
    if (lmax > max) max = lmax;
    error += local*local;
  }
  *err = sqrt(error/static_cast<double>(nl));
  *maxval = max;
}

void TestErrorMultiGPU(ScalarType *ref, ScalarType *eval, IntType nl, double *err, double *maxval) {
  double error = 0;
  double max = 0;
  for (IntType i = 0; i < nl; ++i) {
    double local = static_cast<double>(ref[i]) - static_cast<double>(eval[i]);
    //PetscPrintf(PETSC_COMM_WORLD, "ref[i] = %f, eval[i] = %f\n", ref[i], eval[i]);
    double lmax = std::abs(eval[i]);
    if (lmax > max) max = lmax;
    error += local*local;
  }
  *err = error;
  *maxval = max;
}


namespace reg {
namespace UnitTest {

PetscErrorCode TestInterpolation(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << "starting interpolation unit test" << std::endl;
  
  srand(time(0));
  
  IntType nx[3], istart[3], isize[3], nl;
  ScalarType hx[3], scale;
  
  nl = 1;
  for (int i = 0; i < 3; ++i) {
    nx[i]     = m_Opt->m_Domain.nx[i];
    isize[i]  = m_Opt->m_Domain.isize[i];
    istart[i] = m_Opt->m_Domain.istart[i];
    hx[i]     = m_Opt->m_Domain.hx[i];
    nl *= isize[i];
  }
    
  ScalarType *grid = new ScalarType[nl];
  ScalarType *q    = new ScalarType[nl*3];
  ScalarType *q1   = &q[0];
  ScalarType *q2   = &q[nl];
  ScalarType *q3   = &q[nl*2];
  ScalarType *ref  = new ScalarType[nl];
  ScalarType *eval = new ScalarType[nl];
  
  int iq = 10;

  for (int i1 = 0; i1 < isize[0]; ++i1) {  // x1
    for (int i2 = 0; i2 < isize[1]; ++i2) {  // x2
      for (int i3 = 0; i3 < isize[2]; ++i3) {  // x3
        // compute coordinates (nodal grid)
        double x1 = hx[0]*static_cast<double>(i1 + istart[0]);
        double x2 = hx[1]*static_cast<double>(i2 + istart[1]);
        double x3 = hx[2]*static_cast<double>(i3 + istart[2]);

        // compute linear / flat index
        int i = reg::GetLinearIndex(i1, i2, i3, m_Opt->m_Domain.isize);
        
        TestFunction(grid[i], x1, x2, x3);
        
        x1 += hx[0]*0.5;
        x2 += hx[1]*0.5;
        x3 += hx[2]*0.5;
        
        if (x1 > 2.*M_PI) x1 -= 2.*M_PI;
        if (x2 > 2.*M_PI) x2 -= 2.*M_PI;
        if (x3 > 2.*M_PI) x3 -= 2.*M_PI;
        
        TestFunction(ref[i], x1, x2, x3);
        eval[i] = 0.;
    
#ifdef REG_HAS_CUDA
        q1[i] = x1*nx[0]/(2.*M_PI);
        q2[i] = x2*nx[1]/(2.*M_PI);
        q3[i] = x3*nx[2]/(2.*M_PI);
#else
        q[3*i + 0] = x1;
        q[3*i + 1] = x2;
        q[3*i + 2] = x3;
#endif
      }  // i1
    }  // i2
  }  // i3


#ifdef REG_HAS_CUDA
  ScalarType *pg, *pq1, *pq2, *pq3, *pe;
  float *tmp1, *tmp2;
  
  cudaMalloc((void**)(&pg), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&pq1), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&pq2), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&pq3), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&pe), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&tmp1), sizeof(ScalarType)*nl);
  cudaMalloc((void**)(&tmp2), sizeof(ScalarType)*nl);
  
  cudaMemcpy(pg, grid, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  cudaMemcpy(pq1, q1, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  cudaMemcpy(pq2, q2, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  cudaMemcpy(pq3, q3, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  
  cudaTextureObject_t tex = gpuInitEmptyTexture(nx);
  
  const ScalarType * pq[3] = {pq1, pq2, pq3};
  float timer = 0;
  for (int i=0; i<10; ++i) {
    gpuInterp3D(pg, pq, pe, tmp1, tmp2, nx, nl, tex, m_Opt->m_PDESolver.iporder, &timer);
    cudaDeviceSynchronize();
  }
  
  double error, max;
  error = 0;
  max = 0;
  
  ZeitGeist_define(INTERPOL);
  for (int i=0; i<100; ++i) {
    double err, m;
    ZeitGeist_tick(INTERPOL);
    gpuInterp3D(pg, pq, pe, tmp1, tmp2, nx, nl, tex, m_Opt->m_PDESolver.iporder, &timer);
    //interp0(pg,pq1,pq2,pq3,pe,nx);
    cudaDeviceSynchronize();
    ZeitGeist_tock(INTERPOL);
    cudaMemcpy(eval, pe, sizeof(ScalarType)*nl, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    TestError(ref, eval, nl, &err, &m);
    error += err;
    max += m;
  }
  
  error /= 100;
  max /= 100;
  
  printf("On Grid\n");
  std::cout << "\tthe interpolation has a RMS of " << error << std::endl;
  std::cout << "\tabs max of interpolation is " << max << std::endl;
  
  for (int i1 = 0; i1 < isize[0]; ++i1) {  // x1
      for (int i2 = 0; i2 < isize[1]; ++i2) {  // x2
        for (int i3 = 0; i3 < isize[2]; ++i3) {  // x3
          // compute linear / flat index
          int i = reg::GetLinearIndex(i1, i2, i3, m_Opt->m_Domain.isize);
          
          double x1 = hx[0]*static_cast<double>(i1 + istart[0]);
          double x2 = hx[1]*static_cast<double>(i2 + istart[1]);
          double x3 = hx[2]*static_cast<double>(i3 + istart[2]);
      
          x1 += hx[0]*static_cast<double>(rand()-RAND_MAX/2)/static_cast<double>(RAND_MAX/2);
          x2 += hx[1]*static_cast<double>(rand()-RAND_MAX/2)/static_cast<double>(RAND_MAX/2);
          x3 += hx[2]*static_cast<double>(rand()-RAND_MAX/2)/static_cast<double>(RAND_MAX/2);
          //double x1 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
          //double x2 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
          //double x3 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
          
          if (x1 > 2.*M_PI) x1 -= 2.*M_PI;
          if (x2 > 2.*M_PI) x2 -= 2.*M_PI;
          if (x3 > 2.*M_PI) x3 -= 2.*M_PI;
          if (x1 < 0.) x1 += 2.*M_PI;
          if (x2 < 0.) x2 += 2.*M_PI;
          if (x3 < 0.) x3 += 2.*M_PI;
          
          TestFunction(ref[i], x1, x2, x3);
          
          q1[i] = x1*nx[0]/(2.*M_PI);
          q2[i] = x2*nx[1]/(2.*M_PI);
          q3[i] = x3*nx[2]/(2.*M_PI);
        }  // i1
      }  // i2
    }  // i3
  cudaMemcpy(pq1, q1, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  cudaMemcpy(pq2, q2, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);
  cudaMemcpy(pq3, q3, sizeof(ScalarType)*nl, cudaMemcpyHostToDevice);

  pq[0] = pq1;
  pq[1] = pq2;
  pq[2] = pq3;
  
  for (int i=0; i<10; ++i) {
    gpuInterp3D(pg, pq, pe, tmp1, tmp2, nx, nl, tex, m_Opt->m_PDESolver.iporder, &timer);
    cudaDeviceSynchronize();
  }
  
  error = 0;
  max = 0;
  
  ZeitGeist_define(INTERPOL_RAND);
  for (int i=0; i<100; ++i) {
    double err, m;
    ZeitGeist_tick(INTERPOL_RAND);
    gpuInterp3D(pg, pq, pe, tmp1, tmp2, nx, nl, tex, m_Opt->m_PDESolver.iporder, &timer);
    //interp0(pg,pq1,pq2,pq3,pe,nx);
    cudaDeviceSynchronize();
    ZeitGeist_tock(INTERPOL_RAND);
    cudaMemcpy(eval, pe, sizeof(ScalarType)*nl, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    TestError(ref, eval, nl, &err, &m);
    error += err;
    max += m;
  }
  
  error /= 100;
  max /= 100;
  
  printf("Random Eval\n");
  std::cout << "\tthe interpolation has a RMS of " << error << std::endl;
  std::cout << "\tabs max of interpolation is " << max << std::endl;

  cudaFree(pg);
  cudaFree(q1);
  cudaFree(q2);
  cudaFree(q3);
  cudaFree(pe);
  cudaFree(tmp1);
  cudaFree(tmp2);
  cudaDestroyTextureObject(tex);
#else
  std::cout << "unit test not implemented" << std::endl;
#endif

  printf("Timings\n");
  for (auto zg : ZeitGeist::zgMap()) {
      char txt[120];
      sprintf(txt, "  %16s: %5lix, %10lf s",zg.first.c_str(), zg.second.Count(), zg.second.Total_s());
      reg::Msg(txt);
  }
  
  delete[] grid;
  delete[] q;
  delete[] ref;
  delete[] eval;
  
  PetscFunctionReturn(ierr);
}

/********************************************************************************************************************/
PetscErrorCode TestInterpolationMultiGPU(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  printOnMaster("starting interpolation unit test\n");
  
  srand(time(0));
  
  IntType nx[3], isize[3], istart[3];
  IntType nl = 1, ng = 1;
  ScalarType hx[3], scale;
  double error=0, max=0;
  ScalarType m_GPUtime;
  int nghost = 3;
  IntType isize_g[3], istart_g[3];
  IntType nlghost = 1;
  double timers[4] = {0,0,0,0};
  double global_error, global_max;
  ScalarType* m_tmpInterpol1 = nullptr, *m_tmpInterpol2 = nullptr;
  double global_timers[4] = {0, 0, 0, 0};
  int nprocs, procid;
  PetscRandom rctx;
    
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  Vec q = nullptr; ScalarType* p_q; // query points
  Vec xq = nullptr; ScalarType* p_xq; // query points
  Vec yq = nullptr; ScalarType* p_yq; // query points
  Vec zq = nullptr; ScalarType* p_zq; // query points
  Vec f = nullptr; ScalarType* p_f; // on grid function
  Vec ref = nullptr; ScalarType* p_ref; // on grid function
  Vec fout = nullptr; ScalarType *p_fout;  // output function values 
  Interp3_Plan_GPU* interp_plan = nullptr;
  GhostPlan* ghost_plan = nullptr;
  ScalarType* p_fghost = nullptr;  // ghost padded data
  
  int device;
  char pciBusId[512];
  cudaGetDevice ( &device );
  cudaDeviceGetPCIBusId ( pciBusId, 512, device );

  for (int i = 0; i < 3; ++i) {
    nx[i]     = m_Opt->m_Domain.nx[i];
    isize[i]  = m_Opt->m_Domain.isize[i];
    istart[i] = m_Opt->m_Domain.istart[i];
    hx[i]     = m_Opt->m_Domain.hx[i];
    nl       *= isize[i];
    ng       *= nx[i];
  }
  
  if (m_Opt->m_Verbosity > 2) {
    printVector(nx, 3, "nx");
    printVector(isize, 3, "isize");
    printVector(istart, 3, "istart");
    printVector(hx, 3, "hx");
    printScalar(nl, "nl");
    printScalar(ng, "ng");
  }
  
  ierr = VecCreate(q, 3*nl, 3*ng); CHKERRQ(ierr);
  ierr = VecCreate(xq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(yq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(zq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(f, nl, ng);     CHKERRQ(ierr);
  ierr = VecCreate(ref, nl, ng);     CHKERRQ(ierr);
  ierr = VecCreate(fout, nl, ng);  CHKERRQ(ierr);

  // create ghost plan
  AllocateOnce(ghost_plan, m_Opt, m_Opt->m_PDESolver.iporder);
  size_t g_alloc_max = ghost_plan->get_ghost_local_size_x(isize_g, istart_g);
  cudaMalloc((void**)&p_fghost, g_alloc_max);
  nlghost = isize_g[0]*isize_g[1]*isize_g[2];
    
  if (m_Opt->m_Verbosity > 2) {
    printScalar(g_alloc_max, "g_alloc_max");
    printVector(isize_g, 3, "isize_g");
    printVector(istart_g, 3, "istart_g");
    printScalar(nlghost, "nlghost");
  }
  
  cudaTextureObject_t tex = gpuInitEmptyTexture(isize_g);   

  int data_dofs[2] = {1,3};
  
  // create interpolation plan
  AllocateOnce(interp_plan, g_alloc_max, true); 
  interp_plan->allocate(nl, data_dofs, 2, nghost, isize_g);  

  ZeitGeist_define(overall);
  ZeitGeist_define(interp_overall);
  ZeitGeist_define(scatter_overall);
  ZeitGeist_define(ghost_overall);
  ZeitGeist_define(interp_kernel);
  
  std::string flag = "state";
  
  // initialize the grid points
  ierr = VecCUDAGetArray(xq, &p_xq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(f, &p_f); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(ref, &p_ref); CHKERRQ(ierr);
  initializeGrid(p_xq, p_yq, p_zq, p_f, p_ref, hx, isize, istart, nx, 0);
  ierr = VecCUDARestoreArray(ref, &p_ref); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(f,  &p_f); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(xq, &p_xq); CHKERRQ(ierr);
  
  ierr = VecCUDAGetArray(xq, &p_xq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(zq, &p_zq); CHKERRQ(ierr);
  interp_plan->scatter(nx, isize, istart, nl, nghost, p_xq, p_yq, p_zq, m_Opt->m_CartGridDims, m_Opt->m_Domain.mpicomm, timers, flag);
  ierr = VecCUDARestoreArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(xq, &p_xq); CHKERRQ(ierr);
  if (m_Opt->m_Verbosity > 2)
    reg::DbgMsgCall("Scatter done");
  
  // warmup runs
  MPI_Barrier(MPI_COMM_WORLD);
  for (int rep=0; rep<10; rep++) { 
    ierr = VecCUDAGetArray(f, &p_f);  CHKERRQ(ierr);
    ghost_plan->share_ghost_x(p_f, p_fghost);
    ierr = VecCUDARestoreArray(f, &p_f); CHKERRQ(ierr);

    reg::Assert(p_fghost != nullptr, "nullptr");
      
    if (m_Opt->m_Verbosity > 2)
      reg::DbgMsgCall("Ghost shared");
      
    ierr = VecCUDAGetArray(fout, &p_fout); CHKERRQ(ierr);
    ScalarType *wout[3] = {&p_fout[0]};
    interp_plan->interpolate( p_fghost, 
                              isize_g, 
                              nlghost,
                              nl, 
                              wout,
                              m_Opt->m_Domain.mpicomm, 
                              m_tmpInterpol1, 
                              m_tmpInterpol2, 
                              tex, 
                              m_Opt->m_PDESolver.iporder, 
                              &m_GPUtime, 0, flag);
    ierr = VecCUDARestoreArray(fout, &p_fout); CHKERRQ(ierr);
    if (rep == 0) {
      // compute error in interpolation
      ierr = VecGetArray(ref, &p_ref);  CHKERRQ(ierr);
      ierr = VecGetArray(fout, &p_fout);  CHKERRQ(ierr);
      TestErrorMultiGPU(p_ref, p_fout, nl, &error, &max);
      ierr = VecRestoreArray(fout, &p_fout); CHKERRQ(ierr);
      ierr = VecRestoreArray(ref, &p_ref); CHKERRQ(ierr);
    }
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  
  int repcount = 100;
  ZeitGeist_define(INTERPOL);
  // performance runs
  for (int rep=0; rep<20; rep++) { 
    ZeitGeist_tock(INTERPOL);
    ierr = VecCUDAGetArray(f, &p_f);  CHKERRQ(ierr);
    ghost_plan->share_ghost_x(p_f, p_fghost);
    ierr = VecCUDARestoreArray(f, &p_f); CHKERRQ(ierr);

    reg::Assert(p_fghost != nullptr, "nullptr");
      
    if (m_Opt->m_Verbosity > 2)
      reg::DbgMsgCall("Ghost shared");
      
    ierr = VecCUDAGetArray(fout, &p_fout); CHKERRQ(ierr);
    ScalarType * wout[3] = {&p_fout[0]};
    interp_plan->interpolate( p_fghost, 
                              isize_g, 
                              nlghost,
                              nl, 
                              wout,
                              m_Opt->m_Domain.mpicomm, 
                              m_tmpInterpol1, 
                              m_tmpInterpol2, 
                              tex, 
                              m_Opt->m_PDESolver.iporder, 
                              &m_GPUtime, 0, flag);
    ierr = VecCUDARestoreArray(fout, &p_fout); CHKERRQ(ierr);
    ZeitGeist_tock(INTERPOL);
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  
  // get global runtimes and errors
  std::ostringstream ss;
  MPI_Reduce(&error, &global_error, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
  global_error = sqrt(global_error/static_cast<double>(ng));
  ss << "INTERP ERROR = " << global_error << "\n" << "MAX VALUE = " << global_max << "\n";
  printOnMaster(ss.str());
  ss.str("");
  
  bool write_to_file = false;
  
  if (write_to_file) {
    std::ofstream myfile;
    if (procid == 0) {
      std::string filename = "/home/04716/naveen15/claire-dev/scripts/weak_scaling/multi_node/interp_p" + std::to_string(nprocs) + "_runtimes.csv";
      myfile.open(filename);
      myfile << "task,count,time" << std::endl;
    }

    double global_runtime, local_runtime;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Timings\n");
    for (auto zg : ZeitGeist::zgMap()) {
        char txt[120];
        local_runtime = zg.second.Total_s();
        MPI_Reduce(&local_runtime, &global_runtime, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD); 
        PetscPrintf(PETSC_COMM_WORLD, "%16s: %5lix, %0.4e s\n", zg.first.c_str(), zg.second.Count()/repcount, global_runtime/repcount);
        PetscPrintf(PETSC_COMM_WORLD, "============================================================================\n");
        if (procid == 0) 
            myfile << zg.first.c_str() << "," << zg.second.Count()/repcount << "," << global_runtime/repcount << std::endl;
    }

    if (procid == 0) 
      myfile.close();
  }

  cudaDestroyTextureObject(tex);
  ierr = VecDestroy(&q); CHKERRQ(ierr); q = nullptr;
  ierr = VecDestroy(&f); CHKERRQ(ierr); f = nullptr;
  ierr = VecDestroy(&ref); CHKERRQ(ierr); ref = nullptr;
  ierr = VecDestroy(&fout); CHKERRQ(ierr); fout = nullptr;
  
  //Free(interp_plan); 
  Free(ghost_plan);
  if (p_fghost != nullptr) {
    cudaFree(p_fghost);
    p_fghost = nullptr;
  }
  
  PetscFunctionReturn(ierr);
}

#if 0
PetscErrorCode TestVectorFieldInterpolationMultiGPU(RegOpt *m_Opt) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  printOnMaster("starting vector field interpolation unit test\n");
  
  srand(time(0));
  
  int nx[3], istart[3], isize[3];
  ScalarType hx[3], scale;
  double error[3]={0,0,0}, max[3]={0,0,0};
  ScalarType m_GPUtime;
  int nghost = 3;
  int isize_g[3], istart_g[3];
  int nlghost = 1;
  double timers[4] = {0,0,0,0};
  double global_error[3], global_max[3];
  ScalarType* m_tmpInterpol1 = nullptr, *m_tmpInterpol2 = nullptr;
  double global_timers[4] = {0, 0, 0, 0};
  int nprocs, procid;
  int nl=1, ng=1;
    
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  Vec q = nullptr; ScalarType* p_q; // query points
  Vec xq = nullptr; ScalarType* p_xq; // query points
  Vec yq = nullptr; ScalarType* p_yq; // query points
  Vec zq = nullptr; ScalarType* p_zq; // query points
  reg::VecField* f = nullptr; ScalarType* p_f[3]={nullptr, nullptr, nullptr}; // on grid function
  reg::VecField* ref = nullptr; ScalarType* p_ref[3]={nullptr, nullptr, nullptr}; // on grid function
  Vec fout = nullptr; ScalarType* p_fout = nullptr;  // output function values 
  Interp3_Plan_GPU* interp_plan = nullptr;
  GhostPlan* ghost_plan = nullptr;
  ScalarType* p_fghost = nullptr;  // ghost padded data
  
  int device;
  char pciBusId[512];
  cudaGetDevice ( &device );
  cudaDeviceGetPCIBusId ( pciBusId, 512, device );

  for (int i = 0; i < 3; ++i) {
    nx[i]     = m_Opt->m_Domain.nx[i];
    isize[i]  = m_Opt->m_Domain.isize[i];
    istart[i] = m_Opt->m_Domain.istart[i];
    hx[i]     = m_Opt->m_Domain.hx[i];
    nl       *= isize[i];
    ng       *= nx[i];
  }
  
  if (m_Opt->m_Verbosity > 2) {
    printVector(nx, 3, "nx");
    printVector(isize, 3, "isize");
    printVector(istart, 3, "istart");
    printVector(hx, 3, "hx");
    printScalar(nl, "nl");
    printScalar(ng, "ng");
  }
  
  ierr = VecCreate(q, 3*nl, 3*ng); CHKERRQ(ierr);
  ierr = VecCreate(xq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(yq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(zq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(fout, 3*nl, 3*ng); CHKERRQ(ierr);
  ierr = AllocateOnce(f, m_Opt); CHKERRQ(ierr);
  ierr = AllocateOnce(ref, m_Opt); CHKERRQ(ierr);
  ierr = VecSet(fout, 0); CHKERRQ(ierr);

  // initialize the grid points
  ierr = VecCUDAGetArray(xq, &p_xq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = f->GetArrays(p_f); CHKERRQ(ierr);
  ierr = ref->GetArrays(p_ref); CHKERRQ(ierr);
  initializeGrid(p_xq, p_yq, p_zq, p_f[0], p_ref[0], hx, isize, istart, nx, 0);
  initializeGrid(p_xq, p_yq, p_zq, p_f[1], p_ref[1], hx, isize, istart, nx, 1);
  initializeGrid(p_xq, p_yq, p_zq, p_f[2], p_ref[2], hx, isize, istart, nx, 2);
  ierr = ref->RestoreArrays(p_ref); CHKERRQ(ierr);
  ierr = f->RestoreArrays(p_f); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(xq, &p_xq); CHKERRQ(ierr);
    
  // get ghost dimensions
  if (ghost_plan == nullptr) {
    try {
      ghost_plan = new GhostPlan(m_Opt, m_Opt->m_PDESolver.iporder);
    } catch (std::bad_alloc&) {
      ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
  }
  size_t g_alloc_max = ghost_plan->get_ghost_local_size_x(isize_g, istart_g);

  cudaMalloc((void**)&p_fghost, 3*g_alloc_max);

  for (int i=0; i<3; ++i) {
      nlghost *= isize_g[i];
  }
    
  if (m_Opt->m_Verbosity > 1) {
    printScalar(g_alloc_max, "g_alloc_max");
    printVector(isize_g, 3, "isize_g");
    printVector(istart_g, 3, "istart_g");
    printScalar(nlghost, "nlghost");
  }
  
  cudaTextureObject_t tex = gpuInitEmptyTexture(isize_g);   

  int data_dofs[2] = {1,3};
  
  if (interp_plan == nullptr) {
    try { 
      interp_plan = new Interp3_Plan_GPU(g_alloc_max, true); 
    }
    catch (std::bad_alloc&) {
      ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    interp_plan->allocate(nl, data_dofs, 2, nghost, isize_g);
  }
  
  ZeitGeist_define(overall);
  ZeitGeist_define(interp_overall);
  ZeitGeist_define(scatter_overall);
  ZeitGeist_define(ghost_overall);
  
  MPI_Barrier(PETSC_COMM_WORLD);
  ZeitGeist_tick(overall);
  
  ZeitGeist_tick(scatter_overall);
  ierr = VecCUDAGetArray(xq, &p_xq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(zq, &p_zq); CHKERRQ(ierr);
  interp_plan->scatter(nx, isize, istart, nl, nghost, p_xq, p_yq, p_zq, m_Opt->m_CartGridDims, m_Opt->m_Domain.mpicomm, timers, "state");
  ierr = VecCUDARestoreArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(xq, &p_xq); CHKERRQ(ierr);
  ZeitGeist_tock(scatter_overall);

  ZeitGeist_tick(ghost_overall);
  ierr = f->GetArrays(p_f);  CHKERRQ(ierr);
  ghost_plan->share_ghost_x(p_f[0], p_fghost);
  ghost_plan->share_ghost_x(p_f[1], &p_fghost[nlghost]);
  ghost_plan->share_ghost_x(p_f[2], &p_fghost[2*nlghost]);
  ierr = f->RestoreArrays(p_f); CHKERRQ(ierr);
  ZeitGeist_tock(ghost_overall);

  reg::Assert(p_fghost != nullptr, "nullptr");
      
  ZeitGeist_tick(interp_overall);
  ierr = VecCUDAGetArray(fout, &p_fout); CHKERRQ(ierr);
  ScalarType *wout[3] = {&p_fout[0], &p_fout[1*nl], &p_fout[2*nl]};
  interp_plan->interpolate( p_fghost, 
                            isize_g, 
                            nlghost,
                            nl, 
                            wout,
                            m_Opt->m_Domain.mpicomm, 
                            m_tmpInterpol1, 
                            m_tmpInterpol2, 
                            tex, 
                            m_Opt->m_PDESolver.iporder, 
                            &m_GPUtime, 1, "state");
  ierr = VecCUDARestoreArray(fout, &p_fout); CHKERRQ(ierr);
  ZeitGeist_tock(interp_overall);
  ZeitGeist_tock(overall);
  MPI_Barrier(PETSC_COMM_WORLD);

  // compute error in interpolation
  int k = 0;
  ierr = VecGetArray(fout, &p_fout);  CHKERRQ(ierr);
  ierr = VecGetArray(ref->m_X1, &p_ref[k]);  CHKERRQ(ierr);
  TestError(p_ref[k], &p_fout[0], nl, &error[k], &max[k]);
  ierr = VecRestoreArray(ref->m_X1,  &p_ref[k]); CHKERRQ(ierr);
  
  k = 1;
  ierr = VecGetArray(ref->m_X2, &p_ref[k]);  CHKERRQ(ierr);
  TestError(p_ref[k], &p_fout[nl], nl, &error[k], &max[k]);
  ierr = VecRestoreArray(ref->m_X2,  &p_ref[k]); CHKERRQ(ierr);

  k = 2;
  ierr = VecGetArray(ref->m_X3, &p_ref[k]);  CHKERRQ(ierr);
  TestError(p_ref[k], &p_fout[2*nl], nl, &error[k], &max[k]);
  ierr = VecRestoreArray(ref->m_X3,  &p_ref[k]); CHKERRQ(ierr);
  ierr = VecRestoreArray(fout, &p_fout); CHKERRQ(ierr);
  
  // get global runtimes and errors
  std::ostringstream ss;
  MPI_Reduce((const void*)error, (void*)global_error, 3, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce((const void*)max, (void*)global_max, 3, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
  for (int i=0; i<3; i++) {
    global_error[i] = sqrt(global_error[i]/static_cast<double>(ng));
    ss << "INTERP ERROR[" << i << "] = " << global_error[i] << "\n" << "MAX VALUE[" << i << "] = " << global_max[i] << "\n";
    printOnMaster(ss.str());
    ss.str("");
  }
  
  cudaDestroyTextureObject(tex);
  ierr = VecDestroy(&q); CHKERRQ(ierr); q = nullptr;
  ierr = Free(f); CHKERRQ(ierr); f = nullptr;
  ierr = Free(ref); CHKERRQ(ierr); ref = nullptr;
  ierr = VecDestroy(&fout); CHKERRQ(ierr); fout = nullptr;
  
  
  if (interp_plan != nullptr) {
    delete interp_plan;
    interp_plan = nullptr;
  }

  if (ghost_plan != nullptr) {
    delete ghost_plan;
    ghost_plan = nullptr;
  }

  if (p_fghost != nullptr) {
#if defined(REG_HAS_MPICUDA)
    cudaFree(p_fghost);
#else
    accfft_free(p_fghost); 
#endif
    p_fghost = nullptr;
  }

  PetscFunctionReturn(ierr);
}
#endif


} // namespace UnitTest


} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

