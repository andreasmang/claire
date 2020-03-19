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
#if defined(REG_HAS_MPI_CUDA) || defined(REG_HAS_CUDA)
#include "interp3_gpu_mpi.hpp"
#endif
#include "zeitgeist.hpp"
#define CONSTANT

void TestFunction(ScalarType &val, const ScalarType x1, const ScalarType x2, const ScalarType x3) {
  val =  ( PetscSinReal(8*x1)*PetscSinReal(8*x1)
        + PetscSinReal(2*x2)*PetscSinReal(2*x2)
        + PetscSinReal(4*x3)*PetscSinReal(4*x3) )/3.0;
}

inline void printOnMaster(std::string msg) {
  PetscPrintf(PETSC_COMM_WORLD, msg.c_str());
}



void printGPUMemory() {
    size_t free, used;
    cudaMemGetInfo(&free, &used);
    used -= free;
    std::string msg = "Used mem = " + std::to_string(used/1E9) + " GB, Free mem = " + std::to_string(free/1E9) + " GB\n";
    PetscPrintf(PETSC_COMM_WORLD, msg.c_str());
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
  for (int i = 0; i < nl; ++i) {
    double local = static_cast<double>(ref[i]) - static_cast<double>(eval[i]);
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
  
  int nx[3], istart[3], isize[3], nl;
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
        
      // x1 += hx[0]*0.5;
      // x2 += hx[1]*0.5;
      // x3 += hx[2]*0.5;
        
        x1 +=0* hx[0];
        x2 +=0* hx[1];
        x3 +=0* hx[2];
        
        if (x1 > 2.*M_PI) x1 -= 2.*M_PI;
        if (x2 > 2.*M_PI) x2 -= 2.*M_PI;
        if (x3 > 2.*M_PI) x3 -= 2.*M_PI;
    
        //x1 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
        //x2 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
        //x3 = static_cast<double>(rand())*2.*M_PI/static_cast<double>(RAND_MAX);
        
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
  float timer = 0;
  for (int i=0; i<10; ++i) {
    gpuInterp3D(pg, pq1, pq2, pq3, pe, tmp1, tmp2, nx, nl, tex, m_Opt->m_PDESolver.iporder, &timer);
    cudaDeviceSynchronize();
  }
  
  double error, max;
  error = 0;
  max = 0;
  
  ZeitGeist_define(INTERPOL);
  for (int i=0; i<100; ++i) {
    double err, m;
    ZeitGeist_tick(INTERPOL);
    gpuInterp3D(pg, pq1, pq2, pq3, pe, tmp1, tmp2, nx, nl, tex, m_Opt->m_PDESolver.iporder, &timer);
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
  
  for (int i=0; i<10; ++i) {
    gpuInterp3D(pg, pq1, pq2, pq3, pe, tmp1, tmp2, nx, nl, tex, m_Opt->m_PDESolver.iporder, &timer);
    cudaDeviceSynchronize();
  }
  
  error = 0;
  max = 0;
  
  ZeitGeist_define(INTERPOL_RAND);
  for (int i=0; i<100; ++i) {
    double err, m;
    ZeitGeist_tick(INTERPOL_RAND);
    gpuInterp3D(pg, pq1, pq2, pq3, pe, tmp1, tmp2, nx, nl, tex, m_Opt->m_PDESolver.iporder, &timer);
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
  
  int nx[3], istart[3], isize[3];
  int nl = 1, ng = 1;
  ScalarType hx[3], scale;
  double error=0, max=0;
  ScalarType m_GPUtime;
  int nghost = 3;
  int isize_g[3], istart_g[3];
  int nlghost = 1;
  double timers[4] = {0,0,0,0};
  double global_error, global_max;
  ScalarType* m_tmpInterpol1 = NULL, *m_tmpInterpol2 = NULL;
  double global_timers[4] = {0, 0, 0, 0};
  int nprocs, procid;
  PetscRandom rctx;
    
  MPI_Comm_rank(m_Opt->m_FFT.mpicomm, &procid);
  MPI_Comm_size(m_Opt->m_FFT.mpicomm, &nprocs);

  Vec q = NULL; ScalarType* p_q; // query points
  Vec xq = NULL; ScalarType* p_xq; // query points
  Vec yq = NULL; ScalarType* p_yq; // query points
  Vec zq = NULL; ScalarType* p_zq; // query points
  Vec f = NULL; ScalarType* p_f; // on grid function
  Vec ref = NULL; ScalarType* p_ref; // on grid function
  Vec fout = NULL; ScalarType *p_fout;  // output function values 
  Vec dhx = NULL; ScalarType *p_dhx;  // output function values 
  Interp3_Plan_GPU* interp_plan = NULL;
  ScalarType* p_fghost = NULL;  // ghost padded data
  
  int device;
  char pciBusId[512];
  cudaGetDevice ( &device );
  cudaDeviceGetPCIBusId ( pciBusId, 512, device );
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d][%d] has PCI-id = %s\n", procid, device, pciBusId);
  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

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
  //ierr = VecCreate(dhx, nl/nl, ng/ng); CHKERRQ(ierr);
  ierr = VecCreate(xq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(yq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(zq, nl, ng); CHKERRQ(ierr);
  ierr = VecCreate(f, nl, ng);     CHKERRQ(ierr);
  ierr = VecCreate(ref, nl, ng);     CHKERRQ(ierr);
  ierr = VecCreate(fout, nl, ng);  CHKERRQ(ierr);

  //PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
  //VecSetRandom(dhx,rctx);
  //PetscRandomDestroy(&rctx);
  
  // initialize the grid points
  //ierr = VecCUDAGetArray(dhx, &p_dhx); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(xq, &p_xq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(f, &p_f); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(ref, &p_ref); CHKERRQ(ierr);
  initializeGrid(p_xq, p_yq, p_zq, p_f, p_ref, hx, isize, istart, nx, p_dhx);
  ierr = VecCUDARestoreArray(ref, &p_ref); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(f,  &p_f); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(xq, &p_xq); CHKERRQ(ierr);
  //ierr = VecCUDARestoreArray(dhx, &p_dhx); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "grid initialized\n");

  
  // no need of dhx now
  //if (dhx != NULL) { ierr = VecDestroy(&dhx); CHKERRQ(ierr); dhx = NULL; }
 
    
  // get ghost dimensions
  size_t g_alloc_max = accfft_ghost_xyz_local_size_dft_r2c(m_Opt->m_FFT.fft->m_plan, nghost, isize_g, istart_g);
  for (int i=0; i<3; ++i) {
      nlghost *= isize_g[i];
  }
    
  if (m_Opt->m_Verbosity > 2) {
    printScalar(g_alloc_max, "g_alloc_max");
    printVector(isize_g, 3, "isize_g");
    printVector(istart_g, 3, "istart_g");
    printScalar(nlghost, "nlghost");
  }

  cudaTextureObject_t tex = gpuInitEmptyTexture(isize_g);   
  
  if (interp_plan == NULL) {
    try { 
      interp_plan = new Interp3_Plan_GPU(g_alloc_max); 
    }
    catch (std::bad_alloc&) {
      ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    interp_plan->allocate(nl,1);
  }

  ierr = VecCUDAGetArray(xq, &p_xq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDAGetArray(zq, &p_zq); CHKERRQ(ierr);
  interp_plan->scatter(1, nx, isize, istart, nl, nghost, p_xq, p_yq, p_zq, m_Opt->m_CartGridDims, m_Opt->m_FFT.mpicomm, timers);
  ierr = VecCUDARestoreArray(zq, &p_zq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(yq, &p_yq); CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(xq, &p_xq); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "scatter done\n");

#if defined(REG_HAS_MPICUDA)
  cudaMalloc((void**)&p_fghost, g_alloc_max);
#else
  p_fghost = reinterpret_cast<ScalarType*>(accfft_alloc(g_alloc_max));
#endif
  // get ghost padding 
  ierr = VecCUDAGetArray(f, &p_f);  CHKERRQ(ierr);
  accfft_get_ghost_xyz(m_Opt->m_FFT.fft->m_plan, nghost, isize_g, p_f, p_fghost, timers); 
  ierr = VecCUDARestoreArray(f, &p_f); CHKERRQ(ierr);

  ierr = VecCUDAGetArray(fout, &p_fout); CHKERRQ(ierr);
  interp_plan->interpolate( p_fghost, 
                            1, 
                            nx,
                            isize,
                            istart,
                            isize_g, 
                            nlghost,
                            nl, 
                            nghost, 
                            p_fout,
                            m_Opt->m_CartGridDims,
                            m_Opt->m_FFT.mpicomm, 
                            timers, 
                            m_tmpInterpol1, 
                            m_tmpInterpol2, 
                            tex, 
                            3, 
                            &m_GPUtime);
  ierr = VecCUDARestoreArray(fout, &p_fout); CHKERRQ(ierr);
  
  printGPUMemory();

  // compute error in interpolation
  ierr = VecGetArray(fout, &p_fout);  CHKERRQ(ierr);
  ierr = VecGetArray(ref, &p_ref);  CHKERRQ(ierr);
  TestError(p_ref, p_fout, nl, &error, &max);
  ierr = VecRestoreArray(ref, &p_ref); CHKERRQ(ierr);
  ierr = VecRestoreArray(fout, &p_fout); CHKERRQ(ierr);
  
  // get global runtimes and errors
  std::ostringstream ss;
  
  MPI_Reduce(timers, global_timers, 4, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
  ss << "MAX COMM = " << std::scientific << global_timers[0] << std::scientific << "s, INTERP = " << global_timers[1] << std::scientific << "s, MEMCPY = " << global_timers[2] << std::scientific << "s, MISC = " << global_timers[3] << std::scientific << "s\n"; 
  std::string s = ss.str();
  if (procid==0) {
    std::cout << s << std::scientific << std::endl;
  }
  printOnMaster(s);
  ss.str("");
  
  MPI_Reduce(timers, global_timers, 4, MPI_DOUBLE, MPI_MIN, 0, PETSC_COMM_WORLD);
  ss << "MIN COMM = " << global_timers[0] << "s, INTERP = " << global_timers[1] << "s, MEMCPY = " << global_timers[2] << "s, MISC = " << global_timers[3] << "s\n"; 
  printOnMaster(ss.str());
  ss.str("");
  
  MPI_Reduce(timers, global_timers, 4, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  for (int i=0; i<4; i++) global_timers[i] /= nprocs;
  ss << "AVG COMM = " << global_timers[0] << "s, INTERP = " << global_timers[1] << "s, MEMCPY = " << global_timers[2] << "s, MISC = " << global_timers[3] << "s\n"; 
  printOnMaster(ss.str());
  ss.str("");

  MPI_Reduce(&error, &global_error, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
  global_error = sqrt(global_error/static_cast<double>(ng));
  ss << "INTERP ERROR = " << global_error << "\n" << "MAX VALUE = " << global_max << "\n";
  printOnMaster(ss.str());
  ss.str("");


  cudaDestroyTextureObject(tex);
  ierr = VecDestroy(&q); CHKERRQ(ierr); q = NULL;
  ierr = VecDestroy(&f); CHKERRQ(ierr); f = NULL;
  ierr = VecDestroy(&ref); CHKERRQ(ierr); ref = NULL;
  ierr = VecDestroy(&fout); CHKERRQ(ierr); fout = NULL;
  
  
  if (interp_plan != NULL) {
    delete interp_plan;
    interp_plan = NULL;
  }

  if (p_fghost != NULL) {
#if defined(REG_HAS_MPICUDA)
    cudaFree(p_fghost);
#else
    accfft_free(p_fghost); 
#endif
    p_fghost = NULL;
  }

  PetscFunctionReturn(ierr);

}


} // namespace UnitTest


} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

