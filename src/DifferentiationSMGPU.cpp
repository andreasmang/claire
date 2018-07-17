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

#ifndef _DIFFERENTIATIONSMGPU_CPP_
#define _DIFFERENTIATIONSMGPU_CPP_




#ifdef REG_FFT_CUDA




#include "DifferentiationSMGPU.hpp"
#include "accfft_gpu.h"
#include "accfft_gpuf.h"
#include "accfft_operators_gpu.h"


void deviceMemory() {
  size_t free_byte ;
  size_t total_byte ;

  cudaMemGetInfo( &free_byte, &total_byte );

  double free_db = (double)free_byte ;
  double total_db = (double)total_byte ;
  double used_db = total_db - free_db ;

  printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
    used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
}


namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DifferentiationSM::DifferentiationSM() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DifferentiationSM::DifferentiationSM(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DifferentiationSM::~DifferentiationSM() {
    printf("DifferentiationSMGPU: %i grad, %i div\n",c_grad,c_div);
    this->ClearMemory();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode DifferentiationSM::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    std::cout << "Diff_GPU create plan\n";
    deviceMemory();
    xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;
    c_grad = 0;
    c_div = 0;
    int nx[3];
    nx[0] = m_Opt->m_Domain.nx[0];
    nx[1] = m_Opt->m_Domain.nx[1];
    nx[2] = m_Opt->m_Domain.nx[2];
    cudaMalloc((void**)&u[0],  sizeof(ScalarType)*this->m_Opt->m_Domain.nl);
    cudaMalloc((void**)&u[1],  sizeof(ScalarType)*this->m_Opt->m_Domain.nl);
    cudaMalloc((void**)&u[2],  sizeof(ScalarType)*this->m_Opt->m_Domain.nl);
    cudaMalloc((void**)&u[3],  sizeof(ScalarType)*this->m_Opt->m_Domain.nl);
    plan = accfft_plan_dft_3d_r2c_gpuf(nx, nullptr, nullptr,
                                       m_Opt->m_FFT.mpicomm, ACCFFT_MEASURE);

    std::cout << "Diff_GPU plan created\n";
    deviceMemory();
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DifferentiationSM::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    accfft_destroy_plan(plan);
    accfft_cleanup_gpuf();
    cudaFree(u[0]);
    cudaFree(u[1]);
    cudaFree(u[2]);
    cudaFree(u[3]);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute gradient of a scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Gradient(ScalarType *g1,
                                           ScalarType *g2,
                                           ScalarType *g3,
                                           ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    //accfft_grad_gpu_t(g1, g2, g3, m, plan, &xyz, timer);
    cudaMemcpy(u[0], m, sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyHostToDevice);
    accfft_grad_gpu_t(u[1], u[2], u[3], u[0], this->plan, &xyz, timer);
    cudaMemcpy(g1, u[1], sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost);
    cudaMemcpy(g2, u[2], sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost);
    cudaMemcpy(g3, u[3], sizeof(ScalarType)*this->m_Opt->m_Domain.nl, cudaMemcpyDeviceToHost);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTGRAD);
    
    this->m_Opt->IncreaseFFTTimers(timer);
    c_grad++;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply laplacian operator to scalar field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Laplacian(ScalarType *l,
                                            ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_laplace_t(l, m, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply laplacian operator to vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Divergence(ScalarType *l,
                                            ScalarType *v1,
                                            ScalarType *v2,
                                            ScalarType *v3) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_divergence_t(l, v1, v2, v3, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    this->m_Opt->IncrementCounter(FFT, FFTDIV);

    c_div++;
    
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply biharmonic operator to vector field
 *******************************************************************/
PetscErrorCode DifferentiationSM::Biharmonic(ScalarType *b,
                                            ScalarType *m) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    for (int i=0; i<NFFTTIMERS; ++i) timer[i] = 0;
    
    this->m_Opt->StartTimer(FFTSELFEXEC);
    accfft_biharmonic_t(b, m, this->m_Opt->m_FFT.plan, timer);
    this->m_Opt->StopTimer(FFTSELFEXEC);
    
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(ierr);
}




}  // end of name space




#endif  // REG_HAS_CUDA




#endif  // _DIFFERENTIATIONSMGPU_CPP_
