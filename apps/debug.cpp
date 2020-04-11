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
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#include "CLAIREUtils.hpp"
#include "RegOpt.hpp"

/********************************************************************
 * @brief main function to run registration
 *******************************************************************/
int main(int argc, char **argv) {
    PetscErrorCode ierr = 0;

    // initialize petsc (user is not allowed to set petsc options)
    ierr = PetscInitialize(0, reinterpret_cast<char***>(NULL),
                              reinterpret_cast<char*>(NULL),
                              reinterpret_cast<char*>(NULL)); CHKERRQ(ierr);
    PetscFunctionBegin;
    
    reg::RegOpt opt;
    
    opt.m_Domain.mpicomm = PETSC_COMM_WORLD;
    opt.m_FFT.fft = nullptr;
    opt.m_FFT.nx[0] = 256;
    opt.m_FFT.nx[1] = 256;
    opt.m_FFT.nx[2] = 256;
    
    ierr = AllocateOnce(opt.m_FFT.fft, &opt, &opt.m_FFT);
    
    opt.m_FFT.fft->InitFFT();
    
    ComplexType *spectral;
    ScalarType *real;
    cudaMalloc(&spectral, opt.m_FFT.nalloc);
    cudaMalloc(&real, opt.m_FFT.nalloc);
    
    MPI_Barrier(PETSC_COMM_WORLD);
    for (int r=0; r<20; ++r) {
      opt.m_FFT.fft->FFT_R2C(real, spectral);
      opt.m_FFT.fft->FFT_C2R(spectral, real);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    
    for (int r=0; r<20; ++r) {
      ZeitGeist_define(FFT_LIB);
      ZeitGeist_tick(FFT_LIB);
      opt.m_FFT.fft->FFT_R2C(real, spectral);
      opt.m_FFT.fft->Scale(spectral, 1./static_cast<ScalarType>(opt.m_FFT.nx[0]*opt.m_FFT.nx[1]*opt.m_FFT.nx[2]));
      opt.m_FFT.fft->FFT_C2R(spectral, real);
      ZeitGeist_tock(FFT_LIB);
    }
    
#ifdef ZEITGEIST
    reg::Msg("-----------------------------------------------------------------------------------------------------");
    reg::Msg("ZeitGeist:");
    for (auto zg : ZeitGeist::zgMap()) {
      char txt[120];
      double global_runtime;
      double local_runtime = zg.second.Total_s();
      MPI_Reduce(&local_runtime, &global_runtime, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
      sprintf(txt, "  %16s: %5lix, %0.10lf",zg.first.c_str(), zg.second.Count(), global_runtime);
      reg::Msg(txt);
    }
    reg::Msg("-----------------------------------------------------------------------------------------------------");
#endif

    ierr = reg::Finalize(); CHKERRQ(ierr);

    return 0;
}


