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

#ifndef _LOGGER_CPP_
#define _LOGGER_CPP_

#include "RegUtils.hpp"


// flags for timers
enum TimerType {
    T2SEXEC = 0,  ///< time to solution (execution time)
    PDEEXEC,      ///< pde solves (execution time)
    HMVEXEC,      ///< hessian mat vec (execution time)
    PMVSETUP,     ///< setup time for preconditioner
    PMVEXEC,      ///< precond mat vec (execution time)
    GRADEXEC,     ///< gradient evaluation (execution time)
    OBJEXEC,      ///< objective evluation (execution time)
    FFTSETUP,     ///< fft setup time
    FFTSELFEXEC,  ///< fft execution time
    IPSELFEXEC,   ///< execution time for interpolation
    NTIMERS,      ///< to allocate the timers
};



enum FFTTimers {
    FFTTRANSPOSE = 0,  ///< contains all communication time
    FFTSHUFFLE   = 1,  ///< contained in FFTTRANSPOSE
    FFTCOMM      = 2,  ///< contained in FFTTRANSPOSE
    FFTRESHUFFLE = 3,  ///< contained in FFTTRANSPOSE
    FFTEXECUTE   = 4,  ///< execution time
    FFTNAN       = 5,  ///< currently empty
    FFTHADAMARD  = 6,  ///< hadamard products for differential operators
    NFFTTIMERS   = 7,  ///< total number of itmers
};


enum FFTCounters {
    FFTGRAD = 4,  ///< number of ffts for gradient opteration
    FFTDIV  = 4,  ///< number of ffts for divergence operation
};



// counters (number of operations)
enum CounterType {
    PDESOLVE = 0,  ///< PDE solves
    HESSMATVEC,    ///< # hessian matvecs
    PCMATVEC,      ///< preconditioner matvecs
    OBJEVAL,       ///< objective evaluations
    GRADEVAL,      ///< gradient evaluations
    IPVEC,         ///< interpolation execution time
    IP,            ///< interpolation execution time
    FFT,           ///< fft evaluations
    ITERATIONS,    ///< number of outer iterations
    NCOUNTERS,     ///< to allocate the counters
};


/*! flags for timers */
enum LogType{
    LOGRES,
    LOGCONV,
    LOGKSPRES,
    LOGJAC,
    LOGLOAD,
    NLOGFLAGS
};


/*! parameter for grid continuation */
struct Log {
    enum TimerValue {LOG = 0, MIN, MAX, AVG, NVALTYPES};
    std::vector<ScalarType> distance;        ///< convergence for residual
    std::vector<ScalarType> regularization;  ///< convergence for regularization
    std::vector<ScalarType> objective;       ///< convergence for objective
    std::vector<int> outeriterations;        ///< iterations of solver
    std::vector<ScalarType> kspresidual;     ///< residual of krylov method
    std::vector<int> kspiterations;          ///< iterations of krylov method
    ScalarType finalresidual[4];
    bool enabled[NLOGFLAGS];

    bool memoryusage;
    double timer[NTIMERS][NVALTYPES];
    double temptimer[NTIMERS];
    bool timerruns[NTIMERS];
    unsigned int counter[NCOUNTERS];
    double ffttimers[NFFTTIMERS][NVALTYPES];
    double iptimers[4][NVALTYPES];
};


class Logger {
 public:
    typedef Logger Self;

}




#endif


