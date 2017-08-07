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


#include "Logger.hpp"


/********************************************************************
 * @brief resets all timers
 *******************************************************************/
PetscErrorCode Logger::ResetTimers() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__func__);

    for (int i = 0; i < NTIMERS; ++i) {
        if (i != FFTSETUP) {
            for (int j = 0; j < NVALTYPES; ++j) {
                this->m_Timer[i][j] = 0.0;
            }
        }
        this->m_TimerIsRunning[i] = false;
        this->m_TempTimer[i] = 0.0;
    }

    for (int i = 0; i < NFFTTIMERS; ++i) {
        for (int j = 0; j < NVALTYPES; ++j) {
            this->m_FFTTimers[i][j] = 0.0;
        }
    }
    this->m_FFTAccumTime = 0.0;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < NVALTYPES; ++j) {
            this->m_InterpTimers[i][j] = 0.0;
        }
    }
    this->m_IPAccumTime = 0.0;



    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resets timer
 *******************************************************************/
PetscErrorCode Logger::ResetTimer(TimerType id) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__func__);

    for (int j = 0; j < NVALTYPES; ++j) {
        this->m_Timer[id][j] = 0.0;
    }
    this->m_TimerIsRunning[id] = false;
    this->m_TempTimer[id] = 0.0;

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief start the timer (checks if running)
 *******************************************************************/
PetscErrorCode Logger::StartTimer(TimerType id) {
    PetscErrorCode ierr = 0;
    std::string msg;

    PetscFunctionBegin;

    this->Enter(__func__);

    msg = "fatal error: timer has already been started";
    ierr = Assert(this->m_TimerIsRunning[id] == false, msg); CHKERRQ(ierr);

    this->m_TempTimer[id] = -MPI_Wtime();
    this->m_TimerIsRunning[id] = true;

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief stop setup timer (checks if running)
 *******************************************************************/
PetscErrorCode Logger::StopTimer(TimerType id) {
    PetscErrorCode ierr = 0;
    std::string msg;

    PetscFunctionBegin;

    this->Enter(__func__);

    msg = "fatal error: timer has not been started";
    ierr = Assert(this->m_TimerIsRunning[id], msg); CHKERRQ(ierr);

    this->m_TempTimer[id] += MPI_Wtime();
    this->m_Timer[id][LOG] += this->m_TempTimer[id];

    // tell the world that we stop the timer
    this->m_TimerIsRunning[id] = false;

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief process the timers
 *******************************************************************/
PetscErrorCode Logger::ProcessTimers() {
    PetscErrorCode ierr = 0;
    int rval, rank, nproc;
    double *fftall = NULL, *interpall = NULL, *ttsall = NULL;
    double ival = 0.0, xval = 0.0, ivalsum = 0.0;

    PetscFunctionBegin;

    this->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

    for (int id = 0; id < NTIMERS; ++id) {
        // remember input value
        ival = this->m_Timer[id][LOG];

        // get maximal execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_MIN, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_Timer[id][MIN] = xval;

        // get maximal execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_Timer[id][MAX] = xval;

        // get mean execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_Timer[id][AVG] = xval;

        if (rank == 0) {
            this->m_Timer[id][AVG] /= static_cast<double>(nproc);
        }
    }


    ivalsum = 0.0;
    for (int i = 0; i < NFFTTIMERS; ++i) {
        // remember input value
        ival = this->m_FFTTimers[i][LOG];

        if ((i == FFTTRANSPOSE) || (i == FFTEXECUTE) || (i == FFTHADAMARD)){
            ivalsum += ival;
        }

        // get maximal execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_MIN, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_FFTTimers[i][MIN] = xval;

        // get maximal execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_FFTTimers[i][MAX] = xval;

        // get mean execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_FFTTimers[i][AVG] = xval;

        if (rank == 0) {
            this->m_FFTTimers[i][AVG] /= static_cast<double>(nproc);
        }
    }

    // get max of accumulated time accross all procs
    rval = MPI_Reduce(&ivalsum, &xval, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
    this->m_FFTAccumTime = xval;


    ivalsum = 0.0;
    for (int i = 0; i < 4; ++i) {
        // remember input value
        ival = this->m_InterpTimers[i][LOG];
        ivalsum += ival;

        // get maximal execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_MIN, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_InterpTimers[i][MIN] = xval;

        // get maximal execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_InterpTimers[i][MAX] = xval;

        // get mean execution time
        rval = MPI_Reduce(&ival, &xval, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
        ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
        this->m_InterpTimers[i][AVG] = xval;

        if (rank == 0) {
            this->m_InterpTimers[i][AVG] /= static_cast<double>(nproc);
        }
    }

    // get the timings that correspond to the slowest proc
    if (rank == 0) {
        try {ttsall = new double[nproc];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        try {fftall = new double[nproc];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        try {interpall = new double[nproc];}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    ival = this->m_Timer[FFTSELFEXEC][LOG];
    rval = MPI_Gather(&ival, 1, MPI_DOUBLE, fftall, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    ierr = MPIERRQ(rval); CHKERRQ(ierr);

    ival = this->m_Timer[IPSELFEXEC][LOG];
    rval = MPI_Gather(&ival, 1, MPI_DOUBLE, interpall, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    ierr = MPIERRQ(rval); CHKERRQ(ierr);

    ival = this->m_Timer[T2SEXEC][LOG];
    rval = MPI_Gather(&ival, 1, MPI_DOUBLE, ttsall, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    ierr = MPIERRQ(rval); CHKERRQ(ierr);

    this->m_TTSSlowest = 0.0;
    if (rank == 0) {
        for (int i = 0; i < nproc; ++i) {
            if (this->m_TTSSlowest < ttsall[i]) {
                this->m_IPSlowest  = interpall[i];
                this->m_FFTSlowest = fftall[i];
                this->m_TTSSlowest = ttsall[i];
                this->m_IDSlowest  = i;
            }
        }
    }


    // get max of accumulated time accross all procs
    rval = MPI_Reduce(&ivalsum, &xval, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);
    this->m_IPAccumTime = xval;

    if (ttsall != NULL) {delete [] ttsall; ttsall = NULL;}
    if (fftall != NULL) {delete [] fftall; fftall = NULL;}
    if (interpall != NULL) {delete [] interpall; interpall = NULL;}


    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resets counters
 *******************************************************************/
PetscErrorCode Logger::ResetCounters() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__func__);

    if (this->m_Verbosity > 2) {
        ierr = DbgMsg("resetting counters"); CHKERRQ(ierr);
    }
    for (int i = 0; i < NCOUNTERS; ++i) {
        this->m_Counter[i] = 0;
    }

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resets counter
 *******************************************************************/
PetscErrorCode Logger::ResetCounter(CounterType id) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__func__);

    this->m_Counter[id] = 0;

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write log results to file
 *******************************************************************/
PetscErrorCode Logger::WriteLogFile() {
    PetscErrorCode ierr = 0;

    if (this->m_Log.enabled[LOGLOAD]) {
        ierr = this->WriteWorkLoadLog(); CHKERRQ(ierr);
    }

    if (this->m_Log.enabled[LOGKSPRES]) {
        ierr = this->WriteKSPLog(); CHKERRQ(ierr);
    }

    if (this->m_Log.enabled[LOGCONV]) {
        ierr = this->WriteConvergenceLog(); CHKERRQ(ierr);
    }

    if (this->m_Log.enabled[LOGRES]) {
        ierr = this->WriteConvergenceLog(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write log file for workload
 *******************************************************************/
PetscErrorCode Logger::WriteWorkLoadLog() {
    PetscErrorCode ierr = 0;
    std::string fn, path;
    std::ofstream logwriter;
    int rank;
    PetscFunctionBegin;

    this->Enter(__func__);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    path = this->m_FileNames.xfolder;

    // write out logfile
    if (rank == 0) {
        fn = path + "registration-performance.log";
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);
        std::cout << std::endl;
        ierr = this->WriteWorkLoadLog(std::cout); CHKERRQ(ierr);
        std::cout << std::endl;
        ierr = this->WriteWorkLoadLog(logwriter); CHKERRQ(ierr);
        logwriter.close();
    }

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief write log file for workload
 *******************************************************************/
PetscErrorCode Logger::WriteWorkLoadLog(std::ostream& logwriter) {
    PetscErrorCode ierr = 0;
    std::string line, path;
    int rank, nproc, count = 0;
    PetscFunctionBegin;

    this->Enter(__func__);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

    line = std::string(this->m_LineLength, '-');
    path = this->m_FileNames.xfolder;

    // write out logfile
    if (rank == 0) {
        std::time_t result = std::time(NULL);
        logwriter << "# run finished on " << std::asctime(std::localtime(&result));
#ifdef GIT_VERSION
        logwriter << "# git version " << GIT_VERSION << std::endl;
#endif
        logwriter << "# problem size (nx1,nx2,nx3,nc,nt,nl,ng)=("
                  << this->m_Domain.nx[0] << ","
                  << this->m_Domain.nx[1] << ","
                  << this->m_Domain.nx[2] << ","
                  << this->m_Domain.nc << ","
                  << this->m_Domain.nt << ","
                  << this->m_Domain.nl << ","
                  << this->m_Domain.ng << ")"
                  << std::endl;

        logwriter << "# processors " << nproc
                  << " " << this->m_CartGridDims[0]
                  << "x" << this->m_CartGridDims[1] << std::endl;
        logwriter << "# eventname count minp maxp avgp maxp_by_count" << std::endl;

        count = 1;
        logwriter << "\"time to solution\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[T2SEXEC][MIN]
                  << " " << this->m_Timer[T2SEXEC][MAX]
                  << " " << this->m_Timer[T2SEXEC][AVG]
                  << " " << this->m_Timer[T2SEXEC][MAX]
                  << std::endl;

        ierr = Assert(this->m_Counter[PDESOLVE] > 0, "bug in counter"); CHKERRQ(ierr);
        count = this->m_Counter[PDESOLVE];
        count = count > 0 ? count : 1;
        logwriter << "\"pde solves\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[PDEEXEC][MIN]
                  << " " << this->m_Timer[PDEEXEC][MAX]
                  << " " << this->m_Timer[PDEEXEC][AVG]
                  << " " << this->m_Timer[PDEEXEC][MAX] / static_cast<double>(count)
                  << std::endl;

        count = this->m_Counter[OBJEVAL];
        count = count > 0 ? count : 1;
        logwriter << "\"objective eval\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[OBJEXEC][MIN]
                  << " " << this->m_Timer[OBJEXEC][MAX]
                  << " " << this->m_Timer[OBJEXEC][AVG]
                  << " " << this->m_Timer[OBJEXEC][MAX] / static_cast<double>(count)
                  << std::endl;

        // if time has been logged
        count = this->m_Counter[GRADEVAL];
        count = count > 0 ? count : 1;
        logwriter << "\"gradient eval\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[GRADEXEC][MIN]
                  << " " << this->m_Timer[GRADEXEC][MAX]
                  << " " << this->m_Timer[GRADEXEC][AVG]
                  << " " << this->m_Timer[GRADEXEC][MAX] / static_cast<double>(count)
                  << std::endl;

        // if time has been logged
        count = this->m_Counter[HESSMATVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"hessian matvec\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[HMVEXEC][MIN]
                  << " " << this->m_Timer[HMVEXEC][MAX]
                  << " " << this->m_Timer[HMVEXEC][AVG]
                  << " " << this->m_Timer[HMVEXEC][MAX] / static_cast<double>(count)
                  << std::endl;


        // if time has been logged
        count = this->m_Counter[PCMATVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"precond matvec (setup)\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[PMVSETUP][MIN]
                  << " " << this->m_Timer[PMVSETUP][MAX]
                  << " " << this->m_Timer[PMVSETUP][AVG]
                  << " " << this->m_Timer[PMVSETUP][MAX] / static_cast<double>(count)
                  << std::endl;


        // if time has been logged
        count = this->m_Counter[PCMATVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"precond matvec (exec)\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[PMVEXEC][MIN]
                  << " " << this->m_Timer[PMVEXEC][MAX]
                  << " " << this->m_Timer[PMVEXEC][AVG]
                  << " " << this->m_Timer[PMVEXEC][MAX] / static_cast<double>(count)
                  << std::endl;

        count = this->m_Counter[FFT];
        count = count > 0 ? count : 1;
        logwriter << "\"fft selfexec\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[FFTSELFEXEC][MIN]
                  << " " << this->m_Timer[FFTSELFEXEC][MAX]
                  << " " << this->m_Timer[FFTSELFEXEC][AVG]
                  << " " << this->m_Timer[FFTSELFEXEC][MAX] / static_cast<double>(count)
                  << std::endl;

        // if time has been logged
        count = this->m_Counter[FFT];
        count = count > 0 ? count : 1;
        logwriter << "\"fft accumulated\""
                  << " " << count << std::scientific
                  << " " << 0.0
                  << " " << this->m_FFTAccumTime
                  << " " << 0.0
                  << " " << this->m_FFTAccumTime / static_cast<double>(count)
                  << std::endl;

        count = 1;
        logwriter << "\"fft setup\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[FFTSETUP][MIN]
                  << " " << this->m_Timer[FFTSETUP][MAX]
                  << " " << this->m_Timer[FFTSETUP][AVG]
                  << " " << this->m_Timer[FFTSETUP][MAX] / static_cast<double>(count)
                  << std::endl;

        count = this->m_Counter[FFT];
        count = count > 0 ? count : 1;
        logwriter << "\"fft communication\""
                  << " " << count << std::scientific
                  << " " << this->m_FFTTimers[FFTCOMM][MIN]
                  << " " << this->m_FFTTimers[FFTCOMM][MAX]
                  << " " << this->m_FFTTimers[FFTCOMM][AVG]
                  << " " << this->m_FFTTimers[FFTCOMM][MAX] / static_cast<double>(count)
                  << std::endl;

        count = this->m_Counter[FFT];
        count = count > 0 ? count : 1;
        logwriter << "\"fft execution\""
                  << " " << count << std::scientific
                  << " " << this->m_FFTTimers[FFTEXECUTE][MIN]
                  << " " << this->m_FFTTimers[FFTEXECUTE][MAX]
                  << " " << this->m_FFTTimers[FFTEXECUTE][AVG]
                  << " " << this->m_FFTTimers[FFTEXECUTE][MAX] / static_cast<double>(count)
                  << std::endl;

        // if time has been logged
        count = this->m_Counter[IP] + 3*this->m_Counter[IPVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"interp selfexec\""
                  << " " << count << std::scientific
                  << " " << this->m_Timer[IPSELFEXEC][MIN]
                  << " " << this->m_Timer[IPSELFEXEC][MAX]
                  << " " << this->m_Timer[IPSELFEXEC][AVG]
                  << " " << this->m_Timer[IPSELFEXEC][MAX] / static_cast<double>(count)
                  << std::endl;


        count = this->m_Counter[IP] + 3*this->m_Counter[IPVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"interp accumulated\""
                  << " " << count << std::scientific
                  << " " << 0.0
                  << " " << this->m_IPAccumTime
                  << " " << 0.0
                  << " " << 0.0
                  << std::endl;


        // if time has been logged
        count = this->m_Counter[IP] + 3*this->m_Counter[IPVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"interp comm\""
                  << " " << count << std::scientific
                  << " " << this->m_InterpTimers[0][MIN]
                  << " " << this->m_InterpTimers[0][MAX]
                  << " " << this->m_InterpTimers[0][AVG]
                  << " " << this->m_InterpTimers[0][MAX] / static_cast<double>(count)
                  << std::endl;

        // if time has been logged
        count = this->m_Counter[IP] + 3*this->m_Counter[IPVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"interp exec\""
                  << " " << count << std::scientific
                  << " " << this->m_InterpTimers[1][MIN]
                  << " " << this->m_InterpTimers[1][MAX]
                  << " " << this->m_InterpTimers[1][AVG]
                  << " " << this->m_InterpTimers[1][MAX] / static_cast<double>(count)
                  << std::endl;

        // if time has been logged
        count = this->m_Counter[IP] + 3*this->m_Counter[IPVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"interp alloc\""
                  << " " << count << std::scientific
                  << " " << this->m_InterpTimers[2][MIN]
                  << " " << this->m_InterpTimers[2][MAX]
                  << " " << this->m_InterpTimers[2][AVG]
                  << " " << this->m_InterpTimers[2][MAX] / static_cast<double>(count)
                  << std::endl;

        // if time has been logged
        count = this->m_Counter[IP] + 3*this->m_Counter[IPVEC];
        count = count > 0 ? count : 1;
        logwriter << "\"interp sort\""
                  << " " << count << std::scientific
                  << " " << this->m_InterpTimers[3][MIN]
                  << " " << this->m_InterpTimers[3][MAX]
                  << " " << this->m_InterpTimers[3][AVG]
                  << " " << this->m_InterpTimers[3][MAX] / static_cast<double>(count)
                  << std::endl;

        logwriter << "\"slowest proc tts\""
                  << " " << 0 << std::scientific
                  << " " << 0
                  << " " << this->m_TTSSlowest
                  << " " << 0
                  << " " << 0
                  << std::endl;

        logwriter << "\"slowest proc id\""
                  << " " << 0 << std::scientific
                  << " " << 0
                  << " " << this->m_IDSlowest
                  << " " << 0
                  << " " << 0
                  << std::endl;

        logwriter << "\"slowest proc fft\""
                  << " " << 0 << std::scientific
                  << " " << 0
                  << " " << this->m_FFTSlowest
                  << " " << 0
                  << " " << 0
                  << std::endl;

        logwriter << "\"slowest proc iterp\""
                  << " " << 0 << std::scientific
                  << " " << 0
                  << " " << this->m_IPSlowest
                  << " " << 0
                  << " " << 0
                  << std::endl;
    }

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write log file for workload
 *******************************************************************/
PetscErrorCode Logger::WriteWorkLoadLogReadable() {
    PetscErrorCode ierr = 0;
    std::string fn, path;
    std::ofstream logwriter;
    std::stringstream ss, ssnum;
    int rank;

    PetscFunctionBegin;

    this->Enter(__func__);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    path = this->m_FileNames.xfolder;

    // write out logfile
    if (rank == 0) {
        fn = path + "registration-performance.log";
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);
        std::cout << std::endl;
        ierr = this->WriteWorkLoadLogReadable(std::cout); CHKERRQ(ierr);
        std::cout << std::endl;
        ierr = this->WriteWorkLoadLogReadable(logwriter); CHKERRQ(ierr);
        logwriter.close();
    }

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}


/********************************************************************
 * @brief write log file for workload
 *******************************************************************/
PetscErrorCode Logger::WriteWorkLoadLogReadable(std::ostream& logwriter) {
    PetscErrorCode ierr = 0;
    std::string fn, line, path;
    std::stringstream ss, ssnum;
    int rank, nnum, nstr, nproc;

    PetscFunctionBegin;

    this->Enter(__func__);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

    // write out logfile
    if (rank == 0) {
        line = std::string(this->m_LineLength, '-');
        nnum = 20; nstr = 20;

        logwriter << line << std::endl;
        std::time_t result = std::time(NULL);
        logwriter << "# run finished on " << std::asctime(std::localtime(&result));
#ifdef GIT_VERSION
        logwriter << "# git version " << GIT_VERSION << std::endl;
#endif
        logwriter << std::scientific;
        logwriter << line << std::endl;
        logwriter << "# problem setup" << std::endl;
        logwriter << line << std::endl;
        ss  << this->m_Domain.nx[0] << " x "
            << this->m_Domain.nx[1] << " x "
            << this->m_Domain.nx[2];

        logwriter << std::left
                  << std::setw(nstr) << " nx" << std::right
                  << std::setw(nnum) << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        logwriter << std::left
                  << std::setw(nstr) << " nt" << std::right
                  << std::setw(nnum) << this->m_Domain.nt << std::endl;

        logwriter << std::left
                  << std::setw(nstr) << " n" << std::right
                  << std::setw(nnum) << this->m_Domain.ng << std::endl;

        logwriter << std::left
                  << std::setw(nstr) << " nl" << std::right
                  << std::setw(nnum) << this->m_Domain.nl << std::endl;

        logwriter << std::left
                  << std::setw(nstr) << " nmpi" << std::right
                  << std::setw(nnum) << nproc << std::endl;

        ss << this->m_CartGridDims[0] << " x " << this->m_CartGridDims[1];
        logwriter << std::left
                  << std::setw(nstr) << " proc grid" << std::right
                  << std::setw(nnum) << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        logwriter << std::left
                  << std::setw(nstr) << " num threads" << std::right
                  //<< this->m_NumThreads << std::endl;
                  << omp_get_max_threads() << std::endl;

        logwriter << std::endl;
        logwriter << line << std::endl;
        logwriter << "# timers" << std::endl;
        logwriter << line << std::endl;

        // write heading
        ss  << std::scientific << std::left
            << std::setw(nstr) << " " << std::right
            << std::setw(nnum) << "min(p)"
            << std::setw(nnum) << "max(p)"
            << std::setw(nnum) << "mean(p)"
            << std::setw(nnum) << "max(p)/numeval";
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());


        // if time has been logged
        if (this->m_Timer[T2SEXEC][LOG] > 0.0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " time to solution" << std::right
                << std::setw(nnum) << this->m_Timer[T2SEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[T2SEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[T2SEXEC][AVG]
                << std::setw(nnum) << "n/a";
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }


        // if time has been logged
        if (this->m_Timer[PDEEXEC][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[PDESOLVE] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " pde solves" << std::right
                << std::setw(nnum) << this->m_Timer[PDEEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[PDEEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[PDEEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[PDEEXEC][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[PDESOLVE]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[OBJEXEC][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[OBJEVAL] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " obj eval" << std::right
                << std::setw(nnum) << this->m_Timer[OBJEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[OBJEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[OBJEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[PDEEXEC][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[OBJEVAL]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[GRADEXEC][LOG] > 0.0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " grad eval" << std::right
                << std::setw(nnum) << this->m_Timer[GRADEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[GRADEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[GRADEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[GRADEXEC][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[GRADEVAL]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[HMVEXEC][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[HESSMATVEC] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " hess mat vec" << std::right
                << std::setw(nnum) << this->m_Timer[HMVEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[HMVEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[HMVEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[HMVEXEC][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[HESSMATVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }


        // if time has been logged
        if (this->m_Timer[PMVSETUP][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[PCMATVEC] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " pc mat vec" << std::right
                << std::setw(nnum) << this->m_Timer[PMVSETUP][MIN]
                << std::setw(nnum) << this->m_Timer[PMVSETUP][MAX]
                << std::setw(nnum) << this->m_Timer[PMVSETUP][AVG]
                << std::setw(nnum) << this->m_Timer[PMVSETUP][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[PCMATVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }


        // if time has been logged
        if (this->m_Timer[PMVEXEC][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[PCMATVEC] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " pc mat vec" << std::right
                << std::setw(nnum) << this->m_Timer[PMVEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[PMVEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[PMVEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[PMVEXEC][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[PCMATVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[FFTSELFEXEC][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[FFT] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " fft selfexec" << std::right
                << std::setw(nnum) << this->m_Timer[FFTSELFEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[FFTSELFEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[FFTSELFEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[FFTSELFEXEC][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[FFT]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_FFTAccumTime > 0.0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " FFT accumulated" << std::right
                << std::setw(nnum) << " n/a"
                << std::setw(nnum) << this->m_FFTAccumTime
                << std::setw(nnum) << "n/a"
                << std::setw(nnum) << "n/a";
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[FFTSETUP][LOG] > 0.0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " FFT setup" << std::right
                << std::setw(nnum) << this->m_Timer[FFTSETUP][MIN]
                << std::setw(nnum) << this->m_Timer[FFTSETUP][MAX]
                << std::setw(nnum) << this->m_Timer[FFTSETUP][AVG]
                << std::setw(nnum) << "n/a";
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_FFTTimers[FFTCOMM][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[FFT] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " FFT communication" << std::right
                << std::setw(nnum) << this->m_FFTTimers[FFTCOMM][MIN]
                << std::setw(nnum) << this->m_FFTTimers[FFTCOMM][MAX]
                << std::setw(nnum) << this->m_FFTTimers[FFTCOMM][AVG]
                << std::setw(nnum) << this->m_FFTTimers[FFTCOMM][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[FFT]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_FFTTimers[FFTEXECUTE][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[FFT] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " FFT execution" << std::right
                << std::setw(nnum) << this->m_FFTTimers[FFTEXECUTE][MIN]
                << std::setw(nnum) << this->m_FFTTimers[FFTEXECUTE][MAX]
                << std::setw(nnum) << this->m_FFTTimers[FFTEXECUTE][AVG]
                << std::setw(nnum) << this->m_FFTTimers[FFTEXECUTE][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[FFT]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[IPSELFEXEC][LOG] > 0.0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp selfexec" << std::right
                << std::setw(nnum) << this->m_Timer[IPSELFEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[IPSELFEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[IPSELFEXEC][AVG]
                << std::setw(nnum) << "n/a";
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_IPAccumTime > 0.0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp accumulated" << std::right
                << std::setw(nnum) << " n/a"
                << std::setw(nnum) << this->m_IPAccumTime
                << std::setw(nnum) << "n/a"
                << std::setw(nnum) << "n/a";
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }


        // if time has been logged
        if (this->m_InterpTimers[0][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[IP] > 0 && this->m_Counter[IPVEC] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp comm" << std::right
                << std::setw(nnum) << this->m_InterpTimers[0][MIN]
                << std::setw(nnum) << this->m_InterpTimers[0][MAX]
                << std::setw(nnum) << this->m_InterpTimers[0][AVG]
                << std::setw(nnum) << this->m_InterpTimers[0][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[IP]
                                    + this->m_Counter[IPVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_InterpTimers[1][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[IP] > 0 && this->m_Counter[IPVEC] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp exec" << std::right
                << std::setw(nnum) << this->m_InterpTimers[1][MIN]
                << std::setw(nnum) << this->m_InterpTimers[1][MAX]
                << std::setw(nnum) << this->m_InterpTimers[1][AVG]
                << std::setw(nnum) << this->m_InterpTimers[1][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[IP]
                                    + this->m_Counter[IPVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_InterpTimers[2][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[IP] > 0 && this->m_Counter[IPVEC] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp alloc" << std::right
                << std::setw(nnum) << this->m_InterpTimers[2][MIN]
                << std::setw(nnum) << this->m_InterpTimers[2][MAX]
                << std::setw(nnum) << this->m_InterpTimers[2][AVG]
                << std::setw(nnum) << this->m_InterpTimers[2][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[IP]
                                    + this->m_Counter[IPVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_InterpTimers[3][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[IP] > 0 && this->m_Counter[IPVEC] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp sort" << std::right
                << std::setw(nnum) << this->m_InterpTimers[3][MIN]
                << std::setw(nnum) << this->m_InterpTimers[3][MAX]
                << std::setw(nnum) << this->m_InterpTimers[3][AVG]
                << std::setw(nnum) << this->m_InterpTimers[3][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[IP]
                                    + this->m_Counter[IPVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        logwriter << std::endl;
        logwriter << line << std::endl;
        logwriter << "# counters" << std::endl;
        logwriter << line << std::endl;

        if (this->m_Counter[OBJEVAL] > 0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " objective evals" << std::right
                << std::setw(nnum) << this->m_Counter[OBJEVAL];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[GRADEVAL] > 0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " grad evals" << std::right
                << std::setw(nnum) << this->m_Counter[GRADEVAL];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[PDESOLVE] > 0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " PDE solves" << std::right
                << std::setw(nnum) << this->m_Counter[PDESOLVE];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[HESSMATVEC] > 0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " hess mat vecs" << std::right
                << std::setw(nnum) << this->m_Counter[HESSMATVEC];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[PCMATVEC] > 0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " pc mat vecs" << std::right
                << std::setw(nnum) << this->m_Counter[PCMATVEC];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[IP] > 0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " ips sca" << std::right
                << std::setw(nnum) << this->m_Counter[IP];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[IPVEC] > 0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " ips vec" << std::right
                << std::setw(nnum) << this->m_Counter[IPVEC];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[FFT] > 0) {
            ss  << std::scientific << std::left
                << std::setw(nstr) << " ffts" << std::right
                << std::setw(nnum) << this->m_Counter[FFT];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }
    }

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write residual to file
 *******************************************************************/
PetscErrorCode Logger::WriteFinalResidualLog() {
    PetscErrorCode ierr = 0;
    int rank, nnum;
    std::ofstream logwriter;
    std::stringstream ss;
    std::string path, fn;
    PetscFunctionBegin;

    this->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    nnum = 20;
    path = this->m_FileNames.xfolder;

    if (rank == 0) {
        // create output file
        fn = path + "claire-residual.log";
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-mT||_2" << std::right
            << std::setw(nnum) << this->m_Log.finalresidual[0];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-mT||_infty" << std::right
            << std::setw(nnum) << this->m_Log.finalresidual[1];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-m1||_2" << std::right
            << std::setw(nnum) << this->m_Log.finalresidual[2];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-m1||_infty" << std::right
            << std::setw(nnum) << this->m_Log.finalresidual[3];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        this->m_Log.finalresidual[0] = (this->m_Log.finalresidual[0] > 0.0) ? this->m_Log.finalresidual[0] : 1.0;
        this->m_Log.finalresidual[1] = (this->m_Log.finalresidual[1] > 0.0) ? this->m_Log.finalresidual[1] : 1.0;

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-m1||_2,rel" << std::right
            << std::setw(nnum) <<  this->m_Log.finalresidual[2]/this->m_Log.finalresidual[0];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-m1||_infty,rel" << std::right
            << std::setw(nnum) <<  this->m_Log.finalresidual[3]/this->m_Log.finalresidual[1];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        // close logger
        logwriter.close();
    }

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write out logging information for krylov method
 *******************************************************************/
PetscErrorCode Logger::WriteConvergenceLog() {
    PetscErrorCode ierr = 0;
    int rank, n;
    std::ofstream logwriter;
    std::stringstream ss;
    std::string path, fn;
    PetscFunctionBegin;

    this->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    path = this->m_FileNames.xfolder;

    if (rank == 0) {
        // create output file
        fn = path + "claire-distance-measure-trend.log";
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

        n = static_cast<int>(this->m_Log.distance.size());
        for (int i = 0; i < n; ++i) {
            ss << std::scientific << std::right
               << std::setw(2) << this->m_Log.outeriterations[i]
               << std::setw(20) << this->m_Log.distance[i];
            logwriter << ss.str() << std::endl;
            ss.str(std::string()); ss.clear();
        }
        logwriter.close();  // close logger

        // create output file
        fn = path + "claire-regularization-trend.log";
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

        n = static_cast<int>(this->m_Log.regularization.size());
        for (int i = 0; i < n; ++i) {
            ss << std::scientific << std::right
               << std::setw(2) << this->m_Log.outeriterations[i]
               << std::setw(20) << this->m_Log.regularization[i];
            logwriter << ss.str() << std::endl;
            ss.str(std::string()); ss.clear();
        }
        logwriter.close();  // close logger


        // create output file
        fn = path + "claire-objective-trend.log";
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

        n = static_cast<int>(this->m_Log.objective.size());
        for (int i = 0; i < n; ++i) {
            ss << std::scientific << std::right
               << std::setw(2) << this->m_Log.outeriterations[i]
               << std::setw(20) << this->m_Log.objective[i];
            logwriter << ss.str() << std::endl;
            ss.str(std::string()); ss.clear();
        }
        logwriter.close();  // close logger



    }

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write out logging information for krylov method
 *******************************************************************/
PetscErrorCode Logger::WriteKSPLog() {
    PetscErrorCode ierr = 0;
    int rank, n;
    std::ofstream logwriter;
    std::stringstream ss;
    std::string path, fn;
    PetscFunctionBegin;

    this->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    path = this->m_FileNames.xfolder;

    if (rank == 0) {
        // create output file
        fn = path + "claire-krylov-method-residual.log";
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

        n = static_cast<int>(this->m_Log.kspresidual.size());
        for (int i = 0; i < n; ++i) {
            ss << std::scientific << std::right
               << std::setw(2) << this->m_Log.kspiterations[i]
               << std::setw(20) << this->m_Log.kspresidual[i];
            logwriter << ss.str() << std::endl;
            ss.str(std::string()); ss.clear();
        }
        logwriter.close();  // close logger
    }

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief displays the global exection time
 *******************************************************************/
PetscErrorCode Logger::DisplayTimeToSolution() {
    PetscErrorCode ierr = 0;
    double hours, minutes, seconds, millisec, time;
    int rank;
    std::stringstream ss;
    std::string line;
    PetscFunctionBegin;

    this->Enter(__func__);

    time = this->m_Timer[T2SEXEC][MAX];

    hours    = time / 3600.0;
    minutes  = (hours   - floor(hours))   *   60.0;
    seconds  = (minutes - floor(minutes)) *   60.0;
    millisec = (seconds - floor(seconds)) * 1000.0;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    line = std::string(this->m_LineLength, '-');

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%s\n", line.c_str()); CHKERRQ(ierr);

    ss  << "computation finished (elapsed cpu time "
        << floor(hours) << " h "
        << floor(minutes) << " m "
        << floor(seconds) << " s "
        << floor(millisec) << " ms; "
        << std::scientific << time << ")";

    ierr = Msg(ss.str()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%s\n", line.c_str()); CHKERRQ(ierr);

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}



#endif  // _LOGGER_CPP_
