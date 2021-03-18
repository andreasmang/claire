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

#ifndef _CLPBENCHMARK_CPP_
#define _CLPBENCHMARK_CPP_

#include "BenchmarkOpt.hpp"

#define _TO_STR(s) #s
#define TO_STR(s) _TO_STR(s)


namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
BenchmarkOpt::BenchmarkOpt() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
BenchmarkOpt::BenchmarkOpt(int argc, char** argv) {
    this->Initialize();
    this->ParseArguments(argc,argv);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
BenchmarkOpt::BenchmarkOpt(const BenchmarkOpt& opt) {
    this->Initialize();
    this->Copy(opt);
}




/********************************************************************
 * @brief parse user arguments
 *******************************************************************/
PetscErrorCode BenchmarkOpt::ParseArguments(int argc, char** argv) {
    PetscErrorCode ierr = 0;
    std::string msg;
    std::vector<int> nx;
    std::vector<int> nxr;
    std::vector<int> np;
    PetscFunctionBegin;

    if (argc == 1) {
        ierr = this->Usage(); CHKERRQ(ierr);
    }

    while (argc > 1) {
        if ( (strcmp(argv[1], "-help") == 0)
            || (strcmp(argv[1], "-h") == 0)
            || (strcmp(argv[1], "-HELP") == 0) ) {
            ierr = this->Usage(); CHKERRQ(ierr);
        } else if (strcmp(argv[1], "-advanced") == 0) {
            ierr = this->Usage(true); CHKERRQ(ierr);
        } else if (strcmp(argv[1], "-x") == 0) {
            argc--; argv++;
            this->m_FileNames.xfolder = argv[1];
        } else if (strcmp(argv[1], "-nt") == 0) {
            argc--; argv++;
            this->m_Domain.nt = static_cast<IntType>(atoi(argv[1]));
        } else if (strcmp(argv[1], "-nx") == 0) {
            argc--; argv++;
            const std::string nxinput = argv[1];

            // strip the "x" in the string to get the numbers
            nx = String2Vec( nxinput );

            if (nx.size() == 1) {
                for(int i=0; i < 3; ++i) {
                    this->m_Domain.nx[i] = static_cast<IntType>(nx[0]);
                }
            } else if (nx.size() == 3) {
                for(int i=0; i < 3; ++i) {
                    this->m_Domain.nx[i] = static_cast<IntType>(nx[i]);
                }
            } else {
                msg = "\n\x1b[31m error in grid size argument: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-pdesolver") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "rk2") == 0) {
                this->m_PDESolver.type = RK2;
            } else if (strcmp(argv[1], "rk2a") == 0) {
                this->m_PDESolver.type = RK2A;
            } else if (strcmp(argv[1], "sl") == 0) {
                this->m_PDESolver.type = SL;
            } else {
                msg = "\n\x1b[31m pde solver not implemented: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-adapttimestep") == 0) {
            this->m_PDESolver.adapttimestep = true;
        } else if (strcmp(argv[1], "-cflnumber") == 0) {
            argc--; argv++;
            this->m_PDESolver.cflnumber = atof(argv[1]);
        } else if (strcmp(argv[1], "-iporder") == 0) {
            argc--; argv++;
            this->m_PDESolver.iporder = atoi(argv[1]);
        } else if (strcmp(argv[1], "-nthreads") == 0) {
            argc--; argv++;
            this->m_NumThreads = atoi(argv[1]);
        } else if (strcmp(argv[1], "-np") == 0) {
            argc--; argv++;
            const std::string npinput = argv[1];

            // strip the "x" in the string to get the numbers
            np = String2Vec(npinput);
            if (np.size() == 1) {
                for(int i=0; i < 2; ++i) {
                    this->m_CartGridDims[i] = static_cast<unsigned int>(np[0]);
                }
            } else if (np.size() == 2) {
                for(int i=0; i < 2; ++i) {
                    this->m_CartGridDims[i] = static_cast<unsigned int>(np[i]);
                }
            } else {
                msg="\n\x1b[31m error in number of procs: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-forward") == 0) {
            this->m_BenchmarkID = 0;
        } else if (strcmp(argv[1], "-gradient") == 0) {
            this->m_BenchmarkID = 1;
        } else if (strcmp(argv[1], "-hessmatvec") == 0) {
            this->m_BenchmarkID = 2;
        } else if (strcmp(argv[1], "-terror") == 0) {
            this->m_BenchmarkID = 3;
        } else if (strcmp(argv[1], "-repeats") == 0) {
            argc--; argv++;
            this->m_NumRepeats = atoi(argv[1]);
        } else if (strcmp(argv[1], "-verbosity") == 0) {
            argc--; argv++;
            this->m_Verbosity = std::min(atoi(argv[1]),2);
        } else if (strcmp(argv[1], "-logwork") == 0) {
            this->m_Log.enabled[LOGLOAD] = true;
        } else if (strcmp(argv[1], "-debug") == 0) {
            this->m_Verbosity = 3;
        } else {
            msg = "\n\x1b[31m argument not valid: %s\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
            ierr = this->Usage(); CHKERRQ(ierr);
        }
        argc--; argv++;
    }

    // check the arguments/parameters set by the user
    ierr = this->CheckArguments(); CHKERRQ(ierr);

    // set number of threads
    ierr = InitializeDataDistribution(this->m_NumThreads, this->m_CartGridDims,
                                      this->m_Domain.mpicomm, (this->m_Domain.mpicomm != MPI_COMM_NULL)); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief default constructor
 *******************************************************************/
BenchmarkOpt::~BenchmarkOpt() {
    this->ClearMemory();
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode BenchmarkOpt::ClearMemory() {
    PetscFunctionBegin;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief initialize class variables
 *******************************************************************/
PetscErrorCode BenchmarkOpt::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->SuperClass::Initialize(); CHKERRQ(ierr);

    this->m_BenchmarkID = -1;
    this->m_NumRepeats = 1;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief display usage message for binary
 *******************************************************************/
PetscErrorCode BenchmarkOpt::Usage(bool advanced) {
    PetscErrorCode ierr = 0;
    int rank;
    std::string line;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    line = std::string(this->m_LineLength,'-');

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << line << std::endl;
        std::cout << " usage: regbenchmark [options] " <<std::endl;
        std::cout << line << std::endl;
        std::cout << " where [options] is one or more of the following"<<std::endl;
        // ####################### advanced options #######################
        std::cout << " -forward                    benchmark forward solver"<<std::endl;
        std::cout << " -gradient                   benchmark gradient evaluation"<<std::endl;
        std::cout << " -hessmatvec                 benchmark hessian matvec evaluation"<<std::endl;
        std::cout << " -repeats <int>              set number of repeats"<<std::endl;
        std::cout << " -terror                     compute numerical error for solution of transport equation"<<std::endl;
        std::cout << " -logwork                    log work load (requires -x option)"<<std::endl;
        if (advanced) {
        std::cout << line << std::endl;
        std::cout << " memory distribution and parallelism"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -nthreads <int>             number of threads (default: 1)"<<std::endl;
        std::cout << " -np <int>x<int>             distribution of mpi tasks (cartesian grid) (example: -np 2x4 results"<<std::endl;
        std::cout << "                             results in MPI distribution of size (nx1/2,nx2/4,nx3) for each mpi task)"<<std::endl;
        }
        // ####################### advanced options #######################
        std::cout << line << std::endl;
        std::cout << " ### input parameters"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -x <path>                   output path (by default only deformed template image and velocity"<<std::endl;
        std::cout << "                             field will be written; for more output options, see flags;"<<std::endl;
        std::cout << "                             a prefix can be added by, e.g., doing '-x </path/prefix_>"<<std::endl;
        std::cout << " -nx <int>x<int>x<int>       grid size (e.g., 32x64x32); allows user to control grid size for synthetic" << std::endl;
        std::cout << "                             problems; assumed to be uniform if single integer is provided" << std::endl;
        std::cout << line << std::endl;
        std::cout << " solver specific parameters (numerics)"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -pdesolver <type>           numerical time integrator for transport equations"<<std::endl;
        std::cout << "                             <type> is one of the following"<<std::endl;
        std::cout << "                                 sl           semi-Lagrangian method (default; unconditionally stable)"<<std::endl;
        std::cout << "                                 rk2          rk2 time integrator (conditionally stable)"<<std::endl;
        std::cout << " -nt <int>                   number of time points (for time integration; default: 4)"<<std::endl;
        std::cout << " -adapttimestep              vary number of time steps according to defined number"<<std::endl;
        std::cout << " -cflnumber <dbl>            set cfl number"<<std::endl;
        std::cout << " -iporder <int>              order of interpolation model (default is 3)" << std::endl;
        std::cout << line << std::endl;
        // ####################### advanced options #######################
        std::cout << line << std::endl;
        std::cout << " -usenc                      use netcdf format os output (*.nc; default is *.nii.gz)"<<std::endl;
        std::cout << " -verbosity <int>            verbosity level (ranges from 0 to 3; default: 1)"<<std::endl;
        std::cout << " -help                       display a brief version of the user message"<<std::endl;
        std::cout << " -advanced                   display this message"<<std::endl;
        std::cout << line << std::endl;
        std::cout << line << std::endl;
    }

    ierr = PetscFinalize(); CHKERRQ(ierr);
    exit(0);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief display options
 *******************************************************************/
PetscErrorCode BenchmarkOpt::DisplayOptions() {
    PetscErrorCode ierr = 0;
    int rank, indent;
    std::string msg, line;

    PetscFunctionBegin;

    this->Enter(__func__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    indent = 40;
    line = std::string(this->m_LineLength, '-');

    // display the parameters (only on rank 0)
    if (rank == 0) {
        std::cout << std::endl;
        std::cout << line << std::endl;
        std::cout << " CLAIRE: Constrained Large Deformation Diffeomorphic Registration" << std::endl;
        std::cout << line << std::endl;
        std::time_t result = std::time(NULL);
#ifdef GIT_VERSION
        std::cout << " Version " << TO_STR(GIT_VERSION) << " " << std::asctime(std::localtime(&result));
#else
        std::cout << " " << std::asctime(std::localtime(&result));
#endif
        std::cout << line << std::endl;
        std::cout << " problem setup" << std::endl;
        std::cout << line << std::endl;
        std::cout << std::left << std::setw(indent) << " problem dimensions"
                  << "(nx1,nx2,nx3,nt)=("
                  << this->m_Domain.nx[0] << ", "
                  << this->m_Domain.nx[1] << ", "
                  << this->m_Domain.nx[2] << ", "
                  << this->m_Domain.nt << ")" << std::endl;
        std::cout << std::left << std::setw(indent) << " network dimensions"
                  << this->m_CartGridDims[0] << "x"
                  << this->m_CartGridDims[1] << std::endl;
        std::cout << std::left << std::setw(indent) << " threads"
                  //<< this->m_NumThreads << std::endl;
                  << omp_get_max_threads() << std::endl;
        std::cout << std::left << std::setw(indent) << " (ng,nl)"
                  << "(" << this->m_Domain.ng << ", "
                  << this->m_Domain.nl << ")" << std::endl;
        std::cout << line << std::endl;
    }  // rank

    this->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief check the arguments set by user
 *******************************************************************/
PetscErrorCode BenchmarkOpt::CheckArguments() {
    PetscErrorCode ierr = 0;
    std::string msg;
    bool log = false;
    PetscFunctionBegin;

//    ierr = Assert(this->m_NumThreads > 0, "omp threads < 0"); CHKERRQ(ierr);

    if (this->m_BenchmarkID == -1) {
        msg = "\x1b[31m you need to define a benchmark test\x1b[0m\n";
        ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
        ierr = this->Usage(); CHKERRQ(ierr);
    }

    for (int i = 0; i < NLOGFLAGS; ++i) {
        if (this->m_Log.enabled[i]) {
            log = true;
        }
    }

    if (log) {
        if (this->m_FileNames.xfolder.empty()) {
            msg = "\x1b[31m output folder needs to be set (-x option) \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
    }
    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _REGBENCHMARKOPT_CPP_
