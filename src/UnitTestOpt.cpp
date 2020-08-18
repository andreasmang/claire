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

#ifndef _UNITTESTOPT_CPP_
#define _UNITTESTOPT_CPP_

#include "UnitTestOpt.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
UnitTestOpt::UnitTestOpt() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
UnitTestOpt::UnitTestOpt(int argc, char** argv) {
    this->Initialize();
    this->ParseArguments(argc,argv);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
UnitTestOpt::UnitTestOpt(const UnitTestOpt& opt) {
    this->Initialize();
    this->Copy(opt);
}




/********************************************************************
 * @brief parse user arguments
 *******************************************************************/
PetscErrorCode UnitTestOpt::ParseArguments(int argc, char** argv) {
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
        } else if (strcmp(argv[1], "-rep") == 0) {
            argc--; argv++;
            this->rep = atoi(argv[1]);
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
        } else if (strcmp(argv[1], "-verbosity") == 0) {
            argc--; argv++;
            this->m_Verbosity = std::min(atoi(argv[1]),2);
        } else if (strcmp(argv[1], "-debug") == 0) {
            this->m_Verbosity = 3;
        } else if (strcmp(argv[1], "-regnorm") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "h1s") == 0) {
                this->m_RegNorm.type = H1SN;
                this->m_RegModel = COMPRESSIBLE;
            } else if (strcmp(argv[1], "h2s") == 0) {
                this->m_RegNorm.type = H2SN;
                this->m_RegModel = COMPRESSIBLE;
            } else if (strcmp(argv[1], "h3s") == 0) {
                this->m_RegNorm.type = H3SN;
                this->m_RegModel = COMPRESSIBLE;
            } else if (strcmp(argv[1], "h1") == 0) {
                this->m_RegNorm.type = H1;
                this->m_RegModel = COMPRESSIBLE;
            } else if (strcmp(argv[1], "h1s-div") == 0) {
                this->m_RegNorm.type = H1SN;
                this->m_RegModel = RELAXEDSTOKES;
            } else if (strcmp(argv[1], "h1s-stokes") == 0) {
                this->m_RegNorm.type = H1SN;
                this->m_RegModel = STOKES;
            } else if (strcmp(argv[1], "h2") == 0) {
                this->m_RegNorm.type = H2;
                this->m_RegModel = COMPRESSIBLE;
            } else if (strcmp(argv[1], "h3") == 0) {
                this->m_RegNorm.type = H3;
                this->m_RegModel = COMPRESSIBLE;
            } else if (strcmp(argv[1], "l2") == 0) {
                this->m_RegNorm.type = L2;
                this->m_RegModel = COMPRESSIBLE;
            } else {
                msg = "\n\x1b[31m regularization norm not available: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-beta") == 0) {
            argc--; argv++;
            this->m_RegNorm.beta[0] = atof(argv[1]);
//            this->m_RegNorm.beta[1] = atof(argv[1]);
        } else if (strcmp(argv[1], "-beta-div") == 0) {
            argc--; argv++;
            this->m_RegNorm.beta[2] = atof(argv[1]);
        } else {
            bool found;
            ierr = this->CheckTests(argv[1], found); CHKERRQ(ierr);
            if (!found) {
              msg = "\n\x1b[31m argument not valid: %s\x1b[0m\n";
              ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
              ierr = this->Usage(); CHKERRQ(ierr);
            }
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
UnitTestOpt::~UnitTestOpt() {
    this->ClearMemory();
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode UnitTestOpt::ClearMemory() {
    PetscFunctionBegin;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief initialize class variables
 *******************************************************************/
PetscErrorCode UnitTestOpt::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->SuperClass::Initialize(); CHKERRQ(ierr);

    this->m_TestType = TestType::None;
    this->rep = 1;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief display usage message for binary
 *******************************************************************/
PetscErrorCode UnitTestOpt::Usage(bool advanced) {
    PetscErrorCode ierr = 0;
    int rank;
    std::string line;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    line = std::string(this->m_LineLength,'-');

    if (rank == 0) {
        //this->m_Parser.help("test");
        //std::cout << std::endl;
        //std::cout << line << std::endl;
        std::cout << std::endl;
        std::cout << line << std::endl;
        std::cout << " usage: regbenchmark [options] " <<std::endl;
        std::cout << line << std::endl;
        std::cout << " where [options] is one or more of the following"<<std::endl;
        std::cout << line << std::endl;
        // ####################### advanced options #######################
        ierr = this->PrintTests(); CHKERRQ(ierr);
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
        std::cout << " -regnorm <type>             regularization norm for velocity field" << std::endl;
        std::cout << "                             <type> is one of the following" << std::endl;
        std::cout << "                                 h1s          H1-seminorm" << std::endl;
        std::cout << "                                 h1s-div      H1-seminorm with penalty on div(v); penalty is" << std::endl;
        std::cout << "                                              controlled by regulariztion parameter '-beta-div' (see" << std::endl;
        std::cout << "                                              below); default value for regularization parameter" << std::endl;
        std::cout << "                                              is 1E-4;" << std::endl;
        std::cout << "                                 h1s-stokes   H1-seminorm with incompressiblity constraint (i.e.," << std::endl;
        std::cout << "                                              the jacobian is 1; divergence free velocity)" << std::endl;
        std::cout << "                                 h2s          H2-seminorm" << std::endl;
//        std::cout << "                                 h3s          H3-seminorm" << std::endl;
//        std::cout << "                                 h1           H1-norm" << std::endl;
//        std::cout << "                                 h2           H2-norm" << std::endl;
//        std::cout << "                                 h3           H3-norm" << std::endl;
        std::cout << "                                 l2           l2-norm (discouraged)" << std::endl;
        std::cout << " -beta  <dbl>                set constant regularization parameter for velocity field (default: 1E-2)" << std::endl;
        std::cout << " -beta-div <dbl>             set constant regularization parameter for mass source map (default: 1E-4);" << std::endl;
        std::cout << "                             this parameter controls a penalty on divergence of the velocity, i.e.," << std::endl;
        std::cout << "                             the incompressibility; the penalty is enabled via '-regnorm h1-div'" << std::endl;
        std::cout << "                             option; details can be found above;" << std::endl;
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
 * @brief check the arguments set by user
 *******************************************************************/
PetscErrorCode UnitTestOpt::CheckArguments() {
    PetscErrorCode ierr = 0;
    std::string msg;
    PetscFunctionBegin;

//    ierr = Assert(this->m_NumThreads > 0, "omp threads < 0"); CHKERRQ(ierr);

    if (this->m_TestType == TestType::None) {
        msg = "\x1b[31m you need to define a unit test\x1b[0m\n";
        ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
        ierr = this->Usage(); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(ierr);
}

PetscErrorCode UnitTestOpt::PrintTests() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  std::cout << " -interp                     interpolation unit test"<<std::endl;
  std::cout << " -forward                    forward solver unit test"<<std::endl;
  std::cout << " -trajectory                 trajectory solver unit test"<<std::endl;
  std::cout << " -gradient                   gradient unit test"<<std::endl;
  std::cout << " -hessian                    hessian unit test"<<std::endl;
  std::cout << " -reg                        regularization unit test"<<std::endl;
  std::cout << " -diff                       differentiation unit test"<<std::endl;
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode UnitTestOpt::CheckTests(char* argv, bool &found) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  if (strcmp(argv, "-interp") == 0) {
    this->m_TestType = TestType::Interpolate;
  } else if (strcmp(argv, "-forward") == 0) {
    this->m_TestType = TestType::ForwardSolver;
  } else if (strcmp(argv, "-trajectory") == 0) {
    this->m_TestType = TestType::Trajectory;
  } else if (strcmp(argv, "-gradient") == 0) {
    this->m_TestType = TestType::Gradient;
  } else if (strcmp(argv, "-hessian") == 0) {
    this->m_TestType = TestType::Hessian;
  } else if (strcmp(argv, "-reg") == 0) {
    this->m_TestType = TestType::Reg;
  } else if (strcmp(argv, "-diff") == 0) {
    this->m_TestType = TestType::Diff;
  } else {
    found = false;
  }
  found = true;
  PetscFunctionReturn(ierr);
}


PetscErrorCode UnitTestOpt::Run() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  switch(this->m_TestType) {
  case TestType::All:
    ierr = UnitTest::TestInterpolation(this); CHKERRQ(ierr);
    ierr = UnitTest::TestInterpolationMultiGPU(this); CHKERRQ(ierr);
    ierr = UnitTest::TestForwardSolver(this); CHKERRQ(ierr);
    ierr = UnitTest::TestTrajectory(this); CHKERRQ(ierr);
    ierr = UnitTest::TestGradient(this); CHKERRQ(ierr);
    ierr = UnitTest::TestHessian(this); CHKERRQ(ierr);
    ierr = UnitTest::TestRegularization(this); CHKERRQ(ierr);
    break;
  case TestType::Interpolate:
    if (nprocs == 1) {
      ierr = UnitTest::TestInterpolation(this); CHKERRQ(ierr);
    } else {
      ierr = UnitTest::TestInterpolationMultiGPU(this); CHKERRQ(ierr);
    }
    break;
  case TestType::ForwardSolver:
    ierr = UnitTest::TestForwardSolver(this); CHKERRQ(ierr);
    break;
  case TestType::Trajectory:
//    ierr = UnitTest::TestTrajectory(this); CHKERRQ(ierr);
    ierr = UnitTest::TestTrajectoryMultiGPU(this); CHKERRQ(ierr);
    break;
  case TestType::Gradient:
    ierr = UnitTest::TestGradient(this); CHKERRQ(ierr);
    break;
  case TestType::Hessian:
    ierr = UnitTest::TestHessian(this); CHKERRQ(ierr);
    break;
  case TestType::Reg:
    ierr = UnitTest::TestRegularization(this); CHKERRQ(ierr);
    break;
  case TestType::Diff:
    ierr = UnitTest::TestDifferentiation(this); CHKERRQ(ierr);
    //ierr = UnitTest::TestDifferentiationMultiGPU(this); CHKERRQ(ierr);
    break;
  default:
    break;
  };
  
#ifdef ZEITGEIST
    if (this->m_Domain.level == 0) {
      Msg("-----------------------------------------------------------------------------------------------------");
      Msg("ZeitGeist:");
      for (auto zg : ZeitGeist::zgMap()) {
        char txt[120];
        double global_runtime;
        double local_runtime = zg.second.Total_s();
        MPI_Reduce(&local_runtime, &global_runtime, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
        sprintf(txt, "  %16s: %5lix, %0.10lf",zg.first.c_str(), zg.second.Count(), global_runtime);
        Msg(txt);
      }
      Msg("-----------------------------------------------------------------------------------------------------");
    }
#endif
  
  PetscFunctionReturn(ierr);
}


}  // namespace reg




#endif  // _REGBENCHMARKOPT_CPP_
