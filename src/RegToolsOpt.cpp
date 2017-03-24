#ifndef _REGTOOLSOPT_CPP_
#define _REGTOOLSOPT_CPP_

#include "RegToolsOpt.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegToolsOpt::RegToolsOpt() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
RegToolsOpt::RegToolsOpt(int argc, char** argv) {
    this->Initialize();
    this->ParseArguments(argc,argv);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
RegToolsOpt::RegToolsOpt(const RegToolsOpt& opt) {
    this->Initialize();
    this->Copy(opt);
}




/********************************************************************
 * @brief parse user arguments
 *******************************************************************/
PetscErrorCode RegToolsOpt::ParseArguments(int argc, char** argv) {
    PetscErrorCode ierr = 0;
    std::string msg;
    std::vector<int> nx;
    std::vector<int> nxr;
    std::vector<int> np;
    std::vector<int> sigma;
    PetscFunctionBegin;

    if (argc == 1) {
        ierr = this->Usage(); CHKERRQ(ierr);
    }

    while(argc > 1) {
        if ( (strcmp(argv[1], "-help") == 0)
            || (strcmp(argv[1], "-h") == 0)
            || (strcmp(argv[1], "-HELP") == 0) ) {
            ierr = this->Usage(); CHKERRQ(ierr);
        } else if (strcmp(argv[1], "-advanced") == 0) {
            ierr = this->Usage(true); CHKERRQ(ierr);
        } else if (strcmp(argv[1], "-mr") == 0) {
            argc--; argv++;
            this->m_RFN = argv[1];
        } else if (strcmp(argv[1], "-mt") == 0) {
            argc--; argv++;
            this->m_TFN = argv[1];
        } else if (strcmp(argv[1], "-ifile") == 0) {
            argc--; argv++;
            this->m_iScaFieldFN = argv[1];
        } else if (strcmp(argv[1], "-ivecx1") == 0) {
            argc--; argv++;
            this->m_iVecFieldX1FN = argv[1];
        } else if (strcmp(argv[1], "-ivecx2") == 0) {
            argc--; argv++;
            this->m_iVecFieldX2FN = argv[1];
        } else if (strcmp(argv[1], "-ivecx3") == 0) {
            argc--; argv++;
            this->m_iVecFieldX3FN = argv[1];
        } else if (strcmp(argv[1], "-x") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.xfolder = argv[1];
        } else if (strcmp(argv[1], "-i") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.ifolder = argv[1];
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
        } else if (strcmp(argv[1], "-monitorcflnumber") == 0) {
            this->m_PDESolver.monitorcflnumber = true;
        } else if (strcmp(argv[1], "-interpolationorder") == 0) {
            argc--; argv++;
            this->m_PDESolver.interpolationorder = atoi(argv[1]);
        } else if (strcmp(argv[1], "-sigma") == 0) {
            argc--; argv++;
            const std::string sigmainput = argv[1];
            // strip the "x" in the string to get the numbers
            sigma = String2Vec( sigmainput );

            if (sigma.size() == 1) {
                for (int i=0; i < 3; ++i) {
                    this->m_Sigma[i] = static_cast<ScalarType>(sigma[0]);
                }
            } else if (sigma.size() == 3) {
                for (int i=0; i < 3; ++i) {
                    this->m_Sigma[i] = static_cast<IntType>(sigma[i]);
                }
            } else {
                msg = "\n\x1b[31m error in smoothing kernel size: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-disablesmoothing") == 0) {
            this->m_RegFlags.applysmoothing = false;
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
        } else if (strcmp(argv[1], "-convert") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "2nii") == 0) {
                this->m_ReadWriteFlags.extension = ".nii.gz";
            } else if (strcmp(argv[1], "2nc") == 0) {
                this->m_ReadWriteFlags.extension = ".nc";
            }
            this->m_RegToolFlags.convert = true;
        } else if (strcmp(argv[1], "-usenc") == 0) {
            this->m_ReadWriteFlags.extension = ".nc";
        } else if (strcmp(argv[1], "-xresults") == 0) {
            this->m_ReadWriteFlags.results = true;
        } else if (strcmp(argv[1], "-defgrad") == 0) {
            this->m_ReadWriteFlags.defgrad = true;
        } else if (strcmp(argv[1], "-detdefgrad") == 0) {
            this->m_ReadWriteFlags.detdefgrad = true;
        } else if (strcmp(argv[1], "-invdetdefgrad") == 0) {
            this->m_RegFlags.invdefgrad = true;
            this->m_ReadWriteFlags.detdefgrad = true;
        } else if (strcmp(argv[1], "-defmap") == 0) {
            this->m_ReadWriteFlags.defmap = true;
        } else if (strcmp(argv[1], "-deffield") == 0) {
            this->m_ReadWriteFlags.deffield = true;
        } else if (strcmp(argv[1], "-timeseries") == 0) {
            this->m_ReadWriteFlags.timeseries = true;
        } else if (strcmp(argv[1], "-detdefgradfromdeffield") == 0) {
            this->m_RegFlags.detdefgradfromdeffield = true;
        } else if (strcmp(argv[1], "-residual") == 0) {
            this->m_RegToolFlags.computeresidual = true;
        } else if (strcmp(argv[1], "-error") == 0) {
            this->m_RegToolFlags.computeerror = true;
        } else if (strcmp(argv[1], "-grad") == 0) {
            this->m_RegToolFlags.computegrad = true;
        } else if (strcmp(argv[1], "-tscafield") == 0) {
            this->m_RegToolFlags.tscafield = true;
        } else if (strcmp(argv[1], "-smooth") == 0) {
            this->m_RegToolFlags.applysmoothing = true;
            argc--; argv++;
            const std::string sigmainput = argv[1];

            // strip the "x" in the string to get the numbers
            sigma = String2Vec(sigmainput);

            if (sigma.size() == 1) {
                for (int i = 0; i < 3; ++i) {
                    this->m_Sigma[i] = static_cast<ScalarType>(sigma[0]);
                }
            } else if (sigma.size() == 3) {
                for (int i = 0; i < 3; ++i) {
                    this->m_Sigma[i] = static_cast<IntType>(sigma[i]);
                }
            } else {
                msg = "\n\x1b[31m error in smoothing kernel size: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-tlabelmap") == 0) {
            this->m_RegToolFlags.tlabelmap = true;
        } else if (strcmp(argv[1], "-csynvel") == 0) {
            this->m_RegToolFlags.computesynvel = true;
        } else if (strcmp(argv[1], "-checkfwdsolveerr") == 0) {
            this->m_RegToolFlags.checkfwdsolveerr = true;
        } else if (strcmp(argv[1], "-checkfwdsolvetts") == 0) {
            argc--; argv++;
            this->m_RegToolFlags.numrepeat = atoi(argv[1]);
            this->m_RegToolFlags.checkfwdsolvetts = true;
        } else if (strcmp(argv[1], "-checkadjsolve") == 0) {
            this->m_RegToolFlags.checkadjsolve = true;
        } else if (strcmp(argv[1], "-checkdetdefgradsolve") == 0) {
            this->m_RegToolFlags.checkdetdefgradsolve = true;
        } else if (strcmp(argv[1], "-checkdefmapsolve") == 0) {
            this->m_RegFlags.checkdefmapsolve = true;
        } else if (strcmp(argv[1], "-analyze") == 0) {
            this->m_RegToolFlags.computeanalytics = true;
        } else if (strcmp(argv[1], "-problemid") == 0) {
            argc--; argv++;
            this->m_RegToolFlags.problemid = atoi(argv[1]);
        } else if (strcmp(argv[1], "-rscale") == 0) {
            argc--; argv++;
            this->m_ResamplingPara.gridscale = atof(argv[1]);
        } else if (strcmp(argv[1], "-nxr") == 0) {
            argc--; argv++;
            const std::string nxinput = argv[1];

            // strip the "x" in the string to get the numbers
            nxr = String2Vec(nxinput);

            if (nxr.size() == 1) {
                for(int i = 0; i < 3; ++i) {
                    this->m_ResamplingPara.nx[i] = static_cast<IntType>(nxr[0]);
                }
            } else if (nxr.size() == 3) {
                for(int i = 0; i < 3; ++i) {
                    this->m_ResamplingPara.nx[i] = static_cast<IntType>(nxr[i]);
                }
            } else {
                msg = "\n\x1b[31m error in grid size argument: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-resample") == 0) {
            this->m_RegToolFlags.resample = true;
        } else if (strcmp(argv[1], "-verbosity") == 0) {
            argc--; argv++;
            this->m_Verbosity = std::min(atoi(argv[1]),2);
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
                                      this->m_FFT.mpicomm, this->m_FFT.mpicommexists,
                                      this->m_Domain.nx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief default constructor
 *******************************************************************/
RegToolsOpt::~RegToolsOpt() {
    this->ClearMemory();
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode RegToolsOpt::ClearMemory() {
    PetscFunctionBegin;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief initialize class variables
 *******************************************************************/
PetscErrorCode RegToolsOpt::Initialize() {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr = this->SuperClass::Initialize(); CHKERRQ(ierr);

    this->m_RegToolFlags.readvecfield = false;
    this->m_RegToolFlags.readscafield = false;
    this->m_RegToolFlags.computedeffields = false;
    this->m_RegToolFlags.computegrad = false;
    this->m_RegToolFlags.tlabelmap = false;
    this->m_RegToolFlags.tscafield = false;
    this->m_RegToolFlags.computesynvel = false;
    this->m_RegToolFlags.resample = false;
    this->m_RegToolFlags.checkfwdsolveerr = false;
    this->m_RegToolFlags.checkfwdsolvetts = false;
    this->m_RegToolFlags.checkadjsolve = false;
    this->m_RegToolFlags.convert = false;
    this->m_RegToolFlags.checkdetdefgradsolve = false;
    this->m_RegToolFlags.computeerror = false;
    this->m_RegToolFlags.computeanalytics = false;
    this->m_RegToolFlags.computeresidual = false;
    this->m_RegToolFlags.problemid = 0;
    this->m_RegToolFlags.numrepeat = 1;

    this->m_ResamplingPara.gridscale = -1.0;
    this->m_ResamplingPara.nx[0] = -1.0;
    this->m_ResamplingPara.nx[1] = -1.0;
    this->m_ResamplingPara.nx[2] = -1.0;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief display usage message for binary
 *******************************************************************/
PetscErrorCode RegToolsOpt::Usage(bool advanced) {
    PetscErrorCode ierr = 0;
    int rank;
    std::string line;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    line = std::string(this->m_LineLength,'-');

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << line << std::endl;
        std::cout << " usage: regtools [options] " <<std::endl;
        std::cout << line << std::endl;
        std::cout << " where [options] is one or more of the following"<<std::endl;
        // ####################### advanced options #######################
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
        std::cout << " -mr <file>                  reference image (*.nii, *.nii.gz, *.hdr, *.nc)"<<std::endl;
        std::cout << " -mt <file>                  template image (*.nii, *.nii.gz, *.hdr, *.nc)"<<std::endl;
        std::cout << " -ivecx1 <file>              x1 component of vector field (*.nii, *.nii.gz, *.hdr, *.nc)"<<std::endl;
        std::cout << " -ivecx2 <file>              x2 component of vector field (*.nii, *.nii.gz, *.hdr, *.nc)"<<std::endl;
        std::cout << " -ivecx3 <file>              x3 component of vector field (*.nii, *.nii.gz, *.hdr, *.nc)"<<std::endl;
        std::cout << " -ifile <filename>           input file (scalar field/image)"<<std::endl;
        std::cout << " -i <path>                   input path (defines where registration results (i.e., velocity field, "<<std::endl;
        std::cout << "                             template image, and reference image) are stored; a prefix can be"<<std::endl;
        std::cout << "                             added by, e.g., doing '-i </path/prefix_>"<<std::endl;
        std::cout << " -x <path>                   output path (by default only deformed template image and velocity"<<std::endl;
        std::cout << "                             field will be written; for more output options, see flags;"<<std::endl;
        std::cout << "                             a prefix can be added by, e.g., doing '-x </path/prefix_>"<<std::endl;
        std::cout << " -nx <int>x<int>x<int>       grid size (e.g., 32x64x32); allows user to control grid size for synthetic" << std::endl;
        std::cout << "                             problems; assumed to be uniform if single integer is provided" << std::endl;
        std::cout << line << std::endl;
        std::cout << " ### postprocessing for registration (requires input fields and/or an input folder)"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -defgrad                    compute deformation gradient F = grad(inv(y)) (input: velocity field)"<<std::endl;
        std::cout << " -detdefgrad                 compute determinant of deformation gradient (input: velocity field)"<<std::endl;
        std::cout << " -invdetdefgrad              compute inverse of determinant of deformation gradient (input: velocity field)"<<std::endl;
        std::cout << " -deffield                   compute displacement field u (input: velocity field)"<<std::endl;
        std::cout << " -defmap                     compute deformation map y (input: velocity field)"<<std::endl;
        std::cout << " -tscafield                  transport scalar field (input: velocity field and scalar field)"<<std::endl;
        std::cout << " -tlabelmap                  transport label map (input: velocity field and scalar field)"<<std::endl;
        std::cout << " -residual                   compute residual between scalar fields ('-mr' and '-mt' options)"<<std::endl;
        std::cout << " -error                      compute error between scalar fields ('-mr' and '-mt' options)"<<std::endl;
        std::cout << " -analyze                    compute analytics for scalar field (-ifile option)"<<std::endl;
        // ####################### advanced options #######################
        if (advanced) {
        std::cout << " -detdefgradfromdeffield     compute gradient of some input scalar field ('-ifile' option)"<<std::endl;
        std::cout << " -grad                       compute gradient of some input scalar field ('-ifile' option)"<<std::endl;
        std::cout << " -xtimeseries                store time series (use with caution)"<<std::endl;
        std::cout << "                             problems; assumed to be uniform if single integer is provided"<<std::endl;
        }
        // ####################### advanced options #######################
        // ####################### advanced options #######################
        if (advanced) {
        std::cout << line << std::endl;
        std::cout << " -sigma <int>x<int>x<int>    size of gaussian smoothing kernel applied to input images (e.g., 1x2x1;"<<std::endl;
        std::cout << "                             units: voxel size; if only one parameter is set"<<std::endl;
        std::cout << "                             uniform smoothing is assumed: default: 1x1x1)"<<std::endl;
        std::cout << " -disablesmoothing           disable smoothing"<<std::endl;
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
        std::cout << " -interpolationorder <int>   order of interpolation model (default is 3)" << std::endl;
        }
        // ####################### advanced options #######################
        std::cout << line << std::endl;
        std::cout << " ### resampling"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -resample                   resample data (requires input scalar or vector field;"<<std::endl;
        std::cout << "                             output is input_resampled.ext)"<<std::endl;
        std::cout << " -rscale                     scale for resampling (multiplier applied to number of grid points)"<<std::endl;
        std::cout << " -nxr                        number of grid points for output"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " other parameters/debugging"<<std::endl;
        std::cout << line << std::endl;
        // ####################### advanced options #######################
        if (advanced) {
        std::cout << " -csynvel                    compute synthetic velocity field (use '-nx' to control size)"<<std::endl;
        std::cout << " -checkfwdsolveerr           check numerical error of forward solver"<<std::endl;
        std::cout << " -problemid <int>            problem id for error check"<<std::endl;
        std::cout << " -checkfwdsolvetts           check time-to-solution of forward solver"<<std::endl;
        std::cout << " -checkdetdefgradsolve       check solve for det(grad(y))"<<std::endl;
        std::cout << " -checkdefmapsolve           check solve for y"<<std::endl;
        }
        // ####################### advanced options #######################
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
PetscErrorCode RegToolsOpt::DisplayOptions() {
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
        std::cout << " Constrained Large Deformation Diffeomorphic Registration" << std::endl;
        std::cout << line << std::endl;
        std::cout << " Parallel Algorithms for Data Analysis and Simulation Group" << std::endl;
        std::cout << " The Institute of Computational Engineering and Sciences" << std::endl;
        std::cout << " The University of Texas at Austin" << std::endl;
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
                  << this->m_NumThreads << std::endl;
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
std::string RegToolsOpt::GetVecFieldFN(int i, int flag) {
    if (flag == 0) {
        if (i == 0) {
            return this->m_iVecFieldX1FN;
        } else if (i == 1) {
            return this->m_iVecFieldX2FN;
        } else if (i == 2) {
            return this->m_iVecFieldX3FN;
        } else return "";
    } else if (flag == 1) {
        if (i == 0) {
            return this->m_xVecFieldX1FN;
        } else if (i == 1) {
            return this->m_xVecFieldX2FN;
        } else if (i == 2) {
            return this->m_xVecFieldX3FN;
        } else return "";
    }
    return "";
}




/********************************************************************
 * @brief check the arguments set by user
 *******************************************************************/
std::string RegToolsOpt::GetScaFieldFN(int flag) {
    if (flag == 0) {
        return this->m_iScaFieldFN;
    } else if (flag == 1) {
        return this->m_xScaFieldFN;
    } else if (flag == 2) {
        return this->m_RFN;
    } else if (flag == 3) {
        return this->m_TFN;
    }
    return "";
}




/********************************************************************
 * @brief check the arguments set by user
 *******************************************************************/
PetscErrorCode RegToolsOpt::CheckArguments() {
    PetscErrorCode ierr;
    std::string msg, path, filename, extension;
    PetscFunctionBegin;

    // check output arguments
    if (   this->m_ReadWriteFlags.defgrad
        || this->m_ReadWriteFlags.defmap
        || this->m_ReadWriteFlags.detdefgrad
        || this->m_ReadWriteFlags.deffield ) {
        if (this->m_ReadWriteFlags.xfolder.empty()) {
            msg = "\x1b[31m output folder needs to be set (-x option) \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(); CHKERRQ(ierr);
        }

        if (this->m_ReadWriteFlags.ifolder.empty()) {
            msg = "\x1b[31m input folder needs to be set (-i option) \x1b[0m\n";
//            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
//            ierr = this->Usage(); CHKERRQ(ierr);
        }

        this->m_RegToolFlags.computedeffields = true;

        // set this flag to true, so that containers for reference and
        // template image are not to be deleted in registration class
        this->m_ReadWriteFlags.readfiles = true;
    }

    if ( !this->m_iVecFieldX1FN.empty()
      && !this->m_iVecFieldX2FN.empty()
      && !this->m_iVecFieldX3FN.empty() ) {
        this->m_RegToolFlags.readvecfield = true;
    }

    if (!this->m_iScaFieldFN.empty()) {
        this->m_RegToolFlags.readscafield = true;
    }

    if (this->m_RegToolFlags.resample) {
        if ( (this->m_ResamplingPara.gridscale == -1.0)
           && ( (this->m_ResamplingPara.nx[0] == -1.0)
             || (this->m_ResamplingPara.nx[1] == -1.0)
             || (this->m_ResamplingPara.nx[2] == -1.0))) {
            msg = "\x1b[31m number of grid points/scale for grid points needs to be set \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        if (!this->m_RegToolFlags.readvecfield && !this->m_RegToolFlags.readscafield) {
            msg = "\x1b[31m resampling requires input vector or scalar field \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        if (this->m_RegToolFlags.readvecfield) {
            ierr = GetFileName(path, filename, extension, this->m_iVecFieldX1FN); CHKERRQ(ierr);
            if (this->m_ReadWriteFlags.extension != ".nii.gz") {
                extension = this->m_ReadWriteFlags.extension;
            }
            this->m_xVecFieldX1FN = path + "/" + "resampled_" + filename + extension;

            ierr = GetFileName(path, filename, extension, this->m_iVecFieldX2FN); CHKERRQ(ierr);
            if (this->m_ReadWriteFlags.extension != ".nii.gz") {
                extension = this->m_ReadWriteFlags.extension;
            }
            this->m_xVecFieldX2FN = path + "/" + "resampled_" + filename + extension;

            ierr = GetFileName(path, filename, extension, this->m_iVecFieldX3FN); CHKERRQ(ierr);
            if (this->m_ReadWriteFlags.extension != ".nii.gz") {
                extension = this->m_ReadWriteFlags.extension;
            }
            this->m_xVecFieldX3FN = path + "/" + "resampled_" + filename + extension;
        }

        if (this->m_RegToolFlags.readscafield) {
            ierr = GetFileName(path, filename, extension, this->m_iScaFieldFN); CHKERRQ(ierr);
            if (this->m_ReadWriteFlags.extension != ".nii.gz") {
                extension = this->m_ReadWriteFlags.extension;
            }
            this->m_xScaFieldFN = path + "/" + "resampled_" + filename + extension;
        }
    }


    if (this->m_RegToolFlags.computegrad) {
        if ( !this->m_RegToolFlags.readvecfield && !this->m_RegToolFlags.readscafield ) {
            msg = "\x1b[31m computation of gradient requires input vector or scalar field \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        if (this->m_RegToolFlags.readvecfield) {
/*
            ierr = GetFileName(path,filename,extension,this->m_iVecFieldX1FN); CHKERRQ(ierr);
            if (this->m_ReadWriteFlags.extension != ".nii.gz") {
                extension = this->m_ReadWriteFlags.extension;
            }
            this->m_xVecFieldX1FN = path + "/" + filename + "-gradx1" + extension;

            ierr = GetFileName(path,filename,extension,this->m_iVecFieldX2FN); CHKERRQ(ierr);
            if (this->m_ReadWriteFlags.extension != ".nii.gz") {
                extension = this->m_ReadWriteFlags.extension;
            }
            this->m_xVecFieldX2FN = path + "/" + filename + "-gradx2" + extension;

            ierr = GetFileName(path,filename,extension,this->m_iVecFieldX3FN); CHKERRQ(ierr);
            if (this->m_ReadWriteFlags.extension != ".nii.gz") {
                extension = this->m_ReadWriteFlags.extension;
            }
            this->m_xVecFieldX3FN = path + "/" + filename + "-gradx3" + extension;
*/
        }

        if (this->m_RegToolFlags.readscafield) {
            ierr = GetFileName(path, filename, extension, this->m_iScaFieldFN); CHKERRQ(ierr);
            if (this->m_ReadWriteFlags.extension != ".nii.gz") {
                extension = this->m_ReadWriteFlags.extension;
            }
            this->m_xVecFieldX1FN = path + "/" + filename + "-gradx1" + extension;
            this->m_xVecFieldX2FN = path + "/" + filename + "-gradx2" + extension;
            this->m_xVecFieldX3FN = path + "/" + filename + "-gradx3" + extension;
       }
    }

    if (this->m_RegToolFlags.tscafield || this->m_RegToolFlags.tlabelmap) {
        // transport scalar field
        if (!this->m_RegToolFlags.readvecfield && !this->m_RegToolFlags.readscafield) {
            msg = "\x1b[31m solution of forward problem requires a velocity field and a scalar field \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        ierr = GetFileName(path, filename, extension, this->m_iScaFieldFN); CHKERRQ(ierr);
        if (this->m_ReadWriteFlags.extension != ".nii.gz") {
            extension = this->m_ReadWriteFlags.extension;
        }
        this->m_xScaFieldFN = path + "/" + filename + "-transported" + extension;
    }

    if (this->m_RegToolFlags.computeerror) {
        // transport scalar field
        if (this->m_TFN.empty() || this->m_RFN.empty()) {
            msg = "\x1b[31m reference and template images need to be set\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
    }

    if (this->m_RegToolFlags.computeanalytics) {
        if (this->m_iScaFieldFN.empty()) {
            msg = "\x1b[31m input scalarfield needs to be set\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
    }

    if (this->m_RegToolFlags.computeresidual) {
        // transport scalar field
        if (this->m_TFN.empty() || this->m_RFN.empty()) {
            msg = "\x1b[31m reference and template images need to be set\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        ierr = GetFileName(path, filename, extension, this->m_RFN); CHKERRQ(ierr);
        if (this->m_ReadWriteFlags.extension != ".nii.gz") {
            extension = this->m_ReadWriteFlags.extension;
        }
        this->m_xScaFieldFN = path + "/residual" + extension;
    }

    ierr = Assert(this->m_NumThreads > 0, "omp threads < 0"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




}  // namespace reg




#endif  // _REGTOOLSOPT_CPP_
