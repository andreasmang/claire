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

#ifndef _REGTOOLSOPT_CPP_
#define _REGTOOLSOPT_CPP_

#include "RegToolsOpt.hpp"
#include <algorithm>

#define _TO_STR(s) #s
#define TO_STR(s) _TO_STR(s)

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}


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
PetscErrorCode RegToolsOpt::ParseArguments(int argc, char** av) {
    PetscErrorCode ierr = 0;
    std::string msg;
    std::vector<int> nx;
    std::vector<int> nxr;
    std::vector<int> np;
    std::vector<int> sigma;
    std::vector<char*> av_list(av, av+argc);
    std::vector<char*> config_args;
    PetscFunctionBegin;

    if (argc == 1) {
        ierr = this->Usage(); CHKERRQ(ierr);
    }
    
    std::vector<char*>::iterator argv = av_list.begin();

    while(argc > 1) {
        if (   (strcmp(argv[1], "-h")    == 0)
            || (strcmp(argv[1], "-help") == 0)
            || (strcmp(argv[1], "-HELP") == 0) ) {
            ierr = this->Usage(); CHKERRQ(ierr);
        } else if (strcmp(argv[1], "-advanced") == 0) {
            ierr = this->Usage(true); CHKERRQ(ierr);
        } else if (strcmp(argv[1], "-config") == 0) {
          size_t pos = argv - av_list.begin();
          std::ifstream handle(argv[2]);
          std::string line;
          while (std::getline(handle, line)) {
            std::istringstream lines(line);
            std::string token;
            while(std::getline(lines, token, ' ')) {
              trim(token);
              if (token.size() > 0) {
                config_args.push_back(nullptr);
                config_args.back() = new char[token.size()+1];
                strcpy(config_args.back(), token.c_str());
                av_list.push_back(config_args.back());
                argc++;
              }
            }
          }
          handle.close();
          argv = av_list.begin() + pos;
          argc-=1; argv+=1;
        } else if (strcmp(argv[1], "-mr") == 0) {
            argc--; argv++;
            this->m_FileNames.mr.push_back(argv[1]);
        } else if (strcmp(argv[1], "-mt") == 0) {
            argc--; argv++;
            this->m_FileNames.mt.push_back(argv[1]);
        } else if (strcmp(argv[1], "-ifile") == 0) {
            argc--; argv++;
            this->m_FileNames.isc = argv[1];
        } else if (strcmp(argv[1], "-iref") == 0) {
            argc--; argv++;
            this->m_FileNames.iref = argv[1];
        } else if (strcmp(argv[1], "-ifilevec") == 0) {
            argc--; argv++;
            this->m_Domain.nc = static_cast<IntType>(atoi(argv[1]));
            if (this->m_Domain.nc < 1) {
                msg = "\n\x1b[31m number of components has to be larger than 1: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
            for (IntType i = 0; i < this->m_Domain.nc; ++i) {
                argc--; argv++;
                this->m_FileNames.ivec.push_back(argv[1]);
            }
        } else if (strcmp(argv[1], "-xfile") == 0) {
            argc--; argv++;
            this->m_FileNames.xsc = argv[1];
        } else if (strcmp(argv[1], "-v1") == 0) {
            argc--; argv++;
            this->m_FileNames.iv1 = argv[1];
        } else if (strcmp(argv[1], "-v2") == 0) {
            argc--; argv++;
            this->m_FileNames.iv2 = argv[1];
        } else if (strcmp(argv[1], "-v3") == 0) {
            argc--; argv++;
            this->m_FileNames.iv3 = argv[1];
        } else if (strcmp(argv[1], "-x") == 0) {
            argc--; argv++;
            this->m_FileNames.xfolder = argv[1];
        } else if (strcmp(argv[1], "-i") == 0) {
            argc--; argv++;
            this->m_FileNames.ifolder = argv[1];
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
//        } else if (strcmp(argv[1], "-nlabels") == 0) {
//            argc--; argv++;
//            this->m_NumLabels = static_cast<IntType>(atoi(argv[1]));
        } else if (strcmp(argv[1], "-labels") == 0) {
            argc--; argv++;
            const std::string labels = argv[1];
            const std::string sep = ",";
            // strip the "," in the string to get the numbers
            this->m_LabelIDs = String2Vec(labels, sep);
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
        } else if (strcmp(argv[1], "-iporder") == 0) {
            argc--; argv++;
            this->m_PDESolver.iporder = atoi(argv[1]);
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
                this->m_FileNames.extension = ".nii.gz";
            } else if (strcmp(argv[1], "2nc") == 0) {
                this->m_FileNames.extension = ".nc";
            }
            this->m_RegToolFlags.convert = true;
        } else if (strcmp(argv[1], "-usenc") == 0) {
            this->m_FileNames.extension = ".nc";
        } else if (strcmp(argv[1], "-velocity") == 0) {
            this->m_ReadWriteFlags.velocity = true;
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
        } else if (strcmp(argv[1], "-invdeffield") == 0) {
            this->m_ReadWriteFlags.invdeffield = true;
            this->m_ReadWriteFlags.deffield = true;
        } else if (strcmp(argv[1], "-timeseries") == 0) {
            this->m_ReadWriteFlags.timeseries = true;
        } else if (strcmp(argv[1], "-detdefgradfromdeffield") == 0) {
            this->m_RegFlags.detdefgradfromdeffield = true;
        } else if (strcmp(argv[1], "-residual") == 0) {
            this->m_RegToolFlags.computeresidual = true;
        } else if (strcmp(argv[1], "-dice") == 0) {
            this->m_RegToolFlags.computedice = true;
        } else if (strcmp(argv[1], "-r2t") == 0) {
            this->m_RegToolFlags.reference2template = true;
        } else if (strcmp(argv[1], "-error") == 0) {
            this->m_RegToolFlags.computeerror = true;
        } else if (strcmp(argv[1], "-deformimage") == 0) {
            this->m_RegToolFlags.deformimage = true;
        } else if (strcmp(argv[1], "-tprobmap") == 0) {
            this->m_RegToolFlags.tprobmaps = true;
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
        } else if (strcmp(argv[1], "-saveprob") == 0) {
            this->m_RegToolFlags.saveprob = true;
        } else if (strcmp(argv[1], "-computeravensmap") == 0) {
            this->m_RegToolFlags.computeravensmap = true;
        } else if (strcmp(argv[1], "-csynvel") == 0) {
            this->m_RegToolFlags.computesynvel = true;
        } else if (strcmp(argv[1], "-analyze") == 0) {
            this->m_RegToolFlags.computeanalytics = true;
        } else if (strcmp(argv[1], "-scale") == 0) {
            argc--; argv++;
            this->m_ResamplingPara.gridscale = atof(argv[1]);
        } else if (strcmp(argv[1], "-shift") == 0) {
            argc--; argv++;
            this->m_ResamplingPara.shift = atof(argv[1]);
        } else if (strcmp(argv[1], "-valuescale") == 0) {
            argc--; argv++;
            this->m_ResamplingPara.scale = atof(argv[1]);
        } else if (strcmp(argv[1], "-norm") == 0) {
            this->m_ResamplingPara.normalize = true;
        } else if (strcmp(argv[1], "-clip") == 0) {
            this->m_ResamplingPara.clip = true;
        } else if (strcmp(argv[1], "-nxnew") == 0) {
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
    
    for (size_t i = 0; i < config_args.size(); ++i) {
      delete[] config_args[i];
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
    this->m_RegToolFlags.tlabelmap = false;
    this->m_RegToolFlags.deformimage = false;
    this->m_RegToolFlags.computesynvel = false;
    this->m_RegToolFlags.resample = false;
    this->m_RegToolFlags.convert = false;
    this->m_RegToolFlags.computeerror = false;
    this->m_RegToolFlags.computeanalytics = false;
    this->m_RegToolFlags.computeresidual = false;
    this->m_RegToolFlags.reference2template = false;
    this->m_RegToolFlags.saveprob = false;
    this->m_RegToolFlags.computedice = false;
    this->m_RegToolFlags.tprobmaps = false;
    this->m_RegToolFlags.computeravensmap = false;
    
    
    this->m_ResamplingPara.normalize = false;
    this->m_ResamplingPara.clip = false;
    this->m_ResamplingPara.gridscale = -1.0;
    this->m_ResamplingPara.shift = 0.0;
    this->m_ResamplingPara.scale = 1.0;
    this->m_ResamplingPara.nx[0] = -1.0;
    this->m_ResamplingPara.nx[1] = -1.0;
    this->m_ResamplingPara.nx[2] = -1.0;

//    this->m_NumLabels = -1;

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
        std::cout << " usage: regtools [options] " << std::endl;
        std::cout << line << std::endl;
        std::cout << " where [options] is one or more of the following" << std::endl;
        // ####################### advanced options #######################
        if (advanced) {
        std::cout << line << std::endl;
        std::cout << " ### memory distribution and parallelism" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -nthreads <int>             number of threads (default: 1)" << std::endl;
        std::cout << " -np <int>x<int>             distribution of mpi tasks (cartesian grid) (example: -np 2x4 results" << std::endl;
        std::cout << "                             results in MPI distribution of size (nx1/2,nx2/4,nx3) for each mpi task)" << std::endl;
        }
        // ####################### advanced options #######################
        std::cout << line << std::endl;
        std::cout << " ### input parameters"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -mr <file>                  reference image (*.nii, *.nii.gz, *.hdr, *.nc)" << std::endl;
        std::cout << " -mt <file>                  template image (*.nii, *.nii.gz, *.hdr, *.nc)" << std::endl;
        std::cout << " -v1 <file>                  x1-component of vector field (*.nii, *.nii.gz, *.hdr, *.nc)" << std::endl;
        std::cout << " -v2 <file>                  x2-component of vector field (*.nii, *.nii.gz, *.hdr, *.nc)" << std::endl;
        std::cout << " -v3 <file>                  x3-component of vector field (*.nii, *.nii.gz, *.hdr, *.nc)" << std::endl;
        std::cout << " -ifile <filename>           input file (scalar field/image)" << std::endl;
        std::cout << " -iref <filename>            input reference file (scalar field/image)" << std::endl;
        std::cout << " -ifilevec <filenames>       list of input files  (fields/images)" << std::endl;
        std::cout << " -xfile <filename>           output file (scalar field/image)" << std::endl;
        std::cout << " -i <path>                   input path (defines where registration results (i.e., velocity field, " << std::endl;
        std::cout << "                             template image, and reference image) are stored; a prefix can be" << std::endl;
        std::cout << "                             added by, e.g., doing '-i </path/prefix_>" << std::endl;
        std::cout << " -x <path>                   output path (using this will generate default outputs)" << std::endl;
        std::cout << "                             a prefix can be added by, e.g., doing '-x </path/prefix_>" << std::endl;
        std::cout << " -nx <int>x<int>x<int>       grid size (e.g., 32x64x32); allows user to control grid size for synthetic" << std::endl;
        std::cout << "                             problems; assumed to be uniform if single integer is provided" << std::endl;
        std::cout << line << std::endl;
        std::cout << " ### postprocessing for registration (requires input fields and/or an input folder)" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -defgrad                    compute deformation gradient F = grad(inv(y)) (input: velocity field)" << std::endl;
        std::cout << " -detdefgrad                 compute determinant of deformation gradient (input: velocity field)" << std::endl;
        std::cout << " -invdetdefgrad              compute inverse of determinant of deformation gradient (input: velocity field)" << std::endl;
        std::cout << " -deffield                   compute displacement field u (input: velocity field)" << std::endl;
        std::cout << " -invdeffield                compute displacement field u in template image space (input: velocity field)" << std::endl;
        std::cout << " -defmap                     compute deformation map y (input: velocity field)" << std::endl;
        std::cout << " -residual                   compute residual between scalar fields ('-mr' and '-mt' options)" << std::endl;
        std::cout << " -error                      compute error between scalar fields ('-mr' and '-mt' options)" << std::endl;
        std::cout << " -analyze                    compute analytics for scalar field (-ifile option)" << std::endl;
        std::cout << " -deformimage                transport image (input: velocity field components and image" << std::endl;
        std::cout << "                             image to be deformed)" << std::endl;
        std::cout << line << std::endl;
        std::cout << " ### label maps" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -tlabelmap                  transport label map; the inputs are the components velocity field set via" << std::endl;
        std::cout << "                             '-v1', '-v2', -'v3' and the label map (a scalar field) set via '-ifile';" << std::endl;
        std::cout << "                             output file name (transported label map) needs to be set via '-xfile' option;" << std::endl;
        std::cout << "                             this will be a file that contains the individual label ids of the input;" << std::endl;
        std::cout << "                             the user needs to specify the labels to be transported via the '-labels' option" << std::endl;
        std::cout << "                             option (see below)" << std::endl;
        std::cout << " -labels <l1,l2,...>         labels to be transported (ids/numbers of labels)" << std::endl;
        std::cout << " -tprobmap                   transport multi-component probability maps" << std::endl;
        std::cout << " -saveprob                   enable this flag to write probability maps for individual labels to file" << std::endl;
        std::cout << " -r2t                        map (transport) from reference to template space by enabling this flag" << std::endl;
        // ####################### advanced options #######################
        if (advanced) {
        std::cout << " -detdefgradfromdeffield     compute gradient of some input scalar field ('-ifile' option)" << std::endl;
        std::cout << " -grad                       compute gradient of some input scalar field ('-ifile' option)" << std::endl;
        std::cout << " -xtimeseries                store time series (use with caution)" << std::endl;
        std::cout << "                             problems; assumed to be uniform if single integer is provided" << std::endl;
        }
        // ####################### advanced options #######################
        // ####################### advanced options #######################
        if (advanced) {
        std::cout << line << std::endl;
        std::cout << " -sigma <int>x<int>x<int>    size of gaussian smoothing kernel applied to input images (e.g., 1x2x1;" << std::endl;
        std::cout << "                             units: voxel size; if only one parameter is set" << std::endl;
        std::cout << "                             uniform smoothing is assumed: default: 1x1x1)" << std::endl;
        std::cout << " -disablesmoothing           disable smoothing" << std::endl;
        std::cout << line << std::endl;
        std::cout << " ### solver specific parameters (numerics)" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -pdesolver <type>           numerical time integrator for transport equations" << std::endl;
        std::cout << "                             <type> is one of the following" << std::endl;
        std::cout << "                                 sl           semi-Lagrangian method (default; unconditionally stable)" << std::endl;
        std::cout << "                                 rk2          rk2 time integrator (conditionally stable)" << std::endl;
        std::cout << " -convert <type>             convert data to type" << std::endl;
        std::cout << "                             <type> is one of the following" << std::endl;
        std::cout << "                                 2nii         convert to nifti" << std::endl;
        std::cout << "                                 2nc          convert to netcdf" << std::endl;
        std::cout << " -nt <int>                   number of time points (for time integration; default: 4)" << std::endl;
        std::cout << " -adapttimestep              vary number of time steps according to defined number" << std::endl;
        std::cout << " -cflnumber <dbl>            set cfl number" << std::endl;
        std::cout << " -interpolationorder <int>   order of interpolation model (default is 3)" << std::endl;
        }
        // ####################### advanced options #######################
        std::cout << line << std::endl;
        std::cout << " ### resampling" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -resample                   resample data (requires input scalar or vector field;" << std::endl;
        std::cout << "                             output is resampled_input.ext)" << std::endl;
        std::cout << " -scale                      scale for resampling (multiplier applied to number of grid points)" << std::endl;
        std::cout << " -nxnew                      number of grid points for output" << std::endl;
        std::cout << " -norm                       normalize output to [0, 1]" << std::endl;
        std::cout << " -clip                       clip values below 0" << std::endl;
        std::cout << " -shift                      shift values by value before norm/clip" << std::endl;
        std::cout << " -valuescale                 scale value after shifting, before norm/clip" << std::endl;
        std::cout << line << std::endl;
        std::cout << " ### other parameters/debugging" << std::endl;
        std::cout << line << std::endl;
        // ####################### advanced options #######################
        if (advanced) {
        std::cout << " -csynvel                    compute synthetic velocity field (use '-nx' to control size)" << std::endl;
        }
        // ####################### advanced options #######################
        std::cout << " -usenc                      use netcdf format os output (*.nc; default is *.nii.gz)" << std::endl;
        std::cout << " -verbosity <int>            verbosity level (ranges from 0 to 2; default: 1)" << std::endl;
        std::cout << " -help                       display a brief version of the user message" << std::endl;
        std::cout << " -advanced                   display this message" << std::endl;
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
std::string RegToolsOpt::GetVecFieldFN(int i, int flag) {
    if (flag == 0) {
        if (i == 0) {
            return this->m_FileNames.iv1;
        } else if (i == 1) {
            return this->m_FileNames.iv2;
        } else if (i == 2) {
            return this->m_FileNames.iv3;
        } else return "";
    } else if (flag == 1) {
        if (i == 0) {
            return this->m_FileNames.xv1;
        } else if (i == 1) {
            return this->m_FileNames.xv2;
        } else if (i == 2) {
            return this->m_FileNames.xv3;
        } else return "";
    }
    return "";
}




/********************************************************************
 * @brief check the arguments set by user
 *******************************************************************/
std::string RegToolsOpt::GetScaFieldFN(int flag) {
    if (flag == 0) {
        return this->m_FileNames.isc;
    } else if (flag == 1) {
        return this->m_FileNames.xsc;
    } else if (flag == 2) {
        return this->m_FileNames.mr[0];
    } else if (flag == 3) {
        return this->m_FileNames.mt[0];
    }
    return "";
}




/********************************************************************
 * @brief check the arguments set by user
 *******************************************************************/
PetscErrorCode RegToolsOpt::CheckArguments() {
    PetscErrorCode ierr;
    std::string msg, path, filename, extension;
    bool cdeffield;
    int flag;
    PetscFunctionBegin;


    if ( !this->m_FileNames.iv1.empty()
      && !this->m_FileNames.iv2.empty()
      && !this->m_FileNames.iv3.empty() ) {
        this->m_RegToolFlags.readvecfield = true;
    }

    if (!this->m_FileNames.isc.empty()) {
        this->m_RegToolFlags.readscafield = true;
    }

    cdeffield =  this->m_ReadWriteFlags.defgrad
              || this->m_ReadWriteFlags.detdefgrad
              || this->m_ReadWriteFlags.defmap
              || this->m_ReadWriteFlags.deffield;

    // check output arguments
    if (cdeffield) {
        if (!this->m_RegToolFlags.readvecfield) {
            msg =  "\n\x1b[31m computation of determinant of deformation\n";
            msg += " requires a velocity field as input:\n";
            msg += " use the -v1 <file> -v2 <file> -v3 <file> options\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
        if (this->m_FileNames.xfolder.empty()) {
            msg = "\x1b[31m define an output folder:\n";
            msg += " use the -x <folder> option\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
        this->m_RegToolFlags.computedeffields = true;

        // set this flag to true, so that containers for reference and
        // template image are not to be deleted in registration class
        this->m_ReadWriteFlags.readfiles = true;
    }

    if (this->m_RegToolFlags.resample) {
        if ( (this->m_ResamplingPara.gridscale == -1.0)
           && ( (this->m_ResamplingPara.nx[0] == -1.0)
             || (this->m_ResamplingPara.nx[1] == -1.0)
             || (this->m_ResamplingPara.nx[2] == -1.0)) ) {
            msg = "\x1b[31m number of grid points/scale for grid points needs to be set \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        flag = this->m_FileNames.isc.empty() ? 0 : 1;
        if (flag != 1) {
            flag = this->m_FileNames.iv1.empty()
                || this->m_FileNames.iv2.empty()
                || this->m_FileNames.iv3.empty() ? 0 : 2;
        } else {
          flag = this->m_FileNames.xsc.empty() ? 1 : -1;
        }


        // construct output name for resampling
        if (flag == 1) {
            ierr = GetFileName(path, filename, extension, this->m_FileNames.isc); CHKERRQ(ierr);
            if (this->m_FileNames.extension != ".nii.gz") {
                extension = this->m_FileNames.extension;
            }
            this->m_FileNames.xsc = path + "/" + "resampled_" + filename + extension;
        } else if (flag == 2) {
            ierr = GetFileName(path, filename, extension, this->m_FileNames.iv1); CHKERRQ(ierr);
            if (this->m_FileNames.extension != ".nii.gz") {
                extension = this->m_FileNames.extension;
            }
            this->m_FileNames.xv1 = path + "/" + "resampled_" + filename + extension;

            ierr = GetFileName(path, filename, extension, this->m_FileNames.iv2); CHKERRQ(ierr);
            if (this->m_FileNames.extension != ".nii.gz") {
                extension = this->m_FileNames.extension;
            }
            this->m_FileNames.xv2 = path + "/" + "resampled_" + filename + extension;

            ierr = GetFileName(path, filename, extension, this->m_FileNames.iv3); CHKERRQ(ierr);
            if (this->m_FileNames.extension != ".nii.gz") {
                extension = this->m_FileNames.extension;
            }
            this->m_FileNames.xv3 = path + "/" + "resampled_" + filename + extension;
        } else if (flag == -1) {
          // filename already set
        } else {
            msg = "\x1b[31m resampling requires input vector or scalar field \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
    }

    if (this->m_RegToolFlags.deformimage
        || this->m_RegToolFlags.tlabelmap
        || this->m_RegToolFlags.computeravensmap) {
        // check if flags are set correctly
        if (!this->m_RegToolFlags.readvecfield) {
            msg =  "\n\x1b[31m solution of forward problem requires a velocity field as input:\n";
            msg += " use the -v1 <file> -v2 <file> -v3 <file> options\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        if (!this->m_RegToolFlags.readscafield) {
            msg =  "\n\x1b[31m solution of forward problem requires scalar field to be transported as input:\n";
            msg += " use the -ifile <file> option\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        if (this->m_FileNames.xsc.empty() && !this->m_RegToolFlags.computedice) {
            msg = "\x1b[31m set file name for transported scalar field / ouptut:\n";
            msg += " use the -xfile <file> option\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

/*
        // construct output file name
        ierr = GetFileName(path, filename, extension, this->m_FileNames.isc); CHKERRQ(ierr);
        if (this->m_FileNames.extension != ".nii.gz") {
            extension = this->m_FileNames.extension;
        }
        if (this->m_FileNames.xfolder.empty()) {
            this->m_FileNames.xsc = path + "/" + filename + "-transported" + extension;
        } else {
            this->m_FileNames.xsc = this->m_FileNames.xfolder + "/" + filename + "-transported" + extension;
        }
*/
    }
    if (this->m_RegToolFlags.tlabelmap || this->m_RegToolFlags.computeravensmap) {
//        if (this->m_NumLabels == -1) {
//            msg = "\x1b[31m number of labels to be transported needs to be set\x1b[0m\n";
//            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
//            ierr = this->Usage(true); CHKERRQ(ierr);
//        }
        if (this->m_LabelIDs.size() == 0) {
            msg = "\x1b[31m user needs to set label ids ('-lables' option)\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
    }
    if (this->m_RegToolFlags.computeerror) {
        // transport scalar field
        if (this->m_FileNames.mt[0].empty() || this->m_FileNames.mr[0].empty()) {
            msg = "\x1b[31m reference and template images need to be set\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
    }

    if (this->m_RegToolFlags.computeanalytics) {
        if (this->m_FileNames.isc.empty()) {
            msg = "\x1b[31m input scalarfield needs to be set\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
    }

    if (this->m_RegToolFlags.computeresidual) {
        // transport scalar field
        if (this->m_FileNames.mt[0].empty() || this->m_FileNames.mr[0].empty()) {
            msg = "\x1b[31m reference and template images need to be set\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }

        ierr = GetFileName(path, filename, extension, this->m_FileNames.mr[0]); CHKERRQ(ierr);
        if (this->m_FileNames.extension != ".nii.gz") {
            extension = this->m_FileNames.extension;
        }
        this->m_FileNames.xsc = path + "/residual" + extension;
    }

    //ierr = Assert(this->m_NumThreads > 0, "omp threads < 0"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




}  // namespace reg




#endif  // _REGTOOLSOPT_CPP_
