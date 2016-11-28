/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGOPT_CPP_
#define _REGOPT_CPP_

#include <string>
#include <vector>

#include "RegOpt.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegOpt"
RegOpt::RegOpt() {
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegOpt"
RegOpt::RegOpt(int argc, char** argv) {
    this->Initialize();
    this->ParseArguments(argc, argv);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegOpt"
RegOpt::RegOpt(const RegOpt& opt) {
    this->Initialize();
    this->Copy(opt);
}




/********************************************************************
 * @brief copy entries of input options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Copy"
void RegOpt::Copy(const RegOpt& opt) {
    this->m_SetupDone = false;
    this->m_StoreCheckPoints = opt.m_StoreCheckPoints;
    this->m_FFT.plan = NULL;
    this->m_FFT.mpicomm = NULL;

    this->m_FFT.osize[0] = opt.m_FFT.osize[0];
    this->m_FFT.osize[1] = opt.m_FFT.osize[1];
    this->m_FFT.osize[2] = opt.m_FFT.osize[2];
    this->m_FFT.ostart[0] = opt.m_FFT.ostart[0];
    this->m_FFT.ostart[1] = opt.m_FFT.ostart[1];
    this->m_FFT.ostart[2] = opt.m_FFT.ostart[2];

    this->m_Domain.nlocal = opt.m_Domain.nlocal;
    this->m_Domain.nglobal = opt.m_Domain.nglobal;
    this->m_Domain.isize[0] = opt.m_Domain.isize[0];
    this->m_Domain.isize[1] = opt.m_Domain.isize[1];
    this->m_Domain.isize[2] = opt.m_Domain.isize[2];
    this->m_Domain.istart[0] = opt.m_Domain.istart[0];
    this->m_Domain.istart[1] = opt.m_Domain.istart[1];
    this->m_Domain.istart[2] = opt.m_Domain.istart[2];

    this->m_Domain.nx[0] = opt.m_Domain.nx[0];
    this->m_Domain.nx[1] = opt.m_Domain.nx[0];
    this->m_Domain.nx[2] = opt.m_Domain.nx[0];

    this->m_Domain.nt = opt.m_Domain.nt;
    this->m_Domain.nc = opt.m_Domain.nc;

    this->m_Domain.timehorizon[0] = opt.m_Domain.timehorizon[0];
    this->m_Domain.timehorizon[1] = opt.m_Domain.timehorizon[1];

    this->m_RegNorm.type = opt.m_RegNorm.type;

    this->m_RegNorm.beta[0] = opt.m_RegNorm.beta[0];
    this->m_RegNorm.beta[1] = opt.m_RegNorm.beta[1];
    this->m_RegNorm.beta[2] = opt.m_RegNorm.beta[2];

    this->m_PDESolver = opt.m_PDESolver;
    this->m_RegModel = opt.m_RegModel;

    // smoothing
    this->m_Sigma[0] = opt.m_Sigma[0];
    this->m_Sigma[1] = opt.m_Sigma[1];
    this->m_Sigma[2] = opt.m_Sigma[2];

    this->m_KrylovSolverPara.tol[0] = opt.m_KrylovSolverPara.tol[0];
    this->m_KrylovSolverPara.tol[1] = opt.m_KrylovSolverPara.tol[1];
    this->m_KrylovSolverPara.tol[2] = opt.m_KrylovSolverPara.tol[2];
    this->m_KrylovSolverPara.pcmaxit = opt.m_KrylovSolverPara.pcmaxit;
    this->m_KrylovSolverPara.reesteigvals = opt.m_KrylovSolverPara.reesteigvals;
    this->m_KrylovSolverPara.usepetsceigest = opt.m_KrylovSolverPara.usepetsceigest;
    this->m_KrylovSolverPara.pctol[0] = opt.m_KrylovSolverPara.pctol[0];
    this->m_KrylovSolverPara.pctol[1] = opt.m_KrylovSolverPara.pctol[1];
    this->m_KrylovSolverPara.pctol[2] = opt.m_KrylovSolverPara.pctol[2];
    this->m_KrylovSolverPara.matvectype = opt.m_KrylovSolverPara.matvectype;
    this->m_KrylovSolverPara.checkhesssymmetry = opt.m_KrylovSolverPara.checkhesssymmetry;

    this->m_OptPara.maxit = opt.m_OptPara.maxit;
    this->m_OptPara.minit = opt.m_OptPara.minit;
    this->m_OptPara.tol[0] = opt.m_OptPara.tol[0];
    this->m_OptPara.tol[1] = opt.m_OptPara.tol[1];
    this->m_OptPara.tol[2] = opt.m_OptPara.tol[2];
    this->m_OptPara.presolvemaxit = opt.m_OptPara.presolvemaxit;
    this->m_OptPara.presolvetol[0] = opt.m_OptPara.presolvetol[0];
    this->m_OptPara.presolvetol[1] = opt.m_OptPara.presolvetol[1];
    this->m_OptPara.presolvetol[2] = opt.m_OptPara.presolvetol[2];
    this->m_OptPara.fastsolve = opt.m_OptPara.fastsolve;
    this->m_OptPara.fastpresolve = opt.m_OptPara.fastpresolve;
    this->m_OptPara.method = opt.m_OptPara.method;
    this->m_OptPara.usezeroinitialguess = opt.m_OptPara.usezeroinitialguess;

    this->m_SolveType = opt.m_SolveType;

    // flags
    this->m_ReadWriteFlags.readfiles = opt.m_ReadWriteFlags.readfiles;
    this->m_ReadWriteFlags.readvelocity = opt.m_ReadWriteFlags.readvelocity;

    this->m_ReadWriteFlags.templateim = opt.m_ReadWriteFlags.templateim;
    this->m_ReadWriteFlags.referenceim = opt.m_ReadWriteFlags.referenceim;
    this->m_ReadWriteFlags.timeseries = opt.m_ReadWriteFlags.timeseries;
    this->m_ReadWriteFlags.iterates = opt.m_ReadWriteFlags.iterates;
    this->m_ReadWriteFlags.results = opt.m_ReadWriteFlags.results;
    this->m_ReadWriteFlags.defmap = opt.m_ReadWriteFlags.defmap;
    this->m_ReadWriteFlags.deffield = opt.m_ReadWriteFlags.deffield;
    this->m_ReadWriteFlags.defgrad = opt.m_ReadWriteFlags.defgrad;
    this->m_ReadWriteFlags.detdefgrad = opt.m_ReadWriteFlags.detdefgrad;
    this->m_ReadWriteFlags.residual = opt.m_ReadWriteFlags.residual;
    this->m_ReadWriteFlags.velnorm = opt.m_ReadWriteFlags.velnorm;
    this->m_ReadWriteFlags.deftemplate = opt.m_ReadWriteFlags.deftemplate;

    this->m_RegFlags.applysmoothing = opt.m_RegFlags.applysmoothing;
    this->m_RegFlags.applyrescaling = opt.m_RegFlags.applyrescaling;
    this->m_RegFlags.detdefgradfromdeffield = opt.m_RegFlags.detdefgradfromdeffield;
    this->m_RegFlags.invdefgrad = opt.m_RegFlags.invdefgrad;

    // parameter continuation
    this->m_ParaCont.strategy = opt.m_ParaCont.strategy;
    this->m_ParaCont.enabled = opt.m_ParaCont.enabled;
    this->m_ParaCont.targetbeta = opt.m_ParaCont.targetbeta;
    this->m_ParaCont.beta0 = opt.m_ParaCont.beta0;

    // grid continuation
    this->m_GridCont.enabled = opt.m_GridCont.enabled;
    this->m_GridCont.nlevels = opt.m_GridCont.nlevels;

    // scale continuation
    this->m_ScaleCont.enabled = opt.m_ScaleCont.enabled;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 6; ++j) {
            this->m_ScaleCont.sigma[i][j] = opt.m_ScaleCont.sigma[i][j];
        }
    }

    // monitor for registration
    this->m_RegMonitor.JAC = opt.m_RegMonitor.JAC;
    this->m_RegMonitor.CFL = opt.m_RegMonitor.CFL;
    this->m_RegMonitor.jacmin = opt.m_RegMonitor.jacmin;
    this->m_RegMonitor.jacmax = opt.m_RegMonitor.jacmax;
    this->m_RegMonitor.jacmean = opt.m_RegMonitor.jacmean;
    this->m_RegMonitor.jacbound = opt.m_RegMonitor.jacbound;

    for (int i = 0; i < NLOGFLAGS; ++i) {
        this->m_Log.enabled[i] = opt.m_Log.enabled[i];
    }
    this->m_Log.residual[0] = 0;
    this->m_Log.residual[1] = 0;
    this->m_Log.residual[2] = 0;
    this->m_Log.residual[3] = 0;

    this->m_NumThreads = opt.m_NumThreads;
    this->m_FFT.mpicomm = opt.m_FFT.mpicomm;

    this->m_CartGridDims[0] = opt.m_CartGridDims[0];
    this->m_CartGridDims[1] = opt.m_CartGridDims[1];

    this->m_Verbosity = opt.m_Verbosity;
    this->m_Indent = opt.m_Indent;
}




/********************************************************************
 * @brief parse user arguments
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ParseArguments"
PetscErrorCode RegOpt::ParseArguments(int argc, char** argv) {
    PetscErrorCode ierr = 0;
    std::string msg;
    std::vector<unsigned int> nx;
    std::vector<unsigned int> np;
    std::vector<unsigned int> sigma;
    PetscFunctionBegin;

    while (argc > 1) {
        if (   (strcmp(argv[1], "-help") == 0)
            || (strcmp(argv[1], "-h")    == 0)
            || (strcmp(argv[1], "-HELP") == 0) ) {
            ierr = this->Usage(); CHKERRQ(ierr);
        } else if (strcmp(argv[1], "-advanced") == 0) {
            ierr = this->Usage(true); CHKERRQ(ierr);
        } else if (strcmp(argv[1], "-nx") == 0) {
            argc--; argv++;
            const std::string nxinput = argv[1];

            // strip the "x" in the string to get the numbers
            nx = String2Vec(nxinput);

            if (nx.size() == 1) {
                for (int i = 0; i < 3; ++i) {
                    this->m_Domain.nx[i] = static_cast<IntType>(nx[0]);
                }
            } else if (nx.size() == 3) {
                for (int i = 0; i < 3; ++i) {
                    this->m_Domain.nx[i] = static_cast<IntType>(nx[i]);
                }
            } else {
                msg = "\n\x1b[31m error in grid size argument: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-nt") == 0) {
            argc--; argv++;
            this->m_Domain.nt = static_cast<IntType>(atoi(argv[1]));
        } else if (strcmp(argv[1], "-sigma") == 0) {
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
        } else if (strcmp(argv[1], "-disablesmoothing") == 0) {
            this->m_RegFlags.applysmoothing = false;
        } else if (strcmp(argv[1], "-disablerescaling") == 0) {
            this->m_RegFlags.applyrescaling = false;
        } else if (strcmp(argv[1], "-nthreads") == 0) {
            argc--; argv++;
            this->m_NumThreads = atoi(argv[1]);
        } else if (strcmp(argv[1], "-np") == 0) {
            argc--; argv++;
            const std::string npinput = argv[1];

            // strip the "x" in the string to get the numbers
            np = String2Vec(npinput);

            if (np.size() == 1) {
                for (int i = 0; i < 2; ++i) {
                    this->m_CartGridDims[i] = static_cast<int>(np[0]);
                }
            } else if (np.size() == 2) {
                for (int i = 0; i < 2; ++i) {
                    this->m_CartGridDims[i] = static_cast<int>(np[i]);
                }
            } else {
                msg = "\n\x1b[31m error in number of procs: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-mr") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.mr = argv[1];
        } else if (strcmp(argv[1], "-mt") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.mt = argv[1];
        } else if (strcmp(argv[1], "-vx1") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.vx1 = argv[1];
        } else if (strcmp(argv[1], "-vx2") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.vx2 = argv[1];
        } else if (strcmp(argv[1], "-vx3") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.vx3 = argv[1];
        } else if (strcmp(argv[1], "-x") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.xfolder = argv[1];
        } else if (strcmp(argv[1], "-i") == 0) {
            argc--; argv++;
            this->m_ReadWriteFlags.ifolder = argv[1];
        } else if (strcmp(argv[1], "-usenc") == 0) {
            this->m_ReadWriteFlags.extension = ".nc";
        } else if (strcmp(argv[1], "-usebin") == 0) {
            this->m_ReadWriteFlags.extension = ".bin";
        } else if (strcmp(argv[1], "-usehdf5") == 0) {
            this->m_ReadWriteFlags.extension = ".hdf5";
        } else if (strcmp(argv[1], "-xresult") == 0) {
            this->m_ReadWriteFlags.results = true;
        } else if (strcmp(argv[1], "-xdeftemplate") == 0) {
            this->m_ReadWriteFlags.deftemplate = true;
        } else if (strcmp(argv[1], "-xmt") == 0) {
            this->m_ReadWriteFlags.templateim = true;
        } else if (strcmp(argv[1], "-xmr") == 0) {
            this->m_ReadWriteFlags.referenceim = true;
        } else if (strcmp(argv[1], "-xresidual") == 0) {
            this->m_ReadWriteFlags.residual = true;
        } else if (strcmp(argv[1], "-xvelnorm") == 0) {
            this->m_ReadWriteFlags.velnorm = true;
        } else if (strcmp(argv[1], "-xdefgrad") == 0) {
            this->m_ReadWriteFlags.defgrad = true;
        } else if (strcmp(argv[1], "-xdetdefgrad") == 0) {
            this->m_ReadWriteFlags.detdefgrad = true;
        } else if (strcmp(argv[1], "-xdeffield") == 0) {
            this->m_ReadWriteFlags.deffield = true;
        } else if (strcmp(argv[1], "-xdefmap") == 0) {
            this->m_ReadWriteFlags.defmap = true;
        } else if (strcmp(argv[1], "-xiterates") == 0) {
            this->m_ReadWriteFlags.iterates = true;
        } else if (strcmp(argv[1], "-xtimeseries") == 0) {
            this->m_ReadWriteFlags.timeseries = true;
        } else if (strcmp(argv[1], "-storecheckpoints") == 0) {
            this->m_StoreCheckPoints = true;
        } else if (strcmp(argv[1], "-logjacobian") == 0) {
            this->m_Log.enabled[LOGJAC] = true;
        } else if (strcmp(argv[1], "-logkrylovres") == 0) {
            this->m_Log.enabled[LOGKSPRES] = true;
        } else if (strcmp(argv[1], "-logworkload") == 0) {
            this->m_Log.enabled[LOGLOAD] = true;
        } else if (strcmp(argv[1], "-logresidual") == 0) {
            this->m_Log.enabled[LOGRES] = true;
        } else if (strcmp(argv[1], "-detdefgradfromdeffield") == 0) {
            this->m_RegFlags.detdefgradfromdeffield = true;
        } else if (strcmp(argv[1], "-preset") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "fast-aggressive") == 0) {
                this->m_SolveType = FAST_AGG;
            } else if (strcmp(argv[1], "fast-smooth") == 0) {
                this->m_SolveType = FAST_SMOOTH;
            } else if (strcmp(argv[1], "accurate-aggressive") == 0) {
                this->m_SolveType = ACC_AGG;
            } else if (strcmp(argv[1], "accurate-smooth") == 0) {
                this->m_SolveType = ACC_SMOOTH;
            } else {
                msg = "\n\x1b[31m high level solver flag not available: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-fastsolve") == 0) {
            this->m_OptPara.fastsolve = true;
        } else if (strcmp(argv[1], "-ic") == 0) {
            this->m_RegModel = STOKES;
        } else if (strcmp(argv[1], "-ric") == 0) {
            this->m_RegModel = RELAXEDSTOKES;
        } else if (strcmp(argv[1], "-optmeth") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "fn") == 0) {
                this->m_OptPara.method = FULLNEWTON;
            } else if (strcmp(argv[1], "gn") == 0) {
                this->m_OptPara.method = GAUSSNEWTON;
            } else {
                msg = "\n\x1b[31m optimization method not defined: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-maxit") == 0) {
            argc--; argv++;
            this->m_OptPara.maxit = atoi(argv[1]);
        } else if (strcmp(argv[1], "-gabs") == 0) {
            argc--; argv++;
            this->m_OptPara.tol[0] = atof(argv[1]);
        } else if (strcmp(argv[1], "-grel") == 0) {
            argc--; argv++;
            this->m_OptPara.tol[2] = atof(argv[1]);
        } else if (strcmp(argv[1], "-nonzeroinitialguess") == 0) {
            this->m_OptPara.usezeroinitialguess = false;
        } else if (strcmp(argv[1], "-jbound") == 0) {
            argc--; argv++;
            this->m_RegMonitor.jacbound = atof(argv[1]);
        } else if (strcmp(argv[1], "-krylovsolver") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "pcg") == 0) {
                this->m_KrylovSolverPara.solver = PCG;
                this->m_KrylovSolverPara.name = "PCG";
            } else if (strcmp(argv[1], "gmres") == 0) {
                this->m_KrylovSolverPara.solver = GMRES;
                this->m_KrylovSolverPara.name = "GMRES";
            } else {
                msg = "\n\x1b[31m optimization method not defined: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-krylovmaxit") == 0) {
            argc--; argv++;
            this->m_KrylovSolverPara.maxit = atoi(argv[1]);
        } else if (strcmp(argv[1], "-krylovtol") == 0) {
            argc--; argv++;
            this->m_KrylovSolverPara.tol[0] = atof(argv[1]);
        } else if (strcmp(argv[1], "-krylovfseq") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "none") == 0) {
                this->m_KrylovSolverPara.fseqtype = NOFS;
            } else if (strcmp(argv[1], "quadratic") == 0) {
                this->m_KrylovSolverPara.fseqtype = QDFS;
            } else if (strcmp(argv[1], "suplinear") == 0) {
                this->m_KrylovSolverPara.fseqtype = SLFS;
            } else {
                msg = "\n\x1b[31m forcing sequence not defined: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-precond") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "none") == 0) {
                this->m_KrylovSolverPara.pctype = NOPC;
                this->m_KrylovSolverPara.matvectype = DEFAULTMATVEC;
            } else if (strcmp(argv[1], "invreg") == 0) {
                this->m_KrylovSolverPara.pctype = INVREG;
                this->m_KrylovSolverPara.matvectype = DEFAULTMATVEC;
//                this->m_KrylovSolverPara.pctype = NOPC;
//                this->m_KrylovSolverPara.matvectype = PRECONDMATVECSYM;
            } else if (strcmp(argv[1], "2level") == 0) {
                this->m_KrylovSolverPara.pctype = TWOLEVEL;
                this->m_KrylovSolverPara.matvectype = PRECONDMATVECSYM;
//                 this->m_KrylovSolverPara.matvectype = PRECONDMATVEC;
            } else {
                msg = "\n\x1b[31m preconditioner not defined: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-pcsolver") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "pcg") == 0) {
                this->m_KrylovSolverPara.pcsolver = PCG;
                //this->m_KrylovSolverPara.matvectype = PRECONDMATVECSYM;
            } else if (strcmp(argv[1], "fpcg") == 0) {
                this->m_KrylovSolverPara.pcsolver = FCG;
                //this->m_KrylovSolverPara.matvectype = PRECONDMATVECSYM;
            } else if (strcmp(argv[1], "gmres") == 0) {
                this->m_KrylovSolverPara.pcsolver = GMRES;
                this->m_KrylovSolverPara.name = "GMRES";
                //this->m_KrylovSolverPara.matvectype = PRECONDMATVEC;
            } else if (strcmp(argv[1], "fgmres") == 0) {
                this->m_KrylovSolverPara.pcsolver = FGMRES;
                this->m_KrylovSolverPara.name = "FGMRES";
                //this->m_KrylovSolverPara.matvectype = PRECONDMATVEC;
            } else if (strcmp(argv[1], "cheb") == 0) {
                this->m_KrylovSolverPara.pcsolver = CHEB;
                this->m_KrylovSolverPara.name = "CHEB";
                //this->m_KrylovSolverPara.matvectype = PRECONDMATVECSYM;
            } else {
                msg = "\n\x1b[31m optimization method not defined: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-reesteigvals") == 0) {
            this->m_KrylovSolverPara.reesteigvals = true;
        } else if (strcmp(argv[1], "-pctolscale") == 0) {
            argc--; argv++;
            this->m_KrylovSolverPara.pctolscale = atof(argv[1]);
        } else if (strcmp(argv[1], "-pcsolvermaxit") == 0) {
            argc--; argv++;
            this->m_KrylovSolverPara.pcmaxit = atoi(argv[1]);
        } else if (strcmp(argv[1], "-checksymmetry") == 0) {
            this->m_KrylovSolverPara.checkhesssymmetry = true;
        } else if (strcmp(argv[1], "-pdesolver") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "rk2") == 0) {
                this->m_PDESolver.type = RK2;
            } else if (strcmp(argv[1], "sl") == 0) {
                this->m_PDESolver.type = SL;
            } else {
                msg = "\n\x1b[31m pde solver not implemented: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-regnorm") == 0) {
            argc--; argv++;
            if (strcmp(argv[1], "h1s") == 0) {
                this->m_RegNorm.type = H1SN;
            } else if (strcmp(argv[1], "h2s") == 0) {
                this->m_RegNorm.type = H2SN;
            } else if (strcmp(argv[1], "h3s") == 0) {
                this->m_RegNorm.type = H3SN;
            } else if (strcmp(argv[1], "h1") == 0) {
                this->m_RegNorm.type = H1;
            } else if (strcmp(argv[1], "h2") == 0) {
                this->m_RegNorm.type = H2;
            } else if (strcmp(argv[1], "h3") == 0) {
                this->m_RegNorm.type = H3;
            } else if (strcmp(argv[1], "l2") == 0) {
                this->m_RegNorm.type = L2;
            } else {
                msg = "\n\x1b[31m regularization norm not available: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-betav") == 0) {
            argc--; argv++;
            this->m_RegNorm.beta[0] = atof(argv[1]);
            this->m_RegNorm.beta[1] = atof(argv[1]);
        } else if (strcmp(argv[1], "-betaw") == 0) {
            argc--; argv++;
            this->m_RegNorm.beta[2] = atof(argv[1]);
        } else if (strcmp(argv[1], "-train") == 0) {
            if (this->m_ParaCont.enabled) {
                msg = "\n\x1b[31m you can't do training and continuation simultaneously\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(); CHKERRQ(ierr);
            }

            argc--; argv++;
            if (strcmp(argv[1], "binary") == 0) {
                this->m_ParaCont.strategy = PCONTBINSEARCH;
                this->m_ParaCont.enabled = true;
            } else if (strcmp(argv[1], "reduce") == 0) {
                this->m_ParaCont.strategy = PCONTREDUCESEARCH;
                this->m_ParaCont.enabled = true;
            } else {
                msg = "\n\x1b[31m training method not implemented: %s\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(); CHKERRQ(ierr);
            }
        } else if (strcmp(argv[1], "-betavcont") == 0) {
           if (this->m_ParaCont.enabled) {
                msg = "\n\x1b[31m you can't do training and continuation simultaneously\x1b[0m\n";
                ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
                ierr = this->Usage(true); CHKERRQ(ierr);
            }

            this->m_ParaCont.strategy = PCONTINUATION;
            this->m_ParaCont.enabled = true;

            argc--; argv++;
            this->m_ParaCont.targetbeta = atof(argv[1]);
        } else if (strcmp(argv[1], "-betavinit") == 0) {
            argc--; argv++;
            this->m_ParaCont.beta0 = atof(argv[1]);
        } else if (strcmp(argv[1], "-scalecont") == 0) {
            this->m_ScaleCont.enabled = true;
        } else if (strcmp(argv[1], "-gridcont") == 0) {
            this->m_GridCont.enabled = true;
        } else if (strcmp(argv[1], "-verbosity") == 0) {
            argc--; argv++;
            this->m_Verbosity = atoi(argv[1]);
        } else if (strcmp(argv[1], "-mdefgrad") == 0) {
            this->m_RegMonitor.JAC = true;
        } else if (strcmp(argv[1], "-invdefgrad") == 0) {
            this->m_RegFlags.invdefgrad = true;
        } else {
            msg = "\n\x1b[31m argument not valid: %s\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str(), argv[1]); CHKERRQ(ierr);
            ierr = this->Usage(true); CHKERRQ(ierr);
        }
        argc--; argv++;
    }

    // check the arguments/parameters set by the user
    ierr = this->CheckArguments(); CHKERRQ(ierr);

    if (this->m_SolveType !=  NOTSET) {
        ierr = this->SetPresetParameters(); CHKERRQ(ierr);
    } else {
        if (this->m_OptPara.fastsolve) {
            if (this->m_RegNorm.type == H2SN || this->m_RegNorm.type == H2) {
                this->m_SolveType = FAST_SMOOTH;
            }
            if (this->m_RegNorm.type == H1SN || this->m_RegNorm.type == H1) {
                this->m_SolveType = FAST_AGG;
            }
            ierr = this->SetPresetParameters(); CHKERRQ(ierr);
        }
    }

    // set number of threads
    ierr = InitializeDataDistribution(this->m_NumThreads,
                                      this->m_CartGridDims,
                                      this->m_FFT.mpicomm); CHKERRQ(ierr);
    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegOpt"
RegOpt::~RegOpt() {
    this->ClearMemory();
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode RegOpt::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    ierr = this->DestroyFFT(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief clean up fourier transform and mpi communicator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DestroyFFT"
PetscErrorCode RegOpt::DestroyFFT() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_FFT.plan != NULL) {
        accfft_destroy_plan(this->m_FFT.plan);
        accfft_cleanup();
        this->m_FFT.plan = NULL;
    }

    if (this->m_FFT.mpicomm != NULL) {
        MPI_Comm_free(&this->m_FFT.mpicomm);
        this->m_FFT.mpicomm = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief setup fourier transform and mpi communicator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "InitializeFFT"
PetscErrorCode RegOpt::InitializeFFT() {
    PetscErrorCode ierr = 0;
    int nx[3], isize[3], istart[3], osize[3], ostart[3];
    int nalloc;
    ScalarType *u = NULL, fftsetuptime;
    Complex *uk = NULL;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    // if communicator is not set up
    if (this->m_FFT.mpicomm == NULL) {
        ierr = InitializeDataDistribution(this->m_NumThreads,
                                          this->m_CartGridDims,
                                          this->m_FFT.mpicomm); CHKERRQ(ierr);
    }

    // parse grid size for setup
    for (int i = 0; i < 3; ++i) {
        nx[i] = static_cast<int>(this->m_Domain.nx[i]);
        this->m_Domain.hx[i] = PETSC_PI*2.0/static_cast<ScalarType>(nx[i]);
    }

    // get sizes
    nalloc = accfft_local_size_dft_r2c(nx, isize, istart, osize, ostart, this->m_FFT.mpicomm);
    ierr = Assert(nalloc != 0, "alloc problem"); CHKERRQ(ierr);
    this->m_FFT.nalloc = static_cast<IntType>(nalloc);

    // set up the fft
    u = reinterpret_cast<ScalarType*>(accfft_alloc(nalloc));
    ierr = Assert(u != NULL, "alloc failed"); CHKERRQ(ierr);

    uk = reinterpret_cast<Complex*>(accfft_alloc(nalloc));
    ierr = Assert(uk != NULL, "alloc failed"); CHKERRQ(ierr);

    if (this->m_FFT.plan != NULL) {
        accfft_destroy_plan(this->m_FFT.plan);
        this->m_FFT.plan = NULL;
    }

    fftsetuptime = -MPI_Wtime();
    this->m_FFT.plan = accfft_plan_dft_3d_r2c(nx, u, reinterpret_cast<double*>(uk),
                                              this->m_FFT.mpicomm, ACCFFT_MEASURE);
    fftsetuptime += MPI_Wtime();

    // set the fft setup time
    this->m_Timer[FFTSETUP][LOG] += fftsetuptime;

    // compute global and local size
    this->m_Domain.nlocal = 1;
    this->m_Domain.nglobal = 1;
    for (int i = 0; i < 3; ++i) {
        this->m_Domain.isize[i] = static_cast<IntType>(isize[i]);
        this->m_Domain.istart[i] = static_cast<IntType>(istart[i]);

        this->m_FFT.osize[i] = static_cast<IntType>(osize[i]);
        this->m_FFT.ostart[i] = static_cast<IntType>(ostart[i]);

        this->m_Domain.nlocal *= static_cast<IntType>(isize[i]);
        this->m_Domain.nglobal *= this->m_Domain.nx[i];
    }

    // check if sizes are ok
    ierr = reg::Assert(this->m_Domain.nlocal > 0, "bug in setup"); CHKERRQ(ierr);
    ierr = reg::Assert(this->m_Domain.nglobal > 0, "bug in setup"); CHKERRQ(ierr);

    // clean up
    if (u != NULL) { accfft_free(u); u = NULL; }
    if (uk != NULL) { accfft_free(uk); uk = NULL; }

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief initialize class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode RegOpt::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_SetupDone = false;
    this->m_FFT.plan = NULL;
    this->m_FFT.mpicomm = NULL;
    this->m_FFT.osize[0] = 0;
    this->m_FFT.osize[1] = 0;
    this->m_FFT.osize[2] = 0;
    this->m_FFT.ostart[0] = 0;
    this->m_FFT.ostart[1] = 0;
    this->m_FFT.ostart[2] = 0;

    this->m_Domain.nlocal = 0;
    this->m_Domain.nglobal = 0;
    this->m_Domain.isize[0] = 0;
    this->m_Domain.isize[1] = 0;
    this->m_Domain.isize[2] = 0;
    this->m_Domain.istart[0] = 0;
    this->m_Domain.istart[1] = 0;
    this->m_Domain.istart[2] = 0;
    this->m_Domain.nt = 4;      ///< number of time points
    this->m_Domain.nc = 2;      ///< number of components (vector valued images)
    this->m_Domain.nx[0] = 32;
    this->m_Domain.nx[1] = 32;
    this->m_Domain.nx[2] = 32;
    this->m_Domain.timehorizon[0] = 0.0;
    this->m_Domain.timehorizon[1] = 1.0;

    this->m_RegModel = COMPRESSIBLE;

    this->m_RegNorm.type = H2SN;
    this->m_RegNorm.beta[0] = 1E-2;
    this->m_RegNorm.beta[1] = 1E-2;
    this->m_RegNorm.beta[2] = 1E-4;

    this->m_Verbosity = 0;

    this->m_PDESolver.type = SL;
    this->m_PDESolver.cflnumber = 0;
    this->m_PDESolver.order = 2;

    // smoothing
    this->m_Sigma[0] = 1.0;
    this->m_Sigma[1] = 1.0;
    this->m_Sigma[2] = 1.0;

    this->m_KrylovSolverPara.tol[0] = 1E-12;     ///< relative tolerance
    this->m_KrylovSolverPara.tol[1] = 1E-12;     ///< absolute tolerance
    this->m_KrylovSolverPara.tol[2] = 1E+06;     ///< divergence tolerance
    this->m_KrylovSolverPara.maxit = 1000;       ///< max number of iterations
    this->m_KrylovSolverPara.reltol = 1E-12;     ///< relative tolerance (actually computed in solver)
    this->m_KrylovSolverPara.fseqtype = QDFS;
    this->m_KrylovSolverPara.pctype = INVREG;
    this->m_KrylovSolverPara.solver = PCG;
    this->m_KrylovSolverPara.pcsetupdone = false;
    this->m_KrylovSolverPara.g0norm = 0;
    this->m_KrylovSolverPara.g0normset = false;
    this->m_KrylovSolverPara.iter = 0;           ///< divergence tolerance
    this->m_KrylovSolverPara.pcsolver = PCG;
    this->m_KrylovSolverPara.pctolscale = 1E-1;
    this->m_KrylovSolverPara.pcmaxit = 10;
    this->m_KrylovSolverPara.pcgridscale = 2;
    this->m_KrylovSolverPara.pctol[0] = 1E-12;   ///< relative tolerance
    this->m_KrylovSolverPara.pctol[1] = 1E-12;   ///< absolute tolerance
    this->m_KrylovSolverPara.pctol[2] = 1E+06;   ///< divergence tolerance
    this->m_KrylovSolverPara.usepetsceigest = true;
    this->m_KrylovSolverPara.matvectype = DEFAULTMATVEC;
//    this->m_KrylovSolverPara.matvectype = PRECONDMATVEC;
//    this->m_KrylovSolverPara.matvectype = PRECONDMATVECSYM;
    this->m_KrylovSolverPara.reesteigvals = false;
    this->m_KrylovSolverPara.eigvalsestimated = false;
    this->m_KrylovSolverPara.checkhesssymmetry = false;

    // tolerances for optimization
    this->m_OptPara.tol[0] = 1E-6;                  ///< grad abs tol ||g(x)|| < tol
    this->m_OptPara.tol[1] = 1E-16;                 ///< grad rel tol ||g(x)||/J(x) < tol
    this->m_OptPara.tol[2] = 1E-2;                  ///< grad rel tol ||g(x)||/||g(x0)|| < tol
    this->m_OptPara.maxit = 1E3;                    ///< max number of iterations
    this->m_OptPara.minit = 0;                      ///< min number of iterations
    this->m_OptPara.method = GAUSSNEWTON;           ///< optmization method
    this->m_OptPara.fastsolve = false;              ///< switch on fast solver (less accurate)
    this->m_OptPara.fastpresolve = true;            ///< enable fast (inaccurate) solve for first steps
    this->m_OptPara.usezeroinitialguess = true;     ///< use zero initial guess for optimization

    // tolerances for presolve
    this->m_OptPara.presolvetol[0] = this->m_OptPara.tol[0];    ///< grad abs tol ||g(x)|| < tol
    this->m_OptPara.presolvetol[1] = this->m_OptPara.tol[1];    ///< grad rel tol ||g(x)||/J(x) < tol
    this->m_OptPara.presolvetol[2] = 1E-1;                      ///< grad rel tol ||g(x)||/||g(x0)|| < tol
    this->m_OptPara.presolvemaxit = this->m_OptPara.maxit;      ///< max number of iterations (multilevel; parameter continuation)

    this->m_SolveType = NOTSET;

    // flags
    this->m_ReadWriteFlags.templateim = false;
    this->m_ReadWriteFlags.referenceim = false;
    this->m_ReadWriteFlags.readfiles = false;       ///< read images to file
    this->m_ReadWriteFlags.timeseries = false;      ///< write time series to file (time dependent variables; use with caution) to file
    this->m_ReadWriteFlags.iterates = false;        ///< write iterates (velocity field; use with caution) to file
    this->m_ReadWriteFlags.results = false;         ///< write results (velocity field) to file
    this->m_ReadWriteFlags.defgrad = false;         ///< write deformation gradient (entire tensor) to file
    this->m_ReadWriteFlags.detdefgrad = false;      ///< write deformation gradient (determinant of deformation gradient) to file
    this->m_ReadWriteFlags.residual = false;        ///< write residual images to file
    this->m_ReadWriteFlags.defmap = false;          ///< write deformation map to file
    this->m_ReadWriteFlags.deffield = false;        ///< write deformation field / displacement field to file
    this->m_ReadWriteFlags.velnorm = false;         ///< write norm of velocity field to file
    this->m_ReadWriteFlags.deftemplate = false;     ///< write deformed template image to file
    this->m_ReadWriteFlags.extension = ".nii.gz";   ///< file extension for output

    this->m_RegFlags.applysmoothing = true;           ///< enable/disable image smoothing
    this->m_RegFlags.applyrescaling = true;           ///< enable/disable image rescaling
    this->m_RegFlags.detdefgradfromdeffield = false;    ///< compute det(grad(y)) via displacement field u
    this->m_RegFlags.invdefgrad = false;                ///< compute inverse of det(grad(y))^{-1}

    // parameter continuation
    this->m_ParaCont.strategy = PCONTOFF;       ///< no continuation
    this->m_ParaCont.enabled = false;           ///< flag for parameter continuation
    this->m_ParaCont.targetbeta = 0.0;          ///< has to be set by user
    //this->m_ParaCont.beta0 = 1.0;               ///< default initial parameter for parameter continuation
    this->m_ParaCont.beta0 = 1E-1;               ///< default initial parameter for parameter continuation

    // grid continuation
    this->m_GridCont.enabled = false;
    this->m_GridCont.nlevels = 0;

    // scale continuation
    this->m_ScaleCont.enabled = false;
    for (int i = 0; i < 3; ++i) {
        this->m_ScaleCont.sigma[i][0] = 32.0;
        this->m_ScaleCont.sigma[i][1] = 16.0;
        this->m_ScaleCont.sigma[i][2] =  8.0;
        this->m_ScaleCont.sigma[i][3] =  4.0;
        this->m_ScaleCont.sigma[i][4] =  2.0;
        this->m_ScaleCont.sigma[i][5] =  1.0;
    }

    // monitor for registration
    this->m_RegMonitor.JAC = false;
    this->m_RegMonitor.CFL = false;
    this->m_RegMonitor.jacmin = 0.0;
    this->m_RegMonitor.jacmax = 0.0;
    this->m_RegMonitor.jacmean = 0.0;
    this->m_RegMonitor.jacbound = 2E-1;
    this->m_RegMonitor.boundreached = 2E-1;

    for (int i = 0; i < NLOGFLAGS; ++i) {
        this->m_Log.enabled[i] = false;
    }

    this->m_NumThreads = 1;
    this->m_CartGridDims[0] = 1;
    this->m_CartGridDims[1] = 1;

    this->m_Indent = 0;
    this->m_LineLength = 101;
    this->m_StoreCheckPoints = false;

    ierr = this->ResetTimers(); CHKERRQ(ierr);
    ierr = this->ResetCounters(); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief display usage message for binary
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Usage"
PetscErrorCode RegOpt::Usage(bool advanced) {
    PetscErrorCode ierr = 0;
    int rank;
    std::string line;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    line = std::string(this->m_LineLength, '-');

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << line << std::endl;
        std::cout << " usage: runcoldreg [options] " << std::endl;
        std::cout << line << std::endl;
        std::cout << " where [options] is one or more of the following" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -mr <file>                reference image (*.nii, *.nii.gz, *.hdr)" << std::endl;
        std::cout << " -mt <file>                template image (*.nii, *.nii.gz, *.hdr)" << std::endl;

        // ####################### advanced options #######################
        if (advanced) {
        std::cout << " -vx1 <file>               x1 component of velocity field (*.nii, *.nii.gz, *.hdr, *.nc)" << std::endl;
        std::cout << " -vx2 <file>               x2 component of velocity field (*.nii, *.nii.gz, *.hdr, *.nc)" << std::endl;
        std::cout << " -vx3 <file>               x3 component of velocity field (*.nii, *.nii.gz, *.hdr, *.nc)" << std::endl;
        std::cout << " -sigma <int>x<int>x<int>  size of gaussian smoothing kernel applied to input images" << std::endl;
        std::cout << "                           (e.g., 1x2x1; units: voxel size; if only one parameter is set" << std::endl;
        std::cout << "                           uniform smoothing is assumed: default: 1x1x1)" << std::endl;
        std::cout << " -disablesmoothing         flag: switch off smoothing of image data" << std::endl;
        std::cout << " -disablerescaling         flag: switch off rescaling of intensities of image data to [0,1]" << std::endl;
        }
        // ####################### advanced options #######################

        std::cout << line << std::endl;
        std::cout << " -xresult                  output inversion variable (by default only velocity field will" << std::endl;
        std::cout << "                           be written to file; for more output options, see flags)" << std::endl;
        std::cout << " -x <path>                 output path (a prefix can be added by doing" << std::endl;
        std::cout << "                           '-x /output/path/prefix_')" << std::endl;

        // ####################### advanced options #######################
        if (advanced) {
        std::cout << " -xdefgrad                 write deformation gradient to file" << std::endl;
        std::cout << " -xdetdefgrad              write determinant of deformation gradient to file" << std::endl;
        std::cout << " -xdefmap                  write deformation map to file" << std::endl;
        std::cout << " -xdeffield                write deformation field/displacement field to file" << std::endl;
        std::cout << " -xdeftemplate             write deformed/transported template image to file" << std::endl;
        std::cout << " -xresidual                write pointwise residual (before and after registration) to file" << std::endl;
        std::cout << line << std::endl;
        std::cout << " optimization specific parameters" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -optmeth <type>           optimization method" << std::endl;
        std::cout << "                           <type> is one of the following" << std::endl;
        std::cout << "                               gn           Gauss-Newton (default)" << std::endl;
        std::cout << "                               fn           full Newton" << std::endl;
        std::cout << " -grel <dbl>               tolerance for optimization (default: 1E-2)" << std::endl;
        std::cout << "                               relative change of gradient" << std::endl;
        std::cout << "                               optimization stops if ||g_k||/||g_0|| <= tol" << std::endl;
        std::cout << " -gabs <dbl>               tolerance for optimization (default: 1E-6)" << std::endl;
        std::cout << "                               lower bound for gradient" << std::endl;
        std::cout << "                               optimization stops if ||g_k|| <= tol" << std::endl;
        std::cout << "                               tol <= ||g||/||g_init||" << std::endl;
        std::cout << " -maxit <int>              maximum number of (outer) Newton iterations (default: 50)" << std::endl;
        std::cout << " -krylovsolver <type>      solver for reduced space hessian system H[vtilde]=-g" << std::endl;
        std::cout << "                           <type> is one of the following" << std::endl;
        std::cout << "                               pcg          preconditioned conjugate gradient method" << std::endl;
        std::cout << "                               gmres        generalized minimal residual method" << std::endl;
        std::cout << " -krylovmaxit <int>        maximum number of (inner) Krylov iterations (default: 50)" << std::endl;
        std::cout << " -krylovfseq <type>        forcing sequence for Krylov solver (tolerance for inner iterations)" << std::endl;
        std::cout << "                           <type> is one of the following" << std::endl;
        std::cout << "                               quadratic     quadratic (default)" << std::endl;
        std::cout << "                               suplinear     super-linear" << std::endl;
        std::cout << "                               none          exact solve (expensive)" << std::endl;
        std::cout << " -krylovtol <dbl>          relative tolerance for krylov method (default: 1E-12); forcing sequence" << std::endl;
        std::cout << "                           needs to be switched off (i.g., use with '-krylovfseq none')" << std::endl;
        std::cout << " -precond <type>           preconditioner" << std::endl;
        std::cout << "                           <type> is one of the following" << std::endl;
        std::cout << "                               none         no preconditioner (not recommended)" << std::endl;
        std::cout << "                               invreg       inverse regularization operator (default)" << std::endl;
        std::cout << "                               2level       2-level preconditioner" << std::endl;
        std::cout << " -pcsolver <type>          solver for inversion of preconditioner (in case" << std::endl;
        std::cout << "                           the 2-level preconditioner is used)" << std::endl;
        std::cout << "                           <type> is one of the following" << std::endl;
        std::cout << "                               cheb         chebyshev method (default)" << std::endl;
        std::cout << "                               pcg          preconditioned conjugate gradient method" << std::endl;
        std::cout << "                               fpcg         flexible pcg" << std::endl;
        std::cout << "                               gmres        generalized minimal residual method" << std::endl;
        std::cout << "                               fgmres       flexible gmres" << std::endl;
        std::cout << " -pcsolvermaxit <int>      maximum number of iterations for inverting preconditioner; is" << std::endl;
        std::cout << "                           used for cheb, fgmres and fpcg; default: 10" << std::endl;
        std::cout << " -checksymmetry            check symmetry of hessian operator" << std::endl;
        std::cout << " -reesteigvals             re-estimate eigenvalues of hessian operator at every outer iteration" << std::endl;
        std::cout << "                           (in case a chebyshev method is used to invert preconditioner)" << std::endl;
        std::cout << " -pctolscale <dbl>         scale for tolerance (preconditioner needs to be inverted more" << std::endl;
        std::cout << "                           accurately then hessian; used for gmres and pcg; default: 1E-1)" << std::endl;
        std::cout << " -nonzeroinitialguess      use a non-zero velocity field to compute the initial gradient" << std::endl;
        std::cout << "                           this is only recommended in case one want to solve more accurately" << std::endl;
        std::cout << "                           after a warm start (in general for debugging purposes only)" << std::endl;
        std::cout << line << std::endl;
        std::cout << " regularization/constraints" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -regnorm <type>           regularization norm for velocity field" << std::endl;
        std::cout << "                           <type> is one of the following" << std::endl;
        std::cout << "                               h1s          H1-seminorm" << std::endl;
        std::cout << "                               h2s          H2-seminorm (default)" << std::endl;
        std::cout << "                               h3s          H3-seminorm" << std::endl;
        std::cout << "                               h1           H1-norm" << std::endl;
        std::cout << "                               h2           H2-norm" << std::endl;
        std::cout << "                               h3           H3-norm" << std::endl;
        std::cout << "                               l2           l2-norm (discouraged)" << std::endl;
        std::cout << " -betav <dbl>              regularization parameter (velocity field; default: 1E-2)" << std::endl;
        std::cout << " -betaw <dbl>              regularization parameter (mass source map; default: 1E-4; enable relaxed" << std::endl;
        std::cout << "                           incompressibility to use this parameter via '-ric' option; see below)" << std::endl;
        std::cout << " -ic                       enable incompressibility constraint (det(grad(y))=1)" << std::endl;
        std::cout << " -ric                      enable relaxed incompressibility (control jacobians; det(grad(y)) ~ 1)" << std::endl;
        std::cout << " -scalecont                enable scale continuation (continuation in smoothness of images;" << std::endl;
        std::cout << "                           i.e., use a multi-scale scheme to solve optimization problem)" << std::endl;
        std::cout << " -gridcont                 enable grid continuation (continuation in resolution of images;" << std::endl;
        std::cout << "                           i.e., use multi-resultion scheme to solve optimization probelm)" << std::endl;
        }
        // ####################### advanced options #######################

        std::cout << " -fastsolve                switch on fast solve (preset number of iterations and tolerances to" << std::endl;
        std::cout << "                           reduce the time to solution; inaccurate solve)" << std::endl;
        std::cout << " -train <type>             estimate regularization parameter (use 'jbound' to set bound" << std::endl;
        std::cout << "                           for det(grad(y)) used during estimation)" << std::endl;
        std::cout << "                           <type> is one of the following" << std::endl;
        std::cout << "                               binary       perform binary search (recommended)" << std::endl;
        std::cout << "                               reduce       reduce parameter by one order until bound is breached" << std::endl;
        std::cout << " -jbound <dbl>             lower bound on determinant of deformation gradient (default: 2E-1)" << std::endl;
        std::cout << " -betavcont <dbl>          do parameter continuation in betav until target regularization" << std::endl;
        std::cout << "                           parameter betav=<dbl> is reached (betav must be in (0,1))" << std::endl;
        std::cout << " -betavinit <dbl>          initial regularization weight for continuation" << std::endl;

        // ####################### advanced options #######################
        if (advanced) {
        std::cout << " -mdefgrad                 enable monitor for det(grad(y))" << std::endl;
        std::cout << line << std::endl;
        std::cout << " solver specific parameters (numerics)" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -pdesolver <type>         numerical time integrator for transport equations" << std::endl;
        std::cout << "                           <type> is one of the following" << std::endl;
        std::cout << "                               sl           semi-Lagrangian method (default; unconditionally stable)" << std::endl;
        std::cout << "                               rk2          rk2 time integrator (conditionally stable)" << std::endl;
        std::cout << " -nt <int>                 number of time points (for time integration; default: 4)" << std::endl;
        std::cout << line << std::endl;
        std::cout << " memory distribution and parallelism" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -nthreads <int>           number of threads (default: 1)" << std::endl;
        std::cout << " -np <int>x<int>           distribution of mpi tasks (cartesian grid) (example: -np 2x4 results" << std::endl;
        std::cout << "                           results in MPI distribution of size (nx1/2,nx2/4,nx3) for each mpi task)" << std::endl;
        std::cout << line << std::endl;
        std::cout << " logging" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -logresidual              log residual (user needs to set '-x' option)" << std::endl;
        std::cout << " -logkrylovres             log residual of krylov subpsace method (user needs to set '-x' option)" << std::endl;
        std::cout << " -logworkload              log cpu time and counters (user needs to set '-x' option)" << std::endl;
        std::cout << " -storecheckpoints         store iterates after each iteration (files will be overwritten); this is" << std::endl;
        std::cout << "                           a safeguard for large scale runs in case the code crashes" << std::endl;
        std::cout << line << std::endl;
        std::cout << " other parameters/debugging" << std::endl;
        std::cout << line << std::endl;
        std::cout << " -xiterates                store/write out iterates (deformed template image and velocity field)" << std::endl;
        std::cout << " -xiresults                store intermediate results/data (for scale, grid, and para continuation)" << std::endl;
        std::cout << " -xtimeseries              store time series (use with caution)" << std::endl;
        std::cout << " -nx <int>x<int>x<int>     grid size (e.g., 32x64x32); allows user to control grid size for synthetic" << std::endl;
        std::cout << "                           problems; assumed to be uniform if single integer is provided" << std::endl;
        std::cout << " -usenc                    use netcdf as output format (*.nc); default is NIFTI (*.nii.gz)" << std::endl;
//        std::cout << " -usebin                   use binary files as output format (*.bin)" << std::endl;
//        std::cout << " -usehdf5                  use hdf files as output format (*.hdf5)" << std::endl;
        std::cout << " -verbosity <int>          verbosity level (ranges from 0 to 2; default: 0)" << std::endl;
        }
        // ####################### advanced options #######################

        std::cout << line << std::endl;
        std::cout << " -help                     display a brief version of the user message" << std::endl;
        std::cout << " -advanced                 display this message" << std::endl;
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
#undef __FUNCT__
#define __FUNCT__ "CheckArguments"
PetscErrorCode RegOpt::CheckArguments() {
    PetscErrorCode ierr = 0;
    bool readmR = false, readmT = false, loggingenabled = false,
         readvx1 = false, readvx2 = false, readvx3 = false;
    ScalarType betav;

    std::string msg;
    PetscFunctionBegin;

    if (!this->m_ReadWriteFlags.mt.empty()) {readmT = true;}
    if (!this->m_ReadWriteFlags.mr.empty()) {readmR = true;}

    if (!this->m_ReadWriteFlags.vx1.empty()) {readvx1 = true;}
    if (!this->m_ReadWriteFlags.vx2.empty()) {readvx2 = true;}
    if (!this->m_ReadWriteFlags.vx3.empty()) {readvx3 = true;}

    if (readmT && readmR) {
        // check if files exist
        msg = "file " + this->m_ReadWriteFlags.mt + " does not exist";
        ierr = Assert(FileExists(this->m_ReadWriteFlags.mt), msg); CHKERRQ(ierr);
        msg = "file " + this->m_ReadWriteFlags.mr + " does not exist";
        ierr = Assert(FileExists(this->m_ReadWriteFlags.mr), msg); CHKERRQ(ierr);
        this->m_ReadWriteFlags.readfiles = true;
    } else if ( (readmT == false) && readmR ) {
        msg = "\x1b[31m you need to also assign a template image\x1b[0m\n";
        ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
        ierr = this->Usage(); CHKERRQ(ierr);
    } else if ( readmT && (readmR == false) ) {
        msg = "\x1b[31m you need to also assign a reference image\x1b[0m\n";
        ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
        ierr = this->Usage(); CHKERRQ(ierr);
    } else if ( (readmT == false) && (readmR == false) ) {
        this->m_ReadWriteFlags.readfiles = false;
    }

    if (readvx1 && readvx2 && readvx3) {
        // check if files exist
        msg = "file " + this->m_ReadWriteFlags.vx1 + " does not exist";
        ierr = Assert(FileExists(this->m_ReadWriteFlags.vx1), msg); CHKERRQ(ierr);
        msg = "file " + this->m_ReadWriteFlags.vx2 + " does not exist";
        ierr = Assert(FileExists(this->m_ReadWriteFlags.vx2), msg); CHKERRQ(ierr);
        msg = "file " + this->m_ReadWriteFlags.vx3 + " does not exist";
        ierr = Assert(FileExists(this->m_ReadWriteFlags.vx3), msg); CHKERRQ(ierr);
        this->m_ReadWriteFlags.readvelocity = true;
    }

    if (this->m_ParaCont.strategy == PCONTINUATION) {
        betav = this->m_ParaCont.targetbeta;
        if (betav <= 0.0 || betav > 1.0) {
            msg = "\x1b[31m target betav not in (0.0,1.0]\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(); CHKERRQ(ierr);
        }
    }
    if (   this->m_ParaCont.strategy == PCONTBINSEARCH
        || this->m_ParaCont.strategy == PCONTINUATION  ) {
        betav = this->m_ParaCont.beta0;
        if (betav <= 0.0 || betav > 1.0) {
            msg = "\x1b[31m initial guess for betav not in (0.0,1.0]\x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(); CHKERRQ(ierr);
        }
    }

    if (this->m_ScaleCont.enabled && this->m_ParaCont.enabled) {
        msg = "\x1b[31m combined parameter and scale continuation not available \x1b[0m\n";
        ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
        ierr = this->Usage(); CHKERRQ(ierr);
    }

    // check output arguments
    if (   this->m_ReadWriteFlags.results
        || this->m_ReadWriteFlags.defgrad
        || this->m_ReadWriteFlags.detdefgrad
        || this->m_ReadWriteFlags.defmap
        || this->m_ReadWriteFlags.residual
        || this->m_ReadWriteFlags.timeseries
        || this->m_ReadWriteFlags.iterates ) {
        if (this->m_ReadWriteFlags.xfolder.empty()) {
            msg = "\x1b[31m output folder needs to be set (-x option) \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(); CHKERRQ(ierr);
        }
    }

    for (int i = 0; i < NLOGFLAGS; ++i) {
        if (this->m_Log.enabled[i]) {
            loggingenabled = true;
        }
    }

    // check output arguments
    if (loggingenabled) {
        if (this->m_ReadWriteFlags.xfolder.empty()) {
            msg = "\x1b[31m output folder needs to be set (-x option) \x1b[0m\n";
            ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
            ierr = this->Usage(); CHKERRQ(ierr);
        }
    }

    if (   this->m_KrylovSolverPara.pctolscale < 0.0
        || this->m_KrylovSolverPara.pctolscale >= 1.0 ) {
        msg = "\x1b[31m tolerance for precond solver out of bounds; not in (0,1)\x1b[0m\n";
        ierr = PetscPrintf(PETSC_COMM_WORLD, msg.c_str()); CHKERRQ(ierr);
        ierr = this->Usage(); CHKERRQ(ierr);
    }

    ierr = Assert(this->m_NumThreads > 0, "omp threads < 0"); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief setup options and accfft
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DoSetup"
PetscErrorCode RegOpt::DoSetup(bool dispteaser) {
    PetscErrorCode ierr = 0;
    std::stringstream ss;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    ierr = this->InitializeFFT(); CHKERRQ(ierr);

    // display the options to the user
    if (dispteaser) {
        ierr = this->DisplayOptions(); CHKERRQ(ierr);
    }

    this->m_SetupDone = true;

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief setup of accfft has been done externally; this function
 * uses the fft class to and communcation layout as an input
 * and sets the associated parameters
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CouplingSetup"
PetscErrorCode RegOpt::CouplingSetup(IntType nx[3]) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    // parse grid size for setup
    for (int i = 0; i < 3; ++i) {
        this->m_Domain.nx[i] = nx[i];
    }

    this->m_Verbosity = 2;

    // compute solution faster

    // optimization parameters
    this->m_OptPara.tol[2] = 5E-2;
    this->m_OptPara.maxit = 10;

    // parameters for solver of reduced space KKT system
    this->m_KrylovSolverPara.fseqtype = QDFS;
    this->m_KrylovSolverPara.pctype = INVREG;
    this->m_KrylovSolverPara.maxit = 10;

    ierr = this->DoSetup(true); CHKERRQ(ierr);

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set preset parameters (maybe reduce number of krylov
 * iterations)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EnableFastSolve"
PetscErrorCode RegOpt::EnableFastSolve() {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    this->m_SolveType = FAST_SMOOTH;
    ierr = this->SetPresetParameters(); CHKERRQ(ierr);

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set preset parameters / provide a crude estimate for users
 * either reduce the time to solution or compute high-fidelity
 * results
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetPresetParameters"
PetscErrorCode RegOpt::SetPresetParameters() {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    if (this->m_SolveType == FAST_AGG) {
        // use fast and aggressive method
        this->m_RegNorm.type = H1SN;
        this->m_KrylovSolverPara.fseqtype = QDFS;
        this->m_KrylovSolverPara.maxit = 5;
        this->m_OptPara.maxit = 20;
        this->m_OptPara.presolvemaxit = 10;
        this->m_OptPara.tol[2] = 5E-2;  // use 5E-2 if results are not acceptable reduce; 1E-1 and 2.5E-1
        this->m_OptPara.presolvetol[2] = 5E-1;
    } else if (this->m_SolveType == ACC_AGG) {
        // use slow and aggressive method
        this->m_RegNorm.type = H1SN;
        this->m_KrylovSolverPara.fseqtype = QDFS;
        this->m_KrylovSolverPara.maxit = 50;
        this->m_OptPara.maxit = 50;
        this->m_OptPara.presolvemaxit = 20;
        this->m_OptPara.tol[2] = 1E-2;
        this->m_OptPara.presolvetol[2] = 1E-1;
    } else if (this->m_SolveType == FAST_SMOOTH) {
        // use fast and smooth method
        //this->m_RegNorm.type = H2SN;
        this->m_RegNorm.type = H2;
        this->m_KrylovSolverPara.fseqtype = QDFS;
        this->m_KrylovSolverPara.maxit = 5;
        this->m_OptPara.maxit = 20;
        this->m_OptPara.presolvemaxit = 10;
        this->m_OptPara.tol[2] = 5E-2;
        this->m_OptPara.presolvetol[2] = 5E-1;
    } else if (this->m_SolveType == ACC_SMOOTH) {
        // use slow and smooth method
        //this->m_RegNorm.type = H2SN;
        this->m_RegNorm.type = H2;
        this->m_KrylovSolverPara.fseqtype = QDFS;
        this->m_KrylovSolverPara.maxit = 50;
        this->m_OptPara.maxit = 50;
        this->m_OptPara.presolvemaxit = 20;
        this->m_OptPara.tol[2] = 1E-2;
        this->m_OptPara.presolvetol[2] = 1E-1;
    } else {
        ierr = ThrowError("flag not defined"); CHKERRQ(ierr);
    }

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief display registration options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetBetaMinParaCont"
ScalarType RegOpt::GetBetaMinParaCont() {
    if (this->m_RegNorm.type == H1) {
        return this->m_ParaCont.betavminh1;
    } else if (this->m_RegNorm.type == H2) {
        return this->m_ParaCont.betavminh2;
    } else if (this->m_RegNorm.type == H2SN) {
        return this->m_ParaCont.betavminh2;
    } else if (this->m_RegNorm.type == H1SN) {
        return this->m_ParaCont.betavminh1;
    } else { return 1E-9; }
}




/********************************************************************
 * @brief set up grid continuation
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupGridCont"
PetscErrorCode RegOpt::SetupGridCont() {
    PetscErrorCode ierr = 0;
    IntType nxmin, nxi, nl, ng, nalloc;
    int nx[3], isize[3], istart[3], ostart[3], osize[3];
    int nlevels, level, j;
    ScalarType value;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    // compute number of levels
    nxmin = this->m_Domain.nx[0];
    for (int i = 1; i < 3; ++i) {
        nxi = this->m_Domain.nx[i];
        nxmin = nxmin < nxi ? nxmin : nxi;
    }

    nlevels  = static_cast<int>(std::ceil(PetscLog2Real(static_cast<ScalarType>(nxmin))));
    nlevels -= static_cast<int>(this->m_GridCont.minlevels);
    ierr = Assert(nlevels > 0, "error in size"); CHKERRQ(ierr);
    this->m_GridCont.nlevels = nlevels;

    // allocate arrays for sizes
    this->m_GridCont.nx.resize(nlevels);        ///< grid size per level
    this->m_GridCont.isize.resize(nlevels);     ///< grid size per level (spatial domain)
    this->m_GridCont.istart.resize(nlevels);    ///< start index per level (spatial domain)
    this->m_GridCont.osize.resize(nlevels);     ///< grid size per level (spectral domain)
    this->m_GridCont.ostart.resize(nlevels);    ///< start index per level (spectral domain)

    for (int i = 0; i < nlevels; ++i) {
        this->m_GridCont.nx[i].resize(3);
        this->m_GridCont.istart[i].resize(3);
        this->m_GridCont.isize[i].resize(3);
        this->m_GridCont.ostart[i].resize(3);
        this->m_GridCont.osize[i].resize(3);
    }

    this->m_GridCont.nlocal.resize(nlevels);    ///< local points (MPI task) per level
    this->m_GridCont.nglobal.resize(nlevels);   ///< global points per level
    this->m_GridCont.nalloc.resize(nlevels);    ///< alloc size in fourier domain

    level = 0;
    while (level < nlevels) {
        j = nlevels-(level+1);
        nl = 1;  // reset local size
        ng = 1;  // reset global size

        // compute number of grid points for current level
        for (int i = 0; i < 3; ++i) {
            if (level == 0) {
                this->m_GridCont.nx[j][i] = this->m_Domain.nx[i];
            } else {
                value = static_cast<ScalarType>(this->m_GridCont.nx[j+1][i]);
                this->m_GridCont.nx[j][i] = static_cast<IntType>(std::ceil(value/2.0));
            }

            // compute global size
            ng *= this->m_GridCont.nx[j][i];
            nx[i] = static_cast<int>(this->m_GridCont.nx[j][i]);
        }

        this->m_GridCont.nglobal[j] = ng;  // set global size

        // get the local sizes
        nalloc = accfft_local_size_dft_r2c(nx, isize, istart, osize, ostart, this->m_FFT.mpicomm);
        this->m_GridCont.nalloc[j] = nalloc;

        // compute local sizes
        for (int i = 0; i < 3; ++i) {
            // compute number of local points
            nl *= static_cast<IntType>(isize[i]);

            // sizes in spatial domain
            this->m_GridCont.isize[j][i] = static_cast<IntType>(isize[i]);
            this->m_GridCont.istart[j][i] = static_cast<IntType>(istart[i]);

            // sizes in spectral domain
            this->m_GridCont.osize[j][i] = static_cast<IntType>(osize[i]);
            this->m_GridCont.ostart[j][i] = static_cast<IntType>(ostart[i]);
        }

        this->m_GridCont.nlocal[j] = nl;  // set local size

        ++level;  // increment
    }

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief display registration options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DisplayOptions"
PetscErrorCode RegOpt::DisplayOptions() {
    PetscErrorCode ierr;
    int rank, indent, align;
    std::string msg, line;
    bool newtontype;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    indent = 40;
    align  = 30;
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

        if (this->m_Domain.nc == 1) {
            std::cout << std::left << std::setw(indent) << " problem dimensions"
                        << "(nx1,nx2,nx3,nt)=(" << this->m_Domain.nx[0] << ","
                        <<  this->m_Domain.nx[1] << ","
                        <<  this->m_Domain.nx[2] << ","
                        <<  this->m_Domain.nt << ")" << std::endl;
        } else {
            std::cout << std::left << std::setw(indent) << " problem dimensions"
                        << "(nx1,nx2,nx3,nc,nt)=(" << this->m_Domain.nx[0] << ","
                        <<  this->m_Domain.nx[1] << ","
                        <<  this->m_Domain.nx[2] << ","
                        <<  this->m_Domain.nc << ","
                        <<  this->m_Domain.nt << ")" << std::endl;
        }
        std::cout << std::left << std::setw(indent) <<" network dimensions"
                    << this->m_CartGridDims[0] << "x"
                    << this->m_CartGridDims[1] << std::endl;
        std::cout << std::left << std::setw(indent) << " threads"
                    << this->m_NumThreads<< std::endl;
        std::cout << std::left << std::setw(indent) << " (ng,nl)"
                    << "(" << this->m_Domain.nglobal << ","
                    <<  this->m_Domain.nlocal << ")" << std::endl;

        std::cout << line << std::endl;
        std::cout << " parameters" << std::endl;
        std::cout << line << std::endl;

        // display regularization model
        std::cout << std::left << std::setw(indent) <<" regularization model v";

        if ((this->m_ParaCont.strategy == PCONTBINSEARCH) || (this->m_ParaCont.strategy == PCONTREDUCESEARCH)) {
            switch (this->m_RegNorm.type) {
                case L2:
                {
                    std::cout << "l2-norm (betav estimated)" << std::endl;
                    break;
                }
                case H1:
                {
                    std::cout << "h1-norm (betav estimated)" << std::endl;
                    break;
                }
                case H2:
                {
                    std::cout << "h2-norm (betav estimated)" << std::endl;
                    break;
                }
                case H3:
                {
                    std::cout << "h3-norm (betav estimated)" << std::endl;
                    break;
                }
                case H1SN:
                {
                    std::cout << "h1-seminorm (betav estimated)" << std::endl;
                    break;
                }
                case H2SN:
                {
                    std::cout << "h2-seminorm (betav estimated)" << std::endl;
                    break;
                }
                case H3SN:
                {
                    std::cout << "h3-seminorm (betav estimated)" << std::endl;
                    break;
                }
                default:
                {
                    ierr = ThrowError("regularization model not implemented"); CHKERRQ(ierr);
                    break;
                }
            }

            // display parameters and tolerances
            std::cout << std::left << std::setw(indent) << " parameter continuation";
            if (this->m_ParaCont.strategy == PCONTBINSEARCH) {
                std::cout << "binary search" << std::endl;
            } else if (this->m_ParaCont.strategy == PCONTREDUCESEARCH) {
                std::cout << "search by reduction" << std::endl;
            }
            std::cout << std::left << std::setw(indent) << " "
                      << std::setw(align) << "bound det(grad(y))"
                      << this->m_RegMonitor.jacbound << std::endl;
        } else {
            switch (this->m_RegNorm.type) {
                case L2:
                {
                    std::cout   << std::scientific << "l2-norm (betav="
                                << this->m_RegNorm.beta[0]
                                << ")" << std::endl;
                    break;
                }
                case H1:
                {
                    std::cout   << std::scientific << "h1-norm (betav="
                                << this->m_RegNorm.beta[0]
                                << ", " << this->m_RegNorm.beta[1]
                                << ") " << std::endl;
                    break;
                }
                case H2:
                {
                    std::cout   << std::scientific << "h2-norm (betav="
                                << this->m_RegNorm.beta[0]
                                << ", " << this->m_RegNorm.beta[1]
                                << ")" << std::endl;
                    break;
                }
                case H3:
                {
                    std::cout   << std::scientific << "h3-norm (betav="
                                << this->m_RegNorm.beta[0]
                                << ", " << this->m_RegNorm.beta[1]
                                << ")" << std::endl;
                    break;
                }
                case H1SN:
                {
                    std::cout   << std::scientific <<  "h1-seminorm (betav="
                                <<  this->m_RegNorm.beta[0]
                                << ")" << std::endl;
                    break;
                }
                case H2SN:
                {
                    std::cout   << std::scientific << "h2-seminorm (betav="
                                << this->m_RegNorm.beta[0]
                                << ")" << std::endl;
                    break;
                }
                case H3SN:
                {
                    std::cout   << std::scientific << "h3-seminorm (betav="
                                << this->m_RegNorm.beta[0]
                                << ")" << std::endl;
                    break;
                }
                default: {ierr = ThrowError("regularization model not implemented"); CHKERRQ(ierr); break;}
            }

            // display parameters and tolerances
            std::cout << std::left << std::setw(indent) <<" parameter continuation";
            if (this->m_ParaCont.strategy == PCONTINUATION) {
                std::cout << "enabled" << std::endl;
            } else if (this->m_ParaCont.strategy == PCONTOFF) {
                std::cout << "disabled" << std::endl;
            }
        }

        if (this->m_RegModel == reg::RELAXEDSTOKES) {
            // display regularization model
            std::cout << std::left << std::setw(indent) << " regularization model w";
            std::cout << "h1-seminorm (betaw="
                      << this->m_RegNorm.beta[2]<< ")" << std::endl;
        }

        // display regularization model
        std::cout << std::left << std::setw(indent) << " pde solver (hyperbolic)";
        switch (this->m_PDESolver.type) {
            case RK2:
            {
                std::cout << "second order rk method" << std::endl;
                break;
            }
            case RK2A:
            {
                std::cout << "antisymmetric rk2" << std::endl;
                break;
            }
            case SL:
            {
                std::cout << "semi-lagrangian method" << std::endl;
                break;
            }
            default:
            {
                ierr = ThrowError("solver not implemented"); CHKERRQ(ierr);
                break;
            }
        }

        // display type of optimization method
        newtontype = false;
        std::cout << std::left << std::setw(indent) << " optimization method";
        switch (this->m_OptPara.method) {
            case GRADDESCENT:
            {
                std::cout << "preconditioned gradient descent method" << std::endl;
                break;
            }
            case FULLNEWTON:
            {
                std::cout << "newton method" << std::endl;
                newtontype = true;
                break;
            }
            case GAUSSNEWTON:
            {
                std::cout << "gauss-newton method" << std::endl;
                newtontype = true;
                break;
            }
            default:
            {
                ierr = ThrowError("optimization method not implemented"); CHKERRQ(ierr);
                break;
            }
        }

        std::cout << std::left << std::setw(indent) << " maximal # iterations"
                  << this->m_OptPara.maxit << std::endl;

        // display optimization tolerances
        std::cout << std::left << std::setw(indent) << " convergence tolerances"
                  << std::setw(align) << "||g(v)|| <= tol"
                  << this->m_OptPara.tol[0] << std::endl;
//        std::cout << std::left << std::setw(indent) <<" "
//                  << std::setw(align) <<"||g(v)||/|J(v)| <= tol"
//                  << this->m_OptPara.tol[1] << std::endl;
        std::cout << std::left << std::setw(indent) << " "
                  << std::setw(align) << "||g(v)||/||g(v0)|| <= tol"
                  << this->m_OptPara.tol[2] << std::endl;

        // display parameters for newton type optimization methods
        if ( newtontype ) {
            std::cout << std::left << std::setw(indent)
                      << " hessian system"
                      << std::setw(align) << "solver";

            switch (this->m_KrylovSolverPara.solver) {
                case PCG:
                    this->m_KrylovSolverPara.name = "PCG";
                    break;
                case GMRES:
                    this->m_KrylovSolverPara.name = "GMRES";
                    break;
                default:
                    ierr = ThrowError("solver not defined"); CHKERRQ(ierr);
                    break;
            }

            std::cout << this->m_KrylovSolverPara.name << std::endl;

            std::cout << std::left << std::setw(indent) <<" "
                      << std::setw(align) << "forcing sequence";

            switch (this->m_KrylovSolverPara.fseqtype) {
                case NOFS:
                {
                    std::cout << std::setw(align) << "disabled" << std::endl;
                    std::cout << std::left << std::setw(indent) << " "
                              << std::setw(align) << "absolute"
                              << this->m_KrylovSolverPara.tol[0] << std::endl;
                    std::cout << std::left << std::setw(indent) << " "
                              << std::setw(align) << "relative"
                              << this->m_KrylovSolverPara.tol[1] << std::endl;
                    break;
                }
                case QDFS:
                {
                    std::cout << "quadratic" << std::endl;
                    break;
                }
                case SLFS:
                {
                    std::cout << "superlinear" << std::endl;
                    break;
                }
                default:
                {
                    ierr = ThrowError("forcing sequence not implemented"); CHKERRQ(ierr);
                    break;
                }
            }

            std::cout << std::left << std::setw(indent) << " "
                      << std::setw(align) << "maxit"
                      << this->m_KrylovSolverPara.maxit << std::endl;

            bool twolevel = false;
            std::cout << std::left << std::setw(indent) << " "
                      << std::setw(align) << "preconditioner";

            switch (this->m_KrylovSolverPara.pctype) {
                case INVREG:
                {
                    std::cout << "regularization operator" << std::endl;
                    break;
                }
                case TWOLEVEL:
                {
                    std::cout << "2-level multigrid" << std::endl;
                    twolevel = true;
                    break;
                }
                case NOPC:
                {
                    std::cout << "none" << std::endl;
                    break;
                }
                default:
                {
                    ierr = ThrowError("preconditioner not defined"); CHKERRQ(ierr);
                    break;
                }
            }

            if (twolevel) {
                std::cout << std::left << std::setw(indent) << " "
                          << std::setw(align) << "solver (preconditioner)";

                switch (this->m_KrylovSolverPara.pcsolver) {
                    case PCG:
                    {
                        this->m_KrylovSolverPara.pcname = "PCG";
                        break;
                    }
                    case FCG:
                    {
                        this->m_KrylovSolverPara.pcname = "FPCG";
                        break;
                    }
                    case GMRES:
                    {
                        this->m_KrylovSolverPara.pcname = "GMRES";
                        break;
                    }
                    case FGMRES:
                    {
                        this->m_KrylovSolverPara.pcname = "FGMRES";
                        break;
                    }
                    case CHEB:
                    {
                        this->m_KrylovSolverPara.pcname = "CHEB";
                        break;
                    }
                    default:
                    {
                        ierr = ThrowError("solver not defined"); CHKERRQ(ierr);
                        break;
                    }
                }
                std::cout << this->m_KrylovSolverPara.pcname << std::endl;
            }
/*          std::cout << std::left << std::setw(indent) <<" "
                      << std::setw(align) <<"divergence"
                      << this->m_KrylovSolverPara.tol[2] << std::endl;*/
        }
        std::cout << line << std::endl;
    }

    this->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief compute sizes
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetSizes"
PetscErrorCode RegOpt::GetSizes(IntType* nx, IntType& nl, IntType& ng) {
    PetscErrorCode ierr = 0;
    int nxi[3], isize[3], istart[3], osize[3], ostart[3];

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    ierr = Assert(nx[0] > 0, "error in size"); CHKERRQ(ierr);
    ierr = Assert(nx[1] > 0, "error in size"); CHKERRQ(ierr);
    ierr = Assert(nx[2] > 0, "error in size"); CHKERRQ(ierr);

    nxi[0] = static_cast<int>(nx[0]);
    nxi[1] = static_cast<int>(nx[1]);
    nxi[2] = static_cast<int>(nx[2]);

    accfft_local_size_dft_r2c(nxi, isize, istart, osize, ostart, this->m_FFT.mpicomm);

    nl = 1; ng = 1;
    for (int i = 0; i < 3; ++i) {
        ng *= static_cast<IntType>(nx[i]);
        nl *= static_cast<IntType>(isize[i]);
    }

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute sizes
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetSizes"
PetscErrorCode RegOpt::GetSizes(IntType* nx, IntType* istart, IntType* isize) {
    PetscErrorCode ierr = 0;
    int nxi[3], isizei[3], istarti[3], osize[3], ostart[3];

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    ierr = Assert(nx[0] > 0, "error in size"); CHKERRQ(ierr);
    ierr = Assert(nx[1] > 0, "error in size"); CHKERRQ(ierr);
    ierr = Assert(nx[2] > 0, "error in size"); CHKERRQ(ierr);

    nxi[0] = static_cast<int>(nx[0]);
    nxi[1] = static_cast<int>(nx[1]);
    nxi[2] = static_cast<int>(nx[2]);

    accfft_local_size_dft_r2c(nxi, isizei, istarti, osize, ostart, this->m_FFT.mpicomm);

    for (int i = 0; i < 3; ++i) {
        isize[i] = static_cast<IntType>(isizei[i]);
        istart[i] = static_cast<IntType>(istarti[i]);
    }

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute weight for FFT
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeFFTScale"
ScalarType RegOpt::ComputeFFTScale() {
    ScalarType scale = 1.0;
    for (int i = 0; i < 3; ++i) {
        scale *= static_cast<ScalarType>(this->m_Domain.nx[i]);
    }
    return 1.0/scale;
}




/********************************************************************
 * @brief resets all timers
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResetTimers"
PetscErrorCode RegOpt::ResetTimers() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    for (int i = 0; i < NTIMERS; ++i) {
        if (i != FFTSETUP) {
            for (int j = 0; j < NVALTYPES; ++j) {
                this->m_Timer[i][j] = 0.0;
            }
        }
        this->m_TimerIsRunning[i] = false;
        this->m_TempTimer[i] = 0.0;
    }

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < NVALTYPES; ++j) {
            this->m_FFTTimers[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < NVALTYPES; ++j) {
            this->m_InterpTimers[i][j] = 0.0;
        }
    }

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resets timer
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResetTimer"
PetscErrorCode RegOpt::ResetTimer(TimerType id) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    for (int j = 0; j < NVALTYPES; ++j) {
        this->m_Timer[id][j] = 0.0;
    }
    this->m_TimerIsRunning[id] = false;
    this->m_TempTimer[id] = 0.0;

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief start the timer (checks if running)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "StartTimer"
PetscErrorCode RegOpt::StartTimer(TimerType id) {
    PetscErrorCode ierr = 0;
    std::string msg;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    msg = "fatal error: timer has already been started";
    ierr = Assert(this->m_TimerIsRunning[id] == false, msg); CHKERRQ(ierr);

    this->m_TempTimer[id] = -MPI_Wtime();
    this->m_TimerIsRunning[id] = true;

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief stop setup timer (checks if running)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "StopTimer"
PetscErrorCode RegOpt::StopTimer(TimerType id) {
    PetscErrorCode ierr = 0;
    std::string msg;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    msg = "fatal error: timer has not been started";
    ierr = Assert(this->m_TimerIsRunning[id], msg); CHKERRQ(ierr);

    this->m_TempTimer[id] += MPI_Wtime();
    this->m_Timer[id][LOG] += this->m_TempTimer[id];

    // tell the world that we stop the timer
    this->m_TimerIsRunning[id] = false;

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief process the timers
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ProcessTimers"
PetscErrorCode RegOpt::ProcessTimers() {
    PetscErrorCode ierr = 0;
    int rval, rank, nproc;
    double ival = 0.0, xval = 0.0;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

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

    for (int i = 0; i < 5; ++i) {
        // remember input value
        ival = this->m_FFTTimers[i][LOG];

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

    for (int i = 0; i < 4; ++i) {
        // remember input value
        ival = this->m_InterpTimers[i][LOG];

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

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resets counters
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResetCounters"
PetscErrorCode RegOpt::ResetCounters() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    for (int i = 0; i < NCOUNTERS; ++i) {this->m_Counter[i] = 0;}

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resets counter
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResetCounter"
PetscErrorCode RegOpt::ResetCounter(CounterType id) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    this->m_Counter[id] = 0;

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write log results to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteLogFile"
PetscErrorCode RegOpt::WriteLogFile() {
    PetscErrorCode ierr = 0;

    if (this->m_Log.enabled[LOGLOAD]) {
        ierr = this->WriteWorkLoadLog(); CHKERRQ(ierr);
    }

    if (this->m_Log.enabled[LOGKSPRES]) {
        ierr = this->WriteKSPLog(); CHKERRQ(ierr);
    }

    if (this->m_Log.enabled[LOGRES]) {
        ierr = this->WriteResidualLog(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief displays the global exection time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteWorkLoadLog"
PetscErrorCode RegOpt::WriteWorkLoadLog() {
    PetscErrorCode ierr = 0;
    std::string fn, line, path;
    std::ofstream logwriter;
    std::stringstream ss, ssnum;
    int rank, nnum, nstr, nproc;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

    line = std::string(this->m_LineLength, '-');
    path = this->m_ReadWriteFlags.xfolder;

    // write out logfile
    if (rank == 0) {
        nnum = 20; nstr = 20;
        fn = path + "registration-performance.log";

        // create output file
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

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
                  << std::setw(nnum) << this->GetDomainPara().nglobal << std::endl;

        logwriter << std::left
                  << std::setw(nstr) << " nl" << std::right
                  << std::setw(nnum) << this->GetDomainPara().nlocal << std::endl;

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
                  << std::setw(nnum) << this->m_NumThreads << std::endl;

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
        if (this->m_FFTTimers[2][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[FFT] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " FFT communication" << std::right
                << std::setw(nnum) << this->m_FFTTimers[2][MIN]
                << std::setw(nnum) << this->m_FFTTimers[2][MAX]
                << std::setw(nnum) << this->m_FFTTimers[2][AVG]
                << std::setw(nnum) << this->m_FFTTimers[2][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[FFT]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_FFTTimers[4][LOG] > 0.0) {
            ierr = Assert(this->m_Counter[FFT] > 0, "bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " FFT execution" << std::right
                << std::setw(nnum) << this->m_FFTTimers[4][MIN]
                << std::setw(nnum) << this->m_FFTTimers[4][MAX]
                << std::setw(nnum) << this->m_FFTTimers[4][AVG]
                << std::setw(nnum) << this->m_FFTTimers[4][MAX]
                                    /static_cast<ScalarType>(this->m_Counter[FFT]);
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

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write residual to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteResidualLog"
PetscErrorCode RegOpt::WriteResidualLog() {
    PetscErrorCode ierr = 0;
    int rank, nnum;
    std::ofstream logwriter;
    std::stringstream ss;
    std::string path, fn;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    nnum = 20;
    path = this->m_ReadWriteFlags.xfolder;

    if (rank == 0) {
        // create output file
        fn = path + "cold-residual.log";
        logwriter.open(fn.c_str());
        ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-mT||_2" << std::right
            << std::setw(nnum) << this->m_Log.residual[0];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-mT||_infty" << std::right
            << std::setw(nnum) << this->m_Log.residual[1];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-m1||_2" << std::right
            << std::setw(nnum) << this->m_Log.residual[2];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-m1||_infty" << std::right
            << std::setw(nnum) << this->m_Log.residual[3];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        this->m_Log.residual[0] = (this->m_Log.residual[0] > 0.0) ? this->m_Log.residual[0] : 1.0;
        this->m_Log.residual[1] = (this->m_Log.residual[1] > 0.0) ? this->m_Log.residual[1] : 1.0;

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-m1||_2,rel" << std::right
            << std::setw(nnum) <<  this->m_Log.residual[2]/this->m_Log.residual[0];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        ss  << std::scientific << std::left
            << std::setw(20) << "||mR-m1||_infty,rel" << std::right
            << std::setw(nnum) <<  this->m_Log.residual[3]/this->m_Log.residual[1];
        logwriter << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        // close logger
        logwriter.close();
    }

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write out logging information for krylov method
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteKSPLog"
PetscErrorCode RegOpt::WriteKSPLog() {
    PetscErrorCode ierr = 0;
    int rank, n;
    std::ofstream logwriter;
    std::stringstream ss;
    std::string path, fn;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    path = this->m_ReadWriteFlags.xfolder;

    if (rank == 0) {
        // create output file
        fn = path + "cold-krylov-method-residual.log";
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

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief displays the global exection time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DisplayTimeToSolution"
PetscErrorCode RegOpt::DisplayTimeToSolution() {
    PetscErrorCode ierr = 0;
    double hours, minutes, seconds, millisec, time;
    int rank;
    std::stringstream ss;
    std::string line;
    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    time = this->m_Timer[T2SEXEC][MAX];

    hours    = time / 3600.0;
    minutes  = (hours   - floor(hours))   *   60.0;
    seconds  = (minutes - floor(minutes)) *   60.0;
    millisec = (seconds - floor(seconds)) * 1000.0;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    line = std::string(this->m_LineLength, '-');

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%s\n", line.c_str()); CHKERRQ(ierr);

    ss  << "computation finished ( elapsed cpu time "
        << floor(hours) << " h "
        << floor(minutes) << " m "
        << floor(seconds) << " s "
        << floor(millisec) << " ms )";

    ierr = Msg(ss.str()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%s\n", line.c_str()); CHKERRQ(ierr);

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _REGOPT_CPP_
