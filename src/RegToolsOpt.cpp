#ifndef _REGTOOLSOPT_CPP_
#define _REGTOOLSOPT_CPP_

#include "RegToolsOpt.hpp"




namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegOpt"
RegToolsOpt::RegToolsOpt()
{
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegOpt"
RegToolsOpt::RegToolsOpt(int argc, char** argv)
{
    this->Initialize();
    this->ParseArguments(argc,argv);
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegToolsOpt"
RegToolsOpt::RegToolsOpt(const RegToolsOpt& opt)
{
    this->Initialize();
    this->Copy(opt);
}




/********************************************************************
 * @brief parse user arguments
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ParseArguments"
PetscErrorCode RegToolsOpt::ParseArguments(int argc, char** argv)
{
    PetscErrorCode ierr=0;
    std::string msg;
    std::vector<unsigned int> np;
    std::vector<unsigned int> sigma;
    PetscFunctionBegin;

    if (argc == 1){ ierr=this->Usage(); CHKERRQ(ierr); }

    while(argc > 1){

        if ( (strcmp(argv[1],"-help") == 0)
            || (strcmp(argv[1],"-h") == 0)
            || (strcmp(argv[1],"-HELP") == 0) ){
            ierr=this->Usage(); CHKERRQ(ierr);
        }
        else if (strcmp(argv[1],"-advanced") == 0){
            ierr=this->Usage(true); CHKERRQ(ierr);
        }
        else if(strcmp(argv[1],"-nt") == 0){
            argc--; argv++;
            this->m_Domain.nt = static_cast<IntType>(atoi(argv[1]));
        }
        else if (strcmp(argv[1],"-pdesolver") == 0){
            argc--; argv++;
            if (strcmp(argv[1],"rk2") == 0){
                this->m_PDESolver.type = RK2;
            }
            if (strcmp(argv[1],"rk2a") == 0){
                this->m_PDESolver.type = RK2A;
            }
            else if (strcmp(argv[1],"sl") == 0){
                this->m_PDESolver.type = SL;
            }
            else{
                msg="\n\x1b[31m pde solver not implemented: %s\x1b[0m\n";
                ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
                ierr=this->Usage(true); CHKERRQ(ierr);
            }
        }
        else if(strcmp(argv[1],"-sigma") == 0){

            argc--; argv++;

            const std::string sigmainput = argv[1];

            // strip the "x" in the string to get the numbers
            sigma = String2Vec( sigmainput );

            if (sigma.size() == 1){
                for(unsigned int i=0; i < 3; ++i){
                    this->m_Sigma[i] = static_cast<ScalarType>(sigma[0]);
                }
            }
            else if(sigma.size() == 3){
                for(unsigned int i=0; i < 3; ++i){
                    this->m_Sigma[i] = static_cast<IntType>(sigma[i]);
                }
            }
            else{
                msg="\n\x1b[31m error in smoothing kernel size: %s\x1b[0m\n";
                ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
                ierr=this->Usage(); CHKERRQ(ierr);
            }

        }
        else if(strcmp(argv[1],"-disablesmoothing") == 0){
            this->m_RegFlags.smoothingenabled=false;
        }
        else if(strcmp(argv[1],"-nthreads") == 0){
            argc--; argv++;
            this->m_NumThreads = atoi(argv[1]);
        }
        else if(strcmp(argv[1],"-np") == 0){

            argc--; argv++;
            const std::string npinput = argv[1];

            // strip the "x" in the string to get the numbers
            np = String2Vec( npinput );

            if (np.size() == 1){
                for(unsigned int i=0; i < 2; ++i){
                    this->m_CartGridDims[i] = static_cast<unsigned int>(np[0]);
                }
            }
            else if (np.size() == 2){
                for(unsigned int i=0; i < 2; ++i){
                    this->m_CartGridDims[i] = static_cast<unsigned int>(np[i]);
                }
            }
            else{
                msg="\n\x1b[31m error in number of procs: %s\x1b[0m\n";
                ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
                ierr=this->Usage(); CHKERRQ(ierr);
            }
        }
        else if(strcmp(argv[1],"-x") == 0){
            argc--; argv++;
            this->m_ReadWriteFlags.xfolder = argv[1];
        }
        else if(strcmp(argv[1],"-i") == 0){
            argc--; argv++;
            this->m_ReadWriteFlags.ifolder = argv[1];
        }
        else if(strcmp(argv[1],"-xresults") == 0){
            this->m_ReadWriteFlags.results=true;
        }
        else if(strcmp(argv[1],"-xdefgrad") == 0){
            this->m_ReadWriteFlags.defgrad = true;
        }
        else if(strcmp(argv[1],"-xdefmap") == 0){
            this->m_ReadWriteFlags.defmap = true;
        }
        else if(strcmp(argv[1],"-xdeffield") == 0){
            this->m_ReadWriteFlags.deffield = true;
        }
        else if(strcmp(argv[1],"-xtimeseries") == 0){
            this->m_ReadWriteFlags.timeseries = true;
        }
        else if(strcmp(argv[1],"-detdefgradfromdeffield") == 0){
            this->m_RegFlags.detdefgradfromdeffield = true;
        }
        else if(strcmp(argv[1],"-invdefgrad") == 0){
            this->m_RegFlags.invdefgrad = true;
        }
        else if(strcmp(argv[1],"-ifile") == 0){
            argc--; argv++;
            this->m_iScaFieldFN = argv[1];
        }
        else if(strcmp(argv[1],"-ivecx1") == 0){
            argc--; argv++;
            this->m_iVecFieldX1FN = argv[1];
        }
        else if(strcmp(argv[1],"-ivecx2") == 0){
            argc--; argv++;
            this->m_iVecFieldX2FN = argv[1];
        }
        else if(strcmp(argv[1],"-ivecx3") == 0){
            argc--; argv++;
            this->m_iVecFieldX3FN = argv[1];
        }
        else if(strcmp(argv[1],"-rscale") == 0){
            argc--; argv++;
            this->m_ResamplingPara.gridscale = atof(argv[1]);
        }
        else if(strcmp(argv[1],"-resample") == 0){
            this->m_ResamplingPara.enabled = true;
        }
        else if(strcmp(argv[1],"-verbosity") == 0){
            argc--; argv++;
            this->m_Verbosity = atoi(argv[1]);
        }
        else{
            msg="\n\x1b[31m argument not valid: %s\x1b[0m\n";
            ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
            ierr=this->Usage(); CHKERRQ(ierr);
        }
        argc--; argv++;
    }

    // check the arguments/parameters set by the user
    ierr=this->CheckArguments(); CHKERRQ(ierr);

    // set number of threads
    ierr=Init(this->m_NumThreads,this->m_CartGridDims,this->m_FFT.mpicomm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegOpt"
RegToolsOpt::~RegToolsOpt()
{
    this->ClearMemory();
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode RegToolsOpt::ClearMemory()
{
    PetscFunctionBegin;

    if(this->m_FFT.plan!= NULL){
        accfft_destroy_plan(this->m_FFT.plan);
        accfft_cleanup();
        this->m_FFT.plan = NULL;
    }

    if (this->m_FFT.mpicomm != NULL){
        MPI_Comm_free(&this->m_FFT.mpicomm);
    }

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief initialize class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode RegToolsOpt::Initialize()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=this->SuperClass::Initialize(); CHKERRQ(ierr);

    this->m_RegToolsFlags.readvecfield = false;
    this->m_RegToolsFlags.readscafield = false;

    this->m_PostProcPara.enabled = false;
    this->m_PostProcPara.computedeffields = false;

    this->m_ResamplingPara.enabled = false;
    this->m_ResamplingPara.gridscale = -1.0;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief display usage message for binary
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Usage"
PetscErrorCode RegToolsOpt::Usage(bool advanced)
{

    PetscErrorCode ierr;
    int rank;
    std::string line;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    line = std::string(this->m_LineLength,'-');

    if (rank == 0){

        std::cout << std::endl;
        std::cout << line << std::endl;
        std::cout << " usage: regtools [options] " <<std::endl;
        std::cout << line << std::endl;
        std::cout << " where [options] is one or more of the following"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " ### compute measures from velocity field"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -i <path>                 input path (defines where registration results (i.e., velocity field,"<<std::endl;
        std::cout << "                           template image, and reference image) are stored; a prefix can be"<<std::endl;
        std::cout << "                           added by doing '-x </out/put/path/prefix_>"<<std::endl;
        std::cout << " -x <path>                 output path (by default only deformed template image and velocity"<<std::endl;
        std::cout << "                           field will be written; for more output options, see flags;"<<std::endl;
        std::cout << "                           a prefix can be added by doing '-x </out/put/path/prefix_>"<<std::endl;
        std::cout << " -xdefgrad                 flag: compute deformation gradient and write to file"<<std::endl;
        std::cout << " -xdefmap                  flag: compute deformation map and write to file"<<std::endl;
        std::cout << " -xdeffield                flag: compute displacement field and write to file"<<std::endl;
        std::cout << " -invdefgrad               flag: compute inverse deformation gradient"<<std::endl;


        // ####################### advanced options #######################
        if (advanced)
        {
        std::cout << line << std::endl;
        std::cout << " -sigma <int>x<int>x<int>  size of gaussian smoothing kernel applied to input images (e.g., 1x2x1;"<<std::endl;
        std::cout << "                           units: voxel size; if only one parameter is set"<<std::endl;
        std::cout << "                           uniform smoothing is assumed: default: 1x1x1)"<<std::endl;
        std::cout << " -disablesmoothing         flag: disable smoothing"<<std::endl;
        }
        // ####################### advanced options #######################

        // ####################### advanced options #######################
        if (advanced)
        {
        std::cout << line << std::endl;
        std::cout << " solver specific parameters (numerics)"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -pdesolver <type>         numerical time integrator for transport equations"<<std::endl;
        std::cout << "                           <type> is one of the following"<<std::endl;
        std::cout << "                               sl           semi-Lagrangian method (default; unconditionally stable)"<<std::endl;
        std::cout << "                               rk2          rk2 time integrator (conditionally stable)"<<std::endl;
        std::cout << " -nt <int>                 number of time points (for time integration; default: 4)"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " memory distribution and parallelism"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -nthreads <int>           number of threads (default: 1)"<<std::endl;
        std::cout << " -np <int>x<int>           distribution of mpi tasks (cartesian grid) (example: -np 2x4 results"<<std::endl;
        std::cout << "                           results in MPI distribution of size (nx1/2,nx2/4,nx3) for each mpi task)"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " other parameters/debugging"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -verbosity <int>          verbosity level (ranges from 0 to 3; default: 1)"<<std::endl;
        std::cout << " -xtimeseries              store time series (use with caution)"<<std::endl;
        std::cout << "                           problems; assumed to be uniform if single integer is provided"<<std::endl;
        }
        // ####################### advanced options #######################

        std::cout << line << std::endl;
        std::cout << " ### resampling"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -resample                 flag: resample data"<<std::endl;
        std::cout << " -rscale                   scale for resampling (multiplier applied to number of grid points)"<<std::endl;
        std::cout << " -ivecx1 <file>            x1 component of vector field (*.nii, *.nii.gz, *.hdr)"<<std::endl;
        std::cout << " -ivecx2 <file>            x2 component of vector field (*.nii, *.nii.gz, *.hdr)"<<std::endl;
        std::cout << " -ivecx3 <file>            x3 component of vector field (*.nii, *.nii.gz, *.hdr)"<<std::endl;
        std::cout << " -ifile <filename>         input file (image)"<<std::endl;
        std::cout << line << std::endl;
        std::cout << " -help                     display a brief version of the user message"<<std::endl;
        std::cout << " -advanced                 display this message"<<std::endl;
        std::cout << line << std::endl;
        std::cout << line << std::endl;

    }

    ierr=PetscFinalize(); CHKERRQ(ierr);
    exit(0);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief display options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DisplayOptions"
PetscErrorCode RegToolsOpt::DisplayOptions()
{
    PetscErrorCode ierr=0;
    int rank,indent;
    std::string msg,line;

    PetscFunctionBegin;

    this->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    indent = 40;
    line = std::string(this->m_LineLength,'-');

    // display the parameters (only on rank 0)
    if (rank == 0){

        std::cout<<std::endl;

        std::cout<< line << std::endl;
        std::cout<< " Constrained Large Deformation Diffeomorphic Registration"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " Parallel Algorithms for Data Analysis and Simulation Group"<<std::endl;
        std::cout<< " The Institute of Computational Engineering and Sciences"<<std::endl;
        std::cout<< " The University of Texas at Austin"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " problem setup"<<std::endl;
        std::cout<< line << std::endl;

        std::cout<< std::left << std::setw(indent) <<" problem dimensions"
                    << "(nx1,nx2,nx3,nt)=(" << this->m_Domain.nx[0] <<","
                    <<  this->m_Domain.nx[1] <<","
                    <<  this->m_Domain.nx[2] <<","
                    <<  this->m_Domain.nt <<")" <<std::endl;
        std::cout<< std::left << std::setw(indent) <<" network dimensions"
                    << this->m_CartGridDims[0] <<"x"
                    << this->m_CartGridDims[1]<<std::endl;
        std::cout<< std::left << std::setw(indent) <<" threads"
                    << this->m_NumThreads<<std::endl;
        std::cout<< std::left << std::setw(indent) <<" (ng,nl)"
                    << "(" << this->m_Domain.nglobal <<","
                    <<  this->m_Domain.nlocal <<")" <<std::endl;
        std::cout << line << std::endl;

    } // rank

    this->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief check the arguments set by user
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetVecFieldFN"
std::string RegToolsOpt::GetVecFieldFN(int i, int flag)
{

    if (flag == 0){
        if      (i == 0){ return this->m_iVecFieldX1FN; }
        else if (i == 1){ return this->m_iVecFieldX2FN; }
        else if (i == 2){ return this->m_iVecFieldX3FN; }
        else return "";
    }
    else if (flag == 1){
        if      (i == 0){ return this->m_xVecFieldX1FN; }
        else if (i == 1){ return this->m_xVecFieldX2FN; }
        else if (i == 2){ return this->m_xVecFieldX3FN; }
        else return "";
    }
    return "";
}




/********************************************************************
 * @brief check the arguments set by user
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetScaFieldFN"
std::string RegToolsOpt::GetScaFieldFN(int flag)
{
    if      (flag == 0){ return this->m_iScaFieldFN; }
    else if (flag == 1){ return this->m_xScaFieldFN; }
    return "";
}




/********************************************************************
 * @brief check the arguments set by user
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CheckArguments"
PetscErrorCode RegToolsOpt::CheckArguments()
{
    PetscErrorCode ierr;
    std::string msg,path,filename;
    size_t sep;
    PetscFunctionBegin;


    // check output arguments
    if (   this->m_ReadWriteFlags.defgrad
        || this->m_ReadWriteFlags.defmap
        || this->m_ReadWriteFlags.deffield ){

        if ( this->m_ReadWriteFlags.xfolder.empty() ){
            msg="\x1b[31m output folder needs to be set (-x option) \x1b[0m\n";
            ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr=this->Usage(); CHKERRQ(ierr);
        }

        if ( this->m_ReadWriteFlags.ifolder.empty() ){
            msg="\x1b[31m input folder needs to be set (-i option) \x1b[0m\n";
            ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr=this->Usage(); CHKERRQ(ierr);
        }

        this->m_PostProcPara.computedeffields = true;

        // set this flag to true, so that containers for reference and
        // template image are not to be deleted in registration class
        this->m_ReadWriteFlags.readfiles = true;

    }

    if (    !this->m_iVecFieldX1FN.empty()
         && !this->m_iVecFieldX2FN.empty()
         && !this->m_iVecFieldX3FN.empty() ){
        this->m_RegToolsFlags.readvecfield = true;
    }

    if ( !this->m_iScaFieldFN.empty() ){
        this->m_RegToolsFlags.readscafield = true;
    }

    if ( this->m_ResamplingPara.enabled ){

        if ( this->m_ResamplingPara.gridscale == -1.0 ){
            msg="\x1b[31m scale for rescaling needs to be set \x1b[0m\n";
            ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr=this->Usage(true); CHKERRQ(ierr);
        }

        if ( !this->m_RegToolsFlags.readvecfield && this->m_iScaFieldFN.empty() ){
            msg="\x1b[31m resampling requires input vector or scalar field \x1b[0m\n";
            ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
            ierr=this->Usage(true); CHKERRQ(ierr);
        }

        if ( this->m_RegToolsFlags.readvecfield ){

            sep = this->m_iVecFieldX1FN.find_last_of("\\/");
            if (sep != std::string::npos){
                path=this->m_iVecFieldX1FN.substr(0,sep);
                this->m_xVecFieldX1FN = path + "/resampled_" + this->m_iVecFieldX1FN.substr(sep + 1);
            }

            sep = this->m_iVecFieldX2FN.find_last_of("\\/");
            if (sep != std::string::npos){
                path=this->m_iVecFieldX2FN.substr(0,sep);
                this->m_xVecFieldX2FN = path + "/resampled_" + this->m_iVecFieldX2FN.substr(sep + 1);
            }

            sep = this->m_iVecFieldX3FN.find_last_of("\\/");
            if (sep != std::string::npos){
                path=this->m_iVecFieldX3FN.substr(0,sep);
                this->m_xVecFieldX3FN = path + "/resampled_" + this->m_iVecFieldX3FN.substr(sep + 1);
            }

        }

        if ( this->m_RegToolsFlags.readscafield ){
            sep = this->m_iScaFieldFN.find_last_of("\\/");
            if (sep != std::string::npos){
                path=this->m_iScaFieldFN.substr(0,sep);
                this->m_xScaFieldFN = path + "/resampled_" + this->m_iScaFieldFN.substr(sep + 1);
            }
        }

    }

    ierr=Assert(this->m_NumThreads > 0,"omp threads < 0"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of namespace

#endif //_REGTOOLSOPT_CPP_
