#ifndef _REGCOMMANDOPTIONS_CPP_
#define _REGCOMMANDOPTIONS_CPP_

#include "RegCommandOptions.h"


namespace reg
{


/********************************************************************
 * Name: RegCommandOptions
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegCommandOptions"
RegCommandOptions::RegCommandOptions()
{
    this->Initialize();
}




/********************************************************************
 * Name: OptimalControlRegistration
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegCommandOptions"
RegCommandOptions::~RegCommandOptions(void)
{
    this->ClearMemory();
}




/********************************************************************
 * Name: RegCommandOptions
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegCommandOptions"
RegCommandOptions::RegCommandOptions(int argc, char**argv)
{
    this->Initialize();


    while(argc > 1){

        if(    (strcmp(argv[1],"-help") == 0)
            || (strcmp(argv[1],"-HELP") == 0) ){
            this->Usage();
        }
        else if (strcmp(argv[1],"-h") == 0){
            this->BriefUsage();
        }
        else if (strcmp(argv[1],"-r") == 0){
            argc--; argv++;
            this->m_opt->mRfn = argv[1];
        }
        else if (strcmp(argv[1],"-t") == 0){
            argc--; argv++;
            this->m_opt->mTfn = argv[1];
        }
        else if (strcmp(argv[1],"-x") == 0){
            argc--; argv++;
            this->m_opt->xfolder = argv[1];
        }
        else if(strcmp(argv[1],"-nx") == 0){
            argc--; argv++;
            nx = atoi(argv[1]);
            for(int i=0; i < 3; ++i){
                this->m_opt->nx[i] = nx;
            }
        }
        else {
            MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
            msg="\x1b[31m argument not valid: %s\x1b[0m\n";
            ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
            this->Usage();
        }

        argc--; argv++;
    }

    this->ParseParameters();
}




/********************************************************************
 * Name: Initialize
 * Description: init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode RegCommandOptions::Initialize(void)
{
    PetscFunctionBegin;

    this->m_RegOpt=NULL;
    this->m_MISCOpt=NULL;
    this->m_ACCFFTPlan=NULL;
    this->m_MPIComm=NULL;

    this->m_nx[0] = 32;
    this->m_nx[1] = 32;
    this->m_nx[2] = 32;

    this->m_nt = 4;


    this->m_Beta = 1E-2;
    this->m_JacBound = 2E-1;

    this->m_OptMaxIt = 50;
    this->m_KKTMaxIt = 1E3;
    this->m_Verbosity = 1;

    this->m_NumThreads = 1;
    this->m_CartGridDims[0] = 1;
    this->m_CartGridDims[1] = 1;

    this->m_xfolder = "./results/";
    this->m_mRFN = "/";
    this->m_mTFN = "/";

    this->m_StoreImages = false;
    this->m_DoParameterContinuation = true;
    this->m_Incompressible = false;
    this->m_StoreTimeSeries = false;

    this->m_PDESolverType = SL;
    this->m_RegularizationNorm = H2;
    this->m_OptimizationMethod = GAUSSNEWTON;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode RegCommandOptions::ClearMemory(void)
{
//    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_MISCOpt != NULL){
        delete this->m_MISCOpt;
        this->m_MISCOpt = NULL;
    }

    if (this->m_RegOpt != NULL){
        delete this->m_RegOpt;
        this->m_RegOpt = NULL;
    }

    if(this->m_ACCFFTPlan!= NULL){
        accfft_destroy_plan(this->m_ACCFFTPlan);
        accfft_cleanup();
        this->m_ACCFFTPlan = NULL;
    }

    if (this->m_MPIComm != NULL){
        MPI_Comm_free(&this->m_MPIComm);
    }

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: DoSetup
 * Description: setup options and accfft
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DoSetup"
PetscErrorCode RegCommandOptions::DoSetup()
{
    PetscErrorCode ierr;
    int nx[3],isize[3],istart[3],osize[3],ostart[3],ompthreads;
    IntType alloc_max;
    ScalarType *u, setuptime;
    Complex *uk;
    std::stringstream ss;
    PetscFunctionBegin;

    nx[0] = this->m_nx[0];
    nx[1] = this->m_nx[1];
    nx[2] = this->m_nx[2];

    ierr=reg::Assert(this->m_NumThreads > 0,"omp threads < 0"); CHKERRQ(ierr);

    omp_set_dynamic(0);
    omp_set_num_threads(this->m_NumThreads);

    // check if number of threads is consistent with user options
    ompthreads=omp_get_max_threads();

    ss.str( std::string() );
    ss.clear();
    ss << "max number of openmp threads is not a match (user,set)=("
       << this->m_NumThreads <<"," << ompthreads <<")\n";
    ierr=Assert(ompthreads == this->m_NumThreads,ss.str().c_str()); CHKERRQ(ierr);

    // set up the fft
    fftsetuptime=-MPI_Wtime();
    this->m_ACCFFTPlan = accfft_plan_dft_3d_r2c(nx,u,(double*)uk,this->m_MPIComm,ACCFFT_MEASURE);
    fftsetuptime+=MPI_Wtime();

    // initialize accft
    accfft_create_comm(PETSC_COMM_WORLD,this->m_CartGridDims,&this->m_MPIComm);
    accfft_init(this->m_NumThreads);

    alloc_max = accfft_local_size_dft_r2c(nx,isize,istart,osize,ostart,this->m_MPIComm);
    u  = (ScalarType*)accfft_alloc(alloc_max);
    uk = (Complex*)accfft_alloc(alloc_max);


    try{ this->m_MISCOpt = new N_MISC(nx,isize,istart,this->m_ACCFFTPlan,this->m_MPIComm); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ierr=reg::Assert(this->m_MISCOpt->N_local > 0,"bug in setup"); CHKERRQ(ierr);
    ierr=reg::Assert(this->m_MISCOpt->N_global > 0,"bug in setup"); CHKERRQ(ierr);

    // allocate class for registration options
    try{ *regopt = new reg::RegOpt(*this->m_MISCOpt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ParseParameters
 * Description: parse the set parameters to option class
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ParseArguments"
PetscErrorCode RegCommandOptions::ParseArguments()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=this->DoSetup(); CHKERRQ(ierr);

    this->m_RegOpt->SetNumTimePoints(this->m_nt);
    this->m_RegOpt->SetRegNorm(this->m_RegularizationNorm);
    this->m_RegOpt->SetOptMeth(this->m_OptimizationMethod);
    this->m_RegOpt->SetPDESolver(this->m_PDESolverType);
    this->m_RegOpt->SetOptMaxit(this->m_OptMaxIt);
    this->m_RegOpt->SetKKTMaxit(this->m_KKTMaxIt);
    for(int i=0; i<3; ++i) this->m_RegOpt->SetOptTol(i,this->m_OptTol[i]);

    this->m_RegOpt->SetJacBound(this->m_JacBound);
    this->m_RegOpt->SetRegularizationWeight(this->m_Beta);
    this->m_RegOpt->SetNumThreads(this->m_NumThreads);
    this->m_RegOpt->SetNetworkDims(this->m_CartGridDims);

    this->m_RegOpt->InCompressible(this->m_Incompressible);
    this->m_RegOpt->SetXFolder(this->m_xfolder);
    this->m_RegOpt->WriteImages2File(this->m_StoreImages);
    this->m_RegOpt->StoreTimeSeries(this->m_StoreTimeSeries);
    this->m_RegOpt->DoParameterContinuation(this->m_DoParameterContinuation);

    this->m_RegOpt->SetVerbosity(this->m_Verbosity);

    this->m_RegOpt->DisplayOptions();

    PetscFunctionReturn(0);

}




} // end of name space







#endif // _REGCOMMANDOPTIONS_CPP_
