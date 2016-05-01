
#ifndef _REGOPT_CPP_
#define _REGOPT_CPP_

#include <fstream>

#include "RegOpt.h"

namespace reg
{


/********************************************************************
 * Name: RegOpt
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegOpt"
RegOpt::RegOpt()
{
    this->Initialize();
}



/********************************************************************
 * Name: RegOpt
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegOpt"
RegOpt::RegOpt(N_MISC* opt)
{
    this->Initialize();
    this->m_MiscOpt = opt;
}



/********************************************************************
 * Name: RegOpt
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegOpt"
RegOpt::RegOpt(int argc, char** argv)
{
    this->Initialize();
    this->ParseArguments(argc,argv);
}





/********************************************************************
 * Name: ParseArguments
 * Description: parse user arguments
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ParseArguments"
PetscErrorCode RegOpt::ParseArguments(int argc, char** argv)
{
    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    while(argc > 1){

        if ( (strcmp(argv[1],"-help") == 0)
            ||(strcmp(argv[1],"-h") == 0)
            ||(strcmp(argv[1],"-HELP") == 0) ){
            ierr=this->Usage(); CHKERRQ(ierr);
        }
        else if(strcmp(argv[1],"-nx") == 0){
            argc--; argv++;
            int nx = atoi(argv[1]);
            for(int i=0; i < 3; ++i){
                this->m_nx[i] = static_cast<IntType>(nx);
            }
        }
        else if(strcmp(argv[1],"-nx1") == 0){
            argc--; argv++;
            this->m_nx[0] = static_cast<IntType>(atoi(argv[1]));
        }
        else if(strcmp(argv[1],"-nx2") == 0){
            argc--; argv++;
            this->m_nx[1] = static_cast<IntType>(atoi(argv[1]));
        }
        else if(strcmp(argv[1],"-nx3") == 0){
            argc--; argv++;
            this->m_nx[2] = static_cast<IntType>(atoi(argv[1]));
        }
        else if(strcmp(argv[1],"-nt") == 0){
            argc--; argv++;
            this->m_NT = static_cast<IntType>(atoi(argv[1]));
        }
        else if(strcmp(argv[1],"-nthreads") == 0){
            argc--; argv++;
            this->m_NumThreads = atoi(argv[1]);
        }
        else if(strcmp(argv[1],"-npx1") == 0){
            argc--; argv++;
            this->m_CartGridDims[0] = atoi(argv[1]);
        }
        else if(strcmp(argv[1],"-npx2") == 0){
            argc--; argv++;
            this->m_CartGridDims[0] = atoi(argv[1]);
        }
        else if(strcmp(argv[1],"-mr") == 0){
            argc--; argv++;
            this->m_ReferenceFN = argv[1];
        }
        else if(strcmp(argv[1],"-mt") == 0){
            argc--; argv++;
            this->m_TemplateFN = argv[1];
        }
        else if(strcmp(argv[1],"-x") == 0){
            argc--; argv++;
            this->m_XFolder = argv[1];
        }
        else if(strcmp(argv[1],"-ic") == 0){
            argc--; argv++;
            this->m_InCompressible = true;
        }
        else if(strcmp(argv[1],"-ximages") == 0){
            argc--; argv++;
            this->m_WriteImages = true;
        }
        else if(strcmp(argv[1],"-xlog") == 0){
            argc--; argv++;
            this->m_WriteLogFiles = true;
        }
        else if (strcmp(argv[1],"-optmeth") == 0){
            argc--; argv++;
            if (strcmp(argv[1],"fn") == 0){
                this->m_OptMeth = FULLNEWTON;
            }
            else if (strcmp(argv[1],"gn") == 0){
                this->m_OptMeth = GAUSSNEWTON;
            }
            else {
                msg="\x1b[31m optimization method not defined: %s\x1b[0m\n";
                ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
                ierr=this->Usage(); CHKERRQ(ierr);
            }
        }
        else if(strcmp(argv[1],"-maxit") == 0){
            argc--; argv++;
            this->m_OptPara.maxit = atoi(argv[1]);
        }
        else if(strcmp(argv[1],"-gabs") == 0){
            argc--; argv++;
            this->m_OptPara.tol[0] = atof(argv[1]);
        }
        else if(strcmp(argv[1],"-grel") == 0){
            argc--; argv++;
            this->m_OptPara.tol[2] = atof(argv[1]);
        }
        else if(strcmp(argv[1],"-jbound") == 0){
            argc--; argv++;
            this->m_ParameterCont.jacbound = atof(argv[1]);
        }
        else if(strcmp(argv[1],"-kktmaxit") == 0){
            argc--; argv++;
            this->m_KKTSolverPara.maxit = atoi(argv[1]);
        }
        else if(strcmp(argv[1],"-kktfseq") == 0){
            argc--; argv++;
            if (strcmp(argv[1],"none") == 0){
                this->m_KKTSolverPara.fseqtype = NOFS;
            }
            else if (strcmp(argv[1],"quadratic") == 0){
                this->m_KKTSolverPara.fseqtype = QDFS;
            }
            else if (strcmp(argv[1],"suplinear") == 0){
                this->m_KKTSolverPara.fseqtype = SLFS;
            }
            else {
                msg="\x1b[31m optimization method not defined: %s\x1b[0m\n";
                ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
                ierr=this->Usage(); CHKERRQ(ierr);
            }
        }
        else if(strcmp(argv[1],"-x") == 0){
            argc--; argv++;
            this->m_XFolder = argv[1];
        }
        else if (strcmp(argv[1],"-pdesolver") == 0){
            argc--; argv++;
            if (strcmp(argv[1],"rk2") == 0){
                this->m_PDESolver = RK2;
            }
            else if (strcmp(argv[1],"sl") == 0){
                this->m_PDESolver = SL;
            }
            else {
                msg="\x1b[31m pde solver not implemented: %s\x1b[0m\n";
                ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
                ierr=this->Usage(); CHKERRQ(ierr);
            }
        }
        else if (strcmp(argv[1],"-regnorm") == 0){
            argc--; argv++;
            if (strcmp(argv[1],"h1s") == 0){
                this->m_RegNorm = H1SN;
            }
            else if (strcmp(argv[1],"h2s") == 0){
                this->m_RegNorm = H2SN;
            }
            else if (strcmp(argv[1],"h1") == 0){
                this->m_RegNorm = H1;
            }
            else if (strcmp(argv[1],"h2") == 0){
                this->m_RegNorm = H2;
            }
            else if (strcmp(argv[1],"l2") == 0){
                this->m_RegNorm = L2;
            }
            else {
                msg="\x1b[31m regularization norm not available: %s\x1b[0m\n";
                ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
                ierr=this->Usage(); CHKERRQ(ierr);
            }
        }
        else if(strcmp(argv[1],"-beta") == 0){
            argc--; argv++;
            this->m_Beta[0] = atof(argv[1]);
            this->m_Beta[1] = atof(argv[1]);
        }
        else if(strcmp(argv[1],"-estbeta") == 0){
            argc--; argv++;
            this->m_ParameterCont.enabled=true;
        }
        else if(strcmp(argv[1],"-verbosity") == 0){
            argc--; argv++;
            this->m_Verbosity = atoi(argv[1]);
        }
        else {
            msg="\x1b[31m argument not valid: %s\x1b[0m\n";
            ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str(),argv[1]); CHKERRQ(ierr);
            ierr=this->Usage(); CHKERRQ(ierr);
        }
        argc--; argv++;
    }

    ierr=this->CheckArguments(); CHKERRQ(ierr);

    ierr=this->DoSetup(); CHKERRQ(ierr);
    ierr=this->DisplayOptions(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}





/********************************************************************
 * Name: RegOpt
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegOpt"
RegOpt::~RegOpt()
{
    this->ClearMemory();
}




/********************************************************************
 * Name: ClearMemory
 * Description: clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode RegOpt::ClearMemory()
{
    //PetscErrorCode ierr;

    PetscFunctionBegin;

    if(this->m_MiscOpt!= NULL){
        delete this->m_MiscOpt;
        this->m_MiscOpt = NULL;
    }

    if(this->m_FFTPlan!= NULL){
        accfft_destroy_plan(this->m_FFTPlan);
        accfft_cleanup();
        this->m_FFTPlan = NULL;
    }

    if (this->m_Comm != NULL){
        MPI_Comm_free(&this->m_Comm);
    }

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: Initialize
 * Description: initialize class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode RegOpt::Initialize()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_MiscOpt = NULL;
    this->m_FFTPlan = NULL;
    this->m_Comm = NULL;

    this->m_NT = 4;
    this->m_nx[0] = 32;
    this->m_nx[1] = 32;
    this->m_nx[2] = 32;

    this->m_Beta[0] = 1E-2;
    this->m_Beta[1] = 1E-2;
    this->m_RegNorm = H2SN;
    this->m_Verbosity = 1;
    this->m_TimeHorizon[0] = 0.0;
    this->m_TimeHorizon[1] = 1.0;
    this->m_OptMeth = GAUSSNEWTON;
    this->m_PDESolver = SL;
    this->m_PrecondMeth = INVREG;
    //this->m_PrecondMeth = TWOLEVEL;
    //this->m_PCSolverType = PCPCG;
    this->m_XFolder = "./results/";
    this->m_IFolder = "./results/";
    this->m_CartGridDims[0] = 1;
    this->m_CartGridDims[1] = 1;
    this->m_LineLength = 100;
    this->m_Sigma = 1;
    this->m_WriteImages = false;
    this->m_WriteLogFiles = false;

    this->m_KKTSolverPara.tol[0] = 1E-12; // relative tolerance
    this->m_KKTSolverPara.tol[1] = 1E-12; // absolute tolerance
    this->m_KKTSolverPara.tol[2] = 1E+06; // divergence tolerance
    this->m_KKTSolverPara.maxit  = 1E+03; // maximal iterations
    this->m_KKTSolverPara.fseqtype = QDFS;

    this->m_OptPara.tol[0] = 1E-6;  // grad abs tol
    this->m_OptPara.tol[1] = 1E-16; // grad rel tol
    this->m_OptPara.tol[2] = 1E-2;  // grad rel tol
    this->m_OptPara.maxit = 1E3; // max number of iterations

    this->m_DD.n = 2;

    this->m_MonitorCFLCondition = false;
    this->m_ReadImagesFromFile = false;
    this->m_InCompressible = true;

    this->m_ParameterCont.enabled = false;
    this->m_ParameterCont.betamin = 1E-6;
    this->m_ParameterCont.jacbound = 8E-1;

    this->m_RegMonitor.jacmin = 0.0;
    this->m_RegMonitor.jacmax = 0.0;
    this->m_RegMonitor.jacmean = 0.0;


    this->m_NumThreads=1;
    this->m_CartGridDims[0]=1;
    this->m_CartGridDims[1]=1;

    ierr=this->ResetTimers(); CHKERRQ(ierr);
    ierr=this->ResetCounters(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: Usage
 * Description: display usage message for binary
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Usage"
PetscErrorCode RegOpt::Usage()
{

    PetscErrorCode ierr;
    int rank;
    std::string line;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    line = std::string(this->m_LineLength,'-');

    if (rank == 0){
        std::cout<<std::endl;
        std::cout<< " usage: run [options] " <<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " where [options] is one or more of the following"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< line << std::endl;
        std::cout<< " -mr <file>         reference image (*.nii, *.nc)"<<std::endl;
        std::cout<< " -mt <file>         template image (*.nii, *.nc)"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " -x <folder>        output folder (what's written out is controlled by the flags below)"<<std::endl;
        std::cout<< " -ximages           flag: write images to file (default: not written; requires -x option)"<<std::endl;
        std::cout<< " -xlog              flag: write log files (default: not written; requires -x option)"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " # problem specific parameters"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " -nx <int>          grid size (in all directions uniform; default 32)"<<std::endl;
        std::cout<< " -nx1 <int>         grid size (x1 direction; default 32)"<<std::endl;
        std::cout<< " -nx2 <int>         grid size (x2 direction; default 32)"<<std::endl;
        std::cout<< " -nx3 <int>         grid size (x3 direction; default 32)"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " # optimization specific parameters"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " -optmeth <type>    optimization method (default: gn)"<<std::endl;
        std::cout<< "                    where <types> are"<<std::endl;
        std::cout<< "                      gn           Gauss-Newton (default)"<<std::endl;
        std::cout<< "                      fn           full Newton"<<std::endl;
        std::cout<< " -grel <dbl>        tolerance for optimization (default: 1E-2)"<<std::endl;
        std::cout<<"                       relative change of gradient"<<std::endl;
        std::cout<<"                       optimization stops if ||g_k||/||g_0|| <= tol"<<std::endl;
        std::cout<< " -gabs <dbl>        tolerance for optimization (default: 1E-6)"<<std::endl;
        std::cout<<"                       lower bound for gradient"<<std::endl;
        std::cout<<"                       optimization stops if ||g_k|| <= tol"<<std::endl;
        std::cout<<"                       tol <= ||g||/||g_init||"<<std::endl;
        std::cout<< " -maxit <int>       maximum number of (outer) iterations (default: 50)"<<std::endl;
        std::cout<< " -kktmaxit <int>    maximum number of (inner) iterations (KKT solver; default: 50)"<<std::endl;
        std::cout<< " -kktfseq <type>    forcing sequence for KKT solve (inner iterations; default: quad)"<<std::endl;
        std::cout<< "                    where <types> are"<<std::endl;
        std::cout<< "                      none          exact solve (expensive)"<<std::endl;
        std::cout<< "                      quadratic     quadratic (default)"<<std::endl;
        std::cout<< "                      suplinear     super-linear"<<std::endl;
        std::cout<< " -jbound <dbl>      lower bound on determinant of deformation gradient (default: 2E-1)"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " # regularization parameters"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " -regnorm <type>    regularization norm (default: h2s)"<<std::endl;
        std::cout<< "                    where <types> are"<<std::endl;
        std::cout<< "                      h1s          H1-seminorm"<<std::endl;
        std::cout<< "                      h2s          H2-seminorm (default)"<<std::endl;
        std::cout<< "                      h1           H1-norm"<<std::endl;
        std::cout<< "                      h2           H2-norm"<<std::endl;
        std::cout<< "                      l2           l2-norm (discouraged)"<<std::endl;
        std::cout<< " -beta <dbl>        regularization weight/parameter (default: 1E-2)"<<std::endl;
        std::cout<< " -estbeta           estimate regularization parameter (default: not enabled)"<<std::endl;
        std::cout<< " -ic                flag: enable incompressibility constraint (default: not enabled)"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " # solver specific parameters (numerics)"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " -pdesolver <type>  numerical time integrator for transport equations (default: sl)"<<std::endl;
        std::cout<< "                    where <types> are"<<std::endl;
        std::cout<< "                      sl      semi-Lagrangian method"<<std::endl;
        std::cout<< "                      rk2     rk2 time integrator"<<std::endl;
        std::cout<< " -nt <int>          number of time points (for time integration; default: 4)"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " # memory distribution and parallelism"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " -nthreads <int>    number of threads (default: 1)"<<std::endl;
        std::cout<< " -npx1 <int>        MPI parameter; number of procs in x1 direction (default: 1)"<<std::endl;
        std::cout<< " -npx2 <int>        MPI parameter; number of procs in x2 direction (default: 1)"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " # other parameters"<<std::endl;
        std::cout<< line << std::endl;
        std::cout<< " -verbosity <int>   verbosity level (ranges from 0 to 3; default: 1)"<<std::endl;
        std::cout<< line << std::endl;
    }

    ierr=PetscFinalize(); CHKERRQ(ierr);
    exit(0);

    PetscFunctionReturn(0);


}




/********************************************************************
 * Name: CheckArguments
 * Description: check the arguments set by user
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CheckArguments"
PetscErrorCode RegOpt::CheckArguments()
{
    PetscErrorCode ierr;
    bool readmR,readmT;
    std::string msg;
    PetscFunctionBegin;

    readmR=false;
    if(!this->m_ReferenceFN.empty()){
        readmR=true;
    }
    readmT=false;
    if(!this->m_TemplateFN.empty()){
        readmT=true;
    }

    if (readmT && readmR){
        this->ReadImagesFromFile(true);
    }
    else if( (readmT == false) && readmR ) {
        msg="\x1b[31m you need to also assign a template image\x1b[0m\n";
        ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
        ierr=this->Usage(); CHKERRQ(ierr);
    }
    else if( readmT && (readmR == false) ) {
        msg="\x1b[31m you need to also assign a reference image\x1b[0m\n";
        ierr=PetscPrintf(PETSC_COMM_WORLD,msg.c_str()); CHKERRQ(ierr);
        ierr=this->Usage(); CHKERRQ(ierr);
    }
    else if( (readmT == false) && (readmR == false) ){
        this->ReadImagesFromFile(false);
    }

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: DoSetup
 * Description: setup options and accfft
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DoSetup"
PetscErrorCode RegOpt::DoSetup()
{
    PetscErrorCode ierr;
    int nx[3],isize[3],istart[3],osize[3],ostart[3],ompthreads;
    IntType alloc_max;
    ScalarType *u, fftsetuptime;
    Complex *uk;
    std::stringstream ss;
    PetscFunctionBegin;

    nx[0] = static_cast<int>(this->m_nx[0]);
    nx[1] = static_cast<int>(this->m_nx[1]);
    nx[2] = static_cast<int>(this->m_nx[2]);

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

    // initialize accft
    accfft_create_comm(PETSC_COMM_WORLD,this->m_CartGridDims,&this->m_Comm);
    accfft_init(this->m_NumThreads);

    alloc_max = accfft_local_size_dft_r2c(nx,isize,istart,osize,ostart,this->m_Comm);
    u  = (ScalarType*)accfft_alloc(alloc_max);
    uk = (Complex*)accfft_alloc(alloc_max);

    // set up the fft
    fftsetuptime=-MPI_Wtime();
    this->m_FFTPlan = accfft_plan_dft_3d_r2c(nx,u,(double*)uk,this->m_Comm,ACCFFT_MEASURE);
    fftsetuptime+=MPI_Wtime();

    // set the fft setup time
    this->m_Timer[FFTSETUP][LOG] = fftsetuptime;

    try{ this->m_MiscOpt = new N_MISC(nx,isize,istart,this->m_FFTPlan,this->m_Comm); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ierr=reg::Assert(this->m_MiscOpt->N_local > 0,"bug in setup"); CHKERRQ(ierr);
    ierr=reg::Assert(this->m_MiscOpt->N_global > 0,"bug in setup"); CHKERRQ(ierr);

    accfft_free(u);
    accfft_free(uk);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: DisplayOptions
 * Description: display registration options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DisplayOptions"
PetscErrorCode RegOpt::DisplayOptions()
{
    PetscErrorCode ierr;
    int rank,indent,align;
    std::string msg, line;
    bool newtontype;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    indent = 40;
    align  = 30;
    line = std::string(this->m_LineLength,'-');

    // display the parameters (only on rank 0)
    if (rank == 0){

        std::cout<<std::endl;

        std::cout<< line << std::endl;
        std::cout<< " problem setup"<<std::endl;
        std::cout<< line << std::endl;

        std::cout<< std::left << std::setw(indent) <<" problem dimensions"
                    << "(" << this->m_MiscOpt->N[0] <<","
                    <<  this->m_MiscOpt->N[1] <<","
                    <<  this->m_MiscOpt->N[2] <<","
                    <<  this->m_NT <<")" <<std::endl;
        std::cout<< std::left << std::setw(indent) <<" network dimensions"
                    << this->m_CartGridDims[0] <<"x"
                    << this->m_CartGridDims[1]<<std::endl;
        std::cout<< std::left << std::setw(indent) <<" threads"
                    << this->m_NumThreads<<std::endl;
        std::cout<< std::left << std::setw(indent) <<" (ng,nl)"
                    << "(" << this->GetNGlobal() <<","
                    <<  this->GetNLocal() <<")" <<std::endl;

        std::cout<< line << std::endl;
        std::cout<< " parameters"<<std::endl;
        std::cout<< line << std::endl;

        // display regularization model
        std::cout<< std::left << std::setw(indent) <<" regularization model";

        if(this->m_ParameterCont.enabled){
            switch(this->m_RegNorm){
                case L2: std::cout<<"L2-norm (beta estimated)"<<std::endl; break;
                case H1: std::cout<<"H1-norm (beta estimated)"<<std::endl; break;
                case H2: std::cout<<"H2-norm (beta estimated)"<<std::endl; break;
                case H1SN: std::cout<<"H1-seminorm (beta estimated)"<<std::endl; break;
                case H2SN: std::cout<<"H2-seminorm (beta estimated)"<<std::endl; break;
                default:{ ierr=ThrowError("regularization model not implemented"); CHKERRQ(ierr); break; }
            }

            // display parameters and tolerances
            std::cout<< std::left << std::setw(indent) <<" parameter continuation"
                      << std::setw(align) <<"bound det(grad(y))"
                      << this->m_ParameterCont.jacbound << std::endl;

            std::cout << std::left << std::setw(indent) <<" "
                      << std::setw(align) <<"bound beta"
                      << this->m_ParameterCont.betamin << std::endl;
        }
        else{
            switch(this->m_RegNorm){
                case L2: std::cout<<"L2-norm (beta="<<this->m_Beta[0]<<")"<<std::endl; break;
                case H1: std::cout<<"H1-norm (beta="<<this->m_Beta[0]<< ", "<<this->m_Beta[1]<<")"<<std::endl; break;
                case H2: std::cout<<"H2-norm (beta="<<this->m_Beta[0]<< ", "<<this->m_Beta[1]<<")"<<std::endl; break;
                case H1SN: std::cout<<"H1-seminorm (beta="<<this->m_Beta[0]<<")"<<std::endl; break;
                case H2SN: std::cout<<"H2-seminorm (beta="<<this->m_Beta[0]<<")"<<std::endl; break;
                default:{ ierr=ThrowError("regularization model not implemented"); CHKERRQ(ierr); break; }
            }
        }

        // display regularization model
        std::cout<< std::left << std::setw(indent) <<" PDE solver (hyperbolic)";
        switch(this->m_PDESolver){
            case RK2: std::cout<<"second order RK method"<<std::endl; break;
            case SL: std::cout<<"semi-lagrangian method"<<std::endl; break;
            default:
            {
                ierr=ThrowError("solver not implemented"); CHKERRQ(ierr);
                break;
            }
        }
        // display type of optimization method
        newtontype=false;
        std::cout<< std::left << std::setw(indent) <<" optimization method";
        switch(this->m_OptMeth){
//            case LMVM:        {std::cout<< "LMVM method" << std::endl; break;}
            case GRADDESCENT: {std::cout<< "preconditioned gradient descent method" << std::endl; break;}
            case FULLNEWTON:  {std::cout<< "full newton method" << std::endl; newtontype=true; break;}
            case GAUSSNEWTON: {std::cout<< "gauss newton method" << std::endl; newtontype=true; break;}
            default:
            {
                ierr=ThrowError("optimization method not implemented"); CHKERRQ(ierr);
                break;
            }
        }

        std::cout << std::left << std::setw(indent) <<" maximal # iterations"
                  << this->m_OptPara.maxit << std::endl;

        // display optimization tolerances
        std::cout << std::left << std::setw(indent) <<" convergence tolerances"
                  << std::setw(align) <<"||g(v)|| <= tol"
                  << this->m_OptPara.tol[0] << std::endl;
        std::cout << std::left << std::setw(indent) <<" "
                  << std::setw(align) <<"||g(v)||/|J(v)| <= tol"
                  << this->m_OptPara.tol[1] << std::endl;
        std::cout << std::left << std::setw(indent) <<" "
                  << std::setw(align) <<"||g(v)||/||g(v0)|| <= tol"
                  << this->m_OptPara.tol[2] << std::endl;

        // display parameters for newton type optimization methods
        if ( newtontype ){

            std::cout << std::left << std::setw(indent)
                      <<" KKT sytem"
                      << std::setw(align) << "solver"
                      << "PCG" << std::endl;

            std::cout << std::left << std::setw(indent) <<" "
                      << std::setw(align) << "forcing sequence";

            switch(this->m_KKTSolverPara.fseqtype){
                case NOFS:
                {
                    std::cout << std::setw(align) << "disabled" << std::endl;
                    std::cout << std::left << std::setw(indent) <<" "
                              << std::setw(align) <<"absolute"
                              << this->m_KKTSolverPara.tol[0] << std::endl;
                    std::cout << std::left << std::setw(indent) <<" "
                              << std::setw(align) <<"relative"
                              << this->m_KKTSolverPara.tol[1] << std::endl;
                    break;
                }
                case QDFS: { std::cout << "quadratic" << std::endl; break; }
                case SLFS:{ std::cout << "superlinear" << std::endl; break; }
                default:
                {
                    ierr=ThrowError("forcing sequence not implemented"); CHKERRQ(ierr);
                    break;
                }
            }
            std::cout << std::left << std::setw(indent) <<" "
                      << std::setw(align) <<"divergence"
                      << this->m_KKTSolverPara.tol[2] << std::endl;

            std::cout << std::left << std::setw(indent) <<" "
                      << std::setw(align) <<"maxit"
                      << this->m_KKTSolverPara.maxit << std::endl;

        }

        std::cout<< line << std::endl;

    }

    PetscFunctionReturn(0);
}





/********************************************************************
 * Name: ComputeFFTScale
 * Description: compute weight for FFT
 *******************************************************************/
ScalarType RegOpt::ComputeFFTScale()
{

    ScalarType scale = 1.0;
    for (unsigned int i=0; i < 3; ++i){
        scale *= static_cast<ScalarType>(this->m_MiscOpt->N[i]);
    }
    return 1.0/scale;

};





/********************************************************************
 * Name: ResetTimers
 * Description: resets all timers
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResetTimers"
PetscErrorCode RegOpt::ResetTimers()
{
    PetscFunctionBegin;

    for (int i = 0; i < NTIMERS; ++i){
        if (i != FFTSETUP){
            for (int j = 0; j < NVALTYPES; ++j){
                this->m_Timer[i][j] = 0.0;
            }
            this->m_TimerIsRunning[i] = false;
            this->m_TempTimer[i] = 0.0;
        }
    }

    for(int i = 0; i < 5; ++i){
        for (int j = 0; j < NVALTYPES; ++j){
            this->m_FFTTimers[i][j] = 0.0;
        }
    }

    for(int i = 0; i < 4; ++i){
        for (int j = 0; j < NVALTYPES; ++j){
            this->m_InterpTimers[i][j] = 0.0;
        }
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ResetTimer
 * Description: resets timer
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResetTimer"
PetscErrorCode RegOpt::ResetTimer(TimerType id)
{
    PetscFunctionBegin;

    for (int j = 0; j < NVALTYPES; ++j){
        this->m_Timer[id][j] = 0.0;
    }
    this->m_TimerIsRunning[id] = false;
    this->m_TempTimer[id] = 0.0;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: StartTimer
 * Description: start the timer (checks if running)
 * Author: Andreas Mang
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "StartTimer"
PetscErrorCode RegOpt::StartTimer(TimerType id)
{

    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    msg="fatal error: timer has already been started";
    ierr=Assert(this->m_TimerIsRunning[id]==false,msg); CHKERRQ(ierr);

    this->m_TempTimer[id] = -MPI_Wtime();
    this->m_TimerIsRunning[id] = true;

    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: StopTimer
 * Description: stop setup timer (checks if running)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "StopTimer"
PetscErrorCode RegOpt::StopTimer(TimerType id)
{
    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    msg="fatal error: timer has not been started";
    ierr=Assert(this->m_TimerIsRunning[id],msg); CHKERRQ(ierr);

    this->m_TempTimer[id] += MPI_Wtime();
    this->m_Timer[id][LOG] += this->m_TempTimer[id];

    // tell the world that we stop the timer
    this->m_TimerIsRunning[id] = false;

    PetscFunctionReturn(0);
}





/********************************************************************
 * Name: ProcessTimers
 * Description: process the timers
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ProcessTimers"
PetscErrorCode RegOpt::ProcessTimers()
{
    PetscErrorCode ierr;
    int rval,rank,nproc;
    double ival=0.0,xval=0.0;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nproc);

    for(int id = 0; id < NTIMERS; ++id){

        // remember input value
        ival=this->m_Timer[id][LOG];

        // get maximal execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_MIN,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_Timer[id][MIN]=xval;

        // get maximal execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_Timer[id][MAX] = xval;

        // get mean execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_Timer[id][AVG]=xval;

        if (rank == 0){ this->m_Timer[id][AVG] /= static_cast<double>(nproc); }

    }

    for(int i = 0; i < 5; ++i){

        // remember input value
        ival=this->m_FFTTimers[i][LOG];

        // get maximal execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_MIN,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_FFTTimers[i][MIN] = xval;

        // get maximal execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_FFTTimers[i][MAX] = xval;

        // get mean execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_FFTTimers[i][AVG] = xval;

        if (rank == 0){ this->m_FFTTimers[i][AVG] /= static_cast<double>(nproc); }

    }

    for(int i = 0; i < 4; ++i){

        // remember input value
        ival=this->m_InterpTimers[i][LOG];

        // get maximal execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_MIN,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_InterpTimers[i][MIN] = xval;

        // get maximal execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_InterpTimers[i][MAX] = xval;

        // get mean execution time
        rval=MPI_Reduce(&ival,&xval,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi reduce returned error"); CHKERRQ(ierr);
        this->m_InterpTimers[i][AVG] = xval;

        if (rank == 0){ this->m_InterpTimers[i][AVG] /= static_cast<double>(nproc); }

    }


    PetscFunctionReturn(0);

}




/********************************************************************
 * Name: ResetCounter
 * Description: resets counters
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResetCounters"
PetscErrorCode RegOpt::ResetCounters()
{
    PetscFunctionBegin;

    for(int i = 0; i < NCOUNTERS; ++i){
        this->m_Counter[i] = 0;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ResetCounter
 * Description: resets counter
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResetCounter"
PetscErrorCode RegOpt::ResetCounter(CounterType id)
{
    PetscFunctionBegin;

    this->m_Counter[id] = 0;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: WriteLogFile
 * Description: write log results to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteLogFile"
PetscErrorCode RegOpt::WriteLogFile()
{
    PetscErrorCode ierr;
    std::string filename,fn,line;
    std::ofstream logwriter;
    std::stringstream ss, ssnum;
    int rank, nnum, nstr, nproc;

    PetscFunctionBegin;

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nproc);


    line = std::string(this->m_LineLength,'-');

    // write out logfile
    if (rank == 0){

        nnum = 20; nstr = 20;
        filename = this->m_XFolder + "registration-performance";
        fn = filename + ".log";


        // create output file
        logwriter.open(fn.c_str());
        ierr=Assert(logwriter.is_open(),"could not open file for writing"); CHKERRQ(ierr);

        logwriter << std::scientific;
        logwriter << line <<std::endl;
        logwriter << "# problem setup" <<std::endl;
        logwriter << line <<std::endl;
        ss  << this->m_MiscOpt->N[0] << " x "
            << this->m_MiscOpt->N[1] << " x "
            << this->m_MiscOpt->N[2];

        logwriter << std::left
                  << std::setw(nstr) << " nx" << std::right
                  << std::setw(nnum) << ss.str() << std::endl;
        ss.clear(); ss.str(std::string());

        logwriter << std::left
                  << std::setw(nstr) << " nt" << std::right
                  << std::setw(nnum) << this->m_NT << std::endl;

        logwriter << std::left
                  << std::setw(nstr) << " n" << std::right
                  << std::setw(nnum) << this->GetNGlobal() << std::endl;

        logwriter << std::left
                  << std::setw(nstr) << " nl" << std::right
                  << std::setw(nnum) << this->GetNLocal() << std::endl;

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
        logwriter << line <<std::endl;
        logwriter << "# timers" <<std::endl;
        logwriter << line <<std::endl;

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
        if (this->m_Timer[T2SEXEC][LOG] > 0.0){
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
        if (this->m_Timer[PDEEXEC][LOG] > 0.0){
            ierr=Assert(this->m_Counter[PDESOLVE] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " pde solves" << std::right
                << std::setw(nnum) << this->m_Timer[PDEEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[PDEEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[PDEEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[PDEEXEC][MAX]/static_cast<ScalarType>(this->m_Counter[PDESOLVE]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[OBJEXEC][LOG] > 0.0){
            ierr=Assert(this->m_Counter[OBJEVAL] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " obj eval" << std::right
                << std::setw(nnum) << this->m_Timer[OBJEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[OBJEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[OBJEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[PDEEXEC][MAX]/static_cast<ScalarType>(this->m_Counter[OBJEVAL]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[GRADEXEC][LOG] > 0.0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " grad eval" << std::right
                << std::setw(nnum) << this->m_Timer[GRADEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[GRADEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[GRADEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[GRADEXEC][MAX]/static_cast<ScalarType>(this->m_Counter[GRADEVAL]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[HMVEXEC][LOG] > 0.0){
            ierr=Assert(this->m_Counter[HESSMATVEC] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " hess mat vec" << std::right
                << std::setw(nnum) << this->m_Timer[HMVEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[HMVEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[HMVEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[HMVEXEC][MAX]/static_cast<ScalarType>(this->m_Counter[HESSMATVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_Timer[PMVEXEC][LOG] > 0.0){
            ierr=Assert(this->m_Counter[PCMATVEC] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " pc mat vec" << std::right
                << std::setw(nnum) << this->m_Timer[PMVEXEC][MIN]
                << std::setw(nnum) << this->m_Timer[PMVEXEC][MAX]
                << std::setw(nnum) << this->m_Timer[PMVEXEC][AVG]
                << std::setw(nnum) << this->m_Timer[PMVEXEC][MAX]/static_cast<ScalarType>(this->m_Counter[PCMATVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }


        // if time has been logged
        if (this->m_Timer[FFTSETUP][LOG] > 0.0){
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
        if (this->m_FFTTimers[2][LOG] > 0.0){
            ierr=Assert(this->m_Counter[FFT] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " FFT communication" << std::right
                << std::setw(nnum) << this->m_FFTTimers[2][MIN]
                << std::setw(nnum) << this->m_FFTTimers[2][MAX]
                << std::setw(nnum) << this->m_FFTTimers[2][AVG]
                << std::setw(nnum) << this->m_FFTTimers[2][MAX]/static_cast<ScalarType>(this->m_Counter[FFT]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_FFTTimers[4][LOG] > 0.0){
            ierr=Assert(this->m_Counter[FFT] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " FFT execution" << std::right
                << std::setw(nnum) << this->m_FFTTimers[4][MIN]
                << std::setw(nnum) << this->m_FFTTimers[4][MAX]
                << std::setw(nnum) << this->m_FFTTimers[4][AVG]
                << std::setw(nnum) << this->m_FFTTimers[4][MAX]/static_cast<ScalarType>(this->m_Counter[FFT]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_InterpTimers[0][LOG] > 0.0){
            ierr=Assert(this->m_Counter[IP] > 0 && this->m_Counter[IPVEC] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp comm" << std::right
                << std::setw(nnum) << this->m_InterpTimers[0][MIN]
                << std::setw(nnum) << this->m_InterpTimers[0][MAX]
                << std::setw(nnum) << this->m_InterpTimers[0][AVG]
                << std::setw(nnum) << this->m_InterpTimers[0][MAX]/static_cast<ScalarType>(this->m_Counter[IP] + this->m_Counter[IPVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_InterpTimers[1][LOG] > 0.0){
            ierr=Assert(this->m_Counter[IP] > 0 && this->m_Counter[IPVEC] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp exec" << std::right
                << std::setw(nnum) << this->m_InterpTimers[1][MIN]
                << std::setw(nnum) << this->m_InterpTimers[1][MAX]
                << std::setw(nnum) << this->m_InterpTimers[1][AVG]
                << std::setw(nnum) << this->m_InterpTimers[1][MAX]/static_cast<ScalarType>(this->m_Counter[IP] + this->m_Counter[IPVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_InterpTimers[2][LOG] > 0.0){
            ierr=Assert(this->m_Counter[IP] > 0 && this->m_Counter[IPVEC] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp alloc" << std::right
                << std::setw(nnum) << this->m_InterpTimers[2][MIN]
                << std::setw(nnum) << this->m_InterpTimers[2][MAX]
                << std::setw(nnum) << this->m_InterpTimers[2][AVG]
                << std::setw(nnum) << this->m_InterpTimers[2][MAX]/static_cast<ScalarType>(this->m_Counter[IP] + this->m_Counter[IPVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        // if time has been logged
        if (this->m_InterpTimers[3][LOG] > 0.0){
            ierr=Assert(this->m_Counter[IP] > 0 && this->m_Counter[IPVEC] > 0,"bug in counter"); CHKERRQ(ierr);
            ss  << std::scientific << std::left
                << std::setw(nstr) << " interp sort" << std::right
                << std::setw(nnum) << this->m_InterpTimers[3][MIN]
                << std::setw(nnum) << this->m_InterpTimers[3][MAX]
                << std::setw(nnum) << this->m_InterpTimers[3][AVG]
                << std::setw(nnum) << this->m_InterpTimers[3][MAX]/static_cast<ScalarType>(this->m_Counter[IP] + this->m_Counter[IPVEC]);
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        logwriter << std::endl;
        logwriter << line <<std::endl;
        logwriter << "# counters" <<std::endl;
        logwriter << line <<std::endl;

        if (this->m_Counter[OBJEVAL] > 0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " objective evals" << std::right
                << std::setw(nnum) << this->m_Counter[OBJEVAL];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[GRADEVAL] > 0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " grad evals" << std::right
                << std::setw(nnum) << this->m_Counter[GRADEVAL];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[PDESOLVE] > 0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " PDE solves" << std::right
                << std::setw(nnum) << this->m_Counter[PDESOLVE];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[HESSMATVEC] > 0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " hess mat vecs" << std::right
                << std::setw(nnum) << this->m_Counter[HESSMATVEC];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[PCMATVEC] > 0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " pc mat vecs" << std::right
                << std::setw(nnum) << this->m_Counter[PCMATVEC];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

        }

        if (this->m_Counter[IP] > 0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " ips sca" << std::right
                << std::setw(nnum) << this->m_Counter[IP];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[IPVEC] > 0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " ips vec" << std::right
                << std::setw(nnum) << this->m_Counter[IPVEC];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }

        if (this->m_Counter[FFT] > 0){
            ss  << std::scientific << std::left
                << std::setw(nstr) << " ffts" << std::right
                << std::setw(nnum) << this->m_Counter[FFT];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());
        }


    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: DisplayTimeToSolution
 * Description: displays the global exection time
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DisplayTimeToSolution"
PetscErrorCode RegOpt::DisplayTimeToSolution()
{

    PetscErrorCode ierr;
    double hours,minutes,seconds,millisec,time;
    int rank;
    std::stringstream ss;
    std::string line;
    PetscFunctionBegin;

    time = this->m_Timer[T2SEXEC][MAX];

    hours    = time / 3600.0;
    minutes  = (hours   - floor(hours)  ) *   60.0;
    seconds  = (minutes - floor(minutes)) *   60.0;
    millisec = (seconds - floor(seconds)) * 1000.0;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    line=std::string(80,'-');

    ierr=PetscPrintf(PETSC_COMM_WORLD,"\n%s\n",line.c_str()); CHKERRQ(ierr);

    ss  << "computation finished ( elapsed cpu time "
        << floor(hours) << " h "
        << floor(minutes) << " m "
        << floor(seconds) << " s "
        << floor(millisec) << " ms )";

    ierr=Msg(ss.str()); CHKERRQ(ierr);
    ierr=PetscPrintf(PETSC_COMM_WORLD,"%s\n",line.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



} // end of namespace

#endif //_REGOPT_CPP_
