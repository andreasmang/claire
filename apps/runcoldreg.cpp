/**
 *  Description: main test function for registration functionality
 *  Copyright (c) 2015-2016.
 *  All rights reserved.
 *  This file is part of PGLISTR library.
 *
 *  PGLISTR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  PGLISTR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PGLISTR.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "RegOpt.h"
#include "RegUtils.h"
#include "PreProcessingRegistration.h"
#include "SynProbRegistration.h"
#include "DataReadWriteRegistration.h"
#include "LargeDeformationRegistration.h"
#include "OptimalControlRegistration.h"
#include "OptimalControlRegistrationIC.h"
#include "TaoInterfaceRegistration.h"


static char help[] =
"\n\
 ====================================================================================================\n\
 optimal control based diffeomorphic image registration\n\
 ====================================================================================================\n\
 usage: runreg <options>\n\
 where the options are:\n\
 ----------------------------------------------------------------------------------------------------\n\
  -nx <int>              number of spatial grid points (x1,x2 and x3 direction; default: 32)\n\
  -nx1 <int>             number of spatial grid points (x1 direction; default: 32)\n\
  -nx2 <int>             number of spatial grid points (x2 direction; default: 32)\n\
  -nx3 <int>             number of spatial grid points (x3 direction; default: 32)\n\
  -nt <int>              number of time points (default: nx)\n\
 ----------------------------------------------------------------------------------------------------\n\
  -x <string>            output folder (default: ./results/)\n\
  -mr <string>           input reference image (default: no image)\n\
  -mt <string>           input tempalte image (default: no image)\n\
 ----------------------------------------------------------------------------------------------------\n\
  -nprocx1 <int>         number of procs in x1 direction (default: 1)\n\
  -nprocx2 <int>         number of procs in x2 direction (default: 1)\n\
  -nthreads <int>        number of threads\n\
 ----------------------------------------------------------------------------------------------------\n\
  -beta <double>         regularization paramameter (default: 1E-2)\n\
 ----------------------------------------------------------------------------------------------------\n\
  -pdesolver <type>      switch between PDE solvers\n\
                         where type is one of the following\n\
                         rk2\n\
                         sl\n\
  -xtimehist             flag to store time series of images\n\
  -ximages               flag to store time series of images\n\
 ----------------------------------------------------------------------------------------------------\n\
  -maxit <int>           maximal number of (outer) iterations (default: 1000)\n\
  -kktmaxit <int>        maximal number of krylov (inner) iterations (default: 1000)\n\
  -gttol <dbl>           reduction of gradient ||g[v]||/||g[v0]|| <= tol; default: 1E-4\n\
  -jabstol <dbl>         absolute convergence tolerance: (J[v] - J[v*]) <= tol; default: 1E-16\n\
  -jreltol <dbl>         relative convergence tolerance: (J[v] - J[v*])/J[v*] <= tol; default: 1E-12\n\
  -gabstol <dbl>         absolute convergence tolerance ||g[v]|| <= tol; default: 1E-12\n\
  -greltol <dbl>         relative convergence tolerance ||g[v]||/|J[v]| <= tol; default: 1E-12\n\
 ----------------------------------------------------------------------------------------------------\n\
  -verbosity <int>       control verbosity level (default: 0)\n\
 ====================================================================================================\n\
\n\
\n\
\n\
";



int ParseArguments(reg::RegOpt**,N_MISC**,accfft_plan**,double&,MPI_Comm&);



/********************************************************************
 * Name: main
 * Description: main function to run registration
 *******************************************************************/
int main(int argc,char **argv)
{
    PetscErrorCode ierr;
    int procid,nprocs;
    IntType nl,ng;
    accfft_plan* plan;
    double setuptime;
    Vec mT = NULL, mR=NULL;
    MPI_Comm c_comm;
    PetscBool flag=PETSC_FALSE;

    N_MISC *miscopt = NULL;
    reg::RegOpt* opt = NULL;

    typedef reg::DataReadWriteRegistration ReadWriteType;
    typedef reg::SynProbRegistration SynProbType;
    typedef reg::LargeDeformationRegistration LDRegType;
    typedef reg::OptimalControlRegistration OCRegType;
    typedef reg::OptimalControlRegistrationIC OCRegTypeIC;

    ReadWriteType* io = NULL;
    LDRegType* registration = NULL;
    SynProbType* synprob = NULL;

    // init petsc
    PetscInitialize(&argc,&argv,NULL,help);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD,&procid);

    // init gpu
    //ScalarType *dummy_d;
    //cudaMallocHost((void**)&dummy_d, sizeof(ScalarType));
    //cudaFreeHost(dummy_d);

    ierr=PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"optimal control registration","");
    ierr=PetscOptionsGetBool(NULL,NULL,"-help",&flag,NULL); CHKERRQ(ierr);
    if (flag==PETSC_FALSE){ ierr=PetscOptionsGetBool(NULL,NULL,"-help",&flag,NULL); CHKERRQ(ierr); }
    if (flag==PETSC_FALSE){ ierr==PetscOptionsGetBool(NULL,NULL,"-HELP",&flag,NULL); CHKERRQ(ierr); }
    if (flag==PETSC_FALSE){ ierr=PetscOptionsGetBool(NULL,NULL,"-h",&flag,NULL); CHKERRQ(ierr); }
    if (flag==PETSC_TRUE){ ierr=PetscFinalize(); CHKERRQ(ierr); return 0; }
    ierr=PetscOptionsEnd(); CHKERRQ(ierr);

    // set the user parameters
    ParseArguments(&opt,&miscopt,&plan,setuptime,c_comm);

    // allocate class for io
    try{ io = new ReadWriteType(opt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // allocate class for registration
    if (opt->InCompressible()){
        try{ registration = new OCRegTypeIC(opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    else{
        try{ registration = new OCRegType(opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr=registration->SetIO(io); CHKERRQ(ierr);

    ierr=reg::Msg("running registration"); CHKERRQ(ierr);

    nl = opt->GetNLocal();
    ng = opt->GetNGlobal();


    if(opt->ReadImagesFromFile()){

        ierr=VecCreate(PETSC_COMM_WORLD,&mR); CHKERRQ(ierr);
        ierr=VecSetSizes(mR,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mR); CHKERRQ(ierr);
        ierr=io->Read(mR,opt->GetReferenceFN()); CHKERRQ(ierr);
        ierr=reg::Assert(mR!=NULL, "input reference image is null pointer"); CHKERRQ(ierr);

        ierr=VecCreate(PETSC_COMM_WORLD,&mT); CHKERRQ(ierr);
        ierr=VecSetSizes(mT,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mT); CHKERRQ(ierr);
        ierr=io->Read(mT,opt->GetTemplateFN()); CHKERRQ(ierr);
        ierr=reg::Assert(mT!=NULL, "input template image is null pointer"); CHKERRQ(ierr);

        // pass to registration
        ierr=registration->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr=registration->SetTemplateImage(mT); CHKERRQ(ierr);

    }
    else{

        ierr=VecCreate(PETSC_COMM_WORLD,&mT); CHKERRQ(ierr);
        ierr=VecSetSizes(mT,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mT); CHKERRQ(ierr);
        ierr=VecSet(mT,0.0); CHKERRQ(ierr);

        // allocate class for registration
        try{ synprob = new SynProbType(opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        // set up a synthetic test problem/compute template image
        ierr=synprob->ComputeSmoothScalarField(mT,0); CHKERRQ(ierr);

        // advect template image to obtain reference image;
        // the variables for the images are set internally; so no need
        // to get the output here
        ierr=registration->SetupSyntheticProb(mT); CHKERRQ(ierr);
    }

    Vec v0;
    IntType nlukwn,ngukwn;
    nlukwn = 3*miscopt->N_local;
    ngukwn = 3*miscopt->N_global;

    // allocate vector fields
    ierr=VecCreate(PETSC_COMM_WORLD,&v0); CHKERRQ(ierr);
    ierr=VecSetSizes(v0,nlukwn,ngukwn); CHKERRQ(ierr);
    ierr=VecSetFromOptions(v0); CHKERRQ(ierr);
    ierr=VecSet(v0,0.0); CHKERRQ(ierr);

    // create tao solver
    Tao tao=NULL;
    std::string method = "nls";
    ierr=TaoCreate(PETSC_COMM_WORLD,&tao); CHKERRQ(ierr);
    ierr=TaoSetType(tao,"nls"); CHKERRQ(ierr);
    ierr=TaoSetInitialVector(tao,v0); CHKERRQ(ierr);

    // set the routine to evaluate the objective and compute the gradient
    ierr=TaoSetObjectiveRoutine(tao,reg::EvaluateObjective,(void*)registration); CHKERRQ(ierr);
    ierr=TaoSetGradientRoutine(tao,reg::EvaluateGradient,(void*)registration); CHKERRQ(ierr);
    ierr=TaoSetObjectiveAndGradientRoutine(tao,reg::EvaluateObjectiveGradient,(void*)registration); CHKERRQ(ierr);

    // set the monitor for the optimization process
    ierr=TaoCancelMonitors(tao); CHKERRQ(ierr);
    ierr=TaoSetMonitor(tao,reg::OptimizationMonitor,registration,NULL); CHKERRQ(ierr);

    TaoLineSearch ls;
    ierr=TaoGetLineSearch(tao,&ls); CHKERRQ(ierr);
    ierr=TaoLineSearchSetType(ls,"armijo"); CHKERRQ(ierr);

    // set tolearances for optimizer
    ScalarType gatol,grtol,gttol;
    gatol = opt->GetOptTol(0);
    grtol = opt->GetOptTol(1);
    gttol = opt->GetOptTol(2);
    ierr=TaoSetTolerances(tao,gatol,grtol,gttol); CHKERRQ(ierr);

    ierr=TaoSetMaximumIterations(tao,opt->GetOptMaxit() - 1); CHKERRQ(ierr);

    Mat HMatVec;
    ierr=MatCreateShell(PETSC_COMM_WORLD,nlukwn,nlukwn,ngukwn,ngukwn,static_cast<void*>(registration),&HMatVec); CHKERRQ(ierr);
    ierr=MatShellSetOperation(HMatVec,MATOP_MULT,(void(*)(void))reg::HessianMatVec); CHKERRQ(ierr);
    ierr=MatSetOption(HMatVec,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);
    ierr=TaoSetHessianRoutine(tao,HMatVec,HMatVec,reg::EvaluateHessian,static_cast<void*>(&registration)); CHKERRQ(ierr);

    KSP taoksp;
    PC taokktpc;

    // get the ksp of the optimizer and set options
    ierr=TaoGetKSP(tao,&taoksp); CHKERRQ(ierr);

    // ksp is only nonzero if we use a newton type method
    if (taoksp != NULL){

        if(opt->GetVerbosity() >= 2){
            ierr=KSPMonitorSet(taoksp,reg::KrylovMonitor,registration,NULL); CHKERRQ(ierr);
        }

        // set the preconditioner
        ierr=KSPGetPC(taoksp,&taokktpc); CHKERRQ(ierr);

        if(opt->GetPrecondMeth() == reg::NOPC){
            if (strcmp(method.c_str(),"nls") == 0){
                ierr=PetscOptionsSetValue(NULL,"-tao_nls_pc_type","none"); CHKERRQ(ierr);
                ierr=PetscOptionsSetValue(NULL,"-tao_nls_ksp_type","cg"); CHKERRQ(ierr);
                ierr=TaoSetFromOptions(tao); CHKERRQ(ierr);
            }
            else if (strcmp(method.c_str(),"ntr") == 0){
                ierr=PetscOptionsSetValue(NULL,"-tao_ntr_pc_type","none"); CHKERRQ(ierr);
                ierr=PetscOptionsSetValue(NULL,"-tao_ntr_ksp_type","stcg"); CHKERRQ(ierr);
                ierr=TaoSetFromOptions(tao); CHKERRQ(ierr);
            }
            ierr=PCSetType(taokktpc,PCNONE); CHKERRQ(ierr);
        }
        else if (  (opt->GetPrecondMeth() == reg::INVREG)
                || (opt->GetPrecondMeth() == reg::TWOLEVEL) ) {
            if (strcmp(method.c_str(),"nls") == 0){
                ierr=PetscOptionsSetValue(NULL,"-tao_nls_pc_type","petsc"); CHKERRQ(ierr);
                ierr=PetscOptionsSetValue(NULL,"-tao_nls_ksp_type","cg"); CHKERRQ(ierr);
                ierr=TaoSetFromOptions(tao); CHKERRQ(ierr);
            }
            else if (strcmp(method.c_str(),"ntr") == 0){
                ierr=PetscOptionsSetValue(NULL,"-tao_ntr_pc_type","petsc"); CHKERRQ(ierr);
                ierr=PetscOptionsSetValue(NULL,"-tao_ntr_ksp_type","stcg"); CHKERRQ(ierr);
                ierr=TaoSetFromOptions(tao); CHKERRQ(ierr);
            }
            ierr=PCSetType(taokktpc,PCSHELL); CHKERRQ(ierr);
            ierr=PCShellSetApply(taokktpc,reg::PrecondMatVec); CHKERRQ(ierr);
            ierr=PCShellSetContext(taokktpc,registration); CHKERRQ(ierr);
            //ierr=PCShellSetName(taokktpc,"kktpc"); CHKERRQ(ierr);
            ierr=PCShellSetSetUp(taokktpc,reg::PrecondSetup); CHKERRQ(ierr);
        }
        else{
            ierr=reg::ThrowError("preconditioner not defined"); CHKERRQ(ierr);
        }

        // set tolerances for krylov subspace method
        ScalarType reltol,abstol,divtol;
        IntType maxit;
        reltol = opt->GetKKTSolverTol(0); // 1E-12;
        abstol = opt->GetKKTSolverTol(1); // 1E-12;
        divtol = opt->GetKKTSolverTol(2); // 1E+06;
        maxit  = opt->GetKKTMaxit(); // 1000;
        ierr=KSPSetTolerances(taoksp,reltol,abstol,divtol,maxit); CHKERRQ(ierr);
        //ierr=KSPSetInitialGuessNonzero(taoksp,PETSC_FALSE); CHKERRQ(ierr);
        ierr=KSPSetInitialGuessNonzero(taoksp,PETSC_TRUE); CHKERRQ(ierr);

        //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
        //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
        //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
        ierr=KSPSetNormType(taoksp,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    }

    // reset all the clocks we have used so far
    ierr=opt->ResetTimers(); CHKERRQ(ierr);
    ierr=opt->ResetCounters(); CHKERRQ(ierr);

    // set the fft setup time
    opt->SetFFTSetupTime(setuptime);

    // do the inversion
    ierr=opt->StartTimer(reg::T2SEXEC); CHKERRQ(ierr);
    ierr=reg::Msg("starting optimization"); CHKERRQ(ierr);
    ierr=TaoSolve(tao); CHKERRQ(ierr);
    ierr=reg::Msg("optimization done"); CHKERRQ(ierr);
    ierr=opt->StopTimer(reg::T2SEXEC); CHKERRQ(ierr);

    // get the solution vector
    Vec vstar;
    ierr=TaoGetSolutionVector(tao,&vstar); CHKERRQ(ierr);

    // finalize the registration (write out all data)
    ierr=registration->Finalize(vstar); CHKERRQ(ierr);

    // display info to user, once we're done
    ierr=TaoView(tao,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    ierr=opt->DisplayTimeToSolution(); CHKERRQ(ierr);

    // clean up
    if (io != NULL){ delete io; io = NULL; }
    if (opt != NULL){ delete opt; opt = NULL; }
    if (miscopt != NULL){ delete miscopt; miscopt = NULL; }
    if (synprob != NULL){ delete synprob; synprob = NULL; }
    if (tao!=NULL){ ierr=TaoDestroy(&tao); CHKERRQ(ierr); tao=NULL; }

    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); mT=NULL; }
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); mR=NULL; }

    if (registration != NULL){ delete registration; registration = NULL; }

    accfft_destroy_plan(plan);
    accfft_cleanup();
    MPI_Comm_free(&c_comm);
    ierr=PetscFinalize(); CHKERRQ(ierr);

   return 0;
}



/********************************************************************
 * Name: ParseArguments
 * Description: parse arguments
 *******************************************************************/
int ParseArguments(reg::RegOpt** regopt, N_MISC** miscopt,
                    accfft_plan** plan,double &fftsetuptime,MPI_Comm& c_comm)
{
    PetscErrorCode ierr;
    int nx[3],isize[3],istart[3],osize[3],ostart[3],c_dims[2];
    IntType nthreads,ivalue,maxit,nt,kktmaxit;
    IntType alloc_max;
    int nprocs,rank,verbosity;
    ScalarType beta,opttol[5];
    PetscBool flag;
    bool storetimeseries=false,storeimages=false,readmT=false,readmR=false,isincompressible=false;
    reg::PDESolver pdesoltype;
    reg::RegNorm regnorm;
    reg::OptMeth optmeth;
    std::string msg,xfolder,mTFN,mRFN;
    char istring[PETSC_MAX_PATH_LEN];
    Complex* uk = NULL;
    std::stringstream ss;
    ScalarType* u = NULL;

    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    /////////////////////////////////////////////////////////////////////////////////
    PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"optimal control registration options","");

    // number of grid points
    ivalue=32;
    ierr=PetscOptionsGetInt(NULL,NULL,"-nx",&ivalue,NULL); CHKERRQ(ierr);
    ierr=PetscOptionsGetInt(NULL,NULL,"-nx1",&ivalue,NULL); CHKERRQ(ierr); nx[0] = ivalue;
    ierr=PetscOptionsGetInt(NULL,NULL,"-nx2",&ivalue,NULL); CHKERRQ(ierr); nx[1] = ivalue;
    ierr=PetscOptionsGetInt(NULL,NULL,"-nx3",&ivalue,NULL); CHKERRQ(ierr); nx[2] = ivalue;

    // number of threads
    nthreads=1;
    ierr=PetscOptionsGetInt(NULL,NULL,"-nthreads",&nthreads,&flag); CHKERRQ(ierr);

    // cartesian grid for mpi
    ivalue=1;
    ierr=PetscOptionsGetInt(NULL,NULL,"-nprocx1",&ivalue,&flag); CHKERRQ(ierr);
    c_dims[0] = ivalue;

    ivalue=1;
    ierr=PetscOptionsGetInt(NULL,NULL,"-nprocx2",&ivalue,&flag); CHKERRQ(ierr);
    c_dims[1] = ivalue;

    ivalue=4;
    ierr=PetscOptionsGetInt(NULL,NULL,"-nt",&ivalue,&flag); CHKERRQ(ierr);
    nt = ivalue;

    beta=1E-2;
    ierr=PetscOptionsGetReal(NULL,NULL,"-beta",&beta,&flag); CHKERRQ(ierr);

    strcpy(istring,"");
    ierr=PetscOptionsGetString(NULL,NULL,"-mr",istring,sizeof(istring),&flag); CHKERRQ(ierr);
    if (flag == PETSC_TRUE){
        mRFN = std::string(istring);
        readmR = true;
    }

    strcpy(istring,"");
    ierr=PetscOptionsGetString(NULL,NULL,"-mt",istring,sizeof(istring),&flag); CHKERRQ(ierr);
    if (flag == PETSC_TRUE){
        mTFN = std::string(istring);
        readmT = true;
    }

    strcpy(istring,"./results/");
    ierr=PetscOptionsGetString(NULL,NULL,"-x",istring,sizeof(istring),&flag); CHKERRQ(ierr);
    xfolder = std::string(istring);

    strcpy(istring,"rk2");
    ierr=PetscOptionsGetString(NULL,NULL,"-pdesolver",istring,sizeof(istring),&flag); CHKERRQ(ierr);
    if      (strcmp(istring,"rk2")==0){ pdesoltype=reg::RK2; }
    else if (strcmp(istring,"sl")==0){ pdesoltype=reg::SL; }
    else{ ierr=reg::ThrowError("pde solver not defined"); CHKERRQ(ierr); }

    strcpy(istring,"h2s");
    ierr=PetscOptionsGetString(NULL,NULL,"-regnorm",istring,sizeof(istring),&flag); CHKERRQ(ierr);
    if      (strcmp(istring,"l2")==0){ regnorm=reg::L2; }
    else if (strcmp(istring,"h1")==0){ regnorm=reg::H1; }
    else if (strcmp(istring,"h2")==0){ regnorm=reg::H2; }
    //else if (strcmp(istring,"h3")==0){ regnorm=reg::H3; }
    else if (strcmp(istring,"h1s")==0){ regnorm=reg::H1SN; }
    else if (strcmp(istring,"h2s")==0){ regnorm=reg::H2SN; }
    //else if (strcmp(istring,"h3s")==0){ regnorm=reg::H3SN; }
    else{ ierr=reg::ThrowError("regularization norm not defined"); CHKERRQ(ierr); }

    strcpy(istring,"gn");
    ierr=PetscOptionsGetString(NULL,NULL,"-optmeth",istring,sizeof(istring),&flag); CHKERRQ(ierr);
    if      (strcmp(istring,"gn")==0){ optmeth=reg::GAUSSNEWTON; }
    else if (strcmp(istring,"fn")==0){ optmeth=reg::FULLNEWTON; }
    else{ ierr=reg::ThrowError("optimization method not defined"); CHKERRQ(ierr); }

    ivalue=1E3;
    ierr=PetscOptionsGetInt(NULL,NULL,"-maxit",&ivalue,&flag); CHKERRQ(ierr);
    maxit=ivalue;

    ivalue=1E3;
    ierr=PetscOptionsGetInt(NULL,NULL,"-kktmaxit",&ivalue,&flag); CHKERRQ(ierr);
    kktmaxit=ivalue;

    opttol[0]=1E-6;
    ierr=PetscOptionsGetReal(NULL,NULL,"-gabstol",&opttol[0],&flag); CHKERRQ(ierr);

    opttol[1]=1E-16;
    ierr=PetscOptionsGetReal(NULL,NULL,"-greltol",&opttol[1],&flag); CHKERRQ(ierr);

    opttol[2]=1E-3;
    ierr=PetscOptionsGetReal(NULL,NULL,"-gttol",&opttol[2],&flag); CHKERRQ(ierr);

    flag = PETSC_FALSE;
    ierr=PetscOptionsGetBool(NULL,NULL,"-xtimehist",&flag,NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) storetimeseries = true;
    else storetimeseries = false;

    flag = PETSC_FALSE;
    ierr=PetscOptionsGetBool(NULL,NULL,"-ic",&flag,NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) isincompressible = true;
    else isincompressible = false;

    flag = PETSC_FALSE;
    ierr=PetscOptionsGetBool(NULL,NULL,"-ximages",&flag,NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) storeimages = true;
    else storeimages = false;


    ivalue=1;
    ierr=PetscOptionsGetInt(NULL,NULL,"-verbosity",&ivalue,&flag); CHKERRQ(ierr);
    verbosity=ivalue;


    PetscOptionsEnd();
    /////////////////////////////////////////////////////////////////////////////////

    // initialize accft
    accfft_create_comm(MPI_COMM_WORLD,c_dims,&c_comm);
    accfft_init(nthreads);
    alloc_max = accfft_local_size_dft_r2c(nx,isize,istart,osize,ostart,c_comm);
    u  = (ScalarType*)accfft_alloc(alloc_max);
    uk = (Complex*)accfft_alloc(alloc_max);

    ierr=reg::Assert(nthreads > 0,"omp threads < 0"); CHKERRQ(ierr);

    omp_set_dynamic(0);
    omp_set_num_threads(nthreads);


    // check if number of threads is consistent with user options
    int ompthreads=omp_get_max_threads();

    ss.str( std::string() );
    ss.clear();
    ss << "max number of openmp threads is not a match (user,set)=("
       << nthreads <<"," << ompthreads <<")\n";
    ierr=reg::Assert(ompthreads == nthreads,ss.str().c_str()); CHKERRQ(ierr);

    // set up the fft
    fftsetuptime=-MPI_Wtime();
    *plan = accfft_plan_dft_3d_r2c(nx,u,(double*)uk,c_comm,ACCFFT_MEASURE);
    fftsetuptime+=MPI_Wtime();

    try{ *miscopt = new N_MISC(nx,isize,istart,*plan,c_comm); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    if (rank == 0) std::cout<<(*miscopt)->N_local<<std::endl;
    if (rank == 0) std::cout<<(*miscopt)->N_global<<std::endl;

    // allocate class for registration options
    try{ *regopt = new reg::RegOpt(*miscopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    if (readmT && readmR){
        (*regopt)->ReadImagesFromFile(true);
        (*regopt)->SetReferenceFN(mRFN);
        (*regopt)->SetTemplateFN(mTFN);
    }
    else if( (readmT == false) && readmR ) {
        ierr=reg::ThrowError("you need to assign two images"); CHKERRQ(ierr);
    }
    else if( readmT && (readmR == false) ) {
        ierr=reg::ThrowError("you need to assign two images"); CHKERRQ(ierr);
    }
    else if( (readmT == false) && (readmR == false) ){
        (*regopt)->ReadImagesFromFile(false);
    }

    // parse options
    (*regopt)->SetNumTimePoints(nt);
    (*regopt)->SetXFolder(xfolder);
    (*regopt)->SetRegNorm(regnorm);
    (*regopt)->SetRegularizationWeight(beta);
    (*regopt)->SetOptMeth(optmeth);
    (*regopt)->SetNetworkDims(c_dims);
    (*regopt)->SetNumThreads(nthreads);
    (*regopt)->SetPDESolver(pdesoltype);
    (*regopt)->SetOptMaxit(maxit);
    (*regopt)->SetKKTMaxit(kktmaxit);
    for(int i=0; i<3; ++i) (*regopt)->SetOptTol(i,opttol[i]);

    (*regopt)->WriteImages2File(storeimages);
    (*regopt)->StoreTimeSeries(storetimeseries);
    (*regopt)->InCompressible(isincompressible);
    (*regopt)->SetVerbosity(verbosity);

    (*regopt)->DisplayOptions();

    accfft_free(u);
    accfft_free(uk);

    return 0;
}
