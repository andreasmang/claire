
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
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegOpt"
RegOpt::~RegOpt()
{

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

    this->m_NT = 16;
    this->m_Beta[0] = 1E-2;
    this->m_Beta[1] = 1E-2;
    this->m_RegNorm = H1SN;
    this->m_Verbosity = 1;
    this->m_TimeHorizon[0] = 0.0;
    this->m_TimeHorizon[1] = 1.0;
    this->m_OptMeth = GAUSSNEWTON;
    this->m_PDESolver = SL;
    this->m_PrecondMeth = INVREG;
    //this->m_PrecondMeth = TWOLEVEL;
    this->m_FSeqType = QDFS;
    //this->m_PCSolverType = PCPCG;
    this->m_XFolder = "./results/";
    this->m_IFolder = "./results/";
    this->m_ReferenceFN = "./";
    this->m_TemplateFN = "./";
    this->m_NetworkDims[0] = 1;
    this->m_NetworkDims[1] = 1;
    this->m_LineLength = 100;
    this->m_Sigma = 1;
    this->m_WriteImages2File = false;

    this->m_KKTSolverPara.tol[0] = 1E-12; // relative tolerance
    this->m_KKTSolverPara.tol[1] = 1E-12; // absolute tolerance
    this->m_KKTSolverPara.tol[2] = 1E+06; // divergence tolerance
    this->m_KKTSolverPara.maxit  = 1E+03; // maximal iterations

    this->m_OptPara.tol[0] = 1E-9;  // grad abs tol
    this->m_OptPara.tol[1] = 1E-16; // grad rel tol
    this->m_OptPara.tol[2] = 1E-6;  // grad rel tol
    this->m_OptPara.maxit  = 1E3;

    this->m_DD.n = 2;

    this->m_MonitorCFLCondition = false;
    this->m_ReadImagesFromFile = false;
    this->m_InCompressible = true;

    ierr=this->ResetTimers(); CHKERRQ(ierr);
    ierr=this->ResetCounters(); CHKERRQ(ierr);

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
                    << this->m_NetworkDims[0] <<"x"
                    << this->m_NetworkDims[1]<<std::endl;
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
        switch(this->m_RegNorm){
            case L2: std::cout<<"L2-norm (beta="<<this->m_Beta[0]<<")"<<std::endl; break;
            case H1: std::cout<<"H1-norm (beta="<<this->m_Beta[0]<< ", "<<this->m_Beta[1]<<")"<<std::endl; break;
            case H2: std::cout<<"H2-norm (beta="<<this->m_Beta[0]<< ", "<<this->m_Beta[1]<<")"<<std::endl; break;
            case H1SN: std::cout<<"H1-seminorm (beta="<<this->m_Beta[0]<<")"<<std::endl; break;
            case H2SN: std::cout<<"H2-seminorm (beta="<<this->m_Beta[0]<<")"<<std::endl; break;
            default:
            {
                ierr=ThrowError("regularization model not implemented"); CHKERRQ(ierr);
                break;
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

            switch(this->m_FSeqType){
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
        for (int j = 0; j < NVALTYPES; ++j){
            this->m_Timer[i][j] = 0.0;
        }
        this->m_TimerIsRunning[i] = false;
        this->m_TempTimer[i] = 0.0;
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

    ierr=this->ProcessTimers(); CHKERRQ(ierr);

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

        ss << this->m_NetworkDims[0] << " x " << this->m_NetworkDims[1];
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
