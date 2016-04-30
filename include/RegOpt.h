/**
 *  Description: class (container) for registration related options
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

#ifndef _REGOPT_H_
#define _REGOPT_H_

//global includes
#include <sstream>


// local includes
#include "RegUtils.h"
#include "utils.hpp"


namespace reg
{




// flags for registration norms
enum PDESolver
{
    RK2,  ///< flag for RK2 solver
    RK2A, ///< flag for stabilized RK2 solver
    SL,   ///< flag for SL
};




// flags for registration norms
enum RegNorm
{
    L2,   ///< flag for L2-norm
    H1,   ///< flag for H1-norm
    H2,   ///< flag for H2-norm
    H3,   ///< flag for H3-norm
    H1SN, ///< flag for H1-seminorm
    H2SN, ///< flag for H2-seminorm
    H3SN, ///< flag for H3-seminorm
};




// flags for optimization methods
enum OptMeth
{
    GAUSSNEWTON, ///< Gauss-Newton approximation
    FULLNEWTON,  ///< full Newton
    GRADDESCENT, ///< gradient descent (gradient in sobolev space)
};




// flags for preconditioner methods
enum PrecondMeth
{
    INVREG,  ///< inverse regularization operator
    TWOLEVEL, ///< 2 level preconditioner
    NOPC,    ///< no preconditioner
};




// flags for optimization
enum FSeqType
{
    QDFS, ///< quadratic forcing sequence
    SLFS, ///< superliner forcing sequence
    NOFS, ///< no forcing sequence
};




// flags for optimization methods
enum PCSolverType
{
    PCPCG,    ///< pcg as preconditioner
    PCFCG,    ///< flexible cg as preconditioner
    PCCHEB,   ///< chebyshev method as preconditioner
    PCGMRES,  ///< gmres solver as preconditioner
    PCFGMRES, ///< flexible gmres as preconditioner
};




/********************************************************************
 * Name: TimerTypes
 * Description: flags for timers
 *******************************************************************/
enum TimerType
{
    T2SEXEC = 0, ///< time to solution (execution time)
    PDEEXEC,     ///< pde solves (execution time)
    HMVEXEC,     ///< hessian mat vec (execution time)
    PMVEXEC,     ///< precond mat vec (execution time)
    GRADEXEC,    ///< gradient evaluation (execution time)
    OBJEXEC,     ///< objective evluation (execution time)
    FFTSETUP,    ///< fft setup time
    FFTEXEC,     ///< fft execution time
    NTIMERS,     ///< to allocate the timers
};




/********************************************************************
 * Name: CounterType
 * Description: flags for timers
 *******************************************************************/
enum CounterType
{
    PDESOLVE = 0, ///< PDE solves
    HESSMATVEC,   ///< # hessian matvecs
    PCMATVEC,     ///< preconditioner matvecs
    OBJEVAL,      ///< objective evaluations
    GRADEVAL,     ///< gradient evaluations
    IPVEC,        ///< interpolation execution time
    IP,           ///< interpolation execution time
    FFT,          ///< fft evaluations
    NCOUNTERS,    ///< to allocate the counters
};




class RegOpt
{
public:
    typedef RegOpt Self;

    RegOpt();
    RegOpt(N_MISC*);
    RegOpt(int,char**);
    ~RegOpt();


    // get functions
    inline IntType GetNLocal(void){return this->m_MiscOpt->N_local;};
    inline IntType GetNGlobal(void){return this->m_MiscOpt->N_global;};
    inline IntType GetNumTimePoints(void){return this->m_NT;};

    inline ScalarType GetSpatialStepSize(int i){return this->m_MiscOpt->h[i];};
    inline ScalarType GetTimeStepSize(void){
        return (this->m_TimeHorizon[1] - this->m_TimeHorizon[0])
               /static_cast<ScalarType>(this->m_NT);
    };

    inline ScalarType GetLebesqueMeasure(void)
    {
        return  this->m_MiscOpt->h[0]
               *this->m_MiscOpt->h[1]
               *this->m_MiscOpt->h[2];
    }

    inline void SetRegularizationWeight(ScalarType beta){
        this->m_Beta[0]=beta; this->m_Beta[1]=beta;
    };
    inline ScalarType GetRegularizationWeight(void){return this->m_Beta[0];};

    inline ScalarType GetSigma(void){return this->m_Sigma;};
    inline ScalarType GetTimeHorizon(int i){return this->m_TimeHorizon[i];};
    inline std::string GetXFolder(void){return this->m_XFolder;};
    inline std::string GetIFolder(void){return this->m_IFolder;};
    inline int GetISize(int i){return this->m_MiscOpt->isize[i];};
    inline int GetIStart(int i){return this->m_MiscOpt->istart[i];};
    inline bool ReadImagesFromFile(){return this->m_ReadImagesFromFile;};
    inline void ReadImagesFromFile(bool flag){this->m_ReadImagesFromFile=flag;};


    inline bool StoreTimeSeries(){return this->m_StoreTimeSeries;};
    inline void StoreTimeSeries(bool flag){this->m_StoreTimeSeries = flag;};
    inline bool MonitorCFLCondition(){return this->m_MonitorCFLCondition;};
    inline void MonitorCFLCondition(bool flag){this->m_MonitorCFLCondition=flag;};

    inline RegNorm GetRegNorm(void){return this->m_RegNorm;};
    inline OptMeth GetOptMeth(void){return this->m_OptMeth;};
    inline PDESolver GetPDESolver(void){return this->m_PDESolver;};
    inline PrecondMeth GetPrecondMeth(void){return this->m_PrecondMeth;};
    inline PCSolverType GetPCSolverType(){return this->m_PCSolverType;};
    inline FSeqType GetFSeqType(void){return this->m_KKTSolverPara.fseqtype;};
    inline void SetVerbosity(int vbs){this->m_Verbosity = vbs;};
    inline int GetVerbosity(){return this->m_Verbosity;};

    inline ScalarType GetJacMin(){return this->m_RegMonitor.jacmin;};
    inline ScalarType GetJacMax(){return this->m_RegMonitor.jacmax;};
    inline ScalarType GetJacMean(){return this->m_RegMonitor.jacmean;};

    inline void SetJacMin(ScalarType value){this->m_RegMonitor.jacmin=value;};
    inline void SetJacMax(ScalarType value){this->m_RegMonitor.jacmax=value;};
    inline void SetJacMean(ScalarType value){this->m_RegMonitor.jacmean=value;};


    inline bool DoParameterContinuation(){return this->m_ParameterCont.enabled;};
    inline void DoParameterContinuation(bool flag){this->m_ParameterCont.enabled = flag;};
    inline ScalarType GetJacBound(){return this->m_ParameterCont.jacbound;};
    inline void SetJacBound(ScalarType value){this->m_ParameterCont.jacbound=value;};


    inline unsigned int GetCounter(CounterType id){return this->m_Counter[id];};
    inline void IncrementCounter(CounterType id){this->m_Counter[id]++;};
    inline void IncrementCounter(CounterType id, unsigned int i){this->m_Counter[id]+= i;};

    inline bool InCompressible(){return this->m_InCompressible;};
    inline void InCompressible(bool flag){this->m_InCompressible=flag;};

    inline void GetTimer(TimerType id,double* wtime){
            wtime[0] = this->m_Timer[id][MIN];
            wtime[1] = this->m_Timer[id][MAX];
            wtime[2] = this->m_Timer[id][AVG];
    };

    inline void SetTemplateFN(std::string s){this->m_TemplateFN = s;};
    inline void SetReferenceFN(std::string s){this->m_ReferenceFN = s;};
    inline std::string GetTemplateFN(){return this->m_TemplateFN;};
    inline std::string GetReferenceFN(){return this->m_ReferenceFN;};

    inline void SetXFolder(std::string s){this->m_XFolder = s;};
    inline void SetIFolder(std::string s){this->m_IFolder = s;};
    inline void SetPDESolver(PDESolver id){this->m_PDESolver=id;};
    inline void SetRegNorm(RegNorm regnorm){this->m_RegNorm=regnorm;};
    inline void SetOptMeth(OptMeth meth){this->m_OptMeth=meth;};

    inline void SetNumThreads(int n){this->m_NumThreads=n;};
    inline void SetNumTimePoints(IntType nt){ this->m_NT=nt; };
    inline int GetNetworkDims(int i){return this->m_CartGridDims[i];};
    inline void SetNetworkDims(int dims[2]){
        this->m_CartGridDims[0]=dims[0];
        this->m_CartGridDims[1]=dims[1];
    };
    inline void IncreaseFFTTimers(double timers[5]){
        for(int i=0; i < 5; ++i) this->m_FFTTimers[i][LOG]+=timers[i];
    };
    inline void IncreaseInterpTimers(double timers[4]){
        for(int i=0; i < 4; ++i) this->m_InterpTimers[i][LOG]+=timers[i];
    };

    inline ScalarType GetOptTol(int i){return this->m_OptPara.tol[i];};
    inline void SetOptTol(int i, ScalarType tol){this->m_OptPara.tol[i] = tol;};
    inline int GetOptMaxit(){return this->m_OptPara.maxit;};
    inline void SetOptMaxit(int i){this->m_OptPara.maxit = i;};

    inline ScalarType GetKKTSolverTol(int i){return this->m_KKTSolverPara.tol[i];};
    inline void SetKKTSolverTol(int i,ScalarType tol){this->m_KKTSolverPara.tol[i]=tol;};
    inline int GetKKTMaxit(){return this->m_KKTSolverPara.maxit;};
    inline void SetKKTMaxit(int i){this->m_KKTSolverPara.maxit = i;};

    inline unsigned int GetDDNumDomains(){return this->m_DD.n;};
    inline unsigned int SetDDNumDomains(unsigned int n){return this->m_DD.n = n;};

    inline bool WriteImages(){return this->m_WriteImages;};
    inline void WriteImages(bool flag){this->m_WriteImages = flag;};
    inline bool LoggingEnabled(){return this->m_WriteLogFiles;};

    ScalarType ComputeFFTScale();

    PetscErrorCode StartTimer(TimerType);
    PetscErrorCode StopTimer(TimerType);
    PetscErrorCode ResetTimers(void);
    PetscErrorCode ResetTimer(TimerType);
    PetscErrorCode ResetCounters(void);
    PetscErrorCode ResetCounter(CounterType);
    PetscErrorCode ProcessTimers();

    PetscErrorCode DisplayOptions();

    PetscErrorCode DisplayTimeToSolution();
    PetscErrorCode WriteLogFile();

    // interface to the global glistr options; are set in constructor
    N_MISC* m_MiscOpt;


private:

    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode ParseArguments(int,char**);
    PetscErrorCode DoSetup(void);
    PetscErrorCode Usage(void);

    enum TimerValue{LOG=0,MIN,MAX,AVG,NVALTYPES};

    IntType m_nx[3];
    IntType m_NT; ///< number of time points
    ScalarType m_TimeHorizon[2];

    // parameters for optimization
    struct Optimization{
        ScalarType tol[3];
        int maxit;
    };

    // parameters for KKT solver
    struct KKTSolver{
        ScalarType tol[3];
        int maxit;
        FSeqType fseqtype; ///<forcing sequence type
    };

    // parameters for KKT solver
    struct DomainDecomposition{
        unsigned int n;
    };

    // parameters for parameter continuation
    struct ParameterContinuation{
        ScalarType betamin; ///< minimal regularization parameter
        ScalarType jacbound; ///< lower bound for jacobian
        bool enabled; ///< flag if parameter continuation is switched on
    };

    // parameters for parameter continuation
    struct RegMonitor{
        ScalarType jacmin; ///< min value of jacobian
        ScalarType jacmax; ///< max value of jacobian
        ScalarType jacmean; ///< mean value of jacobian
    };


    // parameters for KKT solver
    struct MultiScale{
        unsigned int n;
    };

    accfft_plan* m_FFTPlan;
    MPI_Comm m_Comm;
    int m_CartGridDims[2];

    Optimization m_OptPara;
    KKTSolver m_KKTSolverPara;
    ParameterContinuation m_ParameterCont;

    DomainDecomposition m_DD;

    RegMonitor m_RegMonitor;

    RegNorm m_RegNorm; ///< flag for regularization norm
    OptMeth m_OptMeth; ///< flag for optimization method
    PrecondMeth m_PrecondMeth; ///< flag for preconditioner
    PDESolver m_PDESolver; ///< flag for PDE solver
    PCSolverType m_PCSolverType; ///< flag for KSP solver for precond

    std::string m_XFolder; ///< identifier for folder to write results to
    std::string m_IFolder; ///< identifier for folder to read data from
    std::string m_TemplateFN; ///< template image file name
    std::string m_ReferenceFN; ///< reference image file name

    double m_Timer[NTIMERS][NVALTYPES];
    double m_TempTimer[NTIMERS];
    bool m_TimerIsRunning[NTIMERS];
    unsigned int m_Counter[NCOUNTERS];
    double m_FFTTimers[5][NVALTYPES];
    double m_InterpTimers[4][NVALTYPES];

    unsigned int m_NumThreads;

    unsigned int m_LineLength;

    bool m_StoreTimeSeries;
    bool m_MonitorCFLCondition;
    bool m_InCompressible;

    bool m_ReadImagesFromFile;
    bool m_WriteImages;
    bool m_WriteLogFiles;

    ScalarType m_Beta[2]; ///< regularization weight
    ScalarType m_Sigma;

    int m_Verbosity;

};

} // end of namespace


#endif
