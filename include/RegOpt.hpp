/**
 *  Copyright (c) 2015-2016.
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
 *
*/


#ifndef _REGOPT_H_
#define _REGOPT_H_

//global includes
#include <sstream>

// local includes
#include "RegUtils.hpp"


namespace reg
{




// flags for hyperbolic PDE solvers
enum PDESolver
{
    RK2,  ///< flag for RK2 solver
    RK2A, ///< flag for stabilized RK2 solver
    SL,   ///< flag for SL
};


// flags for regularization norms
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


// flags for preconditioners
enum PCSolverType
{
    PCPCG,    ///< pcg as preconditioner
    PCFCG,    ///< flexible cg as preconditioner
    PCCHEB,   ///< chebyshev method as preconditioner
    PCGMRES,  ///< gmres solver as preconditioner
    PCFGMRES, ///< flexible gmres as preconditioner
};


// high level flags for solver
enum SolveType
{
    FAST_AGG, // H1-regularization; few iterations; low accuracy
    FAST_SMOOTH, // H2-regularization; few iterations; low accuracy
    ACC_AGG, // H1-regularization; accurate solution;
    ACC_SMOOTH, // H2-regularization; accurate solution;
    NOTSET, // no predefined solver
};


// flags for timers
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


// counters (number of operations)
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


enum RegModel
{
    COMPRESSIBLE,
    STOKES,
    RELAXEDSTOKES
};




class RegOpt
{

public:

    typedef RegOpt Self;

    RegOpt();
    RegOpt(int,char**);
    ~RegOpt();

    // number of points
    inline IntType GetNLocal(void){return this->m_Domain.nlocal;};
    inline IntType GetNGlobal(void){return this->m_Domain.nglobal;};
    inline IntType GetISize(int i){return this->m_Domain.isize[i];};
    inline IntType GetIStart(int i){return this->m_Domain.istart[i];};

    // spatial grid
    inline void SetNumGridPoints(int i,IntType nx){this->m_Domain.nx[i] = nx;};
    inline IntType GetNumGridPoints(int i){return this->m_Domain.nx[i];};
    inline ScalarType GetSpatialStepSize(int i){return this->m_Domain.hx[i];};
    inline ScalarType GetLebesqueMeasure(void)
    {
        return  this->m_Domain.hx[0]
               *this->m_Domain.hx[1]
               *this->m_Domain.hx[2];
    }


    // time horizon, step size, and ....
    inline ScalarType GetTimeHorizon(int i){return this->m_TimeHorizon[i];};
    inline IntType GetNumTimePoints(void){return this->m_nt;};
    inline void SetNumTimePoints(IntType nt){ this->m_nt=nt; };
    inline ScalarType GetTimeStepSize(void){
        return (this->m_TimeHorizon[1] - this->m_TimeHorizon[0])
               /static_cast<ScalarType>(this->m_nt);
    };

    accfft_plan* GetFFTPlan(){return this->m_FFTPlan;};
    MPI_Comm GetComm(){return this->m_Comm;};

    // control input and output
    inline std::string GetIFolder(void){return this->m_IFolder;};
    inline std::string GetXFolder(void){return this->m_XFolder;};
    inline std::string GetXExtension(void){return this->m_XExtension;};
    inline std::string GetTemplateFN(){return this->m_TemplateFN;};
    inline std::string GetReferenceFN(){return this->m_ReferenceFN;};
    inline bool ReadImagesFromFile(){return this->m_ReadImagesFromFile;};
    inline void ReadImagesFromFile(bool flag){this->m_ReadImagesFromFile=flag;};

    // flags
    inline bool StoreTimeSeries(){return this->m_StoreTimeSeries;};
    inline bool StoreIterates(){return this->m_StoreIterates;};

    // registration model
    inline RegModel GetRegModel(void){return this->m_RegModel;};

    // regularization
    inline RegNorm GetRegNorm(void){return this->m_Regularization.norm;};
    inline void SetRegularizationWeight(ScalarType beta){
        this->m_Regularization.beta[0]=beta;
        this->m_Regularization.beta[1]=beta;
    };
    inline ScalarType GetRegularizationWeight(void){return this->m_Regularization.beta[0];};
    inline ScalarType GetRegularizationWeight(int i){return this->m_Regularization.beta[i];};

    // smoothing
    inline ScalarType GetSigma(void){return this->m_Sigma;};

    // solver flags
    inline PDESolver GetPDESolver(void){return this->m_PDESolver;};
    inline PrecondMeth GetPrecondMeth(void){return this->m_PrecondMeth;};
    inline PCSolverType GetPCSolverType(){return this->m_PCSolverType;};
    inline FSeqType GetFSeqType(void){return this->m_KKTSolverPara.fseqtype;};

    // jacobians
    inline ScalarType GetJacMin(){return this->m_RegMonitor.jacmin;};
    inline ScalarType GetJacMax(){return this->m_RegMonitor.jacmax;};
    inline ScalarType GetJacMean(){return this->m_RegMonitor.jacmean;};
    inline ScalarType GetJacBound(){return this->m_RegMonitor.jacbound;};

    inline void SetJacMin(ScalarType value){this->m_RegMonitor.jacmin=value;};
    inline void SetJacMax(ScalarType value){this->m_RegMonitor.jacmax=value;};
    inline void SetJacMean(ScalarType value){this->m_RegMonitor.jacmean=value;};

    inline bool MonitorJacobian(){return this->m_RegMonitor.monitorJAC;};
    inline bool MonitorCFLCondition(){return this->m_RegMonitor.monitorCFL;};
    inline void MonitorCFLCondition(bool flag){this->m_RegMonitor.monitorCFL=flag;};

    // flag for setup
    inline bool SetupDone(){return this->m_SetupDone;};

    // parameter continuation
    inline bool DoRegParaBinarySearch(){return this->m_ParameterCont.strategy[0];};
    inline bool DoRegParaReductionSearch(){return this->m_ParameterCont.strategy[1];};
    inline bool DoRegParaContinuation(){return this->m_ParameterCont.strategy[2];};

    inline int GetMaxStepsParaCont(){return this->m_ParameterCont.maxsteps;};
    inline ScalarType GetBetaScaleParaCont(){return this->m_ParameterCont.betascale;};
    inline ScalarType GetDeltaBetaScaleParaCont(){return this->m_ParameterCont.dbetascale;};
    ScalarType GetBetaMinParaCont();

    // timers and counters
    inline unsigned int GetCounter(CounterType id){return this->m_Counter[id];};
    inline void IncrementCounter(CounterType id){this->m_Counter[id]++;};
    inline void IncrementCounter(CounterType id, unsigned int i){this->m_Counter[id]+= i;};

    inline void GetTimer(TimerType id,double* wtime){
        wtime[0] = this->m_Timer[id][MIN];
        wtime[1] = this->m_Timer[id][MAX];
        wtime[2] = this->m_Timer[id][AVG];
    };

    inline int GetNetworkDims(int i){return this->m_CartGridDims[i];};
    inline void IncreaseFFTTimers(double timers[5]){
        for(int i=0; i < 5; ++i) this->m_FFTTimers[i][LOG]+=timers[i];
    };
    inline void IncreaseInterpTimers(double timers[4]){
        for(int i=0; i < 4; ++i) this->m_InterpTimers[i][LOG]+=timers[i];
    };

    inline ScalarType GetOptTol(int i){return this->m_OptPara.tol[i];};
    inline int GetOptMaxit(){return this->m_OptPara.maxit;};
    inline OptMeth GetOptMeth(void){return this->m_OptPara.method;};

    inline ScalarType GetKKTSolverTol(int i){return this->m_KKTSolverPara.tol[i];};
    inline int GetKKTMaxit(){return this->m_KKTSolverPara.maxit;};
    inline unsigned int GetDDNumDomains(){return this->m_DD.n;};

    inline bool WriteImages(){return this->m_WriteImages;};
    inline void WriteImages(bool flag){this->m_WriteImages = flag;};
    inline bool LoggingEnabled(){return this->m_WriteLogFiles;};
    inline int GetVerbosity(){return this->m_Verbosity;};

    int GetLineLength(){return this->m_LineLength;};

    ScalarType ComputeFFTScale();

    PetscErrorCode StartTimer(TimerType);
    PetscErrorCode StopTimer(TimerType);
    PetscErrorCode ResetTimers(void);
    PetscErrorCode ResetTimer(TimerType);
    PetscErrorCode ResetCounters(void);
    PetscErrorCode ResetCounter(CounterType);
    PetscErrorCode ProcessTimers(void);

    PetscErrorCode DisplayOptions(void);
    PetscErrorCode DisplayTimeToSolution(void);
    PetscErrorCode WriteLogFile(void);
    PetscErrorCode DoSetup(void);

private:

    PetscErrorCode Usage(bool advanced=false);
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode ParseArguments(int,char**);
    PetscErrorCode CheckArguments(void);
    PetscErrorCode SetPresetParameters();

    enum TimerValue{LOG=0,MIN,MAX,AVG,NVALTYPES};

    IntType m_nt; ///< number of time points
    ScalarType m_TimeHorizon[2];


    // parameters for optimization
    struct Domain{
        IntType isize[3];
        IntType istart[3];
        IntType osize[3];
        IntType ostart[3];
        IntType nlocal;
        IntType nglobal;
        ScalarType hx[3];
        IntType nx[3];
    };


    // parameters for optimization
    struct Optimization{
        int maxit;
        ScalarType tol[3];
        OptMeth method;
    };

    // parameters for KKT solver
    struct KKTSolver{
        int maxit;
        ScalarType tol[3];
        FSeqType fseqtype; ///<forcing sequence type
    };

    // parameters for KKT solver
    struct DomainDecomposition{
        unsigned int n;
    };

    // parameters for parameter continuation
    struct ParameterContinuation{
        static const ScalarType betavminh1=1E-3; ///< minimal regularization parameter for h1 type norm
        static const ScalarType betavminh2=1E-6; ///< minimal regularization parameter for h2 type norm
        static const int maxsteps = 10; ///< max number of steps
        static const ScalarType betascale = 1E-1; ///< default reduction factor (one order of magnitude)
        static const ScalarType dbetascale = 1E-2; ///< default reduction factor (one order of magnitude)
        bool strategy[3]; ///< flag: parameter continuation using different strategies
        ScalarType targetbeta;
    };

    // parameters for parameter continuation
    struct RegMonitor{
        bool monitorJAC;
        bool monitorCFL;
        ScalarType jacmin; ///< min value of jacobian
        ScalarType jacmax; ///< max value of jacobian
        ScalarType jacmean; ///< mean value of jacobian
        ScalarType jacbound; ///< lower bound for jacobian
    };

    struct Regularization{
        ScalarType beta[3]; ///< regularization parameter
        RegNorm norm; ///< flag for regularization norm
    };

    // parameters for KKT solver
    struct MultiScale{
        unsigned int n;
    };

    accfft_plan* m_FFTPlan;
    MPI_Comm m_Comm;
    int m_CartGridDims[2];

    Optimization m_OptPara; ///< optimization parameters
    PDESolver m_PDESolver; ///< flag for PDE solver
    KKTSolver m_KKTSolverPara; ///< parameters for KKT solver
    RegMonitor m_RegMonitor;  ///< monitor for registration
    PrecondMeth m_PrecondMeth; ///< flag for preconditioner
    PCSolverType m_PCSolverType; ///< flag for KSP solver for precond
    Regularization m_Regularization; ///< parameters for regularization model
    ParameterContinuation m_ParameterCont; ///< flags for parameter continuation
    DomainDecomposition m_DD; ///< domain decomposition
    Domain m_Domain; ///< domain decomposition
    RegModel m_RegModel;

    std::string m_XFolder; ///< identifier for folder to write results to
    std::string m_XExtension; ///< identifier for extension of files to be written to file
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
    const unsigned int m_LineLength = 101;
    bool m_StoreTimeSeries;
    bool m_StoreIterates;
    bool m_SetupDone;

    bool m_WriteImages;
    bool m_WriteLogFiles;
    bool m_ReadImagesFromFile;

    ScalarType m_Sigma;
    SolveType m_SolveType;

    int m_Verbosity;

};

} // end of namespace


#endif
