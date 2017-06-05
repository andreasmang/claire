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
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGOPT_H_
#define _REGOPT_H_

// global includes
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

// local includes
#include "RegUtils.hpp"




namespace reg {




// flags for hyperbolic PDE solvers
enum PDESolverType {
    RK2,   ///< flag for RK2 solver
    RK2A,  ///< flag for stabilized RK2 solver
    SL,    ///< flag for SL
};


// flags for regularization norms
enum RegNormType {
    L2,    ///< flag for L2-norm
    H1,    ///< flag for H1-norm
    H2,    ///< flag for H2-norm
    H3,    ///< flag for H3-norm
    H1SN,  ///< flag for H1-seminorm
    H2SN,  ///< flag for H2-seminorm
    H3SN,  ///< flag for H3-seminorm
};


enum ParaContType {
    PCONTOFF,
    PCONTBINSEARCH,
    PCONTREDUCESEARCH,
    PCONTINUATION
};


// flags for regularization norms
enum GlobalMethType{
    NOGM,        ///< no linesearch / globalization
    ARMIJOLS,    ///< armijo line search
    MTLS,        ///< more-thuente line search
    GPCGLS,      ///<
    OWARMIJOLS,  ///<
    IPMLS,       ///<
};




/*! hessian mat vec type */
enum HessianMatVecType {
    DEFAULTMATVEC,
    PRECONDMATVEC,
    PRECONDMATVECSYM,
};


// flags for optimization methods
enum OptMeth {
    GAUSSNEWTON,  ///< Gauss-Newton approximation
    FULLNEWTON,   ///< full Newton
    GRADDESCENT,  ///< gradient descent (gradient in sobolev space)
};


// flags for preconditioner methods
enum PrecondMeth {
    INVREG,    ///< inverse regularization operator
    TWOLEVEL,  ///< 2 level preconditioner
    NOPC,      ///< no preconditioner
};


// flags for optimization
enum FSeqType {
    QDFS,  ///< quadratic forcing sequence
    SLFS,  ///< superliner forcing sequence
    NOFS,  ///< no forcing sequence
};


// flags for preconditioners
enum KrylovSolverType {
    PCG,     ///< pcg
    FCG,     ///< flexible cg
    CHEB,    ///< chebyshev method
    GMRES,   ///< gmres
    FGMRES,  ///< flexible gmres
};


// high level flags for solver
enum SolveType {
    FAST_AGG,     ///< H1-regularization; few iterations; low accuracy
    FAST_SMOOTH,  ///< H2-regularization; few iterations; low accuracy
    ACC_AGG,      ///< H1-regularization; accurate solution;
    ACC_SMOOTH,   ///< H2-regularization; accurate solution;
    NOTSET,       ///< no predefined solver
};


// flags for timers
enum TimerType {
    T2SEXEC = 0,  ///< time to solution (execution time)
    PDEEXEC,      ///< pde solves (execution time)
    HMVEXEC,      ///< hessian mat vec (execution time)
    PMVSETUP,     ///< setup time for preconditioner
    PMVEXEC,      ///< precond mat vec (execution time)
    GRADEXEC,     ///< gradient evaluation (execution time)
    OBJEXEC,      ///< objective evluation (execution time)
    FFTSETUP,     ///< fft setup time
    FFTSELFEXEC,  ///< fft execution time
    IPSELFEXEC,   ///< execution time for interpolation
    NTIMERS,      ///< to allocate the timers
};



enum FFTTimers {
    FFTTRANSPOSE = 0,  ///< contains all communication time
    FFTSHUFFLE   = 1,  ///< contained in FFTTRANSPOSE
    FFTCOMM      = 2,  ///< contained in FFTTRANSPOSE
    FFTRESHUFFLE = 3,  ///< contained in FFTTRANSPOSE
    FFTEXECUTE   = 4,  ///< execution time
    FFTNAN       = 5,  ///< currently empty
    FFTHADAMARD  = 6,  ///< hadamard products for differential operators
    NFFTTIMERS   = 7,  ///< total number of itmers
};


enum FFTCounters {
    FFTGRAD = 4,  ///< number of ffts for gradient opteration
    FFTDIV  = 4,  ///< number of ffts for divergence operation
};



// counters (number of operations)
enum CounterType {
    PDESOLVE = 0,  ///< PDE solves
    HESSMATVEC,    ///< # hessian matvecs
    PCMATVEC,      ///< preconditioner matvecs
    OBJEVAL,       ///< objective evaluations
    GRADEVAL,      ///< gradient evaluations
    IPVEC,         ///< interpolation execution time
    IP,            ///< interpolation execution time
    FFT,           ///< fft evaluations
    ITERATIONS,    ///< number of outer iterations
    NCOUNTERS,     ///< to allocate the counters
};


enum RegModel {
    COMPRESSIBLE,
    RELAXEDSTOKES,
    STOKES,
};


/*! flags for timers */
enum LogType{
    LOGRES,
    LOGCONV,
    LOGKSPRES,
    LOGJAC,
    LOGLOAD,
    NLOGFLAGS
};


/*! flags for timers */
enum StopCondType{
    GRAD,
    GRADOBJ,
};


struct ReadWriteFlags {
    bool readfiles;
    bool readvelocity;
    bool timeseries;
    bool iterates;
    bool defgrad;
    bool detdefgrad;
    bool velnorm;
    bool residual;
    bool defmap;
    bool templateim;
    bool referenceim;
    bool deftemplate;
    bool deffield;
    bool results;
    std::string extension;

    std::string xfolder;                ///< identifier for folder to write results to
    std::string ifolder;                ///< identifier for folder to read in results from
    std::vector < std::string > mt;     ///< template image file name
    std::vector < std::string > mr;     ///< reference image file name
    std::string vx1;                    ///< x1-velocity field file name
    std::string vx2;                    ///< x2-velocity field file name
    std::string vx3;                    ///< x3-velocity field file name
};


/* parameters for domain */
struct Domain {
    IntType isize[3];           ///< size of grid in spatial domain for mpi proc
    IntType istart[3];          ///< start index in spatial domain for mpi proc
    IntType nl;                 ///< number of grid points for each mpi proc
    IntType ng;                 ///< number of grid points (global)
    ScalarType hx[3];           ///< spatial grid cell size
    IntType nx[3];              ///< spatial grid size
    IntType nt;                 ///< number of time points
    IntType nc;                 ///< number of components/images to be registered
    ScalarType timehorizon[2];  ///< time horizon
};


/* parameters for optimization */
struct Optimization {
    int maxit;                          ///< maximal number of (outer) iterations
    int minit;                          ///< minimal number of (outer) iterations (for parameter continuation)
    ScalarType tol[3];                  ///< tolerances for optimization
    OptMeth method;                     ///< optimization method
    GlobalMethType glmethod;            ///< method for globalization (line search; trust region; ...)
    StopCondType stopcond;              ///< stopping conditions
    bool fastpresolve;                  ///< flag to switch on fast presolve
    bool fastsolve;                     ///< flag to switch on fast (inaccurate) solve
    ScalarType presolvetol[3];          ///< tolerances for presolve
    int presolvemaxit;                  ///< maximal iterations for presolve
    bool usezeroinitialguess;           ///< use zero initial guess for computing the initial gradient
    int solutionstatus;

    // debug
    bool derivativecheckenabled;        ///< use zero initial guess for computing the initial gradient
};


/* parameters for krylov solver */
struct KrylovSolver {
    int maxit;                ///< max number of iterations for krylov solver
    IntType iter;             ///< current number of iterations for krylov solver
    ScalarType tol[3];        ///< tolerances for krylov method
    FSeqType fseqtype;        ///< forcing sequence type
    std::string name;         ///< name of krylov solver
    ScalarType reltol;        ///< relative tolerance for krylov solver
    ScalarType g0norm;        ///< initial norm of gradient (to normalize stopping condition)
    bool g0normset;           ///< flag to identify if initial norm of gradient has been set
    KrylovSolverType solver;  ///< flag for krylov solver

    ScalarType pctol[3];            ///< tolerances for krylov method (preconditioner)
    IntType pcmaxit;                ///< tolerances for krylov method (preconditioner)
    PrecondMeth pctype;             ///< flag for type of preconditioner
    HessianMatVecType matvectype;   ///< flag for the type of hessian matvec
    std::string pcname;             ///< name of preconditioner
    KrylovSolverType pcsolver;      ///< solver for preconditioner
    bool pcsetupdone;               ///< flag to indicate if setup of preconditioner is done
    ScalarType pctolscale;          ///< tolerance scaling for preconditioner; default: 1E-1
    ScalarType pcgridscale;         ///< this is for the two level preconditioner; defines scale for grid size change; default: 2
    bool usepetsceigest;            ///< in cheb method we need to estimate eigenvalues; use petsc implementation
    bool reesteigvals;              ///< flag to reestimate eigenvalues every iteration
    bool eigvalsestimated;          ///< flag if eigenvalues have already been estimated
    bool checkhesssymmetry;         ///< check symmetry of hessian operator
    ScalarType hessshift;           ///< perturbation to hessian operator
};


/* parameter for parameter continuation (regularization parameter) */
struct ParCont {
//    static constexpr ScalarType betavminh1 = 1E-4;    ///< minimal regularization parameter for h1 type norm
//    static constexpr ScalarType betavminh2 = 1E-7;    ///< minimal regularization parameter for h2 type norm
    static constexpr ScalarType betavminh1 = 1E-9;      ///< minimal regularization parameter for h1 type norm
    static constexpr ScalarType betavminh2 = 1E-9;      ///< minimal regularization parameter for h2 type norm
    static constexpr ScalarType betascale = 1E-1;       ///< default reduction factor (one order of magnitude)
    static constexpr ScalarType dbetascale = 1E-2;      ///< default reduction factor (one order of magnitude)
    static const int maxsteps = 10;                     ///< max number of steps
    ParaContType strategy;                              ///< flag for parameter continuation strategy
    bool enabled;                                       ///< flag: parameter continuation using different strategies
    ScalarType targetbeta;                              ///< target regularization parameter
    ScalarType beta0;                                   ///< initial regularization parameter
    ScalarType stepsize;
};


/* parameter for scale continuation */
struct ScaleCont {
    bool enabled;
    static const int maxlevel = 6;
    ScalarType sigma[3][maxlevel];
};


/* parameter for grid continuation */
struct GridCont {
    bool enabled;
    static const int minlevels = 3;
    int nlevels;
    int minlevel;
    std::vector<int> maxit;
    IntType nxmin;
    std::vector< std::vector<IntType> > nx;
    std::vector< std::vector<IntType> > isize;
    std::vector< std::vector<IntType> > istart;
    std::vector< std::vector<IntType> > osize;
    std::vector< std::vector<IntType> > ostart;
    std::vector<IntType> nl;
    std::vector<IntType> ng;
    std::vector<IntType> nalloc;
};


/* parameters/containers for monitoring registration */
struct RegMonitor {
    bool detdgradenabled;       ///< flag to monitor jacobian during iterations
    ScalarType detdgradmin;     ///< min value of determinant of deformation gradient det(grad(y))
    ScalarType detdgradmax;     ///< max value of determinant of deformation gradient det(grad(y))
    ScalarType detdgradmean;    ///< mean value of determinant of deformation gradient det(grad(y))
    ScalarType detdgradbound;   ///< lower bound of determinant of deformation gradient det(grad(y))
    ScalarType jval;            ///< value of objective functional
    ScalarType dval;            ///< value of distance measure
    ScalarType rval;            ///< value of regularization functional
};


/* regularization operators */
struct RegNorm {
    ScalarType beta[4];  ///< regularization parameter
    RegNormType type;    ///< flag for regularization norm
};


struct FourierTransform{
    accfft_plan_t<ScalarType, ComplexType, FFTWPlanType>* plan;  ///< accfft plan
    //accfft_plan_t<double, Complex, fftw_plan>* plan;  ///< accfft plan
    //accfft_plan* plan;  ///< accfft plan
    MPI_Comm mpicomm;   ///< communicator for accfft
    bool mpicommexists; ///< communicator for accfft
    IntType nalloc;     ///< size for allocation in fourier domain
    IntType osize[3];   ///< size of grid in fourier domain for mpi proc
    IntType ostart[3];  ///< start index in fourier domain for mpi proc
};


struct RegFlags {
    bool applysmoothing;  ///< apply smoothing to images
    bool applyrescaling;  ///< apply rescaling to images (map the intensity range to [0,1])
    bool detdefgradfromdeffield;
    bool invdefgrad;
    bool checkdefmapsolve;
    bool runninginversion;
};


struct PDESolver {
    PDESolverType type;
    int rungekuttaorder;
    int interpolationorder;
    ScalarType cflnumber;
    bool monitorcflnumber;
    bool adapttimestep;
};


/*! parameter for grid continuation */
struct Logger {
    enum TimerValue {LOG = 0, MIN, MAX, AVG, NVALTYPES};
    std::vector<ScalarType> distance;        ///< convergence for residual
    std::vector<ScalarType> regularization;  ///< convergence for regularization
    std::vector<ScalarType> objective;       ///< convergence for objective
    std::vector<int> outeriterations;        ///< iterations of solver
    std::vector<ScalarType> kspresidual;     ///< residual of krylov method
    std::vector<int> kspiterations;          ///< iterations of krylov method
    ScalarType finalresidual[4];
    bool enabled[NLOGFLAGS];

    bool memoryusage;
    double timer[NTIMERS][NVALTYPES];
    double temptimer[NTIMERS];
    bool timerruns[NTIMERS];
    unsigned int counter[NCOUNTERS];
    double ffttimers[NFFTTIMERS][NVALTYPES];
    double iptimers[4][NVALTYPES];
};




class RegOpt {
 public:
    typedef RegOpt Self;

    RegOpt();
    RegOpt(int, char**);
    RegOpt(const RegOpt&);
    virtual ~RegOpt();
    void Copy(const RegOpt&);

    // spatial grid
    inline void SetNumGridPoints(int i, IntType nx) {
        this->m_Domain.nx[i] = nx;
    }
    inline ScalarType GetLebesqueMeasure(void) {
        return  this->m_Domain.hx[0]
               *this->m_Domain.hx[1]
               *this->m_Domain.hx[2];
    }

    inline Domain GetDomainPara() {return this->m_Domain;}
    inline GridCont GetGridContPara() {return this->m_GridCont;}
    inline ScaleCont GetScaleContPara() {return this->m_ScaleCont;}
    inline ParCont GetParaCont() {return this->m_ParaCont;}
    ScalarType GetBetaMinParaCont();

    inline FourierTransform GetFFT() {return this->m_FFT;}
    inline RegFlags GetRegFlags() {return this->m_RegFlags;}

    // monitor functions
    inline RegMonitor GetRegMonitor() {return this->m_RegMonitor;}
    inline void SetDetDefGradBound(ScalarType value) {this->m_RegMonitor.detdgradbound = value;}
    inline void SetDetDGradMin(ScalarType value) {this->m_RegMonitor.detdgradmin = value;}
    inline void SetDetDGradMax(ScalarType value) {this->m_RegMonitor.detdgradmax = value;}
    inline void SetDetDGradMean(ScalarType value) {this->m_RegMonitor.detdgradmean = value;}
    inline void EnableBoundOnDetDGrad() {this->m_RegMonitor.detdgradenabled = true;}
    inline void DisableBoundOnDetDGrad() {this->m_RegMonitor.detdgradenabled = false;}
    inline void SetJVal(ScalarType value){this->m_RegMonitor.jval = value;}
    inline void SetRVal(ScalarType value){this->m_RegMonitor.rval = value;}
    inline void SetDVal(ScalarType value){this->m_RegMonitor.dval = value;}

    inline void DisableInversion() {this->m_RegFlags.runninginversion = false;}
    inline void EnableInversion() {this->m_RegFlags.runninginversion = true;}

    inline Optimization GetOptPara() {return this->m_OptPara;}
    inline void SetOptTol(int i, ScalarType value) {this->m_OptPara.tol[i] = value;}
    inline void SetOptMaxIter(int value) {this->m_OptPara.maxit = value;}
    inline void SetSolutionStatus(int value) {this->m_OptPara.solutionstatus = value;}

    inline void ComputeInvDetDefGrad(bool flag) {
        this->m_RegFlags.invdefgrad = flag;
    }

    inline void EnableRescaling() {this->m_RegFlags.applyrescaling = true;}
    inline void DisableRescaling() {this->m_RegFlags.applyrescaling = false;}

    inline void EnableSmoothing() {this->m_RegFlags.applysmoothing = true;}
    inline void DisableSmoothing() {this->m_RegFlags.applysmoothing = false;}

    PetscErrorCode GetSizes(IntType*, IntType&, IntType&);
    PetscErrorCode GetSizes(IntType*, IntType*, IntType*);

    // time horizon, step size, and ....
    inline void SetNumTimePoints(IntType nt) {
            this->m_Domain.nt = nt;
    }
    inline ScalarType GetTimeStepSize(void) {
        return (this->m_Domain.timehorizon[1] - this->m_Domain.timehorizon[0])
               /static_cast<ScalarType>(this->m_Domain.nt);
    }
    inline void SetNumImageComponents(IntType n) {
        this->m_Domain.nc = n;
    }

    inline ReadWriteFlags GetReadWriteFlags() {return this->m_ReadWriteFlags;}
    inline RegModel GetRegModel(void) {return this->m_RegModel;}

    /* do setup for grid continuation */
    PetscErrorCode SetupGridCont();

    // regularization
    inline RegNorm GetRegNorm() {return this->m_RegNorm;}
    inline void SetRegNormType(RegNormType flag) {this->m_RegNorm.type = flag;}
    inline void SetBeta(int i, ScalarType beta) {this->m_RegNorm.beta[i] = beta;}
    inline ScalarType GetBeta(int i) {return this->m_RegNorm.beta[i];}
    inline void SetSigma(int i, ScalarType sigma) {this->m_Sigma[i] = sigma;}
    inline ScalarType GetSigma(int i) {return this->m_Sigma[i];}

    // solver flags
    inline PDESolver GetPDESolverPara(void) {return this->m_PDESolver;}
    inline void SetPDESolverType(PDESolverType flag) {this->m_PDESolver.type = flag;}

    // krylov method
    inline KrylovSolver GetKrylovSolverPara() {return this->m_KrylovSolverPara;}
    inline void SetRelTolKrylovMethod(ScalarType value) {this->m_KrylovSolverPara.reltol = value;}
    inline void SetMaxItKrylovMethod(int value) {this->m_KrylovSolverPara.maxit = value;}
    inline void SetKrylovIter(IntType value) {this->m_KrylovSolverPara.iter = value;}
    inline void SetInitialGradNormKrylovMethod(ScalarType value) {
        this->m_KrylovSolverPara.g0norm = value;
        this->m_KrylovSolverPara.g0normset = true;
    }
//    inline void InitialGradNormSet(bool flag) {this->m_KrylovSolverPara.g0normset = flag;}
    inline void KrylovMethodEigValsEstimated(bool flag) {this->m_KrylovSolverPara.eigvalsestimated = flag;}
    inline void PrecondSetupDone(bool flag) {this->m_KrylovSolverPara.pcsetupdone = flag;}

    // flag for setup
    inline bool SetupDone() {return this->m_SetupDone;}
    inline bool StoreCheckPoints() {return this->m_StoreCheckPoints;}

    // timers and counters
    inline unsigned int GetCounter(CounterType id) {
        return this->m_Counter[id];
    }
    inline void IncrementCounter(CounterType id) {
        this->m_Counter[id]++;
    }
    inline void IncrementCounter(const CounterType id, const int i) {
        this->m_Counter[id] += i;
    }

    inline void GetTimer(const TimerType id, double* wtime) {
        wtime[0] = this->m_Timer[id][MIN];
        wtime[1] = this->m_Timer[id][MAX];
        wtime[2] = this->m_Timer[id][AVG];
    }

    unsigned int GetNumThreads() {return this->m_NumThreads;}
    inline int GetNetworkDims(const int i) {return this->m_CartGridDims[i];}
    inline void IncreaseFFTTimers(const double timers[NFFTTIMERS]) {
        for (int i = 0; i < NFFTTIMERS; ++i) {
            this->m_FFTTimers[i][LOG] += timers[i];
        }
    }
    inline void IncreaseInterpTimers(const double timers[4]) {
        for (int i = 0; i < 4; ++i) {
            this->m_InterpTimers[i][LOG] += timers[i];
        }
    }

    inline Logger GetLogger() {
        return this->m_Log;
    }
    inline void LogKSPResidual(const int i, const ScalarType value){
        this->m_Log.kspresidual.push_back(value);
        this->m_Log.kspiterations.push_back(i);
    }
    inline void LogConvergence(const int i, const ScalarType J, const ScalarType D, const ScalarType R){
        this->m_Log.distance.push_back(D);
        this->m_Log.regularization.push_back(R);
        this->m_Log.objective.push_back(J);
        this->m_Log.outeriterations.push_back(i);
    }
    inline void LogFinalResidual(const int i, const ScalarType value){
        this->m_Log.finalresidual[i] = value;
    }

    inline void SetVerbosity(int value) {this->m_Verbosity = value;}
    inline int GetVerbosity() {
        return this->m_Verbosity;
    }
    int GetLineLength() {
        return this->m_LineLength;
    }

    ScalarType ComputeFFTScale();

    PetscErrorCode StartTimer(TimerType);
    PetscErrorCode StopTimer(TimerType);
    PetscErrorCode ResetTimers(void);
    PetscErrorCode ResetTimer(TimerType);
    PetscErrorCode ResetCounters(void);
    PetscErrorCode ResetCounter(CounterType);
    PetscErrorCode ProcessTimers(void);

    virtual PetscErrorCode DisplayOptions(void);
    PetscErrorCode DisplayTimeToSolution(void);
    PetscErrorCode WriteLogFile(void);
    PetscErrorCode DoSetup(bool dispteaser = true);

    inline void Enter(std::string fname) {
        #ifdef _REG_DEBUG_
        std::stringstream ss;
        ss << std::string(this->m_Indent++, ' ') << ">> " << fname << std::endl;
        PetscPrintf(PETSC_COMM_WORLD, ss.str().c_str());
        #endif
    }

    inline void Exit(std::string fname) {
        #ifdef _REG_DEBUG_
        std::stringstream ss;
        ss << std::string(--this->m_Indent, ' ') << "<< " << fname << std::endl;
        PetscPrintf(PETSC_COMM_WORLD, ss.str().c_str());
        #endif
    }

    PetscErrorCode EnableFastSolve();

 protected:
    virtual PetscErrorCode Initialize(void);
    PetscErrorCode InitializeFFT();
    PetscErrorCode DestroyFFT();
    virtual PetscErrorCode ClearMemory(void);
    virtual PetscErrorCode ParseArguments(int, char**);
    virtual PetscErrorCode Usage(bool advanced = false);
    virtual PetscErrorCode CheckArguments(void);
    PetscErrorCode SetPresetParameters();
    PetscErrorCode WriteWorkLoadLog();
    PetscErrorCode WriteWorkLoadLogReadable();
    PetscErrorCode WriteWorkLoadLog(std::ostream&);
    PetscErrorCode WriteWorkLoadLogReadable(std::ostream&);
    PetscErrorCode WriteKSPLog();
    PetscErrorCode WriteConvergenceLog();
    PetscErrorCode WriteFinalResidualLog();

    enum TimerValue {LOG = 0, MIN, MAX, AVG, NVALTYPES};

    PDESolver m_PDESolver;              ///< flag for PDE solver
    Optimization m_OptPara;             ///< optimization parameters
    KrylovSolver m_KrylovSolverPara;    ///< parameters for krylov solver
    ReadWriteFlags m_ReadWriteFlags;    ///< flags for io
    RegMonitor m_RegMonitor;            ///< monitor for registration
    RegNorm m_RegNorm;                  ///< parameters for regularization model
    ParCont m_ParaCont;                 ///< flags for parameter continuation
    GridCont m_GridCont;                ///< flags for grid continuation
    ScaleCont m_ScaleCont;              ///< flags for scale continuation
    Domain m_Domain;                    ///< parameters for spatial domain
    RegModel m_RegModel;                ///< flag for particular registration model
    FourierTransform m_FFT;             ///< parameters for FFT/accfft
    RegFlags m_RegFlags;                ///< flags for registration
    SolveType m_SolveType;              ///< solver
    Logger m_Log;                       ///< log

    double m_Timer[NTIMERS][NVALTYPES];
    double m_TempTimer[NTIMERS];
    bool m_TimerIsRunning[NTIMERS];
    unsigned int m_Counter[NCOUNTERS];
    double m_FFTTimers[NFFTTIMERS][NVALTYPES];
    double m_FFTAccumTime;
    double m_InterpTimers[4][NVALTYPES];
    double m_IPAccumTime;
    double m_IPSlowest;
    double m_TTSSlowest;
    double m_FFTSlowest;
    int m_IDSlowest;

    int m_CartGridDims[2];
    unsigned int m_NumThreads;
//    const unsigned int m_LineLength = 101; //C++ 11 feature
    unsigned int m_LineLength;

    bool m_SetupDone;
    bool m_StoreCheckPoints;
    ScalarType m_Sigma[3];

    int m_Indent;
    int m_Verbosity;
};




}  // namespace reg




#endif  // _REGOPT_H_
