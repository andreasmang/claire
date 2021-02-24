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

#ifndef _REGOPT_H_
#define _REGOPT_H_

// local includes
#include "CLAIREUtils.hpp"
#include "Spectral.hpp"
#include <vector>
#include <string>
//#include "ARGVParser.hpp"

namespace reg {




// flags for hyperbolic PDE solvers
enum PDESolverType {
    RK2,   ///< flag for RK2 solver
    RK2A,  ///< flag for stabilized RK2 solver
    SL,    ///< flag for SL
};


// flags for hyperbolic PDE solvers
enum PDEType {
    CONTINUITYEQ,   ///< identifier for continuity equation
    TRANSPORTEQ,    ///< identifier for transport equation
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



// flags for regularization norms
enum DMType {
    SL2,    ///< flag for squared L2-norm
    NCC,    ///< flag for normalized cross-correlation
    SL2AUX, ///< flag for squared L2-norm (coupling)
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

enum DifferentiationType {
    SPECTRAL,    ///< spectral differentiation
    FINITE,      ///< finite difference differentiation
    MIXED        ///< FD for 1st order, Spectral for higher order
};


/*! hessian mat vec type */
enum HessianMatVecType {
    DEFAULTMATVEC,
    H0MATVEC,
    PRECONDMATVEC,
    PRECONDMATVECSYM,
};


/*! hessian mat vec type */
enum GradientType {
    L2GRAD,
    SGRAD,
    SYMSGRAD,
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
    H0,        ///< H(v=0)^-1 preconditioner
    H0MG,      ///< H(v=0)^-1 preconditioner with multi grid
};


// flags for optimization
enum FSeqType {
    QDFS,  ///< quadratic forcing sequence
    SLFS,  ///< superliner forcing sequence
    NOFS,  ///< no forcing sequence
};


// flags for preconditioners
enum KrylovMethodType {
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
    GPUCOMP,      ///< time spent on gpu
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
    LOGDIST,
    LOGCONV,
    LOGKSPRES,
    LOGJAC,
    LOGLOAD,
    LOGGRAD,
    NLOGFLAGS
};


/*! flags for timers */
enum StopCondType{
    GRAD,
    GRADOBJ,
};


struct ReadWriteFlags {
    bool readfiles;           ///< internal flag to indicate that we read files
    bool readvelocity;        ///< internal flag to indicate that we read velocities
    bool timeseries;          ///< write time series to file (debug only; creates a lot of output)
    bool iterates;            ///< write iterates to file
    bool defgrad;             ///< write deformation gradient to file
    bool detdefgrad;          ///< write determinant of deformation gradient
    bool velnorm;             ///< write dataset with norfm ov veloity field
    bool residual;            ///< write dataset with residual between reference and transported template (i.e., the mismatch; abs(m1 - mR))
    bool invresidual;         ///< write dataset with residual between reference and transported template (i.e., the mismatch; 1 - abs(m1 - mR))
    bool defmap;              ///< write deformation map y
    bool templateim;          ///< write template image (original dataset)
    bool referenceim;         ///< write reference image (original dataset)
    bool deftemplate;         ///< write deformed/transported template
    bool deffield;            ///< write deformation field (displacement field)
    bool invdeffield;            ///< write deformation field (displacement field)
    bool velocity;            ///< write velocity field
};


struct FileNames {
    std::string iv1;                    ///< filename for vector field component x1
    std::string iv2;                    ///< filename for vector field component x2
    std::string iv3;                    ///< filename for vector field component x3
    std::string xv1;                    ///< filename for vector field component x1
    std::string xv2;                    ///< filename for vector field component x2
    std::string xv3;                    ///< filename for vector field component x3
    std::string xfolder;                ///< identifier for folder to write results to
    std::string ifolder;                ///< identifier for folder to read in results from
    std::vector < std::string > mt;     ///< template image file name
    std::vector < std::string > mr;     ///< reference image file name
    std::string mask;                   ///< mask for objective
    std::string isc;                    ///< filename for input scalar field
    std::string iref;                   ///< filename for reference input scalar field
    std::vector < std::string > ivec;   ///< filename for input vector image field
    std::string xsc;                    ///< filename for output scalar field
    std::vector < std::string > xvec;    ///< filename for output vector image field
    std::string extension;              ///< identifier for file extension
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
    IntType level;              ///< level for multi-grid schemes, 0 for single grid
    MPI_Comm mpicomm;           ///< communicator for accfft
    MPI_Comm rowcomm;           ///< row communicator for accfft
    MPI_Comm colcomm;           ///< column communicator for accfft
};


/* parameters for optimization */
struct Optimization {
    int maxiter;                 ///< maximal number of (outer) iterations
    int iterbound;               ///< maximal number of (outer) iterations; hard upper limit (in case algorithm diverges)
    int miniter;                 ///< minimal number of (outer) iterations (for parameter continuation)
    ScalarType tol[3];           ///< tolerances for optimization
    ScalarType gtolbound;        ///< tolerances for gradient (lower bound; if maxiter is small, this is what we at least would like to achieve)
    OptMeth method;              ///< optimization method
    GlobalMethType glmethod;     ///< method for globalization (line search; trust region; ...)
    GradientType gradtype;       ///< flag for type of gradient (sobolev or ell-2)
    StopCondType stopcond;       ///< stopping conditions
    bool fastpresolve;           ///< flag to switch on fast presolve
    bool fastsolve;              ///< flag to switch on fast (inaccurate) solve
    ScalarType presolvetol[3];   ///< tolerances for presolve
    int presolvemaxit;           ///< maximal iterations for presolve
    bool usezeroinitialguess;    ///< use zero initial guess for computing the initial gradient
    int solutionstatus;

    // debug
    bool derivativecheckenabled; ///< use zero initial guess for computing the initial gradient
};


/* parameters for krylov solver */
struct KrylovMethod {
    int maxiter;                    ///< max number of iterations for krylov solver
    IntType iter;                   ///< current number of iterations for krylov solver
    ScalarType tol[3];              ///< tolerances for krylov method
    FSeqType fseqtype;              ///< forcing sequence type
    std::string name;               ///< name of krylov solver
    ScalarType reltol;              ///< relative tolerance for krylov solver (used to define tolerance in preconditioner)
    ScalarType g0norm;              ///< initial norm of gradient (to normalize stopping condition)
    bool g0normset;                 ///< flag to identify if initial norm of gradient has been set
    KrylovMethodType solver;        ///< flag for krylov solver

    ScalarType pctol[3];            ///< tolerances for krylov method (preconditioner)
    IntType pcmaxit;                ///< tolerances for krylov method (preconditioner)
    ScalarType pctolint[3];         ///< tolerances for solver of preconditioner
    PrecondMeth pctype;             ///< flag for type of preconditioner
    HessianMatVecType matvectype;   ///< flag for the type of hessian matvec
    std::string pcname;             ///< name of preconditioner
    KrylovMethodType pcsolver;      ///< solver for preconditioner
    bool pcsetupdone;               ///< flag to indicate if setup of preconditioner is done
    ScalarType pctolscale;          ///< tolerance scaling for preconditioner; default: 1E-1
    ScalarType pcgridscale;         ///< this is for the two level preconditioner; defines scale for grid size change; default: 2
    bool usepetsceigest;            ///< in cheb method we need to estimate eigenvalues; use petsc implementation
    int reesteigvals;               ///< flag to reestimate eigenvalues every Krylov(i=1)- or Newton(i=2)-iteration (default: 0)
    bool monitorpcsolver;           ///< flag to monitor PC solver
    bool eigvalsestimated;          ///< flag if eigenvalues have already been estimated
    bool checkhesssymmetry;         ///< check symmetry of hessian operator
    ScalarType hessshift;           ///< perturbation to hessian operator
};


/* parameter for parameter continuation (regularization parameter) */
struct ParCont {
    static constexpr ScalarType betavminh1 = 1E-4;    ///< minimal regularization parameter for h1 type norm
    static constexpr ScalarType betavminh2 = 1E-7;    ///< minimal regularization parameter for h2 type norm
//    static constexpr ScalarType betavminh1 = 1E-9;      ///< minimal regularization parameter for h1 type norm
//    static constexpr ScalarType betavminh2 = 1E-9;      ///< minimal regularization parameter for h2 type norm
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
struct Monitor {
    bool detdgradenabled;       ///< flag to monitor jacobian during iterations
    ScalarType detdgradmin;     ///< min value of determinant of deformation gradient det(grad(y))
    ScalarType detdgradmax;     ///< max value of determinant of deformation gradient det(grad(y))
    ScalarType detdgradmean;    ///< mean value of determinant of deformation gradient det(grad(y))
    ScalarType detdgradbound;   ///< lower bound of determinant of deformation gradient det(grad(y))
    ScalarType jval;            ///< value of objective functional
    ScalarType jval0;           ///< initial value of objective functional
    ScalarType jvalold;         ///< value of objective functional at last iteration
    ScalarType dval;            ///< value of distance measure
    ScalarType dval0;           ///< initial value of distance measure
    ScalarType qmval;           ///< value for inner product between q and m
    ScalarType rval;            ///< value of regularization functional
    ScalarType gradnorm;        ///< norm of gradient at current iteration
    ScalarType gradnorm0;       ///< initial value of norm of gradient
    std::string solverstatus;   ///< string to hold solver status (used in coupling)
};


/* regularization operators */
struct RegNorm {
    ScalarType beta[4];  ///< regularization parameter
    RegNormType type;    ///< flag for regularization norm
};


struct Distance {
    DMType type;
    bool reset;
    ScalarType scale;
};

struct FourierTransform {
    Spectral* fft;      ///< spectral operator
    //FFTPlanType* plan;  ///< accfft plan
    IntType nalloc;     ///< size for allocation in fourier domain
    IntType nx[3];      ///< spatial grid cell size
    IntType osize[3];   ///< size of grid in fourier domain for mpi proc
    IntType ostart[3];  ///< start index in fourier domain for mpi proc
    IntType isize[3];  ///< start index in real domain for mpi proc
    IntType istart[3];  ///< start index in real domain for mpi proc
    ScalarType threshold; ///< threshold to cut of frequencies to zero
};

struct DifferentiationScheme {
    DifferentiationType diffReg;
    DifferentiationType diffPDE;
};

struct RegFlags {
    bool applysmoothing;         ///< apply smoothing to images
    bool applyrescaling;         ///< apply rescaling to images (map the intensity range to [0,1])
    bool registerprobmaps;       ///< flag to identify that we are performing a registration of probabilty maps
    bool detdefgradfromdeffield; ///< compute determinant fo deformation gradient from deformation field (displacement field)
    bool invdefgrad;             ///< compute inverse of deformation gradient
    bool checkdefmapsolve;       ///< check
    bool runinversion;           ///< flag to identify if we are running an inversion or not (lower memory footprint for fwd solve if not)
    bool runsynprob;             ///< true if we run a synthetic test problem
    int synprobid;               ///< flag/id for synthetic problem to be solved
    bool zeroinit;               ///< use zero initial guess for synthetic problems
};


struct PDESolver {
    PDESolverType type;
    PDEType pdetype;
    int rkorder;
    int iporder;
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
    std::vector<ScalarType> gradnorm;        ///< gradient norm
    std::vector<int> newtoniterations;       ///< iterations of solver
    std::vector<ScalarType> krylovresidual;  ///< residual of krylov method
    std::vector<int> kryloviterations;       ///< iterations of krylov method
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


struct MemoryFlags {
    bool savestategrad;       ///< save gradient of state variable for all timesteps
    bool savestategradx;       ///< save interpolated gradient of state variable for all timesteps
};


class RegOpt {
 public:
    typedef RegOpt Self;

    RegOpt();
    RegOpt(int, char**);
    RegOpt(const std::vector<std::string>& args);
    RegOpt(const RegOpt&);
    virtual ~RegOpt();
    void Copy(const RegOpt&);

    inline ScalarType GetLebesgueMeasure(void) {
        return  this->m_Domain.hx[0]
               *this->m_Domain.hx[1]
               *this->m_Domain.hx[2];
    }

    PetscErrorCode SetParameter(const std::vector<std::string>& args);
  
    ScalarType GetBetaMinParaCont();

    PetscErrorCode GetSizes(IntType*, IntType&, IntType&);
    PetscErrorCode GetSizes(IntType*, IntType*, IntType*);

    inline ScalarType GetTimeStepSize(void) {
        return (this->m_Domain.timehorizon[1] - this->m_Domain.timehorizon[0])
               /static_cast<ScalarType>(this->m_Domain.nt);
    }

    /* do setup for grid continuation */
    PetscErrorCode SetupGridCont();

    // timers and counters
    inline unsigned int GetCounter(CounterType id) {return this->m_Counter[id];}
    inline void IncrementCounter(CounterType id) {this->m_Counter[id]++;}
    inline void IncrementCounter(const CounterType id, const int i) {this->m_Counter[id] += i;}

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

    inline void LogKSPResidual(const int i, const ScalarType value){
        this->m_Log.krylovresidual.push_back(value);
        this->m_Log.kryloviterations.push_back(i);
    }
    inline void LogConvergence(const int i, const ScalarType J, const ScalarType D, const ScalarType R){
        this->m_Log.distance.push_back(D);
        this->m_Log.regularization.push_back(R);
        this->m_Log.objective.push_back(J);
        this->m_Log.newtoniterations.push_back(i);
    }
    inline void LogFinalResidual(const int i, const ScalarType value){
        this->m_Log.finalresidual[i] = value;
    }

    ScalarType ComputeFFTScale();

    PetscErrorCode StartTimer(TimerType);
    PetscErrorCode StopTimer(TimerType);
    PetscErrorCode ResetTimers(void);
    PetscErrorCode ResetTimer(TimerType);
    PetscErrorCode ResetCounters(void);
    PetscErrorCode ResetCounter(CounterType);
    PetscErrorCode ProcessTimers(void);

    PetscScalar m_GPUtime = 0;
    PetscScalar m_CPUtime = 0;

    virtual PetscErrorCode DisplayOptions(void);
    PetscErrorCode DisplayTimeToSolution(void);
    PetscErrorCode GetTimeToSolution(double&);
    PetscErrorCode WriteLogFile(bool coarse = false);
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
    PetscErrorCode ResetDM(DMType type);

    RegModel m_RegModel {};              ///< flag for particular registration model
    Domain m_Domain {};                  ///< parameters for spatial domain
    GridCont m_GridCont {};              ///< flags for grid continuation
    ScaleCont m_ScaleCont {};            ///< flags for scale continuation
    Optimization m_OptPara {};           ///< optimization parameters
    ReadWriteFlags m_ReadWriteFlags {};  ///< flags for input/output of fields (defmap, deffield, defgrad, ...)
    PDESolver m_PDESolver {};            ///< flag for PDE solver
    KrylovMethod m_KrylovMethod {};      ///< parameters for krylov solver
    RegFlags m_RegFlags {};              ///< flags for registration
    Monitor m_Monitor {};                ///< monitor for registration
    RegNorm m_RegNorm {};                ///< parameters for regularization model
    Distance m_Distance {};              ///< parameters for distance measure
    FourierTransform m_FFT {};           ///< parameters for FFT/accfft
    FourierTransform m_FFT_coarse {};    ///< parameters for FFT/accfft
    DifferentiationScheme m_Diff {};     ///< differentiation schemes
    ParCont m_ParaCont {};               ///< flags for parameter continuation
    SolveType m_SolveType {};            ///< solver
    FileNames m_FileNames {};            ///< file names for input/output
    Logger m_Log {};                     ///< log
    MemoryFlags m_Mem {};                ///< memory options
    ScalarType m_Sigma[3];               ///< standard deviation for gaussian smoothing

    bool m_SetupDone;
    bool m_StoreCheckPoints;

    int m_CartGridDims[2];
    unsigned int m_NumThreads;
    unsigned int m_LineLength;
//    const unsigned int m_LineLength = 101; //C++ 11 feature
    unsigned int m_Indent;
    int m_Verbosity;
    std::vector<int> m_LabelIDs;       ///< label ids
    std::vector<ScalarType> m_ObjWts;  ///< Objective function component weights
    std::string m_PostFix;
    
#ifdef REG_HAS_CUDA
    int m_gpu_id;                      ///< id of used GPU
#endif
    int rank;                          ///< process id 
    int rank_cnt;                      ///< number of processes


 protected:
    virtual PetscErrorCode Initialize(void);
    PetscErrorCode InitializeFFT();
    PetscErrorCode DestroyFFT();
    virtual PetscErrorCode ClearMemory(void);
    virtual PetscErrorCode ParseArguments(const std::vector<std::string>& args);
    virtual PetscErrorCode ParseArguments(int, char**);
    virtual PetscErrorCode Usage(bool advanced = false);
    virtual PetscErrorCode CheckArguments(void);
    //virtual PetscErrorCode SetupParser(void);
    PetscErrorCode SetPresetParameters();
    PetscErrorCode WriteWorkLoadLog();
    PetscErrorCode WriteWorkLoadLogReadable();
    PetscErrorCode WriteWorkLoadLog(std::ostream&);
    PetscErrorCode WriteWorkLoadLogReadable(std::ostream&);
    PetscErrorCode WriteKSPLog();
    PetscErrorCode WriteConvergenceLog();
    PetscErrorCode WriteFinalResidualLog();

    enum TimerValue {LOG = 0, MIN, MAX, AVG, NVALTYPES};

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
    
    //ARGVParser m_Parser;
};




}  // namespace reg




#endif  // _REGOPT_H_
