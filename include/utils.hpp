
#ifndef _UTILS_
#define _UTILS_

//#define BRAIN
#define ALIGNMENT 32
#define OMP_NUM_THREADS 16
#include <petsc.h>
#include <accfft.h>
#include <accfft_operators.h>
enum SolveType{Euler_Explicit,RK2_Explicit,RK2_Implicit};
class N_MISC{
  public:
	int N[3];
	int isize[3];
	int istart[3];
  double h[3];
	double dt;
	int Nt;
	double T;
  int Np; // number of initial condition parametrization
  long long int N_global;
  long long int N_local;
  accfft_plan* plan;
  MPI_Comm c_comm;
  // Constructor
  N_MISC(int* N,int* isize, int* istart,accfft_plan* plan,MPI_Comm c_comm){
    memcpy(this->N,N,3*sizeof(int));
    memcpy(this->isize,isize,3*sizeof(int));
    memcpy(this->istart,istart,3*sizeof(int));
    this->Np=Np;
    this->plan=plan;
    this->c_comm=c_comm;
    Nt=NULL;
    dt=NULL;
    T=NULL;
    Np=NULL;
    N_global=static_cast<long long int>(N[0])*static_cast<long long int>(N[1])*static_cast<long long int>(N[2]);
    N_local=static_cast<long long int>(isize[0]*isize[1]*isize[2]);
    h[0]=M_PI*2/N[0];
    h[1]=M_PI*2/N[1];
    h[2]=M_PI*2/N[2];

  };
  N_MISC(){
    this->N[0]=NULL;
    this->N[1]=NULL;
    this->N[2]=NULL;
    this->isize[0]=NULL;
    this->isize[1]=NULL;
    this->isize[2]=NULL;
    this->istart[0]=NULL;
    this->istart[1]=NULL;
    this->istart[2]=NULL;
    this->Np=NULL;
    this->plan=NULL;
    this->c_comm=NULL;
    this->Nt=NULL;
    this->dt=NULL;
    this->T=NULL;
    this->Np=NULL;
    this->N_global=NULL;
    this->N_local=NULL;
    h[0]=NULL;
    h[1]=NULL;
    h[2]=NULL;
    return ;

  };

  void Print(){
    int nprocs, procid;
    MPI_Comm_rank(c_comm, &procid);
    MPI_Comm_size(c_comm, &nprocs);
    PCOUT<<"N[0]="<<N[0]<<" N[1]="<<N[1]<<" N[2]="<<N[2]<<std::endl;
    PCOUT<<"isize[0]="<<isize[0]<<" isize[1]="<<isize[1]<<" isize[2]="<<isize[2]<<std::endl;
    PCOUT<<"istart[0]="<<istart[0]<<" istart[1]="<<istart[1]<<" istart[2]="<<istart[2]<<std::endl;
    PCOUT<<"dt="<<dt<<std::endl;
    PCOUT<<"Nt="<<Nt<<std::endl;
    PCOUT<<"T="<<T<<std::endl;
    PCOUT<<"Np="<<Np<<std::endl;
    PCOUT<<"N_global="<<N_global<<std::endl;
    PCOUT<<"N_local="<<N_local<<std::endl;
  }
  // copy assignment
  N_MISC& operator=(N_MISC other){
    memcpy(this->N,other.N,3*sizeof(int));
    memcpy(this->isize,other.isize,3*sizeof(int));
    memcpy(this->istart,other.istart,3*sizeof(int));
    this->Np=other.Np;
    this->plan=other.plan;
    this->c_comm=other.c_comm;
    this->Nt=other.Nt;
    this->dt=other.dt;
    this->T=other.T;
    this->Np=other.Np;
    this->N_global=other.N_global;
    this->N_local=other.N_local;
    return *this;
  };
};

class DiffCoef{
  public:
    // Constructor
    DiffCoef(N_MISC N_Misc);
    MPI_Comm c_comm;
    double k_a; // Average Diff coeff
    Vec kxx;
    Vec kxy;
    Vec kxz;
    Vec kyy;
    Vec kyz;
    Vec kzz;
    Vec kfxx;
    Vec kfxy;
    Vec kfxz;
    Vec kfyy;
    Vec kfyz;
    Vec kfzz;

    // Misc Data
    Vec White_Matter;
    Vec Gray_Matter;
    Vec Glial_Matter;
    Vec Filter;
    Vec IsoDif;
    double kxx_a;
    double kxy_a;
    double kxz_a;
    double kyy_a;
    double kyz_a;
    double kzz_a;

    PetscErrorCode SetValues(double scale=0,double scale_f=0);
    //PetscErrorCode SetValues(double scale=0);
    PetscErrorCode Smooth();
    PetscErrorCode Print();
    PetscErrorCode ZeroFilter(double * c);
    PetscErrorCode ApplyTensor(Vec Tx,Vec Ty,Vec Tz, Vec x, Vec y, Vec z) const;
    // Deconstructor
    ~DiffCoef();
  private:
    long long int N_global,N_local;
    int N[3];
    N_MISC N_Misc;
};
class ReacCoef{
  public:
    // Constructor
    ReacCoef(N_MISC N);

    // Members
    Vec rho;
    double r_a;
    double rho_scale;
    bool linear;
    MPI_Comm c_comm;

    PetscErrorCode SetValues(double scale=0);
    PetscErrorCode Rescale(double rescale);
    PetscErrorCode Print();
    PetscErrorCode Smooth();
    // Deconstructor
    ~ReacCoef();

  private:
    long long int N_global,N_local;
    int N[3];
    N_MISC N_Misc;
};
class Phi{
  public:
    // Constructor
    Phi(N_MISC N);

    // Members
    Vec * phi; // Vec ptr to Nk Vectors
    PetscErrorCode SetValues();
    PetscErrorCode Apply(double *c_out, double *p_in);
    PetscErrorCode Apply(double *c_out, Vec& p_in);
    PetscErrorCode ApplyT(double *p_out, double *c_in);
    MPI_Comm c_comm;

    // Deconstructor
    ~Phi();

    int p_lsize;
    int p_gindex;
  private:
    long long int N_global,N_local;
    void Initialize(double *out,const int * N,double sigma, double *center);
    PetscErrorCode phimesh(double *center, N_MISC N_Misc, double *cm);
    int N[3];
    N_MISC N_Misc;
    int *p_lsizes;
};
class Obs{
  public:
    // Constructor
    Obs(N_MISC N,double * data,double obs_thr);

    // Members
    double Threshold;
    Vec  Filter;
    PetscErrorCode Apply (double * y, double * x);
    PetscErrorCode ApplyT(double * y, double * x);

    // Deconstructor
    ~Obs();

  private:
    long long int N_global,N_local;
    int N[3];
    N_MISC N_Misc;
    MPI_Comm c_comm;
};
class Time_History{
  public:
    // Constructor
    Time_History(N_MISC N);
    //Time_History& Time_History::operator=( const Time_History& other );

    // Members
    MPI_Comm c_comm;
    Vec * C;
    Vec * C_tilde;
    Vec * Alpha;
    Vec * Alpha_tilde;
    bool write_flag; // if 1 data will be written to C and Alpha
    PetscErrorCode CopyC(int time_counter, double * c_in, int linearized);
    PetscErrorCode CopyA(int time_counter, double * c_in, int linearized);
    PetscErrorCode CopyA_tilde(int time_counter, double * c_in, int linearized);
    PetscErrorCode PrintOut();
    PetscErrorCode PrintNorms();
    PetscErrorCode Reset();
    PetscErrorCode Zero_Adjoint();

    // Deconstructor
    ~Time_History();

  private:
    long long int N_global,N_local;
    int N[3];
    N_MISC N_Misc;
    int Nt;
};
const char* commandline_option(int argc, char** argv, const char* opt, const char* def_val, bool required, const char* err_msg);
void DataOut(double * A, N_MISC N_Misc,const char * fname);
void my_daxpy(const int n,const double alpha,double *x,double *y);
void my_daxpy(const int n,const double alpha,double *x);
void my_daxpy(const int *n,const double*  alpha, double *x,const int* stridex, double *y, const int* stride);
int read_binary(double * data, N_MISC N, const char *,const char* ="./brain_data/64/");

#endif
