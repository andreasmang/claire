
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
  N_MISC(int* iN,int* iisize, int* iistart,accfft_plan* iplan,MPI_Comm ic_comm){
    memcpy(this->N,iN,3*sizeof(int));
    memcpy(this->isize,iisize,3*sizeof(int));
    memcpy(this->istart,iistart,3*sizeof(int));
//    this->Np=Np;
    this->plan=iplan;
    this->c_comm=ic_comm;
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


#endif
