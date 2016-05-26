// debugging with ib 4 bin/step1 16 16 16

#include <stdlib.h>
#include <math.h> // M_PI
#include <mpi.h>
#include <accfft.h>
#include "interp3.hpp"

#include <vector>
//#define VERBOSE2
#define NREP 10



size_t accfft_ghost_local_size_dft_r2c(accfft_plan* plan,int g_size, int * isize_g, int* istart_g);
void accfft_get_ghost(accfft_plan* plan,int g_size,int* isize_g, Real* data,Real* ghost_data);

size_t accfft_ghost_xyz_local_size_dft_r2c(accfft_plan* plan,int g_size, int * isize_g, int* istart_g);
void accfft_get_ghost_xyz(accfft_plan* plan,int g_size,int* isize_g, Real* data,Real* ghost_data);

template <class Real>
Real fn(Real x, Real y, Real z){
  Real a=80;

  return 2+sin(2*M_PI*x)*sin(4*M_PI*y)*cos(2*M_PI*z);
  // if you want to add exponential you have to explicity put an if condition
  // and enforce periodicity so that for example X=-0.5 gives the same result as X=0.5.
  // Otherwise if you naively call this function for error checking X=-0.5 would not
  // give the right result
  //+exp(-a*(x-0.5)*(x-0.5)-2*a*(y-0.5)*(y-0.5)-3*a*(z-.5)*(z-.5));

}


// Ghost Region (x,y,z padding), Row Major, General N,
void test1(int * N){

  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  PCOUT<<"\n\n****************************\n";
  PCOUT<<"    Performing test1        \n";
  PCOUT<<"****************************\n";

  /* Create Cartesian Communicator */
  int c_dims[2]={0};
  MPI_Comm c_comm;
  accfft_create_comm(MPI_COMM_WORLD,c_dims,&c_comm);

  /* Get the local pencil size and the allocation size */
  int isize[3],osize[3],istart[3],ostart[3];
  int alloc_max=accfft_local_size_dft_r2c(N,isize,istart,osize,ostart,c_comm);
  size_t N_local=isize[0]*isize[1]*isize[2];

  Real * data=(Real*)accfft_alloc(alloc_max);
  Complex* data_hat=(Complex*)accfft_alloc(alloc_max);

  accfft_init();

  /* Create FFT plan */
  accfft_plan * plan=accfft_plan_dft_3d_r2c(N,data,(Real*)data_hat,c_comm,ACCFFT_MEASURE);


  Real *f; // The function at regular grid values
  Real *query_points; // The query points (linearized array with xyz coordinates of every point entered back to back)
  Real *f_cubic; // Interpolated values corresponding to query points
  Real *f_ref; // reference value of f at query points (used only to compute error)

  long Q=1;//isize[0]*isize[1]; // The number of query points

  f=(Real*) malloc(N_local*sizeof(Real));
  f_cubic=(Real*) malloc(Q*sizeof(Real));
  f_ref=(Real*) malloc(Q*sizeof(Real));
  query_points=(Real*) malloc(Q*sizeof(Real)*COORD_DIM); // Note that you have to provide 3 coordinates back to back thus the factor coord_dim


  // initialize f
  for(long i0=0;i0<isize[0];i0++)
  for(long i1=0;i1<isize[1];i1++)
  for(long i2=0;i2<isize[2];i2++){
    Real c[COORD_DIM]={(i0+istart[0])*1.0/N[0],
                       (i1+istart[1])*1.0/N[1],
                       (i2+istart[2])*1.0/N[2]};
    f[i2+isize[2]*(i1+isize[1]*i0)]=fn(c[0],c[1],c[2]);
  }

  // Get ghost cells of f, and store the resutls into ghost_f
  int g_size=3;
  int isize_g[3],istart_g[3];
  size_t g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c( plan,g_size,isize_g,istart_g);
  Real *ghost_f=(Real*)accfft_alloc(g_alloc_max);
  accfft_get_ghost_xyz(plan,g_size,isize_g,f,ghost_f);
  //PCOUT<<"ghost_f[0]="<<ghost_f[0]<<std::endl;
  //PCOUT<<"f[0]="<<f[0]<<std::endl;


  // Create query points
  Real h[3];
  h[0]=1./N[0];
  h[1]=1./N[1];
  h[2]=1./N[2];
  int ref=0;
  query_points[ref*COORD_DIM+0]=0.2;
  query_points[ref*COORD_DIM+1]=0.2;
  query_points[ref*COORD_DIM+2]=0.2;
  f_ref[ref]=fn(query_points[ref*COORD_DIM+0],query_points[ref*COORD_DIM+1],query_points[ref*COORD_DIM+2]);

  // First create random query points everywhere in the global domain
  // Note that the points can be outside the locally owned grid
  for(long i=1;i<Q;i++){
    query_points[i*COORD_DIM+0]=2*drand48()-0.6;
    query_points[i*COORD_DIM+1]=2*drand48()-0.6;
    query_points[i*COORD_DIM+2]=2*drand48()-0.6;
    f_ref[i]=fn(query_points[i*COORD_DIM+0],query_points[i*COORD_DIM+1],query_points[i*COORD_DIM+2]);
  }

  // Perform parallell interpolation. Note that we are passing ghost_f instead of f.
  par_interp3_ghost_xyz_p(ghost_f, 1, N, isize,istart,Q,g_size,query_points, f_cubic,c_dims,c_comm);

  // Compute the global error
  double err=0,norm=0;
  for(long i=0;i<Q;i++){
    err=std::max(err,fabs(f_ref[i]-f_cubic[i]));
    norm+=fabs(f_ref[i]);
    if(err>=0.941512){
      std::cout<<procid<<" i= "<<i<<std::endl;
      break;}
  }
  //PCOUT<<"err="<<g_err<<" Rel. Err.= "<<g_err/g_norm<<std::endl;

  if(procid==0){
    int i=0;
    std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
    std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
    std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
    std::cout<<"norm="<<norm<<std::endl;
    std::cout<<"f_true= "<<f_ref[i]<<" f_cubic= "<<f_cubic[i]<<" err="<<fabs(f_cubic[i]-f_ref[i])/f_ref[i]<<'\n';
    std::cout<<"L_inf Error in interpolation of "<<Q<<" points is "<<err<<'\n';
    std::cout<<"Rel. L_inf Error in interpolation of "<<Q<<" points is "<<err/norm<<'\n';
  }
  Real g_err,g_norm;
  MPI_Reduce(&err,&g_err,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm,&g_norm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  PCOUT<<"\n Global Error:\n";
  PCOUT<<"L_inf Error in interpolation of "<<Q<<" points is "<<g_err<<'\n';
  PCOUT<<"Rel. L_inf Error in interpolation of "<<Q<<" points is "<<g_err/g_norm<<'\n';


  accfft_free(data);
  accfft_free(data_hat);
  accfft_free(ghost_f);
  free(f);
  free(f_ref);
  free(f_cubic);
  free(query_points);

}

// Same as test1 but with three fields interpolated seperately
void test2(int * N){

  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  PCOUT<<"\n\n****************************\n";
  PCOUT<<"    Performing test2        \n";
  PCOUT<<"****************************\n";

  /* Create Cartesian Communicator */
  int c_dims[2]={0};
  MPI_Comm c_comm;
  accfft_create_comm(MPI_COMM_WORLD,c_dims,&c_comm);

  /* Get the local pencil size and the allocation size */
  int isize[3],osize[3],istart[3],ostart[3];
  int alloc_max=accfft_local_size_dft_r2c(N,isize,istart,osize,ostart,c_comm);
  size_t N_local=isize[0]*isize[1]*isize[2];

  Real * data=(Real*)accfft_alloc(alloc_max);
  Complex* data_hat=(Complex*)accfft_alloc(alloc_max);

  accfft_init();

  /* Create FFT plan */
  accfft_plan * plan=accfft_plan_dft_3d_r2c(N,data,(Real*)data_hat,c_comm,ACCFFT_MEASURE);


  Real *f_x,*f_y,*f_z; // The function at regular grid values
  Real *query_points; // The query points (linearized array with xyz coordinates of every point entered back to back)
  Real *f_cubic_x,*f_cubic_y,*f_cubic_z; // Interpolated values corresponding to query points
  Real *f_ref; // reference value of f at query points (used only to compute error)

  long Q=16384;//isize[0]*isize[1]; // The number of query points

  f_x=(Real*) malloc(N_local*sizeof(Real));
  f_y=(Real*) malloc(N_local*sizeof(Real));
  f_z=(Real*) malloc(N_local*sizeof(Real));

  f_cubic_x=(Real*) malloc(Q*sizeof(Real));
  f_cubic_y=(Real*) malloc(Q*sizeof(Real));
  f_cubic_z=(Real*) malloc(Q*sizeof(Real));
  f_ref=(Real*) malloc(Q*sizeof(Real));
  query_points=(Real*) malloc(Q*sizeof(Real)*COORD_DIM); // Note that you have to provide 3 coordinates back to back thus the factor coord_dim


  // initialize f
  for(long i0=0;i0<isize[0];i0++)
    for(long i1=0;i1<isize[1];i1++)
      for(long i2=0;i2<isize[2];i2++){
        Real c[COORD_DIM]={(i0+istart[0])*1.0/N[0],
          (i1+istart[1])*1.0/N[1],
          (i2+istart[2])*1.0/N[2]};
        f_x[i2+isize[2]*(i1+isize[1]*i0)]=fn(c[0],c[1],c[2]);
        f_y[i2+isize[2]*(i1+isize[1]*i0)]=fn(c[0],c[1],c[2]);
        f_z[i2+isize[2]*(i1+isize[1]*i0)]=fn(c[0],c[1],c[2]);
      }

  // Get ghost cells of f, and store the resutls into ghost_f
  int g_size=3;
  int isize_g[3],istart_g[3];
  size_t g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c( plan,g_size,isize_g,istart_g);
  Real *ghost_f_x=(Real*)accfft_alloc(g_alloc_max);
  Real *ghost_f_y=(Real*)accfft_alloc(g_alloc_max);
  Real *ghost_f_z=(Real*)accfft_alloc(g_alloc_max);
  accfft_get_ghost_xyz(plan,g_size,isize_g,f_x,ghost_f_x);
  accfft_get_ghost_xyz(plan,g_size,isize_g,f_y,ghost_f_y);
  accfft_get_ghost_xyz(plan,g_size,isize_g,f_z,ghost_f_z);
  //PCOUT<<"ghost_f[0]="<<ghost_f[0]<<std::endl;
  //PCOUT<<"f[0]="<<f[0]<<std::endl;


  // Create query points
  Real h[3];
  h[0]=1./N[0];
  h[1]=1./N[1];
  h[2]=1./N[2];
  int ref;

  // First create random query points everywhere in the global domain
  // Note that the points can be outside the locally owned grid
  for(long i=0;i<Q;i++){
    query_points[i*COORD_DIM+0]=drand48();
    query_points[i*COORD_DIM+1]=drand48();
    query_points[i*COORD_DIM+2]=drand48();
    f_ref[i]=fn(query_points[i*COORD_DIM+0],query_points[i*COORD_DIM+1],query_points[i*COORD_DIM+2]);
  }

  // Perform parallell interpolation. Note that we are passing ghost_f instead of f.
  for(int rep=0;rep<10;rep++){
    par_interp3_ghost_xyz_p(ghost_f_x, 1,N,isize,istart,Q,g_size,query_points, f_cubic_x,c_dims,c_comm);
    par_interp3_ghost_xyz_p(ghost_f_y, 1,N,isize,istart,Q,g_size,query_points, f_cubic_y,c_dims,c_comm);
    par_interp3_ghost_xyz_p(ghost_f_z, 1,N,isize,istart,Q,g_size,query_points, f_cubic_z,c_dims,c_comm);
  }
  // Compute the global error
  double err=0,norm=0;
  for(long i=0;i<Q;i++){
    err=std::max(err,fabs(f_ref[i]-f_cubic_x[i]));
    err=std::max(err,fabs(f_ref[i]-f_cubic_y[i]));
    err=std::max(err,fabs(f_ref[i]-f_cubic_z[i]));
    norm+=fabs(f_ref[i]);
  }
  //PCOUT<<"err="<<g_err<<" Rel. Err.= "<<g_err/g_norm<<std::endl;

  if(procid==0){
    int i=0;
    std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
    std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
    std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
    std::cout<<"norm="<<norm<<std::endl;
    std::cout<<"f_true= "<<f_ref[i]<<" f_cubic= "<<f_cubic_x[i]<<" err="<<fabs(f_cubic_x[i]-f_ref[i])/f_ref[i]<<'\n';
    std::cout<<"L_inf Error in interpolation of "<<Q<<" points is "<<err<<'\n';
    std::cout<<"Rel. L_inf Error in interpolation of "<<Q<<" points is "<<err/norm<<'\n';
  }
  Real g_err,g_norm;
  MPI_Reduce(&err,&g_err,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm,&g_norm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  PCOUT<<"\n Global Error:\n";
  PCOUT<<"L_inf Error in interpolation of "<<Q<<" points is "<<g_err<<'\n';
  PCOUT<<"Rel. L_inf Error in interpolation of "<<Q<<" points is "<<g_err/g_norm<<'\n';


  accfft_free(data);
  accfft_free(data_hat);
  accfft_free(ghost_f_x);
  accfft_free(ghost_f_y);
  accfft_free(ghost_f_z);
  free(f_x);
  free(f_y);
  free(f_z);
  free(f_ref);
  free(f_cubic_x);
  free(f_cubic_y);
  free(f_cubic_z);
  free(query_points);

}

// Same as test1 but with three fields interpolated in one shot
void test3(int * N){

  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  PCOUT<<"\n\n****************************\n";
  PCOUT<<"    Performing test3        \n";
  PCOUT<<"****************************\n";

  /* Create Cartesian Communicator */
  int c_dims[2]={0};
  MPI_Comm c_comm;
  accfft_create_comm(MPI_COMM_WORLD,c_dims,&c_comm);

  /* Get the local pencil size and the allocation size */
  int isize[3],osize[3],istart[3],ostart[3];
  int alloc_max=accfft_local_size_dft_r2c(N,isize,istart,osize,ostart,c_comm);
  size_t N_local=isize[0]*isize[1]*isize[2];

  Real * data=(Real*)accfft_alloc(alloc_max);
  Complex* data_hat=(Complex*)accfft_alloc(alloc_max);

  accfft_init();

  /* Create FFT plan */
  accfft_plan * plan=accfft_plan_dft_3d_r2c(N,data,(Real*)data_hat,c_comm,ACCFFT_MEASURE);


  Real *f; // The function at regular grid values
  Real *query_points; // The query points (linearized array with xyz coordinates of every point entered back to back)
  Real *f_cubic; // Interpolated values corresponding to query points
  Real *f_ref; // reference value of f at query points (used only to compute error)

  long Q=isize[0]*isize[1]; // The number of query points
  int data_dof=3; // Three velocity fields

  f=(Real*) malloc(N_local*sizeof(Real)*data_dof);
  f_cubic=(Real*) malloc(Q*sizeof(Real)*data_dof);
  f_ref=(Real*) malloc(Q*sizeof(Real));
  query_points=(Real*) malloc(Q*sizeof(Real)*COORD_DIM); // Note that you have to provide 3 coordinates back to back thus the factor coord_dim


  // initialize f
  for(long i0=0;i0<isize[0];i0++)
  for(long i1=0;i1<isize[1];i1++)
  for(long i2=0;i2<isize[2];i2++){
    Real c[COORD_DIM]={(i0+istart[0])*1.0/N[0],
                       (i1+istart[1])*1.0/N[1],
                       (i2+istart[2])*1.0/N[2]};
    for(int dof=0;dof<data_dof;++dof){
      f[i2+N[2]*(i1+isize[1]*i0)+dof*N_local]=(dof+1)*fn(c[0],c[1],c[2]);
    }
  }


  // Get ghost cells of f, and store the resutls into ghost_f
  int g_size=3;
  int isize_g[3],istart_g[3];
  size_t g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c( plan,g_size,isize_g,istart_g);
  Real *ghost_f=(Real*)accfft_alloc(g_alloc_max*data_dof);
  for(int dof=0;dof<data_dof;++dof){
    accfft_get_ghost_xyz(plan,g_size,isize_g,&f[dof*isize[0]*isize[1]*isize[2]],&ghost_f[dof*isize_g[0]*isize_g[1]*isize_g[2]]);
  }


  // Create query points
  Real h[3];
  h[0]=1./N[0];
  h[1]=1./N[1];
  h[2]=1./N[2];
  int ref=0;
  //query_points[ref*COORD_DIM+0]=-0.490114;
  //query_points[ref*COORD_DIM+1]=-0.490787;
  //query_points[ref*COORD_DIM+2]=-0.520321;
  //f_ref[ref]=fn(query_points[ref*COORD_DIM+0],query_points[ref*COORD_DIM+1],query_points[ref*COORD_DIM+2]);

  // First create random query points everywhere in the global domain
  // Note that the points can be outside the locally owned grid
  for(long i=0;i<Q;i++){
    query_points[i*COORD_DIM+0]=2*drand48()-0.6;
    query_points[i*COORD_DIM+1]=2*drand48()-0.6;
    query_points[i*COORD_DIM+2]=2*drand48()-0.6;
    f_ref[i]=fn(query_points[i*COORD_DIM+0],query_points[i*COORD_DIM+1],query_points[i*COORD_DIM+2]);
  }

  // Perform parallell interpolation. Note that we are passing ghost_f instead of f.
  double time=0;
  for(int rep=0;rep<NREP;rep++){
    time+=-MPI_Wtime();
    par_interp3_ghost_xyz_p(ghost_f, data_dof, N, isize,istart,Q,g_size,query_points, f_cubic,c_dims,c_comm);
    time+=MPI_Wtime();
    MPI_Barrier(c_comm);
  }

  // Compute the global error
  double err=0,norm=0;
  for(long i=0;i<Q;i++){
    for(int dof=0;dof<data_dof;++dof){
      err=std::max(err,fabs((dof+1)*f_ref[i]-f_cubic[i+dof*Q]));
    }
    norm+=fabs(f_ref[i]);
  }

  if(0)
  if(procid==0){
    int i=1;
    std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
    std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
    std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
    std::cout<<"norm="<<norm<<std::endl;
    for(int dof=0;dof<data_dof;++dof){
      std::cout<<"f_true= "<<(dof+1)*f_ref[i]<<" f_cubic= "<<f_cubic[i+dof*Q]<<" err="<<fabs(f_cubic[i+dof*Q]-(dof+1)*f_ref[i])/f_ref[i]<<'\n';
    }
    std::cout<<"L_inf Error in interpolation of "<<Q<<" points is "<<err<<'\n';
    std::cout<<"Rel. L_inf Error in interpolation of "<<Q<<" points is "<<err/norm<<'\n';
  }
  double g_err,g_norm,g_time;
  MPI_Reduce(&err,&g_err,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&time,&g_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm,&g_norm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  PCOUT<<"\n Global Error:\n";
  PCOUT<<"L_inf Error in interpolation of "<<Q<<" points is "<<g_err<<'\n';
  PCOUT<<"Rel. L_inf Error in interpolation of "<<Q<<" points is "<<g_err/g_norm<<'\n';
  PCOUT<<"Time of "<<NREP<<" iterations ="<<g_time<<std::endl;


  accfft_free(data);
  accfft_free(data_hat);
  accfft_free(ghost_f);
  free(f);
  free(f_ref);
  free(f_cubic);
  free(query_points);

}

// Same as test3 but with a planner
void test4(int * N){

  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  PCOUT<<"\n\n****************************\n";
  PCOUT<<"    Performing test4         \n";
  PCOUT<<"****************************\n";

  /* Create Cartesian Communicator */
  int c_dims[2]={0};
  MPI_Comm c_comm;
  accfft_create_comm(MPI_COMM_WORLD,c_dims,&c_comm);

  /* Get the local pencil size and the allocation size */
  int isize[3],osize[3],istart[3],ostart[3];
  int alloc_max=accfft_local_size_dft_r2c(N,isize,istart,osize,ostart,c_comm);
  size_t N_local=isize[0]*isize[1]*isize[2];

  Real * data=(Real*)accfft_alloc(alloc_max);
  Complex* data_hat=(Complex*)accfft_alloc(alloc_max);

  accfft_init();

  /* Create FFT plan */
  accfft_plan * plan=accfft_plan_dft_3d_r2c(N,data,(Real*)data_hat,c_comm,ACCFFT_MEASURE);


  Real *f; // The function at regular grid values
  Real *query_points; // The query points (linearized array with xyz coordinates of every point entered back to back)
  Real *f_cubic; // Interpolated values corresponding to query points
  Real *f_ref; // reference value of f at query points (used only to compute error)

  long Q=isize[0]*isize[1]; // The number of query points
  int data_dof=3; // Three velocity fields

  f=(Real*) malloc(N_local*sizeof(Real)*data_dof);
  f_cubic=(Real*) malloc(Q*sizeof(Real)*data_dof);
  f_ref=(Real*) malloc(Q*sizeof(Real));
  query_points=(Real*) malloc(Q*sizeof(Real)*COORD_DIM); // Note that you have to provide 3 coordinates back to back thus the factor coord_dim


  // initialize f
  for(long i0=0;i0<isize[0];i0++)
    for(long i1=0;i1<isize[1];i1++)
      for(long i2=0;i2<isize[2];i2++){
        Real c[COORD_DIM]={(i0+istart[0])*1.0/N[0],
          (i1+istart[1])*1.0/N[1],
          (i2+istart[2])*1.0/N[2]};
        for(int dof=0;dof<data_dof;++dof){
          f[i2+N[2]*(i1+isize[1]*i0)+dof*N_local]=(dof+1)*fn(c[0],c[1],c[2]);
        }
      }


  // Get ghost cells of f, and store the resutls into ghost_f
  int g_size=3;
  int isize_g[3],istart_g[3];
  size_t g_alloc_max=accfft_ghost_xyz_local_size_dft_r2c( plan,g_size,isize_g,istart_g);
  Real *ghost_f=(Real*)accfft_alloc(g_alloc_max*data_dof);
  for(int dof=0;dof<data_dof;++dof){
    accfft_get_ghost_xyz(plan,g_size,isize_g,&f[dof*isize[0]*isize[1]*isize[2]],&ghost_f[dof*isize_g[0]*isize_g[1]*isize_g[2]]);
  }


  // Create query points
  Real h[3];
  h[0]=1./N[0];
  h[1]=1./N[1];
  h[2]=1./N[2];
  int ref=0;

  // First create random query points everywhere in the global domain
  // Note that the points can be outside the locally owned grid
  for(long i=0;i<Q;i++){
    query_points[i*COORD_DIM+0]=2*drand48()-0.6;
    query_points[i*COORD_DIM+1]=2*drand48()-0.6;
    query_points[i*COORD_DIM+2]=2*drand48()-0.6;
    f_ref[i]=fn(query_points[i*COORD_DIM+0],query_points[i*COORD_DIM+1],query_points[i*COORD_DIM+2]);
  }

  // Perform parallell interpolation. Note that we are passing ghost_f instead of f.
  double time=0;
  double interpolation_timer[4]={0};


  //Real* dummy_query_points=(Real*) malloc(Q*sizeof(Real)*COORD_DIM); // Note that you have to provide 3 coordinates back to back thus the factor coord_dim
  // let's include the setup cost
  //time+=-MPI_Wtime();
  Interp3_Plan interp3_plan;
  interp3_plan.allocate(Q,data_dof);
  //interp3_plan.scatter(data_dof, N, isize,istart,Q,g_size,dummy_query_points,c_dims,c_comm,interpolation_timer);// has to be called only once
  interp3_plan.scatter(data_dof, N, isize,istart,Q,g_size,query_points,c_dims,c_comm,interpolation_timer);// has to be called only once
  //time+=MPI_Wtime();

  for(int rep=0;rep<NREP;rep++){
    time+=-MPI_Wtime();
    //interp3_plan.slow_run(ghost_f, data_dof, N, isize,istart,Q,g_size,query_points, f_cubic,c_dims,c_comm);
    interp3_plan.interpolate(ghost_f, data_dof, N, isize,istart,Q,g_size, f_cubic,c_dims,c_comm,interpolation_timer);
    time+=MPI_Wtime();
    MPI_Barrier(c_comm);
    //PCOUT<<"rep="<<rep<<std::endl;
  }

  // Compute the global error
  double err=0,norm=0;
  for(long i=0;i<Q;i++){
    for(int dof=0;dof<data_dof;++dof){
      err=std::max(err,fabs((dof+1)*f_ref[i]-f_cubic[i+dof*Q]));
    }
    norm+=fabs(f_ref[i]);
  }

  if(0)
  if(procid==0){
    int i=1;
    std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
    std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
    std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
    std::cout<<"norm="<<norm<<std::endl;
    for(int dof=0;dof<data_dof;++dof){
      std::cout<<"f_true= "<<(dof+1)*f_ref[i]<<" f_cubic= "<<f_cubic[i+dof*Q]<<" err="<<fabs(f_cubic[i+dof*Q]-(dof+1)*f_ref[i])/f_ref[i]<<'\n';
    }
    std::cout<<"L_inf Error in interpolation of "<<Q<<" points is "<<err<<'\n';
    std::cout<<"Rel. L_inf Error in interpolation of "<<Q<<" points is "<<err/norm<<'\n';
  }
  double g_err,g_norm,g_time,g_interpolation_timer[4];
  MPI_Reduce(&err,&g_err,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&time,&g_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&interpolation_timer,&g_interpolation_timer,4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm,&g_norm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  PCOUT<<"\n Global Error:\n";
  PCOUT<<"L_inf Error in interpolation of "<<Q<<" points is "<<g_err<<'\n';
  PCOUT<<"Rel. L_inf Error in interpolation of "<<Q<<" points is "<<g_err/g_norm<<'\n';
  PCOUT<<"Time of "<<NREP<<" iterations ="<<g_time<<std::endl;
  PCOUT<<"Interp_timer[0]= "<<g_interpolation_timer[0]<<std::endl;
  PCOUT<<"Interp_timer[1]= "<<g_interpolation_timer[1]<<std::endl;
  PCOUT<<"Interp_timer[2]= "<<g_interpolation_timer[2]<<std::endl;
  PCOUT<<"Interp_timer[3]= "<<g_interpolation_timer[3]<<std::endl;


  accfft_free(data);
  accfft_free(data_hat);
  accfft_free(ghost_f);
  free(f);
  free(f_ref);
  free(f_cubic);
  free(query_points);

}

int main(int argc, char **argv)
{

  int NX,NY,NZ;
  MPI_Init (&argc, &argv);
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Parsing Inputs  */
  if(argc==1){
    NX=128;NY=128;NZ=128;
  }
  else{
    NX=atoi(argv[1]); NY=atoi(argv[2]); NZ=atoi(argv[3]);
  }
  int N[3]={NX,NY,NZ};

  //test1(N); // interpolate a scalar field
  //test2(N); // interpolated three velocity fields seperately
  //test3(N);   // interpolate three velocity fields at once
  test4(N);   // Same as test3 but with a planner

  MPI_Finalize();
  return 0;
} // end main


