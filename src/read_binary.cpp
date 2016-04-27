
/**
 *  Description: This code reads binary input data and writes into data ptr.
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
#include <accfft.h>
#include <utils.hpp>
#include <iostream>
#include <mpi.h>
#include <sstream>      // std::stringstream
#include <math.h>
#include <cmath>        // std::abs
#include <omp.h>
#include <cstring>
#include <petscksp.h>

int read_binary(double * data, N_MISC N,const char * fname,const char* prefix){

  MPI_Comm c_comm=N.c_comm;
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);


	double self_exec_time=-MPI_Wtime();
	const int Nx=N.N[0], Ny=N.N[1], Nz=N.N[2];


  int istart[3], isize[3], osize[3],ostart[3];
  int alloc_max=accfft_local_size_dft_r2c(N.N,isize,istart,osize,ostart,N.c_comm);

	const int N_local=isize[0]*isize[1]*isize[2];
	const int N_global=Nx*Ny*Nz;

	PetscErrorCode ierr;


  std::stringstream str;

  str<<prefix<<fname;

  double * buffer;
  long lSize;
  lSize=Nx*Ny*Nz;
  buffer = (double*) malloc (sizeof(double)*lSize);
  // First only the root process reads the data
  if(procid==0){
    FILE * pFile;
    size_t result;

    pFile = fopen ( str.str().c_str(), "rb" );
    if (pFile==NULL) {
      std::cout<<"File error"<<str.str().c_str()<<std::endl;
      exit (1);
    }

    // obtain file size:

    // allocate memory to contain the whole file:
    if (buffer == NULL) {
      std::cout<<"Memory error\n"<<std::endl;
      exit (2);
    }

    // copy the file into the buffer:
    result = fread (buffer,sizeof(double),lSize,pFile);
    if (result != lSize) {
      std::cout<<"Error reading file in binary read"<<std::endl;
      exit (3);
    }
    fclose (pFile);
  }

  // Now BCast the data to all processes
  MPI_Bcast((void*) buffer,N_global, MPI_DOUBLE,0,MPI_COMM_WORLD);
  self_exec_time +=  MPI_Wtime();  // Total time of the function
  double max_time=0;
  MPI_Reduce(&self_exec_time,&max_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);//get the maximum difference from all procs
  PCOUT<<"time= "<<max_time<<std::endl;



  // Now only keep parts that you need
#pragma omp parallel
  {
  long long int ptr=0;
#pragma omp for
  for(int z=0;z < isize[2];z++)
    for(int y=0;y < isize[1];y++)
      for(int x=0;x < isize[0];x++){
        ptr=x*isize[1]*isize[2]+y*isize[2]+z;
        data[ptr]=buffer[(x+istart[0])*Nz*Ny+Nz*(y+istart[1])+(z+istart[2])];

        //ptr=x + isize[0]*y + isize[0]*isize[1]*z;
        //data[ptr]=buffer[(x+istart[0]-1)+Nx*(y+istart[1]-1)+Nx*Ny*(z+istart[2]-1)];
      }
}
	Vec dummy;
  PetscReal norm;
  double * dummy_ptr;
  ierr = VecCreate(PETSC_COMM_WORLD,&dummy);CHKERRQ(ierr);
  ierr = VecSetSizes(dummy,N_local,N_global);CHKERRQ(ierr);
  ierr = VecSetFromOptions(dummy);CHKERRQ(ierr);

  ierr = VecGetArray(dummy, &dummy_ptr);CHKERRQ(ierr);
  std::memcpy(dummy_ptr, data, sizeof(double)*N_local);
  ierr = VecRestoreArray(dummy, &dummy_ptr);CHKERRQ(ierr);

  VecSum(dummy, &norm);
  PCOUT<<"norm PETSC= "<<norm/N_global<<std::endl;


  free (buffer);
  VecDestroy(&dummy);

  PCOUT<<"DONE with READING"<<std::endl;
  return 0;
}
