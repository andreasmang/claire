
/**
 *  Description: This file is used to dump out input data A, into a
 *  NetCDF file format. The output can be visualized using paraview.
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
#include <mpi.h>
#include <utils.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>  // cout width
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <accfft.h>
#include <accfft_operators.h>
static bool isLittleEndian(){
  uint16_t number = 0x1;
  uint8_t *numPtr = (uint8_t*)&number;
  return (numPtr[0] == 1);
}

void DataOut(double * A, N_MISC N_Misc,const char * fname){

  MPI_Comm c_comm=N_Misc.c_comm;
  int nprocs, procid;
  MPI_Comm_rank(c_comm, &procid);
  MPI_Comm_size(c_comm, &nprocs);

  if(A==NULL){
    PCOUT<<"Error in DataOut ---> Input data is null"<<std::endl;
    return;
  }
  /* Write the output */
  int* istart=N_Misc.istart;
  int* isize=N_Misc.isize;

  std::string filename;
  MPI_Offset istart_mpi[3] = { istart[0], istart[1], istart[2] };
  MPI_Offset isize_mpi[3]  = { isize[0],  isize[1],  isize[2] };
  filename =fname;
  write_pnetcdf(filename,istart_mpi,isize_mpi,c_comm,N_Misc.N,A);
  return;

}
/*
void DataOut(double * A, const int * N,const char * fname){
  MPI_Comm comm=MPI_COMM_WORLD;
  int nproc, proc_id;
  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&proc_id);
  //const char * fname="sample";
  //Open files for writing.
  int binary=1;
    typedef float VTKtype;
  std::stringstream vtifname;
  vtifname<<"./results/"<<fname<<"_"<<std::setfill('0')<<std::setw(2)<<proc_id<<".vti";
  std::ofstream vtifile;
  vtifile.open(vtifname.str().c_str());
  if(vtifile.fail()){ std::cout<<"Failed to open file for writing output"<<std::endl;return;}

  int istart[3],isize[3],iend[3];
  p3dfft_get_dims(istart,iend,isize,1);
  long long int Nlocal=isize[0]*isize[1]*isize[2];

  vtifile.precision(3);
  vtifile<<"<?xml version=\"1.0\"?>\n";
  if(isLittleEndian()) vtifile<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  else                 vtifile<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"BigEndian\">\n";
  vtifile<<"<ImageData WholeExtent=\""<<istart[0]-1<<" "<<iend[0]<<" "<<istart[1]-1<<" "<<iend[1]<<" "<<istart[2]-1<<" "<<iend[2]<<"\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
  vtifile<<"<Piece Extent=\""<<istart[0]-1<<" "<<iend[0]<<" "<<istart[1]-1<<" "<<iend[1]<<" "<<istart[2]-1<<" "<<iend[2]<<"\" GhostLevel=\"0\" >\n";
  vtifile<<"<CellData Scalars=\"value\">\n";
  if(binary==0){
    vtifile<<"<DataArray Name=\"value\" type=\"Float"<<sizeof(VTKtype)*8<<"\" format=\"ascii\">\n";
    for(int z=0;z < isize[2];z++){
      for(int y=0;y < isize[1];y++)
      {
        for(int x=0;x < isize[0];x++)
        {
          vtifile<<" "<<A[x+isize[0]*y+isize[0]*isize[1]*z];
        }
        vtifile<<"\n";
      }
      vtifile<<"\n";
    }



    vtifile<<"\n</DataArray>\n";
    vtifile<<"</CellData>\n";
    vtifile<<"</Piece>\n";
    vtifile<<"</ImageData>\n";
    vtifile<<"</VTKFile>\n";
    vtifile.close();
  }
  else if(binary==1){
    uint32_t binarysize=Nlocal*sizeof(VTKtype);
    VTKtype *A_map=(VTKtype *) malloc(sizeof(VTKtype)*Nlocal);
#pragma omp parallel for
    for (int i=0; i<Nlocal; i++)
      A_map[i]=(VTKtype)A[i];

    vtifile<<"<DataArray Name=\"value\" type=\"Float"<<sizeof(VTKtype)*8<<"\" format=\"appended\" offset=\"0\">\n";
    vtifile<<"\n</DataArray>\n";
    vtifile<<"</CellData>\n";
    vtifile<<"</Piece>\n";
    vtifile<<"</ImageData>\n";
    vtifile<<"<AppendedData encoding=\"raw\">\n";
    vtifile<<"_";
    vtifile.close();
    vtifile.open(vtifname.str().c_str(), std::ios::out | std::ios::app | std::ios::binary | std::ios::ate);



    vtifile.write((char*)&binarysize,1*sizeof(uint32_t));
    vtifile.write((char*)A_map,Nlocal*sizeof(VTKtype));
    vtifile.close();

    vtifile.open(vtifname.str().c_str(), std::ios::out | std::ios::app);
    vtifile<<'\n';
    vtifile<<"</AppendedData>\n";
    vtifile<<"</VTKFile>\n";
    vtifile.close();
    free(A_map);
  }
  //Write parent file.
  // To do so first have to gather start and end of all procs
  int * ISTART, * IEND;
  if ( proc_id==0) {
    ISTART = (int *)malloc(nproc*3*sizeof(int));
    IEND = (int *)malloc(nproc*3*sizeof(int));
  }
  int error =MPI_Gather( istart, 3, MPI_INT, ISTART, 3, MPI_INT, 0, comm);
  if(error!=MPI_SUCCESS) { std::cout<<"error in MPI_Gather= "<<error<<std::endl; exit(-1);}
  error =MPI_Gather( iend, 3, MPI_INT, IEND, 3, MPI_INT, 0, comm);
  if(error!=MPI_SUCCESS) { std::cout<<"error in MPI_Gather= "<<error<<std::endl; exit(-1);}
  if(proc_id==0){
    std::stringstream pvtifname;
    pvtifname<<"./results/"<<fname<<".pvti";

    std::ofstream pvtifile;
    pvtifile.open(pvtifname.str().c_str());
    if(pvtifile.fail()) {std::cout<<"Failed to open file for writing output"<<std::endl;return;}

    pvtifile<<"<?xml version=\"1.0\"?>\n";
    pvtifile<<"<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvtifile<<"<PImageData WholeExtent=\""<<0<<" "<<N[0]<<" "<<0<<" "<<N[1]<<" "<<0<<" "<<N[2]<<"\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    pvtifile<<"<PCellData Scalars=\"value\">\n";
    pvtifile<<"<PDataArray Name=\"value\" type=\"Float"<<sizeof(VTKtype)*8<<"\"/>\n";
    pvtifile<<"</PCellData>\n";
    for (int i=0;i<nproc;i++){
      pvtifile<<"<Piece Extent=\""<<ISTART[3*i]-1<<" "<<IEND[3*i]<<" "<<ISTART[3*i+1]-1<<" "<<IEND[3*i+1]<<" "<<ISTART[3*i+2]-1<<" "<<IEND[3*i+2]<<"\" Source=\""<<fname<<"_"<<std::setfill('0')<<std::setw(2)<<i<<".vti\""<<"/>\n";
    }
    vtifname<<fname<<"_"<<std::setfill('0')<<std::setw(2)<<proc_id<<".vti";
    pvtifile<<"</PImageData>\n";
    pvtifile<<"</VTKFile>\n";


    pvtifile.close();
  }

}
*/
