#ifndef _DATAREADWRITEREGISTRATION_CPP_
#define _DATAREADWRITEREGISTRATION_CPP_

#include "DataReadWriteRegistration.h"




namespace reg
{




/********************************************************************
 * Name: DataReadWriteRegistration
 * Description: default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DataReadWriteRegistration"
DataReadWriteRegistration::DataReadWriteRegistration()
{
    this->Initialize();
}




/********************************************************************
 * Name: DataReadWriteRegistration
 * Description: constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DataReadWriteRegistration"
DataReadWriteRegistration::DataReadWriteRegistration(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * Name: DataReadWriteRegistration
 * Description: default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~DataReadWriteRegistration"
DataReadWriteRegistration::~DataReadWriteRegistration()
{
    this->ClearMemory();
}




/********************************************************************
 * Name: Initialize
 * Description: init class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode DataReadWriteRegistration::Initialize()
{

    this->m_Opt = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ClearMemory
 * Description: clear class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode DataReadWriteRegistration::ClearMemory()
{
    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Read
 * Description: read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Read"
PetscErrorCode DataReadWriteRegistration::Read(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // has to be allocated elsewhere (TODO: needs to be fixed)
    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    if (filename.find(".nc") != std::string::npos){
        ierr=this->ReadNetCDF(x,filename); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("can not write data type to file"); CHKERRQ(ierr); }

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Read
 * Description: read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Read"
PetscErrorCode DataReadWriteRegistration::Read(VecField* v,
                                               std::string fnx1,
                                               std::string fnx2,
                                               std::string fnx3)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->Read(v->m_X1,fnx1); CHKERRQ(ierr);
    ierr=this->Read(v->m_X1,fnx2); CHKERRQ(ierr);
    ierr=this->Read(v->m_X1,fnx3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: WriteTimeSeries
 * Description: write time series data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteTimeSeries"
PetscErrorCode DataReadWriteRegistration::WriteTimeSeries(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    filename = this->m_Opt->GetXFolder() + filename;
    if (filename.find(".nc") != std::string::npos){
        ierr=this->WriteTimeSeriesNetCDF(x,filename); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("can not write data type to file"); CHKERRQ(ierr); }


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ReadTimeSeries
 * Description:read temporal data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadTimeSeries"
PetscErrorCode DataReadWriteRegistration::ReadTimeSeries(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    if (filename.find(".nc") != std::string::npos){
        ierr=this->ReadTimeSeriesNetCDF(x,filename); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("can not write data type to file"); CHKERRQ(ierr); }


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ReadBlock
 * Description: read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadBlock"
PetscErrorCode DataReadWriteRegistration::ReadBlock(Vec x, int isize[3], std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    if (filename.find(".nc") != std::string::npos){
        ierr=this->ReadBlockNetCDF(x,isize,filename); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("can not write data type to file"); CHKERRQ(ierr); }


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: WriteBlock
 * Description: write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteBlock"
PetscErrorCode DataReadWriteRegistration::WriteBlock(Vec x, int isize[3], std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    filename = this->m_Opt->GetXFolder() + filename;
    if (filename.find(".nc") != std::string::npos){
        ierr=this->WriteBlockNetCDF(x,isize,filename); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("can not write data type to file"); CHKERRQ(ierr); }


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Write
 * Description: write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Write"
PetscErrorCode DataReadWriteRegistration::Write(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    filename = this->m_Opt->GetXFolder() + filename;
    if (filename.find(".nc") != std::string::npos){
        ierr=this->WriteNetCDF(x,filename); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("can not write data type to file"); CHKERRQ(ierr); }


    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: Write
 * Description: write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Write"
PetscErrorCode DataReadWriteRegistration::Write(VecField* v,
                                                std::string fnx1,
                                                std::string fnx2,
                                                std::string fnx3)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->Write(v->m_X1,fnx1); CHKERRQ(ierr);
    ierr=this->Write(v->m_X2,fnx2); CHKERRQ(ierr);
    ierr=this->Write(v->m_X3,fnx3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ReadNetCDF
 * Description: read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadNetCDF"
PetscErrorCode DataReadWriteRegistration::ReadNetCDF(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    ScalarType *p_x=NULL;
    MPI_Offset isize[3],istart[3];
    MPI_Comm c;
    int nx[3];
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    isize[0] = this->m_Opt->m_MiscOpt->isize[0];
    isize[1] = this->m_Opt->m_MiscOpt->isize[1];
    isize[2] = this->m_Opt->m_MiscOpt->isize[2];

    istart[0] = this->m_Opt->m_MiscOpt->istart[0];
    istart[1] = this->m_Opt->m_MiscOpt->istart[1];
    istart[2] = this->m_Opt->m_MiscOpt->istart[2];

    c = this->m_Opt->m_MiscOpt->c_comm;

    nx[0] = this->m_Opt->m_MiscOpt->N[0];
    nx[1] = this->m_Opt->m_MiscOpt->N[1];
    nx[2] = this->m_Opt->m_MiscOpt->N[2];

    // call io of accfft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    read_pnetcdf(filename,istart,isize,c,nx,p_x);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ReadBlockNetCDF
 * Description: read data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadBlockNetCDF"
PetscErrorCode DataReadWriteRegistration::ReadBlockNetCDF(Vec x, int bsize[3], std::string filename)
{
    PetscErrorCode ierr;
    ScalarType *p_x=NULL;
    MPI_Offset isize[3],istart[3];
    MPI_Comm c;
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    isize[0] = bsize[0];
    isize[1] = bsize[1];
    isize[2] = bsize[2];

    istart[0] = 0;
    istart[1] = 0;
    istart[2] = 0;

    c = this->m_Opt->m_MiscOpt->c_comm;

    // call io of accfft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    read_pnetcdf(filename,istart,isize,c,bsize,p_x);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: WriteBlockNetCDF
 * Description: write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteBlockNetCDF"
PetscErrorCode DataReadWriteRegistration::WriteBlockNetCDF(Vec x, int bsize[3], std::string filename)
{
    PetscErrorCode ierr;
    ScalarType *p_x=NULL;
    MPI_Offset isize[3],istart[3];
    MPI_Comm c;
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    isize[0] = bsize[0];
    isize[1] = bsize[1];
    isize[2] = bsize[2];

    istart[0] = 0;
    istart[1] = 0;
    istart[2] = 0;

    c = this->m_Opt->m_MiscOpt->c_comm;

    // call io of accfft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    write_pnetcdf(filename,istart,isize,c,bsize,p_x);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: WriteNetCDF
 * Description: write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteNetCDF"
PetscErrorCode DataReadWriteRegistration::WriteNetCDF(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    ScalarType *p_x=NULL;
    MPI_Offset isize[3],istart[3];
    MPI_Comm c;
    int nx[3];
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    isize[0] = this->m_Opt->m_MiscOpt->isize[0];
    isize[1] = this->m_Opt->m_MiscOpt->isize[1];
    isize[2] = this->m_Opt->m_MiscOpt->isize[2];

    istart[0] = this->m_Opt->m_MiscOpt->istart[0];
    istart[1] = this->m_Opt->m_MiscOpt->istart[1];
    istart[2] = this->m_Opt->m_MiscOpt->istart[2];

    c = this->m_Opt->m_MiscOpt->c_comm;

    nx[0] = this->m_Opt->m_MiscOpt->N[0];
    nx[1] = this->m_Opt->m_MiscOpt->N[1];
    nx[2] = this->m_Opt->m_MiscOpt->N[2];

    // call io of accfft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    write_pnetcdf(filename,istart,isize,c,nx,p_x);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}





/********************************************************************
 * Name: WriteTimeSeriesNetCDF
 * Description: write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteTimeSeriesNetCDF"
PetscErrorCode DataReadWriteRegistration::WriteTimeSeriesNetCDF(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    ScalarType *p_x=NULL, *p_xj=NULL;
    Vec xj;
    std::string::size_type pos;
    std::ostringstream ss;;
    std::string fn;
    MPI_Offset is[3],in[3];
    MPI_Comm c;
    IntType nl,ng,nt;
    int nx[3];
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    for (unsigned int i=0; i < 3; ++i){
        nx[i] = this->m_Opt->m_MiscOpt->N[i];
        in[i] = this->m_Opt->m_MiscOpt->isize[i];
        is[i] = this->m_Opt->m_MiscOpt->istart[i];
    }

    c  = this->m_Opt->m_MiscOpt->c_comm;
    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ng = this->m_Opt->GetNGlobal();

    ierr=VecCreate(PETSC_COMM_WORLD,&xj); CHKERRQ(ierr);
    ierr=VecSetSizes(xj,nl,ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(xj); CHKERRQ(ierr);

    // call io of accfft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    ierr=VecGetArray(xj,&p_xj); CHKERRQ(ierr);
    pos = filename.find_last_of(".");
    fn = (std::string::npos == pos) ? filename : filename.substr(0, pos);

    for (unsigned int j = 0; j <= nt; ++j){

        try{ std::copy(p_x+j*nl,p_x+(j+1)*nl,p_xj); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }

        // construct file name
        ss << std::setw(3) << std::setfill('0') << j;
        filename = fn + "-j-" + ss.str() + ".nc"; ss.str("");
        // write to file
        write_pnetcdf(filename,is,in,c,nx,p_xj);
    }

    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);
    ierr=VecRestoreArray(xj,&p_xj); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * Name: ReadTimeSeriesNetCDF
 * Description: write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadTimeSeriesNetCDF"
PetscErrorCode DataReadWriteRegistration::ReadTimeSeriesNetCDF(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    ScalarType *p_x=NULL, *p_xj=NULL;
    Vec xj;
    std::string::size_type pos;
    std::ostringstream ss;;
    std::string fn;
    MPI_Offset is[3],in[3];
    MPI_Comm c;
    IntType nl,ng,nt;
    int nx[3];

    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    for (unsigned int i=0; i < 3; ++i){
        nx[i] = this->m_Opt->m_MiscOpt->N[i];
        in[i] = this->m_Opt->m_MiscOpt->isize[i];
        is[i] = this->m_Opt->m_MiscOpt->istart[i];
    }

    c  = this->m_Opt->m_MiscOpt->c_comm;
    nt = this->m_Opt->GetNumTimePoints();
    nl = this->m_Opt->GetNLocal();
    ng = this->m_Opt->GetNGlobal();

    ierr=VecCreate(PETSC_COMM_WORLD,&xj); CHKERRQ(ierr);
    ierr=VecSetSizes(xj,nl,ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(xj); CHKERRQ(ierr);

    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    ierr=VecGetArray(xj,&p_xj); CHKERRQ(ierr);
    pos = filename.find_last_of(".");
    fn = (std::string::npos == pos) ? filename : filename.substr(0, pos);

    for (unsigned int j = 0; j <= nt; ++j){

        // construct file name
        ss << std::setw(3) << std::setfill('0') << j;
        filename = fn + "-j-" + ss.str() + ".nc"; ss.str("");
        // write to file
        read_pnetcdf(filename,is,in,c,nx,p_xj);
        try{ std::copy(p_xj,p_xj+nl,p_x+j*nl); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }

    }

    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);
    ierr=VecRestoreArray(xj,&p_xj); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}





} // end of name space







#endif // _DATAREADWRITEREGISTRATION_CPP_
