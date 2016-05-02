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
    this->m_NIIImage = NULL;

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
    else if (filename.find(".nii") != std::string::npos){
        ierr=this->ReadNII(&x,filename); CHKERRQ(ierr);
    }
    else if (filename.find(".nii.gz") != std::string::npos){
        ierr=this->ReadNII(&x,filename); CHKERRQ(ierr);
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
    else if (filename.find(".nii") != std::string::npos){
        ierr=this->WriteNII(x,filename); CHKERRQ(ierr);
    }
    else if (filename.find(".nii.gz") != std::string::npos){
        ierr=this->WriteNII(x,filename); CHKERRQ(ierr);
    }
    else if (filename.find(".hdr") != std::string::npos){
        ierr=this->WriteNII(x,filename); CHKERRQ(ierr);
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
 * Description: read netcdf time series
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


/********************************************************************
 * Name: GetComponentType
 * Description: get component type of NII images
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetComponentTypeNII"
PetscErrorCode DataReadWriteRegistration::GetComponentTypeNII(nifti_image* niiimage)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    switch (niiimage->datatype){
        case NIFTI_TYPE_UINT8:
        {
            this->m_ComponentType=UCHAR;
            break;
        }
        case NIFTI_TYPE_INT8:
        {
            this->m_ComponentType=CHAR;
            break;
        }
        case NIFTI_TYPE_UINT16:
        {
            this->m_ComponentType=USHORT;
            break;
        }
        case NIFTI_TYPE_INT16:
        {
            this->m_ComponentType=SHORT;
            break;
        }
        case NIFTI_TYPE_UINT32:
        {
            this->m_ComponentType=UINT;
            break;
        }
        case NIFTI_TYPE_INT32:
        {
            this->m_ComponentType=INT;
            break;
        }
        case NIFTI_TYPE_FLOAT32:
        {
            this->m_ComponentType=FLOAT;
            break;
        }
        case NIFTI_TYPE_FLOAT64:
        {
            this->m_ComponentType=DOUBLE;
            break;
        }
        default:
        {
            ierr=ThrowError("image data not supported"); CHKERRQ(ierr);
            break;
        }
    }

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: ReadNII
 * Description: read nifty image
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadNII"
PetscErrorCode DataReadWriteRegistration::ReadNII(Vec* x, std::string filename)
{
    PetscErrorCode ierr;
    std::string msg, file;

    PetscFunctionBegin;

    ierr=GetFileName(file,filename); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() >= 1){
        msg = "reading " + file;
        ierr=DbgMsg(msg); CHKERRQ(ierr);
    }

    // deallocate nii image
    if(this->m_NIIImage != NULL){
        nifti_image_free(this->m_NIIImage);
    }

    msg = "file " + file + "does not exist";
    ierr=Assert(FileExists(filename),msg); CHKERRQ(ierr);

    // read header file
    this->m_NIIImage = nifti_image_read(filename.c_str(),false);

    msg="could not read nifti image " + file;
    ierr=Assert(this->m_NIIImage != NULL,msg); CHKERRQ(ierr);

    ierr=this->ReadNII(x,this->m_NIIImage,filename); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}






/********************************************************************
 * Name: ReadNII
 * Description: read nifty image with right component type
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadNII"
PetscErrorCode DataReadWriteRegistration::ReadNII(Vec* x,nifti_image* niiimage,std::string filename)
{
    PetscErrorCode ierr;
    std::string msg;
    PetscFunctionBegin;

    switch (niiimage->datatype){
        case NIFTI_TYPE_UINT8:
        {
            this->m_ComponentType=UCHAR;
            ierr=this->ReadNII<unsigned char>(x,niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT8:
        {
            this->m_ComponentType=CHAR;
            ierr=this->ReadNII<char>(x,niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT16:
        {
            this->m_ComponentType=USHORT;
            ierr=this->ReadNII<unsigned short>(x,niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT16:
        {
            this->m_ComponentType=SHORT;
            ierr=this->ReadNII<short>(x,niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT32:
        {
            this->m_ComponentType=UINT;
            ierr=this->ReadNII<unsigned int>(x,niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT32:
        {
            this->m_ComponentType=INT;
            ierr=this->ReadNII<int>(x,niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT32:
        {
            this->m_ComponentType=FLOAT;
            ierr=this->ReadNII<float>(x,niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT64:
        {
            this->m_ComponentType=DOUBLE;
            ierr=this->ReadNII<double>(x,niiimage,filename); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("image data not supported"); CHKERRQ(ierr);
            break;
        }
    }

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: GetComponentType
 * Description: get component type of NII images
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadImageBuffer"
template <typename T>
PetscErrorCode DataReadWriteRegistration::ReadNII(Vec* x,nifti_image* niiimage,std::string filename)
{
    PetscErrorCode ierr;
    T* p_niibuffer;
    ScalarType *p_x;
    unsigned int nx[3],isize[3],istart[3];
    int nprocs, rank;
    std::string msg;
    std::stringstream ss;

    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    nx[0] = static_cast<unsigned int>(niiimage->nx);
    nx[1] = static_cast<unsigned int>(niiimage->ny);
    nx[2] = static_cast<unsigned int>(niiimage->nz);

    istart[0] = static_cast<unsigned int>(this->m_Opt->m_MiscOpt->istart[0]);
    istart[1] = static_cast<unsigned int>(this->m_Opt->m_MiscOpt->istart[1]);
    istart[2] = static_cast<unsigned int>(this->m_Opt->m_MiscOpt->istart[2]);

    isize[0] = static_cast<unsigned int>(this->m_Opt->m_MiscOpt->isize[0]);
    isize[1] = static_cast<unsigned int>(this->m_Opt->m_MiscOpt->isize[1]);
    isize[2] = static_cast<unsigned int>(this->m_Opt->m_MiscOpt->isize[2]);

    // TODO
    //if (!this->m_opt->memorydistset){
    //    ierr=SetupMemoryDistribution(this->m_opt); CHKERRQ(ierr);
    //}

    // TODO
//    if( *x != NULL ){ ierr=VecDestroy(x); CHKERRQ(ierr); }

//    ierr=VecCreate(PETSC_COMM_WORLD,x); CHKERRQ(ierr);
//    ierr=VecSetSizes(*x,this->m_Opt->GetNLocal(),this->m_Opt->GetNGlobal()); CHKERRQ(ierr);
//    ierr=VecSetFromOptions(*x); CHKERRQ(ierr);

    if (nifti_image_load(niiimage) == -1){
        msg="could not read image " + filename;
        ierr=ThrowError(msg); CHKERRQ(ierr);
    }

    p_niibuffer = static_cast<T*>(niiimage->data);

    if (p_niibuffer == NULL){
        msg="image buffer is null pointer";
        ierr=ThrowError(msg); CHKERRQ(ierr);
    }

    ierr=VecGetArray(*x,&p_x); CHKERRQ(ierr);

    for (unsigned int i1 = 0; i1 < isize[0]; ++i1){ // x1
        for (unsigned int i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (unsigned int i3 = 0; i3 < isize[2]; ++i3){ // x3

                unsigned int j1 = i1 + istart[0];
                unsigned int j2 = i2 + istart[1];
                unsigned int j3 = i3 + istart[2];

                IntType j = GetLinearIndex(j1,j2,j3,nx);
                IntType k = GetLinearIndex(i1,i2,i3,isize);

                p_x[k] = static_cast<ScalarType>(p_niibuffer[j]);

            }
        }
    }

    ierr=VecRestoreArray(*x,&p_x); CHKERRQ(ierr);

//    ierr=VecMin(*x,NULL,&this->m_MinIValue); CHKERRQ(ierr);
//    ierr=VecMax(*x,NULL,&this->m_MaxIValue); CHKERRQ(ierr);

    // rescale image intensities to [0,1]
    ierr=Rescale(*x,0.0,1.0); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}



/********************************************************************
 * Name: WriteNII
 * Description: write buffer to nii files
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
PetscErrorCode DataReadWriteRegistration::WriteNII(Vec x,std::string filename)
{
    PetscErrorCode ierr;
    int rank;
    nifti_image* image=NULL;
    PetscFunctionBegin;

    // TODO
    // we want to rescale the image to whatever
    // the input intensity range was
    //rescale = this->m_RescaleImage;
    //this->m_RescaleImage = true;

    // write x buffer to nifti image
    //ierr=this->WriteNII(&this->m_Image,x,filename); CHKERRQ(ierr);
    ierr=this->WriteNII(&image,x,filename); CHKERRQ(ierr);

    //this->m_RescaleImage = rescale;

    // at rank zero write out
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    if(rank == 0){ nifti_image_write(image); }

    PetscFunctionReturn(0);
}





/********************************************************************
 * Name: WriteNII
 * Description: write buffer to nii files
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
PetscErrorCode DataReadWriteRegistration::WriteNII(nifti_image** niiimage,Vec x,std::string filename)
{
    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    // if nifty image is null pointer default to double
    if ((*niiimage) == NULL){
        ierr=this->WriteNII<ScalarType>(niiimage,x,filename); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    switch ((*niiimage)->datatype){
        case NIFTI_TYPE_UINT8:
        {
            ierr=this->WriteNII<unsigned char>(niiimage,x,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT8:
        {
            ierr=this->WriteNII<char>(niiimage,x,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT16:
        {
            ierr=this->WriteNII<unsigned short>(niiimage,x,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT16:
        {
            ierr=this->WriteNII<short>(niiimage,x,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT32:
        {
            ierr=this->WriteNII<unsigned int>(niiimage,x,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT32:
        {
            ierr=this->WriteNII<int>(niiimage,x,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT32:
        {
            ierr=this->WriteNII<float>(niiimage,x,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT64:
        {
            ierr=this->WriteNII<double>(niiimage,x,filename); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("image data not supported"); CHKERRQ(ierr);
            break;
        }
    }

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: WriteNII
 * Description: write buffer to nii files
 * Author: Andreas Mang
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
template <typename T>
PetscErrorCode DataReadWriteRegistration::WriteNII(nifti_image** niiimage,Vec x,std::string filename)
{
    PetscErrorCode ierr;
    IntType nglobal;
    T* p_niibuffer;
    ScalarType *p_xl=NULL,*p_xg=NULL;
    int nprocs,rank,rval,*istart,*isize;
    unsigned int nx[3];
    Vec xg=NULL,xl=NULL;
    VecScatter vecscat;
    std::string msg;

    PetscFunctionBegin;

    // get number of ranks
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // we have to copy the buffer to prevent that the scaling
    // is applied to the input data
    ierr=VecDuplicate(x,&xl); CHKERRQ(ierr);
    ierr=VecCopy(x,xl); CHKERRQ(ierr);

    // TODO
//    if(this->m_RescaleImage){
//        ierr=Rescale(xl,this->m_MinIValue,this->m_MaxIValue); CHKERRQ(ierr);
//    }

    // get local size
    nglobal = this->m_Opt->GetNGlobal();

    // allocate the index buffers on master rank
    if (rank == 0){

        // we need to allocate the image if it's a zero pointer; this
        // will also create a standard header file; not tested (might need
        // to parse the dimensions of the data)
        if( (*niiimage) == NULL){
            if (this->m_Opt->GetVerbosity() >= 3){
                msg="initializing empty nifty image";
                ierr=DbgMsg(msg); CHKERRQ(ierr);
            }
            ierr=this->AllocateNII(niiimage,x); CHKERRQ(ierr);
        }

        msg="nifty image is null pointer";
        ierr=Assert((*niiimage)!=NULL,msg); CHKERRQ(ierr);

        // construct file name
        std::string file(filename);

        // get temp extension
        const char* exttemp = nifti_find_file_extension(file.c_str());
        if (exttemp == NULL){ exttemp=".nii"; }

        // set extension
        const std::string ext(exttemp);

        // get base file name
        char* bnametemp = nifti_makebasename(file.c_str());
        const std::string bname(bnametemp);
        free(bnametemp);

        // is file compressed
        const std::string::size_type sep = ext.rfind(".gz");
        const bool iscompressed = (sep == std::string::npos) ? false : true;

        if ( ext==".nii" || ext==".nii.gz" ){
            (*niiimage)->nifti_type = NIFTI_FTYPE_NIFTI1_1;
        }
        else if ( ext==".nia" ){
            (*niiimage)->nifti_type = NIFTI_FTYPE_ASCII;
        }
        else if ( ext==".hdr" || ext==".img" || ext==".hdr.gz" || ext==".img.gz" ){
            (*niiimage)->nifti_type = NIFTI_FTYPE_NIFTI1_2;
        }
        else{ ierr=ThrowError("file extension not supported"); CHKERRQ(ierr); }

        (*niiimage)->fname = nifti_makehdrname(bname.c_str(),
                                            (*niiimage)->nifti_type,
                                            false, iscompressed);

        (*niiimage)->iname = nifti_makeimgname(bname.c_str(),
                                            (*niiimage)->nifti_type,
                                            false, iscompressed);
    }



    // get array
    if( nprocs == 1 ){

        ierr=VecGetArray(xl,&p_xl); CHKERRQ(ierr);

        // cast pointer of nifti image data
        p_niibuffer = static_cast<T*>((*niiimage)->data);

#pragma omp parallel
{
#pragma omp for
        for (IntType i=0; i < nglobal; ++i){
            p_niibuffer[i] = static_cast<T>(p_xl[i]);
        }
} /// pragma omp parallel

        ierr=VecRestoreArray(xl,&p_xl); CHKERRQ(ierr);

    }
    else{

        // allocate the index buffers on master rank
        if (rank == 0){

            try{ istart = new int[3*nprocs]; }
            catch(std::bad_alloc&){
                ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
            }

            try{ isize = new int[3*nprocs]; }
            catch(std::bad_alloc&){
                ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
            }

        }

        // gather the indices
        rval=MPI_Gather(this->m_Opt->m_MiscOpt->istart,3,MPI_INT,istart,3,MPI_INT,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi gather returned an error"); CHKERRQ(ierr);

        rval=MPI_Gather(this->m_Opt->m_MiscOpt->isize,3,MPI_INT,isize,3,MPI_INT,0,PETSC_COMM_WORLD);
        ierr=Assert(rval==MPI_SUCCESS,"mpi gather returned an error"); CHKERRQ(ierr);

        // create scatter object
        if(xg == NULL){
            ierr=VecScatterCreateToZero(xl,&vecscat,&xg); CHKERRQ(ierr);
        }

        // gather the data
        ierr=VecScatterBegin(vecscat,xl,xg,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr=VecScatterEnd(vecscat,xl,xg,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

        nx[0] = this->m_Opt->m_MiscOpt->N[0];
        nx[1] = this->m_Opt->m_MiscOpt->N[1];
        nx[2] = this->m_Opt->m_MiscOpt->N[2];

        // if we are on master rank
        if (rank == 0){

            p_niibuffer = static_cast<T*>( (*niiimage)->data);
            ierr=VecGetArray(xg,&p_xg); CHKERRQ(ierr);

            unsigned int k = 0;
            for(unsigned int pid = 0; pid < nprocs; ++pid){

                for (unsigned int i1 = 0; i1 < isize[3*pid]; ++i1){ // x1
                    for (unsigned int i2 = 0; i2 < isize[3*pid+1]; ++i2){ // x2
                        for (unsigned int i3 = 0; i3 < isize[3*pid+2]; ++i3){ // x3

                            unsigned int j1 = i1 + istart[3*pid  ];
                            unsigned int j2 = i2 + istart[3*pid+1];
                            unsigned int j3 = i3 + istart[3*pid+2];

                            IntType j = GetLinearIndex(j1,j2,j3,nx);
                            p_niibuffer[j] = static_cast<T>(p_xg[k++]);

                        } // for i1
                    } // for i2
                } // for i3

            } // for all procs

          ierr=VecRestoreArray(xg,&p_xg); CHKERRQ(ierr);

            // clear memory
            delete istart; istart=NULL;
            delete isize; isize=NULL;

        } // if on master

    }// else

    // clean that up
    ierr=VecDestroy(&xl); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}



/********************************************************************
 * Name: AllocateNII
 * Description: allocate buffer for nifty image
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "AllocateNII"
PetscErrorCode DataReadWriteRegistration::AllocateNII(nifti_image** niiimage, Vec x)
{
    PetscErrorCode ierr;
    PetscInt n;
    int rank;
    ScalarType hx[3];

    PetscFunctionBegin;

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // init nifty image
    *niiimage = nifti_simple_init_nim();

    // dimensionalty of data: default is 5 (space, time, components)
    (*niiimage)->dim[0] = (*niiimage)->ndim = 5;

    ierr=VecGetLocalSize(x,&n); CHKERRQ(ierr);
    (*niiimage)->dim[1] = (*niiimage)->nx = this->m_Opt->m_MiscOpt->N[0];
    (*niiimage)->dim[2] = (*niiimage)->ny = this->m_Opt->m_MiscOpt->N[1];
    (*niiimage)->dim[3] = (*niiimage)->nz = this->m_Opt->m_MiscOpt->N[2];

    (*niiimage)->pixdim[1] = static_cast<float>(this->m_Opt->m_MiscOpt->h[0]); // x direction
    (*niiimage)->pixdim[2] = static_cast<float>(this->m_Opt->m_MiscOpt->h[1]); // y direction
    (*niiimage)->pixdim[3] = static_cast<float>(this->m_Opt->m_MiscOpt->h[2]); // z direction

    // TODO: add temporal support
    if (n == this->m_Opt->GetNLocal()){ // scalar field

        (*niiimage)->dim[4] = (*niiimage)->nt = 1;
        (*niiimage)->dim[5] = (*niiimage)->nu = 1;

        // temporal step size
        (*niiimage)->pixdim[4] = 1.0;

    }
    else if (n == 2*this->m_Opt->GetNLocal()){ // 2D vector field

        (*niiimage)->dim[4] = (*niiimage)->nt = 1;
        (*niiimage)->dim[5] = (*niiimage)->nu = 2;

        // temporal step size
        (*niiimage)->pixdim[4] = 1.0;

        // step size (vector field)
        (*niiimage)->pixdim[5] = (*niiimage)->du = static_cast<float>(this->m_Opt->m_MiscOpt->h[0]);
        (*niiimage)->pixdim[6] = (*niiimage)->dv = static_cast<float>(this->m_Opt->m_MiscOpt->h[1]);

    }
    else if (n == 3*this->m_Opt->GetNLocal()){ // 3D vector field

        (*niiimage)->dim[4] = (*niiimage)->nt = 1;
        (*niiimage)->dim[5] = (*niiimage)->nu = 3;

        // temporal step size
        (*niiimage)->pixdim[4] = 1.0;

        // step size (vector field)
        (*niiimage)->pixdim[5] = (*niiimage)->du = static_cast<float>(this->m_Opt->m_MiscOpt->h[0]);
        (*niiimage)->pixdim[6] = (*niiimage)->dv = static_cast<float>(this->m_Opt->m_MiscOpt->h[1]);
        (*niiimage)->pixdim[7] = (*niiimage)->dw = static_cast<float>(this->m_Opt->m_MiscOpt->h[2]);
    }

    // currently fixed to double precision: TODO flexible...
    (*niiimage)->datatype=NIFTI_TYPE_FLOAT64;
    (*niiimage)->nbyper=sizeof(ScalarType);

    (*niiimage)->nvox = 1;

    // compute number of voxels (space and time)
    for(int i = 1; i <= 4; ++i){
        (*niiimage)->nvox *= (*niiimage)->dim[i];
    }
    // add components
    (*niiimage)->nvox *= (*niiimage)->nu;


    // only allocate image buffer on root
    if (rank == 0){
        try { (*niiimage)->data = new ScalarType [(*niiimage)->nvox]; }
        catch (std::bad_alloc&){
            ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
        }

    }

    PetscFunctionReturn(0);
}



} // end of name space







#endif // _DATAREADWRITEREGISTRATION_CPP_
