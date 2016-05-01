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

//    msg = "file " + file + "does not exist";
//    ierr=Assert(FileExists(this->m_FileName),msg); CHKERRQ(ierr);

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



} // end of name space







#endif // _DATAREADWRITEREGISTRATION_CPP_
