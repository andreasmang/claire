#ifndef _READWRITEREG_CPP_
#define _READWRITEREG_CPP_



#include "ReadWriteReg.hpp"




namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadWriteReg"
ReadWriteReg::ReadWriteReg()
{
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadWriteReg"
ReadWriteReg::ReadWriteReg(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~ReadWriteReg"
ReadWriteReg::~ReadWriteReg()
{
    this->ClearMemory();
}




/********************************************************************
 * @brief init class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode ReadWriteReg::Initialize()
{

    this->m_Opt = NULL;
    this->m_Data = NULL;
    this->m_nx[0] = -1;
    this->m_nx[1] = -1;
    this->m_nx[2] = -1;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clear class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode ReadWriteReg::ClearMemory()
{
    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Read"
PetscErrorCode ReadWriteReg::Read(Vec* x, std::string filename)
{
    PetscErrorCode ierr;
    std::string file, msg;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // has to be allocated elsewhere (TODO: needs to be fixed)
//    if (*x != NULL){ ierr=VecDestroy(x); CHKERRQ(ierr); }

    // get file name without path
    ierr=GetFileName(file,filename); CHKERRQ(ierr);

    // display what we are doing
    if (this->m_Opt->GetVerbosity() > 1){
        msg = "reading " + file;
        ierr=DbgMsg(msg); CHKERRQ(ierr);
    }

    if (filename.find(".nii") != std::string::npos){
        ierr=this->ReadNII(x,filename); CHKERRQ(ierr);
    }
    else if (filename.find(".nii.gz") != std::string::npos){
        ierr=this->ReadNII(x,filename); CHKERRQ(ierr);
    }
    else if (filename.find(".hdr") != std::string::npos){
        ierr=this->ReadNII(x,filename); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("could not read: data type not supported"); CHKERRQ(ierr); }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Read"
PetscErrorCode ReadWriteReg::Read(VecField* v,
                                  std::string fnx1,
                                  std::string fnx2,
                                  std::string fnx3)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->Read(&v->m_X1,fnx1); CHKERRQ(ierr);
    ierr=this->Read(&v->m_X1,fnx2); CHKERRQ(ierr);
    ierr=this->Read(&v->m_X1,fnx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief write time series data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteTimeSeries"
PetscErrorCode ReadWriteReg::WriteTimeSeries(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    filename = this->m_Opt->GetXFolder() + filename;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief read temporal data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadTimeSeries"
PetscErrorCode ReadWriteReg::ReadTimeSeries(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief read data from file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadBlock"
PetscErrorCode ReadWriteReg::ReadBlock(Vec x, int isize[3], std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteBlock"
PetscErrorCode ReadWriteReg::WriteBlock(Vec x, int isize[3], std::string filename)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    filename = this->m_Opt->GetXFolder() + filename;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Write"
PetscErrorCode ReadWriteReg::Write(Vec x, std::string filename)
{
    PetscErrorCode ierr;
    std::string file,msg;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    // get file name without path
    ierr=GetFileName(file,filename); CHKERRQ(ierr);

    // display what we are doing
    if (this->m_Opt->GetVerbosity() > 2){
        msg = "writing " + file;
        ierr=DbgMsg(msg); CHKERRQ(ierr);
    }

    filename = this->m_Opt->GetXFolder() + filename;
    if (filename.find(".nii") != std::string::npos){
        ierr=this->WriteNII(x,filename); CHKERRQ(ierr);
    }
    else if (filename.find(".nii.gz") != std::string::npos){
        ierr=this->WriteNII(x,filename); CHKERRQ(ierr);
    }
    else if (filename.find(".hdr") != std::string::npos){
        ierr=this->WriteNII(x,filename); CHKERRQ(ierr);
    }
    else{ ierr=ThrowError("could not write: data type not supported"); CHKERRQ(ierr); }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief write data to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Write"
PetscErrorCode ReadWriteReg::Write(VecField* v,
                                   std::string fnx1,
                                   std::string fnx2,
                                   std::string fnx3)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(v != NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->Write(v->m_X1,fnx1); CHKERRQ(ierr);
    ierr=this->Write(v->m_X2,fnx2); CHKERRQ(ierr);
    ierr=this->Write(v->m_X3,fnx3); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief get component type of NII images
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "GetComponentTypeNII"
PetscErrorCode ReadWriteReg::GetComponentTypeNII(nifti_image* niiimage)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}



/********************************************************************
 * @brief read nifty image
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadNII"
PetscErrorCode ReadWriteReg::ReadNII(Vec* x, std::string filename)
{
    PetscErrorCode ierr;
    std::string msg,file;
    std::stringstream ss;
    int nprocs,rank,rval;
    IntType nx[3],isize[3],istart[3],*isize_c=NULL,*istart_c=NULL,ng,ngx,nl;
    ScalarType *p_x=NULL;
    ScalarType *values=NULL;
    nifti_image *niiimage=NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    // get file name without path
    ierr=GetFileName(file,filename); CHKERRQ(ierr);

    // check if file exists
    msg = "file " + file + "does not exist";
    ierr=Assert(FileExists(filename),msg); CHKERRQ(ierr);

    // read header file
    niiimage = nifti_image_read(filename.c_str(),false);
    msg="could not read nifti image " + file;
    ierr=Assert(niiimage != NULL,msg); CHKERRQ(ierr);

    // get number of grid points
    nx[0] = static_cast<IntType>(niiimage->nx);
    nx[1] = static_cast<IntType>(niiimage->ny);
    nx[2] = static_cast<IntType>(niiimage->nz);

    // if we read images, we want to make sure that they have the same size
    if (   (this->m_nx[0] == -1)
        && (this->m_nx[1] == -1)
        && (this->m_nx[2] == -1) ){

        for (int i = 0; i < 3; ++i){ this->m_nx[i] = nx[i]; }

        if(this->m_Opt->GetVerbosity() > 2){
            ss << "grid size ("<<nx[0]<<","<<nx[1]<<","<<nx[2]<<")";
            ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
        }

    }
    else{

        msg="grid size of images varies: perform affine registration first";
        for (int i = 0; i < 3; ++i){
            ierr=Assert(this->m_nx[i] == nx[i],msg); CHKERRQ(ierr);
        }

    }

    // pass number of grid points to options
    for (int i = 0; i < 3; ++i){ this->m_Opt->SetNumGridPoints(i,nx[i]); }

    // do the setup before running the code (this essentially
    // concerns the memory distribution/the setup of accfft
    if ( !this->m_Opt->SetupDone() ){
        ierr=this->m_Opt->DoSetup(); CHKERRQ(ierr);
    }

    // get local size
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;
    for (int i = 0; i < 3; ++i){
        isize[i] = this->m_Opt->GetDomainPara().isize[i];
        istart[i] = this->m_Opt->GetDomainPara().istart[i];
    }

    //check global size
    ngx=1; for(int i = 0; i < 3; ++i){ ngx*=nx[i]; }
    ierr=Assert(ng==ngx,"global size mismatch"); CHKERRQ(ierr);

    // allocate vector
    if ( *x != NULL ){ ierr=VecDestroy(x); CHKERRQ(ierr); }
    ierr=VecCreate(PETSC_COMM_WORLD,x); CHKERRQ(ierr);
    ierr=VecSetSizes(*x,nl,ng); CHKERRQ(ierr);
    ierr=VecSetFromOptions(*x); CHKERRQ(ierr);


    // read data only on master rank
    if(rank==0){

        // allocate data buffer
        if (this->m_Data != NULL){
            delete this->m_Data;
            this->m_Data=NULL;
        }

        // allocate data buffer
        try{ this->m_Data = new ScalarType[ng]; }
            catch(std::bad_alloc&){
            ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        ierr=this->ReadNII(niiimage,filename); CHKERRQ(ierr);

        // get all the sizes to read and assign data correctly
        if (isize_c==NULL){
            try{ isize_c = new IntType[3*nprocs]; }
            catch(std::bad_alloc&){
                ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        if (istart_c==NULL){
            try{ istart_c = new IntType[3*nprocs]; }
            catch(std::bad_alloc&){
                ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

    }

    int *numsend=NULL,*numoffset=NULL;
    try{ numsend = new int[3*nprocs]; }
    catch(std::bad_alloc&){
         ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try{ numoffset = new int[3*nprocs]; }
    catch(std::bad_alloc&){
         ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // gather isize and istart on master rank
    rval=MPI_Gather(isize,3,MPIU_INT,isize_c,3,MPIU_INT,0,PETSC_COMM_WORLD);
    ierr=MPIERRQ(rval); CHKERRQ(ierr);
    rval=MPI_Gather(istart,3,MPIU_INT,istart_c,3,MPIU_INT,0,PETSC_COMM_WORLD);
    ierr=MPIERRQ(rval); CHKERRQ(ierr);


    // compute offset and number of entries to send
    if (rank == 0){

        IntType offset = 0;
        for (int j = 0; j < nprocs; ++j){
            IntType nsend = 1;
            for (int i = 0; i < 3; ++i){
               nsend *= isize_c[j*3+i];
            }
            numsend[j] = static_cast<int>(nsend);
            numoffset[j] = offset;
            offset += nsend;
        }

        // allocate data buffer
        try{ values = new ScalarType[ng]; }
            catch(std::bad_alloc&){
            ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        IntType i1,i2,i3,j1,j2,j3,j,k=0;
            for(int p = 0; p < nprocs; ++p){

                for (i1 = 0; i1 < isize_c[3*p+0]; ++i1){ // x1
                    for (i2 = 0; i2 < isize_c[3*p+1]; ++i2){ // x2
                        for (i3 = 0; i3 < isize_c[3*p+2]; ++i3){ // x3

                            j1 = i1 + istart_c[3*p+0];
                            j2 = i2 + istart_c[3*p+1];
                            j3 = i3 + istart_c[3*p+2];

                            j = GetLinearIndex(j1,j2,j3,nx);
                            //p_niibuffer[j] = static_cast<T>(p_xg[k++]);
                            values[k++] = this->m_Data[j];

                        } // for i1
                    } // for i2
                } // for i3

            } // for all procs

    }

    int nrecv = static_cast<int>(nl);

//    Vec y; ScalarType*p_y=NULL;

//    ierr=VecCreate(PETSC_COMM_WORLD,&y); CHKERRQ(ierr);
//    ierr=VecSetSizes(y,nl,ng); CHKERRQ(ierr);
//    ierr=VecSetFromOptions(y); CHKERRQ(ierr);

    ierr=VecGetArray(*x,&p_x); CHKERRQ(ierr);

//    ierr=VecGetArray(y,&p_y); CHKERRQ(ierr);
    // TODO: switch between float and double
//    rval=MPI_Bcast(&this->m_Data[0], ng, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
//    ierr=Assert(rval==MPI_SUCCESS,"mpi gather returned an error"); CHKERRQ(ierr);

//    rval=MPI_Scatterv(this->m_Data,numsend,numoffset,
//                        MPI_DOUBLE,p_x,nrecv,MPI_DOUBLE,0,PETSC_COMM_WORLD);
    rval=MPI_Scatterv(values,numsend,numoffset,
                        MPI_DOUBLE,p_x,nrecv,MPI_DOUBLE,0,PETSC_COMM_WORLD);
    //rval=MPI_Scatterv(values,numsend,numoffset,
    //                    MPI_DOUBLE,p_y,nrecv,MPI_DOUBLE,0,PETSC_COMM_WORLD);
    //ierr=MPIERRQ(rval); CHKERRQ(ierr);
/*
    IntType sst=0;
    for (IntType i1 = 0; i1 < isize[0]; ++i1){ // x1
        for (IntType i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < isize[2]; ++i3){ // x3

                IntType l = GetLinearIndex(i1,i2,i3,isize);

                p_x[l] = p_y[sst++];

            }
        }
    }
    ierr=VecRestoreArray(y,&p_y); CHKERRQ(ierr);
    ierr=VecDestroy(&y); CHKERRQ(ierr);
*/
/*
    for (IntType i1 = 0; i1 < isize[0]; ++i1){ // x1
        for (IntType i2 = 0; i2 < isize[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < isize[2]; ++i3){ // x3

                IntType j1 = i1 + istart[0];
                IntType j2 = i2 + istart[1];
                IntType j3 = i3 + istart[2];

                IntType j = GetLinearIndex(j1,j2,j3,nx);
                IntType k = GetLinearIndex(i1,i2,i3,isize);

                p_x[k] = static_cast<ScalarType>(this->m_Data[j]);

            }
        }
    }
*/
    ierr=VecRestoreArray(*x,&p_x); CHKERRQ(ierr);

    // rescale image intensities to [0,1]
//    ierr=Rescale(*x,0.0,1.0); CHKERRQ(ierr);

    if (this->m_Data != NULL){
        delete this->m_Data;
        this->m_Data=NULL;
    }

    if (isize_c!=NULL){ delete [] isize_c; isize_c=NULL; }
    if (istart_c!=NULL){ delete [] istart_c; istart_c=NULL; }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}






/********************************************************************
 * @brief read nifty image with right component type
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadNII"
PetscErrorCode ReadWriteReg::ReadNII(nifti_image* niiimage,std::string filename)
{
    PetscErrorCode ierr;
    std::string msg;
    int rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    msg = "should only be called on master/root rank";
    ierr=Assert(rank==0,msg); CHKERRQ(ierr);

    switch (niiimage->datatype){
        case NIFTI_TYPE_UINT8:
        {
            this->m_ComponentType=UCHAR;
            ierr=this->ReadNII<unsigned char>(niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT8:
        {
            this->m_ComponentType=CHAR;
            ierr=this->ReadNII<char>(niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT16:
        {
            this->m_ComponentType=USHORT;
            ierr=this->ReadNII<unsigned short>(niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT16:
        {
            this->m_ComponentType=SHORT;
            ierr=this->ReadNII<short>(niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_UINT32:
        {
            this->m_ComponentType=UINT;
            ierr=this->ReadNII<unsigned int>(niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_INT32:
        {
            this->m_ComponentType=INT;
            ierr=this->ReadNII<int>(niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT32:
        {
            this->m_ComponentType=FLOAT;
            ierr=this->ReadNII<float>(niiimage,filename); CHKERRQ(ierr);
            break;
        }
        case NIFTI_TYPE_FLOAT64:
        {
            this->m_ComponentType=DOUBLE;
            ierr=this->ReadNII<double>(niiimage,filename); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("image data not supported"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}



/********************************************************************
 * @brief get component type of NII images
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ReadNII"
template <typename T> PetscErrorCode ReadWriteReg::ReadNII(nifti_image* niiimage,std::string filename)
{
    PetscErrorCode ierr;
    T *data=NULL;
    std::string msg;
    IntType ng;
    int rank;
    std::stringstream ss;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(this->m_Data != NULL,"null pointer"); CHKERRQ(ierr);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    msg = "should only be called on master/root rank";
    ierr=Assert(rank==0,msg); CHKERRQ(ierr);

    if (nifti_image_load(niiimage) == -1){
        msg="could not read image " + filename;
        ierr=ThrowError(msg); CHKERRQ(ierr);
    }

    // assign data
    data = static_cast<T*>(niiimage->data);
    ierr=Assert(data!=NULL,"null pointer"); CHKERRQ(ierr);

    // get global number of points
    ng = static_cast<IntType>(this->m_Opt->GetDomainPara().nglobal);

    // copy buffer
    for (IntType i = 0; i < ng; ++i){
        this->m_Data[i] = static_cast<ScalarType>(data[i]);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}



/********************************************************************
 * @brief write buffer to nii files
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
PetscErrorCode ReadWriteReg::WriteNII(Vec x,std::string filename)
{
    PetscErrorCode ierr;
    int rank;
    nifti_image* image=NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief write buffer to nii files
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
PetscErrorCode ReadWriteReg::WriteNII(nifti_image** niiimage,Vec x,std::string filename)
{
    PetscErrorCode ierr;
    std::string msg;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // if nifty image is null pointer default to double
    if ((*niiimage) == NULL){
        ierr=this->WriteNII<ScalarType>(niiimage,x,filename); CHKERRQ(ierr);
        this->m_Opt->Exit(__FUNCT__);
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

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}



/********************************************************************
 * @brief write buffer to nii files
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "WriteNII"
template <typename T> PetscErrorCode ReadWriteReg::WriteNII(nifti_image** niiimage,Vec x,std::string filename)
{
    PetscErrorCode ierr;
    IntType nglobal;
    T* p_niibuffer;
    ScalarType *p_xl=NULL,*p_xg=NULL;
    int nprocs,rank,rval;
    IntType *istart_c,*isize_c,istart[3],isize[3],nx[3];
    IntType i1,i2,i3,j1,j2,j3,j;
    Vec xg=NULL,xl=NULL;
    VecScatter vecscat;
    std::string msg;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // get number of ranks
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    // we have to copy the buffer to prevent that the scaling
    // is applied to the input data
    ierr=VecDuplicate(x,&xl); CHKERRQ(ierr);
    ierr=VecCopy(x,xl); CHKERRQ(ierr);

    // TODO
//    if(this->m_RescaleImage){
//        ierr=Rescale(xl,this->m_MinIValue,this->m_MaxIValue); CHKERRQ(ierr);
//    }

    // get local size
    nglobal = static_cast<IntType>(this->m_Opt->GetDomainPara().nglobal);

    // allocate the index buffers on master rank
    if (rank == 0){

        // we need to allocate the image if it's a zero pointer; this
        // will also create a standard header file; not tested (might need
        // to parse the dimensions of the data)
        if( (*niiimage) == NULL){
            if (this->m_Opt->GetVerbosity() >= 4){
                msg="allocating buffer for nifti image";
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

            try{ istart_c = new IntType[3*nprocs]; }
            catch(std::bad_alloc&){
                ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
            }

            try{ isize_c = new IntType[3*nprocs]; }
            catch(std::bad_alloc&){
                ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
            }

        }

        isize[0] = this->m_Opt->GetDomainPara().isize[0];
        isize[1] = this->m_Opt->GetDomainPara().isize[1];
        isize[2] = this->m_Opt->GetDomainPara().isize[2];

        istart[0] = this->m_Opt->GetDomainPara().istart[0];
        istart[1] = this->m_Opt->GetDomainPara().istart[1];
        istart[2] = this->m_Opt->GetDomainPara().istart[2];

        // gather the indices
        rval=MPI_Gather(istart,3,MPIU_INT,istart_c,3,MPIU_INT,0,PETSC_COMM_WORLD);
        ierr=MPIERRQ(rval); CHKERRQ(ierr);

        rval=MPI_Gather(isize,3,MPIU_INT,isize_c,3,MPIU_INT,0,PETSC_COMM_WORLD);
        ierr=MPIERRQ(rval); CHKERRQ(ierr);

        // create scatter object
        if(xg == NULL){
            ierr=VecScatterCreateToZero(xl,&vecscat,&xg); CHKERRQ(ierr);
        }

        // gather the data
        ierr=VecScatterBegin(vecscat,xl,xg,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr=VecScatterEnd(vecscat,xl,xg,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

        nx[0] = this->m_Opt->GetNumGridPoints(0);
        nx[1] = this->m_Opt->GetNumGridPoints(1);
        nx[2] = this->m_Opt->GetNumGridPoints(2);

        // if we are on master rank
        if (rank == 0){

            p_niibuffer = static_cast<T*>( (*niiimage)->data);
            ierr=VecGetArray(xg,&p_xg); CHKERRQ(ierr);

            int k = 0;
            for(int p = 0; p < nprocs; ++p){

                for (i1 = 0; i1 < isize_c[3*p+0]; ++i1){ // x1
                    for (i2 = 0; i2 < isize_c[3*p+1]; ++i2){ // x2
                        for (i3 = 0; i3 < isize_c[3*p+2]; ++i3){ // x3

                            j1 = i1 + istart_c[3*p+0];
                            j2 = i2 + istart_c[3*p+1];
                            j3 = i3 + istart_c[3*p+2];

                            j = GetLinearIndex(j1,j2,j3,nx);
                            p_niibuffer[j] = static_cast<T>(p_xg[k++]);
                            //p_niibuffer[j] = static_cast<T>(p_xg[j]);

                        } // for i1
                    } // for i2
                } // for i3

            } // for all procs

            ierr=VecRestoreArray(xg,&p_xg); CHKERRQ(ierr);

            // clear memory
            delete istart_c; istart_c=NULL;
            delete isize_c; isize_c=NULL;

        } // if on master

    }// else

    // clean that up
    ierr=VecDestroy(&xl); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief allocate buffer for nifty image
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "AllocateNII"
PetscErrorCode ReadWriteReg::AllocateNII(nifti_image** niiimage, Vec x)
{
    PetscErrorCode ierr;
    PetscInt n;
    int rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // get rank
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // init nifty image
    *niiimage = nifti_simple_init_nim();

    // dimensionalty of data: default is 5 (space, time, components)
    (*niiimage)->dim[0] = (*niiimage)->ndim = 5;

    ierr=VecGetLocalSize(x,&n); CHKERRQ(ierr);
    (*niiimage)->dim[1] = (*niiimage)->nx = this->m_Opt->GetNumGridPoints(0);
    (*niiimage)->dim[2] = (*niiimage)->ny = this->m_Opt->GetNumGridPoints(1);
    (*niiimage)->dim[3] = (*niiimage)->nz = this->m_Opt->GetNumGridPoints(2);

    (*niiimage)->pixdim[1] = static_cast<float>(this->m_Opt->GetDomainPara().hx[0]); // x direction
    (*niiimage)->pixdim[2] = static_cast<float>(this->m_Opt->GetDomainPara().hx[1]); // y direction
    (*niiimage)->pixdim[3] = static_cast<float>(this->m_Opt->GetDomainPara().hx[2]); // z direction

    // TODO: add temporal support
    if (n == this->m_Opt->GetDomainPara().nlocal){ // scalar field

        (*niiimage)->dim[4] = (*niiimage)->nt = 1;
        (*niiimage)->dim[5] = (*niiimage)->nu = 1;

        // temporal step size
        (*niiimage)->pixdim[4] = 1.0;

    }
    else if (n == 2*this->m_Opt->GetDomainPara().nlocal){ // 2D vector field

        (*niiimage)->dim[4] = (*niiimage)->nt = 1;
        (*niiimage)->dim[5] = (*niiimage)->nu = 2;

        // temporal step size
        (*niiimage)->pixdim[4] = 1.0;

        // step size (vector field)
        (*niiimage)->pixdim[5] = (*niiimage)->du = static_cast<float>(this->m_Opt->GetDomainPara().hx[0]);
        (*niiimage)->pixdim[6] = (*niiimage)->dv = static_cast<float>(this->m_Opt->GetDomainPara().hx[1]);

    }
    else if (n == 3*this->m_Opt->GetDomainPara().nlocal){ // 3D vector field

        (*niiimage)->dim[4] = (*niiimage)->nt = 1;
        (*niiimage)->dim[5] = (*niiimage)->nu = 3;

        // temporal step size
        (*niiimage)->pixdim[4] = 1.0;

        // step size (vector field)
        (*niiimage)->pixdim[5] = (*niiimage)->du = static_cast<float>(this->m_Opt->GetDomainPara().hx[0]);
        (*niiimage)->pixdim[6] = (*niiimage)->dv = static_cast<float>(this->m_Opt->GetDomainPara().hx[1]);
        (*niiimage)->pixdim[7] = (*niiimage)->dw = static_cast<float>(this->m_Opt->GetDomainPara().hx[2]);
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
        try { (*niiimage)->data = new ScalarType[(*niiimage)->nvox]; }
        catch (std::bad_alloc&){
            ierr=ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




} // end of name space



#endif //  _READWRITEREG_CPP_
