/*
 *  Copyright (c) 2015-2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 *
 */




#include "RegToolsOpt.hpp"
#include "RegUtils.hpp"
#include "PreProcReg.hpp"
#include "VecField.hpp"
#include "ReadWriteReg.hpp"
#include "SynProbRegistration.hpp"
#include "RegistrationInterface.hpp"

PetscErrorCode RunPostProcessing(reg::RegToolsOpt*);
PetscErrorCode ResampleVecField(reg::RegToolsOpt*);
PetscErrorCode ResampleScaField(reg::RegToolsOpt*);


/********************************************************************
 * @brief main function for registration tools
 *******************************************************************/
int main(int argc,char **argv)
{
    PetscErrorCode ierr;

    reg::RegToolsOpt* regopt = NULL;

    // initialize petsc (user is not allowed to set petsc options)
    ierr=PetscInitialize(0,(char***)NULL,(char*)NULL,(char*)NULL); CHKERRQ(ierr);

    PetscFunctionBegin;

    // allocate class for controlling everything
    try{ regopt = new reg::RegToolsOpt(argc,argv); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    if ( regopt->GetPostProcPara().enabled ){
        ierr=RunPostProcessing(regopt); CHKERRQ(ierr);
    }
    else if (regopt->GetRegFlags().storedefgrad ){
    }
    else if( regopt->GetResamplingPara().enabled ){

        if ( regopt->GetFlags().readvecfield ){
            ierr=ResampleVecField(regopt); CHKERRQ(ierr);
        }

        if ( regopt->GetFlags().readscafield ){
            ierr=ResampleScaField(regopt); CHKERRQ(ierr);
        }

    }

    // clean up
    if (regopt != NULL){ delete regopt; regopt = NULL; }

    // clean up petsc
    ierr=PetscFinalize(); CHKERRQ(ierr);

    return 0;
}




/********************************************************************
 * @brief post process image registration results
 *******************************************************************/
PetscErrorCode RunPostProcessing(reg::RegToolsOpt* regopt)
{
    PetscErrorCode ierr=0;
    std::string ifolder,xfolder,filename;
    Vec mT=NULL,mR=NULL,vx1=NULL,vx2=NULL,vx3=NULL;
    reg::VecField *v=NULL;

    reg::ReadWriteReg* readwrite = NULL;
    reg::RegistrationInterface* registration = NULL;

    PetscFunctionBegin;

    // allocate class for io
    try{ readwrite = new reg::ReadWriteReg(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ifolder=regopt->GetIFolder();
    ierr=reg::Assert(ifolder.empty()!=true,"input folder needs to be provided"); CHKERRQ(ierr);

    // read template image
    filename = ifolder + "template-image.nii.gz";
    ierr=readwrite->Read(&mT,filename); CHKERRQ(ierr);
    ierr=reg::Assert(mT!=NULL,"null pointer"); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ){ ierr=regopt->DoSetup(); CHKERRQ(ierr); }

    // read reference image
    filename = ifolder + "reference-image.nii.gz";
    ierr=readwrite->Read(&mR,filename); CHKERRQ(ierr);
    ierr=reg::Assert(mR!=NULL,"null pointer"); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ){ ierr=regopt->DoSetup(); CHKERRQ(ierr); }

    // allocate container for velocity field
    try{ v = new reg::VecField(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // read velocity components
    filename = ifolder + "velocity-field-x1.nii.gz";
    ierr=readwrite->Read(&vx1,filename); CHKERRQ(ierr);
    ierr=VecCopy(vx1,v->m_X1); CHKERRQ(ierr);

    filename = ifolder + "velocity-field-x2.nii.gz";
    ierr=readwrite->Read(&vx2,filename); CHKERRQ(ierr);
    ierr=VecCopy(vx2,v->m_X2); CHKERRQ(ierr);

    filename = ifolder + "velocity-field-x3.nii.gz";
    ierr=readwrite->Read(&vx3,filename); CHKERRQ(ierr);
    ierr=VecCopy(vx3,v->m_X3); CHKERRQ(ierr);

    // allocate class for registration interface
    try{ registration = new reg::RegistrationInterface(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // set all we need
    ierr=registration->SetReadWrite(readwrite); CHKERRQ(ierr);
    ierr=registration->SetReferenceImage(mR); CHKERRQ(ierr);
    ierr=registration->SetTemplateImage(mT); CHKERRQ(ierr);
    ierr=registration->SetInitialGuess(v); CHKERRQ(ierr);

    // run post processing
    ierr=registration->RunPostProcessing(); CHKERRQ(ierr);

    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); mT=NULL; }
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); mR=NULL; }

    if (vx1!=NULL){ ierr=VecDestroy(&vx1); CHKERRQ(ierr); vx1=NULL; }
    if (vx2!=NULL){ ierr=VecDestroy(&vx2); CHKERRQ(ierr); vx2=NULL; }
    if (vx3!=NULL){ ierr=VecDestroy(&vx3); CHKERRQ(ierr); vx3=NULL; }

    if (v!=NULL){ delete v; v=NULL; }

    if (readwrite!=NULL){ delete readwrite; readwrite=NULL; }
    if (registration!=NULL){ delete registration; registration=NULL; }

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief post process image registration results
 *******************************************************************/
PetscErrorCode ComputeDeformationGradient(reg::RegToolsOpt* regopt)
{
    PetscErrorCode ierr=0;
    std::string ifolder,xfolder,filename;
    Vec vx1=NULL,vx2=NULL,vx3=NULL;
    reg::VecField *v=NULL;

    reg::ReadWriteReg* readwrite = NULL;
    reg::RegistrationInterface* registration = NULL;

    PetscFunctionBegin;

    // allocate class for io
    try{ readwrite = new reg::ReadWriteReg(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ifolder=regopt->GetIFolder();
    ierr=reg::Assert(ifolder.empty()!=true,"input folder needs to be provided"); CHKERRQ(ierr);

    // read velocity components
    filename = ifolder + "velocity-field-x1.nii.gz";
    ierr=readwrite->Read(&vx1,filename); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ){ ierr=regopt->DoSetup(); CHKERRQ(ierr); }

    // allocate container for velocity field
    try{ v = new reg::VecField(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ierr=VecCopy(vx1,v->m_X1); CHKERRQ(ierr);
    if (vx1!=NULL){ ierr=VecDestroy(&vx1); CHKERRQ(ierr); vx1=NULL; }

    filename = ifolder + "velocity-field-x2.nii.gz";
    ierr=readwrite->Read(&vx2,filename); CHKERRQ(ierr);
    ierr=VecCopy(vx2,v->m_X2); CHKERRQ(ierr);
    if (vx2!=NULL){ ierr=VecDestroy(&vx2); CHKERRQ(ierr); vx2=NULL; }

    filename = ifolder + "velocity-field-x3.nii.gz";
    ierr=readwrite->Read(&vx3,filename); CHKERRQ(ierr);
    ierr=VecCopy(vx3,v->m_X3); CHKERRQ(ierr);
    if (vx3!=NULL){ ierr=VecDestroy(&vx3); CHKERRQ(ierr); vx3=NULL; }

    // allocate class for registration interface
    try{ registration = new reg::RegistrationInterface(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // set all we need
    ierr=registration->SetReadWrite(readwrite); CHKERRQ(ierr);
    ierr=registration->SetInitialGuess(v); CHKERRQ(ierr);

    // run post processing
    ierr=registration->ComputeDetDefGrad(); CHKERRQ(ierr);

    if (v!=NULL){ delete v; v=NULL; }
    if (readwrite!=NULL){ delete readwrite; readwrite=NULL; }
    if (registration!=NULL){ delete registration; registration=NULL; }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resample scalar field
 *******************************************************************/
PetscErrorCode ResampleScaField(reg::RegToolsOpt* regopt)
{
    PetscErrorCode ierr=0;
    std::string filename;
    IntType nl,ng,nx[3],nxl[2];
    ScalarType scale;
    Vec m=NULL,ml=NULL;
    reg::PreProcReg* preproc = NULL;
    reg::ReadWriteReg* readwrite = NULL;

    PetscFunctionBegin;

    // allocate class for io
    try{ readwrite = new reg::ReadWriteReg(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // read velocity components
    filename = regopt->GetScaFieldFN(0);
    ierr=readwrite->Read(&m,filename); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ){ ierr=regopt->DoSetup(); CHKERRQ(ierr); }

    // allocate container for velocity field
    try{ preproc = new reg::PreProcReg(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    preproc->ResetGridChangeOps(true);

    // compute grid size
    scale=regopt->GetResamplingPara().scale;
    for(int i=0; i < 3; ++i){
        nx[i]  = regopt->GetDomainPara().nx[i];
        nxl[i] = static_cast<IntType>(ceil(scale*regopt->GetDomainPara().nx[i]));
    }
    ierr=regopt->GetSizes(nxl,nl,ng); CHKERRQ(ierr);

    // allocate array
    ierr=reg::VecCreate(ml,nl,ng); CHKERRQ(ierr);

    // restrict of prolong the vector field
    if (scale > 1.0){ ierr=preproc->Prolong(&ml,m,nxl,nx); CHKERRQ(ierr); }
    else            { ierr=preproc->Restrict(&ml,m,nxl,nx); CHKERRQ(ierr); }

    // reset io
    if (readwrite!=NULL){ delete readwrite; readwrite = NULL; }
    for (int i=0; i<3; ++i){
        regopt->SetNumGridPoints(i,nxl[i]);
    }
    ierr=regopt->DoSetup(false); CHKERRQ(ierr);
    try{ readwrite = new reg::ReadWriteReg(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // read velocity components
    filename = regopt->GetScaFieldFN(1);
    ierr=readwrite->Write(ml,filename); CHKERRQ(ierr);

    if (m!=NULL){ ierr=VecDestroy(&m); CHKERRQ(ierr); }
    if (ml!=NULL){ ierr=VecDestroy(&ml); CHKERRQ(ierr); }

    if (preproc != NULL) { delete preproc; preproc=NULL; }
    if (readwrite != NULL) { delete readwrite; readwrite=NULL; }

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief resample vector field
 *******************************************************************/
PetscErrorCode ResampleVecField(reg::RegToolsOpt* regopt)
{
    PetscErrorCode ierr=0;
    std::string filename,fnx1,fnx2,fnx3;
    IntType nl,ng,nx[3],nxl[3];
    ScalarType scale;
    Vec vx1=NULL,vx2=NULL,vx3=NULL;
    reg::VecField *v=NULL,*vl=NULL;
    reg::PreProcReg* preproc = NULL;
    reg::ReadWriteReg* readwrite = NULL;

    PetscFunctionBegin;

    // allocate class for io
    try{ readwrite = new reg::ReadWriteReg(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // read velocity components
    filename = regopt->GetVecFieldFN(0,0);
    ierr=readwrite->Read(&vx1,filename); CHKERRQ(ierr);

    filename = regopt->GetVecFieldFN(1,0);
    ierr=readwrite->Read(&vx2,filename); CHKERRQ(ierr);

    filename = regopt->GetVecFieldFN(2,0);
    ierr=readwrite->Read(&vx3,filename); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ){ ierr=regopt->DoSetup(); CHKERRQ(ierr); }

    // allocate container for velocity field
    try{ v = new reg::VecField(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    ierr=VecCopy(vx1,v->m_X1); CHKERRQ(ierr);
    ierr=VecCopy(vx2,v->m_X2); CHKERRQ(ierr);
    ierr=VecCopy(vx3,v->m_X3); CHKERRQ(ierr);

    // allocate container for velocity field
    try{ preproc = new reg::PreProcReg(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    preproc->ResetGridChangeOps(true);

    // compute grid size
    scale=regopt->GetResamplingPara().scale;
    for(int i=0; i < 3; ++i){
        nx[i]  = regopt->GetDomainPara().nx[i];
        nxl[i] = scale*regopt->GetDomainPara().nx[i];
    }
    ierr=regopt->GetSizes(nxl,nl,ng); CHKERRQ(ierr);

    // allocate container for velocity field
    try{ vl = new reg::VecField(nl,ng); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // restrict of prolong the vector field
    if (scale > 1.0){ ierr=preproc->Prolong(vl,v,nxl,nx); CHKERRQ(ierr); }
    else            { ierr=preproc->Restrict(vl,v,nxl,nx); CHKERRQ(ierr); }

    // reset io
    if (readwrite!=NULL){ delete readwrite; readwrite = NULL; }
    for (int i=0; i<3; ++i){
        regopt->SetNumGridPoints(i,nxl[i]);
    }
    ierr=regopt->DoSetup(false); CHKERRQ(ierr);
    try{ readwrite = new reg::ReadWriteReg(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // get output file name (based on input file name)
    fnx1 = regopt->GetVecFieldFN(0,1);
    fnx2 = regopt->GetVecFieldFN(1,1);
    fnx3 = regopt->GetVecFieldFN(2,1);

    // write to file
    ierr=readwrite->Write(vl,fnx1,fnx2,fnx3); CHKERRQ(ierr);

    if (v != NULL) { delete v; v=NULL; }
    if (vl != NULL) { delete vl; vl=NULL; }

    if (vx1!=NULL){ ierr=VecDestroy(&vx1); CHKERRQ(ierr); vx1=NULL; }
    if (vx2!=NULL){ ierr=VecDestroy(&vx2); CHKERRQ(ierr); vx2=NULL; }
    if (vx3!=NULL){ ierr=VecDestroy(&vx3); CHKERRQ(ierr); vx3=NULL; }

    if (preproc != NULL) { delete preproc; preproc=NULL; }
    if (readwrite != NULL) { delete readwrite; readwrite=NULL; }

    PetscFunctionReturn(ierr);

}
