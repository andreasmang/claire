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




#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "PreProcReg.hpp"
#include "VecField.hpp"
#include "ReadWriteReg.hpp"
#include "SynProbRegistration.hpp"
#include "RegistrationInterface.hpp"




/********************************************************************
 * @brief post process image registration results
 *******************************************************************/
int main(int argc,char **argv)
{
    PetscErrorCode ierr;
    IntType nl,ng,nx[3],nxres[3],nxpro[3];
    int procid,nprocs;
    std::string ifolder,xfolder,filename;
    Vec mT=NULL,mR=NULL,m=NULL,mres=NULL,mpro=NULL,
        vx1=NULL,vx2=NULL,vx3=NULL;
    reg::VecField *v=NULL,*vres=NULL,*vpro=NULL;

    reg::RegOpt* regopt = NULL;
    reg::PreProcReg* preproc = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    reg::RegistrationInterface* registration = NULL;

    // initialize petsc (user is not allowed to set petsc options)
    ierr=PetscInitialize(0,(char***)NULL,(char*)NULL,(char*)NULL); CHKERRQ(ierr);

    PetscFunctionBegin;

    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD,&procid);

    // allocate class for controlling everything
    try{ regopt = new reg::RegOpt(argc,argv,1); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    if (regopt->GetRegFlags().runpostproc){

        ifolder=regopt->GetIFolder();
        ierr=reg::Assert(ifolder.empty()!=true,"input folder needs to be provided"); CHKERRQ(ierr);

        // allocate class for io
        try{ readwrite = new reg::ReadWriteReg(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        // read template image
        filename = ifolder + "template-image.nii.gz";
        ierr=readwrite->Read(&mT,filename); CHKERRQ(ierr);
        ierr=reg::Assert(mT!=NULL,"null pointer"); CHKERRQ(ierr);
        if ( !regopt->SetupDone() ){
            ierr=regopt->DoSetup(); CHKERRQ(ierr);
        }

        // read reference image
        filename = ifolder + "reference-image.nii.gz";
        ierr=readwrite->Read(&mR,filename); CHKERRQ(ierr);
        ierr=reg::Assert(mR!=NULL,"null pointer"); CHKERRQ(ierr);
        if ( !regopt->SetupDone() ){
            ierr=regopt->DoSetup(); CHKERRQ(ierr);
        }

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

    }
    else if(regopt->GetRegFlags().resampledata){
/*
        regopt->SetNumGridPoints(0,128);
        regopt->SetNumGridPoints(1,128);
        regopt->SetNumGridPoints(2,128);
*/
        if ( !regopt->SetupDone() ){
            ierr=regopt->DoSetup(); CHKERRQ(ierr);
        }

        // allocate container for velocity field
        try{ preproc = new reg::PreProcReg(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        preproc->ResetGridChangeOperators(true);

        reg::SynProbRegistration* synprob=NULL;
        // allocate class for io
        try{ synprob = new reg::SynProbRegistration(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        nl = regopt->GetDomainPara().nlocal;
        ng = regopt->GetDomainPara().nglobal;

        ierr=VecCreate(PETSC_COMM_WORLD,&m); CHKERRQ(ierr);
        ierr=VecSetSizes(m,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(m); CHKERRQ(ierr);

        //ierr=synprob->ComputeExpSin(m); CHKERRQ(ierr);
        ierr=synprob->ComputeSmoothScalarField(m,0); CHKERRQ(ierr);

        // allocate container for velocity field
        try{ v = new reg::VecField(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        ierr=VecCopy(m,v->m_X1); CHKERRQ(ierr);
        ierr=VecCopy(m,v->m_X2); CHKERRQ(ierr);
        ierr=VecCopy(m,v->m_X3); CHKERRQ(ierr);

        // if not yet allocted, do so
        // allocate class for io
        try{ readwrite = new reg::ReadWriteReg(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        // read reference image
        xfolder = "";
        filename = xfolder + "image.nii.gz";
        ierr=readwrite->Write(m,filename); CHKERRQ(ierr);

        delete readwrite; readwrite = NULL;

        nx[0] = regopt->GetDomainPara().nx[0];
        nx[1] = regopt->GetDomainPara().nx[1];
        nx[2] = regopt->GetDomainPara().nx[2];

        nxres[0] = nx[0]/2;
        nxres[1] = nx[1]/2;
        nxres[2] = nx[2]/2;

        ierr=regopt->GetSizes(nxres,nl,ng); CHKERRQ(ierr);

        ierr=VecCreate(PETSC_COMM_WORLD,&mres); CHKERRQ(ierr);
        ierr=VecSetSizes(mres,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mres); CHKERRQ(ierr);
        ierr=VecSet(mres,0.0); CHKERRQ(ierr);

        // allocate container for velocity field
        try{ vres = new reg::VecField(nl,ng); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=vres->SetValue(0.0); CHKERRQ(ierr);

        ierr=preproc->Restrict(&mres,m,nxres,nx); CHKERRQ(ierr);
        ierr=preproc->Restrict(vres,v,nxres,nx); CHKERRQ(ierr);

        nxpro[0] = 2*nx[0];
        nxpro[1] = 2*nx[1];
        nxpro[2] = 2*nx[2];

        ierr=regopt->GetSizes(nxpro,nl,ng); CHKERRQ(ierr);

        ierr=VecCreate(PETSC_COMM_WORLD,&mpro); CHKERRQ(ierr);
        ierr=VecSetSizes(mpro,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mpro); CHKERRQ(ierr);
        ierr=VecSet(mpro,0.0); CHKERRQ(ierr);

        // allocate container for velocity field
        try{ vpro = new reg::VecField(nl,ng); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=vpro->SetValue(0.0); CHKERRQ(ierr);

        ierr=preproc->Prolong(&mpro,m,nxpro,nx); CHKERRQ(ierr);
        ierr=preproc->Prolong(vpro,v,nxpro,nx); CHKERRQ(ierr);

        // initialize
        for (int i=0; i<3; ++i){
            regopt->SetNumGridPoints(i,nxres[i]);
        }
        ierr=regopt->DoSetup(false); CHKERRQ(ierr);

        // allocate class for io
        try{ readwrite = new reg::ReadWriteReg(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        // read reference image
        xfolder = "";
        filename = xfolder + "restricted-image.nii.gz";
        ierr=readwrite->Write(mres,filename); CHKERRQ(ierr);

        ierr=readwrite->Write(vres, xfolder+"restricted-vx1.nii.gz",
                                    xfolder+"restricted-vx2.nii.gz",
                                    xfolder+"restricted-vx3.nii.gz"); CHKERRQ(ierr);

        delete readwrite; readwrite = NULL;

        // initialize
        for (int i=0; i<3; ++i){
            regopt->SetNumGridPoints(i,nxpro[i]);
        }
        ierr=regopt->DoSetup(false); CHKERRQ(ierr);

        // allocate class for io
        try{ readwrite = new reg::ReadWriteReg(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        filename = xfolder + "prolonged-image.nii.gz";
        ierr=readwrite->Write(mpro,filename); CHKERRQ(ierr);

        ierr=readwrite->Write(vpro, xfolder+"prolonged-vx1.nii.gz",
                                    xfolder+"prolonged-vx2.nii.gz",
                                    xfolder+"prolonged-vx3.nii.gz"); CHKERRQ(ierr);

        delete readwrite; readwrite = NULL;

        delete synprob; synprob=NULL;

    }

    // clean up
    if (regopt != NULL){ delete regopt; regopt = NULL; }
    if (readwrite != NULL){ delete readwrite; readwrite = NULL; }
    if (preproc != NULL){ delete preproc; preproc = NULL; }
    if (registration != NULL){ delete registration; registration = NULL; }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); mT=NULL; }
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); mR=NULL; }
    if (m!=NULL){ ierr=VecDestroy(&m); CHKERRQ(ierr); m=NULL; }
    if (mpro!=NULL){ ierr=VecDestroy(&mpro); CHKERRQ(ierr); mpro=NULL; }
    if (mres!=NULL){ ierr=VecDestroy(&mres); CHKERRQ(ierr); mres=NULL; }
    if (v!=NULL){ delete v; v=NULL; }
    if (vres!=NULL){ delete vres; vres=NULL; }
    if (vpro!=NULL){ delete vpro; vpro=NULL; }

    // clean up petsc
    ierr=PetscFinalize(); CHKERRQ(ierr);

    return 0;
}

