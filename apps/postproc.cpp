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
#include "VecField.hpp"
#include "ReadWriteReg.hpp"
#include "RegistrationInterface.hpp"



/********************************************************************
 * @brief post process image registration results
 *******************************************************************/
int main(int argc,char **argv)
{
    PetscErrorCode ierr;
    int procid,nprocs;
    std::string ifolder,filename;
    Vec mT=NULL,mR=NULL,vx1=NULL,vx2=NULL,vx3=NULL;
    reg::VecField *v=NULL;

    reg::RegOpt* regopt = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    reg::RegistrationInterface* registration = NULL;

    // initialize petsc (user is not allowed to set petsc options)
    ierr=PetscInitialize(0,(char***)NULL,(char*)NULL,(char*)NULL); CHKERRQ(ierr);

    PetscFunctionBegin;

    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD,&procid);

    // allocate class for controlling everything
    try{ regopt = new reg::RegOpt(argc,argv); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

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

    // clean up
    if (regopt != NULL){ delete regopt; regopt = NULL; }
    if (readwrite != NULL){ delete readwrite; readwrite = NULL; }
    if (registration != NULL){ delete registration; registration = NULL; }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); mT=NULL; }
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); mR=NULL; }
    if (v!=NULL){ delete v; v=NULL; }

    // clean up petsc
    ierr=PetscFinalize(); CHKERRQ(ierr);

    return 0;
}
