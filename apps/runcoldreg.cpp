/**
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
#include "Optimizer.hpp"
#include "PreProcessingRegistration.hpp"
#include "SynProbRegistration.hpp"
#include "ReadWriteReg.hpp"
#include "LargeDeformationRegistration.hpp"
#include "OptimalControlRegistration.hpp"
#include "OptimalControlRegistrationIC.hpp"
#include "TaoInterfaceRegistration.hpp"



/********************************************************************
 * Name: main
 * Description: main function to run registration
 *******************************************************************/
int main(int argc,char **argv)
{
    PetscErrorCode ierr;
    int procid,nprocs;
    IntType nl,ng;
    Vec mT = NULL, mR=NULL;

    reg::RegOpt* regopt = NULL;
    reg::SynProbRegistration* synprob = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    reg::Optimizer* optimizer = NULL;
    reg::LargeDeformationRegistration* registration = NULL;

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

    // allocate class for io
    try{ optimizer = new reg::Optimizer(regopt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // allocate class for registration
    if (regopt->InCompressible()==false){

        try{ registration = new reg::OptimalControlRegistration(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

    }
    else{

        try{ registration = new reg::OptimalControlRegistrationIC(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

    }
    ierr=registration->SetIO(readwrite); CHKERRQ(ierr);

    // get sizes
    nl = regopt->GetNLocal();
    ng = regopt->GetNGlobal();
    if(regopt->ReadImagesFromFile()){

        ierr=VecCreate(PETSC_COMM_WORLD,&mR); CHKERRQ(ierr);
        ierr=VecSetSizes(mR,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mR); CHKERRQ(ierr);
        ierr=readwrite->Read(mR,regopt->GetReferenceFN()); CHKERRQ(ierr);
        ierr=reg::Assert(mR!=NULL, "input reference image is null pointer"); CHKERRQ(ierr);

        ierr=VecCreate(PETSC_COMM_WORLD,&mT); CHKERRQ(ierr);
        ierr=VecSetSizes(mT,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mT); CHKERRQ(ierr);
        ierr=readwrite->Read(mT,regopt->GetTemplateFN()); CHKERRQ(ierr);
        ierr=reg::Assert(mT!=NULL, "input template image is null pointer"); CHKERRQ(ierr);

        // pass to registration
        ierr=registration->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr=registration->SetTemplateImage(mT); CHKERRQ(ierr);

    }
    else{

        ierr=VecCreate(PETSC_COMM_WORLD,&mT); CHKERRQ(ierr);
        ierr=VecSetSizes(mT,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mT); CHKERRQ(ierr);
        ierr=VecSet(mT,0.0); CHKERRQ(ierr);

        // allocate class for registration
        try{ synprob = new reg::SynProbRegistration(regopt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        // set up a synthetic test problem/compute template image
        ierr=synprob->ComputeSmoothScalarField(mT,0); CHKERRQ(ierr);

        // advect template image to obtain reference image;
        // the variables for the images are set internally; so no need
        // to get the output here
        ierr=registration->SetupSyntheticProb(mT); CHKERRQ(ierr);
    }

    // reset all the clocks we have used so far
    ierr=regopt->ResetTimers(); CHKERRQ(ierr);
    ierr=regopt->ResetCounters(); CHKERRQ(ierr);

    // init solver
    ierr=optimizer->SetProblem(registration); CHKERRQ(ierr);

    // run the solver
    if ( regopt->DoParameterContinuation() ){
        ierr=optimizer->RunBetaCont(); CHKERRQ(ierr);
    }
    else{ ierr=optimizer->Run(); CHKERRQ(ierr); }

    // run the optimizer
    ierr=optimizer->Finalize(); CHKERRQ(ierr);


    // clean up
    if (regopt != NULL){ delete regopt; regopt = NULL; }
    if (synprob != NULL){ delete synprob; synprob = NULL; }
    if (readwrite != NULL){ delete readwrite; readwrite = NULL; }
    if (optimizer != NULL){ delete optimizer; optimizer = NULL; }
    if (registration != NULL){ delete registration; registration = NULL; }

    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); mT=NULL; }
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); mR=NULL; }

    ierr=PetscFinalize(); CHKERRQ(ierr);

   return 0;
}




