/**
 *  Description: main test function for registration functionality
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
 *  along with PGLISTR. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "RegOpt.h"
#include "RegUtils.h"
#include "Optimizer.h"
#include "PreProcessingRegistration.h"
#include "SynProbRegistration.h"
#include "DataReadWriteRegistration.h"
#include "LargeDeformationRegistration.h"
#include "OptimalControlRegistration.h"
#include "OptimalControlRegistrationIC.h"
#include "TaoInterfaceRegistration.h"



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

    reg::RegOpt* opt = NULL;

    reg::LargeDeformationRegistration* registration = NULL;
    reg::DataReadWriteRegistration* io = NULL;
    reg::SynProbRegistration* synprob = NULL;
    reg::Optimizer* optimizer = NULL;

    // initialize petsc (user is not allowed to set petsc options)
    ierr=PetscInitialize(0,(char***)NULL,(char*)NULL,(char*)NULL); CHKERRQ(ierr);

    PetscFunctionBegin;

    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD,&procid);

    // allocate class for controlling everything
    try{ opt = new reg::RegOpt(argc,argv); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // allocate class for io
    try{ io = new reg::DataReadWriteRegistration(opt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // allocate class for io
    try{ optimizer = new reg::Optimizer(opt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // allocate class for registration
    if (opt->InCompressible()==false){

        try{ registration = new reg::OptimalControlRegistration(opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

    }
    else{

        try{ registration = new reg::OptimalControlRegistrationIC(opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

    }

    ierr=registration->SetIO(io); CHKERRQ(ierr);

    // get sizes
    nl = opt->GetNLocal();
    ng = opt->GetNGlobal();
    if(opt->ReadImagesFromFile()){

        ierr=VecCreate(PETSC_COMM_WORLD,&mR); CHKERRQ(ierr);
        ierr=VecSetSizes(mR,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mR); CHKERRQ(ierr);
        ierr=io->Read(mR,opt->GetReferenceFN()); CHKERRQ(ierr);
        ierr=reg::Assert(mR!=NULL, "input reference image is null pointer"); CHKERRQ(ierr);

        ierr=VecCreate(PETSC_COMM_WORLD,&mT); CHKERRQ(ierr);
        ierr=VecSetSizes(mT,nl,ng); CHKERRQ(ierr);
        ierr=VecSetFromOptions(mT); CHKERRQ(ierr);
        ierr=io->Read(mT,opt->GetTemplateFN()); CHKERRQ(ierr);
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
        try{ synprob = new reg::SynProbRegistration(opt); }
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
    ierr=opt->ResetTimers(); CHKERRQ(ierr);
    ierr=opt->ResetCounters(); CHKERRQ(ierr);

    // init solver
    ierr=optimizer->SetProblem(registration); CHKERRQ(ierr);

    // run the solver
    if (opt->DoParameterContinuation()){
        ierr=optimizer->RunBetaCont(); CHKERRQ(ierr);
    }
    else{ ierr=optimizer->Run(); CHKERRQ(ierr); }

    ierr=optimizer->Finalize(); CHKERRQ(ierr);

    ierr=opt->DisplayTimeToSolution(); CHKERRQ(ierr);

    // clean up
    if (io != NULL){ delete io; io = NULL; }
    if (opt != NULL){ delete opt; opt = NULL; }
    if (synprob != NULL){ delete synprob; synprob = NULL; }

    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); mT=NULL; }
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); mR=NULL; }
    if (optimizer != NULL){ delete optimizer; optimizer = NULL; }
    if (registration != NULL){ delete registration; registration = NULL; }

    ierr=PetscFinalize(); CHKERRQ(ierr);

   return 0;
}




