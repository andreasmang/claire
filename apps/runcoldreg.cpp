/*************************************************************************
 *  Copyright (c) 2016.
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
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#include <memory>
#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "ReadWriteReg.hpp"
#include "RegistrationInterface.hpp"



/********************************************************************
 * @brief main function to run registration
 *******************************************************************/
int main(int argc, char **argv) {
    PetscErrorCode ierr = 0;
    int procid, nprocs;
    Vec mT = NULL, mR = NULL, vxi = NULL;
    reg::VecField* v = NULL;
    reg::RegOpt* regopt = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    reg::RegistrationInterface* registration = NULL;

    // initialize petsc (user is not allowed to set petsc options)
    ierr = PetscInitialize(0, reinterpret_cast<char***>(NULL),
                              reinterpret_cast<char*>(NULL),
                              reinterpret_cast<char*>(NULL)); CHKERRQ(ierr);
    PetscFunctionBegin;

    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &procid);

    // allocate class for controlling everything
    try {regopt = new reg::RegOpt(argc, argv);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // allocate class for io
    try {registration = new reg::RegistrationInterface(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    if (regopt->GetReadWriteFlags().readfiles) {
        ierr = readwrite->Read(&mR, regopt->GetReadWriteFlags().mr); CHKERRQ(ierr);
        ierr = reg::Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = readwrite->Read(&mT, regopt->GetReadWriteFlags().mt); CHKERRQ(ierr);
        ierr = reg::Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);

        // pass to registration
        ierr = registration->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr = registration->SetTemplateImage(mT); CHKERRQ(ierr);
    } else {
        ierr = regopt->DoSetup(); CHKERRQ(ierr);
    }

    if (regopt->GetReadWriteFlags().readvelocity) {
        try {v = new reg::VecField(regopt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        ierr = readwrite->Read(&vxi, regopt->GetReadWriteFlags().vx1); CHKERRQ(ierr);
        ierr = reg::Assert(vxi != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = VecCopy(vxi, v->m_X1); CHKERRQ(ierr);
        if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

        ierr = readwrite->Read(&vxi, regopt->GetReadWriteFlags().vx2); CHKERRQ(ierr);
        ierr = reg::Assert(vxi != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = VecCopy(vxi, v->m_X2); CHKERRQ(ierr);
        if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

        ierr = readwrite->Read(&vxi, regopt->GetReadWriteFlags().vx3); CHKERRQ(ierr);
        ierr = reg::Assert(vxi != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = VecCopy(vxi, v->m_X3); CHKERRQ(ierr);
        if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

        ierr = registration->SetInitialGuess(v); CHKERRQ(ierr);
    }

    ierr = registration->SetReadWrite(readwrite); CHKERRQ(ierr);
    ierr = registration->Run(); CHKERRQ(ierr);

    // clean up
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (registration != NULL) {delete registration; registration = NULL;}
    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr); mT = NULL;}
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr); mR = NULL;}
    if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}
    if (regopt != NULL) {delete regopt; regopt = NULL;}
    if (v != NULL) {delete v; v = NULL;}

    ierr = reg::Finalize(); CHKERRQ(ierr);

    return 0;
}




