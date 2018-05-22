/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CLAIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#include <memory>
#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "ReadWriteReg.hpp"
#include "CLAIREInterface.hpp"




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
    reg::CLAIREInterface* registration = NULL;
    std::stringstream ss;

    // initialize petsc (user is not allowed to set petsc options)
    ierr = PetscInitialize(0, reinterpret_cast<char***>(NULL),
                              reinterpret_cast<char*>(NULL),
                              reinterpret_cast<char*>(NULL)); CHKERRQ(ierr);
    PetscFunctionBegin;

    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &procid);

    // allocate class for controlling everything
    try {regopt = new reg::RegOpt(argc, argv);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }

    if (regopt->m_Log.memoryusage) {
        ierr = PetscMemorySetGetMaximumUsage(); CHKERRQ(ierr);
    }

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }

    // allocate class for io
    try {registration = new reg::CLAIREInterface(regopt);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }

    // read reference and template image
    if (regopt->m_ReadWriteFlags.readfiles) {
        if (regopt->m_Verbosity > 1) {
            ierr = reg::DbgMsg("reading reference image"); CHKERRQ(ierr);
        }
        ierr = readwrite->ReadR(&mR, regopt->m_FileNames.mr); CHKERRQ(ierr);
        ierr = reg::Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);
        if (regopt->m_Verbosity > 1) {
            ierr = reg::DbgMsg("reading template image"); CHKERRQ(ierr);
        }
        ierr = readwrite->ReadT(&mT, regopt->m_FileNames.mt); CHKERRQ(ierr);
        ierr = reg::Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);

        // pass to registration
        ierr = registration->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr = registration->SetTemplateImage(mT); CHKERRQ(ierr);
    } else {
        ierr = regopt->DoSetup(); CHKERRQ(ierr);
    }

    // read velocity field
    if (regopt->m_ReadWriteFlags.readvelocity) {
        try {v = new reg::VecField(regopt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        if (regopt->m_Verbosity > 1) {
            ierr = reg::DbgMsg("reading velocity field"); CHKERRQ(ierr);
        }
        ierr = readwrite->Read(&vxi, regopt->m_FileNames.iv1); CHKERRQ(ierr);
        ierr = reg::Assert(vxi != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = VecCopy(vxi, v->m_X1); CHKERRQ(ierr);
        if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}
        if (regopt->m_Verbosity > 2) {
            ierr = reg::ShowValues(v->m_X1); CHKERRQ(ierr);
        }

        //std::cout << regopt->m_ReadWriteFlags.vx2 << std::endl;
        ierr = readwrite->Read(&vxi, regopt->m_FileNames.iv2); CHKERRQ(ierr);
        ierr = reg::Assert(vxi != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = VecCopy(vxi, v->m_X2); CHKERRQ(ierr);
        if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}
        if (regopt->m_Verbosity > 2) {
            ierr = reg::ShowValues(v->m_X2); CHKERRQ(ierr);
        }

        //std::cout << regopt->m_ReadWriteFlags.vx3 << std::endl;
        ierr = readwrite->Read(&vxi, regopt->m_FileNames.iv3); CHKERRQ(ierr);
        ierr = reg::Assert(vxi != NULL, "null pointer"); CHKERRQ(ierr);
        ierr = VecCopy(vxi, v->m_X3); CHKERRQ(ierr);
        if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}
        if (regopt->m_Verbosity > 2) {
            ierr = reg::ShowValues(v->m_X3); CHKERRQ(ierr);
        }
        ierr = registration->SetInitialGuess(v); CHKERRQ(ierr);
    }

    if (!regopt->m_FileNames.mask.empty()) {
        if (regopt->m_Verbosity > 1) {
            ierr = reg::DbgMsg("reading mask image"); CHKERRQ(ierr);
        }
    }

    ierr = registration->SetReadWrite(readwrite); CHKERRQ(ierr);

    ierr = registration->Run(); CHKERRQ(ierr);

    if (regopt->m_Log.memoryusage) {
        PetscLogDouble mem;
        ierr = PetscMemoryGetMaximumUsage(&mem); CHKERRQ(ierr);
        ss << "memory usage (estimate) " << std::scientific << mem/1E9 << " GB";
        ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }
    // clean up
    if (v != NULL) {delete v; v = NULL;}
    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr); mT = NULL;}
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr); mR = NULL;}
    if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}
    if (regopt != NULL) {delete regopt; regopt = NULL;}
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (registration != NULL) {delete registration; registration = NULL;}

    ierr = reg::Finalize(); CHKERRQ(ierr);

    return 0;
}




