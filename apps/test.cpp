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

#include "CLAIREUtils.hpp"
#include "UnitTestOpt.hpp"

/********************************************************************
 * @brief main function to run registration
 *******************************************************************/
int main(int argc, char **argv) {
    PetscErrorCode ierr = 0;
    reg::UnitTestOpt* opt = nullptr;

    // initialize petsc (user is not allowed to set petsc options)
    ierr = PetscInitialize(0, reinterpret_cast<char***>(NULL),
                              reinterpret_cast<char*>(NULL),
                              reinterpret_cast<char*>(NULL)); CHKERRQ(ierr);
    PetscFunctionBegin;

    // allocate class for controlling everything
    ierr = reg::AllocateOnce(opt, argc, argv); CHKERRQ(ierr);
    ierr = opt->DoSetup(); CHKERRQ(ierr);

    ierr = reg::DbgMsg("start unit test"); CHKERRQ(ierr);
    ierr = opt->Run(); CHKERRQ(ierr);

    if (opt != nullptr) {delete opt; opt = nullptr;}

    ierr = reg::Finalize(); CHKERRQ(ierr);

    return 0;
}


