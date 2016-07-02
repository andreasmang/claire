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

#ifndef _TWOLEVELPRECONDREG_H_
#define _TWOLEVELPRECONDREG_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "KrylovInterfaceReg.hpp"


namespace reg
{


class TwoLevelPrecondReg
{

public:

    TwoLevelPrecondReg();
    TwoLevelPrecondReg(RegOpt*);
    ~TwoLevelPrecondReg();


protected:

    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);

    /*! apply two level preconditioner */
    PetscErrorCode Apply(Vec,Vec);


private:

    /*! setup two level preconditioner */
    PetscErrorCode DoSetup();

    /*! apply two level preconditioner */
    PetscErrorCode MatVec(Vec,Vec);

    /*! estimate eigenvalues for hessian */
    PetscErrorCode EstimateEigenValues();

    /*! mat vec for preconditioner */
    PetscErrorCode PCMatVec(Vec,Vec);

    RegOpt* m_Opt;
    KSP m_PCKSP; ///< pointer for KSP method (PETSc)
    Mat m_PCMatVec;

};

} // end of namespace


#endif
