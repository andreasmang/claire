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

#ifndef _PRECONDREG_H_
#define _PRECONDREG_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "KrylovInterfaceReg.hpp"
#include "OptimizationProblem.hpp"


namespace reg
{


class PrecondReg
{

public:

    typedef OptimizationProblem OptProbType;

    PrecondReg();
    PrecondReg(RegOpt*);
    ~PrecondReg();

    /*! set optimization problem */
    PetscErrorCode SetProblem(OptProbType*);

    /*! apply preconditioner */
    PetscErrorCode MatVec(Vec,Vec);

    /*! apply hessian (for inversion) */
    PetscErrorCode HessianMatVec(Vec,Vec);

    inline RegOpt* GetOptions(){ return this->m_Opt; };

protected:

    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);


private:

    /*! setup two level preconditioner */
    PetscErrorCode Setup2LevelPC();

    /*! setup two level preconditioner */
    PetscErrorCode SetupKrylovMethod();

    /*! setup two level preconditioner */
    PetscErrorCode SetTolerancesKrylovMethod();


    PetscErrorCode ApplyInvRegPC(Vec,Vec);
    PetscErrorCode Apply2LevelPC(Vec,Vec);

    RegOpt* m_Opt; ///< registration options
    KSP m_KrylovMethod; ///< pointer for krylov subspace method method (PETSc)
    Mat m_MatVec; ///< mat vec object (PETSc)
    OptProbType* m_OptimizationProblem; ///< pointer to optimization problem

};

} // end of namespace


#endif  //_PRECONDREG_H_
