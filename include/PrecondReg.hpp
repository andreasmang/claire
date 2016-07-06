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
#include "OptimalControlRegistration.hpp"
#include "OptimalControlRegistrationIC.hpp"
#include "OptimalControlRegistrationRelaxedIC.hpp"


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

    /*! get parameters */
    inline RegOpt* GetOptions(){ return this->m_Opt; };

    /*! apply preconditioner */
    PetscErrorCode MatVec(Vec,Vec);

    /*! apply hessian (for inversion) */
    PetscErrorCode HessianMatVec(Vec,Vec);

protected:

    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);

private:

    /*! setup two level preconditioner */
    PetscErrorCode SetUp2LevelPC();

    /*! setup two level preconditioner */
    PetscErrorCode SetupKrylovMethod();

    /*! setup two level preconditioner */
    PetscErrorCode SetTolerancesKrylovMethod();

    PetscErrorCode ApplyInvRegPC(Vec,Vec);
    PetscErrorCode Apply2LevelPC(Vec,Vec);


    RegOpt* m_Opt; ///< registration options
    RegOpt* m_OptCoarse; ///< registration options
    OptProbType* m_OptProb; ///< pointer to optimization problem
    OptProbType* m_OptProbCoarse; ///< pointer to optimization problem (coarse level)

    VecField* m_VelocityField; ///< pointer to velocity field

    Vec m_StateVariableCoarse; ///< pointer to state variable (coarse level)
    Vec m_AdjointVariableCoarse; ///< pointer to adjoint variable (coarse level)
    VecField* m_VelocityFieldCoarse; ///< pointer to velocity field (coarse level)

    Mat m_MatVec; ///< mat vec object (PETSc)



    KSP m_KrylovMethod; ///< pointer for krylov subspace method method (PETSc)

};

} // end of namespace


#endif  //_PRECONDREG_H_
