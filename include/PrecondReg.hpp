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




namespace reg {




class PrecondReg {
 public:
    typedef OptimizationProblem OptProbType;

    PrecondReg();
    PrecondReg(RegOpt*);
    ~PrecondReg();

    /*! get parameters */
    inline RegOpt* GetOptions(){ return this->m_Opt; };

    /*! set optimization problem */
    PetscErrorCode SetProblem(OptProbType*);

    /*! set preprocessing interface */
    PetscErrorCode SetPreProc(PreProcReg*);

    /*! setup preconditioner */
    PetscErrorCode DoSetup();

    /*! setup preconditioner */
    PetscErrorCode Reset();

    /*! apply preconditioner */
    PetscErrorCode MatVec(Vec, Vec);

    /*! apply hessian (for inversion) */
    PetscErrorCode HessianMatVec(Vec,Vec);

    /*! estimate eigenvalues */
    PetscErrorCode EstimateEigenValues();

 protected:
    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);

 private:
    /*! setup two level preconditioner */
    PetscErrorCode ApplyRestriction();
    PetscErrorCode SetupCoarseGrid();
    PetscErrorCode HessianMatVecProRes(Vec, Vec);
    PetscErrorCode HessianMatVecCoarse(Vec, Vec);

    /*! setup krylov method for inversion of preconditioner */
    PetscErrorCode SetupKrylovMethod();

    /*! setup krylov method for estimating eigenvalues */
    PetscErrorCode SetupKrylovMethodEigEst();

    /*! apply inverse regularization operator as preconditioner */
    PetscErrorCode ApplyInvRegPrecond(Vec, Vec);

    /*! apply 2Level PC as preconditioner */
    PetscErrorCode Apply2LevelPrecond(Vec, Vec);

    /*! apply 2Level PC as preconditioner */
    PetscErrorCode Apply2LevelPrecondResPro(Vec, Vec);


    RegOpt* m_Opt;                      ///< registration options
    RegOpt* m_OptCoarse;                ///< registration options (on coarse grid)
    OptProbType* m_OptProb;             ///< pointer to optimization problem
    OptProbType* m_OptProbCoarse;       ///< pointer to optimization problem (coarse level)

    VecField* m_ControlVariable;        ///< pointer to velocity field
    VecField* m_IncControlVariable;     ///< pointer to velocity field

    Vec m_WorkScaField1;                ///< temporary scalar field
    Vec m_WorkScaField2;                ///< temprary scalar field
    VecField* m_WorkVecField;           ///< temporary vector field
    Vec m_WorkScaFieldCoarse1;          ///< temporary scalar field (coarse level)
    Vec m_WorkScaFieldCoarse2;          ///< temporary scalar field (coarse level)

    Vec m_xCoarse;      ///< array for input to hessian mat vec (on coarse level)
    Vec m_HxCoarse;     ///< array for hessian mat vec (on coarse level)

    Vec m_StateVariableCoarse;              ///< pointer to state variable (coarse level)
    Vec m_AdjointVariableCoarse;            ///< pointer to adjoint variable (coarse level)
    VecField* m_ControlVariableCoarse;      ///< pointer to velocity field (coarse level)
    VecField* m_IncControlVariableCoarse;   ///< pointer to velocity field (coarse level)

    Mat m_MatVec;           ///< mat vec object (PETSc)
    Mat m_MatVecEigEst;     ///< mat vec object (PETSc)

    PreProcReg* m_PreProc;  ///< pointer to preprocessing
    KSP m_KrylovMethod;     ///< pointer for krylov subspace method method (PETSc)

    PetscRandom m_RandomNumGen;     ///< random number generated
    KSP m_KrylovMethodEigEst;

    bool m_SetupKSPOnCoarseGrid;
    bool m_CoarseGridSetupDone;
};




}   // namespace reg




#endif  //_PRECONDREG_H_
