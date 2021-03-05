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

#ifndef _PRECONDITIONER_HPP_
#define _PRECONDITIONER_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "KrylovInterface.hpp"
#include "OptimizationProblem.hpp"
#include "CLAIRE.hpp"
#include "CLAIREStokes.hpp"
#include "CLAIREDivReg.hpp"




namespace reg {




class Preconditioner {
 public:
    typedef OptimizationProblem OptProbType;

    Preconditioner();
    Preconditioner(RegOpt*);
    virtual ~Preconditioner();

    /*! get parameters */
    inline RegOpt* GetOptions(){ return this->m_Opt; };

    /*! set optimization problem */
    PetscErrorCode SetProblem(OptProbType*);

    /*! set preprocessing interface */
    PetscErrorCode SetPreProc(Preprocessing*);

    /*! setup preconditioner */
    PetscErrorCode DoSetup();

    /*! setup preconditioner */
    PetscErrorCode Reset();

    /*! apply preconditioner */
    PetscErrorCode MatVec(Vec, Vec);

    /*! apply hessian (for inversion) */
    PetscErrorCode HessianMatVec(Vec, Vec);

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

    /*! setup krylov method for inversion of preconditioner */
    PetscErrorCode SetupKrylovMethod(IntType, IntType);

    /*! setup krylov method for estimating eigenvalues */
    PetscErrorCode SetupKrylovMethodEigEst();

    /*! apply inverse regularization operator as preconditioner */
    PetscErrorCode ApplySpectralPrecond(Vec, Vec);

    /*! apply 2Level PC as preconditioner */
    PetscErrorCode Apply2LevelPrecond(Vec, Vec);
    
    /*! apply inverse of H(v=0) as preconditioner */
    PetscErrorCode ApplyH0Precond(Vec, Vec, bool);

    struct CoarseGrid {
        RegOpt* m_Opt;                        ///< registration options (on coarse grid)
        OptProbType* m_OptimizationProblem;   ///< pointer to optimization problem (on coarse level)
        Vec x;                                ///< array for input to hessian mat vec (on coarse level)
        Vec y;                                ///< array for hessian mat vec (on coarse level)
        Vec m_StateVariable;                  ///< pointer to state variable (on coarse level)
        Vec m_AdjointVariable;                ///< pointer to adjoint variable (on coarse level)
        Vec m_ReferenceImage;                 ///< pointer to adjoint variable (on coarse level)
        Vec m_Mask;                           ///< on coarse level
        VecField* m_ControlVariable;          ///< pointer to velocity field (on coarse level)
        VecField* m_IncControlVariable;       ///< pointer to velocity field (on coarse level)
        Vec m_WorkScaField1;                  ///< temporary scalar field
        Vec m_WorkScaField2;                  ///< temprary scalar field

        inline IntType nl(){return this->m_Opt->m_Domain.nl;};
        inline IntType ng(){return this->m_Opt->m_Domain.ng;};
        bool setupdone;
    };
    
    bool firstrun;                          ///< flag to indicate first run of preconditioner

    CoarseGrid* m_CoarseGrid;

    RegOpt* m_Opt;                          ///< registration options
    OptProbType* m_OptimizationProblem;     ///< pointer to optimization problem

    VecField* m_ControlVariable;            ///< pointer to velocity field
    VecField* m_IncControlVariable;         ///< pointer to velocity field

    Vec m_Mask;                             ///< mask (objective masking)
    Vec m_ReferenceImage;                   ///< reference image
    Vec m_WorkScaField1;                    ///< temporary scalar field
    Vec m_WorkScaField2;                    ///< temprary scalar field
    VecField* m_WorkVecField;               ///< temporary vector field
    
    VecField* m_GradState;                  ///< gradient of state Variable

    Mat m_MatVec;                           ///< mat vec object (PETSc)
    Mat m_MatVecEigEst;                     ///< mat vec object (PETSc)

    Preprocessing* m_PreProc;               ///< pointer to preprocessing
    KSP m_KrylovMethod;                     ///< pointer for krylov subspace method method (PETSc)

    PetscRandom m_RandomNumGen;             ///< random number generated
    KSP m_KrylovMethodEigEst;

};




}   // namespace reg




#endif  //_PRECONDREG_H_
