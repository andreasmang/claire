/*************************************************************************
 *  Copyright (c) 2015-2016.
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


#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"
#include "Preconditioner.hpp"
#include "Preprocessing.hpp"
#include "OptimizationProblem.hpp"

namespace reg {


class Optimizer {
 public:
    typedef OptimizationProblem OptProbType;

    Optimizer();
    virtual ~Optimizer();
    Optimizer(RegOpt*);

    PetscErrorCode SetProblem(OptProbType*);

    PetscErrorCode Run(bool presolve=false);
    PetscErrorCode SetInitialGuess(VecField*);
    PetscErrorCode SetPreconditioner(Preconditioner*);
    PetscErrorCode GetSolution(Vec&);
    PetscErrorCode GetSolutionStatus(bool&);

    PetscErrorCode Finalize();

 private:
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode SetupTao(void);
    PetscErrorCode SetInitialGuess(void);

    RegOpt* m_Opt;
    OptProbType* m_OptimizationProblem;

    Tao m_Tao;
//    TaoLineSearch m_LineSearch; ///< line search type
    KSP m_KrylovMethod; ///< KSP object
//    PC m_KrylovMethodPC; ///< KSP preconditioner object
    Preconditioner* m_Precond;
    Preprocessing* m_PreProc;
    Mat m_MatVec;
    Vec m_Solution; ///< solution vector
    bool m_SolutionAllocated;
};




}  // namespace reg




#endif
