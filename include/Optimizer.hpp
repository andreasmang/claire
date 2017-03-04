/**
 * @file Optimizer.hpp
 *
 * Copyright (c) 2015-2016.
 * All rights reserved.
 * This file is part of the XXX library.
 *
 * XXX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * XXX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @brief Generic interface for optimizer to solve registration problem.
 * This class essentially sets up TAO and provides generic functions to
 * solve the optimization problem.
 *
 * @author Andreas Mang
 *
 */



#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "VecField.hpp"
#include "PrecondReg.hpp"
#include "Preprocessing.hpp"
#include "OptimizationProblem.hpp"

namespace reg
{


class Optimizer
{

public:

    typedef OptimizationProblem OptProbType;

    Optimizer();
    ~Optimizer();
    Optimizer(RegOpt*);

    PetscErrorCode SetProblem(OptProbType*);

    PetscErrorCode Run(bool presolve=false);
    PetscErrorCode SetInitialGuess(VecField*);
    PetscErrorCode SetPreconditioner(PrecondReg*);
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
    PrecondReg* m_Precond;
    Preprocessing* m_PreProc;

    Vec m_Solution; ///< solution vector


};

} // end of name space


#endif
