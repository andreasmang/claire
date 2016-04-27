/**
 *  Description: base class for Optimal Control Registration
 *  Copyright (c) 2015-2016.
 *  All rights reserved.
 *  This file is part of PGLISTR library.
 *
 *  PGLISTR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  PGLISTR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PGLISTR.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "RegOpt.h"
#include "RegUtils.h"
#include "VecField.h"
#include "OptimizationProblemRegistration.h"

namespace reg
{

class Optimizer
{

public:
    typedef OptimizationProblemRegistration OptProbType;

    Optimizer();
    ~Optimizer();
    Optimizer(RegOpt*);

    PetscErrorCode SetProblem(OptProbType*);


    PetscErrorCode Run();
    PetscErrorCode RunBetaCont();
    PetscErrorCode GetSolution(Vec&);
    PetscErrorCode Finalize();


private:

    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode DoSetup(void);

    RegOpt* m_Opt;
    OptProbType* m_OptimizationProblem;

    Tao m_Tao;
    TaoLineSearch m_LineSearch; ///< line search type
    KSP m_KSP; ///< KSP
    PC m_KKTPC; ///< KSP preconditioner



};

} // end of name space

#endif
