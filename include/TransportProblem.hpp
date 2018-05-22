/*************************************************************************
 *  Copyright (c) 2018.
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

#ifndef _TRANSPORTPROBLEM_HPP_
#define _TRANSPORTPROBLEM_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"




namespace reg {




class TransportProblem {
 public:
    typedef TransportProblem Self;

    TransportProblem();
    TransportProblem(RegOpt*);
    virtual ~TransportProblem();

    PetscErrorCode SetReferenceImage(Vec);
    PetscErrorCode SetTemplateImage(Vec);

    PetscErrorCode SetStateVariable(Vec);
    PetscErrorCode SetAdjointVariable(Vec);
    PetscErrorCode SetIncStateVariable(Vec);
    PetscErrorCode SetIncAdjointVariable(Vec);

    virtual PetscErrorCode SolveForwardProblem() = 0;
    virtual PetscErrorCode SolveAdjointProblem() = 0;
    virtual PetscErrorCode SolveIncForwardProblem() = 0;
    virtual PetscErrorCode SolveIncAdjointProblem() = 0;

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    Vec m_ReferenceImage;
    Vec m_TemplateImage;
    Vec m_StateVariable;
    Vec m_AdjointVariable;
    Vec m_IncStateVariable;
    Vec m_IncAdjointVariable;

    RegOpt* m_Opt;

 private:
};




}  // namespace reg




#endif  // _TRANSPORTPROBLEM_HPP_
