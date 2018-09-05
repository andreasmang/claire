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
#include "VecField.hpp"
#include "TransportKernel.hpp"
#include "Differentiation.hpp"


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
    
    PetscErrorCode SetControlVariable(VecField*);
    PetscErrorCode SetIncControlVariable(VecField*);
    
    PetscErrorCode SetWorkScaField(Vec, IntType);
    PetscErrorCode SetWorkVecField(VecField*, IntType);

    virtual PetscErrorCode SolveForwardProblem() = 0;
    virtual PetscErrorCode SolveAdjointProblem() = 0;
    virtual PetscErrorCode SolveIncForwardProblem() = 0;
    virtual PetscErrorCode SolveIncAdjointProblem();
    
    PetscErrorCode SetDifferentiation(Differentiation*);

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    Vec m_ReferenceImage;
    Vec m_TemplateImage;
    Vec m_StateVariable;
    Vec m_AdjointVariable;
    Vec m_IncStateVariable;
    Vec m_IncAdjointVariable;
    
    VecField* m_VelocityField;
    VecField* m_IncVelocityField;
    
    Vec m_WorkScaField[5];  ///< work scalar field
    VecField* m_WorkVecField[5];  ///< data container for vector field (temporary variable)

    RegOpt* m_Opt;
    
    Differentiation* m_Differentiation;

 private:
};




}  // namespace reg




#endif  // _TRANSPORTPROBLEM_HPP_
