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
#include "ScaField.hpp"
#include "TransportKernel.hpp"
#include "Differentiation.hpp"


namespace reg {




class TransportProblem {
 public:
    typedef TransportProblem Self;

    TransportProblem();
    TransportProblem(RegOpt*);
    virtual ~TransportProblem();

    PetscErrorCode SetReferenceImage(ScaField*);
    PetscErrorCode SetTemplateImage(ScaField*);

    PetscErrorCode SetStateVariable(ScaField*);
    PetscErrorCode SetAdjointVariable(ScaField*);
    PetscErrorCode SetIncStateVariable(ScaField*);
    PetscErrorCode SetIncAdjointVariable(ScaField*);
    
    PetscErrorCode SetControlVariable(VecField*);
    PetscErrorCode SetIncControlVariable(VecField*);
    
    PetscErrorCode SetWorkScaField(ScaField*, IntType);
    PetscErrorCode SetWorkVecField(VecField*, IntType);
    
    PetscErrorCode SetGradientState(VecField**);
    PetscErrorCode SetGradientXState(VecField**);
    
    virtual PetscErrorCode InitializeControlVariable(VecField*);

    virtual PetscErrorCode SolveForwardProblem();
    virtual PetscErrorCode SolveAdjointProblem();
    virtual PetscErrorCode SolveIncForwardProblem() = 0;
    virtual PetscErrorCode SolveIncAdjointProblem();
    virtual PetscErrorCode SolveInverseProblem();
    
    PetscErrorCode SetDifferentiation(Differentiation::Type);
    PetscErrorCode SetDifferentiation(Differentiation*);
 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    ScaField* m_ReferenceImage;
    ScaField* m_TemplateImage;
    ScaField* m_StateVariable;
    ScaField* m_AdjointVariable;
    ScaField* m_IncStateVariable;
    ScaField* m_IncAdjointVariable;
    
    VecField* m_VelocityField;
    VecField* m_IncVelocityField;
    
    VecField** m_GradientState;
    VecField** m_GradientXState;
    
    ScaField* m_WorkScaField[5];  ///< work scalar field
    VecField* m_WorkVecField[5];  ///< data container for vector field (temporary variable)

    RegOpt* m_Opt;
    
    Differentiation* m_Differentiation;
    bool m_DiffAllocated;

 private:
};




}  // namespace reg




#endif  // _TRANSPORTPROBLEM_HPP_
