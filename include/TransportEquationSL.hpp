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

#ifndef _TRANSPORTEQUATIONSL_HPP_
#define _TRANSPORTEQUATIONSL_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "TransportProblem.hpp"
#include "SemiLagrangian.hpp"




namespace reg {


class TransportEquationSL : public TransportProblem {
 public:
    typedef TransportEquationSL Self;
    typedef TransportProblem SuperClass;

    TransportEquationSL();
    TransportEquationSL(RegOpt*);
    virtual ~TransportEquationSL();

    PetscErrorCode SolveForwardProblem();
    PetscErrorCode SolveAdjointProblem();
    PetscErrorCode SolveIncForwardProblem();
    PetscErrorCode SolveIncAdjointProblem();
    PetscErrorCode SolveInverseProblem();
    
    PetscErrorCode InitializeControlVariable(VecField*);
    
    SemiLagrangian* GetSemiLagrangian() { return this->m_SemiLagrangianMethod; }

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

 private:
    SemiLagrangian* m_SemiLagrangianMethod;
    
    PetscErrorCode SolveIncAdjointEquationGN();
    PetscErrorCode SolveAdjointEquation();
    PetscErrorCode SolveStateEquation(VecField*);
    
    PetscErrorCode ComputeGradientState(bool, bool, bool);
    
    bool update_grad;
    bool update_gradx;
};




}  // namespace reg




#endif  // _TRANSPORTEQUATION_HPP_
