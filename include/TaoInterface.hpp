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

#ifndef _TAOINTERFACE_HPP_
#define _TAOINTERFACE_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "Preconditioner.hpp"
#include "KrylovInterface.hpp"
#include "OptimizationProblem.hpp"




namespace reg {




// the following functions interface the tao solver
PetscErrorCode EvaluateObjective(Tao, Vec, ScalarType*, void*);
PetscErrorCode EvaluateGradient(Tao, Vec, Vec, void*);
PetscErrorCode EvaluateObjectiveGradient(Tao, Vec, ScalarType*, Vec, void*);

PetscErrorCode ConstructHessian(Tao, Vec, Mat*, Mat*, MatStructure*, void*);
PetscErrorCode EvaluateHessian(Tao, Vec, Mat, Mat, void*);
PetscErrorCode HessianMatVec(Mat, Vec, Vec);
PetscErrorCode PrecondMatVec(PC, Vec, Vec);

PetscErrorCode CheckConvergenceGrad(Tao, void*);
PetscErrorCode CheckConvergenceGradObj(Tao, void*);
PetscErrorCode CheckConvergenceGradObjHess(Tao, void*);

PetscErrorCode OptimizationMonitor(Tao, void*);
PetscErrorCode GetLineSearchStatus(Tao, void*);
PetscErrorCode GetSolverStatus(TaoConvergedReason, std::string&);




}  // namespace reg




#endif  // _TAOINTERFACE_HPP_
