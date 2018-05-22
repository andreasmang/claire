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


#ifndef _KRYLOVINTERFACE_H_
#define _KRYLOVINTERFACE_H_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "OptimizationProblem.hpp"
#include "Preconditioner.hpp"




namespace reg {




// mat vec for two level preconditioner
PetscErrorCode KrylovMonitor(KSP,PetscInt,PetscReal,void*);
PetscErrorCode DispKSPConvReason(KSPConvergedReason);

PetscErrorCode InvertPrecondKrylovMonitor(KSP,PetscInt,PetscReal,void*);
PetscErrorCode InvertPrecondMatVec(Mat,Vec,Vec);
PetscErrorCode InvertPrecondPreKrylovSolve(KSP,Vec,Vec,void*);

PetscErrorCode ProjectGradient(KSP,Vec,void*);
PetscErrorCode PreKrylovSolve(KSP,Vec,Vec,void*);
PetscErrorCode PostKrylovSolve(KSP,Vec,Vec,void*);




}  // namespace reg




#endif  // _KRYLOVINTERFACE_HPP_
