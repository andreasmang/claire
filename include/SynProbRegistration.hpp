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

#ifndef _SYNPROBREGISTRATION_H_
#define _SYNPROBREGISTRATION_H_

#include "RegOpt.hpp"
#include "VecField.hpp"




namespace reg {




class SynProbRegistration {
 public:
    typedef SynProbRegistration Self;
    SynProbRegistration();
    SynProbRegistration(RegOpt*);
    virtual ~SynProbRegistration();

    PetscErrorCode ComputeSquare(Vec);
    PetscErrorCode ComputeExpSin(Vec);
    PetscErrorCode ComputeSphere(Vec);
    PetscErrorCode ComputeHollowSphere(Vec);
    PetscErrorCode ComputeDiamond(Vec, int);
    PetscErrorCode ComputeSmoothScalarField(Vec, int);
    PetscErrorCode ComputeSmoothVectorField(VecField*, int);

 private:
    PetscErrorCode ClearMemory();
    PetscErrorCode Initialize();
    RegOpt* m_Opt;
};




}  // namespace reg




#endif  // _SYNPROBREGISTRATION_H_
