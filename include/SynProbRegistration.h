/**
 *  Description: main class for optimal control based registration;
 *  we consider a PDE constrained formulation, where the constraints
 *  are the transport equations for the image intensities; this class
 *  implements a model that is related to the LDDMM formulation; we
 *  invert for a stationary velocity field; the regularization model is
 *  an H2-seminorm
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
 *
*/

#ifndef _SYNPROBREGISTRATION_H_
#define _SYNPROBREGISTRATION_H_

#include "RegOpt.h"




namespace reg
{


class SynProbRegistration
{

public:
    typedef SynProbRegistration Self;
    SynProbRegistration();
    ~SynProbRegistration();
    SynProbRegistration(RegOpt*);

    PetscErrorCode ComputeSquare(Vec);
    PetscErrorCode ComputeDiamond(Vec,const unsigned  int);
    PetscErrorCode ComputeExpSin(Vec);
    PetscErrorCode ComputeSmoothScalarField(Vec,const unsigned int);


private:
    PetscErrorCode ClearMemory();
    PetscErrorCode Initialize();


    RegOpt* m_Opt;

};



} // end of name space




#endif// _SYNPROBREGISTRATION_H_
