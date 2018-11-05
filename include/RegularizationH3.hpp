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

#ifndef _REGULARIZATIONH3_H_
#define _REGULARIZATIONH3_H_

#include "Regularization.hpp"




namespace reg {




class RegularizationH3 : public Regularization {
 public:
    typedef Regularization SuperClass;
    typedef RegularizationH3 Self;

    RegularizationH3(void);
    RegularizationH3(RegOpt*);
    virtual ~RegularizationH3(void);

    virtual PetscErrorCode EvaluateFunctional(ScalarType*, VecField*);
    virtual PetscErrorCode EvaluateGradient(VecField*, VecField*);
    virtual PetscErrorCode HessianMatVec(VecField*, VecField*);
    virtual PetscErrorCode ApplyInverse(VecField*, VecField*, bool applysqrt = false);
    virtual PetscErrorCode GetExtremeEigValsInvOp(ScalarType&, ScalarType&);
};




}  // end of namespace




#endif
