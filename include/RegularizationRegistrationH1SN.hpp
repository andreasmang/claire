/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGULARIZATIONREGISTRATIONH1SN_H_
#define _REGULARIZATIONREGISTRATIONH1SN_H_

#include "RegularizationRegistration.hpp"




namespace reg {




class RegularizationRegistrationH1SN : public RegularizationRegistration {
 public:
    typedef RegularizationRegistration SuperClass;
    typedef RegularizationRegistrationH1SN Self;

    RegularizationRegistrationH1SN(void);
    RegularizationRegistrationH1SN(RegOpt*);
    ~RegularizationRegistrationH1SN(void);

    virtual PetscErrorCode EvaluateFunctional(ScalarType*, VecField*);
    virtual PetscErrorCode EvaluateGradient(VecField*, VecField*);
    virtual PetscErrorCode HessianMatVec(VecField*, VecField*);
    virtual PetscErrorCode ApplyInverse(VecField*, VecField*, bool applysqrt = false);
    virtual PetscErrorCode GetExtremeEigValsInvOp(ScalarType&, ScalarType&);
};




}  // end of namespace




#endif
