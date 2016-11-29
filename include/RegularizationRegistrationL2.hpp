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
 *  along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGULARIZATIONREGISTRATIONL2_H_
#define _REGULARIZATIONREGISTRATIONL2_H_




#include "RegularizationRegistration.hpp"




namespace reg {




class RegularizationRegistrationL2 : public RegularizationRegistration {
 public:
    typedef RegularizationRegistration SuperClass;
    typedef RegularizationRegistrationL2 Self;

    RegularizationRegistrationL2(void);
    RegularizationRegistrationL2(RegOpt*);
    ~RegularizationRegistrationL2(void);

    PetscErrorCode EvaluateFunctional(ScalarType*, VecField*);
    PetscErrorCode EvaluateGradient(VecField*, VecField*);
    PetscErrorCode HessianMatVec(VecField*, VecField*);
    PetscErrorCode ApplyInvOp(VecField*, VecField*, bool applysqrt = false);
    PetscErrorCode GetExtremeEigValsInvOp(ScalarType&, ScalarType&);
 protected:
 private:
};




}  // namespace reg


#endif  // _REGULARIZATIONREGISTRATIONH2_H_

