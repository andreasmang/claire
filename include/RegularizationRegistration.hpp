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

#ifndef _REGULARIZATIONREGISTRATION_H_
#define _REGULARIZATIONREGISTRATION_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "VecField.hpp"




namespace reg {




class RegularizationRegistration {
 public:
    typedef RegularizationRegistration Self;
    //typedef ScalarType FFTScaType[2];

    RegularizationRegistration(void);
    RegularizationRegistration(RegOpt*);
    virtual ~RegularizationRegistration(void);

    PetscErrorCode SetWorkVecField(VecField*);

    virtual PetscErrorCode EvaluateFunctional(ScalarType*, VecField*) = 0;
    virtual PetscErrorCode EvaluateGradient(VecField*, VecField*) = 0;
    virtual PetscErrorCode HessianMatVec(VecField*, VecField*) = 0;
    virtual PetscErrorCode ApplyInvOp(VecField*, VecField*, bool applysqrt = false) = 0;
    virtual PetscErrorCode GetExtremeEigValsInvOp(ScalarType&, ScalarType&) = 0;


 protected:
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode Allocate();

    RegOpt* m_Opt;
    VecField* m_WorkVecField;

    ComplexType *m_v1hat;
    ComplexType *m_v2hat;
    ComplexType *m_v3hat;
};




}  // namespace reg




#endif
