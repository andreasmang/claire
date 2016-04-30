/**
 *  Copyright (c) 2015-2016.
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
 *
*/


#ifndef _REGULARIZATIONREGISTRATION_H_
#define _REGULARIZATIONREGISTRATION_H_

#include "RegOpt.h"
#include "RegUtils.h"
#include "VecField.h"

namespace reg
{

class RegularizationRegistration
{
public:
    typedef RegularizationRegistration Self;
    typedef ScalarType FFTScaType[2];

    RegularizationRegistration(void);
    RegularizationRegistration(RegOpt*);
    ~RegularizationRegistration(void);

    PetscErrorCode EvaluateFunctional(ScalarType*,VecField*);
    PetscErrorCode EvaluateGradient(VecField*,VecField*);
    PetscErrorCode HessianMatVec(VecField*,VecField*);
    PetscErrorCode ApplyInverseOperator(VecField*,VecField*);

protected:

    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);

    RegOpt* m_Opt;
    VecField* m_WorkVecField;

private:

    PetscErrorCode Allocate();
    PetscErrorCode Deallocate();

    // functions for L2-norm
    PetscErrorCode EvaluateL2(ScalarType*,VecField*);
    PetscErrorCode EvaluateGradL2(VecField*,VecField*);
    PetscErrorCode HessianMatVecL2(VecField*,VecField*);
    PetscErrorCode ApplyInvOpL2(VecField*,VecField*);


    // functions for H1-seminorm
    PetscErrorCode EvaluateH1S(ScalarType*,VecField*);
    PetscErrorCode EvaluateGradH1S(VecField*,VecField*);
    PetscErrorCode HessianMatVecH1S(VecField*,VecField*);
    PetscErrorCode ApplyInvOpH1S(VecField*,VecField*);


    // functions for H2-seminorm
    PetscErrorCode EvaluateH2S(ScalarType*,VecField*);
    PetscErrorCode EvaluateGradH2S(VecField*,VecField*);
    PetscErrorCode HessianMatVecH2S(VecField*,VecField*);
    PetscErrorCode ApplyInvOpH2S(VecField*,VecField*);


    // functions for H1-norm
    PetscErrorCode EvaluateH1(ScalarType*,VecField*);
    PetscErrorCode EvaluateGradH1(VecField*,VecField*);
    PetscErrorCode HessianMatVecH1(VecField*,VecField*);
    PetscErrorCode ApplyInvOpH1(VecField*,VecField*);


    // functions for H2-norm
    PetscErrorCode EvaluateH2(ScalarType*,VecField*);
    PetscErrorCode EvaluateGradH2(VecField*,VecField*);
    PetscErrorCode HessianMatVecH2(VecField*,VecField*);
    PetscErrorCode ApplyInvOpH2(VecField*,VecField*);


    FFTScaType *m_v1hat;
    FFTScaType *m_v2hat;
    FFTScaType *m_v3hat;

    FFTScaType *m_Lv1hat;
    FFTScaType *m_Lv2hat;
    FFTScaType *m_Lv3hat;

};

} // end of namespace


#endif
