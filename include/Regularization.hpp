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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGULARIZATION_HPP_
#define _REGULARIZATION_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"
#include "ScaField.hpp"
#include "Differentiation.hpp"




namespace reg {

class Regularization {
 public:
    typedef Regularization Self;

    Regularization(void);
    Regularization(RegOpt*);
    virtual ~Regularization(void);

    PetscErrorCode SetWorkVecField(VecField*);
    PetscErrorCode SetWorkScaField(ScaField*);
    PetscErrorCode SetDifferentiation(Differentiation::Type);
    PetscErrorCode SetDifferentiation(Differentiation*);
    PetscErrorCode SetSpectralData(ComplexType* =nullptr, ComplexType* =nullptr, ComplexType* =nullptr);

    /// TODO Note: The inner product can be computed in spectral space as $\int |Op(v)|^2 dx = \int |\hat{Op(v)}|^2 dx$, thus no inverse FFT needed here
    virtual PetscErrorCode EvaluateFunctional(ScalarType*, VecField*) = 0;
    
    virtual PetscErrorCode EvaluateGradient(VecField*, VecField*) = 0;
    virtual PetscErrorCode HessianMatVec(VecField*, VecField*) = 0;
    virtual PetscErrorCode ApplyInverse(VecField*, VecField*, bool applysqrt = false) = 0;
    virtual PetscErrorCode GetExtremeEigValsInvOp(ScalarType&, ScalarType&) = 0;

 protected:
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    
    RegOpt* m_Opt;
    VecField* m_WorkVecField;
    ScaField* m_WorkScaField;
    
    Differentiation *m_Differentiation;
    bool m_DiffAllocated;

    ComplexType *m_v1hat;
    ComplexType *m_v2hat;
    ComplexType *m_v3hat;
};




}  // namespace reg




#endif // _REGULARIZATION_HPP_

