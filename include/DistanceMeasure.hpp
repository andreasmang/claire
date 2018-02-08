/*************************************************************************
 *  Copyright (c) 2017.
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

#ifndef _DISTANCEMEASURE_HPP_
#define _DISTANCEMEASURE_HPP_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "VecField.hpp"




namespace reg {




class DistanceMeasure {
 public:
    typedef DistanceMeasure Self;

    DistanceMeasure();
    DistanceMeasure(RegOpt*);
    virtual ~DistanceMeasure();

    virtual PetscErrorCode EvaluateFunctional(ScalarType*) = 0;
    virtual PetscErrorCode SetFinalCondition(Vec) = 0;

    PetscErrorCode SetReferenceImage(Vec);
    PetscErrorCode SetTemplateImage(Vec);
    PetscErrorCode SetStateVariable(Vec);
    PetscErrorCode SetMask(Vec);

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    Vec m_Mask;
    Vec m_ReferenceImage;
    Vec m_TemplateImage;
    Vec m_StateVariable;

    RegOpt* m_Opt;
};




}  // namespace reg




#endif
