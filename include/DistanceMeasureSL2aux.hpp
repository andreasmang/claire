/*************************************************************************
 *  Copyright (c) 2018.
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

#ifndef _DISTANCEMEASURESL2AUX_HPP_
#define _DISTANCEMEASURESL2AUX_HPP_

#include "DistanceMeasure.hpp"




namespace reg {




class DistanceMeasureSL2aux : public DistanceMeasure {
 public:
    typedef DistanceMeasure SuperClass;
    typedef DistanceMeasureSL2aux Self;

    DistanceMeasureSL2aux(void);
    DistanceMeasureSL2aux(RegOpt*);
    virtual ~DistanceMeasureSL2aux(void);

    PetscErrorCode EvaluateFunctional(ScalarType*);
    PetscErrorCode SetFinalCondition();

 protected:
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
};




}  // namespace reg




#endif  // _DISTANCEMEASURESL2AUX_HPP_
