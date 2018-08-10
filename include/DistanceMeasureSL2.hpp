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

#ifndef _DISTANCEMEASURESL2_HPP_
#define _DISTANCEMEASURESL2_HPP_

#include "DistanceMeasure.hpp"




namespace reg {




class DistanceMeasureSL2 : public DistanceMeasure {
 public:
    typedef DistanceMeasure SuperClass;
    typedef DistanceMeasureSL2 Self;

    DistanceMeasureSL2(void);
    DistanceMeasureSL2(RegOpt*);
    virtual ~DistanceMeasureSL2(void);

    PetscErrorCode SetupScale();
    PetscErrorCode EvaluateFunctional(ScalarType*);
    PetscErrorCode SetFinalConditionAE();
    PetscErrorCode SetFinalConditionIAE();

 protected:
    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
};




}  // namespace reg




#endif  // _DISTANCEMEASURESL2_HPP_
