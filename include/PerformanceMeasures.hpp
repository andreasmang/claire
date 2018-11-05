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
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _PERFORMANCEMEASURE_HPP_
#define _PERFORMANCEMEASURE_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"




namespace reg {


class PerformanceMeasures {
 public:
    typedef PerformanceMeasures Self;
    PerformanceMeasures();
    PerformanceMeasures(RegOpt*);
    virtual ~PerformanceMeasures();

    PetscErrorCode ComputeOverlapMeasures(Vec, Vec);

 private:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    double* m_OverlapMeasures;
    RegOpt* m_Opt;
};




}  // namespace reg




#endif  // _PERFORMANCEMEASURE_HPP_
