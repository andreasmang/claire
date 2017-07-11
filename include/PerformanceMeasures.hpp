/*************************************************************************
 *  Copyright (c) 2017.
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

#ifndef _PERFORMANCEMEASURE_H_
#define _PERFORMANCEMEASURE_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"




namespace reg {


class PerformanceMeasures {
 public:
    typedef PerformanceMeasures Self;
    PerformanceMeasures();
    PerformanceMeasures(RegOpt*);
    ~PerformanceMeasures();

    PetscErrorCode ComputeOverlapMeasures(Vec, Vec);

 private:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    double* m_OverlapMeasures;
    RegOpt* m_Opt;
};




}  // namespace reg




#endif  // _PERFORMANCEMEASURE_H_
