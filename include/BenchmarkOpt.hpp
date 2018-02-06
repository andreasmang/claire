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

#ifndef _BENCHMARKOPT_H_
#define _BENCHMARKOPT_H_

#include "RegOpt.hpp"




namespace reg {




class BenchmarkOpt : public RegOpt {
 public:
    typedef BenchmarkOpt Self;
    typedef RegOpt SuperClass;

    BenchmarkOpt();
    BenchmarkOpt(int, char**);
    BenchmarkOpt(const BenchmarkOpt&);
    virtual ~BenchmarkOpt();

    virtual PetscErrorCode DisplayOptions(void);
    inline int Benchmark() const {return this->m_BenchmarkID;};
    inline int NumRepeats() const {return this->m_NumRepeats;};

    inline double GetRunTime() const {return this->m_RunTime;};
    inline void SetRunTime(double value){this->m_RunTime = value;};

 protected:
    virtual PetscErrorCode Initialize(void);
    virtual PetscErrorCode ClearMemory(void);
    virtual PetscErrorCode ParseArguments(int, char**);
    virtual PetscErrorCode Usage(bool advanced = false);
    virtual PetscErrorCode CheckArguments(void);

    int m_BenchmarkID;
    int m_NumRepeats;
    double m_RunTime;
};




}  // namespace reg




#endif  // _REGBENCHMARKOPT_H_
