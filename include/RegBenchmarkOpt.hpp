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

#ifndef _REGBENCHMARKOPT_H_
#define _REGBENCHMARKOPT_H_

#include "RegOpt.hpp"




namespace reg {




class RegBenchmarkOpt : public RegOpt {
 public:
    typedef RegBenchmarkOpt Self;
    typedef RegOpt SuperClass;

    RegBenchmarkOpt();
    RegBenchmarkOpt(int, char**);
    RegBenchmarkOpt(const RegBenchmarkOpt&);
    virtual ~RegBenchmarkOpt();

    virtual PetscErrorCode DisplayOptions(void);

 protected:
    virtual PetscErrorCode Initialize(void);
    virtual PetscErrorCode ClearMemory(void);
    virtual PetscErrorCode ParseArguments(int, char**);
    virtual PetscErrorCode Usage(bool advanced = false);
    virtual PetscErrorCode CheckArguments(void);
};




}  // namespace reg




#endif  // _REGBENCHMARKOPT_H_
