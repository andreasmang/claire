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


#ifndef _PREPROCREG_H_
#define _PREPROCREG_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "ReadWriteReg.hpp"

namespace reg
{

class PreProcReg
{

public:

    typedef PreProcReg Self;
    typedef ScalarType FFTScaType[2];
    typedef ReadWriteReg ReadWriteType;

    PreProcReg();
    PreProcReg(RegOpt*);
    ~PreProcReg();

    PetscErrorCode SetIO(ReadWriteType*);
    PetscErrorCode ApplyGaussianSmoothing(Vec,Vec);

    PetscErrorCode Prolong(Vec,Vec,IntType*);
    PetscErrorCode Restrict(Vec,Vec,IntType*);
    PetscErrorCode Prolong(VecField*,VecField*,IntType*);
    PetscErrorCode Restrict(VecField*,VecField*,IntType*);

private:

    PetscErrorCode ClearMemory();
    PetscErrorCode Initialize();

    RegOpt* m_Opt;
    FFTScaType* m_xhat;
    FFTScaType* m_Kxhat;

    ReadWriteType* m_ReadWrite;


};


} // end of name space


#endif // _PREPROCREG_H_
