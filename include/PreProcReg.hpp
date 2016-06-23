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
    typedef ScalarType ScalarTypeFD[2];
    typedef ReadWriteReg ReadWriteType;

    PreProcReg();
    PreProcReg(RegOpt*);
    ~PreProcReg();

    PetscErrorCode SetReadWrite(ReadWriteType*);
    PetscErrorCode ApplySmoothing(Vec,Vec);

    PetscErrorCode Prolong(Vec*,Vec,IntType*,bool flag=true);
    PetscErrorCode Prolong(VecField*,VecField*,IntType*);

    PetscErrorCode Restrict(Vec*,Vec,IntType*,bool flag=true);
    PetscErrorCode Restrict(VecField*,VecField*,IntType*);


private:

    PetscErrorCode ClearMemory();
    PetscErrorCode Initialize();

    PetscErrorCode SetupRestriction(IntType*);
    PetscErrorCode SetupProlongation(IntType*);

    RegOpt* m_Opt;
    ScalarTypeFD* m_xhat;
    ScalarTypeFD* m_yhat;

    ReadWriteType* m_ReadWrite;

    std::vector<std::vector<IntType>> m_IndicesF;
    std::vector<std::vector<IntType>> m_IndicesC;

};


} // end of name space


#endif // _PREPROCREG_H_
