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
    PetscErrorCode ApplyRectFreqFilter(Vec,Vec,ScalarType,bool flag=true);
    PetscErrorCode ApplyRectFreqFilter(VecField*,VecField*,ScalarType,bool flag=true);

    PetscErrorCode Prolong(Vec*,Vec,IntType*,IntType*);
    PetscErrorCode Prolong(VecField*,VecField*,IntType*,IntType*);

    PetscErrorCode Restrict(Vec*,Vec,IntType*,IntType*);
    PetscErrorCode Restrict(VecField*,VecField*,IntType*,IntType*);

    inline void ResetGridChangeOperators(bool flag){this->m_ResetGridChangeOperators=flag;};

    PetscErrorCode ComputeIndices(IntType*,IntType*);

private:

    PetscErrorCode ClearMemory();
    PetscErrorCode Initialize();

    PetscErrorCode SetupGridChangeOperators(IntType*,IntType*);
    PetscErrorCode CommunicateData();

    RegOpt* m_Opt;
    ScalarTypeFD* m_xhat;
    ScalarTypeFD* m_yhat;
    ReadWriteType* m_ReadWrite;

    std::vector< std::vector<IntType> > m_IndicesF;
    std::vector< std::vector<IntType> > m_IndicesC;
    std::vector< std::vector<ScalarTypeFD> > m_ValuesF;
    std::vector< std::vector<ScalarTypeFD> > m_ValuesC;

    accfft_plan* m_FFTFinePlan;
    accfft_plan* m_FFTCoarsePlan;

    ScalarTypeFD* m_XHatFine;
    ScalarTypeFD* m_XHatCoarse;

    IntType m_osize_c[3];
    IntType m_osize_f[3];
    IntType m_ostart_c[3];
    IntType m_ostart_f[3];

    ScalarType m_FFTFineScale;
    ScalarType m_FFTCoarseScale;

    bool m_ResetGridChangeOperators;
    bool m_IndicesComputed;
};


} // end of name space


#endif // _PREPROCREG_H_
