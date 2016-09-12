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

    PetscErrorCode ComputeOverlapMeasures(Vec,Vec);

    inline void ResetGridChangeOps(bool flag){this->m_ResetGridChangeOps=flag;};
    PetscErrorCode ComputeGridChangeIndices(IntType*,IntType*);

private:

    PetscErrorCode ClearMemory();
    PetscErrorCode Initialize();

    PetscErrorCode GridChangeCommDataRestrict();
    PetscErrorCode GridChangeCommDataProlong();
    PetscErrorCode GridChangeCommIndices();
    PetscErrorCode SetupGridChangeOps(IntType*,IntType*);

    RegOpt* m_Opt;
    ScalarTypeFD* m_xhat;
    ScalarTypeFD* m_yhat;
    ReadWriteType* m_ReadWrite;

    std::vector< std::vector<IntType> > m_IndicesF;
    std::vector< std::vector<IntType> > m_IndicesC;

    ScalarType* m_FourierCoeffSendF;
    ScalarType* m_FourierCoeffSendC;
    ScalarType* m_FourierCoeffRecvF;
    ScalarType* m_FourierCoeffRecvC;

    IntType* m_FourierIndicesSendF;
    IntType* m_FourierIndicesSendC;
    IntType* m_FourierIndicesRecvF;
    IntType* m_FourierIndicesRecvC;

    accfft_plan* m_FFTFinePlan;
    accfft_plan* m_FFTCoarsePlan;

    ScalarTypeFD* m_XHatFine;
    ScalarTypeFD* m_XHatCoarse;

    IntType m_nxC[3];
    IntType m_nxF[3];
    IntType m_osizeC[3];
    IntType m_osizeF[3];
    IntType m_ostartC[3];
    IntType m_ostartF[3];

/*
    IntType *m_NumSend;
    IntType *m_NumRecv;
*/
    IntType *m_NumSend;
    IntType *m_NumRecv;
    IntType *m_OffsetSend;
    IntType *m_OffsetRecv;
    IntType m_nAllocSend;
    IntType m_nAllocRecv;

    ScalarType m_FFTFineScale;
    ScalarType m_FFTCoarseScale;
    MPI_Request *m_SendRequest;
    MPI_Request *m_RecvRequest;

    bool m_GridChangeOpsSet;
    bool m_ResetGridChangeOps;
    bool m_IndicesCommunicated;
    bool m_GridChangeIndicesComputed;


    double *m_OverlapMeasures;

};


} // end of name space


#endif // _PREPROCREG_H_
