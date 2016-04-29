/**
 *  Description: class to handle input options for registration
 *  Copyright (c) 2015-2016.
 *  All rights reserved.
 *  This file is part of PGLISTR library.
 *
 *  PGLISTR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  PGLISTR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PGLISTR.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#ifndef _REGCOMMANDOPTIONS_H_
#define _REGCOMMANDOPTIONS_H_

#include "RegOpt.h"
#include "RegUtils.h"

namespace reg
{

class RegCommandOptions
{

public:
    RegCommandOptions();
    ~RegCommandOptions();
    RegCommandOptions(int, char**);

    RegOpt* GetRegOpt(){return this->m_RegOpt;};


private:

    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();
    PetscErrorCode ParseArguments();
    PetscErrorCode DoSetup();
    PetscErrorCode Usage();

    RegOpt* m_RegOpt;
    N_MISC* m_MISCOpt;
    accfft_plan* m_ACCFFTPlan;
    MPI_Comm m_MPIComm;

    IntType m_nx[3];
    IntType m_nt;

    int m_CartGridDims[2];
    int m_NumThreads;
    int m_OptMaxIt;
    int m_Verbosity;
    int m_KKTMaxIt;


    std::string m_xfolder;

    ScalarType m_Beta;
    ScalarType m_JacBound;
    ScalarType m_OptTol[3];

    bool m_StoreImages;
    bool m_DoParameterContinuation;
    bool m_Incompressible;
    bool m_StoreTimeSeries;

    PDESolver m_PDESolverType;
    RegNorm m_RegularizationNorm;
    OptMeth m_OptimizationMethod;
    std::string m_mRFileName;
    std::string m_mTFileName;


};

} // end of namespace


#endif // _REGCOMMANDOPTIONS_H_




