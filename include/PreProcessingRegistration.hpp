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


#ifndef _PREPROCESSINGREGISTRATION_H_
#define _PREPROCESSINGREGISTRATION_H_

#include "RegOpt.hpp"
#include "RegUtils.hpp"
#include "ReadWriteReg.hpp"

namespace reg
{

class PreProcessingRegistration
{

public:

    typedef PreProcessingRegistration Self;
    typedef ScalarType FFTScaType[2];
    typedef ReadWriteReg ReadWriteType;

    PreProcessingRegistration();
    PreProcessingRegistration(RegOpt*);
    ~PreProcessingRegistration();

    PetscErrorCode ApplyGaussianSmoothing(Vec,Vec);
    PetscErrorCode ApplyMollifier(Vec,Vec);
    PetscErrorCode CompositionData(Vec,unsigned int,std::string);
    PetscErrorCode CompositionTimeDependentData(Vec,unsigned int,std::string);
    PetscErrorCode DecompositionData(Vec,unsigned int,std::string);
    PetscErrorCode DecompositionTimeDependentData(Vec,unsigned int,std::string);

    PetscErrorCode SetIO(ReadWriteType*);

private:

    PetscErrorCode ClearMemory();
    PetscErrorCode Initialize();
    PetscErrorCode ResetDDData(int);

    PetscErrorCode SetupDomainComposition(unsigned int);
    PetscErrorCode SetupDomainDecomposition(unsigned int);

    inline unsigned long GetIndex(int i1, int i2, int i3){

        // enforce periodic boundary conditions
        i1 = (i1 + static_cast<int>(this->m_nx[0])) % static_cast<int>(this->m_nx[0]);
        i2 = (i2 + static_cast<int>(this->m_nx[1])) % static_cast<int>(this->m_nx[1]);
        i3 = (i3 + static_cast<int>(this->m_nx[2])) % static_cast<int>(this->m_nx[2]);

        return GetLinearIndex(i1,i2,i3,this->m_nx);

    }

    RegOpt* m_Opt;
    FFTScaType* m_xhat;
    FFTScaType* m_Kxhat;

    // parameters for domain decomposition
    struct DomainDec{
        unsigned int* isize;
        int* istart;
        unsigned int* iend;
        unsigned int nsubdom; // number of sub domains
        unsigned int nshared; // number of shared grid points
        unsigned int nzeropad; // number of shared grid points
        IntType nglobal; // global number of grid points
        unsigned long* nlocal; // local number of grid ponts (per subdomain)
    };

    DomainDec m_DDPara;
    ReadWriteType* m_IO;
    unsigned int m_nx[3];

};


} // end of name space


#endif // _PREPROCESSINGREGISTRATION_H_
