/*
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




#ifndef _MULTILEVELPYRAMID_H_
#define _MULTILEVELPYRAMID_H_

#include "RegOpt.hpp"

namespace reg
{


struct DataPyramidNode{
    int m_Level;
    Vec m_Data;
    IntType ng;
    IntType nl;

    struct DataPyramidNode* next;
};



class MultiLevelPyramid
{

public:

    MultiLevelPyramid();
    MultiLevelPyramid(RegOpt*);
    ~MultiLevelPyramid();

    typedef struct DataPyramidNode DataType;
    PetscErrorCode GetLevel(Vec*,int);
    PetscErrorCode SetUp(Vec);
    int GetNumLevels(){return this->m_NumLevels;};
    IntType GetNumGridPoints(int level,int i){return this->m_Domain.nx[level][i];};
    IntType GetNLocal(int level){return this->m_Domain.nlocal[level];};
    IntType GetNGlobal(int level){return this->m_Domain.nglobal[level];};

private:

    PetscErrorCode Initialize();
    PetscErrorCode AllocatePyramid();
    PetscErrorCode Allocate(Vec*,IntType,IntType);
    PetscErrorCode ClearMemory();
    PetscErrorCode SetData(Vec,int);
    PetscErrorCode GetDataPointer(Vec**,int);
    PetscErrorCode ComputeGridSize();

    struct DomainPara
    {
        std::vector<std::vector<IntType>> nx;
        std::vector<std::vector<IntType>> isize;
        std::vector<std::vector<IntType>> istart;
        std::vector<std::vector<IntType>> osize;
        std::vector<std::vector<IntType>> ostart;
        std::vector<IntType> nlocal;
        std::vector<IntType> nallocfd;
        std::vector<IntType> nglobal;
    };

    Vec m_DataL01;
    Vec m_DataL02;
    Vec m_DataL03;
    Vec m_DataL04;
    Vec m_DataL05;
    Vec m_DataL06;
    Vec m_DataL07;
    Vec m_DataL08;
    Vec m_DataL09;
    Vec m_DataL10;
    Vec m_DataL11;
    Vec m_DataL12;
    Vec m_DataL13;
    Vec m_DataL14;
    Vec m_DataL15;

    DomainPara m_Domain;

    RegOpt* m_Opt;

    int m_NumLevels;

};




} // end of namespace




#endif // _MULTILEVELPYRAMID_H_
