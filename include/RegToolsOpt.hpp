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


#ifndef _REGTOOLSOPT_H_
#define _REGTOOLSOPT_H_


#include "RegOpt.hpp"

struct ResamplingPara{
    ScalarType gridscale;
    bool enabled;
};


struct PostProcPara{
    bool enabled;
    bool computedeffields;
};


struct RegToolsFlags{
    bool readvecfield;
    bool readscafield;
};


namespace reg
{

class RegToolsOpt : public RegOpt
{

public:

    typedef RegToolsOpt Self;
    typedef RegOpt SuperClass;

    RegToolsOpt();
    RegToolsOpt(int,char**);
    RegToolsOpt(const RegToolsOpt&);
    ~RegToolsOpt();

    std::string GetVecFieldFN(int,int);
    std::string GetScaFieldFN(int);
    inline RegToolsFlags GetFlags(){return this->m_RegToolsFlags;};
    inline ResamplingPara GetResamplingPara(){return this->m_ResamplingPara;};
    inline PostProcPara GetPostProcPara(){return this->m_PostProcPara;};

protected:

    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode ParseArguments(int,char**);
    PetscErrorCode Usage(bool advanced=false);
    PetscErrorCode CheckArguments(void);

    RegToolsFlags m_RegToolsFlags;
    ResamplingPara m_ResamplingPara;
    PostProcPara m_PostProcPara;

    std::string m_iVecFieldX1FN; ///< x1 vector field file name
    std::string m_iVecFieldX2FN; ///< x2 vector field file name
    std::string m_iVecFieldX3FN; ///< x3 vector field file name
    std::string m_iScaFieldFN; ///< input file name

    std::string m_xVecFieldX1FN; ///< x1 vector field file name
    std::string m_xVecFieldX2FN; ///< x2 vector field file name
    std::string m_xVecFieldX3FN; ///< x3 vector field file name
    std::string m_xScaFieldFN; ///< input file name

};


}


#endif // _REGTOOLSOPT_H_
