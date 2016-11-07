/*************************************************************************
 *  Copyright (c) 2016.
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

#ifndef _REGTOOLSOPT_H_
#define _REGTOOLSOPT_H_

#include "RegOpt.hpp"




namespace reg {




struct ResamplingPara {
    ScalarType gridscale;
};


struct RegToolFlags {
    bool computesynprob;     ///< compute synthetic test problem
    bool checkfwdsolve;      ///< perform test for forward solve
    bool readvecfield;       ///< read vector field
    bool resample;           ///< resample scalar / vector field
    bool readscafield;       ///< read scalar field
    bool computedeffields;   ///< compute deformation fields (deformation gradient, displacement field, ...)
    bool computegrad;        ///< compute gradient of scalar field
    bool computesynvel;      ///< compute synthetic velocity field
    bool computeresidual;    ///< compute residual between two images
    bool tscafield;          ///< transport scalar field (forward problem)
    bool tlabelmap;          ///< transport label map (solve forward problem)
    bool convert;            ///< convert image data
    bool applysmoothing;     ///< convert image data
};




class RegToolsOpt : public RegOpt {
 public:
    typedef RegToolsOpt Self;
    typedef RegOpt SuperClass;

    RegToolsOpt();
    RegToolsOpt(int, char**);
    RegToolsOpt(const RegToolsOpt&);
    ~RegToolsOpt();

    std::string GetVecFieldFN(int, int);
    std::string GetScaFieldFN(int);
    virtual PetscErrorCode DisplayOptions(void);
    inline RegToolFlags GetFlags() {return this->m_RegToolFlags;}
    inline ResamplingPara GetResamplingPara() {return this->m_ResamplingPara;}

 protected:
    virtual PetscErrorCode Initialize(void);
    virtual PetscErrorCode ClearMemory(void);
    virtual PetscErrorCode ParseArguments(int, char**);
    virtual PetscErrorCode Usage(bool advanced = false);
    virtual PetscErrorCode CheckArguments(void);

    RegToolFlags m_RegToolFlags;
    ResamplingPara m_ResamplingPara;

    std::string m_RFN;
    std::string m_TFN;

    std::string m_iVecFieldX1FN;    ///< x1 vector field file name
    std::string m_iVecFieldX2FN;    ///< x2 vector field file name
    std::string m_iVecFieldX3FN;    ///< x3 vector field file name
    std::string m_iScaFieldFN;      ///< input file name

    std::string m_xVecFieldX1FN;    ///< x1 vector field file name
    std::string m_xVecFieldX2FN;    ///< x2 vector field file name
    std::string m_xVecFieldX3FN;    ///< x3 vector field file name
    std::string m_xScaFieldFN;      ///< input file name
};




}  // namespace reg




#endif  // _REGTOOLSOPT_H_
