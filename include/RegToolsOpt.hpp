/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CLAIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGTOOLSOPT_H_
#define _REGTOOLSOPT_H_

#include "RegOpt.hpp"




namespace reg {


struct ResamplingPara {
    ScalarType gridscale;
    IntType nx[3];
    bool normalize;
    bool clip;
    ScalarType shift;
    ScalarType scale;
};


struct RegToolFlags {
    bool computedice;           ///< compute DICE
    bool computesynprob;        ///< compute synthetic test problem
    bool readvecfield;          ///< read vector field
    bool resample;              ///< resample scalar / vector field
    bool readscafield;          ///< read scalar field
    bool computedeffields;      ///< compute deformation fields (deformation gradient, displacement field, ...)
    bool computesynvel;         ///< compute synthetic velocity field
    bool computeresidual;       ///< compute residual between two images
    bool computeerror;          ///< compute difference / error between two scalar fields
    bool deformimage;           ///< transport scalar field (forward problem)
    bool tlabelmap;             ///< transport label map (solve forward problem)
    bool tprobmaps;             ///< transport probability map (solve forward problem)
    bool saveprob;              ///< store transported probability maps
    bool convert;               ///< convert image data
    bool computeanalytics;      ///< compute analytics of scalar field
    bool computeravensmap;      ///< compute ravens map
    bool applysmoothing;        ///< convert image data
    bool reference2template;    ///< flag: define direction of deformation map to be inverse (reference space to template space)
};




class RegToolsOpt : public RegOpt {
 public:
    typedef RegToolsOpt Self;
    typedef RegOpt SuperClass;

    RegToolsOpt();
    RegToolsOpt(int, char**);
    RegToolsOpt(const RegToolsOpt&);
    virtual ~RegToolsOpt();

    std::string GetVecFieldFN(int, int);
    std::string GetScaFieldFN(int);
    virtual PetscErrorCode DisplayOptions(void);

//    int m_NumLabels;
    RegToolFlags m_RegToolFlags;
    ResamplingPara m_ResamplingPara;

 protected:
    virtual PetscErrorCode Initialize(void);
    virtual PetscErrorCode ClearMemory(void);
    virtual PetscErrorCode ParseArguments(int, char**);
    virtual PetscErrorCode Usage(bool advanced = false);
    virtual PetscErrorCode CheckArguments(void);
};




}  // namespace reg




#endif  // _REGTOOLSOPT_H_
