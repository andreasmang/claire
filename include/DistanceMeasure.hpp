/*************************************************************************
 *  Copyright (c) 2017.
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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _DISTANCEMEASURE_HPP_
#define _DISTANCEMEASURE_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "VecField.hpp"




namespace reg {




class DistanceMeasure {
 public:
    typedef DistanceMeasure Self;

    DistanceMeasure();
    DistanceMeasure(RegOpt*);
    virtual ~DistanceMeasure();

    /*! evaluate the distance measure (functional) */
    virtual PetscErrorCode EvaluateFunctional(ScalarType*) = 0;

    /*! set the final condition for the adjoint equation */
    virtual PetscErrorCode SetFinalConditionAE() = 0;

    /*! set the final condition for the adjoint equation */
    virtual PetscErrorCode SetFinalConditionIAE() = 0;

    /*! set the reference image */
    PetscErrorCode SetReferenceImage(Vec);

    /*! set the template image */
    PetscErrorCode SetTemplateImage(Vec);

    /*! set the state variable */
    PetscErrorCode SetStateVariable(Vec);

    /*! set the incremental state variable */
    PetscErrorCode SetIncStateVariable(Vec);

    /*! set the adjoint variable */
    PetscErrorCode SetAdjointVariable(Vec);

    /*! set the incremental adjoint variable */
    PetscErrorCode SetIncAdjointVariable(Vec);

    /*! set the auxilary variables */
    PetscErrorCode SetAuxVariable(Vec, int);

    /*! set the mask */
    PetscErrorCode SetMask(Vec);

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    Vec m_Mask;
    Vec m_ReferenceImage;
    Vec m_TemplateImage;

    Vec m_StateVariable;
    Vec m_AdjointVariable;
    Vec m_IncStateVariable;
    Vec m_IncAdjointVariable;

    Vec m_AuxVar1;
    Vec m_AuxVar2;

    RegOpt* m_Opt;
};




}  // namespace reg




#endif
