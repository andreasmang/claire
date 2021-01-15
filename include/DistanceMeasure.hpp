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
#include "ScaField.hpp"




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
    PetscErrorCode SetReferenceImage(ScaField*);

    /*! set the template image */
    PetscErrorCode SetTemplateImage(ScaField*);

    /*! set the state variable */
    PetscErrorCode SetStateVariable(ScaField*);

    /*! set the incremental state variable */
    PetscErrorCode SetIncStateVariable(ScaField*);

    /*! set the adjoint variable */
    PetscErrorCode SetAdjointVariable(ScaField*);

    /*! set the incremental adjoint variable */
    PetscErrorCode SetIncAdjointVariable(ScaField*);

    /*! set the auxilary variables */
    PetscErrorCode SetAuxVariable(ScaField*, int);

    PetscErrorCode SetWorkVecField(VecField*,int);
    PetscErrorCode SetWorkScaField(ScaField*);

    /*! set the mask */
    PetscErrorCode SetMask(ScaField*);

    /* set objective function weights */
    PetscErrorCode SetObjectiveFunctionalWeights();

    /*! set the scale */
    virtual PetscErrorCode SetupScale();

 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    Vec m_ObjWts;
    ScaField *m_Mask;
    ScaField *m_ReferenceImage;
    ScaField *m_TemplateImage;

    ScaField *m_StateVariable;
    ScaField *m_AdjointVariable;
    ScaField *m_IncStateVariable;
    ScaField *m_IncAdjointVariable;

    ScaField *m_AuxVar1;
    ScaField *m_AuxVar2;

    VecField *m_WorkVecField1;
    VecField *m_WorkVecField2;
    VecField *m_WorkVecField3;
    
    ScaField *m_WorkScaField;

    RegOpt* m_Opt;
};




}  // namespace reg




#endif
