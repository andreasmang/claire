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

#ifndef _CLAIREDIVREG_HPP_
#define _CLAIREDIVREG_HPP_

#include "CLAIRE.hpp"




namespace reg {




class CLAIREDivReg : public CLAIRE {
 public:
    typedef CLAIREDivReg Self;
    typedef CLAIRE SuperClass;

    CLAIREDivReg();
    CLAIREDivReg(RegOpt*);
    virtual ~CLAIREDivReg();

 protected:
    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);

    /*! compute body force */
    PetscErrorCode EvaluateObjective(ScalarType*,Vec);

    /*! apply the projection operator to the
        body force and the incremental body force */
    virtual PetscErrorCode ApplyProjection();
    
    virtual PetscErrorCode CreateCoarseReg();
    
 private:
    PetscErrorCode EvaluteRegularizationDIV(ScalarType*);
};




}  // namespace reg




#endif  // _OPTIMALCONTROLREGISTRATIONRELAXEDIC_H_
