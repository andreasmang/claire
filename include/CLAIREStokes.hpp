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


#ifndef _CLAIRESTOKES_HPP_
#define _CLAIRESTOKES_HPP_

#include "CLAIRE.hpp"




namespace reg {




class CLAIREStokes : public CLAIRE {
 public:
    typedef CLAIREStokes Self;
    typedef CLAIRE SuperClass;

    CLAIREStokes();
    CLAIREStokes(RegOpt*);
    virtual ~CLAIREStokes();

 protected:
    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);

    /*! reimplementation (accounting for incompressiblity) */
    PetscErrorCode SolveAdjointEquationSL(void);

    /*! reimplementation (accounting for incompressiblity) */
    PetscErrorCode SolveIncAdjointEquationGNSL(void);

    /*! apply the projection operator to the
        body force and the incremental body force */
    virtual PetscErrorCode ApplyProjection();

 private:
};




}  // namespace reg




#endif  // _CLAIRESTOKES_HPP_

