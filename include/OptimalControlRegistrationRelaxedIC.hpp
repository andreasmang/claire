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

#ifndef _OPTIMALCONTROLREGISTRATIONRELAXEDIC_H_
#define _OPTIMALCONTROLREGISTRATIONRELAXEDIC_H_

#include "OptimalControlRegistration.hpp"




namespace reg {




class OptimalControlRegistrationRelaxedIC : public OptimalControlRegistration {
 public:
    typedef OptimalControlRegistrationRelaxedIC Self;
    typedef OptimalControlRegistration SuperClass;
    typedef ScalarType FFTScaType[2];

    OptimalControlRegistrationRelaxedIC();
    OptimalControlRegistrationRelaxedIC(RegOpt*);
    ~OptimalControlRegistrationRelaxedIC();

 protected:
    /*! init class variables (called by constructor) */
    PetscErrorCode Initialize(void);

    /*! clear memory (called by destructor) */
    PetscErrorCode ClearMemory(void);

    /*! compute body force */
    PetscErrorCode EvaluateObjective(ScalarType*,Vec);

    /*! compute body force */
    PetscErrorCode ComputeBodyForce(void);

    /*! compute incremental body force */
    PetscErrorCode ComputeIncBodyForce(void);

    /*! apply the projection operator to the
        body force and the incremental body force */
    virtual PetscErrorCode ApplyProjection();

 private:
    PetscErrorCode EvaluteRegFunctionalW(ScalarType*);

    FFTScaType *m_x1hat;
    FFTScaType *m_x2hat;
    FFTScaType *m_x3hat;
};




}  // namespace reg




#endif  // _OPTIMALCONTROLREGISTRATIONRELAXEDIC_H_
