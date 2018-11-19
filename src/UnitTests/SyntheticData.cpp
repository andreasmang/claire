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
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _TESTINTERPOLATION_CPP_
#define _TESTINTERPOLATION_CPP_

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "UnitTestOpt.hpp"
#include "CLAIRE.hpp"
#include "VecField.hpp"

namespace reg {
namespace UnitTest {
  
/********************************************************************
 * @brief compute synthetic image
 *******************************************************************/
PetscErrorCode ComputeSyntheticData(Vec& m, reg::RegOpt* opt) {
    PetscErrorCode ierr = 0;
    ScalarType *p_m = NULL;
    ScalarType hx[3], x1, x2, x3;
    IntType i, nl, ng;
    PetscFunctionBegin;

    opt->Enter(__func__);

    // get local and global size
    nl = opt->m_Domain.nl;
    ng = opt->m_Domain.ng;

    // get grid size
    hx[0] = opt->m_Domain.hx[0];
    hx[1] = opt->m_Domain.hx[1];
    hx[2] = opt->m_Domain.hx[2];

    // allocate data
    if (m == NULL) {
        ierr = reg::VecCreate(m, nl, ng); CHKERRQ(ierr);
    }

    ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
    for (IntType i1 = 0; i1 < opt->m_Domain.isize[0]; ++i1) {  // x1
        for (IntType i2 = 0; i2 < opt->m_Domain.isize[1]; ++i2) {  // x2
            for (IntType i3 = 0; i3 < opt->m_Domain.isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + opt->m_Domain.istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + opt->m_Domain.istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + opt->m_Domain.istart[2]);

                // compute linear / flat index
                i = reg::GetLinearIndex(i1, i2, i3, opt->m_Domain.isize);

                p_m[i] =  (PetscSinReal(x1)*PetscSinReal(x1)
                          + PetscSinReal(x2)*PetscSinReal(x2)
                          + PetscSinReal(x3)*PetscSinReal(x3))/3.0;


            }  // i1
        }  // i2
    }  // i3
    ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);

    opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute synthetic velocity field
 *******************************************************************/
PetscErrorCode ComputeSyntheticData(reg::VecField*& v, reg::RegOpt* opt, IntType vcase) {
    PetscErrorCode ierr = 0;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL;
    ScalarType hx[3], x1, x2, x3;
    IntType i;
    PetscFunctionBegin;

    opt->Enter(__func__);

    // get grid size
    hx[0] = opt->m_Domain.hx[0];
    hx[1] = opt->m_Domain.hx[1];
    hx[2] = opt->m_Domain.hx[2];

    // allocate velocity field
    if (v == NULL) {
        try {v = new reg::VecField(opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = VecGetArray(v->m_X1, &p_v1); CHKERRQ(ierr);
    ierr = VecGetArray(v->m_X2, &p_v2); CHKERRQ(ierr);
    ierr = VecGetArray(v->m_X3, &p_v3); CHKERRQ(ierr);
    
    for (IntType i1 = 0; i1 < opt->m_Domain.isize[0]; ++i1) {  // x1
        for (IntType i2 = 0; i2 < opt->m_Domain.isize[1]; ++i2) {  // x2
            for (IntType i3 = 0; i3 < opt->m_Domain.isize[2]; ++i3) {  // x3

                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + opt->m_Domain.istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + opt->m_Domain.istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + opt->m_Domain.istart[2]);

                // compute linear / flat index
                i = reg::GetLinearIndex(i1, i2, i3, opt->m_Domain.isize);

                if (vcase == 0) {
                    // compute the velocity field
                    p_v1[i] = 0.5*PetscSinReal(x3)*PetscCosReal(x2)*PetscSinReal(x2);
                    p_v2[i] = 0.5*PetscSinReal(x1)*PetscCosReal(x3)*PetscSinReal(x3);
                    p_v3[i] = 0.5*PetscSinReal(x2)*PetscCosReal(x1)*PetscSinReal(x1);
                } else if (vcase == 1) {
                    // compute the velocity field
                    p_v1[i] = PetscSinReal(x3)*PetscCosReal(x2)*PetscSinReal(x2);
                    p_v2[i] = PetscSinReal(x1)*PetscCosReal(x3)*PetscSinReal(x3);
                    p_v3[i] = PetscSinReal(x2)*PetscCosReal(x1)*PetscSinReal(x1);
                } else if (vcase == 2) {
                    // compute divergence freee velocity field
                    p_v1[i] = PetscCosReal(x2)*PetscCosReal(x3);
                    p_v2[i] = PetscSinReal(x3)*PetscSinReal(x1);
                    p_v3[i] = PetscCosReal(x1)*PetscCosReal(x2);
                } else if (vcase == 3) {
                    p_v1[i] = 0.5;
                    p_v2[i] = 1;
                    p_v3[i] = 0.5;
                } else if (vcase == 4) {
                    p_v1[i] = 0.0;
                    p_v2[i] = 0.0;
                    p_v3[i] = 0.0;
                } else if (vcase == 5) {
                    p_v1[i] = 1.0;
                    p_v2[i] = 1.0;
                    p_v3[i] = 1.0;
                }
            }  // i1
        }  // i2
    }  // i3

    // ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = VecRestoreArray(v->m_X1, &p_v1); CHKERRQ(ierr);
    ierr = VecRestoreArray(v->m_X2, &p_v2); CHKERRQ(ierr);
    ierr = VecRestoreArray(v->m_X3, &p_v3); CHKERRQ(ierr);
    
    opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief compute synthetic image
 *******************************************************************/
PetscErrorCode ComputeDiffFunction(VecField *v, VecField *dv, int type, reg::RegOpt* opt) {
    PetscErrorCode ierr = 0;
    ScalarType *pv[3];
    ScalarType *pdv[3];
    ScalarType hx[3], x1, x2, x3;
    IntType i;
    PetscFunctionBegin;

    opt->Enter(__func__);

    // get grid size
    hx[0] = opt->m_Domain.hx[0];
    hx[1] = opt->m_Domain.hx[1];
    hx[2] = opt->m_Domain.hx[2];
    
    ierr = VecGetArray(v->m_X1, &pv[0]); CHKERRQ(ierr);
    ierr = VecGetArray(v->m_X2, &pv[1]); CHKERRQ(ierr);
    ierr = VecGetArray(v->m_X3, &pv[2]); CHKERRQ(ierr);
    
    ierr = VecGetArray(dv->m_X1, &pdv[0]); CHKERRQ(ierr);
    ierr = VecGetArray(dv->m_X2, &pdv[1]); CHKERRQ(ierr);
    ierr = VecGetArray(dv->m_X3, &pdv[2]); CHKERRQ(ierr);

    for (IntType i1 = 0; i1 < opt->m_Domain.isize[0]; ++i1) {  // x1
        for (IntType i2 = 0; i2 < opt->m_Domain.isize[1]; ++i2) {  // x2
            for (IntType i3 = 0; i3 < opt->m_Domain.isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + opt->m_Domain.istart[0]) - M_PI;
                x2 = hx[1]*static_cast<ScalarType>(i2 + opt->m_Domain.istart[1]) - M_PI;
                x3 = hx[2]*static_cast<ScalarType>(i3 + opt->m_Domain.istart[2]) - M_PI;

                // compute linear / flat index
                i = reg::GetLinearIndex(i1, i2, i3, opt->m_Domain.isize);
                
                switch (type) {
                case 0: // grad : sca -> vec
                  pv[0][i]  = exp(-x1*x1)*cos(x1*0.5) *
                              exp(-x2*x2)*cos(x2*0.5) *
                              exp(-x3*x3)*cos(x3*0.5);
                  pdv[0][i] = -0.5*exp(-x1*x1)*(4*x1*cos(x1*0.5) + sin(x1*0.5)) *
                              exp(-x2*x2)*cos(x2*0.5) *
                              exp(-x3*x3)*cos(x3*0.5);
                  pdv[1][i] = exp(-x1*x1)*cos(x1*0.5) *
                              -0.5*exp(-x2*x2)*(4*x2*cos(x2*0.5) + sin(x2*0.5)) *
                              exp(-x3*x3)*cos(x3*0.5);
                  pdv[2][i] = exp(-x1*x1)*cos(x1*0.5) *
                              exp(-x2*x2)*cos(x2*0.5) *
                              -0.5*exp(-x3*x3)*(4*x3*cos(x3*0.5) + sin(x3*0.5));
                  break;
                case 1: // div : vec -> sca
                  pv[0][i]  = exp(-x1*x1)*cos(x1*0.5) * sin(x2) * cos(x3);
                  pv[1][i]  = exp(-x2*x2)*cos(x2*0.5) * sin(x3) * cos(x1);
                  pv[2][i]  = exp(-x3*x3)*cos(x3*0.5) * sin(x1) * cos(x2);
                  pdv[0][i] = -0.5*exp(-x1*x1)*(4*x1*cos(x1*0.5) + sin(x1*0.5)) * sin(x2) * cos(x3)
                              -0.5*exp(-x2*x2)*(4*x2*cos(x2*0.5) + sin(x2*0.5)) * sin(x3) * cos(x1)
                              -0.5*exp(-x3*x3)*(4*x3*cos(x3*0.5) + sin(x3*0.5)) * sin(x1) * cos(x2);
                  break;
                case 2: // lap : vec -> vec
                  pv[0][i]  = exp(-x1*x1)*cos(x1*0.5) * sin(x2) * cos(x3);
                  pv[1][i]  = pv[0][i];
                  pv[2][i]  = pv[0][i];
                  pdv[0][i] = 0.25*exp(-x1*x1)*((16.*x1*x1-9.)*cos(x1*0.5) + 8.*x1*sin(x1*0.5)) * sin(x2) * cos(x3)
                            - 2.*exp(-x1*x1)*cos(x1*0.5) * sin(x2) * cos(x3);
                  pdv[1][i] = pdv[0][i];
                  pdv[2][i] = pdv[0][i];
                  break;
                };

            }  // i1
        }  // i2
    }  // i3
    
    ierr = VecRestoreArray(v->m_X1, &pv[0]); CHKERRQ(ierr);
    ierr = VecRestoreArray(v->m_X2, &pv[1]); CHKERRQ(ierr);
    ierr = VecRestoreArray(v->m_X3, &pv[2]); CHKERRQ(ierr);
    
    ierr = VecRestoreArray(dv->m_X1, &pdv[0]); CHKERRQ(ierr);
    ierr = VecRestoreArray(dv->m_X2, &pdv[1]); CHKERRQ(ierr);
    ierr = VecRestoreArray(dv->m_X3, &pdv[2]); CHKERRQ(ierr);

    opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

}} // namespace reg

#endif // _TESTINTERPOLATION_CPP_

