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

#ifndef _SYNTHETICDATA_CPP_
#define _SYNTHETICDATA_CPP_

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
 * @brief compute synthetic characteristic point
 * @description Notes that this is a redundant function, need to merge with ComputeSyntheticData()
 *******************************************************************/
PetscErrorCode ComputeSyntheticVelocity(ScalarType* v, ScalarType* x, int vcase) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    ScalarType x1, x2, x3;
    x1 = x[0];
    x2 = x[1];
    x3 = x[2];

    if (vcase == 0) {
        // compute the velocity field
        v[0] = 0.5*PetscSinReal(x3)*PetscCosReal(x2)*PetscSinReal(x2);
        v[1] = 0.5*PetscSinReal(x1)*PetscCosReal(x3)*PetscSinReal(x3);
        v[2] = 0.5*PetscSinReal(x2)*PetscCosReal(x1)*PetscSinReal(x1);
    } else if (vcase == 1) {
        // compute the velocity field
        v[0] = PetscSinReal(x3)*PetscCosReal(x2)*PetscSinReal(x2);
        v[1] = PetscSinReal(x1)*PetscCosReal(x3)*PetscSinReal(x3);
        v[2] = PetscSinReal(x2)*PetscCosReal(x1)*PetscSinReal(x1);
    } else if (vcase == 2) {
        // compute divergence freee velocity field
        v[0] = PetscCosReal(x2)*PetscCosReal(x3);
        v[1] = PetscSinReal(x3)*PetscSinReal(x1);
        v[2] = PetscCosReal(x1)*PetscCosReal(x2);
    } else if (vcase == 3) {
        v[0] = sin(x1);
        v[1] = cos(x2);
        v[2] = sin(x3);
    } else {
      ierr = DebugNotImplemented(); CHKERRQ(ierr);
    }

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
                    p_v1[i] = sin(x1);
                    p_v2[i] = cos(x2);
                    p_v3[i] = sin(x3);
                } else if (vcase == 4) {
                    p_v1[i] = i;
                    p_v2[i] = 0.0;
                    p_v3[i] = 0.0;
                } else if (vcase == 5) {
                    p_v1[i] = i1+i2+i3;
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
    
    using std::exp;
    using std::sin;
    using std::cos;

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
    
    int s = 1;
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
                  pv[0][i]  = sin(M_PI*sin(x1)) *
                              sin(M_PI*sin(x2)) *
                              sin(M_PI*sin(x3));
                  pdv[0][i] = M_PI*cos(x1)*cos(M_PI*sin(x1)) *
                              sin(M_PI*sin(x2)) *
                              sin(M_PI*sin(x3));
                  pdv[1][i] = sin(M_PI*sin(x1)) *
                              M_PI*cos(x2)*cos(M_PI*sin(x2)) *
                              sin(M_PI*sin(x3));
                  pdv[2][i] = sin(M_PI*sin(x1)) *
                              sin(M_PI*sin(x2)) *
                              M_PI*cos(x3)*cos(M_PI*sin(x3));
                  break;
                case 1: // div : vec -> sca
                  pv[0][i]  = sin(M_PI*sin(x1)) * sin(x2) * cos(x3);
                  pv[1][i]  = sin(M_PI*sin(x2)) * sin(x3) * cos(x1);
                  pv[2][i]  = sin(M_PI*sin(x3)) * sin(x1) * cos(x2);
                  pdv[0][i] = M_PI*cos(x1)*cos(M_PI*sin(x1)) * sin(x2) * cos(x3)
                            + M_PI*cos(x2)*cos(M_PI*sin(x2)) * sin(x3) * cos(x1)
                            + M_PI*cos(x3)*cos(M_PI*sin(x3)) * sin(x1) * cos(x2);
                  break;
                case 2: // lap : vec -> vec
                  pv[0][i]  = sin(M_PI*sin(x1)) * sin(10*x2) * cos(4*x3);
                  pv[1][i]  = pv[0][i];
                  pv[2][i]  = pv[0][i];
                  pdv[0][i] = (-M_PI*M_PI*cos(x1)*cos(x1)*sin(M_PI*sin(x1))-M_PI*sin(x1)*cos(M_PI*sin(x1)))*sin(10*x2)*cos(4*x3)
                            - 100*sin(M_PI*sin(x1))*sin(10*x2)*cos(4*x3)
                            - 16*sin(M_PI*sin(x1))*sin(10*x2)*cos(4*x3);
                  pdv[1][i] = pdv[0][i];
                  pdv[2][i] = pdv[0][i];
                  break;
                case 3: // simple function for grad
                  pv[0][i] = i;
                  pdv[0][i] = 0;
                  pdv[1][i] = 0;
                  pdv[2][i] = 0;
                  break;
                case 4:
                  pv[0][i]  = sin(s*x1)*sin(s*x2)*sin(s*x3);
                  pdv[0][i] = s*cos(s*x1)*sin(s*x2)*sin(s*x3);
                  pdv[1][i] = s*sin(s*x1)*cos(s*x2)*sin(s*x3);
                  pdv[2][i] = s*sin(s*x1)*sin(s*x2)*cos(s*x3);
                  break;
                case 5:
                  pv[0][i]  = sin(s*x3);
                  pdv[0][i] = 0;
                  pdv[1][i] = 0;
                  pdv[2][i] = s*cos(s*x3);
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

/********************************************************************
 * @brief compute synthetic image
 *******************************************************************/
PetscErrorCode ComputeGradSpectral(ScalarType w, VecField *v, VecField *dv, reg::RegOpt* opt) {
    PetscErrorCode ierr = 0;
    ScalarType *pv[3];
    ScalarType *pdv[3];
    ScalarType hx[3], x1, x2, x3;
    IntType i;
    PetscFunctionBegin;
    
    using std::exp;
    using std::sin;
    using std::cos;

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
                
                pv[0][i]  = sin(w*x3) + cos(w*x3);
                pv[1][i] = 0;
                pv[2][i] = 0;
                pdv[0][i] = -w*w*sin(w*x3) - w*w*cos(w*x3);
                pdv[1][i] = 0;
                pdv[2][i] = w*cos(w*x3) - w*sin(w*x3);
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

