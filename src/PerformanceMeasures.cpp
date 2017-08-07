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

#ifndef _PERFORMANCEMEASURES_CPP_
#define _PERFORMANCEMEASURES_CPP_

#include "PerformanceMeasures.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
PerformanceMeasures::PerformanceMeasures() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
PerformanceMeasures::~PerformanceMeasures() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
PerformanceMeasures::PerformanceMeasures(RegOpt* opt) {
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode PerformanceMeasures::Initialize(void) {
    PetscFunctionBegin;

    this->m_Opt = NULL;
    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode PerformanceMeasures::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute overlap between label maps
 * @param mRl label map for reference image
 * @param mTl label map for template image
 *******************************************************************/
PetscErrorCode PerformanceMeasures::ComputeOverlapMeasures(Vec mRl, Vec mTl) {
    PetscErrorCode ierr = 0;
    int nlabels, lR, lT;
    double cj, uj, nlabelsR, nlabelsT, n;
    IntType *icommon = NULL, *iunion = NULL, *nlR = NULL, *nlT = NULL;
    IntType nl;
    ScalarType *p_mrl = NULL, *p_mtl = NULL;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(mRl != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(mTl != NULL, "null pointer"); CHKERRQ(ierr);

    nl = this->m_Opt->m_Domain.nl;

    nlabels = 32;

    if (this->m_OverlapMeasures == NULL) {
        try{this->m_OverlapMeasures = new double [nlabels*4];}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    try {icommon = new IntType[nlabels];}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {iunion = new IntType[nlabels];}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {nlR = new IntType[nlabels];}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {nlT = new IntType[nlabels];}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    for (int i = 0; i < nlabels; ++i) {
        iunion[i] = 0;
        icommon[i] = 0;
        nlR[i] = 0;
        nlT[i] = 0;
    }

    ierr = VecGetArray(mRl, &p_mrl); CHKERRQ(ierr);
    ierr = VecGetArray(mTl, &p_mtl); CHKERRQ(ierr);

    for (IntType i = 0; i < nl; ++i) {
        lR = static_cast<int>(p_mrl[i]);
        lT = static_cast<int>(p_mtl[i]);

        // compute intersection
        if (lR == lT) {
            icommon[lR]++;
        }
        nlR[lR]++;
        nlT[lT]++;

        // compute union
        for (int lj = 0; lj < nlabels; ++lj) {
            if ((lR == lj) || (lT == lj)) iunion[lj]++;
        }
    }

    ierr = VecRestoreArray(mRl, &p_mrl); CHKERRQ(ierr);
    ierr = VecRestoreArray(mTl, &p_mtl); CHKERRQ(ierr);

    for (int lj = 0; lj < nlabels; ++lj) {
        cj = static_cast<double>(icommon[lj]);
        uj = static_cast<double>(iunion[lj]);
        nlabelsR = static_cast<double>(nlR[lj]);
        nlabelsT = static_cast<double>(nlT[lj]);
        n = nlabelsT + nlabelsR;

        // compute jaccard per label
        if (uj != 0.0) {
            this->m_OverlapMeasures[(lj*nlabels)+0] = cj/uj;
        } else {
            this->m_OverlapMeasures[(lj*nlabels)+0] = 0;
        }

        // compute dice per label
        if (n != 0.0) {
            this->m_OverlapMeasures[(lj*nlabels)+1] = 2.0*cj/n;
        } else {
            this->m_OverlapMeasures[(lj*nlabels)+1] = 0.0;
        }

        // compute false positive and false negative per label
        if (nlabelsR != 0.0) {
            this->m_OverlapMeasures[(lj*nlabels)+2] = (nlabelsT-cj)/nlabelsR;
            this->m_OverlapMeasures[(lj*nlabels)+3] = (nlabelsR-cj)/nlabelsR;
        } else {
            this->m_OverlapMeasures[(lj*nlabels)+2] = 0.0;
            this->m_OverlapMeasures[(lj*nlabels)+3] = 0.0;
        }
    }

    if (icommon != NULL) {delete [] icommon; icommon = NULL;}
    if (iunion != NULL) {delete [] iunion; iunion = NULL;}
    if (nlR != NULL) {delete [] nlR; nlR = NULL;}
    if (nlT != NULL) {delete [] nlT; nlT = NULL;}

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg





#endif  //  _PERFORMANCEMEASURES_CPP_
