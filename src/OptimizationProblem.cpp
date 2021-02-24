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

#ifndef _OPTIMIZATIONPROBLEM_CPP_
#define _OPTIMIZATIONPROBLEM_CPP_

#include <string>
#include "OptimizationProblem.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
OptimizationProblem::OptimizationProblem() {
    this->Initialize();
}




/********************************************************************
 * @brief default constructor
 *******************************************************************/
OptimizationProblem::OptimizationProblem(RegOpt* opt) {
    this->Initialize();

    this->m_Opt = opt;
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
OptimizationProblem::~OptimizationProblem(void) {
    this->ClearMemory();
}




/********************************************************************
 * @brief init class variables
 *******************************************************************/
PetscErrorCode OptimizationProblem::Initialize(void) {
    PetscFunctionBegin;

    this->m_Opt = NULL;
    this->m_Iterate = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
PetscErrorCode OptimizationProblem::ClearMemory(void) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    if (this->m_Iterate != NULL) {
        ierr = VecDestroy(&this->m_Iterate); CHKERRQ(ierr);
        this->m_Iterate = NULL;
    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the registration options
 *******************************************************************/
PetscErrorCode OptimizationProblem::SetOptions(RegOpt* opt) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // check for null pointer
    ierr = Assert(opt != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Opt != NULL) {
        delete this->m_Opt;
        this->m_Opt = NULL;
    }

    // overwrite
    this->m_Opt = opt;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the registration options
 *******************************************************************/
PetscErrorCode OptimizationProblem::ComputeUpdateNorm(Vec x, ScalarType& normdx, ScalarType& normx) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    // check for null pointer
    ierr = Assert(x != NULL, "null pointer"); CHKERRQ(ierr);

    if (this->m_Iterate == NULL) {
        ierr = VecDuplicate(x, &this->m_Iterate); CHKERRQ(ierr);
        ierr = VecSet(this->m_Iterate, 0.0); CHKERRQ(ierr);
    }

    // compute norm of new iterate
    ierr = VecNorm(x, NORM_INFINITY, &normx); CHKERRQ(ierr);

    // compute update and it's norm
    ierr = VecAXPY(this->m_Iterate, -1.0, x); CHKERRQ(ierr);
    ierr = VecNorm(this->m_Iterate, NORM_INFINITY, &normdx); CHKERRQ(ierr);

    // overwrite iterate
    ierr = VecCopy(x, this->m_Iterate); CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief check gradient based on a taylor expansion
 *******************************************************************/
PetscErrorCode OptimizationProblem::DerivativeCheckGradient() {
    PetscErrorCode ierr = 0;
    Vec v = NULL, vtilde = NULL, w = NULL, g = NULL;
    PetscRandom rctx;
    IntType nl, ng;
    ScalarType h, htilde, Jv, dvJw, Jvtilde, e[2], normv, normw;
    char buffer[256];

    PetscFunctionBegin;

    ierr = Assert(this->m_Opt != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = DbgMsg("performing derivative check (gradient)"); CHKERRQ(ierr);
    snprintf(buffer, 256, "%-12s %-12s %-12s", "h", "e(h)", "e(h^2)");
    ierr = DbgMsg(buffer); CHKERRQ(ierr);

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // create an extra array for initial guess (has to be flat for optimizer)
    //ierr = VecCreate(PETSC_COMM_WORLD, &v); CHKERRQ(ierr);
    //ierr = VecSetSizes(v, 3*nl, 3*ng); CHKERRQ(ierr);
    //ierr = VecSetFromOptions(v); CHKERRQ(ierr);
    ierr = VecCreate(v, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecSet(v, 0.0); CHKERRQ(ierr);

    ierr = VecDuplicate(v, &vtilde); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &g); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &w); CHKERRQ(ierr);

    // create random vectors
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);
    ierr = VecSetRandom(v, rctx); CHKERRQ(ierr);
    ierr = VecSetRandom(w, rctx); CHKERRQ(ierr);

    // compute norm of random vectors
    ierr = VecNorm(v, NORM_2, &normv); CHKERRQ(ierr);
    ierr = VecNorm(w, NORM_2, &normw); CHKERRQ(ierr);

    // normalize random number
    ierr = VecScale(v, 1.0/normv); CHKERRQ(ierr);
    ierr = VecScale(w, 1.0/normw); CHKERRQ(ierr);

    // compute value of objective functional
    ierr = this->EvaluateObjective(&Jv, v); CHKERRQ(ierr);
    ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);

    // do the derivative check
    h = 1E-8;
    for (int i = 0; i < 10; ++i) {
        // compute step size
        htilde = h*pow(10.0, i);

        // perturb velocity field
        ierr = VecWAXPY(vtilde, htilde, w, v); CHKERRQ(ierr);

        // evaluate objective
        ierr = this->EvaluateObjective(&Jvtilde, vtilde); CHKERRQ(ierr);

        // inner product between perturbation and gradient
        ierr = VecTDot(w, g, &dvJw); CHKERRQ(ierr);

        e[0] = (Jvtilde - Jv);
        e[1] = (Jvtilde - Jv - htilde*dvJw);

        e[0] = std::abs(e[0]);
        e[1] = std::abs(e[1]);

        snprintf(buffer, 256, "%e %e %e", htilde, e[0], e[1]);
        ierr = DbgMsg(buffer); CHKERRQ(ierr);
    }

    // clean up
    if (v != NULL) {ierr = VecDestroy(&v); CHKERRQ(ierr); v = NULL;}
    if (w != NULL) {ierr = VecDestroy(&w); CHKERRQ(ierr); w = NULL;}
    if (g != NULL) {ierr = VecDestroy(&g); CHKERRQ(ierr); g = NULL;}
    if (vtilde != NULL) {ierr = VecDestroy(&vtilde); CHKERRQ(ierr); vtilde = NULL;}
    if (rctx != NULL) {ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr); rctx = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief check gradient based on a taylor expansion
 *******************************************************************/
PetscErrorCode OptimizationProblem::DerivativeCheckHessian() {
    PetscErrorCode ierr = 0;
    Vec v = NULL, vtilde = NULL, w = NULL, g = NULL, hvtilde = NULL;
    PetscRandom rctx;
    IntType nl, ng;
    ScalarType h, htilde, Jv, wTdvJ, Jvtilde, wTdvvJw, e[3], normv, normw;
    char buffer[256];

    PetscFunctionBegin;

    ierr = Assert(this->m_Opt != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = DbgMsg("performing derivative check (gradient)"); CHKERRQ(ierr);
    snprintf(buffer, 256, "%-12s %-12s %-12s %-12s", "h", "e(h)", "e(h^2)", "e(h^3)");
    ierr = DbgMsg(buffer); CHKERRQ(ierr);

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // create an extra array for initial guess (has to be flat for optimizer)
    //ierr = VecCreate(PETSC_COMM_WORLD, &v); CHKERRQ(ierr);
    //ierr = VecSetSizes(v, 3*nl, 3*ng); CHKERRQ(ierr);
    //ierr = VecSetFromOptions(v); CHKERRQ(ierr);
    ierr = VecCreate(v, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecSet(v, 0.0); CHKERRQ(ierr);

    ierr = VecDuplicate(v, &vtilde); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &g); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &w); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &hvtilde); CHKERRQ(ierr);

    // create random vectors
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);
    ierr = VecSetRandom(v, rctx); CHKERRQ(ierr);
    ierr = VecSetRandom(w, rctx); CHKERRQ(ierr);

    // compute norm of random vectors
    ierr = VecNorm(v, NORM_2, &normv); CHKERRQ(ierr);
    ierr = VecNorm(w, NORM_2, &normw); CHKERRQ(ierr);

    // normalize random number
    ierr = VecScale(v, 1.0/normv); CHKERRQ(ierr);
    ierr = VecScale(w, 1.0/normw); CHKERRQ(ierr);

    // compute value of objective functional
    ierr = this->EvaluateObjective(&Jv, v); CHKERRQ(ierr);
    ierr = this->EvaluateGradient(g, v); CHKERRQ(ierr);
    ierr = this->HessianMatVec(hvtilde, w); CHKERRQ(ierr);
    ierr = VecTDot(w, hvtilde, &wTdvvJw); CHKERRQ(ierr);

    // inner product between perturbation and gradient
    ierr = VecTDot(w, g, &wTdvJ); CHKERRQ(ierr);

    // do the derivative check
    h = 1E-4;
    for (int i = 0; i < 10; ++i) {
        // compute step size
        htilde = h*pow(10.0, i);

        // perturb velocity field
        ierr = VecWAXPY(vtilde, htilde, w, v); CHKERRQ(ierr);

        // evaluate objective
        ierr = this->EvaluateObjective(&Jvtilde, vtilde); CHKERRQ(ierr);

        e[0] = (Jvtilde - Jv);
        e[1] = (Jvtilde - (Jv + htilde*wTdvJ));
        e[2] = (Jvtilde - (Jv + htilde*wTdvJ + 0.5*htilde*htilde*wTdvvJw));

        e[0] = std::abs(e[0]);
        e[1] = std::abs(e[1]);
        e[2] = std::abs(e[2]);

        snprintf(buffer, 256, "%e %e %e %e", htilde, e[0], e[1], e[2]);
        ierr = DbgMsg(buffer); CHKERRQ(ierr);
    }

    // clean up
    if (v != NULL) {ierr = VecDestroy(&v); CHKERRQ(ierr); v = NULL;}
    if (w != NULL) {ierr = VecDestroy(&w); CHKERRQ(ierr); w = NULL;}
    if (g != NULL) {ierr = VecDestroy(&g); CHKERRQ(ierr); g = NULL;}
    if (vtilde != NULL) {ierr = VecDestroy(&vtilde); CHKERRQ(ierr); vtilde = NULL;}
    if (hvtilde != NULL) {ierr = VecDestroy(&hvtilde); CHKERRQ(ierr); hvtilde = NULL;}
    if (rctx != NULL) {ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr); rctx = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief check gradient based on a taylor expansion
 *******************************************************************/
PetscErrorCode OptimizationProblem::DerivativeCheckHessianFD() {
    PetscErrorCode ierr = 0;
    Vec v = NULL, vtilde = NULL, w = NULL, gv = NULL, gvtilde = NULL, hess = NULL, delta = NULL;
    PetscRandom rctx;
    IntType nl, ng;
    ScalarType h, htilde, temp, epsilon, normv, normw;
    char buffer[256];

    PetscFunctionBegin;

    ierr = Assert(this->m_Opt != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = DbgMsg("performing derivative check (hessian gk - gk+1 = H(xk - xk+1))"); CHKERRQ(ierr);
    snprintf(buffer, 256, "%-12s %-12s", "h", "err");
    ierr = DbgMsg(buffer); CHKERRQ(ierr);

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // create an extra array for initial guess (has to be flat for optimizer)
    ierr = VecCreate(PETSC_COMM_WORLD, &v); CHKERRQ(ierr);
    ierr = VecSetSizes(v, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecSetFromOptions(v); CHKERRQ(ierr);
    ierr = VecSet(v, 0.0); CHKERRQ(ierr);

    ierr = VecDuplicate(v, &vtilde); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &gv); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &gvtilde); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &w); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &hess); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &delta); CHKERRQ(ierr);

    // create random vectors
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);
    ierr = VecSetRandom(v, rctx); CHKERRQ(ierr);
    ierr = VecSetRandom(w, rctx); CHKERRQ(ierr);

    // compute norm of random vectors
    ierr = VecNorm(v, NORM_2, &normv); CHKERRQ(ierr);
    ierr = VecNorm(w, NORM_2, &normw); CHKERRQ(ierr);

    // normalize random number
    ierr = VecScale(v, 1.0/normv); CHKERRQ(ierr);
    ierr = VecScale(w, 1.0/normw); CHKERRQ(ierr);

    // compute value of objective functional
    ierr = this->EvaluateObjective(&temp, v); CHKERRQ(ierr);
    ierr = this->EvaluateGradient(gv, v); CHKERRQ(ierr);

    // do the derivative check
    h = 1E-4;
    for (int i = 0; i < 10; ++i) {
        // compute step size
        htilde = h*pow(10.0,i);

        // perturb velocity field
        ierr = VecWAXPY(vtilde, htilde, w, v); CHKERRQ(ierr);

        // evaluate objective
        ierr = this->EvaluateObjective(&temp, vtilde); CHKERRQ(ierr);
        ierr = this->EvaluateGradient(gvtilde, vtilde); CHKERRQ(ierr);

        ierr = VecCopy(w, delta); CHKERRQ(ierr);
        ierr = VecScale(delta, htilde); CHKERRQ(ierr);
        ierr = this->HessianMatVec(hess, delta); CHKERRQ(ierr);

        // compute finite difference approximation of hessian
        ierr = VecAXPY(gvtilde, -1.0, gv); CHKERRQ(ierr);
        ierr = VecAXPY(gvtilde, -htilde, hess); CHKERRQ(ierr);

        ierr = VecNorm(gvtilde, NORM_2, &epsilon); CHKERRQ(ierr);

        snprintf(buffer, 256, "%e %e", htilde, epsilon);
        ierr = DbgMsg(buffer); CHKERRQ(ierr);
    }

    // clean up
    if (v != NULL) {ierr = VecDestroy(&v); CHKERRQ(ierr); v = NULL;}
    if (w != NULL) {ierr = VecDestroy(&w); CHKERRQ(ierr); w = NULL;}
    if (gv != NULL) {ierr = VecDestroy(&gv); CHKERRQ(ierr); gv = NULL;}
    if (delta != NULL) {ierr = VecDestroy(&delta); CHKERRQ(ierr); delta = NULL;}
    if (gvtilde != NULL) {ierr = VecDestroy(&gvtilde); CHKERRQ(ierr); gvtilde = NULL;}
    if (vtilde != NULL) {ierr = VecDestroy(&vtilde); CHKERRQ(ierr); vtilde = NULL;}
    if (hess != NULL) {ierr = VecDestroy(&hess); CHKERRQ(ierr); hess = NULL;}
    if (rctx != NULL) {ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr); rctx = NULL;}

    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief check symmetry of hessian
 * the idea is to use the identity
 *   \langle A x, A x \rangle = \langle A^T*Ax, x \rangle
 * for the inner product
 *******************************************************************/
PetscErrorCode OptimizationProblem::HessianSymmetryCheck() {
    PetscErrorCode ierr = 0;
    IntType nl, ng;
    Vec v = NULL, Hv = NULL, HHv = NULL;
    ScalarType HvHv, HHvv, symerr, relsymerr, normHv;
    std::string msg;
    std::stringstream sserr, ssrelerr;
    PetscRandom rctx;

    PetscFunctionBegin;

    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;

    // create arrays
    ierr = VecCreate(v, 3*nl, 3*ng); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &Hv); CHKERRQ(ierr);
    ierr = VecDuplicate(v, &HHv); CHKERRQ(ierr);

    // create random vectors
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);
    ierr = VecSetRandom(v, rctx); CHKERRQ(ierr);

    // apply hessian to vector field
    ierr = this->HessianMatVec(Hv, v); CHKERRQ(ierr);
    ierr = this->HessianMatVec(HHv, Hv); CHKERRQ(ierr);

    ierr = VecTDot(Hv, Hv, &HvHv); CHKERRQ(ierr);
    ierr = VecTDot(HHv, v, &HHvv); CHKERRQ(ierr);
    ierr = VecNorm(Hv, NORM_2, &normHv); CHKERRQ(ierr);

    symerr = std::abs(HvHv-HHvv);
    relsymerr = symerr/normHv;

    sserr << symerr;
    ssrelerr << relsymerr;
    msg = "symmetry error of hessian: " + sserr.str()
        + " (relative " + ssrelerr.str() + ")";

    ierr = DbgMsg(msg); CHKERRQ(ierr);

    if (v != NULL) {ierr = VecDestroy(&v); CHKERRQ(ierr);v = NULL;}
    if (Hv != NULL) {ierr = VecDestroy(&Hv); CHKERRQ(ierr); Hv = NULL;}
    if (HHv != NULL) {ierr = VecDestroy(&HHv); CHKERRQ(ierr); HHv = NULL;}
    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _OPTIMIZATIONPROBLEM_CPP_

