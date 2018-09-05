/*************************************************************************
 *  Copyright (c) 2018.
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

#ifndef _TRANSPORTEQUATIONSL_CPP_
#define _TRANSPORTEQUATIONSL_CPP_

#include "TransportEquationSL.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
TransportEquationSL::TransportEquationSL() : SuperClass() {
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
TransportEquationSL::~TransportEquationSL() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
TransportEquationSL::TransportEquationSL(RegOpt* opt) : SuperClass(opt) {
    this->Initialize();
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
PetscErrorCode TransportEquationSL::Initialize() {
    PetscFunctionBegin;

    this->m_TemplateImage = NULL;
    this->m_ReferenceImage = NULL;
    
    this->m_SemiLagrangianMethod = nullptr;
    this->m_Differentiation = nullptr;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode TransportEquationSL::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;
    
    if (this->m_SemiLagrangianMethod != nullptr) {
        delete this->m_SemiLagrangianMethod;
        this->m_SemiLagrangianMethod = nullptr;
    }

    PetscFunctionReturn(ierr);
}



/********************************************************************
 * @brief solve the forward problem (i.e., the continuity equation)
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveForwardProblem() {
    PetscErrorCode ierr = 0;
    IntType nl, nc, nt, l, lnext;
    ScalarType *p_m = NULL;
    bool store = true;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);

    // flag to identify if we store the time history
    store = this->m_Opt->m_RegFlags.runinversion;

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;

    if (this->m_SemiLagrangianMethod == nullptr) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangian(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // compute trajectory
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);

    // get state variable m
    ierr = GetRawPointerReadWrite(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    for (IntType j = 0; j < nt; ++j) {  // for all time points
        if (store) {
            l = j*nl*nc; lnext = (j+1)*nl*nc;
        } else {
            l = 0; lnext = 0;
        }
        for (IntType k = 0; k < nc; ++k) {  // for all image components
            // compute m(X,t^{j+1}) (interpolate state variable)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_m + lnext + k*nl, p_m + l + k*nl, "state"); CHKERRQ(ierr);
        }
    }
    ierr = RestoreRawPointerReadWrite(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveAdjointProblem() {
    PetscErrorCode ierr = 0;
    ScalarType *p_v1 = nullptr, *p_v2 = nullptr, *p_v3 = nullptr,
                *p_l = nullptr, *p_m = nullptr;
    IntType nl, ng, nc, nt, ll, lm, llnext;
    TransportKernelAdjointSL kernel;
    bool fullnewton = false;

    PetscFunctionBegin;
    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    kernel.scale = kernel.ht/static_cast<ScalarType>(nc);
    kernel.nl = nl;

    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[2] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);
  
    if (this->m_SemiLagrangianMethod == nullptr) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangian(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // compute trajectory for adjoint equations
    ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);

    // for full newton we store the adjoint variable
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        fullnewton = true;
    }
    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField[0], &kernel.pDivV); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField[1], &kernel.pDivVx); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField[2], &kernel.pLx); CHKERRQ(ierr);

    // compute divergence of velocity field
    ierr = this->m_Differentiation->Divergence(kernel.pDivV,  this->m_VelocityField); CHKERRQ(ierr);
        
    // compute divergence of velocity field at X
    ierr = this->m_WorkVecField[0]->GetArrays(kernel.pGm); CHKERRQ(ierr);

    // evaluate div(v) along characteristic X
    ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pDivVx, kernel.pDivV, "adjoint"); CHKERRQ(ierr);
    
    // init body force for numerical integration
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArrays(kernel.pB); CHKERRQ(ierr);

    // perform numerical time integration for adjoint variable and
    // add up body force
    for (IntType j = 0; j < nt; ++j) {
        lm = (nt-j)*nc*nl;
        if (fullnewton) {
            ll = (nt-j)*nc*nl; llnext = (nt-(j+1))*nc*nl;
        } else {
            ll = 0; llnext = 0;
        }

        // scaling for trapezoidal rule (for body force)
        if (j == 0) kernel.scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {
            kernel.pL = p_l + ll + k*nl;
            kernel.pLnext = p_l + llnext + k*nl;
          
            // compute lambda(t^j,X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pLx, kernel.pL, "adjoint"); CHKERRQ(ierr);
            
            // compute gradient of m (for incremental body force)
            ierr = this->m_Differentiation->Gradient(kernel.pGm, p_m + lm + k*nl); CHKERRQ(ierr);
            
            // compute \lambda(x,t^{j+1}) and bodyforce
            ierr = kernel.ComputeBodyForcePart1(); CHKERRQ(ierr);
        }
        // trapezoidal rule (revert scaling; for body force)
        if (j == 0) kernel.scale *= 2.0;
    }
    
    kernel.scale *= 0.5;
    
    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
      lm = k*nl;
      kernel.pL = p_l + k*nl;

      // compute gradient of m (for incremental body force)
      ierr = this->m_Differentiation->Gradient(kernel.pGm, p_m + lm); CHKERRQ(ierr);
      
      ierr = kernel.ComputeBodyForcePart2(); CHKERRQ(ierr);
    }
        
    ierr = this->m_WorkVecField[1]->RestoreArrays(kernel.pB); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->RestoreArrays(kernel.pGm); CHKERRQ(ierr);

    ierr = RestoreRawPointer(this->m_WorkScaField[2], &kernel.pLx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField[1], &kernel.pDivVx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField[0], &kernel.pDivV); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_AdjointVariable, &p_l); CHKERRQ(ierr);
  
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveIncForwardProblem() {
    PetscErrorCode ierr = 0;
    IntType nl, ng, nt, nc, lm, lmnext, lmt, lmtnext;
    ScalarType *p_mtilde = nullptr, *p_m = nullptr, *p_mx = nullptr;
    TransportKernelIncStateSL kernel;
    bool fullnewton = false;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    kernel.hthalf = 0.5*this->m_Opt->GetTimeStepSize();
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncVelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);

    if (this->m_SemiLagrangianMethod == nullptr) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangian(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "state"); CHKERRQ(ierr);
    }

    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {   // gauss newton
        fullnewton = true;
    }

    ierr = GetRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField[0], &p_mx); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArrays(kernel.pGm[0], kernel.pGm[1], kernel.pGm[2]); CHKERRQ(ierr);

    ierr = this->m_SemiLagrangianMethod->Interpolate(this->m_WorkVecField[1], this->m_IncVelocityField, "state"); CHKERRQ(ierr);

    ierr = this->m_WorkVecField[1]->GetArraysRead(kernel.pVtildex); CHKERRQ(ierr);
    ierr = this->m_IncVelocityField->GetArraysRead(kernel.pVtilde); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j) {  // for all time points
        lm = j*nl*nc; lmnext = (j+1)*nl*nc;
        if (fullnewton) {   // full newton
            lmt = j*nl*nc; lmtnext = (j+1)*nl*nc;
        } else {
            lmt = 0; lmtnext = 0;
        }

        for (IntType k = 0; k < nc; ++k) {  // for all image components
            kernel.pMtilde = p_mtilde + lmtnext + k*nl;
            
            // interpolate incremental adjoint variable \tilde{m}^j(X)
            ierr = this->m_SemiLagrangianMethod->Interpolate(p_mtilde + lmtnext + k*nl, p_mtilde + lmt + k*nl, "state"); CHKERRQ(ierr);

            // compute gradient for state variable
            this->m_Differentiation->Gradient(kernel.pGm,p_m+lm+k*nl);

            ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pGm[0], kernel.pGm[1], kernel.pGm[2], kernel.pGm[0], kernel.pGm[1], kernel.pGm[2], "state"); CHKERRQ(ierr);
            
            DBGCHK();
            
            // first part of time integration
            ierr = kernel.TimeIntegrationPart1(); CHKERRQ(ierr);

            DBGCHK();
            
            // compute gradient for state variable at next time time point
            this->m_Differentiation->Gradient(kernel.pGm, p_m+lmnext+k*nl);

            // second part of time integration
            ierr = kernel.TimeIntegrationPart2(); CHKERRQ(ierr);
        }
    }  // for all time points

    ierr = this->m_IncVelocityField->RestoreArraysRead(kernel.pVtilde); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->RestoreArraysRead(kernel.pVtildex); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->RestoreArrays(kernel.pGm); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField[0], &p_mx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_IncStateVariable, &p_mtilde); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &p_m); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve the adjoint problem
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveIncAdjointProblem() {
    PetscErrorCode ierr = 0;
    bool VelocityIsZero = false;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    
    if (this->m_Opt->m_OptPara.method == GAUSSNEWTON) {   // gauss newton
        ierr = this->m_VelocityField->IsZero(VelocityIsZero); CHKERRQ(ierr);
        if (VelocityIsZero) {
            ierr = SuperClass::SolveIncAdjointProblem(); CHKERRQ(ierr);
        } else {
            ierr = this->SolveIncAdjointEquationGN(); CHKERRQ(ierr);
        }
    } else if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ierr = ThrowError("PDE solver not implemented"); CHKERRQ(ierr);
    } else {
        ierr = ThrowError("update method not defined"); CHKERRQ(ierr);
    }
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

/********************************************************************
 * @brief solve the adjoint problem with Gauss-Newton
 *******************************************************************/
PetscErrorCode TransportEquationSL::SolveIncAdjointEquationGN() {
    PetscErrorCode ierr = 0;
    ScalarType *pM = nullptr, *pLtilde = nullptr;
    TransportKernelAdjointSL kernel;
    IntType nt, nc, nl, lm;
    PetscFunctionBegin;
    
    this->m_Opt->Enter(__func__);
    
    ierr = Assert(this->m_StateVariable != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_VelocityField != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkScaField[2] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[0] != nullptr, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_WorkVecField[1] != nullptr, "null pointer"); CHKERRQ(ierr);
    
    ierr = Assert(this->m_Differentiation != nullptr, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    nc = this->m_Opt->m_Domain.nc;
    nl = this->m_Opt->m_Domain.nl;
    kernel.ht = this->m_Opt->GetTimeStepSize();
    kernel.scale = kernel.ht/static_cast<ScalarType>(nc);
    kernel.nl = nl;

    if (this->m_SemiLagrangianMethod == nullptr) {
        try {this->m_SemiLagrangianMethod = new SemiLagrangian(this->m_Opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
        ierr = this->m_SemiLagrangianMethod->SetWorkVecField(this->m_WorkVecField[0]); CHKERRQ(ierr);
        ierr = this->m_SemiLagrangianMethod->ComputeTrajectory(this->m_VelocityField, "adjoint"); CHKERRQ(ierr);
    }

    // compute divergence of velocity field
    ierr = GetRawPointer(this->m_WorkScaField[0], &kernel.pDivV); CHKERRQ(ierr);
    ierr = this->m_Differentiation->Divergence(kernel.pDivV,this->m_VelocityField); CHKERRQ(ierr);

    ierr = GetRawPointer(this->m_WorkScaField[1], &kernel.pDivVx); CHKERRQ(ierr);
    ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pDivVx, kernel.pDivV, "adjoint"); CHKERRQ(ierr);

    ierr = GetRawPointer(this->m_StateVariable, &pM); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_IncAdjointVariable, &pLtilde); CHKERRQ(ierr);
    ierr = GetRawPointer(this->m_WorkScaField[2], &kernel.pLx); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[0]->GetArrays(kernel.pGm); CHKERRQ(ierr);

    // initialize work vec field
    ierr = this->m_WorkVecField[1]->SetValue(0.0); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->GetArrays(kernel.pB); CHKERRQ(ierr);

    for (IntType j = 0; j < nt; ++j) {
        lm = (nt-j)*nc*nl;
        if (j == 0) kernel.scale *= 0.5;
        for (IntType k = 0; k < nc; ++k) {
            kernel.pL = pLtilde + k*nl;
            kernel.pLnext = kernel.pL;
            
            ierr = this->m_SemiLagrangianMethod->Interpolate(kernel.pLx, kernel.pL, "adjoint"); CHKERRQ(ierr);

            // compute gradient of m^j
            ierr = this->m_Differentiation->Gradient(kernel.pGm, pM+lm+k*nl); CHKERRQ(ierr);
            
            // compute body force
            ierr = kernel.ComputeBodyForcePart1();
        }  // for all image components
        if (j == 0) kernel.scale *= 2.0;
    }  // for all time points
    
    kernel.scale *= 0.5;

    // compute body force for last time point t = 0 (i.e., for j = nt)
    for (IntType k = 0; k < nc; ++k) {  // for all image components
        lm = k*nl;
        kernel.pL = pLtilde + k*nl;
        
        // compute gradient of m (for incremental body force)
        ierr = this->m_Differentiation->Gradient(kernel.pGm,pM+lm); CHKERRQ(ierr);
        
        ierr = kernel.ComputeBodyForcePart2(); CHKERRQ(ierr);
    }

    ierr = RestoreRawPointer(this->m_IncAdjointVariable, &pLtilde); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_StateVariable, &pM); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField[2], &kernel.pLx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField[1], &kernel.pDivVx); CHKERRQ(ierr);
    ierr = RestoreRawPointer(this->m_WorkScaField[0], &kernel.pDivV); CHKERRQ(ierr);

    ierr = this->m_WorkVecField[0]->RestoreArrays(kernel.pGm); CHKERRQ(ierr);
    ierr = this->m_WorkVecField[1]->RestoreArrays(kernel.pB); CHKERRQ(ierr);
    
    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




}  // namespace reg


#endif
