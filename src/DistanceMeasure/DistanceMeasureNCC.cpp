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
 *  along with CLAIRE.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _DISTANCEMEASURENCC_CPP_
#define _DISTANCEMEASURENCC_CPP_

#include "DistanceMeasureNCC.hpp"
#include "DistanceMeasureKernel.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
DistanceMeasureNCC::DistanceMeasureNCC() : SuperClass() {
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
DistanceMeasureNCC::~DistanceMeasureNCC() {
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
DistanceMeasureNCC::DistanceMeasureNCC(RegOpt* opt) : SuperClass(opt) {
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
PetscErrorCode DistanceMeasureNCC::ClearMemory() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    PetscFunctionReturn(ierr);
}

/********************************************************************
* @brief Set up scale for NCC measure given by
* scale = D_SL2(mR,mT)/D_NCC(mR,mT)
********************************************************************/
PetscErrorCode DistanceMeasureNCC::SetupScale(){
    PetscErrorCode ierr = 0;	
    DistanceMeasureKernel::EvaluateFunctionalNCC kernel;
    //ScalarType *p_mr = NULL, *p_mt = NULL, *p_w = NULL;
    IntType nc, ng;
    ScalarType norm_l2, norm_mT, norm_mR, inpr_mT_mR, sum_mT, sum_mR;
    int rval;
    ScalarType l2distance, nccdistance, hd;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);
    

    ierr = Assert(this->m_TemplateImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // get sizes
    kernel.nc = this->m_Opt->m_Domain.nc;
    nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    ng = this->m_Opt->m_Domain.ng;
    hd  = this->m_Opt->GetLebesgueMeasure();   

    ierr = this->m_TemplateImage->GetArrayRead(kernel.pMt); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->GetArrayRead(kernel.pMr); CHKERRQ(ierr);
    ierr = GetRawPointerRead(this->m_ObjWts, &kernel.pWts); CHKERRQ(ierr);

    if (this->m_Mask != NULL) {
        // Mask objective functional
        ierr = this->m_Mask->GetArrayRead(kernel.pW); CHKERRQ(ierr);

        ierr = kernel.ComputeScaleMask(); CHKERRQ(ierr);
    
        ierr = this->m_Mask->RestoreArray(); CHKERRQ(ierr);
    } else {
        ierr = kernel.ComputeScale(); CHKERRQ(ierr);
    }
    // All reduce the pieces
    rval = MPI_Allreduce(&kernel.norm_l2_loc, &norm_l2, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&kernel.norm_mT_loc, &norm_mT, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&kernel.norm_mR_loc, &norm_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&kernel.inpr_mT_mR_loc, &inpr_mT_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&kernel.sum_mR_loc, &sum_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.sum_mT_loc, &sum_mT, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    
    ierr = this->m_ReferenceImage->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_TemplateImage->RestoreArray(); CHKERRQ(ierr);
    ierr = RestoreRawPointerRead(this->m_ObjWts, &kernel.pWts); CHKERRQ(ierr);

    // update the inner products and norms to include the image mean differences
    inpr_mT_mR = inpr_mT_mR - sum_mR*sum_mT/ng;
    norm_mT = norm_mT - sum_mT*sum_mT/ng;
    norm_mR = norm_mR - sum_mR*sum_mR/ng;

    // Calculate two distance measures and scale
    l2distance = 0.5*hd*norm_l2/static_cast<ScalarType>(nc);
    nccdistance = 0.5 - 0.5*(inpr_mT_mR*inpr_mT_mR)/(norm_mT*norm_mR);
    this->m_Opt->m_Distance.scale = l2distance/nccdistance;
    //std::cout<<"scale = " << l2distance/nccdistance << std::endl;
    this->m_Opt->Exit(__func__);


    PetscFunctionReturn(ierr);   
}


/********************************************************************
 * @brief evaluate the functional (i.e., the distance measure)
 * D = 0.5 - 0.5* <m1,mR>_L2/(||m1||_L2 * ||mR||_L2)
 *******************************************************************/
PetscErrorCode DistanceMeasureNCC::EvaluateFunctional(ScalarType* D) {
    PetscErrorCode ierr = 0;
    //ScalarType *p_mr = NULL, *p_m = NULL, *p_w = NULL;
    DistanceMeasureKernel::EvaluateFunctionalNCC kernel;
    IntType nt, ng;
    ScalarType norm_m1, norm_mR, inpr_m1_mR, sum_m1, sum_mR, scale;
    int rval;

    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    // Get sizes
    nt = this->m_Opt->m_Domain.nt;
    ng = this->m_Opt->m_Domain.ng;
    kernel.nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    scale = this->m_Opt->m_Distance.scale;
     
    ierr = this->m_StateVariable->GetArrayRead(kernel.pM, 0, nt); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->GetArrayRead(kernel.pMr); CHKERRQ(ierr);

    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = this->m_Mask->GetArrayRead(kernel.pW); CHKERRQ(ierr);

        ierr = kernel.ComputeFunctionalMask(); CHKERRQ(ierr);
    
        ierr = this->m_Mask->RestoreArray(); CHKERRQ(ierr);
    } else {
        ierr = kernel.ComputeFunctional(); CHKERRQ(ierr);
    }
    // All reduce various pieces
    rval = MPI_Allreduce(&kernel.norm_m1_loc, &norm_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&kernel.norm_mR_loc, &norm_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&kernel.inpr_m1_mR_loc, &inpr_m1_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.sum_m1_loc, &sum_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.sum_mR_loc, &sum_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    ierr = this->m_ReferenceImage->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    
    // save the inner products in order to reuse them for computing final condition for AE and IAE
    
    //PetscPrintf(PETSC_COMM_WORLD, "sum_m1 = %0.4e\n", sum_m1);
    //PetscPrintf(PETSC_COMM_WORLD, "sum_mR = %0.4e\n", sum_mR);
    //PetscPrintf(PETSC_COMM_WORLD, "inpr_m1_mR = %0.4e\n", inpr_m1_mR);
    //PetscPrintf(PETSC_COMM_WORLD, "norm_m1 = %0.4e\n", norm_m1);
    //PetscPrintf(PETSC_COMM_WORLD, "norm_mR = %0.4e\n", norm_mR);
    //PetscPrintf(PETSC_COMM_WORLD, "ng diff = %0.4e\n", 128*128*128-ng);
    //PetscPrintf(PETSC_COMM_WORLD, "D_o = %0.4e\n", 0.5 - 0.5*inpr_m1_mR*inpr_m1_mR/(norm_m1*norm_mR));

    // subtract means from inner products and norms
    inpr_m1_mR = inpr_m1_mR - sum_m1*sum_mR/ng;
    norm_m1 = norm_m1 - sum_m1*sum_m1/ng;
    norm_mR = norm_mR - sum_mR*sum_mR/ng;

    //PetscPrintf(PETSC_COMM_WORLD, "inpr_m1_mR_s = %0.4e\n", inpr_m1_mR);
    //PetscPrintf(PETSC_COMM_WORLD, "norm_m1_s = %0.4e\n", norm_m1);
    //PetscPrintf(PETSC_COMM_WORLD, "norm_mR_s = %0.4e\n", norm_mR);
    //PetscPrintf(PETSC_COMM_WORLD, "scale_S = %0.4e\n", scale);
    //PetscPrintf(PETSC_COMM_WORLD, "true = %0.4e\n", 0.5 - 0.5*(inpr_m1_mR*inpr_m1_mR)/(norm_m1*norm_mR));

    // Objective value
    *D = scale*(0.5 - 0.5*(inpr_m1_mR*inpr_m1_mR)/(norm_m1*norm_mR));
    //PetscPrintf(PETSC_COMM_WORLD, "D = %0.4e\n", *D);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for adjoint equation
 * (varies for different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureNCC::SetFinalConditionAE() {
    PetscErrorCode ierr = 0;
    IntType nt, ng;
    int rval;
    //ScalarType *p_mr = NULL, *p_m = NULL, *p_l = NULL, *p_w = NULL;
    ScalarType norm_m1, norm_mR, inpr_m1_mR, sum_m1, sum_mR;
    ScalarType hd, scale;
    DistanceMeasureKernel::FinalConditionNCC kernel;
    PetscFunctionBegin;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_AdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    ng = this->m_Opt->m_Domain.ng;
    kernel.nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    hd  = this->m_Opt->GetLebesgueMeasure();
    scale = this->m_Opt->m_Distance.scale;

    // Index for final condition
    ierr = this->m_StateVariable->GetArrayRead(kernel.pM, 0, nt); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->GetArrayRead(kernel.pMr); CHKERRQ(ierr);

    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
      ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL, 0, nt); CHKERRQ(ierr);
    } else {
      ierr = this->m_AdjointVariable->GetArrayReadWrite(kernel.pL); CHKERRQ(ierr);
    }

    /* compute terminal condition
	lambda = (<m1,mR>/<m1,m1><mR,mR>) * (mR - (<m1,mR>/<m1,m1>)*m1)
    */
    ierr = kernel.ComputeInnerProductsFinalConditionAE(); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.norm_mR_loc, &norm_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.norm_m1_loc, &norm_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.inpr_m1_mR_loc, &inpr_m1_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.sum_m1_loc, &sum_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.sum_mR_loc, &sum_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    // first load the inner products
    //l = nt*nc*nl;
    //norm_m1 = this->inpr_m1;
    //norm_mR = this->inpr_mR;
    //inpr_m1_mR = this->inpr_m1_mR;
    //
    //subtract means from inner products
    inpr_m1_mR = inpr_m1_mR - sum_m1*sum_mR/ng;
    norm_m1 = norm_m1 - sum_m1*sum_m1/ng;
    norm_mR = norm_mR - sum_mR*sum_mR/ng;

    // Now, write the terminal condition to lambda
    kernel.const1 = scale*inpr_m1_mR/(hd*norm_m1*norm_mR);
    kernel.const2 = scale*(inpr_m1_mR*inpr_m1_mR)/(hd*norm_m1*norm_m1*norm_mR);

    //std::cout << "gradient const1 = " << kernel.const1 << std::endl;
    //std::cout << "gradient const2 = " << kernel.const2 << std::endl;
    //

    kernel.mean_m1 = sum_m1/ng;
    kernel.mean_mR = sum_mR/ng;
    

   if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = this->m_Mask->GetArrayRead(kernel.pW); CHKERRQ(ierr);
        
        // not implemented
        //ierr = kernel.ComputeFinalConditionMaskAE(); CHKERRQ(ierr);
    
        ierr = this->m_Mask->RestoreArray(); CHKERRQ(ierr);
    } else {
        ierr = kernel.ComputeFinalConditionAE(); CHKERRQ(ierr);
    }  

    ierr = this->m_AdjointVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set final condition for incremental adjoint equation
 * (varies for different distance measures)
 *******************************************************************/
PetscErrorCode DistanceMeasureNCC::SetFinalConditionIAE() {
    PetscErrorCode ierr = 0;
    IntType nt, ng;
    int rval;
    //ScalarType const1, const2, const3, const4, const5;
    ScalarType hd, scale; 
    ScalarType norm_m1, norm_mR, inpr_m1_mR, inpr_m1_mtilde, inpr_mR_mtilde, sum_m1, sum_mR, sum_mtilde;
    DistanceMeasureKernel::FinalConditionNCC kernel;

    this->m_Opt->Enter(__func__);

    ierr = Assert(this->m_IncAdjointVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncStateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_StateVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_ReferenceImage != NULL, "null pointer"); CHKERRQ(ierr);

    nt = this->m_Opt->m_Domain.nt;
    ng = this->m_Opt->m_Domain.ng;
    kernel.nc = this->m_Opt->m_Domain.nc;
    kernel.nl = this->m_Opt->m_Domain.nl;
    hd  = this->m_Opt->GetLebesgueMeasure();
    scale = this->m_Opt->m_Distance.scale;

    // Index for final condition
    ierr = this->m_StateVariable->GetArrayRead(kernel.pM, 0, nt); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->GetArrayRead(kernel.pMr); CHKERRQ(ierr);
    
    if (this->m_Opt->m_OptPara.method == FULLNEWTON) {
        ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pLtilde, 0, nt); CHKERRQ(ierr);
        ierr = this->m_IncStateVariable->GetArrayRead(kernel.pMtilde, 0, nt); CHKERRQ(ierr);
    } else {
      ierr = this->m_IncAdjointVariable->GetArrayReadWrite(kernel.pLtilde); CHKERRQ(ierr);
      ierr = this->m_IncStateVariable->GetArrayRead(kernel.pMtilde); CHKERRQ(ierr);
    }
    
    ierr = kernel.ComputeInnerProductsFinalConditionIAE(); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.norm_mR_loc, &norm_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.norm_m1_loc, &norm_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.inpr_m1_mR_loc, &inpr_m1_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.inpr_m1_mtilde_loc, &inpr_m1_mtilde, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);

    rval = MPI_Allreduce(&kernel.inpr_mR_mtilde_loc, &inpr_mR_mtilde, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.sum_m1_loc, &sum_m1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.sum_mR_loc, &sum_mR, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    rval = MPI_Allreduce(&kernel.sum_mtilde_loc, &sum_mtilde, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    ierr = Assert(rval == MPI_SUCCESS, "mpi error"); CHKERRQ(ierr);
    
    // Get the rest of the inner products
    //norm_m1 = this->inpr_m1;
    //norm_mR = this->inpr_mR;
    //inpr_m1_mR = this->inpr_m1_mR;

    norm_m1 = norm_m1 - sum_m1*sum_m1/ng;
    norm_mR = norm_mR - sum_mR*sum_mR/ng;
    inpr_m1_mR = inpr_m1_mR - sum_m1*sum_mR/ng;
    inpr_m1_mtilde = inpr_m1_mtilde - sum_m1*sum_mtilde/ng;
    inpr_mR_mtilde = inpr_mR_mtilde - sum_mR*sum_mtilde/ng;
    kernel.mean_m1 = sum_m1/ng;
    kernel.mean_mR = sum_mR/ng;

    
    // Now, write the terminal condition to lambda tilde
    kernel.const1 = scale*inpr_mR_mtilde/(hd*norm_m1*norm_mR);
    Assert(!PetscIsNanReal(kernel.const1), "is nan");
    kernel.const2 = 2.0*scale*(inpr_m1_mR*inpr_m1_mtilde)/(hd*norm_m1*norm_m1*norm_mR);
    Assert(!PetscIsNanReal(kernel.const2), "is nan");
    kernel.const3 = 4.0*scale*(inpr_m1_mR*inpr_m1_mR*inpr_m1_mtilde)/(hd*norm_m1*norm_m1*norm_m1*norm_mR);
    Assert(!PetscIsNanReal(kernel.const3), "is nan");
    kernel.const4 = 2.0*scale*(inpr_m1_mR*inpr_mR_mtilde)/(hd*norm_m1*norm_m1*norm_mR);
    Assert(!PetscIsNanReal(kernel.const4), "is nan");
    kernel.const5 = scale*(inpr_m1_mR*inpr_m1_mR)/(hd*norm_m1*norm_m1*norm_mR);
    Assert(!PetscIsNanReal(kernel.const5), "is nan");
    
    //std::cout << "norm2 mtilde = " << norm_mtilde << std::endl;
    //std::cout << "hessian const1 = " << kernel.const1 << std::endl;
    //std::cout << "hessian const2 = " << kernel.const2 << std::endl;
    //std::cout << "hessian const3 = " << kernel.const3 << std::endl;
    //std::cout << "hessian const4 = " << kernel.const4 << std::endl;
    //std::cout << "hessian const5 = " << kernel.const5 << std::endl;

    if (this->m_Mask != NULL) {
        // mask objective functional
        ierr = this->m_Mask->GetArrayRead(kernel.pW); CHKERRQ(ierr);
        
        // not implemented
        //ierr = kernel.ComputeFinalConditionMaskIAE(); CHKERRQ(ierr);
    
        ierr = this->m_Mask->RestoreArray(); CHKERRQ(ierr);
    } else {
        ierr = kernel.ComputeFinalConditionIAE(); CHKERRQ(ierr);
    }  

    ierr = this->m_IncAdjointVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_IncStateVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_StateVariable->RestoreArray(); CHKERRQ(ierr);
    ierr = this->m_ReferenceImage->RestoreArray(); CHKERRQ(ierr);

    this->m_Opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

}  // namespace reg




#endif  // _DISTANCEMEASURENCC_CPP_
