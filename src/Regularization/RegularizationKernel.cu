#include "RegularizationKernel.hpp"
#include "cuda_helper.hpp"

#include "RegularizationKernel.txx"

namespace reg {
  
using KernelUtils::ReductionKernelCallGPU;

PetscErrorCode RegularizationKernel::LocalNorm (ScalarType &lnorm) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ReductionKernelCallGPU<NormKernel>(lnorm, nl, pX[0], pX[1], pX[2]); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

}
