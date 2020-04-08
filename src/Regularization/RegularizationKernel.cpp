#include "RegularizationKernel.hpp"

#include "RegularizationKernel.txx"

namespace reg {
  
using KernelUtils::ReductionKernelCall;

PetscErrorCode RegularizationKernel::LocalNorm (ScalarType &lnorm) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  ierr = ReductionKernelCall<NormKernel>(lnorm, nl, pX[0], pX[1], pX[2]); CHKERRQ(ierr);
    
  PetscFunctionReturn(ierr);
}

}
