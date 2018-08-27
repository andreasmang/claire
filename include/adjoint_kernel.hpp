#include "petsc.h"
#include "petsccuda.h"

#ifndef __ADJOINT_KERNEL_HPP__
#define __ADJOINT_KERNEL_HPP__

void ComputeAdjointBodyForceGPU(PetscReal *p_lnext, PetscReal *p_l, PetscReal *p_lx, PetscReal *p_divv, PetscReal *p_divvx, PetscReal *p_vec1, PetscReal *p_vec2, PetscReal *p_vec3, PetscReal *p_b1, PetscReal *p_b2, PetscReal *p_b3, PetscReal ht, PetscReal scale, PetscInt nl);
void ComputeAdjointBodyForceGPU(PetscReal *p_l, PetscReal *p_vec1, PetscReal *p_vec2, PetscReal *p_vec3, PetscReal *p_b1, PetscReal *p_b2, PetscReal *p_b3, PetscReal scale, PetscInt nl);

#endif
