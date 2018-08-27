#include "petsc.h"
#include "petsccuda.h"

#ifndef __DISTANCE_KERNEL_HPP__
#define __DISTANCE_KERNEL_HPP__

void DistanceMeasureSetFinalGPU(PetscReal *p_l, const PetscReal *p_m, const PetscReal *p_mr, PetscInt nl);
void DistanceMeasureSetFinalMaskGPU(PetscReal *p_l, const PetscReal *p_m, const PetscReal *p_mr, const PetscReal *p_w, PetscInt nl, PetscInt nc);

#endif
