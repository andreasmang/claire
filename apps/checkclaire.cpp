#include <iostream>
#include "petsc.h"

int main(int argc, char **argv) {
    PetscErrorCode ierr;
    PetscScalar maxval, minval;
    PetscInt posmin, posmax;
    Vec x; PetscScalar* p_x = NULL;

    ierr = PetscInitialize(0, reinterpret_cast<char***>(NULL),
                              reinterpret_cast<char*>(NULL),
                              reinterpret_cast<char*>(NULL)); CHKERRQ(ierr);
    PetscInt n = 16;

    ierr = VecCreate(PETSC_COMM_WORLD, &x); CHKERRQ(ierr);
    ierr = VecSetSizes(x, n, n); CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); CHKERRQ(ierr);

    ierr = VecGetArray(x, &p_x); CHKERRQ(ierr);
    for (PetscInt i = 0; i < n; ++i) {
        p_x[i] = static_cast<PetscScalar>(i);
        std::cout << p_x[i] << " ";
    }
    ierr = VecRestoreArray(x, &p_x); CHKERRQ(ierr);

    std::cout << std::endl;

    ierr = VecMin(x, &posmin, &minval); CHKERRQ(ierr);
    ierr = VecMax(x, &posmax, &maxval); CHKERRQ(ierr);

    std::cout << "min " << minval << " at " << posmin << std::endl;
    std::cout << "max " << maxval << " at " << posmax << std::endl;

    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}
