#!/bin/bash
set -x

# execute setup (in case velocities are not available)
./check4vel.sh

# We provide an example on how to compute the determinant
# of the deformation gradient det \grad y given a
# velocity field 

# define directory for binary (relative path)
bindir=../../bin

# define directory for data (relative path)
datdir=../data

# execute clairetools (backslashes are added for line breaks)
mpirun -np 20 $bindir/clairetools -v1 velocity-field-x1.nii.gz       \
                                  -v2 velocity-field-x2.nii.gz       \
                                  -v3 velocity-field-x3.nii.gz       \
                                  -x ./ -detdefgrad


#### commands explained
# > mpirun -np 20        --- execute with 20 MPI tasks 
# > $bindir/clairetools  --- name of the binary 
# > -v1 <file>           --- x1 component of velocity field (input)
# > -v2 <file>           --- x2 component of velocity field (input)
# > -v3 <file>           --- x3 component of velocity field (input)
# > -x <folder>          --- output path (where to store det-deformation-grad.nii.gz) 
# > -detdefgrad          --- flag to compute determinant of deformation gradient 




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################

