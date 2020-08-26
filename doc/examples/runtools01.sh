#!/bin/bash
set -x

# execute setup (in case velocities are not available)
./check4vel.sh

# We provide an example on how to transport an image.
# Our input are an image and the three components of
# the velocity field (x1, x2, and x3). The output is
# the transported/deformed intput image.
# The dataset is a brain image of size 128x150x128.

# define directory for binary (relative path)
bindir=../../bin

# define directory for data (relative path)
datdir=../data

# execute clairetools (backslashes are added for line breaks)
mpirun -np 20 $bindir/clairetools -v1 velocity-field-x1.nii.gz       \
                                  -v2 velocity-field-x2.nii.gz       \
                                  -v3 velocity-field-x3.nii.gz       \
                                  -ifile $datdir/brain01.nii.gz      \
                                  -xfile brain01-transported.nii.gz -deformimage


#### commands explained
# > mpirun -np 20        --- execute with 20 MPI tasks 
# > $bindir/clairetools  --- name of the binary 
# > -v1 <file>           --- x1 component of velocity field (input)
# > -v2 <file>           --- x2 component of velocity field (input)
# > -v3 <file>           --- x3 component of velocity field (input)
# > -ifile <file>        --- template image (image to be transported; input)
# > -xfile <file>        --- output image (output)
# > -deformimage         --- flag to execute deformation / transport input image 




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################

