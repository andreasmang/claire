#!/bin/bash
set -x

# In this script we showcase how to run CLAIRE in parallel.
# We use a neuroimaging dataset. We perform a parameter 
# continuation until we reach a target regularization 
# parameter. We store the velocity field as output.
# The datasets are two brain images of size 128x150x128.

# define directory for binary (relative path)
bindir=../../bin

# define directory for data (relative path)
datdir=../data

# execute claire with default settings in parallel
mpirun -np 20 $bindir/claire -mr $datdir/brain01.nii.gz -mt $datdir/brain02.nii.gz -betacont 7.750000e-04 -x ./ -velocity


#### commands explained
# > mpirun -np 20   --- execute with 20 MPI tasks 
# > $bindir/claire  --- name of the binary 
# > -mr <file>      --- reference image (fixed image)
# > -mt <file>      --- template image (image to be registered) 
# > -betacont <dbl> --- use continuation in regularization
#                       parameter until target regularization
#                       parameter is reached 
# > -x <folder>     --- output folder
# > -velocity       --- write computed velocities to file: 
#                       velocity-field-x1.nii.gz (x1 component) 
#                       velocity-field-x2.nii.gz (x2 component) 
#                       velocity-field-x3.nii.gz (x3 component)




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################
