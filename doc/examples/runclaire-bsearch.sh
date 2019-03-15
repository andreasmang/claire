#!/bin/bash

# In this script we showcase how to run CLAIRE in parallel.
# We use a neuroimaging dataset. We try to identify an
# adequate regularization parameter using a binary search  


# the datasets are two brain images of size 128x150x128

# define directory for binary (relative path)
bindir=../../bin

# define directory for data (relative path)
datdir=../data

# execute claire with default settings in parallel
mpirun -np 20 $bindir/claire -mr $datdir/brain01.nii.gz -mt $datdir/brain02.nii.gz -train binary




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################
