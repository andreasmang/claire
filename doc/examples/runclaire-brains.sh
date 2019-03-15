#!/bin/bash

# In this script we showcase how to run CLAIRE in parallel
# with the default settings. We execute CLAIRE with 20 MPI
# tasks (the number of MPI tasks be less or equql to the
# number of cores on your system). It is recommended to
# execute CLAIRE in parallel. If MPI is not available on
# your system remove the "mpirun -np 20" command.

# the datasets are two brain images of size 128x150x128

# define directory for binary (relative path)
bindir=../../bin

# define directory for data (relative path)
datdir=../data

# execute claire with default settings in parallel
mpirun -np 20 $bindir/claire -mr $datdir/brain01.nii.gz -mt $datdir/brain02.nii.gz




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################
