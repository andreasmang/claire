#!/bin/bash

# In this script we execute CLAIRE for a synthetic test problem.
# We also show how one can execute CLAIRE in parallel for a
# larger problem size (nx=128x128x128) 

# define directory for binary (relative path)
bindir=../../bin

# execute claire with default settings
$bindir/claire -synthetic 0

# uncomment to execute in parallel 
# mpirun -np 20 $bindir/claire -synthetic 0 -nx 128




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################
