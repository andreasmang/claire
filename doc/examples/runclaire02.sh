#!/bin/bash

# In this script we execute CLAIRE for a synthetic test problem.
# We also show how one can execute CLAIRE in parallel for a
# problem size of nx=128x128x128 with 20 MPI tasks. We recommend
# to execute CLAIRE on multi-core systems to significantly
# reduce the runtime.




# define directory for binary (relative path)
bindir=../../bin

# execute claire with default settings in parallel
# (using 20 MPI tasks (command: "mpirun -np 20"))
mpirun -np 20 $bindir/claire -synthetic 0 -nx 128

#### commands explained
# > mpirun -np 20    --- execute with 20 MPI tasks 
# > $bindir/claire   --- name of the binary 
# > -synthetic <int> --- select synthetic test problem 
# > -nx size         --- define size of problem
#                        (other possibility: -nx 128x128x128) 




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################
