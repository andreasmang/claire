#!/bin/bash

# In this script we execute CLAIRE for a synthetic test problem.

# define directory for binary (relative path)
bindir=../../bin

# execute claire with default settings
$bindir/claire -synthetic 0


#### commands explained
# > $bindir/claire --- name of the binary 
# > -synthetic 0   --- select synthetic test problem 




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################
