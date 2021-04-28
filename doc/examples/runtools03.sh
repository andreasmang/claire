#!/bin/bash

# In this script we show hot to convert a dataset from *.nii.gz to
# PNETCDF data (for parallel IO; important in large-scale applications).
# Notice that the output file will be stored with the same file name
# and in the same folder where the input file is located (i.e., for
# example $datdir/brain.nii.gz is converted to $datdir/brain.nc).




# define directory for binary (relative path)
bindir=../../bin

# define directory for data (relative path)
datdir=../data

# execute claire with default settings
$bindir/clairetools -ifile $datdir/brain01.nii.gz -convert 2nc

#### commands explained
# > $bindir/clairetools  --- name of the binary 
# > -ifile <file>        --- name of input file (nii.gz data in our case) 
# > -convert 2nc         --- operation to be performed (conversion to
#                            a NETCDF format




################################################################
# This script is part of the C++ software
#                        ---  CLAIRE  ---
# For details see https://github.com/andreasmang/claire
################################################################
