#!/bin/bash

compute_system=longhorn

gpu="V100" # TACC Longhorn
#gpu="RTX" # TACC Frontera
#gpu="P100" # TACC Maverick2

./build_libs.sh --gpu=$gpu

if [ $compute_system == "longhorn" ]; then
    ./build_libs.sh --bpetsccudasgl --enableCUDA --gpu=$gpu
fi

if [ $compute_system == "longhorn" ]; then
    ./build_libs.sh --bfftw --enableCUDA --POWER9
else
    ./build_libs.sh --bfftw --enableCUDA
fi

./build_libs.sh --baccfft --enableCUDA --gpu=$gpu

./build_libs.sh --bzlib --enableCUDA

./build_libs.sh --bnifti --enableCUDA

