#!/bin/bash

compute_system=$1
gpu=$2

./build_libs.sh --gpu=$gpu

#if [ $compute_system -neq "frontera" ]; then
    ./build_libs.sh --bpetsccudasgl --enableCUDA --gpu=$gpu
#fi

if [ $compute_system == "longhorn" ]; then
    ./build_libs.sh --bfftw --enableCUDA --POWER9
else
    ./build_libs.sh --bfftw --enableCUDA
fi

./build_libs.sh --baccfft --enableCUDA --gpu=$gpu

./build_libs.sh --bzlib --enableCUDA

./build_libs.sh --bnifti --enableCUDA

./build_libs.sh --bmorton --enableCUDA

