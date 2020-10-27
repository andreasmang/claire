#!/bin/bash

gpu="V100" # TACC Longhorn
#gpu="RTX" # TACC Frontera
#gpu="P100" # TACC Maverick2

./build_libs.sh --gpu=$gpu

./build_libs.sh --bpetsccudasgl --enableCUDA --gpu=$gpu

./build_libs.sh --bfftw --enableCUDA --POWER9

./build_libs.sh --baccfft --enableCUDA --gpu=$gpu

./build_libs.sh --bzlib --enableCUDA

./build_libs.sh --bnifti --enableCUDA

./build_libs.sh --bmorton --enableCUDA

