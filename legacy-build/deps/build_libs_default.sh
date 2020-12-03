#!/bin/bash


# extract all libraries
./build_libs

./build_libs --bfftw

./build_libs --baccfft

./build_libs --petscdbl

./build_libs --petscsgl
