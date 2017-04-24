#!/bin/bash

if [[ ! -e accfft.tar.gz ]]; then
	echo "----------------------------------------------------------------------------------"
	echo "downloading accFFT"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/accfft.tar.gz .
fi

if [[ ! -e fftw-3.3.6-pl1.tar.gz ]]; then
	echo "----------------------------------------------------------------------------------"
	echo "downloading FFTW 3.3.6 PL1"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/fftw-3.3.6-pl1.tar.gz .
fi


if [[ ! -e libmorton.tar.gz ]]; then
	echo "----------------------------------------------------------------------------------"
	echo "downloading morton library"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/libmorton.tar.gz
fi

if [[ ! -e nifticlib-2.0.0.tar.gz ]]; then
	echo "----------------------------------------------------------------------------------"
	echo "downloading nifti library 2.0.0"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/nifticlib-2.0.0.tar.gz
fi

if [[ ! -e parallel-netcdf-1.7.0.tar.gz ]]; then
	echo "----------------------------------------------------------------------------------"
	echo "downloading netcdf library 1.7.0"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/parallel-netcdf-1.7.0.tar.gz
fi

if [[ ! -e petsc-lite-3.7.3.tar.gz ]]; then
	echo "----------------------------------------------------------------------------------"
	echo "downloading petsc library 3.7.3"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/petsc-lite-3.7.3.tar.gz
fi

if [[ ! -e zlib-1.2.8.tar.gz ]]; then
	echo "----------------------------------------------------------------------------------"
	echo "downloading zlib 1.2.8"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/zlib-1.2.8.tar.gz
fi
