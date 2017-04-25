#!/bin/bash

if [[ ! -e accfft.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading accFFT"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/accfft.tar.gz .
fi

if [[ ! -e fftw.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading FFTW"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/fftw.tar.gz .
fi


if [[ ! -e morton.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading morton library"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/morton.tar.gz
fi

if [[ ! -e nifticlib.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading nifti library"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/nifticlib.tar.gz
fi

if [[ ! -e parallel-netcdf.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading netcdf library"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/parallel-netcdf.tar.gz
fi

if [[ ! -e petsc-lite.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading petsc library"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/petsc.tar.gz
fi

if [[ ! -e zlib.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading zlib"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/zlib.tar.gz
fi
