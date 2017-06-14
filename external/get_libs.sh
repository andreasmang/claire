#!/bin/bash

if [[ ! -e accfft.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading accFFT"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/accfft.tar.gz
#	wget -O accfft-master.zip https://github.com/amirgholami/accfft/archive/master.zip
#	tar -czvf master.tar.gz master
fi

if [[ ! -e fftw.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading FFTW"
	echo "----------------------------------------------------------------------------------"
	wget -O fftw.tar.gz http://www.fftw.org/fftw-3.3.6-pl2.tar.gz
#	wget http://users.ices.utexas.edu/~andreas/libs/fftw.tar.gz .
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
	#wget -O nifticlib.tar.gz --no-check-certificate http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz
	wget http://users.ices.utexas.edu/~andreas/libs/nifticlib.tar.gz
fi

if [[ ! -e parallel-netcdf.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading netcdf library"
	echo "----------------------------------------------------------------------------------"
	wget -O parallel-netcdf.tar.gz http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.8.1.tar.gz 
	#wget http://users.ices.utexas.edu/~andreas/libs/parallel-netcdf.tar.gz
fi

if [[ ! -e petsc-lite.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading petsc library"
	echo "----------------------------------------------------------------------------------"
	wget -O petsc-lite.tar.gz http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz
#	wget http://users.ices.utexas.edu/~andreas/libs/petsc.tar.gz
fi

if [[ ! -e zlib.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading zlib"
	echo "----------------------------------------------------------------------------------"
	wget http://users.ices.utexas.edu/~andreas/libs/zlib.tar.gz
fi
