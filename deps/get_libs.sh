#!/bin/bash

if [[ ! -e accfft.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading accFFT"
	echo "----------------------------------------------------------------------------------"
	#wget --no-check-certificate --content-disposition https://github.com/amirgholami/accfft/archive/master.tar.gz -O accfft.tar.gz
	wget --no-check-certificate https://github.com/amirgholami/accfft/archive/master.tar.gz -O accfft.tar.gz
fi

if [[ ! -e fftw.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading FFTW"
	echo "----------------------------------------------------------------------------------"
	wget http://www.fftw.org/fftw-3.3.6-pl2.tar.gz -O fftw.tar.gz 
fi

if [[ ! -e nifticlib-2.0.0.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading nifti library"
	echo "----------------------------------------------------------------------------------"
#	wget http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz -O nifticlib.tar.gz
	wget http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz
fi

if [[ ! -e parallel-netcdf.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading netcdf library"
	echo "----------------------------------------------------------------------------------"
	wget http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.8.1.tar.gz -O parallel-netcdf.tar.gz
fi

#if [[ ! -e petsc.tar.gz ]]; then
if [[ ! -e petsc-lite-3.11.4.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading petsc library"
	echo "----------------------------------------------------------------------------------"
#	wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz -O petsc.tar.gz
#	wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz
#   wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.4.tar.gz # does not work for OpenMPI 4.0.3
	wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.4.tar.gz
fi

if [[ ! -e zlib.tar.gz ]]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "downloading zlib"
	echo "----------------------------------------------------------------------------------"
	wget https://zlib.net/zlib-1.2.11.tar.gz -O zlib.tar.gz
fi
