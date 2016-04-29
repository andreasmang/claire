#!/bin/bash

CODE_DIR=/work/03279/amang/maverick/code/cold
BUILD_DIR=/work/03279/amang/maverick/code/build/libs/build

MPI_CXX=mpicxx
MPI_C=mpicc

builddep=1


LIB_DIR=${CODE_DIR}/external

if [ -e ${BUILD_DIR}/environment_vars.sh ]; then
	rm ${BUILD_DIR}/environment_vars.sh
fi

################################
# FFTW
################################
FFTW_LIB_DIR=${BUILD_DIR}/fftw

if [ ! -d ${FFTW_LIB_DIR} ]; then
	mkdir ${FFTW_LIB_DIR}
	mkdir ${FFTW_LIB_DIR}/fftw-srcs
	mkdir ${FFTW_LIB_DIR}/fftw-build
	echo extracting FFTW lib...
	tar -xzf ${LIB_DIR}/fftw-3.3.4.tar.gz -C ${SRC_DIR} --strip-components=1
fi

SRC_DIR=${FFTW_LIB_DIR}/fftw-srcs
BLD_DIR=${FFTW_LIB_DIR}/fftw-build

if [ ${builddep} -eq 1 ]; then 
	echo "configuring FFTW (double precision)" 
	cd ${SRC_DIR}
#	./configure --prefix=${BLD_DIR} --enable-mpi --enable-threads --enable-sse2 --enable-openmp --enable-avx CFLAGS=-O3 MAKEINFO=missing
#	make
#	make install

	echo "configuring FFTW (float precision)" 
	#./configure --prefix=${BLD_DIR} --enable-mpi --enable-threads --enable-sse2 --enable-openmp --enable-avx --enable-float CFLAGS=-O3 MAKEINFO=missing
	#make
	#make install

fi

echo "export FFTW_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/lib:${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh



################################
# PNETCDF
################################
PNETCDF_LIB_DIR=${BUILD_DIR}/pnetcdf
SRC_DIR=${PNETCDF_LIB_DIR}/pnetcdf-srcs
BLD_DIR=${PNETCDF_LIB_DIR}/pnetcdf-build
if [ ! -d ${PNETCDF_LIB_DIR} ]; then
	mkdir ${PNETCDF_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}
	echo extracting PNETCDF lib...
	tar -xzf ${LIB_DIR}/parallel-netcdf-1.7.0.tar.gz -C ${SRC_DIR} --strip-components=1
fi


if [ ${builddep} -eq 1 ]; then 
	echo "configuring PNETCDF" 
	cd ${SRC_DIR}
	#./configure --prefix=${BLD_DIR}
	echo "building PNETCDF" 
	#make -j
	#make install
fi


echo "export PNETCDF_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh




################################
# ACCFFT
################################
ACCFFT_LIB_DIR=${BUILD_DIR}/accfft
SRC_DIR=${ACCFFT_LIB_DIR}/accfft-srcs
BLD_DIR=${ACCFFT_LIB_DIR}/accfft-build
CMK_DIR=${ACCFFT_LIB_DIR}/accfft-tmp
if [ ! -d ${ACCFFT_LIB_DIR} ]; then
	mkdir ${ACCFFT_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}
	mkdir ${CMK_DIR}
	echo extracting ACCFFT lib...
	tar -xzf ${LIB_DIR}/accfft.tar.gz -C ${SRC_DIR} --strip-components=1
fi


if [ ${builddep} -eq 1 ]; then 
	echo "configuring ACCFFT" 
	cd ${CMK_DIR}
	#cmake ${SRC_DIR} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/fftw-build -DFFTW_USE_STATIC_LIBS=true -DBUILD_GPU=false -DBUILD_STEPS=true -DCXX_FLAGS="-O3" -DBUILD_SHARED=false -DPNETCDF_DIR=${PNETCDF_LIB_DIR}/pnetcdf-build
	echo "building ACCFFT" 
	#make -j
	#make install
fi

echo "export ACCFFT_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh


################################
# PETSC
################################
PETSC_LIB_DIR=${BUILD_DIR}/petsc
SRC_DIR=${PETSC_LIB_DIR}/petsc-srcs
BLD_DIR=${PETSC_LIB_DIR}/petsc-build
if [ ! -d ${PETSC_LIB_DIR} ]; then
	mkdir ${PETSC_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}
	echo extracting PETSC lib...
	tar -xzf ${LIB_DIR}/petsc-lite-3.7.0.tar.gz -C ${SRC_DIR} --strip-components=1
fi

PETSC_ARCH=cxx_opt
if [ ${builddep} -eq 1 ]; then 
	echo "configuring PETSC" 
	cd ${SRC_DIR}
#	./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR} --with-clanguage=C++ --with-debugging=0 --download-f2cblaslapack --with-shared=1 --with-x=0 --with-64-bit-indices --with-c++-support=1 --with-pthread=1 --with-cc=${MPI_C} --with-cxx=${MPI_CXX} --with-fc=0
	echo "building PETSC" 
#	make PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH}
#	make install
fi

echo "export PETSC_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh
echo "export PETSC_ARCH=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh

echo "export LD_LIBRARY_PATH=${BLD_DIR}/lib:${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh

################################
# NIFTICLIB
################################
NIFTI_LIB_DIR=${BUILD_DIR}/nifticlib
SRC_DIR=${NIFTI_LIB_DIR}/nifticlib-srcs
BLD_DIR=${NIFTI_LIB_DIR}/nifticlib-build
if [ ! -d ${NIFTI_LIB_DIR} ]; then
	mkdir ${NIFTI_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}
	echo extracting PETSC lib...
	tar -xzf ${LIB_DIR}/nifticlib-2.0.0.tar.gz -C ${SRC_DIR} --strip-components=1
fi

PETSC_ARCH=cxx_opt
if [ ${builddep} -eq 1 ]; then 
	echo "configuring NIFTICLIB" 
	cd ${SRC_DIR}
	#cmake ${SRC_DIR} -DBUILD_SHARED_LIBS:BOOL=OFF -DCMAKE_INSTALL_PREFIX:PATH=${BLD_DIR}
	echo "building NIFTICLIB" 
	#make
	#make install
fi

echo "export NIFTI_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh

