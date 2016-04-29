#!/bin/bash

# define the MPI compilers
MPI_CXX=mpicxx
MPI_C=mpicc

builddep=0		# set to 1 if you wanna build all libraries
enableavx=0		# enable AVX	# (Advanced Vector Extensions (AVX)
								# are extensions to the x86 instruction
								# set architecture for microprocessors
								# from Intel and AMD)

for i in "$@"
do
case $i in
    -mpicxx=*)
    MPI_CXX="${i#*=}"
    shift # past argument=value
    ;;
    -mpic=*)
    MPI_C="${i#*=}"
    ;;
    --build)
    builddep=1
    ;;
    *)
    # unknown option
    ;;
esac
shift
done


##############################################################
##############################################################
# DO NOT MODIFY
##############################################################
##############################################################

LIB_DIR=${PWD}

# go up one level
BUILD_DIR=${LIB_DIR}/libs
if [ ! -d ${BUILD_DIR} ]; then
	mkdir ${BUILD_DIR}
fi

cd ${LIB_DIR}

MPI_DIR=$(which ${MPI_CXX})
MPI_DIR=$(dirname "${MPI_DIR}")
echo detected MPI directory: $MPI_DIR




if [ -e ${BUILD_DIR}/environment_vars.sh ]; then
	rm ${BUILD_DIR}/environment_vars.sh
fi




################################
# FFTW
################################
FFTW_LIB_DIR=${BUILD_DIR}/fftw
SRC_DIR=${FFTW_LIB_DIR}/src
BLD_DIR=${FFTW_LIB_DIR}/build

if [ ! -d ${FFTW_LIB_DIR} ]; then
	mkdir ${FFTW_LIB_DIR}
	mkdir ${FFTW_LIB_DIR}/src
	mkdir ${FFTW_LIB_DIR}/build
	echo extracting FFTW lib...
	tar -xzf ${LIB_DIR}/fftw-3.3.4.tar.gz -C ${SRC_DIR} --strip-components=1
fi

if [ ${builddep} -eq 1 ]; then 
	echo "configuring FFTW (double precision)" 
	cd ${SRC_DIR}
	if [ ${enableavx} -eq 1 ]; then
		./configure --prefix=${BLD_DIR} --enable-mpi --enable-threads --enable-sse2 --enable-openmp --enable-avx CFLAGS=-O3 MAKEINFO=missing
	else
		./configure --prefix=${BLD_DIR} --enable-mpi --enable-threads --enable-sse2 --enable-openmp CFLAGS=-O3 MAKEINFO=missing
	fi
	make
	make install

	echo "configuring FFTW (float precision)" 
	if [ ${enableavx} -eq 1 ]; then
		./configure --prefix=${BLD_DIR} --enable-mpi --enable-threads --enable-sse2 --enable-openmp --enable-avx --enable-float CFLAGS=-O3 MAKEINFO=missing
	else 
		./configure --prefix=${BLD_DIR} --enable-mpi --enable-threads --enable-sse2 --enable-openmp --enable-float CFLAGS=-O3 MAKEINFO=missing
	fi
	make
	make install

fi

echo "export FFTW_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh




################################
# PNETCDF
################################
PNETCDF_LIB_DIR=${BUILD_DIR}/pnetcdf
SRC_DIR=${PNETCDF_LIB_DIR}/src
BLD_DIR=${PNETCDF_LIB_DIR}/build
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
	./configure --prefix=${BLD_DIR} --with-mpi=${MPI_DIR} MPICC=${MPI_C} MPICXX=${MPI_CXX} CXXFLAGS='-O3' CPPFLAGS='-O3' --disable-fortran #FFLAGS='-O3' FCFLAGS='-O3'  
	echo "building PNETCDF" 
	make -j
	make install
fi


echo "export PNETCDF_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh




################################
# ACCFFT
################################
ACCFFT_LIB_DIR=${BUILD_DIR}/accfft
SRC_DIR=${ACCFFT_LIB_DIR}/src
BLD_DIR=${ACCFFT_LIB_DIR}/build
CMK_DIR=${ACCFFT_LIB_DIR}/build/config
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
	cmake ${SRC_DIR} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/build -DFFTW_USE_STATIC_LIBS=true -DBUILD_GPU=false -DBUILD_STEPS=true -DCXX_FLAGS="-O3" -DBUILD_SHARED=false -DPNETCDF_DIR=${PNETCDF_LIB_DIR}/build
	echo "building ACCFFT" 
	make -j
	make install
fi

echo "export ACCFFT_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh




################################
# PETSC
################################
PETSC_LIB_DIR=${BUILD_DIR}/petsc
SRC_DIR=${PETSC_LIB_DIR}/src
BLD_DIR=${PETSC_LIB_DIR}/build
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
	./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR} --with-clanguage=C++ --with-debugging=0 --download-f2cblaslapack --with-shared=1 --with-x=0 --with-64-bit-indices --with-c++-support=1 --with-pthread=1 COPTFLAGS='-O3' CXXOPTFLAGS='-O3' --with-cc=${MPI_C} --with-cxx=${MPI_CXX} --with-fc=0
	echo "building PETSC" 
	make PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH}
	make install
fi

echo "export PETSC_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh
echo "export PETSC_ARCH=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh

echo "export LD_LIBRARY_PATH=${BLD_DIR}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh





################################
# NIFTICLIB
################################
NIFTI_LIB_DIR=${BUILD_DIR}/nifticlib
SRC_DIR=${NIFTI_LIB_DIR}/srcs
BLD_DIR=${NIFTI_LIB_DIR}/build
if [ ! -d ${NIFTI_LIB_DIR} ]; then
	mkdir ${NIFTI_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}
	echo extracting NIFTI lib...
	tar -xzf ${LIB_DIR}/nifticlib-2.0.0.tar.gz -C ${SRC_DIR} --strip-components=1
fi

if [ ${builddep} -eq 1 ]; then 
	echo "configuring NIFTICLIB" 
	cd ${SRC_DIR}
	cmake ${SRC_DIR} -DBUILD_SHARED_LIBS:BOOL=OFF -DCMAKE_INSTALL_PREFIX:PATH=${BLD_DIR} -Wno-dev
	echo "building NIFTICLIB" 
	make
	make install
fi

echo "export NIFTI_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh


