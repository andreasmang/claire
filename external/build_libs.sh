#!/bin/bash

# define the MPI compilers
MPI_CXX=mpicxx
MPI_C=mpicc
#MPI_F77=mpif77
#MPI_F99=mpif90

builddep=0		# set to 1 if you wanna build all libraries
enableavx=0		# enable AVX	# (Advanced Vector Extensions (AVX)
								# are extensions to the x86 instruction
								# set architecture for microprocessors
								# from Intel and AMD)

buildfftw=0
buildaccfft=0
buildpnetcdf=0
buildnifticlib=0
buildpetsc=0
cleanup=0

for i in "$@"
do
case $i in
    --cxx=*)
    MPI_CXX="${i#*=}"
    shift # past argument=value
    ;;
    --c=*)
    MPI_C="${i#*=}"
    ;;
    --build)
    builddep=1
    ;;
    --bfftw)
    buildfftw=1
    ;;
    --baccfft)
    buildaccfft=1
    ;;
    --bpnetcdf)
    buildpnetcdf=1
    ;;
    --bnifti)
    buildnifticlib=1
    ;;
    --bpetsc)
    buildpetsc=1
    ;;
    --clean)
    cleanup=1
    echo ""
    echo "----------------------------------------------------------------------------------"
    echo " cleaning up"
    echo "----------------------------------------------------------------------------------"
    ;;
    --help)
    echo "script to build dependencies (libraries)" 
    echo "----------------------------------------------------------------------------------"
    echo " the libraries are: FFTW; ACCFFT; PNETCDF; NIFTICLIB; PETSc;"
    echo "----------------------------------------------------------------------------------"
    echo " options for this script are"
    echo "----------------------------------------------------------------------------------"
    echo "     --help         print this message"
    echo "     --build        build all libraries"
    echo "----------------------------------------------------------------------------------"
    echo "     --cxx          MPI C++ compiler (typically mpicxx; mpicxx is used if not set)"
    echo "     --c            MPI C compiler (typically mpicc; mpicc is used if not set)"
    echo "----------------------------------------------------------------------------------"
    echo "     --bpnetcdf     build PNETCDF library"
    echo "     --bfftw        build FFTW library"
    echo "     --baccfft      build ACCFFT library (depends on FFTW & PNETCDF)"
    echo "     --bnifti       build NIFTI library"
    echo "     --bpetsc       build PETSc library"
    echo "----------------------------------------------------------------------------------"
    echo "     --clean        remove all libraries (deletes all subfolders)"
    echo "----------------------------------------------------------------------------------"
    echo ""
    ;;
    *)
    # unknown option
    ;;
esac
shift
done


PETSC_OPTIONS="
--with-cc=${MPI_C}
COPTFLAGS='-O3'
--with-cxx=${MPI_CXX}
CXXOPTFLAGS='-O3'
--download-f2cblaslapack
--with-debugging=0
--with-x=0
--with-c++-support=1
--with-clanguage=C++
--with-fc=0
--with-64-bit-indices"


FFTW_OPTIONS="
 --enable-mpi
--enable-threads
--enable-sse2
--enable-openmp
CFLAGS=-O3
MAKEINFO=missing"


PNETCDF_OPTIONS="
--disable-fortran"


ACCFFT_OPTIONS="
-DFFTW_USE_STATIC_LIBS=true
-DBUILD_GPU=false
-DBUILD_SHARED=false"


NIFTICLIB_OPTIONS="
-DBUILD_SHARED_LIBS:BOOL=OFF
-Wno-dev"

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
cd ${MPI_DIR}
cd ..
MPI_DIR=${PWD}
echo " detected MPI directory: ${MPI_DIR}"




if [ -e ${BUILD_DIR}/environment_vars.sh ]; then
	rm ${BUILD_DIR}/environment_vars.sh
fi


################################
# FFTW
################################
FFTW_LIB_DIR=${BUILD_DIR}/fftw
SRC_DIR=${FFTW_LIB_DIR}/src
BLD_DIR=${FFTW_LIB_DIR}/build

if [ ! -d ${FFTW_LIB_DIR} -a ! ${cleanup} -eq 1 ]; then
	mkdir ${FFTW_LIB_DIR}
	mkdir ${FFTW_LIB_DIR}/src
	mkdir ${FFTW_LIB_DIR}/build
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo extracting FFTW lib...
	echo "----------------------------------------------------------------------------------"
	tar -xzf ${LIB_DIR}/fftw-3.3.4.tar.gz -C ${SRC_DIR} --strip-components=1
else
	if [ ${cleanup} -eq 1 -a ! ${FFTW_LIB_DIR} == ${HOME} ]; then
		rm -rf ${FFTW_LIB_DIR}
	fi
fi

if [ ${builddep} -eq 1 -o ${buildfftw} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring FFTW (double precision)" 
	echo "----------------------------------------------------------------------------------"
	cd ${SRC_DIR}
	if [ ${enableavx} -eq 1 ]; then
		./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-avx
	else
		./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS}
	fi

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building FFTW (double precision)" 
	echo "----------------------------------------------------------------------------------"
	make
	make install

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring FFTW (float precision)" 
	echo "----------------------------------------------------------------------------------"
	if [ ${enableavx} -eq 1 ]; then
		./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-avx --enable-float
	else 
		./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-float
	fi

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building FFTW (float precision)" 
	echo "----------------------------------------------------------------------------------"
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
if [ ! -d ${PNETCDF_LIB_DIR} -a ! ${cleanup} -eq 1 ]; then
	mkdir ${PNETCDF_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo extracting PNETCDF lib...
	echo "----------------------------------------------------------------------------------"
#	tar -xzf ${LIB_DIR}/parallel-netcdf-1.7.0.tar.gz -C ${SRC_DIR} --strip-components=1
	tar -xzf ${LIB_DIR}/parallel-netcdf-1.6.1.tar.gz -C ${SRC_DIR} --strip-components=1
else
	if [ ${cleanup} -eq 1 -a ! ${PNETCDF_LIB_DIR} == ${HOME} ]; then
		rm -rf ${PNETCDF_LIB_DIR}
	fi
fi


if [ ${builddep} -eq 1 -o ${buildpnetcdf} -eq 1 ]; then 

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring PNETCDF" 
	echo "----------------------------------------------------------------------------------"
	cd ${SRC_DIR}
	./configure --prefix=${BLD_DIR} ${PNETCDF_OPTIONS}

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building PNETCDF" 
	echo "----------------------------------------------------------------------------------"
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
if [ ! -d ${ACCFFT_LIB_DIR} -a ! ${cleanup} -eq 1 ]; then
	mkdir ${ACCFFT_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}
	mkdir ${CMK_DIR}

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo extracting ACCFFT lib...
	echo "----------------------------------------------------------------------------------"
	tar -xzf ${LIB_DIR}/accfft.tar.gz -C ${SRC_DIR} --strip-components=1
else
	if [  ${cleanup} -eq 1  -a  ! ${ACCFFT_LIB_DIR} == ${HOME}  ]; then
		rm -rf ${ACCFFT_LIB_DIR}
	fi
fi


if [ ${builddep} -eq 1 -o ${buildaccfft} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring ACCFFT" 
	echo "----------------------------------------------------------------------------------"
	cd ${CMK_DIR}
	cmake ${SRC_DIR} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/build -DPNETCDF_DIR=${PNETCDF_LIB_DIR}/build ${ACCFFT_OPTIONS}
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building ACCFFT"
	echo "----------------------------------------------------------------------------------"
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
if [ ! -d ${PETSC_LIB_DIR} -a ! ${cleanup} -eq 1 ]; then
	mkdir ${PETSC_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo extracting PETSC lib...
	echo "----------------------------------------------------------------------------------"
	tar -xzf ${LIB_DIR}/petsc-lite-3.7.0.tar.gz -C ${SRC_DIR} --strip-components=1
fi

PETSC_ARCH=cxx_opt
if [ ${builddep} -eq 1 -o ${buildpetsc} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring PETSC" 
	echo "----------------------------------------------------------------------------------"
	cd ${SRC_DIR}
	./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR} ${PETSC_OPTIONS}

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building PETSC" 
	echo "----------------------------------------------------------------------------------"
	make PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH}
	make PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} install
else
	if [ ${cleanup} -eq 1 -a ! ${PETSC_LIB_DIR} == ${HOME} ]; then
		rm -rf ${PETSC_LIB_DIR}
	fi
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
if [ ! -d ${NIFTI_LIB_DIR} -a ! ${cleanup} -eq 1 ]; then
	mkdir ${NIFTI_LIB_DIR}
	mkdir ${BLD_DIR}
	mkdir ${SRC_DIR}
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo extracting NIFTI lib...
	echo "----------------------------------------------------------------------------------"
	tar -xzf ${LIB_DIR}/nifticlib-2.0.0.tar.gz -C ${SRC_DIR} --strip-components=1
else
	if [ ${cleanup} -eq 1 -a ! ${NIFTI_LIB_DIR} == "~/" ]; then
		rm -rf ${NIFTI_LIB_DIR}
	fi
fi

if [ ${builddep} -eq 1 -o ${buildnifticlib} -eq 1 ]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring NIFTICLIB"
	echo "----------------------------------------------------------------------------------"
	cd ${SRC_DIR}
	cmake ${SRC_DIR} -DCMAKE_INSTALL_PREFIX:PATH=${BLD_DIR} ${NIFTICLIB_OPTIONS}
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building NIFTICLIB"
	echo "----------------------------------------------------------------------------------"
	make
	make install
fi

echo "export NIFTI_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh

if [ ${cleanup} -eq 1 ]; then
	rm -f ${BUILD_DIR}/environment_vars.sh
fi
