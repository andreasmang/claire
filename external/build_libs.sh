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
buildnifticlib=0
buildpetsc=0
uselocalmkl=0
cleanup=0

for i in "$@"
do
case $i in
    --bldir=*)
    BL_DIR="${i#*=}"
    uselocalmkl=1
    shift # past argument=value
    ;;
    --cxx=*)
    MPI_CXX="${i#*=}"
    shift # past argument=value
    ;;
    --c=*)
    MPI_C="${i#*=}"
    shift # past argument=value
    ;;
    --build)
    builddep=1
    shift # past argument=value
    ;;
    --bfftw)
    buildfftw=1
    shift # past argument=value
    ;;
    --baccfft)
    buildaccfft=1
    shift # past argument=value
    ;;
    --bnifti)
    buildnifticlib=1
    shift # past argument=value
    ;;
    --bpetsc)
    buildpetsc=1
    shift # past argument=value
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
    echo " the libraries are: FFTW; ACCFFT; NIFTICLIB; PETSc;"
    echo "----------------------------------------------------------------------------------"
    echo " options for this script are"
    echo "----------------------------------------------------------------------------------"
    echo "     --help         print this message"
    echo "     --build        build all libraries"
    echo "----------------------------------------------------------------------------------"
    echo "     --cxx=<CXX>    MPI C++ compiler (typically mpicxx; mpicxx is used if not set)"
    echo "     --c=<CC>       MPI C compiler (typically mpicc; mpicc is used if not set)"
    echo "     --bldir=<DIR>  path to your local lapack and blas installation (PETSc)"
    echo "----------------------------------------------------------------------------------"
    echo "     --bfftw        build FFTW library"
    echo "     --baccfft      build ACCFFT library (depends on FFTW & PNETCDF)"
    echo "     --bnifti       build NIFTI library"
    echo "     --bpetsc       build PETSc library"
    echo "----------------------------------------------------------------------------------"
    echo "     --clean        remove all libraries (deletes all subfolders)"
    echo "----------------------------------------------------------------------------------"
    echo ""
    exit;
    shift # past argument=value
    ;;
    *)
    # unknown option
    ;;
esac
shift
done


if [ ${uselocalmkl} -eq 1 ]; then
echo " using lapack and blas dir: ${BL_DIR}"
PETSC_OPTIONS="
--with-cc=${MPI_C}
COPTFLAGS='-O3'
--with-cxx=${MPI_CXX}
CXXOPTFLAGS='-O3'
--with-blas-lapack-dir=${BL_DIR}
--with-debugging=0
--with-x=0
--with-c++-support=1
--with-clanguage=C++
--with-fc=0
--with-64-bit-indices"
else
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
fi

#CFALGS='-mt_mpi'

FFTW_OPTIONS="
--enable-mpi
--enable-threads
--enable-sse2
--enable-openmp
MAKEINFO=missing"
#CFLAGS='-O3'

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
		./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-avx CFLAGS='-O3'
	else
		./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} CFLAGS='-O3'
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
		./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-avx --enable-float CFLAGS='-O3'
	else 
		./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-float CFLAGS='-O3'
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
	cmake ${SRC_DIR} ${ACCFFT_OPTIONS} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/build -DCFLAGS='-O3'
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
