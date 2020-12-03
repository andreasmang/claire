#!/bin/bash

# define the MPI compilers
MPI_CXX=mpicxx
MPI_C=mpicc
CUDA_C=nvcc
#MPI_F77=mpif77
#MPI_F99=mpif90

petscvar='lite-3.11.4'


builddep=0		# set to 1 if you wanna build all libraries
enableAVX=0		# enable AVX	# (Advanced Vector Extensions (AVX)
								# are extensions to the x86 instruction
								# set architecture for microprocessors
								# from Intel and AMD)
enableOMP=1 # needed by ACCFFT
enableCUDA=0
useIMPI=0
POWER9=0

GPU="V100" # TACC Longhorn
#GPU="RTX" # TACC Frontera
#GPU="P100" # TACC Maverick2

buildfftw=0
buildaccfft=0
buildnifticlib=0
buildpetscdbl=0
buildpetscsgl=0
buildpetsccudasgl=0
buildpetsccudadbl=0
buildpetsccudasgldbg=0
buildpetsccudadbldbg=0
buildpetscdbgsgl=0
buildpetscdbgdbl=0
buildpnetcdf=0
buildzlib=0
cleanup=0

myline="----------------------------------------------------------------------------------"

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
    --gpu=*)
    GPU="${i#*=}"
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
    --enableCUDA)
    enableCUDA=1
    shift # past argument=value
    ;;
    --bnifti)
    buildnifticlib=1
    shift # past argument=value
    ;;
    --bzlib)
    buildzlib=1
    shift # past argument=value
    ;;
    --bpnetcdf)
    buildpnetcdf=1
    shift # past argument=value
    ;;
    --bpetscdbl)
    buildpetscdbl=1
    shift # past argument=value
    ;;
    --bpetscdbgdbl)
    buildpetscdbgdbl=1
    shift # past argument=value
    ;;
    --bpetscsgl)
    buildpetscsgl=1
    shift # past argument=value
    ;;
    --bpetsccudasgl)
    buildpetsccudasgl=1
    shift # past argument=value
    ;;
    --bpetsccudasgldbg)
    buildpetsccudasgldbg=1
    shift # past argument=value
    ;;
    --bpetsccudadbl)
    buildpetsccudadbl=1
    shift # past argument=value
    ;;
    --bpetsccudadbldbg)
    buildpetsccudadbldbg=1
    shift # past argument=value
    ;;
    --bpetscdbgsgl)
    buildpetscdbgsgl=1
    shift # past argument=value
    ;;
    --useIMPI)
    useIMPI=1
    shift # past argument=value
    ;;
    --enableOMP)
    enableOMP=1
    shift # past argument=value
    ;;
    --enableAVX)
    enableAVX=1
    shift # past argument=value
    ;;
    --POWER9)
    POWER9=1
    shift # past argument=value
    ;;
    --clean)
    cleanup=1
    echo ""
    echo ${myline}
    echo " cleaning up"
    echo ${myline}
    ;;
    --help)
    echo "script to build dependencies (libraries)" 
    echo ${myline}
    echo " the libraries are: FFTW; ACCFFT; NIFTICLIB; PETSc; MORTON, PNETCDF"
    echo ${myline}
    echo " options for this script are"
    echo ${myline}
    echo "     --help          print this message"
    echo "     --build         build all libraries"
    echo ${myline}
    echo "     --cxx=<CXX>     MPI C++ compiler (typically mpicxx; mpicxx is used if not set)"
    echo "     --c=<CC>        MPI C compiler (typically mpicc; mpicc is used if not set)"
    echo "     --bldir=<DIR>   path to your local lapack and blas installation (PETSc)"
    echo "     --useIMPI       flag: use intel MPI (instead ov MVAPICH and OpenMPI)"
    #echo "     --enableOMP     flag: use OpenMP"
    echo "     --enableAVX     flag: use AVX"
    echo "     --enableCUDA    flag: use CUDA for AccFFT"
    echo "     --POWER9        flag: build for IBM POWER9 architecture"
    echo ${myline}
    echo " build libraries"
    echo ${myline}
    echo "     --bfftw         build FFTW library (mandatory)"
    echo "     --baccfft       build ACCFFT library (mandatory; depends on FFTW)"
    echo "     --bnifti        build NIFTI library (default)"
    echo "     --bpetscsgl     build PETSc library (mandatory; default; single precision)"
    echo "     --bpetscdbl     build PETSc library (double precision)"
    echo "     --bpetscdbgsgl  build PETSc library (for developers; single precision; debug mode)"
    echo "     --bpetscdbgdbl  build PETSc library (for developers; double precision; debug mode)"
    echo "     --bpetsccudasgl build PETSc library with nVidia-CUDA library (single precision)"
    echo "     --bpetsccudasgldbg build PETSc library with nVidia-CUDA library (single precision; debug mode)"
    echo "     --bpnetcdf      build pnetCDF library (optional)"
    echo ${myline}
    echo " clean libraries"
    echo ${myline}
    echo "     --clean        remove all libraries (deletes all subfolders)"
    echo ${myline}
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


CFLAGS+=''
if [ ${useIMPI} -eq 1 ]; then
CFLAGS+=-mt_mpi
fi

CXXFLAGS+=''
if [ ${useIMPI} -eq 1 ]; then
CXXFLAGS+=-mt_mpi
fi

#--download-hdf5
#--with-hdf5


#### FFTW OPTIONS
FFTW_OPTIONS="--enable-sse2 MAKEINFO=missing"
if [ ${POWER9} -eq 1 ]; then
  FFTW_OPTIONS="MAKEINFO=missing"
fi
if [ ${enableOMP} -eq 1 ]; then
	FFTW_OPTIONS="${FFTW_OPTIONS} --enable-threads --enable-openmp"
fi
if [ ${enableAVX} -eq 1 ]; then
	FFTW_OPTIONS="${FFTW_OPTIONS} --enable-avx"
fi

### ACCFFT OPTIONS
ACCFFT_OPTIONS="
-DCMAKE_CXX_COMPILER=${MPI_CXX}
-DCMAKE_C_COMPILER=${MPI_C}
-DFFTW_USE_STATIC_LIBS=true
-DBUILD_GPU=false
-DBUILD_SHARED=false"

if [ ${enableCUDA} -eq 1 ]; then
  if [ "$GPU" == "V100" ]; then
    ACCFFT_OPTIONS="${ACCFFT_OPTIONS} -DBUILD_GPU=true -DCUDA_NVCC_FLAGS=-gencode;arch=compute_70,code=sm_70"
  fi
  if [ "$GPU" == "P100" ]; then
    ACCFFT_OPTIONS="${ACCFFT_OPTIONS} -DBUILD_GPU=true -DCUDA_NVCC_FLAGS=-gencode;arch=compute_60,code=sm_60"
  fi
  if [ "$GPU" == "RTX" ]; then
    ACCFFT_OPTIONS="${ACCFFT_OPTIONS} -DBUILD_GPU=true -DCUDA_NVCC_FLAGS=-gencode;arch=compute_75,code=sm_75"
  fi
fi


### NIFTI OPTIONS
NIFTICLIB_OPTIONS="
-DCMAKE_CXX_COMPILER=${MPI_CXX}
-DCMAKE_C_COMPILER=${MPI_C}
-DBUILD_SHARED_LIBS:BOOL=OFF
-Wno-dev"

##############################################################
##############################################################
# DO NOT MODIFY
##############################################################
##############################################################

#LIB_DIR=${HOME}/apps

LIB_DIR=$PWD/

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

echo "export MPI_DIR=${MPI_DIR}" >> ${BUILD_DIR}/environment_vars.sh

if [ ${enableCUDA} -eq 1 ]; then
  CUDA_DIR=$(which ${CUDA_C})
  CUDA_DIR=$(dirname "${CUDA_DIR}")
  CUDA_DIR=$(dirname "${CUDA_DIR}")
  echo " detected CUDA-toolkit directory: ${CUDA_DIR}"
  echo "export CUDA_DIR=${CUDA_DIR}" >> ${BUILD_DIR}/environment_vars.sh
  echo "export LD_LIBRARY_PATH=${CUDA_DIR}/lib64:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh
fi

################################
# FFTW
################################
FFTW_LIB_DIR=${BUILD_DIR}/fftw
SRC_DIR=${FFTW_LIB_DIR}/src
BLD_DIR=${FFTW_LIB_DIR}/build

if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${FFTW_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo ${myline} 
		echo extracting FFTW lib...
		echo ${myline} 
		tar -xzf ${LIB_DIR}/fftw.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [ ${cleanup} -eq 1 -a ! ${FFTW_LIB_DIR} == ${HOME} ]; then
		rm -rf ${FFTW_LIB_DIR}
	fi
fi

if [ ${builddep} -eq 1 -o ${buildfftw} -eq 1 ]; then 
	echo ""
	echo ${myline} 
	echo "configuring FFTW (double precision)"
	echo ${myline} 
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	cd ${SRC_DIR}
	echo ./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} CFLAGS='-O3'
	./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} CFLAGS='-O3'

	echo ""
	echo ${myline} 
	echo "building FFTW (double precision)" 
	echo ${myline} 
	make
	make install

	echo ""
	echo ${myline} 
	echo "configuring FFTW (float precision)"
	echo ${myline} 
	echo ./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-float CFLAGS='-O3'
	./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-float CFLAGS='-O3'

	echo ""
	echo ${myline} 
	echo "building FFTW (float precision)" 
	echo ${myline} 
	make
	make install
fi

echo "export FFTW_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh


################################
# PETSC
################################
function build_petsc() {
  PETSC_ARCH=$1
  MESSAGE=$2
	echo ""
	echo ${myline} 
	echo "${MESSAGE}"
	echo ${myline} 
	if [ -d ${SRC_DIR}/${PETSC_ARCH} -a ! ${SRC_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${SRC_DIR}/${PETSC_ARCH}
	fi
	if [ -d ${BLD_DIR}/${PETSC_ARCH} -a ! ${BLD_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${BLD_DIR}/${PETSC_ARCH}
	fi
	
	cd ${SRC_DIR}
	python2 config/examples/configure_petsc.py $PETSC_ARCH $GPU
	make PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} all
}

# extract the source files
PETSC_LIB_DIR=${BUILD_DIR}/petsc-${petscvar}
SRC_DIR=${PETSC_LIB_DIR}
BLD_DIR=${PETSC_LIB_DIR}
export PETSC_DIR=${SRC_DIR}

echo "export PETSC_DIR=${PETSC_LIB_DIR}" >> ${BUILD_DIR}/environment_vars.sh

if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${PETSC_LIB_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo ${myline} 
		echo extracting PETSC lib...
		echo ${myline} 
		tar -xzf ${LIB_DIR}/petsc-${petscvar}.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [  ${cleanup} -eq 1 -a ! ${PETSC_LIB_DIR} == ${HOME} ]; then
		rm -rf ${PETSC_LIB_DIR}
	fi
fi

if [ -d ${SRC_DIR} ]; then
  cp ${LIB_DIR}/configure_petsc.py ${SRC_DIR}/config/examples/
fi
################################
# PETSC-DOUBLE-PRECISION
################################
PETSC_ARCH=cxx-opt-dbl
if [ ${builddep} -eq 1 -o ${buildpetscdbl} -eq 1 ]; then 
	build_petsc $PETSC_ARCH "configuring PETSC (double precision)"
fi
echo "export PETSC_ARCH_DOUBLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh

################################
# PETSC-SINGLE-PRECISION
################################
PETSC_ARCH=cxx-opt-sgl
if [ ${builddep} -eq 1 -o ${buildpetscsgl} -eq 1 ]; then 
	build_petsc $PETSC_ARCH "configuring PETSC (single precision)"
fi
echo "export PETSC_ARCH_SINGLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh

################################
# PETSC-CUDA-SINGLE-PRECISION
################################
PETSC_ARCH=cuda-opt-sgl-${GPU}
if [ ${builddep} -eq 1 -o ${buildpetsccudasgl} -eq 1 ]; then 
	build_petsc $PETSC_ARCH "configuring PETSC-CUDA (single precision)"
fi
echo "export PETSC_ARCH_CUDA_SINGLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh

################################
# PETSC-DOUBLE-DEBUG
################################
PETSC_ARCH=cxx-opt-dbg-dbl
if [ ${builddep} -eq 1 -o ${buildpetscdbgdbl} -eq 1 ]; then 
	build_petsc $PETSC_ARCH "configuring PETSC (double precision) in debug mode"
fi
echo "export PETSC_ARCH_DBG_DOUBLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh

################################
# PETSC DBG SINGLE PRECISION
################################
PETSC_ARCH=cxx-opt-dbg-sgl
if [ ${builddep} -eq 1 -o ${buildpetscdbgsgl} -eq 1 ]; then 
	build_petsc $PETSC_ARCH "configuring PETSC (single precision) in debug mode"
fi
echo "export PETSC_ARCH_DBG_SINGLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh

################################
# PETSC-CUDA DBG SINGLE
################################
PETSC_ARCH=cuda-opt-sgl-dbg-${GPU}
if [ ${builddep} -eq 1 -o ${buildpetsccudasgldbg} -eq 1 ]; then 
	build_petsc $PETSC_ARCH "configuring PETSC-cuda (single precision) in debug mode"
fi
echo "export PETSC_ARCH_CUDA_SINGLE_DBG=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh

################################
# MORTON LIBRARY
################################
#M_LIB_DIR=${BUILD_DIR}/libmorton
#SRC_DIR=${M_LIB_DIR}
#if [ ! ${cleanup} -eq 1 ]; then
#	if [ ! -d ${M_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
#		mkdir -p ${SRC_DIR}
#		echo ""
#		echo ${myline} 
#		echo extracting morton library...
#		echo ${myline} 
#		tar -xzf ${LIB_DIR}/morton.tar.gz -C ${SRC_DIR} --strip-components=1
#	fi
#else
#	if [  ${cleanup} -eq 1 -a ! ${M_LIB_DIR} == ${HOME} ]; then
#		rm -rf ${M_LIB_DIR}
#	fi
#fi
#echo "export MORTON_DIR=${SRC_DIR}" >> ${BUILD_DIR}/environment_vars.sh




################################
# zlib
################################
Z_LIB_DIR=${BUILD_DIR}/zlib
SRC_DIR=${Z_LIB_DIR}/src
BLD_DIR=${Z_LIB_DIR}/build

if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${Z_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo ${myline} 
		echo extracting zlib...
		echo ${myline} 
		tar -xzf ${LIB_DIR}/zlib.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [  ${cleanup} -eq 1 -a ! ${Z_LIB_DIR} == ${HOME} ]; then
		rm -rf ${Z_LIB_DIR}
	fi
fi

if [ ${builddep} -eq 1 -o ${buildzlib} -eq 1 ]; then 
	echo ""
	echo ${myline} 
	echo "configuring zlib" 
	echo ${myline} 
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	cd ${SRC_DIR}
	echo ./configure --prefix=${BLD_DIR} --static
	./configure --prefix=${BLD_DIR} --static
	echo ""
	echo ${myline} 
	echo "building zlib" 
	echo ${myline} 
	make 
	make install
else
	if [  ${cleanup} -eq 1 -a ! ${Z_LIB_DIR} == ${HOME} ]; then
		rm -rf ${Z_LIB_DIR}
	fi
fi
echo "export ZLIB_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh


################################
# NIFTICLIB
################################
NIFTI_LIB_DIR=${BUILD_DIR}/nifticlib
SRC_DIR=${NIFTI_LIB_DIR}/src
BLD_DIR=${NIFTI_LIB_DIR}/build
if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${NIFTI_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo ${myline} 
		echo extracting NIFTI lib...
		echo ${myline} 
		#tar -xzf ${LIB_DIR}/nifticlib.tar.gz -C ${SRC_DIR} --strip-components=1
		tar -xf ${LIB_DIR}/nifticlib-2.0.0.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [ ${cleanup} -eq 1 -a ! ${NIFTI_LIB_DIR} == ${HOME} ]; then
		rm -rf ${NIFTI_LIB_DIR}
	fi
fi

if [ ${builddep} -eq 1 -o ${buildnifticlib} -eq 1 ]; then
	echo ""
	echo ${myline} 
	echo "configuring NIFTICLIB"
	echo ${myline} 
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	mkdir ${BLD_DIR}/config
	cd ${SRC_DIR}

	echo cmake ${SRC_DIR} -DCMAKE_INSTALL_PREFIX:PATH=${BLD_DIR} ${NIFTICLIB_OPTIONS} -DZLIB_ROOT=${Z_LIB_DIR}/build

	cmake ${SRC_DIR} -DCMAKE_INSTALL_PREFIX:PATH=${BLD_DIR} ${NIFTICLIB_OPTIONS} -DZLIB_ROOT=${Z_LIB_DIR}/build

	echo ""
	echo ${myline} 
	echo "building NIFTICLIB"
	echo ${myline} 
	make
	make install
fi

if [ ${cleanup} -eq 1 ]; then
	rm -f ${BUILD_DIR}/environment_vars.sh
fi
echo "export NIFTI_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh

################################
# PNETCDF
################################
PNETCDF_LIB_DIR=${BUILD_DIR}/pnetcdflib
SRC_DIR=${PNETCDF_LIB_DIR}/src
BLD_DIR=${PNETCDF_LIB_DIR}/build
if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${PNETCDF_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo ${myline} 
		echo extracting PNETCDF lib...
		echo ${myline} 
		tar -xzf ${LIB_DIR}/parallel-netcdf.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [ ${cleanup} -eq 1 -a ! ${PNETCDF_LIB_DIR} == ${HOME} ]; then
		rm -rf ${PNETCDF_LIB_DIR}
	fi
fi

#if [ ${builddep} -eq 1 -o ${buildnifticlib} -eq 1 ]; then
if [ ${buildpnetcdf} -eq 1 ]; then
	echo ""
	echo ${myline} 
	echo "configuring PNETCDF lib..."
	echo ${myline} 
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	cd ${SRC_DIR}

	echo ./configure --prefix=${BLD_DIR}
	./configure --prefix=${BLD_DIR} FFLAGS='-O3' CFLAGS='-O3'

	echo ""
	echo ${myline} 
	echo "building PNETCDF lib"
	echo ${myline} 
	make
	make install
fi

if [ ${cleanup} -eq 1 ]; then
	rm -f ${BUILD_DIR}/environment_vars.sh
fi
echo "export PNETCDF_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh

################################
# ACCFFT
################################
ACCFFT_LIB_DIR=${BUILD_DIR}/accfft
SRC_DIR=${ACCFFT_LIB_DIR}/src
BLD_DIR=${ACCFFT_LIB_DIR}/build

if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${ACCFFT_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${ACCFFT_LIB_DIR}/src
		echo ""
		echo ${myline} 
		echo extracting ACCFFT lib...
		echo ${myline} 
		tar -xzf ${LIB_DIR}/accfft.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [ ${cleanup} -eq 1 -a  ! ${ACCFFT_LIB_DIR} == ${HOME} ]; then
		rm -rf ${ACCFFT_LIB_DIR}
	fi
fi


if [ ${builddep} -eq 1 -o ${buildaccfft} -eq 1 ]; then 
	echo ""
	echo ${myline} 
	echo "configuring ACCFFT"
	echo ${myline} 
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	mkdir ${BLD_DIR}/config
	cd ${BLD_DIR}/config
	echo cmake ${SRC_DIR} ${ACCFFT_OPTIONS} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/build -DCFLAGS='-O3'
	cmake ${SRC_DIR} ${ACCFFT_OPTIONS} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/build -DCFLAGS='-O3'
	echo ""
	echo ${myline} 
	echo "building ACCFFT"
	echo ${myline} 
	make -j
	make install
fi

if [ ${cleanup} -eq 1 ]; then
	rm -f ${BUILD_DIR}/environment_vars.sh
fi
echo "export ACCFFT_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh

