#!/bin/bash

# define the MPI compilers
MPI_CXX=mpicxx
MPI_C=mpicc
#MPI_F77=mpif77
#MPI_F99=mpif90

builddep=0		# set to 1 if you wanna build all libraries
enableAVX=0		# enable AVX	# (Advanced Vector Extensions (AVX)
								# are extensions to the x86 instruction
								# set architecture for microprocessors
								# from Intel and AMD)
enableOMP=1 # needed by ACCFFT
useIMPI=0

buildfftw=0
buildaccfft=0
buildnifticlib=0
buildpetscdbl=0
buildpetscsgl=0
buildpetscdbgsgl=0
buildpetscdgbdbl=0
buildpnetcdf=0
buildzlib=0
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
    buildpetscdgbdbl=1
    shift # past argument=value
    ;;
    --bpetscsgl)
    buildpetscsgl=1
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
    echo "     --help          print this message"
    echo "     --build         build all libraries"
    echo "----------------------------------------------------------------------------------"
    echo "     --cxx=<CXX>     MPI C++ compiler (typically mpicxx; mpicxx is used if not set)"
    echo "     --c=<CC>        MPI C compiler (typically mpicc; mpicc is used if not set)"
    echo "     --bldir=<DIR>   path to your local lapack and blas installation (PETSc)"
    echo "     --useIMPI       flag: use intel MPI (instead ov MVAPICH and OpenMPI)"
    #echo "     --enableOMP     flag: use OpenMP"
    echo "     --enableAVX     flag: use AVX"
    echo "----------------------------------------------------------------------------------"
    echo " build libraries"
    echo "----------------------------------------------------------------------------------"
    echo "     --bfftw         build FFTW library"
    echo "     --baccfft       build ACCFFT library (depends on FFTW)"
    echo "     --bnifti        build NIFTI library"
    echo "     --bpetscsgl     build PETSc library (single precision)"
    echo "     --bpetscdbl     build PETSc library (double precision)"
    echo "     --bpetscdbgsgl  build PETSc library (single precision; debug mode)"
    echo "     --bpetscdbgdbl  build PETSc library (double precision; debug mode)"
    echo "     --bpnetcdf      build pnetCDF library"
    echo "----------------------------------------------------------------------------------"
    echo " clean libraries"
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


CFLAGS+=''
if [ ${useIMPI} -eq 1 ]; then
CFLAGS+=-mt_mpi
fi

CXXFLAGS+=''
if [ ${useIMPI} -eq 1 ]; then
CXXFLAGS+=-mt_mpi
fi


### PETSC OPTIONS
PETSC_OPTIONS="
--with-cc=${MPI_C}
--CFLAGS=${CFLAGS}
COPTFLAGS='-O3'
--with-cxx=${MPI_CXX}
--CXXFLAGS=${CXXFLAGS}
--download-f2cblaslapack
CXXOPTFLAGS='-O3'
--with-ssl=0
--with-debugging=0
--with-64-bit-indices
--with-shared=0
--with-x=0
--with-fc=0"

PETSC_DBG_OPTIONS="
--with-cc=${MPI_C}
--CFLAGS=${CFLAGS}
--with-cxx=${MPI_CXX}
--CXXFLAGS=${CXXFLAGS}
--with-ssl=0
--download-f2cblaslapack
--with-debugging=1
--with-64-bit-indices
--with-shared=0
--with-x=0
--with-fc=0"


#--download-hdf5
#--with-hdf5


#### FFTW OPTIONS
FFTW_OPTIONS="--enable-sse2 MAKEINFO=missing"
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

if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${FFTW_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo "----------------------------------------------------------------------------------"
		echo extracting FFTW lib...
		echo "----------------------------------------------------------------------------------"
		tar -xzf ${LIB_DIR}/fftw-3.3.6-pl1.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
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
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	cd ${SRC_DIR}
	echo ./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} CFLAGS='-O3'
	./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} CFLAGS='-O3'

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
	echo ./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-float CFLAGS='-O3'
	./configure --prefix=${BLD_DIR} ${FFTW_OPTIONS} --enable-float CFLAGS='-O3'

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
# PETSC
################################
PETSC_LIB_DIR=${BUILD_DIR}/petsc
SRC_DIR=${PETSC_LIB_DIR}/src
BLD_DIR=${PETSC_LIB_DIR}/build

echo "export PETSC_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh

if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${PETSC_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo "----------------------------------------------------------------------------------"
		echo extracting PETSC lib...
		echo "----------------------------------------------------------------------------------"
		tar -xzf ${LIB_DIR}/petsc-lite-3.7.3.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [  ${cleanup} -eq 1 -a ! ${PETSC_LIB_DIR} == ${HOME} ]; then
		rm -rf ${PETSC_LIB_DIR}
	fi
fi

PETSC_ARCH=cxx_opt_dbl
if [ ${builddep} -eq 1 -o ${buildpetscdbl} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring PETSC (double precision)"
	echo "----------------------------------------------------------------------------------"
	if [ -d ${SRC_DIR}/${PETSC_ARCH} -a ! ${SRC_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${SRC_DIR}/${PETSC_ARCH}
	fi
	if [ -d ${BLD_DIR}/${PETSC_ARCH} -a ! ${BLD_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${BLD_DIR}/${PETSC_ARCH}
	fi
	cd ${SRC_DIR}
	echo ./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH}/${PETSC_ARCH} --prefix=${BLD_DIR} ${PETSC_OPTIONS}
	./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR}/${PETSC_ARCH} ${PETSC_OPTIONS}
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

echo "export PETSC_ARCH_DOUBLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh



################################
# PETSC
################################
PETSC_ARCH=cxx_opt_sgl
if [ ${builddep} -eq 1 -o ${buildpetscsgl} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring PETSC (single precision)"
	echo "----------------------------------------------------------------------------------"
	if [ -d ${SRC_DIR}/${PETSC_ARCH} -a ! ${SRC_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${SRC_DIR}/${PETSC_ARCH}
	fi
	if [ -d ${BLD_DIR}/${PETSC_ARCH} -a ! ${BLD_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${BLD_DIR}/${PETSC_ARCH}
	fi
	cd ${SRC_DIR}
	echo ./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR}/${PETSC_ARCH} ${PETSC_OPTIONS} --with-precision=single
	./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR}/${PETSC_ARCH} ${PETSC_OPTIONS} --with-precision=single
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

echo "export PETSC_ARCH_SINGLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh




################################
# PETSC DBG
################################
PETSC_ARCH=cxx_dbg_dbl
if [ ${builddep} -eq 1 -o ${buildpetscdgbdbl} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring PETSC (double precision; debug)"
	echo "----------------------------------------------------------------------------------"
	if [ -d ${SRC_DIR}/${PETSC_ARCH} -a ! ${SRC_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${SRC_DIR}/${PETSC_ARCH}
	fi
	if [ -d ${BLD_DIR}/${PETSC_ARCH} -a ! ${BLD_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${BLD_DIR}/${PETSC_ARCH}
	fi
	cd ${SRC_DIR}
	echo ./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR}/${PETSC_ARCH} ${PETSC_DBG_OPTIONS}
	./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR}/${PETSC_ARCH} ${PETSC_DBG_OPTIONS}
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

echo "export PETSC_ARCH_DBG_DOUBLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh


################################
# PETSC DBG SINGLE PRECISION
################################
PETSC_ARCH=cxx_dbg_sgl
if [ ${builddep} -eq 1 -o ${buildpetscdbgsgl} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring PETSC (single precision; debug)"
	echo "----------------------------------------------------------------------------------"
	if [ -d ${SRC_DIR}/${PETSC_ARCH} -a ! ${SRC_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${SRC_DIR}/${PETSC_ARCH}
	fi
	if [ -d ${BLD_DIR}/${PETSC_ARCH} -a ! ${BLD_DIR}/${PETSC_ARCH} == ${HOME} ]; then
		rm -rf ${BLD_DIR}/${PETSC_ARCH}
	fi
	cd ${SRC_DIR}
	echo ./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR}/${PETSC_ARCH} ${PETSC_DBG_OPTIONS} --with-precision=single
	./configure PETSC_DIR=${SRC_DIR} PETSC_ARCH=${PETSC_ARCH} --prefix=${BLD_DIR}/${PETSC_ARCH} ${PETSC_DBG_OPTIONS} --with-precision=single
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

echo "export PETSC_ARCH_DBG_SINGLE=${PETSC_ARCH}" >> ${BUILD_DIR}/environment_vars.sh
echo "export LD_LIBRARY_PATH=${BLD_DIR}/${PETSC_ARCH}/lib:\${LD_LIBRARY_PATH}" >> ${BUILD_DIR}/environment_vars.sh


################################
# MORTON LIBRARY
################################
M_LIB_DIR=${BUILD_DIR}/libmorton
SRC_DIR=${M_LIB_DIR}
if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${M_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo "----------------------------------------------------------------------------------"
		echo extracting morton library...
		echo "----------------------------------------------------------------------------------"
		tar -xzf ${LIB_DIR}/morton.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [  ${cleanup} -eq 1 -a ! ${M_LIB_DIR} == ${HOME} ]; then
		rm -rf ${M_LIB_DIR}
	fi
fi
echo "export MORTON_DIR=${SRC_DIR}" >> ${BUILD_DIR}/environment_vars.sh




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
		echo "----------------------------------------------------------------------------------"
		echo extracting ACCFFT lib...
		echo "----------------------------------------------------------------------------------"
		tar -xzf ${LIB_DIR}/accfft.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [ ${cleanup} -eq 1 -a  ! ${ACCFFT_LIB_DIR} == ${HOME} ]; then
		rm -rf ${ACCFFT_LIB_DIR}
	fi
fi


if [ ${builddep} -eq 1 -o ${buildaccfft} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring ACCFFT"
	echo "----------------------------------------------------------------------------------"
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	mkdir ${BLD_DIR}/config
	cd ${BLD_DIR}/config
	echo cmake ${SRC_DIR} ${ACCFFT_OPTIONS} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/build -DCFLAGS='-O3'
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
# zlib
################################
Z_LIB_DIR=${BUILD_DIR}/zlib
SRC_DIR=${Z_LIB_DIR}/src
BLD_DIR=${Z_LIB_DIR}/build

if [ ! ${cleanup} -eq 1 ]; then
	if [ ! -d ${Z_LIB_DIR} -o ! -d ${SRC_DIR} ]; then
		mkdir -p ${SRC_DIR}
		echo ""
		echo "----------------------------------------------------------------------------------"
		echo extracting zlib...
		echo "----------------------------------------------------------------------------------"
		tar -xzf ${LIB_DIR}/zlib-1.2.8.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [  ${cleanup} -eq 1 -a ! ${Z_LIB_DIR} == ${HOME} ]; then
		rm -rf ${Z_LIB_DIR}
	fi
fi

if [ ${builddep} -eq 1 -o ${buildzlib} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring zlib" 
	echo "----------------------------------------------------------------------------------"
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	cd ${SRC_DIR}
	echo ./configure --prefix=${BLD_DIR} --static
	./configure --prefix=${BLD_DIR} --static
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building zlib" 
	echo "----------------------------------------------------------------------------------"
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
		echo "----------------------------------------------------------------------------------"
		echo extracting NIFTI lib...
		echo "----------------------------------------------------------------------------------"
		tar -xzf ${LIB_DIR}/nifticlib-2.0.0.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [ ${cleanup} -eq 1 -a ! ${NIFTI_LIB_DIR} == ${HOME} ]; then
		rm -rf ${NIFTI_LIB_DIR}
	fi
fi

if [ ${builddep} -eq 1 -o ${buildnifticlib} -eq 1 ]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring NIFTICLIB"
	echo "----------------------------------------------------------------------------------"
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	mkdir ${BLD_DIR}/config
	cd ${SRC_DIR}

	echo cmake ${SRC_DIR} -DCMAKE_INSTALL_PREFIX:PATH=${BLD_DIR} ${NIFTICLIB_OPTIONS} -DZLIB_ROOT=${Z_LIB_DIR}/build

	cmake ${SRC_DIR} -DCMAKE_INSTALL_PREFIX:PATH=${BLD_DIR} ${NIFTICLIB_OPTIONS} -DZLIB_ROOT=${Z_LIB_DIR}/build

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
		echo "----------------------------------------------------------------------------------"
		echo extracting PNETCDF lib...
		echo "----------------------------------------------------------------------------------"
		tar -xzf ${LIB_DIR}/parallel-netcdf-1.7.0.tar.gz -C ${SRC_DIR} --strip-components=1
		#tar -xzf ${LIB_DIR}/parallel-netcdf-1.6.1.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [ ${cleanup} -eq 1 -a ! ${PNETCDF_LIB_DIR} == ${HOME} ]; then
		rm -rf ${PNETCDF_LIB_DIR}
	fi
fi

#if [ ${builddep} -eq 1 -o ${buildnifticlib} -eq 1 ]; then
if [ ${buildpnetcdf} -eq 1 ]; then
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring PNETCDF lib..."
	echo "----------------------------------------------------------------------------------"
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	cd ${SRC_DIR}

	echo ./configure --prefix=${BLD_DIR}
	./configure --prefix=${BLD_DIR} FFLAGS='-O3' CFLAGS='-O3'

	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building PNETCDF lib"
	echo "----------------------------------------------------------------------------------"
	make
	make install
fi

echo "export PNETCDF_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh

if [ ${cleanup} -eq 1 ]; then
	rm -f ${BUILD_DIR}/environment_vars.sh
fi



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
		echo "----------------------------------------------------------------------------------"
		echo extracting ACCFFT lib...
		echo "----------------------------------------------------------------------------------"
		tar -xzf ${LIB_DIR}/accfft.tar.gz -C ${SRC_DIR} --strip-components=1
	fi
else
	if [ ${cleanup} -eq 1 -a  ! ${ACCFFT_LIB_DIR} == ${HOME} ]; then
		rm -rf ${ACCFFT_LIB_DIR}
	fi
fi


if [ ${builddep} -eq 1 -o ${buildaccfft} -eq 1 ]; then 
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "configuring ACCFFT"
	echo "----------------------------------------------------------------------------------"
	if [ -d ${BLD_DIR} -a ! ${BLD_DIR} == ${HOME} ]; then
		rm -rf ${BLD_DIR}
	fi
	mkdir ${BLD_DIR}
	mkdir ${BLD_DIR}/config
	cd ${BLD_DIR}/config
	echo cmake ${SRC_DIR} ${ACCFFT_OPTIONS} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/build -DCFLAGS='-O3'
	cmake ${SRC_DIR} ${ACCFFT_OPTIONS} -DCMAKE_INSTALL_PREFIX=${BLD_DIR} -DFFTW_ROOT=${FFTW_LIB_DIR}/build -DCFLAGS='-O3'
	echo ""
	echo "----------------------------------------------------------------------------------"
	echo "building ACCFFT"
	echo "----------------------------------------------------------------------------------"
	make -j
	make install
fi

echo "export ACCFFT_DIR=${BLD_DIR}" >> ${BUILD_DIR}/environment_vars.sh


