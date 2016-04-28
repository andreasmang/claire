INCLUDE(ExternalProject)


# define path for external libraries
SET(COLD_EXTERNAL_LIBDIR "${CMAKE_CURRENT_BINARY_DIR}/external")




##################################
# FFTW
##################################
SET(COLD_FFTW_DIR "${COLD_EXTERNAL_LIBDIR}/fftw")
SET(ENV{FFTW_DIR} "${COLD_FFTW_DIR}")
SET(COLD_FFTW_VERSION 3.3.4)
SET(COLD_FFTW_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")

MESSAGE( STATUS "COLD_FFTW_DIR: " ${COLD_FFTW_DIR} )

IF(NOT COLD_FFTW_LIBPRJ)

    EXTERNALPROJECT_ADD(COLD_FFTW_LIBPRJ
      URL ftp://ftp.fftw.org/pub/fftw/fftw-${COLD_FFTW_VERSION}.tar.gz
      URL_MD5 ${COLD_FFTW_MD5}
      SOURCE_DIR ${COLD_FFTW_DIR}/src
      BINARY_DIR ${COLD_FFTW_DIR}/src
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/fftw
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/fftw
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/fftw
      INSTALL_DIR ${COLD_FFTW_DIR}/build
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND ./configure
            --prefix=${COLD_FFTW_DIR}/build
            --enable-mpi --enable-threads --enable-shared
            --enable-sse2 --enable-openmp --enable-avx
            CFLAGS=-O3 MAKEINFO=missing
      BUILD_COMMAND make)

    # ACCFFT also needs the float version; we can't make it at once, so
    # we have got to add an additional three steps for the float version
    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIBPRJ lfftw-config
      DEPENDEES install
      COMMAND cd ${COLD_FFTW_DIR}/src && ./configure
            --prefix=${COLD_FFTW_DIR}/build
            --enable-mpi --enable-threads --enable-shared
            --enable-sse2 --enable-openmp --enable-avx
            --enable-float CFLAGS=-O3 MAKEINFO=missing)

    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIBPRJ lfftw-build
      DEPENDEES lfftw-config
      COMMAND cd ${COLD_FFTW_DIR}/src && make)

    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIBPRJ lfftw-install
      DEPENDEES lfftw-build
      COMMAND cd ${COLD_FFTW_DIR}/src && make install)

ENDIF(NOT COLD_FFTW_LIBPRJ)

## find library
FIND_LIBRARY(COLD_FFTW_LIB NAMES fftw3 PATHS ${COLD_FFTW_DIR}/build/lib DOC "FFTW libs")
FIND_LIBRARY(COLD_FFTW_THREADS_LIB NAMES fftw3_threads PATHS ${COLD_FFTW_DIR}/build/lib DOC "FFTW threads libs")
FIND_LIBRARY(COLD_FFTWF_LIB NAMES fftw3f PATHS ${COLD_FFTW_DIR}/build/lib DOC "FFTWF libs")
FIND_LIBRARY(COLD_FFTWF_THREADS_LIB NAMES fftw3f_threads PATHS ${COLD_FFTW_DIR}/build/lib DOC "FFTWF threads libs")

LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_FFTW_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_FFTW_THREADS_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_FFTWF_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_FFTWF_THREADS_LIB})
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_FFTW_DIR}/build/include)

INCLUDE_DIRECTORIES(${COLD_FFTW_DIR}/build/include)




##################################
# PNETCDF
##################################
SET(COLD_PNETCDF_DIR "${COLD_EXTERNAL_LIBDIR}/pnetcdf")
SET(ENV{PNETDCF_DIR} "${COLD_PNETCDF_DIR}")

MESSAGE( STATUS "COLD_PNETCDF_DIR: " ${COLD_PNETCDF_DIR} )

IF(NOT COLD_PNETCDF_LIBPRJ)

    EXTERNALPROJECT_ADD(COLD_PNETCDF_LIBPRJ
      URL http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.7.0.tar.gz
      SOURCE_DIR ${COLD_PNETCDF_DIR}/src
      BINARY_DIR ${COLD_PNETCDF_DIR}/src
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/pnetcdf
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/pnetcdf
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/pnetcdf
      INSTALL_DIR ${COLD_PNETCDF_DIR}/build
      CONFIGURE_COMMAND ./configure --prefix=${COLD_PNETCDF_DIR}/build CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
      BUILD_COMMAND make
      INSTALL_COMMAND make install
      UPDATE_COMMAND "")

ENDIF(NOT COLD_PNETCDF_LIBPRJ)

# get the library
FIND_LIBRARY(COLD_PNETCDF_LIB NAMES pnetcdf PATHS ${COLD_PNETCDF_DIR}/build/lib DOC "PNETCDF libs")

LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_PNETCDF_LIB})
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_PNETCDF_DIR}/build/include)

INCLUDE_DIRECTORIES(${COLD_PNETCDF_DIR}/build/include)




##################################
# ACCFFT
##################################
SET(COLD_ACCFFT_DIR "${COLD_EXTERNAL_LIBDIR}/accfft")
SET(ENV{ACCFFT_DIR} "${COLD_ACCFFT_DIR}")

MESSAGE( STATUS "COLD_ACCFFT_DIR: " ${COLD_ACCFFT_DIR} )

IF(NOT COLD_ACCFFT_LIBPRJ)
    EXTERNALPROJECT_ADD(COLD_ACCFFT_LIBPRJ
      DEPENDS COLD_PNETCDF_LIBPRJ COLD_FFTW_LIBPRJ
      GIT_REPOSITORY git@github.com:amirgholami/accfft.git
      SOURCE_DIR ${COLD_ACCFFT_DIR}/src
      BINARY_DIR ${COLD_ACCFFT_DIR}/build
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      INSTALL_DIR ${COLD_ACCFFT_DIR}/build
      CMAKE_COMMAND cmake
      BUILD_COMMAND make
      INSTALL_COMMAND make install
      UPDATE_COMMAND "")
ENDIF(NOT COLD_ACCFFT_LIBPRJ)

#-DFFTW_USE_STATIC_LIBS=true
#-DBUILD_GPU=false
#-DBUILD_STEPS=false
#-DCXX_FLAGS="-O3"
#-DPNETCDF_DIR=${COLD_PNETCDF_DIR}/build

## find library
ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_PNETCDF_LIB})
ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_FFTW_LIB})
ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_FFTW_THREADS_LIB})
ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_FFTWF_LIB})
ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_FFTWF_THREADS_LIB})
FIND_LIBRARY(COLD_ACCFFT_LIB NAMES accfft PATHS ${COLD_ACCFFT_DIR}/build/lib DOC "ACCFFT libs")
FIND_LIBRARY(COLD_ACCFFTUTILS_LIB NAMES accfft_utils PATHS ${COLD_ACCFFT_DIR}/build/lib DOC "ACCFFT libs")
FIND_PATH(COLD_ACCFFT_INC NAMES "accfft.h" PATHS ${COLD_ACCFFT_DIR}/build/include)


LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_ACCFFT_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_ACCFFTUTILS_LIB})
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_ACCFFT_INC})
INCLUDE_DIRECTORIES(${COLD_ACCFFT_INC})




##################################
# PETSC
##################################
SET(COLD_PETSC_DIR "${COLD_EXTERNAL_LIBDIR}/petsc")
SET(COLD_PETSC_ARCH "cxxopt")
SET(ENV{PETSC_DIR} "${COLD_PETSC_DIR}")
SET(ENV{PETSC_ARCH} "${COLD_PETSC_ARCH}")

SET(PETSC_CONFIG_STD "  " CACHE STRING "PETSc config flags")
SET(PETSC_CONFIG_ADV "  " CACHE STRING "PETSc advanced config flags")

# download and build petsc
IF(NOT COLD_PETSC_LIBPRJ)
    EXTERNALPROJECT_ADD(COLD_PETSC_LIBPRJ
      URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz
      #GIT_REPOSITORY https://bitbucket.org/petsc/petsc
      SOURCE_DIR ${COLD_PETSC_DIR}/${COLD_PETSC_ARC}
      BINARY_DIR ${COLD_PETSC_DIR}/${COLD_PETSC_ARC}
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      INSTALL_DIR ${COLD_PETSC_DIR}/${COLD_PETSC_ARC}
      CONFIGURE_COMMAND ./configure PETSC_DIR=${COLD_PETSC_DIR}
                                    PETSC_ARCH=${COLD_PETSC_ARCH}
                                    --with-clanguage=C++ 
                                    --with-debugging=0 CXXOPTFLAGS='-O3' 
                                    --download-fblaslapack=1 
                                    --with-x=0 --with-64-bit-indices 
                                    --with-c++-support=1 --with-pthread=1 
                                    --with-mpi-dir=/opt/apps/intel14/mvapich2/2.0b
      BUILD_COMMAND make PETSC_DIR=${COLD_PETSC_DIR} PETSC_ARCH=${COLD_PETSC_ARCH} all
      INSTALL_COMMAND "")
ENDIF(NOT COLD_PETSC_LIBPRJ)


## find library
FIND_LIBRARY(COLD_PETSC_LIB NAMES petsc PATHS ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/lib DOC "PETSC libs")

LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_PETSC_LIB})
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_PETSC_DIR}/include)
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/include)

INCLUDE_DIRECTORIES(${COLD_PETSC_DIR}/include)
INCLUDE_DIRECTORIES(${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/include)



##################################
# NIFTICLIB
##################################
SET(COLD_NIFTICLIB_DIR "${COLD_EXTERNAL_LIBDIR}/nifticlib")
SET(ENV{NIFTICLIB_DIR} "${COLD_NIFTICLIB_DIR}")

MESSAGE( STATUS "COLD_NIFTICLIB_DIR: " ${COLD_NIFTICLIB_DIR} )

IF(NOT COLD_NIFTI_LIBPRJ)

    EXTERNALPROJECT_ADD(COLD_NIFTICLIB_LIBPRJ
      URL ${COLD_SOURCE_DIR}/external/nifticlib-2.0.0.tar.gz
      #URL http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz
      SOURCE_DIR ${COLD_NIFTICLIB_DIR}/src
      BINARY_DIR ${COLD_NIFTICLIB_DIR}/src
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/nifticlib
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/nifticlib
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/nifticlib
      INSTALL_DIR ${COLD_NIFTICLIB_DIR}/build
      CMAKE_ARGS -Wno-dev # suppress missing cmake_minimum_required() warning
      CMAKE_CACHE_ARGS "-DBUILD_SHARED_LIBS:BOOL=OFF"
                       "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
                       "-DCMAKE_INSTALL_PREFIX:PATH=${COLD_NIFTICLIB_DIR}/build"
      CMAKE_COMMAND cmake
      BUILD_COMMAND make
      INSTALL_COMMAND make install
)

ENDIF(NOT COLD_NIFTI_LIBPRJ)

# get the library
#FIND_LIBRARY(COLD_NIFTI_LIB NAMES libniftiio PATHS ${COLD_NIFTICLIB_DIR}/build/lib DOC "NIFTICLIB libs")

LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_NIFTI_LIB})
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_NIFTICLIB_DIR}/build/include)

INCLUDE_DIRECTORIES(${COLD_NIFTICLIB_DIR}/build/include)
