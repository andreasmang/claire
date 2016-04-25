INCLUDE(ExternalProject)

# define path for external libraries
SET(COLD_EXTERNAL_LIBS "${CMAKE_CURRENT_BINARY_DIR}/external")



##################################
# FFTW
##################################
SET(COLD_FFTW_DIR "${COLD_EXTERNAL_LIBS}/fftw")
SET(ENV{FFTW_DIR} "${COLD_FFTW_DIR}")
SET(COLD_FFTW_VERSION 3.3.4)
SET(COLD_FFTW_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")

MESSAGE( STATUS "COLD_FFTW_DIR: " ${COLD_FFTW_DIR} )

EXTERNALPROJECT_ADD(FFTW
  URL ftp://ftp.fftw.org/pub/fftw/fftw-${COLD_FFTW_VERSION}.tar.gz
  URL_MD5 ${COLD_FFTW_MD5}
  SOURCE_DIR ${COLD_FFTW_DIR}/src
  BINARY_DIR ${COLD_FFTW_DIR}/src
  TMP_DIR ${COLD_EXTERNAL_LIBS}/tmp/fftw
  STAMP_DIR ${COLD_EXTERNAL_LIBS}/tmp/fftw
  DOWNLOAD_DIR ${COLD_EXTERNAL_LIBS}/tmp/fftw
  INSTALL_DIR ${COLD_FFTW_DIR}/build
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ./configure
        --prefix=${COLD_FFTW_DIR}/build
        --enable-mpi --enable-threads --enable-shared
        --enable-sse2 --enable-openmp --enable-avx
        CFLAGS=-O3 MAKEINFO=missing)


##################################
# ACCFFT
##################################
SET(COLD_ACCFFT_DIR "${COLD_EXTERNAL_LIBS}/accfft")
SET(ENV{ACCFFT_DIR} "${COLD_ACCFFT_DIR}")

MESSAGE( STATUS "COLD_ACCFFT_DIR: " ${COLD_ACCFFT_DIR} )

EXTERNALPROJECT_ADD(ACCFFT
  DEPENDS FFTW
  GIT_REPOSITORY git@github.com:amirgholami/accfft.git #https://github.com/amirgholami/accfft
  SOURCE_DIR ${COLD_ACCFFT_DIR}/src
  BINARY_DIR ${COLD_ACCFFT_DIR}/build
  TMP_DIR ${COLD_EXTERNAL_LIBS}/tmp/accfft
  STAMP_DIR ${COLD_EXTERNAL_LIBS}/tmp/accfft
  DOWNLOAD_DIR ${COLD_EXTERNAL_LIBS}/tmp/accfft
  INSTALL_DIR ${COLD_ACCFFT_DIR}/build
  CMAKE_COMMAND cmake -DFFTW_ROOT=${COLD_FFTW_DIR}/build
                -DFFTW_USE_STATIC_LIBS=false
                -DBUILD_GPU=false
                -DBUILD_STEPS=false
                -DCXX_FLAGS="-O3"
                -DBUILD_SHARED=false
  UPDATE_COMMAND "")



##################################
### PETSC
##################################
SET(COLD_PETSC_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/petsc")
SET(COLD_PETSC_ARCH "cxxopt")
SET(ENV{PETSC_DIR} "${COLD_PETSC_DIR}")
SET(ENV{PETSC_ARCH} "${COLD_PETSC_ARCH}")

SET(PETSC_CONFIG_STD " --with-clanguage=C++ --with-debugging=0 CXXOPTFLAGS='-O3' --download-fblaslapack=1 --with-x=0 --with-64-bit-indices --with-c++-support=1 --with-pthread=1 " CACHE STRING "PETSc config flags")
SET(PETSC_CONFIG_ADV " --with-mpi-dir=/opt/apps/intel14/mvapich2/2.0b " CACHE STRING "PETSc advanced config flags")

# download and build petsc
#EXTERNALPROJECT_ADD(petsclib
#  URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz
  #GIT_REPOSITORY https://bitbucket.org/petsc/petsc
#  BUILD_IN_SOURCE 1
#  SOURCE_DIR ${EXTERNAL_LIBS}/petsc
#  CONFIGURE_COMMAND ./configure PETSC_DIR=${COLD_PETSC_DIR} PETSC_ARCH=${COLD_PETSC_ARCH} ${PETSC_CONFIG_STD} ${PETSC_CONFIG_ADV}
#  BUILD_COMMAND make PETSC_DIR=${COLD_PETSC_DIR} PETSC_ARCH=${COLD_PETSC_ARCH} all
#  INSTALL_COMMAND "")
