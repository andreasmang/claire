INCLUDE(ExternalProject)




# define path for external libraries
SET(COLD_EXTERNAL_LIBDIR "${CMAKE_CURRENT_BINARY_DIR}/external")
SET(COLD_EXTERNAL_LIBS "")
SET(COLD_EXTERNAL_INCLUDES "")



##################################
# FFTW
##################################
SET(COLD_FFTW_DIR "${COLD_EXTERNAL_LIBDIR}/fftw")
SET(ENV{FFTW_DIR} "${COLD_FFTW_DIR}")
SET(COLD_FFTW_VERSION 3.3.4)
SET(COLD_FFTW_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")

MESSAGE( STATUS "COLD_FFTW_DIR: " ${COLD_FFTW_DIR} )

IF(NOT COLD_FFTW_LIB)

    EXTERNALPROJECT_ADD(COLD_FFTW_LIB
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
            CFLAGS=-O3 MAKEINFO=missing)

    # ACCFFT also needs the float version; we can't make it at once, so
    # we have got to add an additional three steps for the float version
    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIB lfftw-config
      DEPENDEES install
      COMMAND cd ${COLD_FFTW_DIR}/src && ./configure
            --prefix=${COLD_FFTW_DIR}/build
            --enable-mpi --enable-threads --enable-shared
            --enable-sse2 --enable-openmp --enable-avx
            --enable-float CFLAGS=-O3 MAKEINFO=missing)

    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIB lfftw-build
      DEPENDEES lfftw-config
      COMMAND cd ${COLD_FFTW_DIR}/src && make)

    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIB lfftw-install
      DEPENDEES lfftw-build
      COMMAND cd ${COLD_FFTW_DIR}/src && make install)

ENDIF(NOT COLD_FFTW_LIB)

LIST(APPEND COLD_EXTERNAL_LIBS COLD_FFTW_LIB) 
LIST(APPEND COLD_EXTERNAL_INCLUDES ${COLD_FFTW_DIR}/build/include)


##################################
# PNETCDF
##################################
SET(COLD_PNETCDF_DIR "${COLD_EXTERNAL_LIBDIR}/pnetcdf")
SET(ENV{PNETDCF_DIR} "${COLD_PNETCDF_DIR}")

MESSAGE( STATUS "COLD_PNETCDF_DIR: " ${COLD_PNETCDF_DIR} )

IF(NOT COLD_PNETCDF_LIB)

    EXTERNALPROJECT_ADD(COLD_PNETCDF_LIB
      URL http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.7.0.tar.gz
      SOURCE_DIR ${COLD_PNETCDF_DIR}/src
      BINARY_DIR ${COLD_PNETCDF_DIR}/src
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/pnetcdf
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/pnetcdf
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/pnetcdf
      INSTALL_DIR ${COLD_PNETCDF_DIR}/build
      CONFIGURE_COMMAND ./configure --prefix=${COLD_PNETCDF_DIR}/build
            CFLAGS=-O3 MAKEINFO=missing
      UPDATE_COMMAND "")

ENDIF(NOT COLD_PNETCDF_LIB)

LIST(APPEND COLD_EXTERNAL_LIBS COLD_PNETCDF_LIB) 
LIST(APPEND COLD_EXTERNAL_INCLUDES ${COLD_PNETCDF_DIR}/build/include)



##################################
# ACCFFT
##################################
SET(COLD_ACCFFT_DIR "${COLD_EXTERNAL_LIBDIR}/accfft")
SET(ENV{ACCFFT_DIR} "${COLD_ACCFFT_DIR}")

MESSAGE( STATUS "COLD_ACCFFT_DIR: " ${COLD_ACCFFT_DIR} )

IF(NOT COLD_ACCFFT_LIB)
    EXTERNALPROJECT_ADD(COLD_ACCFFT_LIB
      DEPENDS COLD_PNETCDF_LIB COLD_FFTW_LIB
      GIT_REPOSITORY git@github.com:amirgholami/accfft.git
      SOURCE_DIR ${COLD_ACCFFT_DIR}/src
      BINARY_DIR ${COLD_ACCFFT_DIR}/build
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      INSTALL_DIR ${COLD_ACCFFT_DIR}/build
      CMAKE_COMMAND cmake -DFFTW_ROOT=${COLD_FFTW_DIR}/build
                    -DFFTW_USE_STATIC_LIBS=false
                    -DBUILD_GPU=false
                    -DBUILD_STEPS=false
                    -DCXX_FLAGS="-O3"
                    -DPNETCDF_DIR=${COLD_PNETCDF_DIR}/build
                    -DBUILD_SHARED=false
      UPDATE_COMMAND "")
ENDIF(NOT COLD_ACCFFT_LIB)

LIST(APPEND COLD_EXTERNAL_LIBS COLD_ACCFFT_LIB)
LIST(APPEND COLD_EXTERNAL_INCLUDES ${COLD_ACCFFT_DIR}/build/include)



##################################
### PETSC
##################################
SET(COLD_PETSC_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/petsc")
SET(COLD_PETSC_ARCH "cxxopt")
SET(ENV{PETSC_DIR} "${COLD_PETSC_DIR}")
SET(ENV{PETSC_ARCH} "${COLD_PETSC_ARCH}")

SET(PETSC_CONFIG_STD "  " CACHE STRING "PETSc config flags")
SET(PETSC_CONFIG_ADV "  " CACHE STRING "PETSc advanced config flags")

# download and build petsc
IF(NOT COLD_PETSC_LIB)
    EXTERNALPROJECT_ADD(COLD_PETSC_LIB
      URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz
      #GIT_REPOSITORY https://bitbucket.org/petsc/petsc
      BUILD_IN_SOURCE 1
      SOURCE_DIR ${COLD_PETSC_DIR}/${PETSC_ARC}
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      CONFIGURE_COMMAND ./configure PETSC_DIR=${COLD_PETSC_DIR}
                                    PETSC_ARCH=${COLD_PETSC_ARCH}
                                    --with-clanguage=C++ 
                                    --with-debugging=0 CXXOPTFLAGS='-O3' 
                                    --download-fblaslapack=1 
                                    --with-x=0 --with-64-bit-indices 
                                    --with-c++-support=1 --with-pthread=1 
                                    --with-mpi-dir=/opt/apps/intel14/mvapich2/2.0b
      BUILD_COMMAND  make PETSC_DIR=${COLD_PETSC_DIR} PETSC_ARCH=${COLD_PETSC_ARCH} all
      INSTALL_COMMAND "")
ENDIF(NOT COLD_PETSC_LIB)

LIST(APPEND COLD_EXTERNAL_LIBS COLD_PETSC_LIB)
LIST(APPEND COLD_EXTERNAL_INCLUDES ${COLD_PETSC_DIR}/include)
