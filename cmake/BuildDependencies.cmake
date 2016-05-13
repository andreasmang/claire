INCLUDE(ExternalProject)
# define path for external libraries
SET(COLD_EXTERNAL_LIBDIR "${PROJECT_BINARY_DIR}/external")
MESSAGE( STATUS "cold external library path " ${COLD_EXTERNAL_LIBDIR} )


##################################
# FFTW
##################################
SET(COLD_FFTW_DIR "${COLD_EXTERNAL_LIBDIR}/fftw")
SET(ENV{FFTW_DIR} "${COLD_FFTW_DIR}")
SET(ENV{FFTW_INCLUDES} "${COLD_FFTW_DIR}/build/include")
SET(ENV{FFTW_LIBRARIES} "${COLD_FFTW_DIR}/build/lib")
SET(COLD_FFTW_VERSION 3.3.4)
SET(COLD_FFTW_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")

MESSAGE( STATUS "COLD_FFTW_DIR: " ${COLD_FFTW_DIR} )

#IF(NOT COLD_FFTW_LIBPRJ)

    EXTERNALPROJECT_ADD(COLD_FFTW_LIBPRJ
      #URL ftp://ftp.fftw.org/pub/fftw/fftw-${COLD_FFTW_VERSION}.tar.gz
      URL ${PROJECT_SOURCE_DIR}/external/fftw-3.3.4.tar.gz
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
            --enable-mpi --enable-threads #--enable-shared
            --enable-sse2 --enable-openmp --enable-avx
            CFLAGS=-O3 MAKEINFO=missing
      BUILD_COMMAND make)

    # ACCFFT also needs the float version; we can't make it at once, so
    # we have got to add an additional three steps for the float version
    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIBPRJ lfftw-config
      DEPENDEES install
      COMMAND cd ${COLD_FFTW_DIR}/src && ./configure
            --prefix=${COLD_FFTW_DIR}/build
            --enable-mpi --enable-threads #--enable-shared
            --enable-sse2 --enable-openmp --enable-avx
            --enable-float CFLAGS=-O3 MAKEINFO=missing)

    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIBPRJ lfftw-build
      DEPENDEES lfftw-config
      COMMAND cd ${COLD_FFTW_DIR}/src && make)

    EXTERNALPROJECT_ADD_STEP(COLD_FFTW_LIBPRJ lfftw-install
      DEPENDEES lfftw-build
      COMMAND cd ${COLD_FFTW_DIR}/src && make install)

#ENDIF(NOT COLD_FFTW_LIBPRJ)

## find library
#FIND_LIBRARY(COLD_FFTW_LIB NAMES fftw3 PATHS ${COLD_FFTW_DIR}/build/lib DOC "FFTW libs")
#FIND_LIBRARY(COLD_FFTW_THREADS_LIB NAMES fftw3_threads PATHS ${COLD_FFTW_DIR}/build/lib DOC "FFTW threads libs")
#FIND_LIBRARY(COLD_FFTWF_LIB NAMES fftw3f PATHS ${COLD_FFTW_DIR}/build/lib DOC "FFTWF libs")
#FIND_LIBRARY(COLD_FFTWF_THREADS_LIB NAMES fftw3f_threads PATHS ${COLD_FFTW_DIR}/build/lib DOC "FFTWF threads libs")

# add library path to serach
LINK_DIRECTORIES(${COLD_FFTW_DIR}/build/lib)

#SET(COLD_FFTW_LIB ${COLD_FFTW_DIR}/build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3${CMAKE_SHARED_LIBRARY_SUFFIX})
#SET(COLD_FFTW_THREADS_LIB ${COLD_FFTW_DIR}/build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3_threads${CMAKE_SHARED_LIBRARY_SUFFIX})
#SET(COLD_FFTWF_LIB ${COLD_FFTW_DIR}/build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f${CMAKE_SHARED_LIBRARY_SUFFIX})
#SET(COLD_FFTWF_THREADS_LIB ${COLD_FFTW_DIR}/build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f_threads${CMAKE_SHARED_LIBRARY_SUFFIX})


#SET(COLD_FFTW_LIB ${COLD_FFTW_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_FFTW_THREADS_LIB ${COLD_FFTW_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_threads${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_FFTWF_LIB ${COLD_FFTW_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_FFTWF_THREADS_LIB ${COLD_FFTW_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX})

SET(COLD_FFTW_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
SET(COLD_FFTW_THREADS_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_threads${CMAKE_STATIC_LIBRARY_SUFFIX})
SET(COLD_FFTWF_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f${CMAKE_STATIC_LIBRARY_SUFFIX})
SET(COLD_FFTWF_THREADS_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX})

#ADD_LIBRARY(COLD_FFTW_LIB STATIC IMPORTED)
#SET_TARGET_PROPERTIES(COLD_FFTW_LIB PROPERTIES IMPORTED_LOCATION ${COLD_FFTW_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})

#ADD_LIBRARY(COLD_FFTW_THREADS_LIB STATIC IMPORTED)
#SET_TARGET_PROPERTIES(COLD_FFTW_THREADS_LIB PROPERTIES IMPORTED_LOCATION ${COLD_FFTW_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_threads${CMAKE_STATIC_LIBRARY_SUFFIX})


#ADD_LIBRARY(COLD_FFTWF_LIB STATIC IMPORTED)
#SET_TARGET_PROPERTIES(COLD_FFTWF_LIB PROPERTIES IMPORTED_LOCATION ${COLD_FFTW_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f${CMAKE_STATIC_LIBRARY_SUFFIX})

#ADD_LIBRARY(COLD_FFTWF_THREADS_LIB STATIC IMPORTED)
#SET_TARGET_PROPERTIES(COLD_FFTWF_THREADS_LIB PROPERTIES IMPORTED_LOCATION ${COLD_FFTW_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX})


LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_FFTW_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_FFTW_THREADS_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_FFTWF_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_FFTWF_THREADS_LIB})

LIST(APPEND COLD_EXTERNAL_INCS ${COLD_FFTW_DIR}/build/include)

INCLUDE_DIRECTORIES(${COLD_FFTW_DIR}/build/include)




##################################
# ACCFFT
##################################
SET(COLD_ACCFFT_DIR "${COLD_EXTERNAL_LIBDIR}/accfft")
SET(ENV{ACCFFT_DIR} "${COLD_ACCFFT_DIR}")

MESSAGE( STATUS "COLD_ACCFFT_DIR: " ${COLD_ACCFFT_DIR} )

#IF(NOT COLD_ACCFFT_LIBPRJ)
    EXTERNALPROJECT_ADD(COLD_ACCFFT_LIBPRJ
      DEPENDS COLD_FFTW_LIBPRJ
      #GIT_REPOSITORY git@github.com:amirgholami/accfft.git
      URL ${PROJECT_SOURCE_DIR}/external/accfft.tar.gz
      SOURCE_DIR ${COLD_ACCFFT_DIR}/src
      BINARY_DIR ${COLD_ACCFFT_DIR}/build
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/accfft
      INSTALL_DIR ${COLD_ACCFFT_DIR}/build
      CMAKE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${COLD_ACCFFT_DIR}/build
                          -DFFTW_ROOT=${COLD_FFTW_DIR}/build -DFFTW_USE_STATIC_LIBS=true
      BUILD_COMMAND make
      INSTALL_COMMAND make install
      UPDATE_COMMAND "")
#ENDIF(NOT COLD_ACCFFT_LIBPRJ)

#ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_PNETCDF_LIB})
#ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_FFTW_LIB})
#ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_FFTW_THREADS_LIB})
#ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_FFTWF_LIB})
#ADD_DEPENDENCIES(COLD_ACCFFT_LIBPRJ ${COLD_FFTWF_THREADS_LIB})
#TARGET_LINK_LIBRARIES(COLD_ACCFFT_LIBPRJ ${COLD_PNETCDF_LIB})
#COLD_PNETCDF_LIBPRJ 

#-DPNETCDF_DIR=${COLD_PNETCDF_DIR}/build -DBUILD_SHARED=false
#-DFFTW_USE_STATIC_LIBS=true
#-DBUILD_GPU=false
#-DBUILD_STEPS=false
#


#FIND_LIBRARY(COLD_ACCFFT_LIB NAMES accfft PATHS ${COLD_ACCFFT_DIR}/build/lib DOC "ACCFFT libs")
#FIND_LIBRARY(COLD_ACCFFTUTILS_LIB NAMES accfft_utils PATHS ${COLD_ACCFFT_DIR}/build/lib DOC "ACCFFT libs")
#FIND_PATH(COLD_ACCFFT_INC NAMES "accfft.h" PATHS ${COLD_ACCFFT_DIR}/build/include)

#SET(COLD_ACCFFT_LIB ${COLD_ACCFFT_DIR}/build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}accfft${CMAKE_SHARED_LIBRARY_SUFFIX})
#SET(COLD_ACCFFTUTILS_LIB ${COLD_ACCFFT_DIR}/build/lib/${CMAKE_SHARED_LIBRARY_PREFIX}accfft_utils${CMAKE_SHARED_LIBRARY_SUFFIX})

#SET(COLD_ACCFFT_LIB ${COLD_ACCFFT_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}accfft${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_ACCFFTUTILS_LIB ${COLD_ACCFFT_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}accfft_utils${CMAKE_STATIC_LIBRARY_SUFFIX})
SET(COLD_ACCFFT_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}accfft${CMAKE_STATIC_LIBRARY_SUFFIX})
SET(COLD_ACCFFTUTILS_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}accfft_utils${CMAKE_STATIC_LIBRARY_SUFFIX})
TARGET_LINK_LIBRARIES( COLD_ACCFFT_LIBPRJ COLD_FFTW_LIB )
TARGET_LINK_LIBRARIES( COLD_ACCFFT_LIBPRJ COLD_FFTW_THREADS_LIB )
TARGET_LINK_LIBRARIES( COLD_ACCFFT_LIBPRJ COLD_FFTWF_LIB )
TARGET_LINK_LIBRARIES( COLD_ACCFFT_LIBPRJ COLD_FFTWF_THREADS_LIB )

#ADD_LIBRARY(COLD_ACCFFT_LIB STATIC IMPORTED)
#SET_TARGET_PROPERTIES(COLD_ACCFFT_LIB PROPERTIES IMPORTED_LOCATION ${COLD_ACCFFT_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}accfft${CMAKE_STATIC_LIBRARY_SUFFIX})

#ADD_LIBRARY(COLD_ACCFFTUTILS_LIB STATIC IMPORTED)
#SET_TARGET_PROPERTIES(COLD_ACCFFTUTILS_LIB PROPERTIES IMPORTED_LOCATION ${COLD_ACCFFT_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}accfft${CMAKE_STATIC_LIBRARY_SUFFIX})


LINK_DIRECTORIES(${COLD_ACCFFT_DIR}/build/lib)
INCLUDE_DIRECTORIES(${COLD_ACCFFT_DIR}/build/include)
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_ACCFFT_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_ACCFFTUTILS_LIB})
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_ACCFFT_DIR}/build/include)




##################################
# PETSC
##################################
SET(COLD_PETSC_DIR "${COLD_EXTERNAL_LIBDIR}/petsc")
SET(COLD_PETSC_ARCH "cxx_opt")
SET(ENV{PETSC_DIR} "${COLD_PETSC_DIR}")
SET(ENV{PETSC_ARCH} "${COLD_PETSC_ARCH}")

# download and build petsc
#IF(NOT COLD_PETSC_LIBPRJ)
    EXTERNALPROJECT_ADD(COLD_PETSC_LIBPRJ
      #URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz
      #GIT_REPOSITORY https://bitbucket.org/petsc/petsc
      URL ${PROJECT_SOURCE_DIR}/external/petsc-lite-3.7.0.tar.gz
      URL_MD5 f04bf6163e735b164159c71ab50af232
      SOURCE_DIR ${COLD_PETSC_DIR}/${COLD_PETSC_ARC}
      BINARY_DIR ${COLD_PETSC_DIR}/${COLD_PETSC_ARC}
      TMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      STAMP_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      DOWNLOAD_DIR ${COLD_EXTERNAL_LIBDIR}/tmp/petsc
      INSTALL_DIR ${COLD_PETSC_DIR}/build
      CONFIGURE_COMMAND ./configure PETSC_DIR=${COLD_PETSC_DIR}
                                    PETSC_ARCH=${COLD_PETSC_ARCH}
                                    --prefix ${COLD_PETSC_DIR}/build
                                    --with-clanguage=C++
                                    --with-debugging=0 CXXOPTFLAGS=${CMAKE_CXX_FLAGS} COPTFLAGS=${CMAKE_C_FLAGS}
                                    --with-shared=0
                                    --with-blas-lapack-dir=$ENV{TACC_MKL_LIB}
                                    --with-x=0 --with-64-bit-indices
                                    --with-c++-support=1 --with-pthread=1
                                    --with-cc=${MPI_C_COMPILER} --with-cxx=${MPI_CXX_COMPILER} --with-fc=0
      BUILD_COMMAND make PETSC_DIR=${COLD_PETSC_DIR} PETSC_ARCH=${COLD_PETSC_ARCH} all
      INSTALL_COMMAND make install)
#ENDIF(NOT COLD_PETSC_LIBPRJ)

#--download-f2cblaslapack 
#--with-mpi-dir=/opt/apps/intel14/mvapich2/2.0b
## find library
#FIND_LIBRARY(COLD_PETSC_LIB NAMES petsc PATHS ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/lib DOC "PETSC libs")
#SET(COLD_PETSC_LIB ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}petsc${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_LAPACK_LIB ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}f2clapack${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_BLAS_LIB ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}f2cblas${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_PETSC_LIB ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}petsc${CMAKE_SHARED_LIBRARY_SUFFIX})
#SET(COLD_PETSC_LIB ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}petsc${CMAKE_STATIC_LIBRARY_SUFFIX})

#SET(COLD_LAPACK_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}f2clapack${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_BLAS_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}f2cblas${CMAKE_STATIC_LIBRARY_SUFFIX})
SET(COLD_PETSC_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}petsc${CMAKE_STATIC_LIBRARY_SUFFIX})

#LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_BLAS_LIB})
#LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_LAPACK_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_PETSC_LIB})

LIST(APPEND COLD_EXTERNAL_INCS ${COLD_PETSC_DIR}/include)
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/include)



#TARGET_INCLUDE_DIRECTORIES (COLD_PETSC_LIB PUBLIC ${COLD_PETSC_DIR}/include)
#TARGET_INCLUDE_DIRECTORIES (COLD_PETSC_LIB PUBLIC ${COLD_PETSC_DIR}/${COLD_PETSC_ARCH}/include)

#SET_PROPERTY(TARGET COLD_PETSC_LIBPRJ PROPERTY IMPORTED_LOCATION ${COLD_PETSC_LIB})
#SET_PROPERTY(TARGET COLD_PETSC_LIBPRJ PROPERTY IMPORTED_LOCATION ${COLD_LAPACK_LIB})
#SET_PROPERTY(TARGET COLD_PETSC_LIBPRJ PROPERTY IMPORTED_LOCATION ${COLD_BLAS_LIB})

# add include path
INCLUDE_DIRECTORIES(SYSTEM ${COLD_PETSC_DIR}/build/include/)
INCLUDE_DIRECTORIES(SYSTEM $ENV{TACC_MKL_DIR}/include)

# add library path
LINK_DIRECTORIES(${COLD_PETSC_DIR}/build/lib)
LINK_DIRECTORIES(SYSTEM $ENV{TACC_MKL_DIR}/lib/intel64)



##################################
# NIFTICLIB
##################################
SET(COLD_NIFTICLIB_DIR "${COLD_EXTERNAL_LIBDIR}/nifticlib")
SET(ENV{NIFTICLIB_DIR} "${COLD_NIFTICLIB_DIR}")

MESSAGE( STATUS "COLD_NIFTICLIB_DIR: " ${COLD_NIFTICLIB_DIR} )

#IF(NOT COLD_NIFTI_LIBPRJ)

    EXTERNALPROJECT_ADD(niftiio
      #URL http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz
      URL ${PROJECT_SOURCE_DIR}/external/nifticlib-2.0.0.tar.gz
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

#ENDIF(NOT COLD_NIFTI_LIBPRJ)

#ExternalProject_Get_Property(COLD_NIFTI_LIBPRJ source_dir)
#ExternalProject_Get_Property(COLD_NIFTI_LIBPRJ binary_dir)
#MESSAGE(STATUS "source directory" "${source_dir}")
#MESSAGE(STATUS "binary directory" "${binary_dir}")


#SET(COLD_NIFTIIO_LIB ${COLD_NIFTICLIB_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}niftiio${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_NIFTICDF_LIB ${COLD_NIFTICLIB_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}nifticdf${CMAKE_STATIC_LIBRARY_SUFFIX})
#SET(COLD_NIFTIZNZ_LIB ${COLD_NIFTICLIB_DIR}/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}znz${CMAKE_STATIC_LIBRARY_SUFFIX})

SET(COLD_NIFTIIO_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}niftiio${CMAKE_STATIC_LIBRARY_SUFFIX})
SET(COLD_NIFTICDF_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}nifticdf${CMAKE_STATIC_LIBRARY_SUFFIX})
SET(COLD_NIFTIZNZ_LIB ${CMAKE_STATIC_LIBRARY_PREFIX}znz${CMAKE_STATIC_LIBRARY_SUFFIX})



LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_NIFTICDF_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_NIFTIIO_LIB})
LIST(APPEND COLD_EXTERNAL_LIBS ${COLD_NIFTIZNZ_LIB})
LIST(APPEND COLD_EXTERNAL_INCS ${COLD_NIFTICLIB_DIR}/build/include/nifti)
LIST(APPEND COLD_EXTERNAL_LIBS $ENV{TACC_MKL_DIR}/lib/intel64)

LINK_DIRECTORIES(${COLD_NIFTICLIB_DIR}/build/lib)
MESSAGE(STATUS "nifti libaries ${COLD_NIFTICLIB_DIR}/build/lib" )
INCLUDE_DIRECTORIES(${COLD_NIFTICLIB_DIR}/build/include/nifti)

