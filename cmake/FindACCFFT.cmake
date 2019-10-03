# - Find the ACCFFT library
#
# Usage:
#   find_package(ACCFFT [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   ACCFFT_FOUND               ... true if accfft is found on the system
#   ACCFFT_LIBRARIES           ... full path to accfft library
#   ACCFFT_INCLUDES            ... accfft include directory
#
# The following variables will be checked by the function
#   ACCFFT_USE_STATIC_LIBS    ... if true, only static libraries are found
#   ACCFFT_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   ACCFFT_LIBRARY            ... accfft library to use
#   ACCFFT_INCLUDE_DIR        ... accfft include directory
#

#If environment variable ACCFFT_ROOT_DIR is specified, it has same effect as ACCFFT_ROOT
if( NOT ACCFFT_ROOT AND DEFINED ENV{ACCFFT_DIR} )
  set( ACCFFT_ROOT $ENV{ACCFFT_DIR} )
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT ACCFFT_ROOT )
  pkg_check_modules( PKG_ACCFFT QUIET "accfft" )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

#if( ${ACCFFT_USE_STATIC_LIBS} )
  set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
#else()
#  set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )
#endif()

if( ACCFFT_ROOT )

  find_library(
    ACCFFT_LIB
    NAMES "accfft"
    PATHS ${ACCFFT_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )

  find_library(
    ACCFFT_UTILS_LIB
    NAMES "accfft_utils"
    PATHS ${ACCFFT_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )


  #find includes
  find_path(
    ACCFFT_INCLUDES
    NAMES "accfft.h"
    PATHS ${ACCFFT_ROOT}
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
  )

else()

  find_library(
    ACCFFT_LIB
    NAMES "accfft"
    PATHS ${PKG_ACCFFT_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
  )
  find_library(
    ACCFFT_UTILS_LIB
    NAMES "accfft_utils"
    PATHS ${PKG_ACCFFT_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
  )

  find_path(
    ACCFFT_INCLUDES
    NAMES "accfft.h"
    PATHS ${PKG_ACCFFT_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
  )

endif( ACCFFT_ROOT )

set(ACCFFT_LIBRARIES ${ACCFFT_LIB} ${ACCFFT_UTILS_LIB})

if(ACCFFT_UTILS_LIB)
  set(ACCFFT_LIBRARIES ${ACCFFT_LIBRARIES} ${ACCFFT_UTILS_LIB})
endif()

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ACCFFT DEFAULT_MSG
                                  ACCFFT_INCLUDES ACCFFT_LIBRARIES)

mark_as_advanced(ACCFFT_INCLUDES ACCFFT_LIBRARIES ACCFFT_UTILS_LIB ACCFFT_LIB)

