# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##############################################################################
# @file  FindNIFTI.cmake
# @brief Find nifticlib package.
#
# @par Input variables:
# <table border="0">
#   <tr>
#     @tp @b NIFTI_DIR @endtp
#     <td>The nifticlib package files are searched under the specified root
#         directory. If they are not found there, the default search paths
#         are considered. This variable can also be set as environment variable.</td>
#   </tr>
#   <tr>
#     @tp @b NIFTICLIB_DIR @endtp
#     <td>Alternative environment variable for @p NIFTI_DIR.</td>
#   </tr>
#   <tr>
#     @tp @b NIFTI_USE_STATIC_LIB @endtp
#     <td>Forces this module to search for the static library. Otherwise,
#         the shared library is preferred.</td>
#   </tr>
# </table>
#
# @par Output variables:
# <table border="0">
#   <tr>
#     @tp  @b NIFTI_FOUND @endtp
#     <td>Whether the nifticlib package was found and the following CMake
#         variables are valid.</td>
#   </tr>
#   <tr>
#     @tp @b NIFTI_INCLUDE_DIR @endtp
#     <td>Cached include directory/ies.</td>
#   </tr>
#   <tr>
#     @tp @b NIFTI_INCLUDE_DIRS @endtp
#     <td>Alias for @p NIFTI_INCLUDE_DIR (not cached).</td>
#   </tr>
#   <tr>
#     @tp @b NIFTI_INCLUDES @endtp
#     <td>Alias for @p NIFTI_INCLUDE_DIR (not cached).</td>
#   </tr>
#   <tr>
#     @tp @b NIFTI_LIBRARY @endtp
#     <td>Path of @c niftiio library.</td>
#   </tr>
#   <tr>
#     @tp @b NIFTI_LIB @endtp
#     <td>Alias for @p NIFTI_LIBRARY (not cached).</td>
#   </tr>
#   <tr>
#     @tp @b NIFTI_LIBRARIES @endtp
#     <td>Path of @c niftiio library and prerequisite libraries.</td>
#   </tr>
# </table>
#
# @par Imported targets:
# <table border="0">
#   <tr>
#     @tp @b niftiio @endtp
#     <td>The library target of the @c nifticlib library.</td>
#   </tr>
# </table>
#
# @ingroup CMakeFindModules
##############################################################################

# ----------------------------------------------------------------------------
# initialize search
if (NOT NIFTI_DIR)
  if (NOT "$ENV{NIFTICLIB_DIR}" STREQUAL "")
    set (NIFTI_DIR "$ENV{NIFTICLIB_DIR}" CACHE PATH "Installation prefix for NIFTI." FORCE)
  else ()
    set (NIFTI_DIR "$ENV{NIFTI_DIR}" CACHE PATH "Installation prefix for NIFTI." FORCE)
  endif ()
endif ()

set (NIFTI_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if (NIFTI_USE_STATIC_LIB)
  if (WIN32)
    set (CMAKE_FIND_LIBRARY_SUFFIXES .lib)
  else ()
    set (CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()
else ()
  if (WIN32)
    set (CMAKE_FIND_LIBRARY_SUFFIXES .dll .lib)
  elseif(APPLE)
    set (CMAKE_FIND_LIBRARY_SUFFIXES .dylib .a)
  else ()
    set (CMAKE_FIND_LIBRARY_SUFFIXES .so .a)
  endif()
endif ()

# ----------------------------------------------------------------------------
# find paths/files
if (NIFTI_DIR)

  find_path (
    NIFTI_INCLUDE_DIR
      NAMES         nifti/nifti1_io.h
      HINTS         ${NIFTI_DIR}
      PATH_SUFFIXES "include"
      DOC           "path to directory containing nifti1_io.h file."
      NO_DEFAULT_PATH
  )

  find_library (
    NIFTI_LIBRARY
      NAMES         niftiio
      HINTS         ${NIFTI_DIR}
      PATH_SUFFIXES lib
      DOC           "Path of niftiio library"
      NO_DEFAULT_PATH
  )
  
  find_library (
    NIFTI_CDF_LIBRARY
      NAMES         nifticdf
      HINTS         ${NIFTI_DIR}
      PATH_SUFFIXES lib
      DOC           "Path of nifticdf library"
      NO_DEFAULT_PATH
  )

  find_library (
    NIFTI_znz_LIBRARY
      NAMES znz
      HINTS ENV LD_LIBRARY_PATH ${NIFTI_DIR}
      PATH_SUFFIXES lib
      DOC   "Path of znz library"
  )

else ()

  find_path (
    NIFTI_INCLUDE_DIR
      NAMES         nifti/nifti1_io.h
      HINTS         ENV C_INCLUDE_PATH ENV CXX_INCLUDE_PATH
      DOC           "path to directory containing nifti1_io.h file."
  )

  find_library (
    NIFTI_LIBRARY
      NAMES niftiio
      HINTS ENV LD_LIBRARY_PATH
      DOC   "Path of niftiio library"
  )
  
  find_library (
    NIFTI_CDF_LIBRARY
      NAMES nifticdf
      HINTS ENV LD_LIBRARY_PATH
      DOC   "Path of nifticdf library"
  )

  find_library (
    NIFTI_znz_LIBRARY
      NAMES znz
      HINTS ENV LD_LIBRARY_PATH
      DOC   "Path of znz library"
  )

endif ()

#set( NIFTI_LIBRARY ${NIFTI_CDF_LIBRARY} ${NIFTI_LIBRARY})

mark_as_advanced (NIFTI_INCLUDE_DIR)
mark_as_advanced (NIFTI_LIBRARY)
mark_as_advanced (NIFTI_CDF_LIBRARY)
mark_as_advanced (NIFTI_znz_LIBRARY)

# ----------------------------------------------------------------------------
# prerequisites
if (NIFTI_USE_STATIC_LIB OR NIFTI_znz_LIBRARY MATCHES "\\.a$")
#  find_package (ZLIB REQUIRED)
endif ()

set(ZLIB_LIBRARIES "")

set (NIFTI_LIBRARIES "${ZLIB_LIBRARIES}")
if (NIFTI_CDF_LIBRARY)
    list (APPEND NIFTI_LIBRARIES "${NIFTI_CDF_LIBRARY}")
endif ()
if (NIFTI_LIBRARY)
  list (APPEND NIFTI_LIBRARIES "${NIFTI_LIBRARY}")
endif ()
if (NIFTI_znz_LIBRARY)
  list (APPEND NIFTI_LIBRARIES "${NIFTI_znz_LIBRARY}")
endif ()

# ----------------------------------------------------------------------------
# import targets
if (NIFTI_znz_LIBRARY)
  if (NIFTI_USE_STATIC_LIB OR NIFTI_znz_LIBRARY MATCHES "\\.a$")
    add_library (niftiznz STATIC IMPORTED)
  else ()
    add_library (niftiznz SHARED IMPORTED)
  endif ()
  set_target_properties (
    niftiznz
    PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION                 "${NIFTI_znz_LIBRARY}"
      IMPORTED_LINK_INTERFACE_LIBRARIES "${ZLIB_LIBRARIES}"
  )
endif ()

if (NIFTI_LIBRARY)
  if (NIFTI_USE_STATIC_LIB OR NIFTI_LIBRARY MATCHES "\\.a$")
    add_library (niftiio STATIC IMPORTED)
  else ()
    add_library (niftiio SHARED IMPORTED)
  endif ()
  set_target_properties (
    niftiio
    PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION                 "${NIFTI_LIBRARY}"
  )
  if (TARGET niftiznz)
    set_target_properties (niftiio PROPERTIES IMPORTED_LINK_INTERFACE_LIBRARIES niftiznz)
  endif ()
endif ()

# ----------------------------------------------------------------------------
# aliases / backwards compatibility
if (NIFTI_INCLUDE_DIR)
  set (NIFTI_INCLUDE_DIRS "${NIFTI_INCLUDE_DIR}")
  if (NOT NIFTI_INCLUDE_DIR MATCHES "/nifti/?$")
    list (APPEND NIFTI_INCLUDE_DIRS "${NIFTI_INCLUDE_DIR}/nifti")
  endif ()
  set (NIFTI_INCLUDES "${NIFTI_INCLUDE_DIRS}")
endif ()

if (NIFTI_LIBRARY)
  set (NIFTI_LIB "${NIFTI_LIBRARY}")
endif ()

# ----------------------------------------------------------------------------
# reset CMake variables
set (CMAKE_FIND_LIBRARY_SUFFIXES ${NIFTI_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

# ----------------------------------------------------------------------------
# handle the QUIETLY and REQUIRED arguments and set *_FOUND to TRUE
# if all listed variables are found or TRUE
include (FindPackageHandleStandardArgs)

find_package_handle_standard_args (
  NIFTI
  REQUIRED_VARS
    NIFTI_INCLUDE_DIR
    NIFTI_LIBRARY
    NIFTI_znz_LIBRARY
)

set (NIFTI_FOUND ${NIFTICLIB_FOUND})

# ----------------------------------------------------------------------------
# set NIFTI_DIR
if (NOT NIFTI_DIR AND NIFTI_FOUND)
  string (REGEX REPLACE "include(/nifti)?/?" "" NIFTI_PREFIX "${NIFTI_INCLUDE_DIR}")
  set (NIFTI_DIR "${NIFTI_PREFIX}" CACHE PATH "Installation prefix for NIFTI." FORCE)
  unset (NIFTI_PREFIX)
endif ()
