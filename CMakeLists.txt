cmake_minimum_required (VERSION 3.8)

project(claire)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
# CMP0075 Include file check macros honor CMAKE_REQUIRED_LIBRARIES
if(POLICY CMP0075)
    cmake_policy(SET CMP0075 NEW)
endif()

# find_package packages
include(basicFind)
include(FindMPI)
find_package(MPI REQUIRED)
find_package(PETSc 3.6 REQUIRED)
set(PETSC_VERSIONS "")
set(PETSC_VERSION_MAJOR "")
set(PETSC_VERSION_MINOR "")
string(REGEX MATCHALL "[0-9]+" PETSC_VERSIONS ${PETSC_VERSION})
list(GET PETSC_VERSIONS 0 PETSC_VERSION_MAJOR)
list(GET PETSC_VERSIONS 1 PETSC_VERSION_MINOR)
find_package(FFTW REQUIRED)
find_package(ACCFFT REQUIRED)
find_package(NIFTI REQUIRED)
find_package(PNETCDF REQUIRED)
find_package(ZLIB REQUIRED)

# add sources
file(GLOB_RECURSE CLAIRE_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM CLAIRE_SOURCES "${PROJECT_SOURCE_DIR}/src/DomainDecomposition.cpp")
list(REMOVE_ITEM CLAIRE_SOURCES "${PROJECT_SOURCE_DIR}/src/Logger.cpp")

# add target
add_library(claireobj OBJECT ${CLAIRE_SOURCES})

# set include directories
target_include_directories(claireobj PUBLIC "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/deps/3rdparty" "${PROJECT_SOURCE_DIR}/deps/3rdparty/libmorton")
target_include_directories(claireobj PUBLIC ${PETSC_INCLUDES} ${FFTW_INCLUDES} ${ACCFFT_INCLUDES} ${NIFTI_INCLUDES} ${PNETCDF_INCLUDE_DIRS})

# set compiler stuff
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_FLAGS "-O3 -ansi -xhost -parallel -qopenmp -std=c++11 -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wwrite-strings")
target_compile_definitions(claireobj PUBLIC REG_HAS_NIFTI REG_HAS_PNETCDF)
#set_target_properties(claireobj PROPERTIES CXX_STANDARD 11)

# install all targets
install(TARGETS claireobj)



