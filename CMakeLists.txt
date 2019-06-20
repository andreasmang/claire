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
find_package(PETSc REQUIRED)
find_package(FFTW REQUIRED)
find_package(ACCFFT REQUIRED)
find_package(NIFTI REQUIRED)
find_package(PNETCDF REQUIRED)
find_package(ZLIB REQUIRED)


