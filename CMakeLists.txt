cmake_minimum_required (VERSION 3.0)

project(claire)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# includes
include(basicFind)
include(CopyTargetProperty)
include(FindMPI)

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
# CMP0075 Include file check macros honor CMAKE_REQUIRED_LIBRARIES
if(POLICY CMP0075)
    cmake_policy(SET CMP0075 NEW)
endif()

# find_package packages
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# set optimization level
set(CMAKE_CXX_FLAGS "-O3 -ansi ${OpenMP_CXX_FLAGS}")

# check for intel compiler
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(USE_INTEL "yes")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xhost -parallel")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# set PEDANTIC flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wwrite-strings")

find_package(PETSc 3.6 REQUIRED)
set(PETSC_VERSIONS "")
set(PETSC_VERSION_MAJOR "")
set(PETSC_VERSION_MINOR "")
string(REGEX MATCHALL "[0-9]+" PETSC_VERSIONS ${PETSC_VERSION})
list(GET PETSC_VERSIONS 0 PETSC_VERSION_MAJOR)
list(GET PETSC_VERSIONS 1 PETSC_VERSION_MINOR)
if (${PETSC_VERSION_MINOR} GREATER "9")
    message(FATAL_ERROR "Claire does not support PETSc version >= 3.10. Use PETSc version < 3.10")
endif()

find_package(FFTW REQUIRED)
# for FFTW, if not using single precision, exclude the coressponding libraries
if (NOT USE_SINGLE_PRECISION)
    list(FILTER FFTW_LIBRARIES EXCLUDE REGEX "fftw3f+")
endif()

find_package(ACCFFT REQUIRED)
find_package(NIFTI REQUIRED)
if (USE_PNETCDF)
    find_package(PNETCDF REQUIRED)
endif()
find_package(ZLIB REQUIRED)

# add sources
file(GLOB_RECURSE CLAIRE_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM CLAIRE_SOURCES "${PROJECT_SOURCE_DIR}/src/DomainDecomposition.cpp")
list(REMOVE_ITEM CLAIRE_SOURCES "${PROJECT_SOURCE_DIR}/src/Logger.cpp")

# add target
if (BUILD_SHARED)
    add_library(registration SHARED ${CLAIRE_SOURCES})
else()
    add_library(registration STATIC ${CLAIRE_SOURCES})
endif()

# set include directories
target_include_directories(registration PUBLIC "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/deps/3rdparty" "${PROJECT_SOURCE_DIR}/deps/3rdparty/libmorton")
target_include_directories(registration PUBLIC ${PETSC_INCLUDES} ${FFTW_INCLUDES} ${ACCFFT_INCLUDES} ${NIFTI_INCLUDES} ${PNETCDF_INCLUDE_DIRS})

target_compile_definitions(registration PUBLIC REG_HAS_NIFTI)
if (${USE_PNETCDF})
    target_compile_definitions(registration PUBLIC  REG_HAS_PNETCDF)
endif()
if (USE_HASWELL AND USE_SINGLE_PRECISION)
    target_compile_definitions(registration PUBLIC HASWELL)
endif()
set_target_properties(registration PROPERTIES CXX_STANDARD 11 ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

target_link_libraries(registration ${PETSC_LIBRARIES} ${NIFTI_LIBRARIES} ${ZLIB_LIBRARIES} ${ACCFFT_LIBRARIES} ${PNETCDF_LIBRARIES} ${FFTW_LIBRARIES})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS registration DESTINATION lib)

# add app executables and install
file( GLOB APP_SOURCES "${PROJECT_SOURCE_DIR}/apps/*.cpp" )
foreach( appsourcefile ${APP_SOURCES} )
    get_filename_component( app ${appsourcefile} NAME_WE "${PROJECT_SOURCE_DIR}/apps")
    add_executable( ${app} ${appsourcefile} )
    add_dependencies( ${app} registration )
    target_link_libraries (${app} ${PETSC_LIBRARIES} ${NIFTI_LIBRARIES} ${ZLIB_LIBRARIES} ${ACCFFT_LIBRARIES} ${PNETCDF_LIBRARIES} ${FFTW_LIBRARIES} registration)
    if (${USE_INTEL})
        target_link_libraries (${app} "imf")
        target_link_libraries (${app} "m")
    endif()
    set_target_properties(${app} PROPERTIES CXX_STANDARD 11)
    copy_target_property(registration ${app} INCLUDE_DIRECTORIES)
    copy_target_property(registration ${app} COMPILE_OPTIONS)
    copy_target_property(registration ${app} COMPILE_DEFINITIONS)
    install( TARGETS ${app} 
             DESTINATION bin )
endforeach( appsourcefile ${APP_SOURCES} )

