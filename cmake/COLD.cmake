# prevent in build installs
STRING(TOLOWER "${CMAKE_INSTALL_PREFIX}" _PREFIX)
STRING(TOLOWER "${${CMAKE_PROJECT_NAME}_BINARY_DIR}" _BUILD)
IF("${_PREFIX}" STREQUAL "${_BUILD}")
  MESSAGE(FATAL_ERROR "Current CMAKE_INSTALL_PREFIX points at build tree:\n ${CMAKE_INSTALL_PREFIX}\n This is not supported.")
ENDIF()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

# prevent in source builds
INCLUDE(${COLD_SOURCE_DIR}/cmake/PreventInSourceBuild.cmake)

SET(CMAKE_BUILD_TYPE)

# set path for cmake modules
SET(CMAKE_MODULE_PATH ${OCREG_SOURCE_DIR}/cmake)

SET(COLD_EXTERNAL_LIBS)
SET(COLD_EXTERNAL_INCS)

# set warning level for compiler
SET( CMAKE_CXX_WARNING_LEVEL 5 )
#SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
#SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3")

FIND_PACKAGE(MPI REQUIRED)
IF (MPI_C_FOUND AND MPI_CXX_FOUND)
    MESSAGE(STATUS "MPI library found")
    LIST( APPEND COLD_EXTERNAL_INCS ${MPI_INCLUDE_PATH} )
    LIST( APPEND COLD_EXTERNAL_LIBS ${MPI_LIBRARIES} )
    SET( CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} )
    SET( CMAKE_C_COMPILER ${MPI_C_COMPILER} )
ELSE( )
    MESSAGE(FATAL_ERROR "FATAL ERROR: MPI library could not be found")
ENDIF( )


FIND_PACKAGE(OpenMP REQUIRED)
IF (OPENMP_FOUND)
    MESSAGE(STATUS "OpenMP library found")
    SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} )
    SET( CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} )
ELSE( )
    MESSAGE(FATAL_ERROR "FATAL ERROR: OpenMP library could not be found")
ENDIF( )


include(FindMPI)

# set cmake compiler to mpi compiler wrapper found by FindMPI (usually mpicxx) 
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )


ADD_DEFINITIONS(
"-Wall
-Wshadow
-Wunused-variable
-Wunused-parameter
-Wunused-function
-Wunused
-Woverloaded-virtual
-Wwrite-strings
-Wno-deprecated
-Wformat
-Woverloaded-virtual
-ansi")



# build external dependencies (PETSc, FFTW, ACCFFT, PNETCDF)
INCLUDE(${COLD_SOURCE_DIR}/cmake/BuildDependencies.cmake)

#IF( COLD_EXTERNAL_INCS )
#    LIST( REMOVE_DUPLICATES COLD_EXTERNAL_INCS)
#    LIST( REMOVE_DUPLICATES COLD_EXTERNAL_LIBS)
#ENDIF( COLD_EXTERNAL_INCS )


INCLUDE_DIRECTORIES( ${COLD_SOURCE_DIR}/include )
INCLUDE_DIRECTORIES( ${COLD_BINARY_DIR}/include )
INCLUDE_DIRECTORIES( ${COLD_BINARY_DIR} )

# add all directories
INCLUDE_DIRECTORIES(${COLD_EXTERNAL_INCS})


MESSAGE(STATUS "CMAKE_CXX_FLAGS are ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "CMAKE_C_FLAGS are ${CMAKE_C_FLAGS}")
MESSAGE(STATUS "cmake compiler: ${CMAKE_CXX_COMPILER}")
MESSAGE(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "external libraries: ${COLD_EXTERNAL_LIBS}")
MESSAGE(STATUS "external includes: ${COLD_EXTERNAL_INCS}")


SET(COLD_SRCS
${COLD_SOURCE_DIR}/src/ghost.cpp
${COLD_SOURCE_DIR}/src/interp3.cpp
${COLD_SOURCE_DIR}/src/Interp3_Plan.cpp
${COLD_SOURCE_DIR}/src/RegOpt.cpp
${COLD_SOURCE_DIR}/src/RegUtils.cpp
${COLD_SOURCE_DIR}/src/DataReadWriteRegistration.cpp
${COLD_SOURCE_DIR}/src/LargeDeformationRegistration.cpp
${COLD_SOURCE_DIR}/src/OptimalControlRegistration.cpp
${COLD_SOURCE_DIR}/src/OptimalControlRegistrationIC.cpp
${COLD_SOURCE_DIR}/src/OptimizationProblemRegistration.cpp
${COLD_SOURCE_DIR}/src/PreProcessingRegistration.cpp
${COLD_SOURCE_DIR}/src/RegularizationRegistration.cpp
${COLD_SOURCE_DIR}/src/SemiLagrangian.cpp
${COLD_SOURCE_DIR}/src/SynProbRegistration.cpp
${COLD_SOURCE_DIR}/src/TaoInterfaceRegistration.cpp
${COLD_SOURCE_DIR}/src/VecField.cpp
)
#${COLD_SOURCE_DIR}/src/gpu_interp3.cpp
#${COLD_SOURCE_DIR}/src/read_binary.cpp
#${COLD_SOURCE_DIR}/src/DataOut.cpp
#${COLD_SOURCE_DIR}/src/SemiLagrangianGPU.cpp

LINK_DIRECTORIES( ${COLD_BINARY_DIR}/lib )

# build library
#IF(BUILD_SHARED_LIBS)
#    ADD_LIBRARY(coldreg SHARED ${COLD_SRCS})
#ELSE(BUILD_SHARED_LIBS)
#    ADD_LIBRARY(coldreg STATIC ${COLD_SRCS})
#ENDIF(BUILD_SHARED_LIBS)

ADD_LIBRARY(coldreg ${COLD_SRCS})

# link libraries
TARGET_LINK_LIBRARIES( coldreg ${COLD_EXTERNAL_LIBS} )

INSTALL(TARGETS coldreg DESTINATION lib)

