IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	SET(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/bin" CACHE PATH "install path prefix" FORCE)
ENDIF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


# create a configure file (used to parse variables/flags to code)
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/ColdConfig.h.in" "${PROJECT_BINARY_DIR}/ColdConfig.h")
INSTALL(FILES "${PROJECT_BINARY_DIR}/ColdConfig.h" DESTINATION include)

# add binary tree to search path for include files
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# make sure include files can be found
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)

# prevent in source builds
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake/PreventInSourceBuild.cmake )

# init external library and include variables
SET(COLD_EXTERNAL_LIBS)
SET(COLD_EXTERNAL_INCS)

# set warning level for compiler
SET( CMAKE_CXX_WARNING_LEVEL 5 )
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3")


# find MPI (required)
FIND_PACKAGE(MPI REQUIRED)
IF (MPI_C_FOUND AND MPI_CXX_FOUND)
    LIST( APPEND COLD_EXTERNAL_INCS ${MPI_INCLUDE_PATH} )
    LIST( APPEND COLD_EXTERNAL_LIBS ${MPI_LIBRARIES} )
    SET( CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} )
    SET( CMAKE_C_COMPILER ${MPI_C_COMPILER} )
ELSE( )
    MESSAGE(FATAL_ERROR "FATAL ERROR: MPI library could not be found")
ENDIF( )


# find OpenMP (required by accfft) 
FIND_PACKAGE(OpenMP REQUIRED)
IF (OPENMP_FOUND)
    SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} )
    SET( CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} )
    SET( CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} )
ELSE()
    MESSAGE(FATAL_ERROR "FATAL ERROR: OpenMP library could not be found")
ENDIF()


# search for zlib (if not found, we compile it)
FIND_PACKAGE( ZLIB )
IF ( ZLIB_FOUND )
    INCLUDE_DIRECTORIES( ${ZLIB_INCLUDE_DIRS} )
    LIST( APPEND COLD_EXTERNAL_LIBS ${ZLIB_LIBRARIES} )
ENDIF( ZLIB_FOUND )



#FIND_PACKAGE( MKL REQUIRED )
#IF ( MKL_FOUND )
#    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
#    LIST( APPEND COLD_EXTERNAL_LIBS ${MKL_LIBRARIES} )
#ENDIF( MKL_FOUND )



# build external dependencies (PETSc, FFTW, ACCFFT, nifti)
INCLUDE(${PROJECT_SOURCE_DIR}/cmake/BuildDependencies.cmake)

LIST( APPEND COLD_EXTERNAL_LIBS -lssl -lcrypto )


#IF( COLD_EXTERNAL_INCS )
#    LIST( REMOVE_DUPLICATES COLD_EXTERNAL_INCS)
#    LIST( REMOVE_DUPLICATES COLD_EXTERNAL_LIBS)
#ENDIF( COLD_EXTERNAL_INCS )


MESSAGE(STATUS "CMAKE_CXX_FLAGS are ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "CMAKE_C_FLAGS are ${CMAKE_C_FLAGS}")
MESSAGE(STATUS "cmake compiler: ${CMAKE_CXX_COMPILER}")
MESSAGE(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "external libraries: ${COLD_EXTERNAL_LIBS}")
MESSAGE(STATUS "external includes: ${COLD_EXTERNAL_INCS}")


SET(COLD_SRCS
${PROJECT_SOURCE_DIR}/src/ghost.cpp
${PROJECT_SOURCE_DIR}/src/interp3.cpp
${PROJECT_SOURCE_DIR}/src/Interp3_Plan.cpp
${PROJECT_SOURCE_DIR}/src/RegOpt.cpp
${PROJECT_SOURCE_DIR}/src/RegUtils.cpp
${PROJECT_SOURCE_DIR}/src/Optimizer.cpp
${PROJECT_SOURCE_DIR}/src/ReadWriteReg.cpp
${PROJECT_SOURCE_DIR}/src/TaoInterfaceRegistration.cpp
${PROJECT_SOURCE_DIR}/src/OptProbRegistration.cpp
${PROJECT_SOURCE_DIR}/src/LargeDeformationRegistration.cpp
${PROJECT_SOURCE_DIR}/src/OptimalControlRegistration.cpp
${PROJECT_SOURCE_DIR}/src/OptimalControlRegistrationIC.cpp
${PROJECT_SOURCE_DIR}/src/OptimalControlRegistrationRIC.cpp
${PROJECT_SOURCE_DIR}/src/PreProcessingRegistration.cpp
${PROJECT_SOURCE_DIR}/src/RegularizationRegistration.cpp
${PROJECT_SOURCE_DIR}/src/RegularizationRegistrationH1.cpp
${PROJECT_SOURCE_DIR}/src/RegularizationRegistrationH1SN.cpp
${PROJECT_SOURCE_DIR}/src/RegularizationRegistrationH2.cpp
${PROJECT_SOURCE_DIR}/src/RegularizationRegistrationH2SN.cpp
${PROJECT_SOURCE_DIR}/src/SemiLagrangian.cpp
${PROJECT_SOURCE_DIR}/src/SynProbRegistration.cpp
${PROJECT_SOURCE_DIR}/src/VecField.cpp
)


SET(COLD_INCS
${PROJECT_SOURCE_DIR}/include/utils.hpp
${PROJECT_SOURCE_DIR}/include/interp3.hpp
${PROJECT_SOURCE_DIR}/include/interp3_common.hpp
${PROJECT_SOURCE_DIR}/include/LargeDeformationRegistration.hpp
${PROJECT_SOURCE_DIR}/include/OptimalControlRegistration.hpp
${PROJECT_SOURCE_DIR}/include/OptimalControlRegistrationIC.hpp
${PROJECT_SOURCE_DIR}/include/OptimalControlRegistrationRIC.hpp
${PROJECT_SOURCE_DIR}/include/OptProbRegistration.hpp
${PROJECT_SOURCE_DIR}/include/Optimizer.hpp
${PROJECT_SOURCE_DIR}/include/PreProcessingRegistration.hpp
${PROJECT_SOURCE_DIR}/include/ReadWriteReg.hpp
${PROJECT_SOURCE_DIR}/include/RegOpt.hpp
${PROJECT_SOURCE_DIR}/include/RegularizationRegistration.hpp
${PROJECT_SOURCE_DIR}/include/RegularizationRegistrationH1.hpp
${PROJECT_SOURCE_DIR}/include/RegularizationRegistrationH1SN.hpp
${PROJECT_SOURCE_DIR}/include/RegularizationRegistrationH2.hpp
${PROJECT_SOURCE_DIR}/include/RegularizationRegistrationH2SN.hpp
${PROJECT_SOURCE_DIR}/include/RegUtils.hpp
${PROJECT_SOURCE_DIR}/include/SemiLagrangian.hpp
${PROJECT_SOURCE_DIR}/include/SemiLagrangianGPU.hpp
${PROJECT_SOURCE_DIR}/include/SynProbRegistration.hpp
${PROJECT_SOURCE_DIR}/include/TaoInterfaceRegistration.hpp
${PROJECT_SOURCE_DIR}/include/VecField.hpp
)


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/includes)
INCLUDE_DIRECTORIES(${COLD_EXTERNAL_INCS})
ADD_LIBRARY(registration STATIC ${COLD_SRCS} ${COLD_INCS})

# link libraries
TARGET_LINK_LIBRARIES( registration ${COLD_EXTERNAL_LIBS} )
INSTALL( TARGETS registration DESTINATION lib )
LINK_DIRECTORIES( ${PROJECT_BINARY_DIR}/bin/lib )

# build the binaries
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/apps)
