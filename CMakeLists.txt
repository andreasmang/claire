CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)

# if user specifies -DCMAKE_BUILD_TYPE in command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set 
# CMAKE_BUILD_TYPE to Debug prior to calling PROJECT()

IF(CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "build type")
ELSE()
	SET(CMAKE_BUILD_TYPE Debug CACHE STRING "build type")
ENDIF()

PROJECT(COLD)
SET(COLD_VERSION_MAJOR 0)
SET(COLD_VERSION_MINOR 0)
SET(COLD_VERSION_PATCH 0)
SET(COLD_VERSION ${COLD_VERSION_MAJOR}.${COLD_VERSION_MINOR}.${COLD_VERSION_PATCH})

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

MESSAGE(STATUS "compiling COLD v${COLD_VERSION_MAJOR}.${COLD_VERSION_MAJOR}.${COLD_VERSION_PATCH}")
MESSAGE(STATUS "source folder: ${COLD_SOURCE_DIR}")


# run main cmake file
INCLUDE(${COLD_SOURCE_DIR}/cmake/COLD.cmake)





#SET( COLD_LIB_LIST core registration numerics )

# process lib directories 
#FOREACH( libpack ${COLD_LIB_LIST} )
#	ADD_SUBDIRECTORY( ${COLD_SOURCE_DIR}/${libpack} lib/${libpack} )
#ENDFOREACH()

# process bin directories 
#ADD_SUBDIRECTORY(${COLD_SOURCE_DIR}/apps apps)
