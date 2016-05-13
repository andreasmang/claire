# add binary tree directory to search path for include files
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# add executable and link libraries
ADD_EXECUTABLE ( runcoldreg runcoldreg.cpp )

MESSAGE(STATUS "libraries in binary ${COLD_EXTERNAL_LIBS}")
TARGET_LINK_LIBRARIES( runcoldreg registration ${COLD_EXTERNAL_LIBS} )
INSTALL(TARGETS runcoldreg DESTINATION bin)
