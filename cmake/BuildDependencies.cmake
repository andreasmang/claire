INCLUDE(ExternalProject)

ExternalProject_Add(petsclib
  URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz
  #GIT_REPOSITORY https://bitbucket.org/petsc/petsc
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/petsc
  SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/petsc
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND "${CMAKE_CURRENT_BINARY_DIR}/external/configure"
  BUILD_COMMAND make
  INSTALL_COMMAND make install
)
