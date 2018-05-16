[0]PETSC ERROR: ------------------------------------------------------------------------
[0]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation, probably memory access out of range
[0]PETSC ERROR: Try option -start_in_debugger or -on_error_attach_debugger
[0]PETSC ERROR: or see http://www.mcs.anl.gov/petsc/documentation/faq.html#valgrind
[0]PETSC ERROR: or try http://valgrind.org on GNU/linux and Apple Mac OS X to find memory corruption errors
[0]PETSC ERROR: configure using --with-debugging=yes, recompile, link, and run 
[0]PETSC ERROR: to get more information on the crash.
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Signal received
[0]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
[0]PETSC ERROR: Petsc Release Version 3.9.1, unknown 
[0]PETSC ERROR: Unknown Name on a arch-cuda-single named c221-303.maverick.tacc.utexas.edu by naveen15 Wed May 16 18:30:03 2018
[0]PETSC ERROR: Configure options --with-cuda=1 --download-cusp=yes --with-precision=single --download-f2cblaslapack=1 --with-mpi=1 --with-mpi-dir=/home/04716/naveen15/claire/external/libs/openmpi-3.0.1 --with-debugging=0 --CUDAFLAGS=-arch=sm_35 --with-cuda-dir=/opt/apps/cuda/8.0 --with-ssl=0 --with-64-bit-indices --with--shared --with-x COPTFLAGS="-O3" CXXOPTFLAGS="-O3" PETSC_ARCH=arch-cuda-single
[0]PETSC ERROR: #1 User provided function() line 0 in  unknown file
--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD
with errorcode 59.

NOTE: invoking MPI_ABORT causes Open MPI to kill all MPI processes.
You may or may not see output from other processes, depending on
exactly when Open MPI kills them.
--------------------------------------------------------------------------
