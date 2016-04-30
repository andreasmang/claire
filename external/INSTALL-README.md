ADD external dependencies (tarballs) to this folder.

To compile the code you have to do the following:
    1) compile external libraries
    2) init or add environment variables
    3) compile registration binary



########################################################
# 1) compile external libraries 
########################################################
The libraries can be compiled by running the

    build_libs.sh

script. For a more detailed descriptions do

    ./build_libs.sh --help

This provides information on what parameters to parse. In general it should be suifficient to do

    ./build_libs.sh --build

This will install all libraries in a local folder called "lib".

If your mpi compilers are NOT

   mpicc     and     mpicxx

you will have to pass the mpi compiler manually via

    ./build_libs.sh --cxx YOURMPICXXCOMPILER --c YOURMPICCOMPILER

The script will figure out the path from the binary you set. An example for "mpicxx" and "mpicc" is

    ./build_libs.sh --cxx mpicxx --c mpicc


Please check the cmake, make  and configure outputs for errors. A simple check if everything worked is to inspect the "build" subfolders for the individual liberaries in "lib". See if folders in "build" were created and the libfiles exist.



########################################################
# 2) init and add external variables
########################################################

After step 1) there also should be a file called

    environment_vars.sh

in the "lib" folder. To run and compile the code (via make) either do

   source environment_vars.sh

or copy its entries to your ~/.bashrc.



########################################################
# 2) compile registration binary
########################################################

Go to the top level directory of the code and run

    make

This should compile the code and create the desired binary in the "bin" folder. To get help about the binary, just do

   ./runcoldreg -help

If you get a message about a missing "petsc" library when running the code, you probably forgot to source the

    environment_vars.sh

file. A quick fix to not have to redo this is to add its content to your ~/.bashrc.



########################################################
# libraries
########################################################

The dependencies are

----------------------------
# NIFTICLIB:
----------------------------
file: nifticlib-2.0.0.tar.gz
url: https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0/
description: lib to read and write nifti images


----------------------------
# FFTW
----------------------------
file: fftw-3.3.4.tar.gz
url: ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz
description: lib for computing FFTs 


----------------------------
# PNETCDF
----------------------------
file: parallel-netcdf-1.7.0.tar.gz
url: http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.7.0.tar.gz
description: lib for reading and writing *.nc files (for paraview)


----------------------------
# ACCFFT
----------------------------
file: accfft.tar.gz
git: git@github.com:amirgholami/accfft.git
note: requires PNETCDF and FFTW
description: lib to compute FFT in parallel


----------------------------
# PETSc
----------------------------
file: petsc-lite-3.7.0.tar.gz
git: https://bitbucket.org/petsc/petsc
url: http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz
description: lib for numerics and optimization

