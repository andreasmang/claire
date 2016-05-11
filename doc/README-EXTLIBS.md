# Installing Dependencies

COLDREG depends on the following libraries:

* [FFTW](http://www.fftw.org) (version 3.3.4)
* [ACCFFT](http://accfft.org) (requires FFTW)
* [PETSc](https://www.mcs.anl.gov/petsc/) (version 3.7; requires [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/))
* [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) (version 2.0.0; requires zlib)

These need to be installed and made available on your system before compiling the code (instruction for compiling COLDREG can be found in [doc/README-INSTALLATION.md](README-INSTALLATION.md)).

## Before Compiling

Make sure the standard MPI wrappers for `mpicc` and `mpicxx` are available on your system (either by loading the right modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions). For instance, add the following definitions to your `~/.bashrc`:

```bash
export PATH=/path/to/mpicxx:/path/to/mpicc:${PATH}
export LD_LIBRARY_PATH=/path/to/mpi/lib/${LD_LIBRARY_PATH}
```

## External Libraries/Dependencies

In general you should have received tarball files for the individual libraries. The **compressed** tarball files (`LIBNAME.tar.gz`) should already be in or be added to the [external](../external) folder.



## Building Dependencies

### General Info

I recommend to add the path to the **FFTW library** and the **PETSc library** to your `LD_LIBRARY_PATH` (lib folders for FFTW and PETSc) as well. If you decide to use PETSc with your local MKL implementation, also add the corresponding path to the `LD_LIBRARY_PATH`.


### Quick Shot

If `mpicc` and `mpicxx` are available, you can install all external dependencies at once as follows:

```bash
cd external
./build_libs.sh --build
```


### More Detailed Instructions

The libraries can be compiled by running the [build_libs.sh](../external/build_libs.sh) script in the [external](../external) subdirectory. For options do

    ./build_libs.sh --help

This will provide information on what parameters to parse. Ideally it should be sufficient to do

    ./build_libs.sh --build

This will install all libraries in a local folder called "lib" ([external/lib](../external/lib/)). You can also build the individual libraries one after another. For more information do `./build_libs.sh --help`. 

If the wrapper for your MPI implementation does **not** include `mpicc` and `mpicxx` you will have to pass the MPI compiler manually (**not tested**)

    ./build_libs.sh --cxx YOURMPICXXCOMPILER --c YOURMPICCOMPILER

The script will figure out the path from the binary you set. An example for `mpicxx` and `mpicc` is

    ./build_libs.sh --cxx mpicxx --c mpicc

Please check the `cmake`, `make` and `automake` outputs for errors. To check if everything worked you can also inspect the "build" subfolders of the individual libraries in the "lib" folder (subdirectories of [external](../external)). See if folders in "build" were created and the library and include files exist.


### Adding Libraries to System

Before you are able to compile and run COLDREG you will have to add some **environment variables** to your system. When building the libraries a file called `environment_vars.sh` is created. This file should be located in [external/libs](../external/libs). To add the corresponding environment variables temporarily (for the current session) to your system, do

   source environment_vars.sh

or copy its content of `environment_vars.sh` to your `~/.bashrc`.


## Compiling the Code

See [doc/README-INSTALLATION.md](README-INSTALLATION.md)



## Additional Info


### NIFTICLIB

* file: [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0/)
* description: library to read and write NIFTI images
* see [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) (requires zlib)


### FFTW

* file: [fftw-3.3.4.tar.gz](ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz)
* description: library for computing FFTs
* see [FFTW](http://www.fftw.org) (version 3.3.4)


### ACCFFT

* file: [accfft.tar.gz](https://github.com/amirgholami/accfft)
* git: `git clone git@github.com:amirgholami/accfft.git`
* description: library to compute FFT in parallel (requires FFTW)


### PETSc

* file: [petsc-lite-3.7.0.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz)
* [bitbucket](https://bitbucket.org/petsc/petsc)
* description: library for numerics and optimization

