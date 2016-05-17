# Installing and Running COLDREG




## Content

* Requirements
* Installation of Dependencies
	* General Overview
	* Installing Dependencies: Quick Shot
	* Installing Dependencies: Details
	* Environment Variables
* Building COLDREG
* Running COLDREG
* FAQ
* More info about the Dependencies




## Requirements

* MPI (Open MPI; MVAPICH; Intel MPI); required by ACCFFT, PETSc, and COLDREG
* cmake [https://cmake.org](https://cmake.org); required by *ACCFFT* and *niftilib*
* python ([https://www.python.org](https://www.python.org)); version 2.7; required by PETSc


Make sure that the standard *wrappers* for *mpicc* and *mpicxx* are available on your system (either by loading the appropriate modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions; see below). The compilation has been tested for *Open MPI*, *MVAPICH*, and *Intel MPI*. An **FAQ** that summarizes common issues encountered when building the code can be found on the bottom of the page. The dependencies are listed below.




## Installing Dependencies


### General Overview

COLDREG depends on the following libraries:

* FFTW [http://www.fftw.org](http://www.fftw.org)
* ACCFFT [http://accfft.org](http://accfft.org) 
* PETSc [https://www.mcs.anl.gov/petsc/](https://www.mcs.anl.gov/petsc/) 
* zlib [http://zlib.net](http://zlib.net) 
* niftilib [https://sourceforge.net/projects/niftilib/files/nifticlib/](https://sourceforge.net/projects/niftilib/files/nifticlib/) 

We provide the following libraries with the code: FFTW version 3.3.4; ACCFFT (downloaded in Mai, 2016), PETSc (version 3.7; requires *python 2.7*); zlib (version 1.2.8); and niftilib (version 2.0.0). These libraries have to be installed and made available on your system before compiling the code. We build all libraries as **static** by default. We provide *tarball* files for these individual libraries. The *compressed* tarball files (i.e, *LIBRARYNAME.tar.gz*) should be located in or be added to the [external](../external) folder.


#### Installing Dependencies: Quick Shot

```bash
cd external
./build_libs.sh --build
```

If you use *Intel MPI* provide the `--useimpi` option to `./build_libs.sh`.


#### Installing Dependencies: Details

The libraries can be compiled by running the [build_libs.sh](../external/build_libs.sh) script in the [external](../external) subdirectory. To see all the options do

```bash
./build_libs.sh --help
```

This will provide information on what parameters you can parse. Ideally it should be sufficient to do `./build_libs.sh --build`. You can also build the individual libraries one after another, via the `--bLIBNAME` option, where `LIBNAME` is the name of the libarary. For precise instructions, do `./build_libs.sh --help`. If you want to clean up the libraries folder, you can do `./build_libs.sh --clean`. The *build* folders will be removed each time you recompile the libraries.

Please check the `cmake`, `make` and `automake` outputs for errors. To check if everything worked you can also take a look at the "build" subdirectories of the individual libraries in the "lib" folder (subdirectories of [external](../external)). See if folders in "build" were created and the library and include files exist.


#### Environment Variables 

Before you are able to compile and run COLDREG you will have to add some *environment variables* to your system. When building the libraries a file called `environment_vars.sh` is created. This file should be located in [external/libs](../external/libs). To add the corresponding environment variables temporarily (for the current session) to your system, do

```bash
source environment_vars.sh
```

To add them permanently, copy the content of `environment_vars.sh` to your `~/.bashrc`. Note that the script will define *absolute paths*.




## Building COLDREG

Before you can build COLDREG you need to

* make sure that you have installed all *dependencies*
* check the [makefile](makefile) before building the code:
	* if you use an *intel compiler* set the `USEINTEL` flag to `yes`
	* if you use a *GNU compiler* set the `USEINTEL` flag to `no`
	* if you use *Intel MPI* (impi) set the `USEINTELMPI` flag to `yes` (if not, set it to `no`)
* make sure all paths needed in the makefile are available on your system (to check, you can do `env` in your bash); to add the paths temporarily `source external/libs/environment_vars.sh` or add the content of `external/libs/environment_vars.sh` to your `~/.bashrc`

To build the code using the `make` system do (in the top level directory):

```bash
make -j
```




## Running COLDREG

If everything compiled correctly, you can run a test example by doing:

```bash
./bin/runcoldreg
```

To get a general idea on how to run the binary do:

```bash
./bin/runcoldreg -help
```

For more advanced options do:

```bash
./bin/runcoldreg -advanced
```

You can also find a list of the available options for the binary in [doc/help.txt](help.txt) and [doc/advanced-help.txt](advanced-help.txt).




## FAQ

1. `libstdc++` or `GLIBCXX_ ...` not found
	* required by *ACCFFT* and *niftilib*
	* `libstdc++` is the *GNU Standard C++ library*
	* you might be able to locate this library via `locate libstdc`
	* **fix**: if `gcc` is on your system, you need to add it to the `LD_LIBRARY_PATH`; add `export LD_LIBRARY_PATH=/path/to/gcc/lib(64):$LD_LIBRARY_PATH` to your `bashrc`
2. PETSc requires at least *python 2.7* for compilation; *python 3* is not supported in PETSc 3.7
3. other dependencies (should in general be available on your system)
	* cmake ([https://cmake.org](https://cmake.org); required by ACCFFT and niftilib; version 2.8 or greater)
	* openMP (required by ACCFFT)
	* `crypt` library (required by PETSc)
	* `ssl` library (required by PETSc)
	* BLAS (required by PETSc; we install it along with PETSc; [http://www.netlib.org/blas/](http://www.netlib.org/blas/))
	* LAPACK (required by PETSc; we install it along with PETSc; [http://www.netlib.org/lapack/](http://www.netlib.org/lapack/))




## More Info about the Dependencies


### FFTW

* file: [fftw-3.3.4.tar.gz](ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz)
* description: library for computing FFTs
* additional details can be found at [http://www.fftw.org](http://www.fftw.org) (version 3.3.4)


### ACCFFT

* file: accfft.tar.gz
* source code also available on gitub: [https://github.com/amirgholami/accfft](https://github.com/amirgholami/accfft)
* description: library to compute FFT in parallel (requires FFTW)
* additional details can be found at [http://www.accfft.org](http://www.accfft.org)
* requires `FFTW` to be installed


### PETSc

* file: [petsc-lite-3.7.0.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz)
* source code also available on bitbucket: [https://bitbucket.org/petsc/petsc](https://bitbucket.org/petsc/petsc)
* description: library for numerics and optimization


### NIFTICLIB

* file: [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0/)
* description: library to read and write NIFTI images
* see [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) 
