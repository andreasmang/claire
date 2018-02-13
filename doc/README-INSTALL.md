# Installing and Running CLAIRE




## Content

* Requirements
* Installation of Dependencies
	* General Overview
	* Installing Dependencies: Quick Shot
	* Installing Dependencies: Details
	* Environment Variables
* Building CLAIRE
* Running CLAIRE
* FAQ
* More info about the Dependencies




## Requirements

* MPI (Open MPI; MVAPICH; Intel MPI); required by ACCFFT, PETSc, and CLAIRE
* cmake [https://cmake.org](https://cmake.org); required by *ACCFFT* and *niftilib*
* python ([https://www.python.org](https://www.python.org)); version 2.7; required by PETSc


Make sure that the standard *wrappers* for *mpicc* and *mpicxx* are available on your system (either by loading the appropriate modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions; see below). The compilation has been tested for *Open MPI*, *MVAPICH*, and *Intel MPI*. An **FAQ** that summarizes common issues encountered when building the code can be found on the bottom of the page. The dependencies are listed below.




## Installing Dependencies


### General Overview

CLAIRE depends on the following libraries:

* FFTW [http://www.fftw.org](http://www.fftw.org)
* ACCFFT [http://accfft.org](http://accfft.org)
* PETSc [https://www.mcs.anl.gov/petsc/](https://www.mcs.anl.gov/petsc/) (requires *python 2.7*)
* niftilib [https://sourceforge.net/projects/niftilib/files/nifticlib/](https://sourceforge.net/projects/niftilib/files/nifticlib/)
* zlib [http://zlib.net](http://zlib.net)
* libmorton [https://github.com/Forceflow/libmorton](https://github.com/Forceflow/libmorton)


To download the libraries we provide a script called 'get_libs.sh' (in the 'external' subfolder). By running this script you can download the libraries: 
```bash
cd external
./get_libs.sh
```

The *tarball* files will be downloaded to this very folder. The *compressed* tarball files (i.e, *LIBRARYNAME.tar.gz*) should remain located in or be added to the [external](../external) folder. Make sure that all libraries are really downloaded (the progress bar of 'wget' should be full). Occasionally, this is not the case. Just delete the corresponding tarball file and run 'get_libs.sh' again. For the second execution it will work. These libraries have to be installed and made available on your system before compiling the code. We build all libraries as **static** by default. 


#### Installing Dependencies: Quick Shot

```bash
cd external
./build_libs.sh --build
```

<If you use *Intel MPI* provide the `--useimpi` option to `./build_libs.sh`.>


#### Installing Dependencies: Details

The libraries can be compiled by running the [build_libs.sh](../external/build_libs.sh) script in the [external](../external) subdirectory. To see all the options do

```bash
./build_libs.sh --help
```

This will provide information on what parameters you can parse. Ideally it should be sufficient to do `./build_libs.sh --build`. You can also build the individual libraries one after another, via the `--bLIBNAME` option, where `LIBNAME` is the name of the library. For precise instructions, do `./build_libs.sh --help`. If you want to clean up the libraries folder, you can do `./build_libs.sh --clean`. The *build* folders will be removed each time you recompile the libraries.

Please check the `cmake`, `make` and `automake` outputs for errors. To check if everything worked you can also take a look at the "build" subdirectories of the individual libraries in the "lib" folder (subdirectories of [external](../external)). See if folders in "build" were created and the library and include files exist.


#### Environment Variables 

Before you are able to compile and run COLDREG you will have to add some *environment variables* to your system. When building the libraries a file called `environment_vars.sh` is created. This file should be located in [external/libs](../external/libs). To add the corresponding environment variables temporarily (for the current session) to your system, do

```bash
source environment_vars.sh
```

To add them permanently, copy the content of `environment_vars.sh` to your `~/.bashrc`. Note that the script will define *absolute paths*.




## Building CLAIRE

Before you can build CLAIRE you need to

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

If you build in parallel using 'make -j', on certain systems to many threads will be used. This will result in compilation errors. To fix this, run 'make -j 12' instead (for quick access, you may want to define an alias in your ~/.bashrc).



## Running CLAIRE

Details on how to run CLAIRE are summarized in [doc/README-RUNME.md](README-RUNME.md). A list of the available options for the binary can be found in [doc/help.txt](help.txt) and [doc/advanced-help.txt](advanced-help.txt).




## FAQ

1. `libstdc++` or `GLIBCXX...` not found
	* required by *ACCFFT* and *niftilib*
	* `libstdc++` is the *GNU Standard C++ library*
	* you might be able to locate this library via `locate libstdc`
	* **fix**: if `gcc` is on your system, you need to add it to the `LD_LIBRARY_PATH`; add `export LD_LIBRARY_PATH=/path/to/gcc/lib(64):$LD_LIBRARY_PATH` to your `bashrc`
2. PETSc requires at least *python 2.7(.11)* for compilation; *python 3* is not supported in PETSc 3.7
3. when compiling PETSc: ERROR:root:code for hash md5 was not found...
	* this is a problem with your python 2.7 installation
	* **fix** upgrade to python 2.7.11
4. g++: error: unrecognized command line option -parallel
	* you probably set `USEINTEL` in the [makefile](../makefile) to `yes`, but are not using an intel compiler
	* **fix**: `USEINTEL=no`
5. ./include/RegUtils.hpp:33:19: fatal error: petsc.h: No such file or directory
	* you probably forgot to set the environment variables
	* **fix**: `source external/libs/environment_vars.sh`
6. definition in ...libmpi_mt.so section .tbss mismatches non-TLS definition in ...libmpi.so.4 section .bss
	* you have used inconsistent MPI libraries during the build (thread safe vs. non-thread safe)
	* **fix**:
		* set `USEINTELMPI` in the [makefile](../makefile) to `yes` 
		* when building the libraries add `--useimpi`: `./build_libs.sh --build --useimpi`
7. libmpi.so: could not read symbols: Bad value (see 6.)
8. libimf.so: warning: warning: feupdateenv is not implemented and will always fail
	* you probably use an intel compiler but did not set `USEINTEL` to `yes`
	* **fix**: set `USEINTEL` in the [makefile](../makefile) to `yes` 
8. other dependencies (should in general be available on your system)
	* cmake (https://cmake.org; required by ACCFFT and niftilib; version 2.8 or greater)
	* openMP (required by ACCFFT)
	* `crypt` library (required by PETSc)
	* `ssl` library (required by PETSc)
	* BLAS (required by PETSc; we install it along with PETSc; http://www.netlib.org/blas/)
	* LAPACK (required by PETSc; we install it along with PETSc; http://www.netlib.org/lapack/)


## Build Systems
CLAIRE has been successfully compiled using the following modules:

### [Opuntia](http://www.uh.edu/cacds/resources/hpc/opuntia/)
```bash
intel/psxe-2016-64-bit
intelmpi/5.1.1
```

```bash
gcc/4.9.2
openmpi/1.8.3-gcc4.9.2
```


Additional packages are:
```bash
python/2.7
cmake/3.0.2
```

## More Info about the Dependencies


### FFTW

* file: [fftw-3.3.6-pl2.tar.gz](http://www.fftw.org/fftw-3.3.6-pl2.tar.gz) 
* description: library for computing FFTs
* additional details can be found at [http://www.fftw.org](http://www.fftw.org) (version 3.3.4)
	* [fftw-3.3.4.tar.gz](http://www.fftw.org/fftw-3.3.4.tar.gz)

### ACCFFT

* file: accfft.tar.gz
* source code also available on gitub: [https://github.com/amirgholami/accfft](https://github.com/amirgholami/accfft)
* description: library to compute FFT in parallel (requires FFTW)
* additional details can be found at [http://www.accfft.org](http://www.accfft.org)
* requires `FFTW` to be installed


### PETSc

* file: [petsc-lite-3.8.3.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.8.3.tar.gz)
* source code also available on bitbucket: [https://bitbucket.org/petsc/petsc](https://bitbucket.org/petsc/petsc)
* description: library for numerics, linear algebra, and optimization
* older versions that have been used successfully:
	* [petsc-lite-3.7.6.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz)
	* [petsc-lite-3.7.0.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz)


### NIFTICLIB

* file: [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0/)
* description: library to read and write NIFTI images
* see [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) 
