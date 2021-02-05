# CLAIRE: Installation and Requirements

Go back to [README.md](../README.md).

## Content
* [Installation Overview](#installation)
* [Detailed Installation Guide](#verboseinstall)
	* [Requirements](#requirements)
	* [Dependencies](#dependencies)
	* [Building CLAIRE](#buildclaire)
* [Instructions for Specific Systems](#systems)
* [Trouble Shooting](#faq)

## Installation Overview <a name="installation"></a>

What follows is a minimal installation guide. We provide scripts to download and install the dependencies using generic settings that have worked on our systems. A more detailed installation guide can be found [below](#verboseinstall).

To install the dependencies go to the top level directory of CLAIRE in your command window and execute the following commands in your command window:
```bash
cd deps
./get_libs.sh
./build_libs.sh --build
```

The first bash script downloads the libraries you need to install. The second bash script compiles these libraries using default settings. To add these dependencies to your environment type the following into your command line and press return:
```bash
source environment_vars.sh
```

This step needs to be done every time you log out of your computer. As an alternative, you can add the content of `environment_vars.sh` to your `.bashrc` or `bash_profile`.

Return to the top level directory of CLAIRE and type
```bash
cd ..
make -j
```

If you would like to verify if CLAIRE has been installed correctly run the following command in your command window:

```bash
./bin/claire -synthetic 0
```

If you have trouble executing CLAIRE take a look at the [examples and tests](README-RUNME.md) described in [README-RUNME.md](README-RUNME.md).

## Detailed Installation Guide <a name="verboseinstall"></a>

### Requirements <a name="requirements"></a>

* MPI (Open MPI; MVAPICH; Intel MPI); required by `AccFFT`, `PETSc`, and `CLAIRE`
* cmake ([https://cmake.org](https://cmake.org)); required by `AccFFT` and `niftilib`
* python ([https://www.python.org](https://www.python.org)); version 2.7; required by `PETSc`

Make sure that the standard *wrappers* for `mpicc` and `mpicxx` are available on your system (either by loading the appropriate modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions below). The compilation has been tested with `Open MPI`, `MVAPICH`, and `Intel MPI`.

### Dependencies <a name="dependencies"></a>

CLAIRE requires the following libraries to be installed on your system:
* `FFTW` [http://www.fftw.org](http://www.fftw.org)
	* file: [fftw-3.3.6-pl2.tar.gz](http://www.fftw.org/fftw-3.3.6-pl2.tar.gz)
	* description: library for computing FFTs
	* additional details can be found at [http://www.fftw.org](http://www.fftw.org)
* `AccFFT` [http://accfft.org](http://accfft.org)
	* file: accfft.tar.gz
	* source code also available on gitub: [https://github.com/amirgholami/accfft](https://github.com/amirgholami/accfft)
	* description: library to compute FFT in parallel (requires FFTW)
	* additional details can be found at [http://www.accfft.org](http://www.accfft.org)
	* requires `FFTW` to be installed
* `PETSc` [https://www.mcs.anl.gov/petsc/](https://www.mcs.anl.gov/petsc/) (requires `python 2.7`)
	* file: [petsc-lite-3.9.1.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.9.1.tar.gz)
	* source code also available on bitbucket: [https://bitbucket.org/petsc/petsc](https://bitbucket.org/petsc/petsc)
	* description: library for numerics, linear algebra, and optimization
	* older versions that have been succesfully used by our group:
		* [petsc-lite-3.8.3.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.8.3.tar.gz)
		* [petsc-lite-3.7.6.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz)
		* [petsc-lite-3.7.0.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz)
* `niftilib` [https://sourceforge.net/projects/niftilib/files/nifticlib/](https://sourceforge.net/projects/niftilib/files/nifticlib/)
	* file: [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0)
	* description: library to read and write NIFTI images
	* see [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/)
* `zlib` [http://zlib.net](http://zlib.net)
* `libmorton` [https://github.com/Forceflow/libmorton](https://github.com/Forceflow/libmorton)

#### Step 1: Downloading Dependencies
To download the libraries we provide a script called `get_libs.sh` (see [deps/get_libs.sh](../deps/get_libs.sh)). Simply run this script in your command window to download *tarball* files of the libraries identified above.

```bash
cd debs
./get_libs.sh
```

The *compressed* tarball files (i.e, `LIBRARY-NAME.tar.gz`) should remain located in or be added to the [deps](../deps) folder. Make sure that all libraries are downloaded (the progress bar of `wget` should be full). If the download failed, delete the tarball file and run `get_libs.sh` again. To view the urls for the libraries you can take a look at [deps/get_libs.sh](https://github.com/andreasmang/claire/tree/master/deps/get_libs.sh). We provide additional information about these libraries [below](#depsinf). This also includes links to versions for these libraries that we have used with CLAIRE before.


#### Step 2: Installing Dependencies

##### Quick Shot
To build all libraries at once, execute the build script in your command window as follows:

```bash
cd deps
./build_libs.sh --build
```

This will decompress the *tarball* files and compile them with standard settings. The compiled libraries will be installed in in a subfolder of your `deps` folder called `libs`.


##### Detailed Instructions
We build all libraries as **static** by default.

The libraries can be compiled by running the [build_libs.sh](https://github.com/andreasmang/claire/tree/master/deps/build_libs.sh) script in the [external](https://github.com/andreasmang/claire/tree/master/deps) subdirectory. To see all the options do

```bash
./build_libs.sh --help
```

This will provide information on what parameters you can parse. Ideally, it is sufficient to do `./build_libs.sh --build`. You can also build the individual libraries one after another, via the `--bLIBNAME` option, where `LIBNAME` is the name of the library. For precise instructions, do `./build_libs.sh --help`. If you want to clean up the libraries folder, you can do `./build_libs.sh --clean`. The *build* folders will be removed each time you recompile the libraries.

Please check the `cmake`, `make` and `automake` outputs for errors. To check if everything worked you can also take a look at the "build" subdirectories of the individual libraries in the "lib" folder (subdirectories of [deps](../deps)). See if folders in "build" were created and the library and include files exist. If not, something went wrong during the build.


#### Step 3: Setting Environment Variables
Before you are able to compile and run CLAIRE you need to add *environment variables* to your system. When building the libraries a file called `environment_vars.sh` is created. This file should be located in the [debs/libs](https://github.com/andreasmang/claire/tree/master/deps) subfolder. To add the environment variables temporarily (for the current session) to your system, do

```bash
source environment_vars.sh
```

To add them permanently, copy the content of `environment_vars.sh` to your `~/.bashrc`. Notice that `environment_vars.sh` defines *absolute paths*.

### Building CLAIRE <a name="buildclaire"></a>

Before you can build CLAIRE you need to

* Make sure that you have installed all *dependencies*.
* Check the [makefile](https://github.com/andreasmang/claire/tree/master/makefile) before building the code:
	* If you use an *intel compiler* set the `USEINTEL` flag to `yes`.
	* If you use a *GNU compiler* set the `USEINTEL` flag to `no`.
	* If you use *Intel MPI* (impi) set the `USEINTELMPI` flag to `yes` (if not, set it to `no`).
* Make sure all paths needed in the makefile are available on your system (to check, you can do `env` in your bash). to add the paths temporarily `source external/libs/environment_vars.sh` or add the content of `external/libs/environment_vars.sh` to your `~/.bashrc`.

To build the code using the `make` system do (in the top level directory):

```bash
make -j
```

If you build in parallel using `make -j`, on certain systems to many threads will be used. This will result in compilation errors. To fix this, run `make -j 12` instead (for quick access, you may want to define an alias in your `~/.bashrc`).


## Instructions for Specific Systems <a name="systems"></a>

Information about the Texas Advanced Computing Center ([TACC](https://www.tacc.utexas.edu)) can be found at [https://www.tacc.utexas.edu](https://www.tacc.utexas.edu).In the following we will mainly provide information on which modules you need to load on the individual TACC systems to compile the code and the libraries. For more information on the systems please visit the [TACC](https://www.tacc.utexas.edu) web pages.


## Maverick Systems

Info about Maverick can be found here: [https://www.tacc.utexas.edu/systems/maverick](https://www.tacc.utexas.edu/systems/maverick).

CLAIRE has been tested with the following module settings:

```bash
module load intel/15.0.3
module load impi/5.0.3
module load git/2.7.0
module load cmake/2.8.12.2
```



## Stampede

```bash
module load intel/14.0.1.106
module load impi/5.0.3 
```


## Troubleshooting <a name="faq"></a>

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
5. ./include/RegOpt.hpp: fatal error: petsc.h: No such file or directory
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
8. Other dependencies that may cause problems (should in general be available on your system)
	* cmake (https://cmake.org; required by ACCFFT and niftilib; version 2.8 or greater)
	* openMP (required by ACCFFT)
	* `crypt` library (required by PETSc)
	* `ssl` library (required by PETSc)
	* BLAS (required by PETSc; we install it along with PETSc; http://www.netlib.org/blas/)
	* LAPACK (required by PETSc; we install it along with PETSc; http://www.netlib.org/lapack/)
