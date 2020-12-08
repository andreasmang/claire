# CLAIRE: Detailed Installation Guide


## Content

* [Requirements](#requirements)
* [Dependencies](#dependencies)
	* Required Dependencies
  * Step 1: Downloading Dependencies
  * Step 2: Installing Dependencies
  	* Quick Shot
  	* Detailed Instructions
  * Step 3: Setting Environment Variables
* [Building CLAIRE](#buildclaire)
* [Additional Info for Dependencies](#depsinf)
* [Troubleshooting](#faq)


## Requirements <a name="requirements"></a>

* MPI (Open MPI; MVAPICH; Intel MPI); required by `PETSc`, and `CLAIRE`
* cmake ([https://cmake.org](https://cmake.org)); required by `niftilib`
* python ([https://www.python.org](https://www.python.org)); required by `PETSc` and the optional `pyclaire` bindings
* zlib ([https://www.zlib.net/](https://www.zlib.net/)); requiered by `niftilib`
* CUDA-API


Make sure that the standard *wrappers* for `mpicc`, `mpicxx`, and `nvcc` are available on your system (either by loading the appropriate modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions below). The compilation has been tested with `Open MPI`, `MVAPICH`, and `Intel MPI`.

## Dependencies <a name="dependencies"></a>

### Required Dependencies
CLAIRE-GPU requires the following libraries to be installed on your system:

* `PETSc` with CUDA support [https://www.mcs.anl.gov/petsc/](https://www.mcs.anl.gov/petsc/)
* `niftilib` [https://sourceforge.net/projects/niftilib/files/nifticlib/](https://sourceforge.net/projects/niftilib/files/nifticlib/)
* `zlib` [http://zlib.net](http://zlib.net)

### Step 1: Downloading and Installing Dependencies

To download and compile the libraries we provide a makefile (see [deps/makefile](../deps/makefile)). Simply run `make` with this script in your command window to download *tarball* files of the libraries identified above.

```bash
cd deps
make
```

The *compressed* tarball files (i.e, `LIBRARY-NAME.tar.gz`) should remain located in or be added to the [deps](../deps) folder. Make sure that all libraries are downloaded (the progress bar of `wget` should be full). To view the urls for the libraries you can take a look at the [deps/makefile](../deps/makefile). We provide additional information about these libraries [below](#depsinf). This also includes links to versions for these libraries that we have used with CLAIRE before.

The makefile has some optional parameters to configure the build. The parameters can be set by `make PARAMETER=value`. Multiple parameter-value pairs can be passed to the make command. The makefile to compile the dependencies has following parameters.

| PARAMETER   | Description                                           | Default | Valid Values  |
| ----------- | ----------------------------------------------------- | ------- | ------        |
| BUILD_PETSC | PETSc version to download and compile; empty for none | 3.12.4  | PETSc Version |
| BUILD_NIFTI | Download and build `niftilib`                         | yes     | yes, no       |
| WITH_BATCH  | Option to build petsc on a batch system, e.g. slurm   | no      | yes, no       |
| CC          | Path to C compiler                                    | mpicc   | file path     |
| CXX         | Path to CXX compiler                                  | mpicxx  | file path     |
| NVCC        | Path to CUDA compiler                                 | nvcc    | file path     |

The libraries will be extraxted and build in the `deps/lib` subfolder.

### Step 2: Setting Environment Variables
Before you are able to compile and run CLAIRE you need to add *environment variables* to your system. When building the libraries a file called `env_source.sh` is created. This file should be located in the [debs](../deps) subfolder. To add the environment variables temporarily (for the current session) to your system, do

```bash
source env_source.sh
```

To add them permanently, copy the content of `env_source.sh` to your `~/.bashrc`. Notice that `env_source.sh` defines *absolute paths*.

## Building CLAIRE for CPU <a name="buildclaire"></a>

Before you can build CLAIRE you need to

* Make sure that you have installed all *dependencies*.
* Make sure all paths and compilers needed in the makefile are available on your system, i.e. mpicxx, nvcc, and the dependencies.

To first see the all options used by the makefile do (in the top level directory):

```bash
make VERBOSE=1 VVERBOSE=1 config
```

To build the code using the `make` system do (in the top level directory):

```bash
make -j
```

If you build in parallel using `make -j`, on certain systems to many threads will be used. This will result in compilation errors. To fix this, run `make -j 12` instead (for quick access, you may want to define an alias in your `~/.bashrc`).

The makefile has some optional parameters to configure the build. The parameters can be set by `make PARAMETER=value`. Multiple parameter-value pairs can be passed to the make command. The makefile to compile the dependencies has following parameters.

| PARAMETER    | Description                                           | Default | Valid Values  |
| -----------  | ----------------------------------------------------- | ------- | ------        |
| BUILD_TEST   | build the unit test applications                      | no      | yes; no       |
| BUILD_PYTHON | build `pyclaire` python bindings                      | no      | yes; no       |
| BUILD_SHARED | build CLAIRE as shared library                        | no      | yes; no       |
| WITH_NIFTI   | build with `niftilib`                                 | yes     | yes; no       |
| WITH_DEBUG   | build with additional debug informations              | no      | yes; no       |
| WITH_DEVELOP | build CLAIRE additional development informations      | no      | yes; no       |
| BUILD_TARGET | target CPU architecture                               | X86     | X86; POWER9   |
| GPU_VERSION  | GPU CUDA version to compile, e.g. 35, 60, 70, 75      |         | Compute Capability |
| CPP_VERSION  | C++ Standard to use                                   | c++11   | c++11; c++14  |
| LD_FLAGS     | additional flags for the linker                       |         |               |
| CXX_FLAGS    | additional flags for the C++ compiler                 |         |               |
| NVCC_FLAGS   | additional flags for the CUDA compiler                |         |               |
| MPI_DIR      | main path to the MPI include and lib directory        |         |               |
| CUDA_DIR     | main path to the CUDA include and lib directory       |         |               |
| PETSC_DIR    | main path to the PETSc include and lib directory      |         |               |
| NIFTI_DIR    | main path to the libnifti include and lib directory   |         |               |
| ZLIB_DIR     | main path to the zlib include and lib directory       |         |               |
| PYTHON_DIR   | main path to the Python3 include and lib directory    |         |               |
| VERBOSE      | if set to any value the make command is verbose       |         |               |
| VVERBOSE     | if set to any value the make command is very verbose  |         |               |

Not that the makefile generates a cache (`make.cache`) to detect if a complete rebuild of CLAIRE is needed. If this file is removed or not exsiting the next build will first do `make clean` automatically.


## Additional Info for Dependencies <a name="depsinf"></a>

### PETSc

* file: [petsc-lite-3.12.4.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.12.4.tar.gz)
* source code also available on bitbucket: [https://bitbucket.org/petsc/petsc](https://bitbucket.org/petsc/petsc)
* description: library for numerics, linear algebra, and optimization
* older versions that have been succesfully used by our group:
	* [petsc-lite-3.11.4.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.4.tar.gz)

### nifticlib

* file: [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0)
* description: library to read and write NIFTI images
* see [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/)


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
