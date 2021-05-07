# CLAIRE: Installation and Requirements

Go back to [README.md](../README.md).

## Content

* [Installation Overview](#installation)
	* [One Shot](#oneshot)
	* [Step by Step](#stepbystep)
* [Detailed Installation Guide](#verboseinstall)
	* [Requirements](#requirements)
	* [Dependencies](#dependencies)
		* Required Dependencies and Compatibility
		* Step 1: Downloading and Installing Dependencies
		* Step 2: Setting Environment Variables
	* [Building CLAIRE](#buildclaire)
	* [Executing CLAIRE](#execclaire)
* [Additional Info for Dependencies](#depsinf)
* [CLAIRE on Specific Systems](#clairesys)
* [Troubleshooting and Known Issues](#faq)


## Installation Overview (Quick Guide)<a name="installation"></a>

In this section, we provide a minimal installation guide. We provide a make environment to download and install the dependencies using generic settings that have worked on most of our systems. If this brief installation guide does not work for you, please consult the [detailed installation guide](#verboseinstall) below.

### One Shot <a name="oneshot"></a>

To use the default settings to build dependencies and CLAIRE itself do the following:

```bash
cd deps
make
source env_source.sh
cd ..
make -j
./bin/claire -synthetic 0
```

The enviroment variables need to be sourced every time you log out of your computer or start a new bash (`source env_source.sh`). As an alternative, you can add the content of `env_source.sh` (for example) to your `.bashrc` or `.bash_profile`.

### Step by Step  <a name="stepbystep"></a>

Next, we go over the steps outlined above step by step. Again, more details are provided [below](#verboseinstall).

#### Step 1) Installing Dependencies

To install the dependencies (the PETSc and NIFTI libraries) go to the top level directory of CLAIRE in your command window and execute the following commands within your command window:

```bash
cd deps
make
```

This makefile downloads and compiles the dependencies for CLAIRE. To add these dependencies to your environment type the following into your command line and press return:

```bash
source env_source.sh
```

Notice that the enviroment variables need to be sourced every time you log out of your computer or start a new bash. As an alternative, you can add the content of `env_source.sh` to your `.bashrc` or `bash_profile`.


#### Step 2) Compiling CLAIRE

Assuming that you are in the top level directory of CLAIRE, all you need to do is to type

```bash
make -j
```

#### Step 3) Executing CLAIRE

If you would like to verify if CLAIRE has been installed correctly run the following command in your command window:

```bash
./bin/claire -synthetic 0
```

Additional examples for executing CLAIRE are described in [doc/README-RUNME.md](README-RUNME.md).


## Detailed Installation Guide <a name="verboseinstall"></a>

In this section we provide a more detailed description of the installation process to help users with troubleshooting. The following table lists system configurations on which we have successfully installed CLAIRE.

|Test   | Compiler  | MPI            | CUDA | PETSc  | CPU    | GPU   | System       |
|---    |---------- |-----           |------|------- |---     |---    |---           |
|b5213fa| GCC 9.3   | OpenMPI 4.0.3  | 11.0 | 3.14.2 | x86_64 | GA102 | Ubuntu 20.04 |
|6f40316| GCC 9.3   | OpenMPI 4.0.3  | 11.1 | 3.14.2 | x86_64 | GK110 | Ubuntu 20.04 |
|4967052| GCC 8.4   | OpenMPI 1.10.2 | 10.1 | 3.12.4 | x86_64 | GK110 | Ubuntu 16.04 |
|4967052| GCC 5.4.0 | OpenMPI 1.10.2 | 10.0 | 3.12.4 | x86_64 | GM200 | Ubuntu 16.04 |
|4967052| GCC 7.4   | OpenMPI 4.0.1  | 10.1 | 3.12.4 | x86_64 | GP100 | Ubuntu 16.04 |
|4967052| GCC 4.8.5 | OpenMPI 3.1.6  | 10.2 | 3.12.4 | Power9 | GV100 | CentOS 7.8   |
|4967052| XLC 16.1  | Spectrum 10.3  | 10.2 | 3.12.4 | Power9 | GV100 | RHEL 7.8     |


### Requirements <a name="requirements"></a>

The minimal requirements for compiling CLAIRE on your system are:
* MPI (Open MPI; MVAPICH; Intel MPI; ...; required by [PETSc](https://www.mcs.anl.gov/petsc), and CLAIRE)
* cmake (see [https://cmake.org](https://cmake.org); required by niftilib)
* python (see [https://www.python.org](https://www.python.org); required by [PETSc](https://www.mcs.anl.gov/petsc) and the optional pyclaire bindings)
* zlib (see [https://www.zlib.net](https://www.zlib.net); required by niftilib)
* CUDA-API

Make sure that the standard *wrappers* for `mpicc`, `mpicxx`, and `nvcc` are available on your system (either by loading the appropriate modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions). The compilation has been tested with Open MPI, MVAPICH, and Intel MPI.



### Dependencies <a name="dependencies"></a>


#### Required Dependencies and Compatibility

The compiler needs `C++11` support. The GPU version of CLAIRE requires the following libraries to be installed on your system:

* MPI (with GPU support (CUDA-aware MPI) for multi-GPU multi-node)
* PETSc with CUDA support (see [https://www.mcs.anl.gov/petsc](https://www.mcs.anl.gov/petsc))
* niftilib (see [https://sourceforge.net/projects/niftilib/files/nifticlib](https://sourceforge.net/projects/niftilib/files/nifticlib))
* zlib (see [http://zlib.net](http://zlib.net))

We provide functionality to build PETSc, niftilib, and zlip on your system (see next section).


#### Step 1: Downloading and Installing Dependencies

To download and compile the libraries we provide a `makefile` (see [deps/makefile](../deps/makefile)). Simply run `make` with this script in your command window to download *tarball* files of the libraries identified above.

```bash
cd deps
make
```

The *compressed* tarball files (i.e, `LIBRARY-NAME.tar.gz`) should remain located in or be added to the [deps](../deps) folder. Make sure that all libraries are downloaded (the progress bar of `wget` should be full). To view the urls for the libraries you can take a look at the [deps/makefile](../deps/makefile). We provide additional information about these libraries [below](#depsinf). This also includes links to versions for these libraries that we have used to compile the GPU version of CLAIRE before.

The [makefile](../deps/makefile) has some optional parameters to configure the build. The parameters can be set by `make PARAMETER=value`. Multiple parameter-value pairs can be passed to the make command. The [makefile](../deps/makefile) to compile the dependencies has the following parameters.

| PARAMETER       | Description                                           | Default | Valid Values  |
| --------------- | ----------------------------------------------------- | ------- | ------        |
| BUILD_PETSC     | PETSc version to download and compile; empty for none | 3.12.4  | PETSc Version |
| BUILD_NIFTI     | Download and build `niftilib`                         | yes     | yes, no       |
| WITH_BATCH      | Option to build petsc on a batch system, e.g. slurm   | no      | yes, no       |
| WITH_CUDA_MPI   | MPI is CUDA-aware                                     | yes     | yes, no       |
| CC              | Path to C compiler                                    | mpicc   | file path     |
| CXX             | Path to CXX compiler                                  | mpicxx  | file path     |
| NVCC            | Path to CUDA compiler                                 | nvcc    | file path     |
| WITH_PETSC_OPTS | additional PETSC compile options                      |         |               |

The libraries will be extracted and build in the `deps/lib` subfolder.


#### Step 2: Setting Environment Variables

Before you are able to compile and run CLAIRE you need to add *environment variables* to your system. When building the libraries a file called `env_source.sh` is created. This file should be located in the [debs](../deps) subfolder. To add the environment variables temporarily (for the current session) to your system, do

```bash
source env_source.sh
```

To add them permanently, copy the content of `env_source.sh` to your `~/.bashrc`. Notice that `env_source.sh` defines *absolute paths*.


## Building CLAIRE <a name="buildclaire"></a>

Before you can build CLAIRE you need to

* Make sure that you have installed all *dependencies* (see prior sections).
* Make sure all paths and compilers needed in the `makefile` are available on your system, i.e. `mpicxx`, `nvcc`, and the dependencies.

To inspect all options used in the `makefile` for CLAIRE (see [makefile](../makefile)) do (in the top level directory):

```bash
make VERBOSE=1 VVERBOSE=1 config
```

To build the code using the `make` system do (in the top level directory):

```bash
make -j
```

If you build in parallel using `make -j`, on certain systems to many threads will be used. This will result in compilation errors. To fix this, run `make -j 12` instead (for quick access, you may want to define an alias in your `~/.bashrc`).

The [makefile](../makefile) also contains optional parameters to configure the build. The parameters can be set by `make PARAMETER=value`. Multiple parameter-value pairs can be passed to the `make` command. The `makefile` to compile the dependencies has following parameters.

| PARAMETER      | Description                                           | Default | Valid Values  |
| -------------- | ----------------------------------------------------- | ------- | ------        |
| BUILD_TEST     | build the unit test applications                      | no      | yes; no       |
| BUILD_PYTHON   | build `pyclaire` python bindings                      | no      | yes; no       |
| BUILD_SHARED   | build CLAIRE as shared library                        | no      | yes; no       |
| WITH_NIFTI     | build with `niftilib`                                 | yes     | yes; no       |
| WITH_DEBUG     | build with additional debug informations              | no      | yes; no       |
| WITH_DEVELOP   | build CLAIRE additional development informations      | no      | yes; no       |
| WITH_CUDA_MPI  | MPI is CUDA-aware                                     | yes     | yes, no       |
| BUILD_TARGET   | target CPU architecture                               | X86     | X86; POWER9   |
| GPU_VERSION    | GPU CUDA version to compile, e.g. 35, 60, 70, 75      |         | Compute Capability |
| CPP_VERSION    | C++ Standard to use                                   | c++11   | c++11; c++14  |
| LD_FLAGS       | additional flags for the linker                       |         |               |
| CXX_FLAGS      | additional flags for the C++ compiler                 |         |               |
| NVCC_FLAGS     | additional flags for the CUDA compiler                |         |               |
| MPI_DIR        | main path to the MPI include and lib directory        |         |               |
| CUDA_DIR       | main path to the CUDA include and lib directory       |         |               |
| PETSC_DIR      | main path to the PETSc include and lib directory      |         |               |
| NIFTI_DIR      | main path to the libnifti include and lib directory   |         |               |
| ZLIB_DIR       | main path to the zlib include and lib directory       |         |               |
| PYTHON_DIR     | main path to the Python3 include and lib directory    |         |               |
| VERBOSE        | if set to any value the make command is verbose       |         |               |
| VVERBOSE       | if set to any value the make command is very verbose  |         |               |

Not that the `makefile` generates a cache (`make.cache`) to detect if a complete rebuild of CLAIRE is needed. If this file is removed or does not exsit the next build will first do `make clean` automatically.


## Executing CLAIRE <a name="execclaire"></a>

If you would like to verify if CLAIRE has been installed correctly, run the following command in your command window:

```bash
./bin/claire -synthetic 0
```

Additional examples for executing CLAIRE are described in [doc/README-RUNME.md](README-RUNME.md)


## Additional Info for Dependencies <a name="depsinf"></a>

### PETSc
* PETSc webpage: [https://www.mcs.anl.gov/petsc](https://www.mcs.anl.gov/petsc)
* description: library for numerics, linear algebra, and optimization
* source code also available on bitbucket: [https://bitbucket.org/petsc/petsc](https://bitbucket.org/petsc/petsc)
* versions that have been succesfully used by our group:
	* [petsc-lite-3.11.4.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.4.tar.gz)

### nifticlib
* NIFTICLIB webpage: [http://niftilib.sourceforge.net](http://niftilib.sourceforge.net)
* description: library to read and write NIFTI images
* see [https://sourceforge.net/projects/niftilib/files/nifticlib](https://sourceforge.net/projects/niftilib/files/nifticlib)
* versions that have been succesfully used by our group:
	* [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0)



## CLAIRE on Specific Systems <a name="clairesys"></a>

### TACC's Longhorn System (03/17/21)

More information about TACC's Longhorn system can be found at [https://www.tacc.utexas.edu/systems/longhorn](https://www.tacc.utexas.edu/systems/longhorn).

Modules loaded:
```bash
1) xl/16.1.1             4) autotools/1.2   7) TACC
2) spectrum_mpi/10.3.0   5) cmake/3.16.1    8) cuda/10.2 (g)
3) git/2.24.1            6) xalt/2.10.2
```

Compilation of CLAIRE and its dependencies:

```bash
cd deps
make WITH_BATCH=yes
source deps/env_source.sh
make BUILD_TARGET=POWER9 GPU_VERSION=70
```

To test if the compilation worked check if the binaries are available in the `bin` folder. To execute CLAIRE using an interactive job do

```bash
cd bin
idev -N1 # launch an interactive session with one node
ibrun ./claire -help
ibrun ./claire -synthetic 0 -nx 128
```

A job submission file for TACC's Longhorn system (for multi-GPU exection) can be found in [doc/examples/longhorn_mgpu.slurm](examples/longhorn_mgpu.slurm).



## Troubleshooting / Known Issues <a name="faq"></a>

* if MPI is not compiled with CUDA-aware options, add the file `.petscrc` to the working directory and add the option `-use_gpu_aware_mpi 0`
* CUDA >= 11.0 is only supported with PETSc >= 3.14.
* Kepler GPUs work with PETSc 3.12.4  (others not tested)
* Compiling PETSc with CUDA support on cluster login nodes without GPUs might fail
* PNETCDF is currently not tested for GPUs
* The GPU version of CLAIRE can currently only be compiled in single precision. This limits the selection of regularization operators to H1-type regularization only. There are issues with the numerical accuracy of H2- and H3-type regularization operators for single precision. Applying these operators requires a compilation in double precision (available on the GPU branch) 
