# Installing Dependencies

COLDREG depends on the following libraries:

* [FFTW](http://www.fftw.org) (version 3.3.4)
* [ACCFFT](http://accfft.org) (requires FFTW)
* [PETSc](https://www.mcs.anl.gov/petsc/) (version 3.7; requires [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/))
* [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) (version 2.0.0; requires `zlib` and `libstdc++6`)

These need to be installed and made available on your system before compiling the code (instruction for compiling COLDREG can be found in [doc/README-INSTALLATION.md](README-INSTALLATION.md)). Additional details for each individual library can be found below.


## Before Compiling

Make sure the standard MPI wrappers for `mpicc` and `mpicxx` are available on your system (either by loading the right modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions). For instance, add the following definitions to your `~/.bashrc`:

```bash
export PATH=/path/to/mpicxx:/path/to/mpicc:${PATH}
export LD_LIBRARY_PATH=/path/to/mpi/lib/${LD_LIBRARY_PATH}
```


## External Libraries/Dependencies

In general you should have received tarball files for the individual libraries. The **compressed** tarball files (`LIBRARYNAME.tar.gz`) should be located in or be added to the [external](../external) folder. If you use your local installations make sure you set the following variables in your `~/.bashrc`:

```bash
export FFTW_DIR=/path/to/fftw
export ACCFFT_DIR=/path/to/accfft
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=your_petsc_arch
export NIFTI_DIR=/path/to/nifticlib
```


## Building Dependencies

### General Info

I recommend that you add the path to the **FFTW** and the **PETSc** library to your `LD_LIBRARY_PATH`:

```bash
export LD_LIBRARY_PATH=/path/to/petsc/lib/${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/path/to/fftw/lib/${LD_LIBRARY_PATH}
```

If you are compiling the code with the shipped `build_libs.sh` (see below) this is eqivalent to

```bash
export LD_LIBRARY_PATH=/path/to/cold/external/libs/fftw/build/lib/${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/path/to/cold/external/libs/petsc/build/lib/${LD_LIBRARY_PATH}
```


### Quick Shot

If `mpicc` and `mpicxx` are available, you can install all external dependencies at once as follows:

```bash
cd external
./build_libs.sh --build
```


### More Detailed Instructions

The libraries can be compiled by running the [build_libs.sh](../external/build_libs.sh) script in the [external](../external) subdirectory. To see all the options do

```bash
./build_libs.sh --help
```

This will provide information on what parameters you can parse. Ideally it should be sufficient to do `./build_libs.sh --build`.  This will install all libraries in a local folder called "lib" in [external](../external/)). You can also build the individual libraries one after another. For instructions do `./build_libs.sh --help`.

If the wrapper for your MPI implementation does **not** provide `mpicc` and `mpicxx` you will have to pass the MPI compiler manually (**not tested**)

```bash
./build_libs.sh --cxx YOURMPICXXCOMPILER --c YOURMPICCOMPILER
```

An example for `mpicxx` and `mpicc` is

```bash
./build_libs.sh --cxx mpicxx --c mpicc
```


If you want to use your own BLAS and LAPACK installations, you can pass the path with the `--bldir` option: 

```bash
./build_libs.sh --bldir /path/to/libraries/
```

Please check the `cmake`, `make` and `automake` outputs for errors. To check if everything worked you can also inspect the "build" subfolders of the individual libraries in the "lib" folder (subdirectories of [external](../external)). See if folders in "build" were created and the library and include files exist.


### Adding Libraries to System

Before you are able to compile and run COLDREG you will have to add some **environment variables** to your system. When building the libraries a file called `environment_vars.sh` is created. This file should be located in [external/libs](../external/libs). To add the corresponding environment variables temporarily (for the current session) to your system, do

```bash
source environment_vars.sh
```

or copy the content of `environment_vars.sh` to your `~/.bashrc`.


## Compiling the Code

See [doc/README-INSTALLATION.md](README-INSTALLATION.md)



## Additional Info


### FFTW

* file: [fftw-3.3.4.tar.gz](ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz)
* description: library for computing FFTs
* see [FFTW](http://www.fftw.org) (version 3.3.4)


### ACCFFT

* file: [accfft.tar.gz](https://github.com/amirgholami/accfft)
* git: `git clone git@github.com:amirgholami/accfft.git`
* description: library to compute FFT in parallel (requires FFTW)
* see [ACCFFT](http://www.accfft.org)
* typical problems
	* `libstdc++6`:
		* if there are errors when compiling `ACCFFT` make sure that this library can be found on your system
		* they are typically located in your `gcc` library path (simply add the `lib(64)` folder to the `LD_LIBRARY_PATH`)
		* you might be able to locate these libraries via `locate libstdc`
	* `FFTW`: accfft depends on FFTW; info on FFTW can be found above


### PETSc

* file: [petsc-lite-3.7.0.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz)
* [bitbucket](https://bitbucket.org/petsc/petsc)
* description: library for numerics and optimization

### NIFTICLIB

* file: [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0/)
* description: library to read and write NIFTI images
* see [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) 
* typical problems
	* `zlib` and `libstdc++6`:
		* if there are errors when compiling `NIFTICLIB` make sure that these libraries can be found on your system
		* they are typically located in your `gcc` library path (simply add the `lib(64)` folder to the `LD_LIBRARY_PATH`)
		* you might be able to locate these libraries via `locate zlib` and `locate libstdc`
