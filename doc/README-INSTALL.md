# Installing and Running COLDREG


## Before Compiling

Make sure that the standard **MPI wrappers** for `mpicc` and `mpicxx` are available on your system (either by loading the appropriate modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions). Add the following to your `~/.bashrc`:

```bash
export PATH=/path/to/mpicxx:/path/to/mpicc:${PATH}
export LD_LIBRARY_PATH=/path/to/mpi/lib:${LD_LIBRARY_PATH}
```

To compile ACCFFT and NIFTICLIB you need to make sure that `libstdc++` and `zlib` are available on your system. To check you can try to run `locate libstdc++` and `locate zlib` in your bash.


## Installing Dependencies

COLDREG depends on the following libraries:

* [FFTW](http://www.fftw.org) (version 3.3.4)
* [ACCFFT](http://accfft.org) (requires *FFTW*, `libstdc++`, *OpenMP* and [cmake](https://cmake.org))
* [PETSc](https://www.mcs.anl.gov/petsc/) (version 3.7; requires *python 2.7* ([https://www.python.org](https://www.python.org)), [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/))
* [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) (version 2.0.0; requires [cmake](https://cmake.org), `zlib` and `libstdc++`)

These libraries to be installed and made available on your system before compiling the code. We build all libraries as **static** by default.

### External Libraries/Dependencies

In general, you should have received tarball files for the individual libraries. The *compressed* tarball files (i.e, *LIBRARYNAME.tar.gz*) should be located in or be added to the [external](../external) folder. If you decide to use your local installations (*not recommended*) make sure you set the following variables in your `~/.bashrc`:

```bash
export FFTW_DIR=/path/to/fftw
export ACCFFT_DIR=/path/to/accfft
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=your_petsc_arch
export NIFTI_DIR=/path/to/nifticlib
```


### Building Dependencies


#### Quick Shot

If `mpicc` and `mpicxx` are available, you can install all external dependencies at once as follows:

```bash
cd external
./build_libs.sh --build
```
This is not recommended, if you install the libraries for the first time. If you use **IntelMPI** add the `--useimpi` option.


#### More Detailed Instructions

The libraries can be compiled by running the [build_libs.sh](../external/build_libs.sh) script in the [external](../external) subdirectory. To see all the options do

```bash
./build_libs.sh --help
```

This will provide information on what parameters you can parse. Ideally it should be sufficient to do `./build_libs.sh --build`.  You can also build the individual libraries one after another. For instructions do `./build_libs.sh --help`. If you want to reconfigure and recompile the libraries you can do `./build_libs.sh --clean`. The *build* folders will be removed each time you recompile the libraries. 

If your MPI implementation does **not** provide wrappers for `mpicc` and `mpicxx`, you will have to pass the MPI compiler manually (**not tested** and probably fails)

```bash
./build_libs.sh --cxx YOURMPICXXCOMPILER --c YOURMPICCOMPILER
```

An example for `mpicxx` and `mpicc` is

```bash
./build_libs.sh --cxx mpicxx --c mpicc
```

Please check the `cmake`, `make` and `automake` outputs for errors. To check if everything worked you can also take a look at the "build" subfolders of the individual libraries in the "lib" folder (subdirectories of [external](../external)). See if folders in "build" were created and the library and include files exist.


#### Adding Libraries to System

Before you are able to compile and run COLDREG you will have to add some **environment variables** to your system. When building the libraries a file called `environment_vars.sh` is created. This file should be located in [external/libs](../external/libs). To add the corresponding environment variables temporarily (for the current session) to your system, do

```bash
source environment_vars.sh
```

or copy the content of `environment_vars.sh` to your `~/.bashrc`. Note that the script will set **absolute paths**. 


#### Additional Info


##### FFTW

* file: [fftw-3.3.4.tar.gz](ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz)
* description: library for computing FFTs
* see [FFTW](http://www.fftw.org) (version 3.3.4)


##### ACCFFT

* file: [accfft.tar.gz](https://github.com/amirgholami/accfft)
* git: `git clone git@github.com:amirgholami/accfft.git`
* description: library to compute FFT in parallel (requires FFTW)
* see [ACCFFT](http://www.accfft.org)
* typical problems
	* cmake needs to be installed [cmake](https://cmake.org)
	* `libstdc++6`:
		* if there are errors when compiling `ACCFFT` make sure that this library can be found on your system
		* they are typically located in your `gcc` library path (simply add the `lib(64)` folder to the `LD_LIBRARY_PATH`)
		* you might be able to locate these libraries via `locate libstdc`
	* `FFTW`: accfft depends on FFTW; info on FFTW can be found above


##### PETSc

* file: [petsc-lite-3.7.0.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.0.tar.gz)
* [bitbucket](https://bitbucket.org/petsc/petsc)
* description: library for numerics and optimization
* typical problems:
	* requires at least `python 2.7`; afaik `python 3` is not supported


##### NIFTICLIB

* file: [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0/)
* description: library to read and write NIFTI images
* see [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) 
* typical problems
	* cmake needs to be installed [cmake](https://cmake.org)
	* `zlib` and `libstdc++6`:
		* if there are errors when compiling `NIFTICLIB` make sure that these libraries can be found on your system
		* they are typically located in your `gcc` library path (simply add the `lib(64)` folder to the `LD_LIBRARY_PATH`)
		* you might be able to locate these libraries via `locate zlib` and `locate libstdc`


## Building COLDREG

Before you can build COLDREG you need to 

* make sure that you have installed the **external dependencies** (visit [doc/README-EXTLIBS.md](README-EXTLIBS.md) to learn more)
* check the [makefile](makefile) before building the code:
	* if you use an **intel compiler** (`icc`) set the `USEINTEL` flag to `yes`
	* if you use a **GNU compiler** (`gcc`) set the `USEINTEL` flag to `no`
	* if you use a **IntelMPI** (`impi`) set the `USEINTELMPI` flag to `yes` (if not, set it to `no`)
* make sure all paths needed in the makefile are available on your system (to check, you can do `env` in your bash); to add the pathes necessary to find and link against the library you can `source libs/environment_vars.sh` or add the content of `libs/environment_vars.sh` to your `~/.bashrc`

To build the code using the `make` system do (in the top level directory):

```bash
make -j
```


## Running COLDREG

If everything compiled correctly, you can run a test example by doing:

```bash
./bin/runcoldreg
```

If you get an error message that indicates that the PETSc library could not be found, you probably forgot to

```bash
source libs/environment_vars.sh
```

To get a general idea on how to run the binary do: 

```bash
./bin/runcoldreg -help
```

For more advanced options do:

```bash
./bin/runcoldreg -advanced
```

You can also find a list of the available options for the binary in [help.txt](help.txt) and [doc/advanced-help.txt](advanced-help.txt).