# COLDREG



## Installation

### Dependencies

* FFTW (version 3.3.4)
* ACCFFT (needs FFTW)
* PETSc (version 3.7; needs BLAS and LAPACK)
* NIFTICLIB (version 2.0.0; needs zlib)

You can find these in the **external** folder.

### Build Code and Dependencies



#### General (Paths)

Make sure `mpicc` and `mpicxx` are in your path. Add the following to your `~/.bashrc`:

```bash
export PATH=/path/to/mpicxx:/path/to/mpicc:${PATH}
export LD_LIBRARY_PATH=/path/to/mpi/lib/${LD_LIBRARY_PATH}
```

I recommend to add the path to the FFTW library and the PETSc library to your `LD_LIBRARY_PATH` (lib folders for FFTW and PETSc). If you decide to use PETSc with your local MKL implementation, also add these to the `LD_LIBRARY_PATH`.



#### Build Dependencies

If `mpicc` and `mpicxx` are available, you can install all external dependencies at once as follows:

```bash
cd external
./build_libs.sh --build
source libs/environment_vars.sh
cd ..
make -j
```

### Run 

```bash
./bin/runcoldreg
```
