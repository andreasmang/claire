# COLDREG


## Installation


### General Setup

Make sure the standard MPI wrappers for `mpicc` and `mpicxx` are available on your system (either by loading the right modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions). Add the following to your `~/.bashrc`:

```bash
export PATH=/path/to/mpicxx:/path/to/mpicc:${PATH}
export LD_LIBRARY_PATH=/path/to/mpi/lib/${LD_LIBRARY_PATH}
```


### Dependencies

COLDREG depends on the following libraries:

* FFTW (version 3.3.4)
* ACCFFT (needs FFTW)
* PETSc (version 3.7; needs BLAS and LAPACK)
* NIFTICLIB (version 2.0.0; needs zlib)

More information on how to add, install, and link these libraries, can be found in [external/README-LIBS.md](external/README-LIBS.md)


### Build COLDREG

To build the code via `make` do:

```bash
source libs/environment_vars.sh
make -j
```

Check the [makefile](makefile) before doing so:

* Set the `USEINTEL` flag to `1` if you use an **intel compiler** (`icc`). If you use a **GNU** compiler (`gcc`) set the `USEINTEL` flag to `0`.

* You can avoid `source libs/environment_vars.sh` by adding the entries in `libs/environment_vars.sh` to your `~/.bashrc`.


## Run COLDREG

To run the code using a test example do:
```bash
./bin/runcoldreg
```

For general options do:
```bash
./bin/runcoldreg -help
```

For more advanced options do:
```bash
./bin/runcoldreg -advanced
```
