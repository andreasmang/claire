# XXX



## Installation

### Dependencies
The code depends on the following libraries:

* FFTW
* ACCFFT
* PETSc
* NIFTICLIB


### Build Code and Dependencies

Make sure `mpicc` and `mpicxx` are in your path. Add the following to your `~/.bashrc`:

```bash
export PATH=/path/to/mpicxx:/path/to/mpicc:${PATH}
export LD_LIBRARY_PATH=/path/to/mpi/lib/${LD_LIBRARY_PATH}
```

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
