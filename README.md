# COLDREG

**COLDREG** implements a parallel algorithm for **Constrained Large Deformation Diffeomorphic Image Registration**. Additional information on the methodology can be found in [doc/README-REFERENCES.md](doc/README-REFERENCES.md).


## Installation

### COLDREG

Instructions on how to install COLDREG can be found in [doc/README-INSTALLATION.md](doc/README-INSTALLATION.md).


### Dependencies

COLDREG depends on the following libraries:

* [FFTW](http://www.fftw.org) (version 3.3.4)
* [ACCFFT](http://accfft.org) (requires FFTW)
* [PETSc](https://www.mcs.anl.gov/petsc/) (version 3.7; requires [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/))
* [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) (version 2.0.0; requires zlib)

More information on how to add, install, and link these libraries, can be found in [doc/README-EXTLIBS.md](doc/README-EXTLIBS.md)


## Run COLDREG

Instructions on how to use COLDREG can be found in [doc/README-INSTALLATION.md](doc/README-INSTALLATION.md). To check if your installation works you can do:

```bash
./bin/runcoldreg
```

This will run a 32x32x32 registration test problem with default parameters.
