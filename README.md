# COLDREG

**COLDREG** implements a parallel algorithm for **Constrained Large Deformation Diffeomorphic Image Registration**. Additional information on the methodology can be found in [doc/README-REFERENCES.md](doc/README-REFERENCES.md).


## Installation

The installation consists of three steps:

* compile and install the libraries/dependencies (see below and [doc/README-EXTLIBS.md](doc/README-EXTLIBS.md))
* set the environment variables to be able to link to the libraries (see [doc/README-EXTLIBS.md](doc/README-EXTLIBS.md))
* compile the code (see below and [doc/README-INSTALLATION.md](doc/README-INSTALLATION.md))


### Install Dependencies

COLDREG depends on the following libraries:

* [FFTW](http://www.fftw.org) (version 3.3.4)
* [ACCFFT](http://accfft.org) (requires `FFTW` and `libstdc++`)
* [PETSc](https://www.mcs.anl.gov/petsc/) (version 3.7; requires [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/))
* [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) (version 2.0.0; requires `zlib` and `libstdc++`)

More information on how to **add**, **install**, and **link** these libraries, can be found in [doc/README-EXTLIBS.md](doc/README-EXTLIBS.md)


### Install COLDREG

Instructions on how to install COLDREG can be found in [doc/README-INSTALLATION.md](doc/README-INSTALLATION.md). Instruction for running and installing the software on the TACC systems are provided in [doc/README-TACC.md](doc/README-TACC.md).


## Run COLDREG

To run COLDREG with a 32x32x32 test example do

```bash
./bin/runcoldreg
```

To run a simple test problem using the provided test images do:

```bash
./bin/runcoldreg -mr ./external/mR.nii.gz -mt ./external/mT.nii.gz -nx 256x256x256 -betav 1E-2 -regnorm h2s -xresults -x ./results
```

Here, `-mr ./external/mR.nii.gz` defines the *reference image*, `-mt ./external/mT.nii.gz` the *template image*, `-nx 256x256x256` the *size* of the images, `-betav 1E-2` the *regularization weight*,  `-regnorm h2s` the *regularization norm* (H2-seminorm in this case), `-x ./results` the *output folder*, and `-xresults` enables the write out of output images. To learn about the registration options you can do 

```bash
./bin/runcoldreg -help
```

For more advanced options do

```bash
./bin/runcoldreg -advanced
```

You can also find the options in [doc/help.txt](doc/help.txt) and [doc/advanced.txt](doc/advanced.txt). More instructions on how to run COLDREG can be found in [doc/README-RUNME.md](doc/README-RUNME.md).
