# COLDREG

**COLDREG** implements a parallel solver for **Constrained Large Deformation Diffeomorphic Image Registration**. Additional information on the methodology can be found in [doc/README-REFERENCES.md](doc/README-REFERENCES.md).

* [Installation](##installation)
* [Run](#run)
* [Advanced Instructions](#advanced-install)
* [License](#license)


If there are any issues or you have any questions do not hesitate to send an email to <andreas@ices.utexas.edu>.

## Installation

The installation consists of three steps:

* **Step 1**: compilation and installation of the libraries/dependencies
* **Step 2**: setting the environment variables to be able to link to the libraries
* **Step 3**: compilation of the code

Instructions for these four steps can be found below.


### Before you Begin

* make sure **cmake** is available
* make sure **python** is available
* make sure wrappers for **mpicc** and **mpicxx** are available


### STEP 1: Installation of Dependencies

COLDREG depends on the following libraries:

* [FFTW](http://www.fftw.org) (version 3.3.4)
* [ACCFFT](http://accfft.org) (requires FFTW, OpenMP, and [cmake](https://cmake.org))
* [PETSc](https://www.mcs.anl.gov/petsc/) (version 3.7; requires python 2.7)
* [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) (version 2.0.0; requires [cmake](https://cmake.org))

The *tarball files* for these libraries are in the *external* subfolder. Assuming you are in the *top level directory* of the code: To **build** all **dependencies at once** run the *build_libs.sh* in the *external* subfolder:

```bash
cd external
./build_libs.sh --build
```


### STEP 2: Set Environment Variables

If you are still in the *external* subfolder do:

```bash
source libs/environment_vars.sh
```


### STEP 3: Compile COLDREG

If your still in the *external* subfolder, change back to the top level directory of the code (i.e., `cd ..`). Then do

```bash
make -j
```

If you are using an *intel compiler* set `USEINTEL` in the makefile to `yes`; if not, set it to `no`. If you are using *IntelMPI* set `USEINTELMPI` in the makefile to `yes`; if not, set it to `no`.



##<a name="run"></a> Run COLDREG

To run COLDREG with a 32x32x32 test example do

```bash
./bin/runcoldreg
```

To **run** a simple test problem using the provided test images do:

```bash
./bin/runcoldreg -mr ./external/mR.nii.gz -mt ./external/mT.nii.gz -nx 256x256x256 -betav 1E-2 -regnorm h2s -xresults -x ./results
```

Here, `-mr ./external/mR.nii.gz` defines the *reference image* (fixed image), `-mt ./external/mT.nii.gz` the *template image* (image to be registered), `-nx 256x256x256` the *size* of the images, `-betav 1E-2` the *regularization weight*,  `-regnorm h2s` the *regularization norm* (H2-seminorm in this case), `-x ./results` the *output folder*, and `-xresults` enables the output of images, the computed velocity field, the deformation map, and derived measures.

To learn about the **options** you can do

```bash
./bin/runcoldreg -help
```

For more advanced options do

```bash
./bin/runcoldreg -advanced
```


##<a name="advanced-install"></a>  Advanced Instructions

More information on how to **add**, **install**, and **link** these libraries, can be found in [doc/README-INSTALL.md](doc/README-INSTALL.md). You can find a list of the available options for the binary in [doc/help.txt](doc/help.txt) and [doc/advanced-help.txt](doc/advanced-help.txt). More instructions on how to run COLDREG can be found in [doc/README-RUNME.md](doc/README-RUNME.md).

##<a name="license"></a>  License

Read the [License](LICENSE) for more details.