# COLDREG

**COLDREG** implements a parallel solver for **Constrained Large Deformation Diffeomorphic Image Registration**. Additional information on the methodology can be found in [doc/README-REFERENCES.md](doc/README-REFERENCES.md).

## Content

* Installation
* Running COLDREG
* Advanced Instructions
* License


If there are any issues or you have any questions send me an email: <andreas@ices.utexas.edu>.

## Installation

The installation consists of three steps:

* **Step 1**: compilation and installation of the libraries/dependencies
* **Step 2**: setting the environment variables to be able to link to the libraries
* **Step 3**: compilation of the code

Instructions for these steps can be found below.


### Before you Begin

* make sure **cmake** is available
* make sure **python** is available
* make sure wrappers for **mpicc** and **mpicxx** are available (the code has been tested with *Open MPI*, *MVAPICH*, and *Intel MPI* on Linux systems)
* make sure **OpenMP** is available on your system


### STEP 1: Installation of Dependencies

COLDREG depends on the following libraries:

* [FFTW](http://www.fftw.org) (version 3.3.4)
* [ACCFFT](http://accfft.org) (requires FFTW, OpenMP, and [cmake](https://cmake.org))
* [PETSc](https://www.mcs.anl.gov/petsc/) (version 3.7; requires python 2.7)
* [NIFTICLIB](https://sourceforge.net/projects/niftilib/files/nifticlib/) (version 2.0.0; requires [cmake](https://cmake.org))

The *tarball files* for these libraries are in the *external* subfolder. I am assuming you are in the *top level directory* of the code. If you install for the first time, it is recommended that you build the libraries one by one. First change to the *external folder*:

```bash
cd external
```

Then build *FFTW*, *ACCFFT*, *PETSc*, and *NIFTICLIB* do doing the following **line by line**:

```bash
./build_libs.sh --bfftw
./build_libs.sh --baccfft
./build_libs.sh --bpetsc
./build_libs.sh --bnifti
```

Please check the *cmake*, *make* and *automake* outputs for errors. To check if everything worked, you can also take a look at the "build" subfolders of the individual libraries in the "libs" subdirectory (subdirectories of [external](external)). See if folders in "build" were created and the libraries and include files exist. If there are problems, consult [doc/README-INSTALL.md](doc/README-INSTALL.md). 

To **build** all **dependencies at once** run the *build_libs.sh* in the *external* subfolder:

```bash
cd external
./build_libs.sh --build
```


### STEP 2: Set Environment Variables

If you are still in the *external* subfolder, do:

```bash
source libs/environment_vars.sh
```

If you are in the top level directory of your code, do:

```bash
source external/libs/environment_vars.sh
```

### STEP 3: Compile COLDREG

If you are using an *intel compiler* set `USEINTEL` in the [makefile](makefile) to `yes`; if not, set it to `no`. If you are using *IntelMPI* set `USEINTELMPI` in the makefile to `yes`; if not, set it to `no`.

In the *top level directory* of the code, do

```bash
make -j
```


## Runing COLDREG

To run COLDREG with a 32x32x32 test example do

```bash
./bin/runcoldreg
```

To **run** a simple test problem using some test images do:

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

You can also find a list of the available options for the binary in [doc/help.txt](doc/help.txt) and [doc/advanced-help.txt](doc/advanced-help.txt).


## Advanced Instructions

More information on how to **add**, **install**, and **link** these libraries, can be found in [doc/README-INSTALL.md](doc/README-INSTALL.md).


## License

Read the [license](LICENSE) file for more details.