# COLDREG

**COLDREG** implements a parallel solver for **Constrained Large Deformation Diffeomorphic Image Registration**. Additional information on the methodology can be found in [doc/README-REFERENCES.md](doc/README-REFERENCES.md).

## Content

* Installation
* Running COLDREG
* Advanced Instructions
* License

If there are any issues or you have any questions send me an email: <andreas@ices.utexas.edu>.

## Installation

### Requirements

* **cmake** [cmake](https://cmake.org)
* **python** version 2.7
* **mpicc**, **mpicxx** should be in the path and should support openMP

### INSTALL

```bash
cd external
./build_libs.sh --build
source libs/environment_vars.sh
cd ..
```

If you are using an *intel compiler* set `USEINTEL` in the [makefile](makefile) to `yes`; if not, set it to `no`. If you are using *IntelMPI* set `USEINTELMPI` in the makefile to `yes`; if not, set it to `no`.

In the *top level directory* of the code, do

```bash
make -j
```


## Runing COLDREG

To run COLDREG with a 32x32x32 image registration test example do

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