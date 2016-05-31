# COLDREG

**COLDREG** implements a parallel solver for *Constrained Large Deformation Diffeomorphic Image Registration*. Additional information on the methodology can be found in [doc/README-REFERENCES.md](doc/README-REFERENCES.md).

## Installation

If there are any issues or you have any questions send an email to <andreas@ices.utexas.edu>.

### Requirements

* **cmake** (at least version 2.8; [https://cmake.org](https://cmake.org))
* **python** (version 2.7)
* **mpicc** and **mpicxx**, and the MPI libraries, should be in the path

If there are issues with compiling the code, take a look at [doc/README-INSTALL.md](doc/README-INSTALL.md). We also provide an FAQ in [doc/README-INSTALL.md](doc/README-INSTALL.md) that lists common problems with compiling the code.


### Installation Instructions

#### Build Dependencies

```bash
cd external
./build_libs.sh --build
source libs/environment_vars.sh
cd ..
```

#### Build COLDREG

* Do you use an *intel compiler*? Set `USEINTEL` in the [makefile](makefile) to `yes` or `no`.
* Are you using *Intel MPI*? Set `USEINTELMPI` in the [makefile](makefile) to `yes` or `no`.

In the *top level directory* of the code, do

```bash
make -j
```


## Run COLDREG


### Test Problem

To run an image registration test example do:

```bash
./bin/runcoldreg
```

To run the code with different grid sizes use the `-nx` option (i.e., for a 128x128x128 problem, use `-nx 128x128x128`).

### Using Input Images

To run an image registration problem with input images do:

```bash
./bin/runcoldreg -mr ./external/mR.nii.gz -mt ./external/mT.nii.gz -betav 1E-2 -regnorm h2s -xresults -x ./results
```

Here, `-mr ./external/mR.nii.gz` defines the *reference image* (fixed image), `-mt ./external/mT.nii.gz` the *template image* (image to be registered), `-betav 1E-2` the *regularization weight*,  `-regnorm h2s` the *regularization norm* (H2-seminorm in this case), `-x ./results` the *output folder*, and `-xresults` enables the output of images, the computed velocity field, the deformation map, and derived measures.

**Important**: We assume that the images have been **affinely pre-registered** (same voxel dimensions and grid size). More details about these options and the output can be found in [doc/README-RUNME.md](doc/README-RUNME.md).

### Options

To see the basic options do:

```bash
./bin/runcoldreg -help
```

For more advanced options do:

```bash
./bin/runcoldreg -advanced
```

You can also find a list of the available options for the binary in [doc/help.txt](doc/help.txt) and [doc/advanced-help.txt](doc/advanced-help.txt).


## Advanced Instructions

More information on how to **add**, **install**, and **link** these libraries, can be found in [doc/README-INSTALL.md](doc/README-INSTALL.md). More details about how to run the software can be found in [doc/README-RUNME.md](doc/README-RUNME.md).

## License

Read the [LICENSE](LICENSE) file for more details.
