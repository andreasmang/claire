# CLAIRE

**CLAIRE** implements a parallel solver for *Constrained Large Deformation Diffeomorphic Image Registration*. Additional information on the methodology can be found in [doc/README-REFERENCES.md](doc/README-REFERENCES.md).

## Installation

If there are any issues or you have any questions send an email to <andreas@math.uh.edu>.

### Requirements

* **cmake** (at least version 2.8; [https://cmake.org](https://cmake.org))
* **python** (version 2.7)
* **mpicc** and **mpicxx**, and the MPI libraries, should be in the path

If there are issues with compiling the code, take a look at [doc/README-INSTALL.md](doc/README-INSTALL.md). We also provide an FAQ in [doc/README-INSTALL.md](doc/README-INSTALL.md) that lists common problems with compiling the code.


### Installation Instructions

#### Build Dependencies

The following explains how to get the libraries, build them, and set the environmental variables.

```bash
cd external
./get_libs.sh
./build_libs.sh --build
source libs/environment_vars.sh
cd ..
```

#### Build CLAIRE


In the *top level directory* of the code, do

```bash
make -j
```

To user can change some options in the makefile:

* Are you going to run the code in single precision? Set `USESINGLE` in the [makefile](makefile) to `yes` or `no`.
* Are you going to use the toolbox (e.g., compute jacobians)? Set `BUILDTOOLS` in the [makefile](makefile) to `yes` or `no`.
* Are your input files netcdf files (.nc)? Set `USEPNETCDF` in the [makefile](makefile) to `yes` or `no`.
* Are your input files nifti files (.nii)? Set `USENIFTI` in the [makefile](makefile) to `yes` or `no`.
* Do you use an *intel compiler*? Set `USEINTEL` in the [makefile](makefile) to `yes` or `no`.
* Are you using *Intel MPI*? Set `USEINTELMPI` in the [makefile](makefile) to `yes` or `no`.

## Run CLAIRE


### Test Problem

To run an image registration test example do:

```bash
./bin/claire -synthetic 0
```

To run the code with different grid sizes use the `-nx` option (i.e., for a 128x128x128 problem, use `-nx 128x128x128`).

### Using Input Images

To run an image registration problem with input images do:

```bash
./bin/claire -mr ./external/mR.nii.gz -mt ./external/mT.nii.gz -beta 1E-2 -regnorm h2s -velocity -x ./results -disablesmoothing
```

Here, `-mr ./external/mR.nii.gz` defines the *reference image* (fixed image), `-mt ./external/mT.nii.gz` the *template image* (image to be registered), `-beta 1E-2` the *regularization weight*,  `-regnorm h2s` the *regularization norm* (H2-seminorm in this case), `-x ./results` the *output folder*, and `-velocity` enables the output of the computed velocity field. These images are smooth; we can disable the default smoothing by adding the `-disablesmoothing` flag to the command line. **Warning**: do not do this for real images.

**Important**: We assume that the images have been **affinely pre-registered** (same voxel dimensions and grid size). More details about these options and the output can be found in [doc/README-RUNME.md](doc/README-RUNME.md).

### Options

To see the basic options do:

```bash
./bin/claire -help
```

For more advanced options do:

```bash
./bin/claire -advanced
```

You can also find a list of the available options for the binary in [doc/help.txt](doc/help.txt) and [doc/advanced-help.txt](doc/advanced-help.txt).


## Advanced Instructions

More information on how to **add**, **install**, and **link** these libraries, can be found in [doc/README-INSTALL.md](doc/README-INSTALL.md). More details about how to run the software can be found in [doc/README-RUNME.md](doc/README-RUNME.md).

## License

Read the [LICENSE](LICENSE) file for more details.
