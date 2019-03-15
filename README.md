# CLAIRE

**CLAIRE** stands for *Constrained Large Deformation Diffeomorphic Image Registration*. It is a C/C++ software package for velocity-based diffeomorphic image registration in three dimensions. Its performance is optimized for multi-core systems. It uses MPI for data parallelism, and has been demonstrated to scale on several supercomputing platforms. CLAIRE can be executed on large-scale state-of-the-art computing systems as well as on local compute systems with limited resources.

<p align="center">
<img src="https://github.com/andreasmang/claire/blob/master/doc/figs/claire4brains.jpg" alt="CLAIRE4Brains"  width="800"/>
</p>

If you are interested in the foundations of CLAIRE we encourage you to read our past work [doc/README-REFERENCES.md](doc/README-REFERENCES.md).

If there are any issues, you have questions, you would like to give us feedback or you have feature requests, do not hesitate to send an email to <andreas@math.uh.edu>.

If you plan on using CLAIRE in your research please cite the following manuscript:

A. Mang, A. Gholami, C. Davatzikos & G. Biros. *CLAIRE: A distributed-memory solver for constrained large deformation diffeomorphic image registration*. arXiv.1808.0448. 2018.

## Installation



### Requirements

* **cmake** (at least version 2.8; [https://cmake.org](https://cmake.org))
* **python** (version 2.7)
* **mpicc** and **mpicxx**, and the MPI libraries, should be in the path

If there are issues with compiling the code, take a look at [doc/README-INSTALL.md](doc/README-INSTALL.md). We also provide an FAQ in [doc/README-INSTALL.md](doc/README-INSTALL.md) that lists common problems with compiling the code.


### Installation Instructions

#### Installing Dependencies

The following explains how to download the libraries, build them, and set the appropriate environmental variables.

```bash
cd external
./get_libs.sh
./build_libs.sh --build
source libs/environment_vars.sh
cd ..
```

#### Installing CLAIRE


In the *top level directory* of the code, do

```bash
make -j
```

To user can change some options in the makefile:

* Are your input files nifti files (.nii)? Set `USENIFTI` in the [makefile](makefile) to `yes` or `no`.
* Do you use an *intel compiler*? Set `USEINTEL` in the [makefile](makefile) to `yes` or `no`.
* Are you using *Intel MPI*? Set `USEINTELMPI` in the [makefile](makefile) to `yes` or `no`.
* Are you going to run the code in single precision? Set `USESINGLE` in the [makefile](makefile) to `yes` or `no`.
* Are you going to use the toolbox (e.g., compute jacobians)? Set `BUILDTOOLS` in the [makefile](makefile) to `yes` or `no`.
* Are your input files netcdf files (.nc)? Set `USEPNETCDF` in the [makefile](makefile) to `yes` or `no`.

## Executing CLAIRE


### Simple Synthetic Test Problem

To run an image registration test example do:

```bash
mpirun -n 20 ./bin/claire -synthetic 0
```
To run the code with different grid sizes use the `-nx` option (i.e., for a 128x128x128 problem, use `-nx 128x128x128`).

### Run CLAIRE with Images

To run an image registration problem with input images do:

```bash
mpirun -n 20 ./bin/claire -mr ./external/mR.nii.gz -mt ./external/mT.nii.gz -beta 1E-2 -regnorm h2s -velocity -x ./results -disablesmoothing
```

Here, `-mr ./external/mR.nii.gz` defines the *reference image* (fixed image), `-mt ./external/mT.nii.gz` the *template image* (image to be registered), `-beta 1E-2` the *regularization weight*,  `-regnorm h2s` the *regularization norm* (H2-seminorm in this case), `-x ./results` the *output folder*, and `-velocity` enables the output of the computed velocity field. These images are smooth; we can disable the default smoothing by adding the `-disablesmoothing` flag to the command line. **Warning**: do not do this for real images.

**Important**: We assume that the images have been **affinely pre-registered** (same voxel dimensions and grid size). More details about these options and the output can be found in [doc/README-RUNME.md](doc/README-RUNME.md).

### Options

To see the basic options do:

```bash
./bin/claire -help
```

## Advanced Instructions

More information on how to **add**, **install**, and **link** these libraries, can be found in [doc/README-INSTALL.md](doc/README-INSTALL.md). More details about how to run the software can be found in [doc/README-RUNME.md](doc/README-RUNME.md).

## License

Read the [LICENSE](LICENSE) file for more details.
