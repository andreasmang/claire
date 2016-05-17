COLDREG

COLDREG implements a parallel solver for Constrained Large Deformation Diffeomorphic Image Registration. Additional information on the methodology can be found in doc/README-REFERENCES.md.

1) Installation

If there are any issues or you have any questions send an email to andreas@ices.utexas.edu.

1.1) Requirements

[] cmake (at least version 2.8; https://cmake.org)
[] python (version 2.7)
[] mpicc and mpicxx should be in the path


1.2) Installation Instructions

1.2.1) Build Dependencies

cd external
./build_libs.sh --build
source libs/environment_vars.sh
cd ..

1.2.2) Build COLDREG

[] Do you use an intel compiler? Set USEINTEL in the makefile to `yes` or `no`.
[] Are you using Intel MPI? Set `USEINTELMPI` in the makefile to `yes` or `no`.

To build the binary, in the top level directory of the code, do

make -j


2) Run COLDREG

[] Run 32x32x32 image registration test example:

./bin/runcoldreg

[] Run image registration problem with input images:

./bin/runcoldreg -mr ./external/mR.nii.gz -mt ./external/mT.nii.gz -nx 256x256x256 -betav 1E-2 -regnorm h2s -xresults -x ./results

Here, `-mr ./external/mR.nii.gz` defines the reference image (fixed image), `-mt ./external/mT.nii.gz` the template image (image to be registered), `-nx 256x256x256` the size of the images, `-betav 1E-2` the regularization weight,  `-regnorm h2s` the regularization norm (H2-seminorm in this case), `-x ./results` the output folder, and `-xresults` enables the output of images, the computed velocity field, the deformation map, and derived measures.

[] General options:

./bin/runcoldreg -help

[] Advanced options

./bin/runcoldreg -advanced


You can also find a list of the available options for the binary in doc/help.txt and advanced-help.txt.


3) Advanced Instructions

More information on how to add, install, and link these libraries, can be found in doc/README-INSTALL.md.


4) License

Read the LICENSE file for more details.