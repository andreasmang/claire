# Running COLDREG




## Content

* Running a Test Problem
* Registering Images
* Controlling the Output
* Parallel Execution




## Running a Test Problem

To run a simple **synthetic registration problem** using some default set of parameters do:

```bash
./bin/runcoldreg
```

To run the problem with different grid size use the `-nx` option; i.e., for a problem of size 128x64x256 do `-nx 128x64x256`.




## Registering Images

To register two NIfTI images `mR.nii.gz` (reference image) and `mT.nii.gz` do 

```bash
./bin/runcoldreg -mr /path/to/image/mR.nii.gz -mt ./path/to/image/mT.nii.gz
```

**Important:** The images have to be affinely preregistered (i.e., have the same grid size).




## Controlling the Output

To output registration results add the `-xresults` option. This also requires you to add the output path `-x /path/to/results/` (**mandatory**):

```bash
./bin/runcoldreg -mr /path/to/image/mR.nii.gz -mt ./path/to/image/mT.nii.gz -xresults -x /path/to/results/
```

You can also add a **prefix** by doing `-x /path/to/results/prefix_`. This will add the specified `prefix_` infront of every file written to disc.

The outputs are:

file                            | explanation
--------------------------------|--------------------------------------------
reference-image.nii.gz          | reference image (fixed image) *mR*
template-image.nii.gz           | template image (image to be registered) *mT*
velocity-field-x1.nii.gz        | x1 component of computed velocity field *v*
velocity-field-x2.nii.gz        | x2 component of computed velocity field *v*
velocity-field-x3.nii.gz        | x3 component of computed velocity field *v*
deformed-template-image.nii.gz  | deformed template image (template image after registration)
deformation-map-x1.nii.gz       | x1 component of computed deformation map *y*
deformation-map-x2.nii.gz       | x2 component of computed deformation map *y*
deformation-map-x3.nii.gz       | x3 component of computed deformation map *y*
det-deformation-grad.nii.gz     | determinant of deformation gradient (jacobian; det(grad(y)))
residual-after.nii.gz           | residual / mismatch after registration
residual-before.nii.gz          | residual / mismatch before registration
velocity-field-2norm.nii.gz     | l2 norm of velocity field




## Parallel Execution

To run the code with 2 MPI tasks do (assuming you use `mpiexec`):

```bash
mpiexec -n 2 ./bin/runcoldreg
```

ACCFFT ([http://accfft.org](http://accfft.org)) will automatically decide on the data distribution. ACCFFT uses a pencil decomposition, i.e., assuming we have p = p1 p2 MPI tasks and n = n1 n2 n3 grid points. Then, the data is going to be distributed so that each MPI task gets (n1 / p1)(n2 / p2) n3 data points. To control the layout you can use the `-np` option. For instance, for 20 MPI tasks, you could use `-np 4x5`, which yields (n1 / 4)(n2 / 5) n3 data points for each MPI task.

```bash
mpiexec -n 20 ./bin/runcoldreg -np 4x5
```
