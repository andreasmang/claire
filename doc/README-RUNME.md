# Running COLDREG




## Content

* Running a Test Problem
* Registering Images 
	* Mandatory Input Arguments 
	* Output Arguments
* Parallel Execution




## Running a Test Problem

To run a simple **synthetic registration problem** using some default set of parameters do:

```bash
./bin/runcoldreg
```



## Registering Images

To register two NIfTI images `mR.nii.gz` (reference image) and `mT.nii.gz` of size 256x256x256 do 

```bash
./bin/runcoldreg -mr /path/to/image/mR.nii.gz -mt ./path/to/image/mT.nii.gz -nx 256x256x256
```

The size argument `-nx 256x256x256` is **mandatory**.



### Controlling the Output

To output registration results add the `-xresults` option. This also requires you to add the output path `-x /path/to/results/` (**mandatory**):

```bash
./bin/runcoldreg -mr /path/to/image/mR.nii.gz -mt ./path/to/image/mT.nii.gz -nx 256x256x256 -xresults -x /path/to/results/
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

