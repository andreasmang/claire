# Running COLDREG




## Content

* Running a Test Problem
* Registering Images
* Controlling the Output
* Parallel Execution




## Running a Test Problem

To run a simple **synthetic registration problem** using some default set of parameters do:

```bash
./bin/claire
```

To run the problem with different grid size use the `-nx` option; i.e., for a problem of size 128x64x256 do `-nx 128x64x256`.




## Registering Images

To register two NIfTI images `mR.nii.gz` (reference image) and `mT.nii.gz` do 

```bash
./bin/claire -mr /path/to/image/mR.nii.gz -mt ./path/to/image/mT.nii.gz
```

**Important:** The images have to be affinely preregistered (i.e., have the same grid size).




## Controlling the Output

To output registration results add the `-x` option. This also requires you to add the output path `-x /path/to/results/` (**mandatory**):

```bash
./bin/claire -mr /path/to/image/mR.nii.gz -mt ./path/to/image/mT.nii.gz -x /path/to/results/ -velocity
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
residual-t=0.nii.gz             | residual / mismatch before registration
residual-t=1.nii.gz             | residual / mismatch after registration
velocity-field-2norm.nii.gz     | l2 norm of velocity field




## Post Processing 

### Post Processing 

Map (transport) template image to reference space (registration performed from template to reference space):

```bash
./bin/clairetools -mt results/template-image.nii.gz -deformimage -v1 results/velocity-field-x1.nii.gz -v2 results/velocity-field-x2.nii.gz -v3 results/velocity-field-x3.nii.gz
```


Map (transport) reference image to template space (registration performed from template to reference space; add `-r2t` option):

```bash
./bin/clairetools -mt results/reference-image.nii.gz -deformimage -v1 results/velocity-field-x1.nii.gz -v2 results/velocity-field-x2.nii.gz -v3 results/velocity-field-x3.nii.gz -r2t
```


Map (transport) label map to reference space (registration performed from template to reference space):

```bash
./bin/clairetools -mt results/template-image.nii.gz -deformimage -v1 results/velocity-field-x1.nii.gz -v2 results/velocity-field-x2.nii.gz -v3 results/velocity-field-x3.nii.gz
```


## Parallel Execution

To run the code with 2 MPI tasks do (assuming you use `mpiexec`):

```bash
mpiexec -n 2 ./bin/claire
```

ACCFFT ([http://accfft.org](http://accfft.org)) will automatically decide on the data distribution. ACCFFT uses a pencil decomposition, i.e., assuming we have p = p1 p2 MPI tasks and n = n1 n2 n3 grid points. Then, the data is going to be distributed so that each MPI task gets (n1 / p1)(n2 / p2) n3 data points. To control the layout you can use the `-np` option. For instance, for 20 MPI tasks, you could use `-np 4x5`, which yields (n1 / 4)(n2 / 5) n3 data points for each MPI task.

```bash
mpiexec -n 20 ./bin/claire -np 4x5
```
