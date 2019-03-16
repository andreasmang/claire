# CLAIRE: The Binaries

CLAIRE has two binaries: `claire` and `clairetools`.

* `claire`: perform registrations
* `clairetools`: post and pre-processing

We provide several examples for executing these binaries in the [doc/examples](examples) subfolder. We briefly explain these examples below.

These binaries can be found in the `bin` folder after CLAIRE has been built successfully. To learn more about building claire take a look at our **quick installation guide** [doc/README-INSTALL-QUICK.md](README-INSTALL-QUICK.md) or our **detailed installation guide** [doc/README-INSTALL.md](README-INSTALL.md).


## Content
* [Simple Examples: `claire`](#claireexamples)
* [Simple Examples: `clairetools`](#toolsexamples)


## Simple Examples: `claire` <a name="claireexamples"></a>

To learn more about the options of

### Example 01: Synthetic Problem

In [runclaire01.sh](examples/runclaire01.sh) we execute CLAIRE for a synthetic test problem of size 32x32. We use default settings for our solver.
```bash
$bindir/claire -synthetic 0
```

The flag `-synthetic` allows one to select several smooth test problems. To change the problem size simply add the `-nx <n1xn2xn3>` flag, where `<n$i$>` represents the problem size in each spatial direction (i.e., -nx <128x128x128> executes CLAIRE with a problem size of `$128^3$`.)




## Simple Examples: `clairetools` <a name="claireexamples"></a>



## Content

* Running a Test Problem
* Registering Images
* Controlling the Output
* Parallel Execution
* Estimating the Regularization Parameter
* Post Processing




## Running a Test Problem

To run a simple **synthetic registration problem** using some default set of parameters do:

```bash
./bin/claire -synthetic 0
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


## Estimating the Regularization Parameter

CLAIRE features a method based on parameter continuation to identify an adequate regularization parameter for the velocity field. The search uses the determinant of the deformation gradient as a "metric". The user defines a lower bound for the Jacobian via the `-jbound <dbl>` option (the upper bound is `1/<dbl>`). To perform the search specify the `-train <type>` option. There are two strategies implemented: A simple reduction of the regularization parameter until the bound is hit (use `-train reduce`) and a binary search (use `-train binary`). An example is

```bash
./bin/claire -synthetic 0 -train binary -jbound 0.8
```

For neuroimaging applications with the standard regularization model (penalty on the divergence of v and H1-regularity for v) a bound of 0.2 is recommended. For H2 regularity, smaller bounds can be used (as volume conservation is not enforced). The test problem above is very smooth. Hence the large bound.

Once the optimal regularization parameter is found, the user can run CLAIRE with this parameter on similar sets of images using a similar parameter continuation scheme (using parameter continuation is recommended since it usually yields a significant speedup for small regularization parameters):

```bash
./bin/claire -synthetic 0 -betacont 1E-3
```

The continuation starts with a regularization parameter of 1.0 and reduces this parameter until the target regularization parameter (1E-3 in the example above) is reached.


## Post Processing

Map (transport) template image to reference space (registration performed from template to reference space):

```bash
./bin/clairetools -ifile inputfolder/reference-image.nii.gz -deformimage -v1 inputfolder/velocity-field-x1.nii.gz -v2 inputfolder/velocity-field-x2.nii.gz -v3 inputfolder/velocity-field-x3.nii.gz -xfile outputfolder/output-file.nii.gz
```


Map (transport) reference image to template space (registration performed from template to reference space; add `-r2t` option):

```bash
./bin/clairetools -ifile inputfolder/reference-image.nii.gz -deformimage -v1 inputfolder/velocity-field-x1.nii.gz -v2 inputfolder/velocity-field-x2.nii.gz -v3 inputfolder/velocity-field-x3.nii.gz -xfile outputfolder/output-file.nii.gz -r2t
```


Map (transport) label map to reference space (registration performed from template to reference space; you need to specify the label IDs of the labels that will be transported):

```bash
./bin/clairetools -ifile inputfolder/template-labels.nii.gz -tlabelmap -v1 inputfolder/velocity-field-x1.nii.gz -v2 inputfolder/velocity-field-x2.nii.gz -v3 inputfolder/velocity-field-x3.nii.gz -labels 1,2,10,40 -xfile outputfolder/output-file.nii.gz
```

Create Ravens map:

```bash
./bin/clairetools -ifile inputfolder/template-labels.nii.gz -computeravensmap -v1 inputfolder/velocity-field-x1.nii.gz -v2 inputfolder/velocity-field-x2.nii.gz -v3 inputfolder/velocity-field-x3.nii.gz -labels 1,2,10,40 -xfile outputfolder/ravens-map.nii.gz
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
