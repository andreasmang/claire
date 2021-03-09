# CLAIRE: The Binaries

Go back to [README.md](../README.md).

## Content
* [Overview](#clairebins)
* [Get Help](#clairehelp)
* [Simple Examples: `claire`](#claireexmp)
	* [Synthetic Problem](#clairexmp1)
	* [Synthetic Problem (Parallel Execution)](#clairexmp2)
	* [Real Data (Parallel Execution)](#clairexmp3)
	* [Regularization Parameter Estimation](#clairexmp4)
	* [Parameter Continuation](#clairexmp5)
	* [Output Velocities](#clairexmp6)
* [Simple Examples: `clairetools`](#toolsxmp)
	* [Transporting Images](#toolsxmp1)
	* [Computing Jacobians](#toolsxmp2)
* [Testing and Benchmarks](#testing)


## Overview <a name="clairebins"></a>

CLAIRE has two main binaries: `claire` and `clairetools`.

  * `claire`: perform registrations
  * `clairetools`: post and pre-processing

We provide **several examples** for executing these binaries in the [doc/examples](https://github.com/andreasmang/claire/tree/gpu/examples) subfolder. We briefly explain these examples below.

In addition to that we have added two binaries named `test` and `benchmark` for developers to test the code. They are described in greater detail in the [Testing and Benchmarks](#testing) section below.

These binaries can be found in the `bin` folder after CLAIRE has been built successfully. To learn more about building claire take a look at our [installation guide](README-INSTALL.md) found in [doc/README-INSTALL.md](README-INSTALL.md).


## Get Help <a name="clairehelp"></a>

To learn more about the options available in `claire` and `clairetools` add the `-help` flag:

```bash
claire -help
```

```bash
clairetools -help
```

## Simple Examples: `claire` <a name="claireexmp"></a>

### Example 01: Synthetic Problem <a name="clairexmp1"></a>

In [runclaire01.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire01.sh) we execute CLAIRE for a synthetic test problem of size 32x32x32. We use default settings for our solver:
```bash
$bindir/claire -synthetic 0
```

`$bindir` is the directory in which the binary is located. The flag `-synthetic` allows one to select several smooth test problems. To change the problem size simply add the `-nx <n1xn2xn3>` flag, where `<n$i$>` represents the problem size in each spatial direction (i.e., `-nx 128x128x128` executes CLAIRE with a problem size of `128x128x128`.) We recommend executing CLAIRE in parallel for larger problem sizes (see [example 2](#clairexmp2))

### Example 02: Synthetic Problem (Parallel Execution) <a name="clairexmp2"></a>

In [runclaire02.sh](examples/runclaire02.sh) we execute CLAIRE for a synthetic test problem of size 128x128x128 in parallel. We use 20 MPI tasks. We use default settings for our solver:

```bash
mpirun -np 20 $bindir/claire -synthetic 0 -nx 128
```

The options used with `claire` are explained in [example 1](#clairexmp1). The key difference is the instruction `mpirun -np 20` infront of the executable. This instructs your compute node to use 20 MPI tasks. CLAIRE will determine the processor layout for you. We recommend executing CLAIRE in parallel. If you use `mpiexec` replace `mpirun -np 20` with `mpiexec -n 20`.


### Example 03: Real Data <a name="clairexmp3"></a>

In [runclaire03.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire03.sh) we execute CLAIRE for real medical images (in NIfTI format) of size 128x150x128. We use 20 MPI tasks. The data can be found in the [docs/data](data) subdirectory. We use default settings for our solver:

```bash
mpirun -np 20 $bindir/claire -mr $datdir/brain01.nii.gz \
                             -mt $datdir/brain02.nii.gz
```

**<span style="color:red">Important:</span>** The images have to be **affinely preregistered** (i.e., have the same grid size).

`$datdir` points to the location where the data is located. The `-mr` flag identifies the image to be used as a **reference image** (alas, *fixed* or *target* image) and the `-mt` flag identifies the image to be used as a **template image** (i.e., the image to be registered; alas *floating* or *moving* image). The line break (backslash `\`) is only added for readability.


### Example 04: Regularization Parameter Estimation <a name="clairexmp4"></a>

In [runclaire04.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire04.sh) we execute CLAIRE to automatically identify an adequate regularization parameter for a given set of images. We use default settings for our solver:

```bash
mpirun -np 20 $bindir/claire -mr $datdir/brain01.nii.gz \
                             -mt $datdir/brain02.nii.gz -train binary
```

Running `claire` on real image data is explained in [example 3](#clairexmp3). We use a method based on parameter continuation to identify an adequate regularization parameter. The search uses the determinant of the deformation gradient as a "metric". The user can define a lower bound for the Jacobian via the `-jbound <dbl>` option (the upper bound is `1/<dbl>`). To perform the search specify the `-train <type>` option. There are two strategies implemented: A simple reduction of the regularization parameter until the bound is hit (use `-train reduce`) and a binary search (use `-train binary`). The line break (backslash `\`) is only added for readability.


### Example 05: Parameter Continuation <a name="clairexmp5"></a>

In [runclaire05.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire05.sh) we show how to execute CLAIRE using a parameter continuation scheme with a target regularization parameter for the velocity. We use default settings for our solver:

```bash
mpirun -np 20 $bindir/claire -mr $datdir/brain01.nii.gz \
                             -mt $datdir/brain02.nii.gz \
                             -betacont 7.75e-04
```

We have observed that a parameter continuation scheme speeds up the rate of convergence of our solver. We recommend using it in practical settings. We show how to estimate an adequate regularization parameter in  [example 4](#clairexmp4). The line breaks (backslashes `\`) are only added for readability.


### Example 06: Output Velocities <a name="clairexmp6"></a>

In [runclaire06.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire06.sh) we show how to store the computed velocity field on file. We use default settings for our solver:

```bash
mpirun -np 20 $bindir/claire -mr $datdir/brain01.nii.gz \
                             -mt $datdir/brain02.nii.gz \
                             -betacont 7.75e-04  -x ./ -velocity
```

This example is a direct extension of [example 4](#clairexmp4). The only difference is that we added an output. We need to provide an output folder. This is done with the `-x <folder>` option. We write the output to the current directory (`./`). **The outputs in CLAIRE will have default names**. If you prefer to store all files for multiple runs in a single folder, we recommend to use a prefix:
```bash
-x /my/output/folder/name/PREFIX_
```

The `-velocity` option tells CLAIRE to write out the velocity field. There are multiple other outputs available, most of which can be computed from the velocity field. This can be done using `clairetools`. To learn more about how to use `clairetools` continue reading. The line breaks (backslashes `\`) are only added for readability.


## Simple Examples: `clairetools` <a name="toolsxmp"></a>

### Example 01: Transporting Images <a name="toolsxmp1"></a>

In [runtools01.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runtools01.sh) we show how to transport an image (i.e., e.g., compute the deformed template image after a velocity has been computed using `claire`.)

```bash
mpirun -np 20 $bindir/clairetools -v1 velocity-field-x1.nii.gz       \
                                  -v2 velocity-field-x2.nii.gz       \
                                  -v3 velocity-field-x3.nii.gz       \
                                  -ifile $datdir/brain01.nii.gz      \
                                  -xfile brain01-transported.nii.gz -deformimage
```

The input are the three components of the computed velocity (`-v$i$ velocity-field-x$i$.nii.gz `) and the image to be transported (`-ifile $datdir/brain01.nii.gz`; `$datdir` points to the folder the data is located in, i.e., [doc/data](https://github.com/andreasmang/claire/tree/gpu/doc/data)). The output is the transported brain image (`-xfile brain01-transported.nii.gz`). The user can add a path as prefix if desired. The command to tell `clairetools` that we are interested in solving the forward problem (i.e., transporting/deforming an image) is `-deformimage`. The line breaks (backslashes `\`) are only added for readability.


### Example 02: Computing Jacobians <a name="toolsxmp2"></a>

In [runtools02.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runtools02.sh) we show how to compute the determinant of the deformation gradient (alas Jacobian) from a velocity field that has been computed using `claire`.

```bash
mpirun -np 20 $bindir/clairetools -v1 velocity-field-x1.nii.gz       \
                                  -v2 velocity-field-x2.nii.gz       \
                                  -v3 velocity-field-x3.nii.gz       \
                                  -x ./ -detdefgrad
```



## Testing and Benchmarks <a name="testing"></a>

We have implemented numerical tests for the main computational kernels available in CLAIRE to study the performance and accuracy of the mathematical operators that appear in the optimality system. For reproducability, we also posted the NIREP data at [https://github.com/andreasmang/nirep](https://github.com/andreasmang/nirep). We have used this data extensively in our most recent [prior work](README-REFERENCES.md).

To build **binaries for testing** set `BUILD_TEST=yes` in the makefile. This will build two binaries: `benchmark` and `test`.  The `test` application allows users to check the main computational kernels:
* the interpolation kernels (`-interp` flag; checks interpolation accuracy)
* the regularization operators (e.g., biharmonic or laplacian regularization operators; `-reg` flag)
* numerical differentiation (`-diff`; checks differentiation accuracy)

These tests are implemented in the `*.cpp` files in the [UnitTests](https://github.com/andreasmang/claire/tree/gpu/src/UnitTests) subfolder.

We have also implemented several **high-level numerical checks** to assess the performance of our methodology and ensure that the mathematical operators are correct.

The `benchmark` binary allows users to (i) check the accuracy of the forward operator and (ii) report the runtime for evaluating several key mathematical operators. Use the `-help` flag to see all options. The main tests are the following:
* check error of forward operator (solve the forward problem for $v$ and $-v$ and check error with respect to initial condition for $t=0$)
* report runtimes for evaluating the forward operator, the gradient operator and the Hessian matvec.

The **tests/debug options** directly available within the `claire` binary are the following (Use the `-help` flag to see all options.):
* The default test in CLAIRE is to consider synthetic test problems. The user can select between several test problems of varying complexity by setting the flag `-synthetic i`, where `i` selects the particular test case (valid values for `i` are `0`, `1`, ..., `5`).
* The user can control the verbosity level of `claire` by using the `-verbose 2` flag (debug mode verbosity). This will, e.g., enable command window outputs such as the residual in each iteration of the Krylov subspace method used to compute the search direction (and much more).
* The Newton--Krylov solver monitors several critical values during the course of the iterations. The user can see outputs such as
	* number of (Gauss--)Newton iterations.
	* trend of the objective value (has to decrease monotonically)
	* trend of the gradient norm (should decrease but not necessarily monotonically)
	* number of line search steps (should be 1 for a Newton method subject to accuracy requirements)
	* and much more...
* The accuracy of the symmetry of the discretized Hessian operator can be monitored by enabling the `-checksymmetry` flag in `claire`. Notice that we consider an optimize-then-discretize approach and numerical schemes that do not preserve symmetry; consequently, in our current implementation, the Hessian is only symmetric up to the discretization error of the adjoint operators (~1e-2).
* The approximation accuracy of the gradient and Hessian can be monitored by enabling the `-derivativecheck` flag in `claire`. We report the assymptotic behavior of the Taylor expansion. The approximation error should decrease with decreasing perburbation (i.e., the error should converge). However, we do, in general not expect to observe quadratic or cubic convergence (as we would if we considered a discretize-then-optimize approach).
