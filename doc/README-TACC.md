# TACC Systems (University of Texas at Austin)

Information about the Texas Advanced Computing Center ([TACC](https://www.tacc.utexas.edu)) can be found at [https://www.tacc.utexas.edu](https://www.tacc.utexas.edu).In the following we will mainly provide information on which modules you need to load on the individual TACC systems to compile the code and the libraries. For more information on the systems please visit the [TACC](https://www.tacc.utexas.edu) web pages.

CLAIRE GPU code has been tested on the following systems:


## [Longhorn](https://portal.tacc.utexas.edu/user-guides/longhorn)

Make sure the following modules are loaded before building external dependencies for CLAIRE.

### Using *MVAPICH2* compiler

```bash
git/2.24.1
xalt/2.8.1
autotools/1.2
TACC
cmake/3.16.1
cuda/10.2
gcc/7.3.0
mvapich2-gdr/2.3.4
```



# RCDC Systems (University of Houston)

Information about the Research Computing Data Core ([RCDC](https://www.uh.edu/rcdc/resources/)) center at the Univeristy of Houston can be found at [https://www.uh.edu/rcdc/resources](https://www.uh.edu/rcdc/resources). In the following we will provide information on which modules you need to load on the individual TACC systems to compile the cod and the libraries.


## [Sabine](https://www.uh.edu/rcdc/resources/hpc/sabine)

GPU version of CLAIRE (GIT version v0.07-363-gb5ed; date: 11/21/2020)

### Modules:

```bash
python/3.6
cmake/3.15.4
CUDA/9.2.88
PSM2/10.3.35-cuda
OpenMPI/gcc-cuda/3.1.2
```

### Dependencies:
* nifticlib
* petsc-lite-3.11.4
