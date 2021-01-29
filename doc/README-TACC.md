# TACC Systems (University of Texas at Austin)

Information about the Texas Advanced Computing Center ([TACC](https://www.tacc.utexas.edu)) can be found at [https://www.tacc.utexas.edu](https://www.tacc.utexas.edu).In the following we will mainly provide information on which modules you need to load on the individual TACC systems to compile the code and the libraries. For more information on the systems please visit the [TACC](https://www.tacc.utexas.edu) web pages.

CLAIRE GPU code has been tested on the following systems:


## [Longhorn](https://portal.tacc.utexas.edu/user-guides/longhorn)

Make sure the following modules are loaded before building external dependencies for CLAIRE.

### Modules:

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

```bash
export LD_LIBRARY_PATH=/opt/apps/gcc7_3/mvapich2-gdr/2.3.4/include/:$LD_LIBRARY_PATH
```

### Dependencies:

* nifticlib
* petsc-lite-3.11.4

```bash
./build_libs.sh --enableCUDA=1 --gpu=V100 --bnifti
./build_libs.sh --enableCUDA=1 --gpu=V100 --bpetsccudasgl
```


# RCDC Systems (University of Houston)

Information about the Research Computing Data Core ([RCDC](https://www.uh.edu/rcdc/resources/)) center at the Univeristy of Houston can be found at [https://www.uh.edu/rcdc/resources](https://www.uh.edu/rcdc/resources). In the following we will provide information on which modules you need to load on the individual TACC systems to compile the cod and the libraries.


## [Sabine](https://www.uh.edu/rcdc/resources/hpc/sabine)

GPU version of CLAIRE (GIT version v0.07-363-gb5ed; date: 11/21/2020)

### Modules

```bash
module load python$
module load CMake$
module load OpenMPI/intel/4.0.1~$
module load CUDA/10.0.130$
```


### Build

```bash
cd deps
make -j
cd ..
source deps/env_source.sh
make -j
```

### Batch Job Submission

```bash
#!/bin/bash
#SBATCH -J claire
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o job%j.out
#SBATCH -e job%j.err
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=64GB
#SBATCH -t 01:30:00
#SBATCH --mail-user=YOUR EMAIL
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

### directory of your code
CDIR=${HOME}/code/claire

nvidia-smi

#### define paths
source ${CDIR}/deps/env_source.sh

mpirun ${CDIR}/bin/claire -synthetic 2 -nx 128
```
