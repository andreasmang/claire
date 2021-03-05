# CLAIRE

* Are you looking for **examples**? Check the [doc/examples](https://github.com/andreasmang/claire/tree/master/doc/examples) folder.
* Are in interested in **how CLAIRE works**? Check the [documentation](#clairedoc).
* Are you interested in **what CLAIRE is**? Read the [about](#claireabout) section.

## About <a name="claireabout"></a>

**CLAIRE** stands for *Constrained Large Deformation Diffeomorphic Image Registration*. It is a C/C++ software package for velocity-based diffeomorphic image registration in three dimensions. Its performance is optimized for multi-core systems. It uses MPI for data parallelism, and has been demonstrated to scale on several supercomputing platforms. CLAIRE can be executed on large-scale state-of-the-art computing systems as well as on local compute systems with limited resources.

<p align="center">
<img src="doc/figs/claire4brains.jpg" alt="CLAIRE4Brains"  width="800"/>
</p>

If there are any issues, you have questions, you would like to give us feedback or you have feature requests, do not hesitate to send an email to <andreas@math.uh.edu>.

If you plan on using CLAIRE in your research please cite the following manuscript:
A. Mang, A. Gholami, C. Davatzikos & G. Biros. *CLAIRE: A distributed-memory solver for constrained large deformation diffeomorphic image registration*. SIAM Journal on Scientific Computing 41(5):C548--C584, 2019 [[arxiv](https://arxiv.org/abs/1808.04487), [sisc](https://epubs.siam.org/doi/abs/10.1137/18M1207818)].

## Documentation <a name="clairedoc"></a>
* [Quick Installation Guide](doc/README-INSTALL-QUICK.md)
* [Detailed Installation Guide](doc/README-INSTALL.md)
* [Examples](doc/README-RUNME.md)
* [News](doc/README-NEWS.md)
* [Publications](doc/README-REFERENCES.md)

The links above point to individual markdown files. These files can be found in the [doc](https://github.com/andreasmang/claire/tree/master/doc) subfolder. Basic examples for how to execute CLAIRE can be found in the [doc/examples](https://github.com/andreasmang/claire/tree/master/doc/examples) folder. The NIREP dataset used to test CLAIRE can be downloaded here: [NIREP Data](https://github.com/andreasmang/nirep).

# Compatibility and Dependencies
The compiler needs C++11 support.

|Test   | Compiler  | MPI            | CUDA | PETSc  | CPU    | GPU   | System       |
|---    |---------- |-----           |------|------- |---     |---    |---           |
|b5213fa| GCC 9.3   | OpenMPI 4.0.3  | 11.0 | 3.14.2 | x86_64 | GA102 | Ubuntu 20.04 |
|6f40316| GCC 9.3   | OpenMPI 4.0.3  | 11.1 | 3.14.2 | x86_64 | GK110 | Ubuntu 20.04 |
|4967052| GCC 8.4   | OpenMPI 1.10.2 | 10.1 | 3.12.4 | x86_64 | GK110 | Ubuntu 16.04 |
|4967052| GCC 5.4.0 | OpenMPI 1.10.2 | 10.0 | 3.12.4 | x86_64 | GM200 | Ubuntu 16.04 |
|4967052| GCC 7.4   | OpenMPI 4.0.1  | 10.1 | 3.12.4 | x86_64 | GP100 | Ubuntu 16.04 |
|4967052| GCC 4.8.5 | OpenMPI 3.1.6  | 10.2 | 3.12.4 | Power9 | GV100 | CentOS 7.8   |
|4967052| XLC 16.1  | Spectrum 10.3  | 10.2 | 3.12.4 | Power9 | GV100 | RHEL 7.8     |

## Dependencies
* MPI (with GPU support (CUDA-aware MPI) for multi-GPU multi-node)
* PETSc (with GPU support)
* Nifti (./deps or libnifti-dev)
* ZLib

## Known issues
* if MPI is not compiled with CUDA-aware options, add the file `.petscrc` to the working directory and add the option `-use_gpu_aware_mpi 0`
* CUDA >= 11.0 is only supported with PETSc >= 3.14.
* Kepler GPUs work with PETSc 3.12.4  (others not tested)
* Compiling PETSc with CUDA support on cluster login nodes without GPUs might fail
* PNETCDF is currently not tested for GPUs

## License
Read the [LICENSE](https://github.com/andreasmang/claire/tree/master/LICENSE) file for more details.


## Contributors
George Biros, Malte Brunn, Amir Gholami, James Herring, Naveen Himthani, Andreas Mang, Miriam Mehl
