# CLAIRE

## About

**CLAIRE** stands for *Constrained Large Deformation Diffeomorphic Image Registration*. It is a C/C++ software package for velocity-based diffeomorphic image registration in three dimensions. Its performance is optimized for multi-core systems. It uses MPI for data parallelism, and has been demonstrated to scale on several supercomputing platforms. CLAIRE can be executed on large-scale state-of-the-art computing systems as well as on local compute systems with limited resources.

<p align="center">
<img src="https://github.com/andreasmang/claire/blob/master/doc/figs/claire4brains.jpg" alt="CLAIRE4Brains"  width="800"/>
</p>

If there are any issues, you have questions, you would like to give us feedback or you have feature requests, do not hesitate to send an email to <andreas@math.uh.edu>.

If you plan on using CLAIRE in your research please cite the following manuscript: A. Mang, A. Gholami, C. Davatzikos & G. Biros. *CLAIRE: A distributed-memory solver for constrained large deformation diffeomorphic image registration*. arXiv.1808.0448. 2018.

## Documentation
* [Installation Guide](doc/README-INSTALL-QUICK.md)
* [Detailed Installation Guide](doc/README-INSTALL.md)
* [Examples](doc/README-RUNME.md)
* [Executables](doc/README-EXEC.md)
* [News](doc/README-NEWS.md)
* [Publications](doc/README-REFERENCES.md)

## Requirements
* **cmake** (at least version 2.8; [https://cmake.org](https://cmake.org))
* **python** (version 2.7)
* **mpicc** and **mpicxx**, and the MPI libraries, should be in the path

If there are issues with compiling the code, take a look at [Installation Guide](doc/README-INSTALL.md) or contact us.

## License
Read the [LICENSE](LICENSE) file for more details.
