---
title: 'CLAIRE: Constrained Large Deformation Diffeomorphic Image Registration on Parallel Computing Architectures'
tags:
  - C++
  - GPUs
  - parallel computing
  - high performance computing
  - large deformation diffeomorphic registration
  - optimal control
  - medical imaging
authors:
  - name: Malte Brunn
    affiliation: "1"
  - name: Naveen Himthani
    affiliation: "2"
  - name: George Biros
    affiliation: "2"
  - name: Miriam Mehl
    affiliation: "1"
  - name: Andreas Mang^[Corresponding Author]
    orcid: 0000-0003-4605-3290
    affiliation: "3"
affiliations:
 - name: Institute for Parallel and Distributed Systems, University Stuttgart
   index: 1
 - name: Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin
   index: 2
 - name: Department of Mathematics, University of Houston
   index: 3
date: 28 May 2021
bibliography: paper.bib
---

# Summary
[`CLAIRE`](https://andreasmang.github.io/claire) [@claire-web] is a computational framework for **C**onstrained **LA**rge deformation diffeomorphic **I**mage **RE**gistration [@Mang:2019a]. It supports highly-optimized, parallel computational kernels for (multi-node) CPU [@Mang:2016a; @Gholami:2017a; @Mang:2019a] and (multi-node multi-)GPU architectures [@Brunn:2020a; @Brunn:2021a]. `CLAIRE` uses MPI for distributed-memory parallelism and can be scaled up to thousands of cores [@Mang:2019a; @Mang:2016a] and GPU devices [@Brunn:2020a]. The multi-GPU implementation uses device direct communication. The computational kernels are interpolation for semi-Lagrangian time integration, and a mixture of high-order finite difference operators and Fast-Fourier-Transforms (FFTs) for differentiation. `CLAIRE` uses a Newton--Krylov solver for numerical optimization [@Mang:2015a; @Mang:2017a]. It features various schemes for regularization of the control problem [@Mang:2016a] and different similarity measures. `CLAIRE` implements different preconditioners for the reduced space Hessian [@Brunn:2020a; @Mang:2019a] to optimize computational throughput and enable fast convergence. It uses `PETSc` [@petsc-web] for scalable and efficient linear algebra operations and solvers and `TAO` [@petsc-web; @Munson:2015a] for numerical optimization. `CLAIRE` can be downloaded at <https://github.com/andreasmang/claire>.

# Statement of Need
Image registration is required whenever images are taken at different points in time, from different viewpoints, and/or using different imaging modalities and these images need to be compared, combined, or integrated [@Fischer:2008a; @Modersitzki:2004a; @Modersitzki:2009a; @Sotiras:2013a]. Image registration is an inverse problem. The inputs to this inverse problem are two (or more) images $m_0(x)$ (the template image) and $m_1(x)$ (the reference image) of the same object. The task of image registration is to find a plausible map $y(x)$ that establishes spatial correspondences between the reference and template image, i.e., $m_0(x) \approx m_1(y(x))$. In `CLAIRE` the set of admissible spatial transformations $y$ is limited to diffeomorphisms, i.e., maps $y$ that are continuous, one-to-one, and have a smooth inverse. `CLAIRE` is related to a prominent class of formulations for these types of problems referred to as <em>large-deformation diffeomorphic metric mapping</em> [@Beg:2005a; @Trouve:1998a; @Younes:2010a].

Diffeomorphic image registration is an indispensable tool in medical image analysis [@Sotiras:2013a]. Computing diffeomorphisms that map one image to another is expensive. Deformable image registration is an infinite-dimensional problem that upon discretization leads to nonlinear optimality systems with millions or even billions of unknowns. For example, registering two typical medical imaging datasets of size $256^3$ necessitates solving for about 50 million unknowns (in our formulation). Additional complications are the ill-posedness and non-linearty of this inverse problem [@Fischer:2008a]. Consequently, image registration can take several minutes on multi-core high-end CPUs. Many of the available methods reduce the number of unknowns by using coarser resolutions either through parameterization or by solving the problem on coarser grids; they use simplified algorithms and deliver subpar registration quality. In the age of big data, clinical population studies that require thousands of registrations are incresingly common, and execution times of individual registrations become more critical. We provide technology that allows solving registration problems for clinical datasets in seconds. In addition, we have made available to the public a software that works on multi-node, multi-GPU architectures [@Brunn:2020a; @Brunn:2021a] that allows the registration of large-scale microscopic imaging data such as CLARITY imaging [@Tomer:2014a; @Kutten:2017a].

# Highlights
`CLAIRE` can be used to register images of $2048^3$ (25 B unknowns) on 64 nodes with 256 GPUs on TACC’s Longhorn system [@Brunn:2020a]. `CLAIRE` has been used for the registration of high resolution CLARITY imaging data [@Brunn:2020a]. The GPU version of `CLAIRE` can solve clinically relevant problems (50 M unknowns) in approximately 5 seconds on a single NVIDIA Tesla V100 [@Brunn:2020a]. `CLAIRE` has also been applied to hundreds of images in brain tumor imaging studies [@Bakas:2018a; @Mang:2017c; @Scheufele:2021a], and has been integrated with models for biophysics inversion [@Mang:2018a; @Mang:2020a; @Scheufele:2020a; @Scheufele:2019a; @Scheufele:2021a; @Subramanian:2020b] and Alzheimer's disease progression [@Scheufele:2020c]. `CLAIRE` uses highly optimized computational kernels and effective, state-of-the-art algorithms for time integration and numerical optimization. Our most recent version of `CLAIRE` features a Python interface to assist users in their applications.

We provide a detailed documentation on how to execute, compile, and install `CLAIRE` on various systems at our deployment page <https://andreasmang.github.io/claire>.

# Mathematics
`CLAIRE` uses an optimal control formulation. The diffeomorphism $y(x)$ is parameterized using a smooth, stationary velocity field $v(x)$. Given the template image $m_0(x)$ and the reference image $m_1(x)$, this velocity is found by solving the partial-differential equation constrained optimization problem of the form
$$
\operatorname{minimize}_{v,m} \operatorname{dist}(m(x,t=1),m_1) + \alpha\operatorname{reg}(v)
$$

subject to
$$
\begin{aligned}
\partial_t  m(x,t) + v(x) \cdot \nabla m(x,t) &= 0 \\
m(x,t=0) & = m_0(x)
\end{aligned}
$$

The first term in the objective functional measures the proximity of the deformed template image $m(x,t=1)$ and the reference image $m_1(x)$. The default option availble in `CLAIRE` is an $L^2$-distance. The second term controls the regularity of $v$. `CLAIRE` features different Sobolev norms. The default option is an $H^1$-seminorm. The constraint models the deformation the template image (i.e., the transport of the intensities of $m_0(x)$). `CLAIRE` also features additional hard constraints for controlling the divergence of $v(x)$ [@Mang:2016a]. For optimization, we use the method of Lagrange multipliers and solve the associated Karush--Kuhn--Tucker optimality system using a Newton--Krylov reduced space method [@Mang:2015a; @Mang:2015a].


# Acknowledgements
This work was partly supported by the National Science Foundation (DMS-1854853, DMS-2009923, DMS-2012825, CCF-1817048, CCF-1725743), the NVIDIA Corporation (NVIDIA GPU Grant Program), the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany’s Excellence Strategy-EXC 2075-390740016, by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research, Applied Mathematics program under Award Number DE-SC0019393; by the U.S. Air Force Office of Scientific Research award FA9550-17-1-0190; by the Portugal Foundation for Science and Technology and the UT Austin-Portugal program, and by NIH award 5R01NS042645-11A1. Any opinions, findings, and conclusions or recommendations expressed herein are those of the authors and do not necessarily reflect the views of the DFG, AFOSR, DOE, NIH, and NSF. Computing time on the Texas Advanced Computing Centers’ (TACC) systems was provided by an allocation from TACC and the NSF. This work was completed in part with resources provided by the Research Computing Data Core at the University of Houston.

# References
