#!/usr/bin/python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  mode=sys.argv[1]
  gpu=sys.argv[2]
  configure_options = [
    '--download-f2cblaslapack=1',
    '--with-cxx=mpicxx',
    '--with-cc=mpicc',
    '--with-mpiexec=mpiexec.hydra',
    '--with-fc=0',
    '--with-ssl=0',
    '--with-shared=0',
    '--with-64-bit-indices',
    '--with-x=0',
    'COPTFLAGS="-O3"',
    'CXXOPTFLAGS="-O3"',
  ]

    #'--with-mpi-dir='+os.environ['TACC_SPECTRUM_MPI_DIR'],
    #'--with-mpiexec='+os.path.join(os.environ['TACC_SPECTRUM_MPI_BIN'], 'mpiexec'),
    #'--with-mpiexec=mpiexec',
    #'--with-mpi-dir=/opt/apps/gcc7_3/mvapich2-gdr/2.3.4',
  #if "cuda" in mode:
  #  mode = mode + "-" + gpu
  configure_options.append('--with-petsc-arch='+mode)

  if "dbg" in mode:
    configure_options.append('--with-debugging=1')

  if "cuda" in mode:
    #configure_options.append('--with-cuda-dir='+os.environ['TACC_CUDA_DIR'])
    configure_options.append('--with-cuda=1')
    configure_options.append('--with-cudac=nvcc')
    configure_options.append('--download-cusp=yes')
    if gpu == "V100":
      # on TACC-longhorn
      configure_options.append('CUDAFLAGS="-arch=sm_70"')
    elif gpu == "P100":
      configure_options.append('CUDAFLAGS="-arch=sm_60"')
    elif gpu == "RTX":
      # on TACC-frontera
      configure_options.append('CUDAFLAGS="-arch=sm_75"')
    else:
      configure_options.append('CUDAFLAGS="-arch=sm_70"')


    configure_options.append('CUDAOPTFLAGS="-O3"')

  if "sgl" in mode:
    configure_options.append('--with-precision=single')

  if "dbl" in mode:
    configure_options.append('--with-precision=double')

  print configure_options
  configure.petsc_configure(configure_options)

