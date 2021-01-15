BUILD_PETSC   = 3.12.4
BUILD_NIFTI   = yes
BUILD_PNETCDF = no
BUILD_ACCFFT  = no

BUILD_GPU     = yes
BUILD_DEBUG   = no
BUILD_DOUBLE  = no

WITH_CUDA_MPI = yes
WITH_PETSC_OPTS =

WITH_BATCH = no

BUILD_DIR = $(PWD)/lib

GPU_VERSION =

CC = mpicc
CXX = mpicxx
NVCC = nvcc

ifneq ($(BUILD_PETSC), no)
	TARGETS += petsc
endif

ifneq ($(BUILD_NIFTI), no)
	TARGETS += nifti
endif

PETSC_OPTIONS += --with-mkl_sparse=0
PETSC_OPTIONS += --download-f2cblaslapack
PETSC_OPTIONS += --with-fortran-bindings=0
PETSC_OPTIONS += --with-valgrind=0
PETSC_OPTIONS += --with-c2html=0
PETSC_OPTIONS += --with-64-bit-indices
PETSC_OPTIONS += --with-x=0
PETSC_OPTIONS += --with-fc=0
PETSC_OPTIONS += --with-ssl=0
PETSC_OPTIONS += --COPTFLAGS='-O3'
PETSC_OPTIONS += --CXXOPTFLAGS='-O3'
PETSC_OPTIONS += --with-precision=single
PETSC_OPTIONS += --with-cc=$(CC)
PETSC_OPTIONS += --with-cxx=$(CXX)
PETSC_OPTIONS += --with-shared=1

ifeq ($(WITH_CUDA_MPI), yes)
	PETSC_OPTIONS += -use-gpu-aware-mpi=1
else
	PETSC_OPTIONS += -use-gpu-aware-mpi=0
endif

ifeq ($(WITH_BATCH), yes)
	PETSC_OPTIONS += --with-batch=1
endif

NIFTI_OPTIONS += -DCMAKE_CXX_COMPILER=$(CXX)
NIFTI_OPTIONS += -DCMAKE_C_COMPILER=$(CC)
NIFTI_OPTIONS += -Wno-dev
NIFTI_OPTIONS += -DBUILD_SHARED_LIBS:BOOL=ON

ifeq ($(BUILD_GPU), yes)
	PETSC_ARCH = gpu
	PETSC_OPTIONS += --with-cuda=1 
#	PETSC_OPTIONS += --download-cusp=yes
	PETSC_OPTIONS += --CUDAOPTFLAGS='-O3'
	PETSC_OPTIONS += --with-cudac='$(NVCC) -ccbin=$(CXX)'
ifdef GPU_VERSION
	PETSC_OPTIONS += --with-cuda-gencodearch=$(GPU_VERSION)
endif	
else
	PETSC_ARCH =cpu
	PETSC_OPTIONS += --with-cuda=0
endif

ifeq ($(BUILD_DOUBLE), yes)
	PETSC_ARCH :=$(PETSC_ARCH)_double
else
	PETSC_ARCH :=$(PETSC_ARCH)_single
endif

ifeq ($(BUILD_DEBUG), yes)
	PETSC_ARCH :=$(PETSC_ARCH)_debug
	PETSC_OPTIONS += --with-debugging=1
else
	PETSC_OPTIONS += --with-debugging=0
endif

PETSC_OPTIONS += $(WITH_PETSC_OPTS)

BASE_DIR=$(PWD)

all: config $(TARGETS)
	@echo "================================================================================"
	@echo "done"
	@echo "================================================================================"

config:
	@echo "export LD_LIBRARY_PATH=$(BUILD_DIR)/lib:\$$LD_LIBRARY_PATH" > env_source.sh
	@echo "================================================================================"
	@echo "options"
	@echo "================================================================================"
	@echo "build PETSc:      $(BUILD_PETSC)"
	@echo "build NIFTI:      $(BUILD_NIFTI)"
	@echo "================================================================================"
	@echo "build with DEBUG:          $(BUILD_DEBUG)"
	@echo "build with DOUBLE:         $(BUILD_DOUBLE)"
	@echo "build with batch system:   $(WITH_BATCH)"
	@echo "build with CUDA-aware MPI: $(WITH_CUDA_MPI)"
	@echo "build for GPU:             $(BUILD_GPU)"
	@echo "build with DEBUG:          $(BUILD_DEBUG)"
	@echo "build directory:           $(BUILD_DIR)"
	@echo "================================================================================"
	@echo "using CC:   $(CC)"
	@echo "using CXX:  $(CXX)"
	@echo "using NVCC: $(NVCC)"
	@echo "================================================================================"

petsc: petsc-lite-$(BUILD_PETSC).tar.gz
	@echo "================================================================================"
	@echo "building PETSc: $(BUILD_PETSC)"
	@echo "configure with: $(PETSC_OPTIONS)"
	@echo "================================================================================"
	cd $(BASE_DIR)
	rm -rf $(BUILD_DIR)/src/petsc
	mkdir -p $(BUILD_DIR)/src/petsc
	tar -xzf petsc-lite-$(BUILD_PETSC).tar.gz -C $(BUILD_DIR)/src/petsc --strip-components=1
	@echo "================================================================================"
	cd $(BUILD_DIR)/src/petsc; ./configure --prefix=$(BUILD_DIR) $(PETSC_OPTIONS) PETSC_DIR=$(BUILD_DIR)/src/petsc PETSC_ARCH=$(PETSC_ARCH)
	@echo "================================================================================"
	cd $(BUILD_DIR)/src/petsc; make PETSC_DIR=$(BUILD_DIR)/src/petsc PETSC_ARCH=$(PETSC_ARCH) all
	cd $(BUILD_DIR)/src/petsc; make PETSC_DIR=$(BUILD_DIR)/src/petsc PETSC_ARCH=$(PETSC_ARCH) install
	@echo "================================================================================"
	@echo "export PETSC_DIR=$(BUILD_DIR)" >> env_source.sh

petsc-lite-$(BUILD_PETSC).tar.gz:
	@echo "================================================================================"
	@echo "download PETSc Lite $(BUILD_PETSC)"
	cd $(BASE_DIR)
	wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${BUILD_PETSC}.tar.gz
	@echo "================================================================================"
	
nifti: nifticlib-2.0.0.tar.gz
	@echo "================================================================================"
	@echo "building NIFTI"
	@echo "configure with: $(NIFTI_OPTIONS)"
	@echo "================================================================================"
	cd $(BASE_DIR)
	rm -rf $(BUILD_DIR)/src/nifti
	mkdir -p $(BUILD_DIR)/src/nifti
	tar -xzf nifticlib-2.0.0.tar.gz -C $(BUILD_DIR)/src/nifti --strip-components=1
	@echo "================================================================================"
	cd $(BUILD_DIR)/src/nifti; cmake -DCMAKE_INSTALL_PREFIX=$(BUILD_DIR) $(NIFTI_OPTIONS)
	@echo "================================================================================"
	cd $(BUILD_DIR)/src/nifti; make
	cd $(BUILD_DIR)/src/nifti; make install
	@echo "================================================================================"
	@echo "export NIFTI_DIR=$(BUILD_DIR)" >> env_source.sh
	
nifticlib-2.0.0.tar.gz:
	@echo "================================================================================"
	@echo "download Nifti C Lib"
	cd $(BASE_DIR)
	wget http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz
	@echo "================================================================================"

.PHONY: config
