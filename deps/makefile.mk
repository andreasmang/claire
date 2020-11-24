BUILD_PETSC   = 3.12.4
BUILD_NIFTI   = yes
BUILD_PNETCDF = no
BUILD_ACCFFT  = no

BUILD_GPU     = yes
BUILD_DEBUG   = no
BUILD_SHARED  = no

BUILD_DIR = $(PWD)/lib

CC = mpicc
CXX = mpicxx
NVCC = nvcc

ifneq ($(BUILD_PETSC), no)
	TARGETS += petsc
endif

ifneq ($(BUILD_NIFTI), no)
	TARGETS += nifti
endif

PETSC_OPTIONS += --download-f2cblaslapack
PETSC_OPTIONS += --with-64-bit-indices
PETSC_OPTIONS += --with-x=0
PETSC_OPTIONS += --with-fc=0
PETSC_OPTIONS += --with-ssl=0
PETSC_OPTIONS += --COPTFLAGS='-O3'
PETSC_OPTIONS += --CXXOPTFLAGS='-O3'
PETSC_OPTIONS += --with-precision=single
PETSC_OPTIONS += --with-cc=$(CC)
PETSC_OPTIONS += --with-cxx=$(CXX)

NIFTI_OPTIONS += -DCMAKE_CXX_COMPILER=$(CXX)
NIFTI_OPTIONS += -DCMAKE_C_COMPILER=$(CC)
NIFTI_OPTIONS += -Wno-dev

ifeq ($(BUILD_SHARED), yes)
	PETSC_OPTIONS += --with-shared=1
	NIFTI_OPTIONS += -DBUILD_SHARED_LIBS:BOOL=ON
else
	PETSC_OPTIONS += --with-shared=0
	NIFTI_OPTIONS += -DBUILD_SHARED_LIBS:BOOL=OFF
endif

ifeq ($(BUILD_GPU), yes)
	PETSC_OPTIONS += --with-cuda=1 
	PETSC_OPTIONS += --download-cusp=yes
	PETSC_OPTIONS += --CUDAOPTFLAGS='-O3'
	PETSC_OPTIONS += --with-cudac=$(NVCC)
else
	PETSC_OPTIONS += --with-cuda=0
endif

ifeq ($(BUILD_DEBUG), yes)
	PETSC_OPTIONS += --with-debugging=1
else
	PETSC_OPTIONS += --with-debugging=0
endif

BASE_DIR=$(PWD)

all: config $(TARGETS)
	@echo "================================================================================"
	@echo "done"
	@echo "================================================================================"

config:
	@echo "================================================================================"
	@echo "options"
	@echo "================================================================================"
	@echo "build PETSc:      $(BUILD_PETSC)"
	@echo "================================================================================"
	@echo "build for GPU:    $(BUILD_GPU)"
	@echo "build with DEBUG: $(BUILD_DEBUG)"
	@echo "build shared:     $(BUILD_SHARED)"
	@echo "build directory:  $(BUILD_DIR)"
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
	cd $(BUILD_DIR)/src/petsc; ./configure --prefix=$(BUILD_DIR) $(PETSC_OPTIONS)
	@echo "================================================================================"
	cd $(BUILD_DIR)/src/petsc; make
	cd $(BUILD_DIR)/src/petsc; make install
	@echo "================================================================================"

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
	
nifticlib-2.0.0.tar.gz:
	@echo "================================================================================"
	@echo "download Nifti C Lib"
	cd $(BASE_DIR)
	wget http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz
	@echo "================================================================================"

.PHONY: config
