#directories
override APP_DIR := ./apps
override SRC_DIR := ./src
override OBJ_DIR := ./obj
override LIB_DIR := ./lib
override EXSRC_DIR := ./3rdparty
INCLUDES = ./include ./deps/3rdparty

BUILD_DIR = ./bin

INCLUDES += $(EXSRC_DIR)

#programs
CXX = mpicxx
NVCC = nvcc
RM = rm -f
MKDIR = mkdir -p

LIBRARIES = $(LIB_DIR)
LDFLAGS += -fopenmp

#locate dependencies
MPI_DIR = $(abspath $(subst bin,,$(dir $(shell which $(CXX)))))
INCLUDES += $(MPI_DIR)/include
INCLUDES += $(MPI_DIR)/include/mpi
LIBRARIES += $(MPI_DIR)/lib

CXXFLAGS += -g $(CXX_FLAGS)
LDFLAGS += -g $(LD_FLAGS)

ifeq ($(WITH_DEVELOP), yes)
	CXXFLAGS += -DZEITGEIST
endif

ifeq ($(WITH_DEVELOP_EXTRA), yes)
	CXXFLAGS += -DZEITGEIST_DEV
endif

ifeq ($(WITH_CUDA_MPI), yes)
	CXXFLAGS += -DREG_HAS_MPICUDA
endif

ifeq ($(BUILD_TARGET), POWER9)
	CXXFLAGS += -DREG_HAS_MPICUDA
	CXXFLAGS += -DPOWER9
endif

ifeq ($(BUILD_GPU), yes)
	ifeq ($(WITH_DOUBLE), yes)
$(error GPU build only supports single precision)
	endif
	CUDA_DIR = $(abspath $(subst bin,,$(dir $(shell which $(NVCC)))))
	INCLUDES += $(CUDA_DIR)/include
	LIBRARIES += $(CUDA_DIR)/lib
	LIBRARIES += $(CUDA_DIR)/lib64
	LDFLAGS += -lcuda -lcudart -lcufft -lcublas -lcusparse -lcusolver
	NVCCFLAGS += -Xcompiler "$(CXXFLAGS)" --std=$(CPP_VERSION) -O3 -c $(NVCC_FLAGS)
	CXXFLAGS += -DREG_HAS_CUDA
	CXXFLAGS += -DREG_FFT_CUDA
	ifeq ($(WITH_DEBUG),yes)
		CXXFLAGS += -DREG_DBG_CUDA
		LDFLAGS += -lnvToolsExt
	endif
	ifdef GPU_VERSION
		NVCCFLAGS += -gencode arch=compute_$(GPU_VERSION),code=sm_$(GPU_VERSION)
	endif 
else
$(error This branch only supports GPU build)
endif

PETSC_DIR ?= ./deps/lib
INCLUDES += $(PETSC_DIR)/include
LIBRARIES += $(PETSC_DIR)/lib
LDFLAGS += -lpetsc -lf2clapack -lf2cblas

#python build needs shared library
PYTHON_DIR = /usr/include/python3.5
BUILD_SHARED ?= no
ifeq ($(BUILD_PYTHON),yes)
	override BUILD_SHARED = yes
	INCLUDES += $(PYTHON_DIR)
endif

GIT_VERSION := $(shell git describe --abbrev=4 --always --tags)
CXXFLAGS += -DGIT_VERSION=$(GIT_VERSION)

CXXFLAGS += -Wall -Wno-unused-function -c -std=$(CPP_VERSION)

ifeq ($(WITH_DEBUG), yes)
	CXXFLAGS += -g
else
	CXXFLAGS += -O3
	#CXXFLAGS += -march=native
	CXXFLAGS += -fopenmp
endif

#link ACCFFT and FFTW if not building for GPUs
ifneq ($(BUILD_GPU), yes)
	FFTW_DIR ?= ./
	ACCFFT_DIR ?= ./
	INCLUDES += $(FFTW_DIR)/include
	INCLUDES += $(ACCFFT_DIR)/include
	LIBRARIES += $(FFTW_DIR)/lib
	LIBRARIES += $(ACCFFT_DIR)/lib
	LDFLAGS += -lfftw3_threads -lfftw3 -laccfft -laccfft_utils
	ifneq ($(WITH_DOUBLE), yes)
		LDFLAGS += -lfftw3f_threads -lfftw3f
	endif
endif

ifeq ($(BUILD_SHARED), yes)
	CXXFLAGS += -fPIC
	#LDFLAGS += -shared
endif

ifeq ($(WITH_NIFTI), yes)
	CXXFLAGS += -DREG_HAS_NIFTI
	NIFTI_DIR ?= ./deps/lib
	ZLIB_DIR ?= ./
	INCLUDES += $(NIFTI_DIR)/include/nifti
	LIBRARIES += $(NIFTI_DIR)/lib
	LIBRARIES += $(ZLIB_DIR)/lib
	LDFLAGS += -lnifticdf -lniftiio -lznz -lz
endif

ifeq ($(WITH_PNETCDF), yes)
	CXXFLAGS += -DREG_HAS_PNETCDF
	PNETCDF_DIR ?= ./
	INCLUDES += $(PNETCDF_DIR)/include
	LIBRARIES += $(PNETCDF_DIR)/lib
	LDFLAGS += -lpnetcdf
endif

ifeq ($(BUILD_SHARED),yes)
	LDFLAGS += -lclaire
endif

CACHE = "$(BUILD_GPU) $(BUILD_TEST) $(BUILD_PYTHON) $(WITH_NIFTI) $(WITH_PNETCDF) $(WITH_DOUBLE) $(WITH_DEBUG) $(BUILD_SHARED) $(WITH_DEBUG) $(BUILD_TARGET) $(GPU_VERSION) $(shell uname -a)"

OLD_CACHE = "$(shell touch make.cache; cat make.cache)"

ifneq ($(strip $(OLD_CACHE)), $(strip $(CACHE)))
	PRESTEPS += clean
endif
