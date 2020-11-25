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
MPI_DIR = /usr/lib/openmpi
INCLUDES += $(MPI_DIR)/include
LIBRARIES += $(MPI_DIR)/lib

ifeq ($(BUILD_GPU), yes)
	CUDA_DIR = $(abspath $(subst bin,,$(dir $(shell which $(NVCC)))))
	INCLUDES += $(CUDA_DIR)/include
	LIBRARIES += $(CUDA_DIR)/lib
	LIBRARIES += $(CUDA_DIR)/lib64
	LDFLAGS += -lcuda -lcudart -lcufft -lcublas -lcusparse -lcusolver
	NVCCFLAGS += -Xcompiler "$(CXXFLAGS)" --std=c++11 -O3 -c
	#NVCCFLAGS += -gencode arch=compute_35,code=sm_35
	CXXFLAGS += -DREG_HAS_CUDA
	CXXFLAGS += -DREG_FFT_CUDA
	ifeq ($(WITH_DEBUG),yes)
		CXXFLAGS += -DREG_DBG_CUDA
		LDFLAGS += -lnvToolsExt
	endif
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

CXXFLAGS += -Wall -c -std=c++11

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
