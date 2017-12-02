CXX=mpicxx

USEINTEL=yes
USEINTELMPI=yes
USESINGLE=no
USEPNETCDF=yes
USENIFTI=yes
USEKNL=no
USEHASWELL=yes
BUILDTOOLS=yes


# developer flags (ignore)
DBGCODE=no
PEDANTIC=yes

RM = rm -f
MKDIRS = mkdir -p

ifeq ($(DBGCODE),yes)
	CXXFLAGS = -g -debug all
else
	CXXFLAGS = -O3 -ansi
endif


ifeq ($(USEINTEL),yes)
	CXXFLAGS += -xhost -parallel
	#CXXFLAGS += -openmp
	CXXFLAGS += -qopenmp
else
	CXXFLAGS += -fopenmp
endif
CXXFLAGS += -std=c++11


ifeq ($(USEKNL),yes)
	CXXFLAGS += -DKNL
endif

ifeq ($(USEHASWELL),yes)
	ifeq ($(USESINGLE),yes)
		CXXFLAGS += -DHASWELL
	endif
endif


ifeq ($(PEDANTIC),yes)
	CXXFLAGS += -Warray-bounds -Wchar-subscripts -Wcomment
	CXXFLAGS += -Wenum-compare -Wformat -Wuninitialized
	CXXFLAGS += -Wmaybe-uninitialized -Wmain -Wnarrowing
	CXXFLAGS += -Wnonnull -Wparentheses -Wpointer-sign
	CXXFLAGS += -Wreorder -Wreturn-type -Wsign-compare
	CXXFLAGS += -Wsequence-point -Wtrigraphs -Wunused-function
	CXXFLAGS += -Wunused-but-set-variable -Wunused-variable -Wwrite-strings
endif

ifeq ($(USEPNETCDF),yes)
	CXXFLAGS += -DREG_HAS_PNETCDF
endif

ifeq ($(USENIFTI),yes)
	CXXFLAGS += -DREG_HAS_NIFTI
endif

BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include
APPDIR = ./apps

#GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
GIT_VERSION := $(shell git describe --abbrev=4 --always --tags)
CXXFLAGS += -DGIT_VERSION=\"$(GIT_VERSION)\"

CLAIRE_INC = -I$(INCDIR)
ifeq ($(DBGCODE),yes)
	ifeq ($(USESINGLE),yes)
		CLAIRE_INC += -isystem$(PETSC_DIR)/include -isystem$(PETSC_DIR)/$(PETSC_ARCH_DBG_SINGLE)/include
	else 
		CLAIRE_INC += -isystem$(PETSC_DIR)/include -isystem$(PETSC_DIR)/$(PETSC_ARCH_DBG_DOUBLE)/include
	endif
else
	ifeq ($(USESINGLE),yes)
		CLAIRE_INC += -isystem$(PETSC_DIR)/include -isystem$(PETSC_DIR)/$(PETSC_ARCH_SINGLE)/include
	else
		CLAIRE_INC += -isystem$(PETSC_DIR)/include -isystem$(PETSC_DIR)/$(PETSC_ARCH_DOUBLE)/include
	endif
endif
CLAIRE_INC += -I$(ACCFFT_DIR)/include
CLAIRE_INC += -I$(FFTW_DIR)/include
CLAIRE_INC += -I$(MORTON_DIR)
CLAIRE_INC += -I./3rdparty
ifeq ($(USENIFTI),yes)
	CLAIRE_INC += -I$(NIFTI_DIR)/include/nifti
endif

ifeq ($(USEPNETCDF),yes)
	CLAIRE_INC += -I$(PNETCDF_DIR)/include
endif


ifeq ($(DBGCODE),yes)
	ifeq ($(USESINGLE),yes)
		LDFLAGS += -L$(PETSC_DIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH_DBG_SINGLE)/lib
	else
		LDFLAGS += -L$(PETSC_DIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH_DBG_DOUBLE)/lib
	endif
else
	ifeq ($(USESINGLE),yes)
		LDFLAGS += -L$(PETSC_DIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH_SINGLE)/lib
	else
		LDFLAGS += -L$(PETSC_DIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH_DOUBLE)/lib
	endif
endif
LDFLAGS += -lpetsc -lf2clapack -lf2cblas 

ifeq ($(USENIFTI),yes)
	LDFLAGS += -L$(NIFTI_DIR)/lib -lnifticdf -lniftiio -lznz -L$(ZLIB_DIR)/lib -lz
endif

#LDFLAGS+= -lcrypto -lssl -ldl
ifeq ($(USEINTEL),yes)
	LDFLAGS += -limf
endif


ifeq ($(USEINTELMPI),yes)
#	LDFLAGS += -lmpi_mt
endif
LDFLAGS += -lm


# FFT LIBRARIES
LDFLAGS += -L$(ACCFFT_DIR)/lib -laccfft -laccfft_utils
ifeq ($(USEPNETCDF),yes)
	LDFLAGS += -L$(PNETCDF_DIR)/lib -lpnetcdf
endif

ifeq ($(USESINGLE),yes)
	LDFLAGS += -L$(FFTW_DIR)/lib
endif
LDFLAGS += -L$(FFTW_DIR)/lib

ifeq ($(USESINGLE),yes)
	LDFLAGS += -lfftw3f_threads -lfftw3f
endif
LDFLAGS += -lfftw3_threads -lfftw3


BIN += $(BINDIR)/claire
ifeq ($(BUILDTOOLS),yes)
	BIN += $(BINDIR)/benchmark
	BIN += $(BINDIR)/clairetools
endif


CPPFILES=$(SRCDIR)/RegOpt.cpp \
		$(SRCDIR)/RegToolsOpt.cpp \
		$(SRCDIR)/CLPBenchmark.cpp \
		$(SRCDIR)/RegUtils.cpp \
		$(SRCDIR)/ghost.cpp \
		$(SRCDIR)/interp3.cpp \
		$(SRCDIR)/Interp3_Plan.cpp \
		$(SRCDIR)/VecField.cpp \
		$(SRCDIR)/TenField.cpp \
		$(SRCDIR)/ReadWriteReg.cpp \
		$(SRCDIR)/SynProbRegistration.cpp \
		$(SRCDIR)/SemiLagrangian.cpp \
		$(SRCDIR)/Optimizer.cpp \
		$(SRCDIR)/KrylovInterfaceReg.cpp \
		$(SRCDIR)/TaoInterfaceRegistration.cpp \
		$(SRCDIR)/RegistrationInterface.cpp \
		$(SRCDIR)/MultiLevelPyramid.cpp \
		$(SRCDIR)/PrecondReg.cpp \
		$(SRCDIR)/RegularizationRegistration.cpp \
		$(SRCDIR)/RegularizationRegistrationL2.cpp \
		$(SRCDIR)/RegularizationRegistrationH1.cpp \
		$(SRCDIR)/RegularizationRegistrationH2.cpp \
		$(SRCDIR)/RegularizationRegistrationH1SN.cpp \
		$(SRCDIR)/RegularizationRegistrationH2SN.cpp \
		$(SRCDIR)/RegularizationRegistrationH3.cpp \
		$(SRCDIR)/RegularizationRegistrationH3SN.cpp \
		$(SRCDIR)/OptimizationProblem.cpp \
		$(SRCDIR)/OptimalControlRegistrationBase.cpp \
		$(SRCDIR)/OptimalControlRegistration.cpp \
		$(SRCDIR)/OptimalControlRegistrationIC.cpp \
		$(SRCDIR)/OptimalControlRegistrationRelaxedIC.cpp \
		$(SRCDIR)/Preprocessing.cpp

OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(CPPFILES))

.SECONDARY: $(OBJS)

all: $(BIN)

$(BINDIR)/%: $(OBJDIR)/%.o $(OBJS)
	-@$(MKDIRS) $(dir $@) # if bin exists dont give an error
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) $^ $(LDFLAGS) $(CLAIRE_LIB) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) -c $^ -o $@

$(OBJDIR)/%.o: $(APPDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) -c $^ -o $@

.PHONY: clean

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) -r $(BINDIR) $(OBJDIR) results
	$(RM) *~ */*~
