CXX=mpicxx

USEINTEL=yes
USEINTELMPI=yes
BUILDTOOLS=yes
DBGCODE=yes
PEDANTIC=yes

RM = rm -f
MKDIRS = mkdir -p

ifeq ($(DBGCODE),yes)
	CXXFLAGS = -g -debug all
else
	CXXFLAGS = -O3 -ansi
endif

ifeq ($(USEINTEL),yes)
	CXXFLAGS+= -std=c++11 -DINVERT_RHO -xhost -parallel
	CXXFLAGS+= -openmp
else
	CXXFLAGS+= -fopenmp
endif

ifeq ($(PEDANTIC),yes)
	CXXFLAGS+= -Warray-bounds -Wchar-subscripts -Wcomment
	CXXFLAGS+= -Wenum-compare -Wformat -Wuninitialized
	CXXFLAGS+= -Wmaybe-uninitialized -Wmain -Wnarrowing
	CXXFLAGS+= -Wnonnull -Wparentheses -Wpointer-sign
	CXXFLAGS+= -Wreorder -Wreturn-type -Wsign-compare
	CXXFLAGS+= -Wsequence-point -Wtrigraphs -Wunused-function
	CXXFLAGS+= -Wunused-but-set-variable -Wunused-variable -Wwrite-strings
endif


BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include
APPDIR = ./apps

COLD_INC = -I$(INCDIR)
ifeq ($(DBGCODE),yes)
	COLD_INC+= -isystem$(PETSC_DBG_DIR)/include -isystem$(PETSC_DBG_DIR)/$(PETSC_DBG_ARCH)/include
else
	COLD_INC+= -isystem$(PETSC_DIR)/include -isystem$(PETSC_DIR)/$(PETSC_ARCH)/include
endif
COLD_INC+= -I$(ACCFFT_DIR)/include
COLD_INC+= -I$(FFTW_DIR)/include
COLD_INC+= -I$(NIFTI_DIR)/include/nifti

LDFLAGS+= -L$(ACCFFT_DIR)/lib -laccfft -laccfft_utils
LDFLAGS+= -L$(FFTW_DIR)/lib -lfftw3 -lfftw3_threads
ifeq ($(DBGCODE),yes)
	LDFLAGS+= -L$(PETSC_DBG_DIR)/lib -L$(PETSC_DBG_DIR)/$(PETSC_DBG_ARCH)/lib -lpetsc -lf2clapack -lf2cblas
else
	LDFLAGS+= -L$(PETSC_DIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lf2clapack -lf2cblas
endif
LDFLAGS+= -L$(NIFTI_DIR)/lib -lnifticdf -lniftiio -lznz
LDFLAGS+= -L$(ZLIB_DIR)/lib -lz
#LDFLAGS+= -lcrypto -lssl -ldl
ifeq ($(USEINTEL),yes)
	LDFLAGS+= -limf
endif


ifeq ($(USEINTELMPI),yes)
#	LDFLAGS+= -lmpi_mt
endif
LDFLAGS+= -lm

BIN+=$(BINDIR)/runcoldreg
BIN+=$(BINDIR)/regtools
ifeq ($(BUILDTOOLS),yes)
#	BIN+=$(BINDIR)/regtools
#	BIN+= $(BINDIR)/par_interp3_driver
endif


CPPFILES=$(SRCDIR)/RegOpt.cpp \
		$(SRCDIR)/RegToolsOpt.cpp \
		$(SRCDIR)/RegUtils.cpp \
		$(SRCDIR)/ghost.cpp \
		$(SRCDIR)/interp3.cpp \
		$(SRCDIR)/Interp3_Plan.cpp \
		$(SRCDIR)/VecField.cpp \
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
		$(SRCDIR)/RegularizationRegistrationH1.cpp \
		$(SRCDIR)/RegularizationRegistrationH2.cpp \
		$(SRCDIR)/RegularizationRegistrationH3.cpp \
		$(SRCDIR)/RegularizationRegistrationH1SN.cpp \
		$(SRCDIR)/RegularizationRegistrationH2SN.cpp \
		$(SRCDIR)/RegularizationRegistrationH3SN.cpp \
		$(SRCDIR)/OptimizationProblem.cpp \
		$(SRCDIR)/OptimalControlRegistrationBase.cpp \
		$(SRCDIR)/OptimalControlRegistration.cpp \
		$(SRCDIR)/OptimalControlRegistrationIC.cpp \
		$(SRCDIR)/OptimalControlRegistrationRelaxedIC.cpp \
		$(SRCDIR)/PreProcReg.cpp

OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(CPPFILES))

.SECONDARY: $(OBJS)

all: $(BIN)

$(BINDIR)/%: $(OBJDIR)/%.o $(OBJS)
	-@$(MKDIRS) $(dir $@) # if bin exists dont give an error
	$(CXX) $(CXXFLAGS) $(COLD_INC) $^ $(LDFLAGS) $(COLD_LIB) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(COLD_INC) -c $^ -o $@

$(OBJDIR)/%.o: $(APPDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(COLD_INC) -c $^ -o $@

.PHONY: clean

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) -r $(BINDIR) $(OBJDIR) results
	$(RM) *~ */*~
