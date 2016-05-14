CXX=mpicxx
USEINTEL=yes
USEINTELMPI=no
RM = rm -f
MKDIRS = mkdir -p

CXXFLAGS = -O3 -ansi
ifeq ($(USEINTEL),yes)
	CXXFLAGS+= -openmp -std=c++11 -DINVERT_RHO -xhost -parallel
else
	CXXFLAGS+= -fopenmp
endif

BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include
APPDIR = ./apps

COLD_INC = -I$(INCDIR)
COLD_INC+= -isystem$(PETSC_DIR)/include -isystem$(PETSC_DIR)/$(PETSC_ARCH)/include
COLD_INC+= -I$(ACCFFT_DIR)/include
COLD_INC+= -I$(PNETCDF_DIR)/include
COLD_INC+= -I$(FFTW_DIR)/include
COLD_INC+= -I$(NIFTI_DIR)/include/nifti

LDFLAGS+= -L$(ACCFFT_DIR)/lib -laccfft -laccfft_utils
LDFLAGS+= -L$(FFTW_DIR)/lib -lfftw3 -lfftw3_threads
LDFLAGS+= -L$(PETSC_DIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
LDFLAGS+= -L$(NIFTI_DIR)/lib -lnifticdf -lniftiio -lznz
ifeq ($(USEINTEL),yes)
	LDFLAGS+= -limf
endif

ifeq ($(USEINTELMPI),yes)
	LDFLAGS+= -lmpi_mt
endif
LDFLAGS+= -lm

BIN+= $(BINDIR)/runcoldreg


INCFILES=RegOpt.h RegUtils.h interp3.hpp utils.hpp interp3_common.hpp VecField.h ReadWriteReg.h SynProbRegistration.h SemiLagrangian.h Optimizer.h TaoInterfaceRegistration.h RegularizationRegistration.h LargeDeformationRegistration.h OptimalControlRegistration.h OptimalControlRegistrationIC.h OptProbRegistration.h PreProcessingRegistration.h
DEPS = $(patsubst %,$(INCDIR)/%.hpp,$(INCFILES))


CPPFILES=$(SRCDIR)/RegOpt.cpp \
		$(SRCDIR)/RegUtils.cpp \
		$(SRCDIR)/ghost.cpp \
		$(SRCDIR)/interp3.cpp \
		$(SRCDIR)/Interp3_Plan.cpp \
		$(SRCDIR)/VecField.cpp \
		$(SRCDIR)/ReadWriteReg.cpp \
		$(SRCDIR)/SynProbRegistration.cpp \
		$(SRCDIR)/SemiLagrangian.cpp \
		$(SRCDIR)/Optimizer.cpp \
		$(SRCDIR)/TaoInterfaceRegistration.cpp \
		$(SRCDIR)/RegularizationRegistration.cpp \
		$(SRCDIR)/RegularizationRegistrationH1.cpp \
		$(SRCDIR)/RegularizationRegistrationH2.cpp \
		$(SRCDIR)/RegularizationRegistrationH1SN.cpp \
		$(SRCDIR)/RegularizationRegistrationH2SN.cpp \
		$(SRCDIR)/LargeDeformationRegistration.cpp \
		$(SRCDIR)/OptimalControlRegistration.cpp \
		$(SRCDIR)/OptimalControlRegistrationIC.cpp \
		$(SRCDIR)/OptimalControlRegistrationRIC.cpp \
		$(SRCDIR)/OptProbRegistration.cpp \
		$(SRCDIR)/PreProcessingRegistration.cpp

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
