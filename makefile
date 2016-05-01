CXX=mpicxx
CXXFLAGS= -O3 -ansi -openmp -xhost -DINVERT_RHO -std=c++11 #-Wfatal-errors -Wall -Wextra -Wconversion -Wshadow
RM = rm -f
MKDIRS = mkdir -p

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
LDFLAGS+= -L$(PNETCDF_DIR)/lib -lpnetcdf
LDFLAGS+= -L$(NIFTI_DIR)/lib -lnifticdf -lniftiio -lznz
LDFLAGS+= -limf -lm

BIN = $(BINDIR)/runcoldreg

SRCFILES=$(SRCDIR)/RegOpt.cpp \
		$(SRCDIR)/RegUtils.cpp \
		$(SRCDIR)/ghost.cpp \
		$(SRCDIR)/interp3.cpp \
		$(SRCDIR)/Interp3_Plan.cpp \
		$(SRCDIR)/VecField.cpp \
		$(SRCDIR)/DataReadWriteRegistration.cpp \
		$(SRCDIR)/SynProbRegistration.cpp \
		$(SRCDIR)/SemiLagrangian.cpp \
		$(SRCDIR)/Optimizer.cpp \
		$(SRCDIR)/TaoInterfaceRegistration.cpp \
		$(SRCDIR)/RegularizationRegistration.cpp \
		$(SRCDIR)/LargeDeformationRegistration.cpp \
		$(SRCDIR)/OptimalControlRegistration.cpp \
		$(SRCDIR)/OptimalControlRegistrationIC.cpp \
		$(SRCDIR)/OptimizationProblemRegistration.cpp \
		$(SRCDIR)/PreProcessingRegistration.cpp



OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES))

.SECONDARY: $(OBJS)

all: $(BIN)

$(BINDIR)/%: $(OBJDIR)/%.o $(OBJS)
	-@$(MKDIRS) $(dir $@) # if bin exists dont give an error
	$(CXX) $(CXXFLAGS) $(COLD_INC) $^ $(LDFLAGS) ${COLD_LIB} -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(COLD_INC) -c $^ -o $@



$(OBJDIR)/%.o: $(APPDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(COLD_INC) -c $^ -o $@

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) *~ */*~
