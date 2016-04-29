CXX=mpicxx
RM = rm -f
MKDIRS = mkdir -p

BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include ./apps -I$(ACCFFT_DIR)/include -I$(PNETCDF_DIR)/include -I$(FFTW_DIR)/include
APPDIR = ./apps

PSC_INC = -isystem$(PETSC_DIR)/include -isystem$(PETSC_DIR)/$(PETSC_ARCH)/include
PSC_LIB = -L$(PETSC_DIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

CXXFLAGS= -O3 -ansi -openmp -xhost -DINVERT_RHO

LDFLAGS=  -L$(FFTW_DIR)/lib  -lfftw3 -lfftw3_threads -L$(ACCFFT_DIR)/lib
LDFLAGS+= -L$(ACCFFT_DIR)/lib -laccfft -laccfft_utils -lfftw3 -lfftw3_threads -L$(PNETCDF_DIR)/lib -limf -lm -lpnetcdf
TARGET_BIN = $(BINDIR)/runcoldreg


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

#		$(SRCDIR)/RegCommandOptions.cpp \

OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES))

.SECONDARY: $(OBJS)

all: $(TARGET_BIN)

$(BINDIR)/%: $(OBJDIR)/%.o $(OBJS)
	-@$(MKDIRS) $(dir $@) # if bin exists dont give an error
	$(CXX) $(CXXFLAGS) ${PSC_INC} $^ $(LDFLAGS) ${PSC_LIB} -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) ${PSC_INC} -I$(INCDIR) -c $^ -o $@

$(OBJDIR)/%.o: $(APPDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) ${PSC_INC} -I$(INCDIR) -c $^ -o $@

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) *~ */*~
