CXX=mpicxx            # select the MPI compiler on your system
USEINTEL=no           # enable if you want plan on using intel compilers
USEINTELMPI=no        # enable if you plan on using INTEL MPI
USESINGLE=yes         # compile CLAIRE in single precision (needs PETSC in single precision)
USEPNETCDF=no         # use PNETCDF as data format (do not forget to build the dependency)
USENIFTI=yes          # use NIFTI as data format; this is the default for IO and should always be enabled
USEHASWELL=no         # build CLAIRE for a HASWELL architecture
BUILDTOOLS=yes        # build the `clairetools` binary

include config/setup.mk
include config/files.mk

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
