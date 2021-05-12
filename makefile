# select MPI compiler
CXX=mpicxx
# use intel compilder (yes,no) 
USEINTEL=no
# use intel MPI (yes,no) 
USEINTELMPI=no
# use single precision (yes,no) 
USESINGLE=yes
# use PNETCDF data format (yes,no; requires PNETCDF lib; see dependencies) 
USEPNETCDF=no
# use NIFTI data format (yes,no; default IO format; see dependencies)
USENIFTI=yes
# build CLAIRE for a HASWELL architecture (yes,no)
USEHASWELL=no
# build `clairetools` binary (yes,no)
BUILDTOOLS=yes

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
