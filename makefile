CXX=$(MPICXX)/mpicxx
CUDAC=$(CUDA_DIR)/bin/nvcc
USEINTEL=no
USEINTELMPI=no
USECUDA=yes
USESINGLE=yes
USEPNETCDF=yes
USENIFTI=no
USEKNL=no
USEHASWELL=yes
BUILDTOOLS=yes


include config/setup.mk
include config/files.mk

OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(CPPFILES))
CUDA_OBJS = $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(CUFILES))


.SECONDARY: $(OBJS)

all: $(BIN)

$(BINDIR)/%: $(OBJDIR)/%.o $(CUDA_OBJS) $(OBJS)
	-@$(MKDIRS) $(dir $@) # if bin exists dont give an error
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) $^ $(LDFLAGS) $(CLAIRE_LIB) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) -c $^ -o $@

$(CUDA_OBJS): $(CUFILES)
	-@$(MKDIRS) $(dir $@)
	$(CUDAC) $(CUDA_FLAGS) -I$(CUDA_INC) -c $^ -o $@

$(OBJDIR)/%.o: $(APPDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) -c $^ -o $@

.PHONY: clean

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) -r $(BINDIR) $(OBJDIR) results
	$(RM) *~ */*~
