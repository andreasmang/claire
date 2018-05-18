CXX=$(MPICXX)/mpicxx
USECUDA=no
ifeq ($(USECUDA),yes)
CUDAC=$(CUDA_DIR)/bin/nvcc
endif
USEINTEL=no
USEINTELMPI=no
USESINGLE=yes
USEPNETCDF=yes
USENIFTI=no
USEKNL=no
USEHASWELL=yes
BUILDTOOLS=yes


include config/setup.mk
include config/files.mk


ifeq ($(USECUDA),yes)
CUDA_OBJS = $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(CUFILES))
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(CPPFILESCUDA))
else
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(CPPFILES))
endif

.SECONDARY: $(OBJS)

all: $(BIN)

ifeq ($(USECUDA),yes)

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

else



$(BINDIR)/%: $(OBJDIR)/%.o $(OBJS)
	-@$(MKDIRS) $(dir $@) # if bin exists dont give an error
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) $^ $(LDFLAGS) $(CLAIRE_LIB) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) -c $^ -o $@

$(OBJDIR)/%.o: $(APPDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CLAIRE_INC) -c $^ -o $@

endif

.PHONY: clean

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) -r $(BINDIR) $(OBJDIR) results
	$(RM) *~ */*~
