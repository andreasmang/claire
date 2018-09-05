CXX=mpicxx
USECUDA=yes
USEINTEL=no
USEINTELMPI=no
USESINGLE=yes
USEPNETCDF=yes
USENIFTI=yes
USEKNL=no
USEHASWELL=no
BUILDTOOLS=yes


include config/setup.mk
include config/files.mk


ifeq ($(USECUDA),yes)
CUDAC=nvcc
CUDA_OBJS = $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/cuda/%.o,$(CUFILES))
CUDA_OBJS += $(patsubst $(EXSRCDIR)/%.cu,$(OBJDIR)/cuda/%.o,$(EXCUFILES))
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

$(OBJDIR)/cuda/%.o: $(SRCDIR)/%.cu
	-@$(MKDIRS) $(dir $@)
	$(CUDAC) $(CUDA_FLAGS) $(CUDA_INC) -c $^ -o $@
	
$(OBJDIR)/cuda/%.o: $(EXSRCDIR)/%.cu
	-@$(MKDIRS) $(dir $@)
	$(CUDAC) $(CUDA_FLAGS) $(CUDA_INC) -c $^ -o $@

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
