CXX=mpicxx

USEINTEL=no
USEINTELMPI=no
USESINGLE=no
USEPNETCDF=yes
USENIFTI=yes
USEKNL=no
USEHASWELL=yes
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
