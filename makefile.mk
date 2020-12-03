#build code for GPUs (yes, no)
BUILD_GPU = yes
#build unit tests (yes, no)
BUILD_TEST = no
#build Python bindings (yes, no)
BUILD_PYTHON = no

#additional options
#enable nifti support (yes, no)
WITH_NIFTI = yes
#enable pnetcdf support (yes, no)
WITH_PNETCDF = no
#enable double precision (yes, no)
WITH_DOUBLE = no
#enable debugging information (yes, no)
WITH_DEBUG = no

WITH_DEBUG = no

WITH_DEVELOP = no

BUILD_TARGET = X86
GPU_VERSION = 35
CPP_VERSION = c++14

include config.mk
include filelist.mk

all: .PRIMARY $(PRESTEPS) $(BINS)
	@echo "================================================================================"
	@echo "writing cache"
	@echo $(CACHE) > make.cache
	@echo "done"
	@echo "================================================================================"

.PRIMARY: config
	@echo "Building"
	@echo "================================================================================"

$(BUILD_DIR)/%: $(OBJ_DIR)/$(APP_DIR)/%.o $(OBJS)
	-@$(MKDIR) $(dir $@)
ifdef VERBOSE
	$(CXX) $(addprefix -L,$(LIBRARIES)) $^ $(LDFLAGS) -o $@
else
	@echo linking $@
	@$(CXX) $(addprefix -L,$(LIBRARIES)) $^ $(LDFLAGS) -o $@
endif

$(LIB_DIR)/libclaire.so: $(CPU_OBJS) $(GPU_OBJS)
	-@$(MKDIR) $(dir $@)
ifdef VERBOSE
	$(CXX) $(addprefix -L,$(LIBRARIES)) $^ $(subst -lclaire,,$(LDFLAGS)) -shared -o $@
else
	@echo linking $@
	@$(CXX) $(addprefix -L,$(LIBRARIES)) $^ $(subst -lclaire,,$(LDFLAGS)) -shared -o $@
endif

$(LIB_DIR)/_pyclaire.so: $(LIB_DIR)/libclaire.so $(SWIG_OBJS)
	-@$(MKDIR) $(dir $@)
ifdef VERBOSE
	$(CXX) $(addprefix -L,$(LIBRARIES)) $^ $(subst -lclaire,,$(LDFLAGS)) -lclaire -shared -o $@
else
	@echo linking $@
	@$(CXX) $(addprefix -L,$(LIBRARIES)) $^ $(subst -lclaire,,$(LDFLAGS)) -lclaire -shared -o $@
endif

$(OBJ_DIR)/%.co: %.cu
	-@$(MKDIR) $(dir $@)
ifdef VERBOSE
	$(NVCC) -ccbin $(CXX) $(NVCCFLAGS) $(addprefix -I,$(INCLUDES)) $^ -o $@
else
	@echo building $@
	@$(NVCC) -ccbin $(CXX) $(NVCCFLAGS) $(addprefix -I,$(INCLUDES)) $^ -o $@
endif

$(OBJ_DIR)/%.o: %.cpp
	-@$(MKDIR) $(dir $@)
ifdef VERBOSE
	$(CXX) $(CXXFLAGS) $(addprefix -I,$(INCLUDES)) $^ -o $@
else
	@echo building $@
	@$(CXX) $(CXXFLAGS) $(addprefix -I,$(INCLUDES)) $^ -o $@
endif

%_wrap.swig.cpp: %.i
ifdef VERBOSE
	swig -c++ -python -cppext swig.cpp $(addprefix -I,$(INCLUDES)) -outdir $(dir $^) $^
	-@$(MKDIR) $(LIB_DIR)
	mv $(basename $^).py $(LIB_DIR)/
else
	@echo building SWIG interface
	@swig -c++ -python -cppext swig.cpp $(addprefix -I,$(INCLUDES)) -outdir $(dir $^) $^
	-@$(MKDIR) $(LIB_DIR)
	@mv $(basename $^).py $(LIB_DIR)/
endif 

config:
	@echo "================================================================================"
	@echo "Options"
	@echo "================================================================================"
	@echo "BUILD_GPU:    $(BUILD_GPU); [yes, no]"
	@echo "BUILD_TEST:   $(BUILD_TEST); [yes, no]"
	@echo "BUILD_PYTHON: $(BUILD_PYTHON); [yes, no]"
	@echo "================================================================================"
	@echo "WITH_NIFTI:   $(WITH_NIFTI); [yes, no]"
	@echo "WITH_PNETCDF: $(WITH_PNETCDF); [yes, no]"
	@echo "WITH_DOUBLE:  $(WITH_DOUBLE); [yes, no]"
	@echo "WITH_DEBUG:   $(WITH_DEBUG); [yes, no]"
	@echo "WITH_DEVELOP: $(WITH_DEVELOP); [yes, no]"
	@echo "================================================================================"
	@echo "BUILD_DIR:    $(BUILD_DIR)"
	@echo "================================================================================"
	@echo "CXX:           $(CXX)"
	@echo "NVCC:          $(NVCC)"
	@echo "================================================================================"
ifdef VERBOSE
	@echo "internal build options"
	@echo "================================================================================"
	@echo "BUILD_SHARED:  $(BUILD_SHARED); [yes, no]"
	@echo "BUILD_TARGET:  $(BUILD_TARGET); [POWER9, X86]"
	@echo "================================================================================"
	@echo "MPI_DIR:       $(MPI_DIR)"
	@echo "CUDA_DIR:      $(CUDA_DIR)"
	@echo "PETSC_DIR:     $(PETSC_DIR)"
	@echo "NIFTI_DIR:     $(NIFTI_DIR)"
	@echo "ZLIB_DIR:      $(ZLIB_DIR)"
	@echo "PNETCDF_DIR:   $(PNETCDF_DIR)"
	@echo "PYTHON_DIR:    $(PYTHON_DIR)"
	@echo "================================================================================"
	@echo "GPU_VERSION:   $(GPU_VERSION)"
	@echo "CPP_VERSION:   $(CPP_VERSION)"
	@echo "================================================================================"
endif
ifdef VVERBOSE
	@echo "APP_DIR:      $(APP_DIR)"
	@echo "SRC_DIR:      $(SRC_DIR)"
	@echo "OBJ_DIR:      $(OBJ_DIR)"
	@echo "LIB_DIR:      $(LIB_DIR)"
	@echo "EXSRC_DIR:    $(EXSRC_DIR)"
	@echo "================================================================================"
	@echo "CXX_FLAGS:    $(CXX_FLAGS)"
	@echo "NVCC_FLAGS:   $(NVCC_FLAGS)"
	@echo "LD_FLAGS:     $(LD_FLAGS)"
	@echo "================================================================================"
endif

clean:
	@echo "================================================================================"
	@echo "Cleaning up build"
	@echo "================================================================================"
	$(RM) -r $(LIB_DIR)
	$(RM) -r $(BUILD_DIR)
	$(RM) -r $(OBJ_DIR)
	$(RM) *~ */*~
	$(RM) make.cache
	@echo "================================================================================"

.PHONY: config clean
