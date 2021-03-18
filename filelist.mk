CPU_FILES=
GPU_FILES=
SWIG_FILES=
APP_FILES=

#general C++ files
CPU_FILES += $(SRC_DIR)/RegOpt.cpp
CPU_FILES += $(SRC_DIR)/RegToolsOpt.cpp
CPU_FILES += $(SRC_DIR)/BenchmarkOpt.cpp
CPU_FILES += $(SRC_DIR)/CLAIREUtils.cpp
CPU_FILES += $(SRC_DIR)/CLAIREUtilsKernel.cpp
CPU_FILES += $(SRC_DIR)/GhostPlan.cpp
CPU_FILES += $(SRC_DIR)/Differentiation/Differentiation.cpp
CPU_FILES += $(SRC_DIR)/Differentiation/DifferentiationFD.cpp
CPU_FILES += $(SRC_DIR)/Differentiation/DifferentiationSM.cpp
CPU_FILES += $(SRC_DIR)/ScaField.cpp
CPU_FILES += $(SRC_DIR)/VecField.cpp
CPU_FILES += $(SRC_DIR)/TenField.cpp
CPU_FILES += $(SRC_DIR)/ReadWriteReg.cpp
CPU_FILES += $(SRC_DIR)/SynProbRegistration.cpp
CPU_FILES += $(SRC_DIR)/Solver/TransportProblem.cpp
CPU_FILES += $(SRC_DIR)/Solver/TransportEquationSL.cpp
CPU_FILES += $(SRC_DIR)/Solver/TransportEquationRK2.cpp
CPU_FILES += $(SRC_DIR)/Solver/ContinuityEquation.cpp
CPU_FILES += $(SRC_DIR)/DeformationFields/DeformationFields.cpp
CPU_FILES += $(SRC_DIR)/DistanceMeasure/DistanceMeasure.cpp
CPU_FILES += $(SRC_DIR)/DistanceMeasure/DistanceMeasureNCC.cpp
CPU_FILES += $(SRC_DIR)/DistanceMeasure/DistanceMeasureSL2.cpp
CPU_FILES += $(SRC_DIR)/DistanceMeasure/DistanceMeasureSL2aux.cpp
CPU_FILES += $(SRC_DIR)/Optimizer.cpp
CPU_FILES += $(SRC_DIR)/KrylovInterface.cpp
CPU_FILES += $(SRC_DIR)/TaoInterface.cpp
CPU_FILES += $(SRC_DIR)/CLAIREInterface.cpp
CPU_FILES += $(SRC_DIR)/MultiLevelPyramid.cpp
CPU_FILES += $(SRC_DIR)/Preconditioner.cpp
CPU_FILES += $(SRC_DIR)/Regularization/Regularization.cpp
CPU_FILES += $(SRC_DIR)/Regularization/RegularizationL2.cpp
CPU_FILES += $(SRC_DIR)/Regularization/RegularizationH1.cpp
CPU_FILES += $(SRC_DIR)/Regularization/RegularizationH2.cpp
CPU_FILES += $(SRC_DIR)/Regularization/RegularizationH1SN.cpp
CPU_FILES += $(SRC_DIR)/Regularization/RegularizationH2SN.cpp
CPU_FILES += $(SRC_DIR)/Regularization/RegularizationH3.cpp
CPU_FILES += $(SRC_DIR)/Regularization/RegularizationH3SN.cpp
CPU_FILES += $(SRC_DIR)/OptimizationProblem.cpp
CPU_FILES += $(SRC_DIR)/CLAIREBase.cpp
CPU_FILES += $(SRC_DIR)/CLAIRE.cpp
CPU_FILES += $(SRC_DIR)/CLAIREStokes.cpp
CPU_FILES += $(SRC_DIR)/CLAIREDivReg.cpp
CPU_FILES += $(SRC_DIR)/Preprocessing.cpp
CPU_FILES += $(SRC_DIR)/Spectral/Spectral.cpp
CPU_FILES += $(SRC_DIR)/TwoLevel/TwoLevelFFT.cpp
CPU_FILES += $(SRC_DIR)/TwoLevel/TwoLevelFinite.cpp

ifneq ($(BUILD_GPU),yes)
	# CPU build specific C++ files
	CPU_FILES += $(SRC_DIR)/ghost.cpp
	CPU_FILES += $(SRC_DIR)/Interpolation/interp3.cpp
	CPU_FILES += $(SRC_DIR)/Interpolation/Interp3_Plan.cpp
	CPU_FILES += $(SRC_DIR)/Differentiation/DifferentiationKernel.cpp
	CPU_FILES += $(SRC_DIR)/Solver/TransportKernel.cpp
	CPU_FILES += $(SRC_DIR)/DeformationFields/DeformationKernel.cpp
	CPU_FILES += $(SRC_DIR)/DistanceMeasure/DistanceMeasureKernel.cpp
	CPU_FILES += $(SRC_DIR)/SemiLagrangian/SemiLagrangian.cpp
	CPU_FILES += $(SRC_DIR)/PreconditionerKernel.cpp
	CPU_FILES += $(SRC_DIR)/Regularization/RegularizationKernel.cpp
	CPU_FILES += $(SRC_DIR)/Spectral/SpectralKernel.cpp
else
	# GPU build specific C++ files
	CPU_FILES += $(SRC_DIR)/SemiLagrangian/SemiLagrangianGPUNew.cpp
	CPU_FILES += $(SRC_DIR)/Spectral/mpicufft.cpp
	# GPU build specific CUDA files
	GPU_FILES += $(SRC_DIR)/Interpolation/interp3_gpu_new.cu
	GPU_FILES += $(SRC_DIR)/Solver/TransportKernel.cu
	GPU_FILES += $(SRC_DIR)/CLAIREUtilsKernel.cu
	GPU_FILES += $(SRC_DIR)/DeformationFields/DeformationKernel.cu
	GPU_FILES += $(SRC_DIR)/DistanceMeasure/DistanceMeasureKernel.cu
	GPU_FILES += $(SRC_DIR)/Differentiation/DifferentiationKernel.cu
	GPU_FILES += $(SRC_DIR)/Differentiation/TextureDifferentiationKernel.cu
	GPU_FILES += $(SRC_DIR)/Spectral/SpectralKernel.cu
	GPU_FILES += $(SRC_DIR)/PreconditionerKernel.cu
	GPU_FILES += $(SRC_DIR)/Regularization/RegularizationKernel.cu
	GPU_FILES += $(SRC_DIR)/Interpolation/Interp3_Plan_GPU.cu
	GPU_FILES += $(SRC_DIR)/Interpolation/Interp3_Plan_GPU_kernel.cu
	GPU_FILES += $(SRC_DIR)/SemiLagrangian/SemiLagrangianKernel.cu
	GPU_FILES += $(SRC_DIR)/TwoLevel/TwoLevelKernel.cu
endif

APP_FILES += $(APP_DIR)/claire.cpp
APP_FILES += $(APP_DIR)/clairetools.cpp

ifeq ($(BUILD_TEST),yes)
	CPU_FILES += $(SRC_DIR)/UnitTestOpt.cpp
	CPU_FILES += $(SRC_DIR)/UnitTests/TestRegularization.cpp
	CPU_FILES += $(SRC_DIR)/UnitTests/TestDifferentiation.cpp
	CPU_FILES += $(SRC_DIR)/UnitTests/SyntheticData.cpp
	CPU_FILES += $(SRC_DIR)/UnitTests/TestClaire.cpp
	CPU_FILES += $(SRC_DIR)/UnitTests/TestInterpolation.cpp
	APP_FILES += $(APP_DIR)/benchmark.cpp
	APP_FILES += $(APP_DIR)/test.cpp
endif

ifeq ($(BUILD_PYTHON),yes)
	# C++ files for Python (depend on swig)
	SWIG_FILES += $(SRC_DIR)/Interface/PythonInterface.cpp
	SWIG_FILES += $(SRC_DIR)/Interface/pyclaire_wrap.swig.cpp
	# SWIG files
	#SWIG_FILES += $(SRC_DIR)/Interface/pyclaire.i
endif

CPU_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(CPU_FILES))
SWIG_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SWIG_FILES))
GPU_OBJS = $(patsubst %.cu,$(OBJ_DIR)/%.co,$(GPU_FILES))
BINS = $(patsubst $(APP_DIR)/%.cpp,$(BUILD_DIR)/%,$(APP_FILES))

ifeq ($(BUILD_SHARED),yes)
	BINS += $(LIB_DIR)/libclaire.so
	OBJS = $(LIB_DIR)/libclaire.so
else
	OBJS = $(CPU_OBJS) $(GPU_OBJS)
endif

ifeq ($(BUILD_PYTHON),yes)
	BINS += $(LIB_DIR)/_pyclaire.so
endif
