CPPFILES=$(SRCDIR)/RegOpt.cpp \
		$(SRCDIR)/RegToolsOpt.cpp \
		$(SRCDIR)/BenchmarkOpt.cpp \
		$(SRCDIR)/CLAIREUtils.cpp \
		$(SRCDIR)/CLAIREUtilsKernel.cpp \
		$(SRCDIR)/ghost.cpp \
		$(SRCDIR)/GhostPlan.cpp \
		$(SRCDIR)/Interpolation/interp3.cpp \
		$(SRCDIR)/Interpolation/Interp3_Plan.cpp \
		$(SRCDIR)/Differentiation/Differentiation.cpp \
		$(SRCDIR)/Differentiation/DifferentiationKernel.cpp \
		$(SRCDIR)/Differentiation/DifferentiationFD.cpp \
		$(SRCDIR)/Differentiation/DifferentiationSM.cpp \
		$(SRCDIR)/ScaField.cpp \
		$(SRCDIR)/VecField.cpp \
		$(SRCDIR)/TenField.cpp \
		$(SRCDIR)/ReadWriteReg.cpp \
		$(SRCDIR)/SynProbRegistration.cpp \
		$(SRCDIR)/Solver/TransportProblem.cpp \
		$(SRCDIR)/Solver/TransportEquationSL.cpp \
		$(SRCDIR)/Solver/TransportEquationRK2.cpp \
		$(SRCDIR)/Solver/TransportKernel.cpp \
		$(SRCDIR)/Solver/ContinuityEquation.cpp \
		$(SRCDIR)/DeformationFields/DeformationFields.cpp \
		$(SRCDIR)/DeformationFields/DeformationKernel.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasure.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasureKernel.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasureNCC.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasureSL2.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasureSL2aux.cpp \
		$(SRCDIR)/SemiLagrangian/SemiLagrangian.cpp \
		$(SRCDIR)/Optimizer.cpp \
		$(SRCDIR)/KrylovInterface.cpp \
		$(SRCDIR)/TaoInterface.cpp \
		$(SRCDIR)/CLAIREInterface.cpp \
		$(SRCDIR)/MultiLevelPyramid.cpp \
		$(SRCDIR)/Preconditioner.cpp \
		$(SRCDIR)/PreconditionerKernel.cpp \
		$(SRCDIR)/Regularization/Regularization.cpp \
		$(SRCDIR)/Regularization/RegularizationL2.cpp \
		$(SRCDIR)/Regularization/RegularizationH1.cpp \
		$(SRCDIR)/Regularization/RegularizationH2.cpp \
		$(SRCDIR)/Regularization/RegularizationH1SN.cpp \
		$(SRCDIR)/Regularization/RegularizationH2SN.cpp \
		$(SRCDIR)/Regularization/RegularizationH3.cpp \
		$(SRCDIR)/Regularization/RegularizationH3SN.cpp \
		$(SRCDIR)/Regularization/RegularizationKernel.cpp \
		$(SRCDIR)/OptimizationProblem.cpp \
		$(SRCDIR)/CLAIREBase.cpp \
		$(SRCDIR)/CLAIRE.cpp \
		$(SRCDIR)/CLAIREStokes.cpp \
		$(SRCDIR)/CLAIREDivReg.cpp \
		$(SRCDIR)/Preprocessing.cpp \
		$(SRCDIR)/Spectral/Spectral.cpp \
		$(SRCDIR)/Spectral/SpectralKernel.cpp

EXCPPFILES=

#EXCUFILES=
EXCUFILES=$(EXSRCDIR)/interp3_gpu_new.cu

CUFILES=$(SRCDIR)/Solver/TransportKernel.cu \
		$(SRCDIR)/CLAIREUtilsKernel.cu \
		$(SRCDIR)/DeformationFields/DeformationKernel.cu \
		$(SRCDIR)/DistanceMeasure/DistanceMeasureKernel.cu \
		$(SRCDIR)/Differentiation/DifferentiationKernel.cu \
		$(SRCDIR)/Differentiation/TextureDifferentiationKernel.cu \
		$(SRCDIR)/Spectral/SpectralKernel.cu \
		$(SRCDIR)/PreconditionerKernel.cu \
		$(SRCDIR)/Regularization/RegularizationKernel.cu \
		$(SRCDIR)/Interpolation/Interp3_Plan_GPU.cu \
		$(SRCDIR)/Interpolation/Interp3_Plan_GPU_kernel.cu \
		$(SRCDIR)/SemiLagrangian/SemiLagrangianKernel.cu

CPPFILESCUDA=$(SRCDIR)/RegOpt.cpp \
		$(SRCDIR)/RegToolsOpt.cpp \
		$(SRCDIR)/BenchmarkOpt.cpp \
		$(SRCDIR)/CLAIREUtils.cpp \
		$(SRCDIR)/GhostPlan.cpp \
		$(SRCDIR)/Differentiation/Differentiation.cpp \
		$(SRCDIR)/CLAIREUtilsKernel.cpp \
		$(SRCDIR)/Differentiation/DifferentiationFD.cpp \
		$(SRCDIR)/Differentiation/DifferentiationSM.cpp \
		$(SRCDIR)/VecField.cpp \
		$(SRCDIR)/ScaField.cpp \
		$(SRCDIR)/TenField.cpp \
		$(SRCDIR)/ReadWriteReg.cpp \
		$(SRCDIR)/SynProbRegistration.cpp \
		$(SRCDIR)/Solver/TransportProblem.cpp \
		$(SRCDIR)/Solver/TransportEquationSL.cpp \
		$(SRCDIR)/Solver/TransportEquationRK2.cpp \
		$(SRCDIR)/Solver/ContinuityEquation.cpp \
		$(SRCDIR)/DeformationFields/DeformationFields.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasure.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasureNCC.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasureSL2.cpp \
		$(SRCDIR)/DistanceMeasure/DistanceMeasureSL2aux.cpp \
		$(SRCDIR)/Optimizer.cpp \
		$(SRCDIR)/KrylovInterface.cpp \
		$(SRCDIR)/TaoInterface.cpp \
		$(SRCDIR)/CLAIREInterface.cpp \
		$(SRCDIR)/MultiLevelPyramid.cpp \
		$(SRCDIR)/Preconditioner.cpp \
		$(SRCDIR)/Regularization/Regularization.cpp \
		$(SRCDIR)/Regularization/RegularizationL2.cpp \
		$(SRCDIR)/Regularization/RegularizationH1.cpp \
		$(SRCDIR)/Regularization/RegularizationH2.cpp \
		$(SRCDIR)/Regularization/RegularizationH1SN.cpp \
		$(SRCDIR)/Regularization/RegularizationH2SN.cpp \
		$(SRCDIR)/Regularization/RegularizationH3.cpp \
		$(SRCDIR)/Regularization/RegularizationH3SN.cpp \
		$(SRCDIR)/OptimizationProblem.cpp \
		$(SRCDIR)/CLAIREBase.cpp \
		$(SRCDIR)/CLAIRE.cpp \
		$(SRCDIR)/CLAIREStokes.cpp \
		$(SRCDIR)/CLAIREDivReg.cpp \
		$(SRCDIR)/Preprocessing.cpp \
		$(SRCDIR)/Spectral/Spectral.cpp \
		$(SRCDIR)/SemiLagrangian/SemiLagrangianGPUNew.cpp \
		$(SRCDIR)/Spectral/mpicufft.cpp
		#$(SRCDIR)/ghost.cpp \
		#$(SRCDIR)/SemiLagrangian/SemiLagrangian.cpp \
		#$(SRCDIR)/Interpolation/interp3.cpp \
		#$(SRCDIR)/Interpolation/Interp3_Plan.cpp

TESTFILES=$(SRCDIR)/UnitTestOpt.cpp \
		$(SRCDIR)/UnitTests/TestRegularization.cpp \
		$(SRCDIR)/UnitTests/TestDifferentiation.cpp \
		$(SRCDIR)/UnitTests/SyntheticData.cpp \
		$(SRCDIR)/UnitTests/TestClaire.cpp \
		$(SRCDIR)/UnitTests/TestInterpolation.cpp

SWIGFILE=$(SRCDIR)/Interface/PythonInterface.cpp \
		$(SRCDIR)/Interface/pyclaire_wrap.cxx
