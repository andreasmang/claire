/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CLAIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef _REGISTRATIONINTERFACE_H_
#define _REGISTRATIONINTERFACE_H_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "ReadWriteReg.hpp"
#include "Optimizer.hpp"
#include "Preconditioner.hpp"
#include "MultiLevelPyramid.hpp"
#include "CLAIREBase.hpp"
#include "CLAIRE.hpp"
#include "CLAIREStokes.hpp"
#include "CLAIREDivReg.hpp"




namespace reg {




class CLAIREInterface {
 public:
    typedef Optimizer OptimizerType;
    typedef ReadWriteReg ReadWriteType;
    typedef Preprocessing PreProcType;
    typedef CLAIREBase RegProblemType;

    CLAIREInterface();
    virtual ~CLAIREInterface();
    CLAIREInterface(RegOpt*);

    PetscErrorCode Run();
    PetscErrorCode Finalize();

    PetscErrorCode SetReadWrite(ReadWriteReg*);
    PetscErrorCode SetTemplateImage(Vec);
    PetscErrorCode SetAuxVariable(Vec);
    PetscErrorCode SetMask(Vec);
    PetscErrorCode SetCellDensity(Vec);
    PetscErrorCode SetReferenceImage(Vec);
    PetscErrorCode SetSolutionVector(VecField*);
    PetscErrorCode SetInitialGuess(VecField*, bool copy = false);
    PetscErrorCode GetSolution(VecField*, bool copy = false);

//    PetscErrorCode EvaluateDistanceMeasure(ScalarType&, VecField*);  // TODO
    PetscErrorCode EvaluateRegularizationFunctional(ScalarType*, VecField*);
    PetscErrorCode EvaluateGradient(ScalarType*, VecField*);

    PetscErrorCode GetFinalState(Vec);
    PetscErrorCode RunPostProcessing();
    PetscErrorCode ComputeDefFields();
    PetscErrorCode ComputeDetDefGrad(Vec);
    PetscErrorCode ComputeDeformationMap(VecField*);

    PetscErrorCode SolveForwardProblem();
    PetscErrorCode SolveForwardProblem(Vec, Vec);
    PetscErrorCode SolveAdjointProblem(Vec, Vec);

 private:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();
    PetscErrorCode SetupData(Vec&, Vec&);
    PetscErrorCode SetupSolver();
    PetscErrorCode SetupRegProblem();

    PetscErrorCode DispLevelMsg(std::string, int);

    PetscErrorCode RunSolver();
    PetscErrorCode RunSolverGridCont();
    PetscErrorCode RunSolverScaleCont();
    PetscErrorCode RunSolverRegParaCont();
    PetscErrorCode RunSolverRegParaContBinarySearch();
    PetscErrorCode RunSolverRegParaContReductSearch();
    PetscErrorCode RunSolverRegParaContReduction();

    PetscErrorCode ProlongVelocityField(VecField*&, int);

    RegOpt* m_Opt;
    PreProcType* m_PreProc;
    Preconditioner* m_Precond;
    ReadWriteType* m_ReadWrite;
    OptimizerType* m_Optimizer;
    RegProblemType* m_RegProblem;

    MultiLevelPyramid *m_TemplatePyramid;
    MultiLevelPyramid *m_ReferencePyramid;

    Vec m_TemplateImage;    ///< original template image (not overwritten)
    Vec m_ReferenceImage;   ///< original reference image (not overwritten)
    Vec m_Mask;             ///< mask image
    VecField* m_Solution;   ///< initial guess
    Vec m_AuxVariable;      ///< auxilariy variable
    Vec m_CellDensity;      ///< cell density

    bool m_IsTemplateSet;   ///< flag: delete the template image (allocated locally)
    bool m_IsReferenceSet;  ///< flag: delete the reference image (allocated locally)
    bool m_IsMaskSet;       ///< flag: delete the mask image (allocated locally)
    bool m_DeleteSolution;  ///< flag: delete the solution vector (allocated locally)
};




}  // namespace reg




#endif  // _REGISTRATIONINTERFACE_H_
