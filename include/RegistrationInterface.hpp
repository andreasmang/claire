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
#include "RegUtils.hpp"
#include "ReadWriteReg.hpp"
#include "Optimizer.hpp"
#include "PrecondReg.hpp"
#include "MultiLevelPyramid.hpp"
#include "OptimalControlRegistrationBase.hpp"
#include "OptimalControlRegistration.hpp"
#include "OptimalControlRegistrationIC.hpp"
#include "OptimalControlRegistrationRelaxedIC.hpp"




namespace reg {




class RegistrationInterface {
 public:
    typedef Optimizer OptimizerType;
    typedef ReadWriteReg ReadWriteType;
    typedef Preprocessing PreProcType;
    typedef OptimalControlRegistrationBase RegProblemType;

    RegistrationInterface();
    ~RegistrationInterface();
    RegistrationInterface(RegOpt*);

    PetscErrorCode Run();
    PetscErrorCode Finalize();

    PetscErrorCode SetReadWrite(ReadWriteReg*);
    PetscErrorCode SetTemplateImage(Vec);
    PetscErrorCode SetAuxVariable(Vec);
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
    PrecondReg* m_Precond;
    ReadWriteType* m_ReadWrite;
    OptimizerType* m_Optimizer;
    RegProblemType* m_RegProblem;

    MultiLevelPyramid *m_TemplatePyramid;
    MultiLevelPyramid *m_ReferencePyramid;

    Vec m_TemplateImage;    ///< original template image (not overwritten)
    Vec m_ReferenceImage;   ///< original reference image (not overwritten)
    VecField* m_Solution;   ///< initial guess

    bool m_IsTemplateSet;   ///< flag: delete the template image (allocated locally)
    bool m_IsReferenceSet;  ///< flag: delete the reference image (allocated locally)
    bool m_DeleteSolution;  ///< flag: delete the solution vector (allocated locally)
};




}  // namespace reg




#endif  // _REGISTRATIONINTERFACE_H_



