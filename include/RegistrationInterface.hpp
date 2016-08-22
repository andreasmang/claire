
/**
 * @file RegistrationInterface.hpp
 *
 * Copyright (c) 2015-2016.
 * All rights reserved.
 * This file is part of the XXX library.
 *
 * XXX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * XXX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @brief Generic interface for registration (allows to control the solver)
 *
 * @author Andreas Mang
 *
 */

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



namespace reg
{


class RegistrationInterface
{

public:

    typedef Optimizer OptimizerType;
    typedef ReadWriteReg ReadWriteType;
    typedef PreProcReg PreProcType;
    typedef OptimalControlRegistrationBase RegProblemType;

    RegistrationInterface();
    ~RegistrationInterface();
    RegistrationInterface(RegOpt*);

    PetscErrorCode Run();
    PetscErrorCode Finalize();

    PetscErrorCode SetReadWrite(ReadWriteReg*);
    PetscErrorCode SetTemplateImage(Vec);
    PetscErrorCode SetReferenceImage(Vec);
    PetscErrorCode SetInitialGuess(VecField*);

    PetscErrorCode RunPostProcessing();
    PetscErrorCode ComputeDefFields();

private:

    PetscErrorCode Initialize(void);
    PetscErrorCode ClearMemory(void);
    PetscErrorCode SetupSolver(void);
    PetscErrorCode SetupRegProblem(void);

    PetscErrorCode DispLevelMsg(std::string,int);

    PetscErrorCode RunSolver(void);
    PetscErrorCode RunSolverGridCont(void);
    PetscErrorCode RunSolverScaleCont(void);
    PetscErrorCode RunSolverRegParaCont(void);
    PetscErrorCode RunSolverRegParaContBinarySearch(void);
    PetscErrorCode RunSolverRegParaContReductSearch(void);
    PetscErrorCode RunSolverRegParaContReduction(void);

    PetscErrorCode ProlongVelocityField(VecField*&,int);

    RegOpt* m_Opt;
    PreProcType* m_PreProc;
    PrecondReg* m_Precond;
    ReadWriteType* m_ReadWrite;
    OptimizerType* m_Optimizer;
    RegProblemType* m_RegProblem;

    MultiLevelPyramid *m_TemplatePyramid;
    MultiLevelPyramid *m_ReferencePyramid;

    Vec m_TemplateImage; ///< original template image (not overwritten)
    Vec m_ReferenceImage; ///< original reference image (not overwritten)
    VecField* m_Solution; ///< initial guess
};



} // end of name space

#endif // _REGISTRATIONINTERFACE_H_



