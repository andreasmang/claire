/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the XXX library.
 *
 *  XXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  XXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XXX. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


#include <string>

#include "RegToolsOpt.hpp"
#include "RegUtils.hpp"
#include "PreProcReg.hpp"
#include "VecField.hpp"
#include "ReadWriteReg.hpp"
#include "SynProbRegistration.hpp"
#include "RegistrationInterface.hpp"

PetscErrorCode RunPostProcessing(reg::RegToolsOpt*);
PetscErrorCode ResampleVecField(reg::RegToolsOpt*);
PetscErrorCode ResampleScaField(reg::RegToolsOpt*);
PetscErrorCode ComputeDefFields(reg::RegToolsOpt*);
PetscErrorCode ComputeGrad(reg::RegToolsOpt*);
PetscErrorCode ComputeResidual(reg::RegToolsOpt*);
PetscErrorCode ComputeSynVel(reg::RegToolsOpt*);
PetscErrorCode SolveForwardProblem(reg::RegToolsOpt*);
PetscErrorCode CheckForwardSolve(reg::RegToolsOpt*);


/********************************************************************
 * @brief main function for registration tools
 *******************************************************************/
int main(int argc, char **argv) {
    PetscErrorCode ierr = 0;

    reg::RegToolsOpt* regopt = NULL;

    // initialize petsc (user is not allowed to set petsc options)
    ierr = PetscInitialize(0, reinterpret_cast<char***>(NULL),
                              reinterpret_cast<char*>(NULL),
                              reinterpret_cast<char*>(NULL)); CHKERRQ(ierr);
    PetscFunctionBegin;

    // allocate class for controlling everything
    try {regopt = new reg::RegToolsOpt(argc, argv);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    if (regopt->GetFlags().computedeffields) {
        ierr = ComputeDefFields(regopt); CHKERRQ(ierr);
    } else if (regopt->GetFlags().computegrad) {
        ierr = ComputeGrad(regopt); CHKERRQ(ierr);
    } else if (regopt->GetFlags().resample) {
        if (regopt->GetFlags().readvecfield) {
            ierr = ResampleVecField(regopt); CHKERRQ(ierr);
        }
        if (regopt->GetFlags().readscafield) {
            ierr = ResampleScaField(regopt); CHKERRQ(ierr);
        }
    } else if (regopt->GetFlags().tscafield) {
        ierr = SolveForwardProblem(regopt); CHKERRQ(ierr);
    } else if (regopt->GetFlags().tlabelmap) {
        ierr = SolveForwardProblem(regopt); CHKERRQ(ierr);
    } else if (regopt->GetFlags().computeresidual) {
        ierr = ComputeResidual(regopt); CHKERRQ(ierr);
    } else if (regopt->GetFlags().computesynvel) {
        ierr = ComputeSynVel(regopt); CHKERRQ(ierr);
    } else if (regopt->GetFlags().checkfwdsolve) {
        ierr = CheckForwardSolve(regopt); CHKERRQ(ierr);
    }

    // clean up
    if (regopt != NULL) {delete regopt; regopt = NULL;}

    // clean up petsc
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}




/********************************************************************
 * @brief compute gradient of scalar field
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeGrad"
PetscErrorCode ComputeGrad(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    std::string filename, fnx1, fnx2, fnx3;
    std::stringstream ss;
    int rank;
    double timers[5] = {0, 0, 0, 0, 0};
    ScalarType *p_m = NULL, *p_gm1 = NULL, *p_gm2 = NULL, *p_gm3 = NULL;
    Vec m = NULL;
    std::bitset<3> XYZ; XYZ[0] = 1; XYZ[1] = 1; XYZ[2] = 1;
    reg::VecField* grad = NULL;
    reg::ReadWriteReg* readwrite = NULL;

    PetscFunctionBegin;

    regopt->Enter(__FUNCT__);

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (regopt->GetFlags().readscafield) {
        // read velocity components
        filename = regopt->GetScaFieldFN(0);
        ierr = readwrite->Read(&m, filename); CHKERRQ(ierr);
        ierr = reg::Assert(m != NULL, "null pointer"); CHKERRQ(ierr);

        if (!regopt->SetupDone()) {ierr = regopt->DoSetup(); CHKERRQ(ierr);}

        try {grad = new reg::VecField(regopt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        // computing gradient of m
        ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
        ierr = grad->GetArrays(p_gm1, p_gm2, p_gm3); CHKERRQ(ierr);
        accfft_grad(p_gm1, p_gm2, p_gm3, p_m, regopt->GetFFT().plan, &XYZ, timers);
        ierr = grad->RestoreArrays(p_gm1, p_gm2, p_gm3); CHKERRQ(ierr);
        ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);

        // write to file
        fnx1 = regopt->GetVecFieldFN(0, 1);
        fnx2 = regopt->GetVecFieldFN(1, 1);
        fnx3 = regopt->GetVecFieldFN(2, 1);
        ierr = readwrite->Write(grad, fnx1, fnx2, fnx3); CHKERRQ(ierr);

    } else if (regopt->GetFlags().readvecfield) {
    }

    if (m != NULL) {ierr = VecDestroy(&m); CHKERRQ(ierr);}
    if (grad != NULL) {delete grad; grad = NULL;}
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}

    regopt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief post process image registration results
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunPostProcessing"
PetscErrorCode RunPostProcessing(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    std::string ifolder, xfolder, filename, ext;
    Vec mT = NULL, mR = NULL, vx1 = NULL, vx2 = NULL, vx3 = NULL;
    reg::VecField *v = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    reg::RegistrationInterface* registration = NULL;

    PetscFunctionBegin;

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    ierr = reg::Msg("processing results"); CHKERRQ(ierr);

    ext = regopt->GetReadWriteFlags().extension;
    ifolder = regopt->GetReadWriteFlags().ifolder;
    ierr = reg::Assert(ifolder.empty() != true, "input folder needs to be provided"); CHKERRQ(ierr);

    // read template image
    filename = ifolder + "template-image" + ext;
    ierr = readwrite->Read(&mT, filename); CHKERRQ(ierr);
    ierr = reg::Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ) { ierr = regopt->DoSetup(); CHKERRQ(ierr); }

    // read reference image
    filename = ifolder + "reference-image" + ext;
    ierr = readwrite->Read(&mR, filename); CHKERRQ(ierr);
    ierr = reg::Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);

    if (!regopt->SetupDone()) {ierr = regopt->DoSetup(); CHKERRQ(ierr);}

    // allocate container for velocity field
    try { v = new reg::VecField(regopt); }
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // read velocity components
    filename = ifolder + "velocity-field-x1" + ext;
    ierr = readwrite->Read(&vx1, filename); CHKERRQ(ierr);
    ierr = VecCopy(vx1, v->m_X1); CHKERRQ(ierr);

    filename = ifolder + "velocity-field-x2" + ext;
    ierr = readwrite->Read(&vx2, filename); CHKERRQ(ierr);
    ierr = VecCopy(vx2, v->m_X2); CHKERRQ(ierr);

    filename = ifolder + "velocity-field-x3" + ext;
    ierr = readwrite->Read(&vx3, filename); CHKERRQ(ierr);
    ierr = VecCopy(vx3, v->m_X3); CHKERRQ(ierr);

    // allocate class for registration interface
    try {registration = new reg::RegistrationInterface(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // set all we need
    ierr = registration->SetReadWrite(readwrite); CHKERRQ(ierr);
    ierr = registration->SetReferenceImage(mR); CHKERRQ(ierr);
    ierr = registration->SetTemplateImage(mT); CHKERRQ(ierr);
    ierr = registration->SetInitialGuess(v); CHKERRQ(ierr);

    // run post processing
    ierr = registration->RunPostProcessing(); CHKERRQ(ierr);

    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr); mT = NULL;}
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr); mR = NULL;}

    if (vx1 != NULL) {ierr = VecDestroy(&vx1); CHKERRQ(ierr); vx1 = NULL;}
    if (vx2 != NULL) {ierr = VecDestroy(&vx2); CHKERRQ(ierr); vx2 = NULL;}
    if (vx3 != NULL) {ierr = VecDestroy(&vx3); CHKERRQ(ierr); vx3 = NULL;}

    if (v != NULL) {delete v; v = NULL;}

    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (registration != NULL) {delete registration; registration = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief post process image registration results
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDefFields"
PetscErrorCode ComputeDefFields(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    std::string ifolder, xfolder, filename, ext;
    Vec vxi = NULL;
    reg::VecField* v = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    reg::RegistrationInterface* registration = NULL;

    PetscFunctionBegin;

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ext = regopt->GetReadWriteFlags().extension;
    ifolder = regopt->GetReadWriteFlags().ifolder;
    ierr = reg::Assert(ifolder.empty() != true, "input folder needs to be provided"); CHKERRQ(ierr);

    // read velocity components
    filename = ifolder + "velocity-field-x1" + ext;
    ierr = readwrite->Read(&vxi, filename); CHKERRQ(ierr);
    if (!regopt->SetupDone()) {ierr = regopt->DoSetup(); CHKERRQ(ierr);}

    // allocate container for velocity field
    try {v = new reg::VecField(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    ierr = VecCopy(vxi, v->m_X1); CHKERRQ(ierr);
    if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

    filename = ifolder + "velocity-field-x2" + ext;
    ierr = readwrite->Read(&vxi, filename); CHKERRQ(ierr);
    ierr = VecCopy(vxi, v->m_X2); CHKERRQ(ierr);
    if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

    filename = ifolder + "velocity-field-x3" + ext;
    ierr = readwrite->Read(&vxi, filename); CHKERRQ(ierr);
    ierr = VecCopy(vxi, v->m_X3); CHKERRQ(ierr);
    if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

    // allocate class for registration interface
    try {registration = new reg::RegistrationInterface(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // set all we need
    ierr = registration->SetReadWrite(readwrite); CHKERRQ(ierr);
    ierr = registration->SetInitialGuess(v); CHKERRQ(ierr);

    // run post processing
    ierr = registration->ComputeDefFields(); CHKERRQ(ierr);

    if (v != NULL) {delete v; v = NULL;}
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (registration != NULL) {delete registration; registration = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resample scalar field
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResampleScaField"
PetscErrorCode ResampleScaField(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    std::string filename;
    std::stringstream ss;
    int rank;
    IntType nl, ng, nx[3], nxl[3];
    ScalarType gridscale, value;
    Vec m = NULL, ml = NULL;
    reg::PreProcReg* preproc = NULL;
    reg::ReadWriteReg* readwrite = NULL;

    PetscFunctionBegin;

    regopt->Enter(__FUNCT__);

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // read velocity components
    filename = regopt->GetScaFieldFN(0);
    ierr = readwrite->Read(&m, filename); CHKERRQ(ierr);
    ierr = reg::Assert(m != NULL, "null pointer"); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ) {ierr = regopt->DoSetup(); CHKERRQ(ierr);}

    // compute grid size
    gridscale = regopt->GetResamplingPara().gridscale;

    if (gridscale != 1.0) {
        // allocate container for velocity field
        try {preproc = new reg::PreProcReg(regopt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        preproc->ResetGridChangeOps(true);

        for (int i=0; i < 3; ++i) {
            nx[i] = regopt->GetDomainPara().nx[i];
            value = gridscale*static_cast<ScalarType>(regopt->GetDomainPara().nx[i]);
            nxl[i] = static_cast<IntType>(ceil(value));
        }
        ierr = regopt->GetSizes(nxl, nl, ng); CHKERRQ(ierr);

        ss << "resampling scalar field  (" << nx[0] << "," << nx[1] << "," << nx[2] << ")"
           << " -> (" << nxl[0] << "," << nxl[1] << "," << nxl[2] << ")";
        ierr = reg::Msg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();

        // allocate array
        ierr = reg::VecCreate(ml, nl, ng); CHKERRQ(ierr);

        // restrict of prolong the vector field
        if (gridscale > 1.0) {
            ierr = preproc->Prolong(&ml, m, nxl, nx); CHKERRQ(ierr);
        } else {
            ierr = preproc->Restrict(&ml, m, nxl, nx); CHKERRQ(ierr);
        }

        // reset io
        if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
        for (int i=0; i < 3; ++i) {
            regopt->SetNumGridPoints(i, nxl[i]);
        }
        ierr = regopt->DoSetup(false); CHKERRQ(ierr);

        try {readwrite = new reg::ReadWriteReg(regopt);}
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        // write resampled scalar field to file
        filename = regopt->GetScaFieldFN(1);
        ierr = readwrite->Write(ml, filename); CHKERRQ(ierr);

    } else {
        // simply write field to file
        filename = regopt->GetScaFieldFN(1);
        ierr = readwrite->Write(m, filename); CHKERRQ(ierr);
    }

    if (m != NULL) {ierr = VecDestroy(&m); CHKERRQ(ierr);}
    if (ml != NULL) {ierr = VecDestroy(&ml); CHKERRQ(ierr);}

    if (preproc != NULL) {delete preproc; preproc = NULL;}
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}

    regopt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief resample vector field
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ResampleVecField"
PetscErrorCode ResampleVecField(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    std::string filename, fnx1, fnx2, fnx3;
    std::stringstream ss;
    IntType nl, ng, nx[3], nxl[3];
    ScalarType scale;
    Vec vx1 = NULL, vx2 = NULL, vx3 = NULL;
    reg::VecField *v = NULL, *vl = NULL;
    reg::PreProcReg* preproc = NULL;
    reg::ReadWriteReg* readwrite = NULL;

    PetscFunctionBegin;

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // read velocity components
    filename = regopt->GetVecFieldFN(0, 0);
    ierr = readwrite->Read(&vx1, filename); CHKERRQ(ierr);

    filename = regopt->GetVecFieldFN(1, 0);
    ierr = readwrite->Read(&vx2, filename); CHKERRQ(ierr);

    filename = regopt->GetVecFieldFN(2, 0);
    ierr = readwrite->Read(&vx3, filename); CHKERRQ(ierr);

    if (!regopt->SetupDone()) {ierr = regopt->DoSetup(); CHKERRQ(ierr);}

    // allocate container for velocity field
    try {v = new reg::VecField(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    ierr = VecCopy(vx1, v->m_X1); CHKERRQ(ierr);
    ierr = VecCopy(vx2, v->m_X2); CHKERRQ(ierr);
    ierr = VecCopy(vx3, v->m_X3); CHKERRQ(ierr);

    // allocate container for velocity field
    try {preproc = new reg::PreProcReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    preproc->ResetGridChangeOps(true);

    // compute grid size
    scale = regopt->GetResamplingPara().gridscale;
    for (int i=0; i < 3; ++i) {
        nx[i]  = regopt->GetDomainPara().nx[i];
        nxl[i] = scale*regopt->GetDomainPara().nx[i];
    }
    ierr = regopt->GetSizes(nxl, nl, ng); CHKERRQ(ierr);

    ss << "resampling vector field  (" << nx[0] << "," << nx[1] << "," << nx[2] << ")"
       << " -> (" << nxl[0] << "," << nxl[1] << "," << nxl[2] << ")";
    ierr = reg::Msg(ss.str()); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();

    // allocate container for velocity field
    try {vl = new reg::VecField(nl, ng);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // restrict of prolong the vector field
    if (scale > 1.0) {
        ierr = preproc->Prolong(vl, v, nxl, nx); CHKERRQ(ierr);
    } else {
        ierr = preproc->Restrict(vl, v, nxl, nx); CHKERRQ(ierr);
    }

    // reset io
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    for (int i=0; i < 3; ++i) {
        regopt->SetNumGridPoints(i, nxl[i]);
    }
    ierr = regopt->DoSetup(false); CHKERRQ(ierr);
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // get output file name (based on input file name)
    fnx1 = regopt->GetVecFieldFN(0, 1);
    fnx2 = regopt->GetVecFieldFN(1, 1);
    fnx3 = regopt->GetVecFieldFN(2, 1);

    // write to file
    ierr = readwrite->Write(vl, fnx1, fnx2, fnx3); CHKERRQ(ierr);

    if (v != NULL) {delete v; v = NULL; }
    if (vl != NULL) {delete vl; vl = NULL;}

    if (vx1 != NULL) {ierr = VecDestroy(&vx1); CHKERRQ(ierr); vx1 = NULL;}
    if (vx2 != NULL) {ierr = VecDestroy(&vx2); CHKERRQ(ierr); vx2 = NULL;}
    if (vx3 != NULL) {ierr = VecDestroy(&vx3); CHKERRQ(ierr); vx3 = NULL;}

    if (preproc != NULL) {delete preproc; preproc = NULL;}
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief solve forward problem
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SolveForwardProblem"
PetscErrorCode SolveForwardProblem(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    IntType nl;
    std::string fn;
    std::stringstream ss;
    reg::VecField* v = NULL;
    Vec m0 = NULL, m1 = NULL, vxi = NULL;
    ScalarType *p_m1 = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    reg::RegistrationInterface* registration = NULL;
    PetscFunctionBegin;

    regopt->Enter(__FUNCT__);

    // allocate class for io
    try { readwrite = new reg::ReadWriteReg(regopt); }
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // read velocity components
    fn = regopt->GetScaFieldFN(0);
    ierr = readwrite->Read(&m0, fn); CHKERRQ(ierr);
    ierr = reg::Assert(m0 != NULL, "null pointer"); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ) {ierr = regopt->DoSetup(); CHKERRQ(ierr);}
    nl = regopt->GetDomainPara().nlocal;

    // do allocation
    ierr = VecDuplicate(m0, &m1); CHKERRQ(ierr);
    try {v = new reg::VecField(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // read the individual components of the vector fields
    fn = regopt->GetVecFieldFN(0, 0);
    ierr = readwrite->Read(&vxi, fn); CHKERRQ(ierr);
    ierr = VecCopy(vxi, v->m_X1); CHKERRQ(ierr);
    if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

    fn = regopt->GetVecFieldFN(1, 0);
    ierr = readwrite->Read(&vxi, fn); CHKERRQ(ierr);
    ierr = VecCopy(vxi, v->m_X2); CHKERRQ(ierr);
    if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

    fn = regopt->GetVecFieldFN(2, 0);
    ierr = readwrite->Read(&vxi, fn); CHKERRQ(ierr);
    ierr = VecCopy(vxi, v->m_X3); CHKERRQ(ierr);
    if (vxi != NULL) {ierr = VecDestroy(&vxi); CHKERRQ(ierr); vxi = NULL;}

    // allocate class for registration interface
    try {registration = new reg::RegistrationInterface(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // set all we need
    ierr = registration->SetReadWrite(readwrite); CHKERRQ(ierr);
    ierr = registration->SetInitialGuess(v); CHKERRQ(ierr);

    // run post processing
    ierr = registration->SolveForwardProblem(m1, m0); CHKERRQ(ierr);

    // if we consider a lable map, we want to truncate the values
    if (regopt->GetFlags().tlabelmap) {
        // map to integer
        ierr = VecGetArray(m1, &p_m1); CHKERRQ(ierr);
        for (IntType i = 0; i < nl; ++i) {
            p_m1[i] = round(p_m1[i]);
        }
        ierr = VecGetArray(m1, &p_m1); CHKERRQ(ierr);
    }

    // write resampled scalar field to file
    fn = regopt->GetScaFieldFN(1);
    ierr = readwrite->Write(m1, fn); CHKERRQ(ierr);

    if (v != NULL) {delete v; v = NULL;}
    if (m0 != NULL) {ierr = VecDestroy(&m0); CHKERRQ(ierr); m0 = NULL;}
    if (m1 != NULL) {ierr = VecDestroy(&m1); CHKERRQ(ierr); m1 = NULL;}
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (registration != NULL) {delete registration; registration = NULL;}

    regopt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute residual
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeResidual"
PetscErrorCode ComputeResidual(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    IntType nl;
    std::string fn;
    std::stringstream ss;
    int rank;
    Vec mR = NULL, mT = NULL;
    ScalarType *p_mr = NULL, *p_mt = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    PetscFunctionBegin;

    regopt->Enter(__FUNCT__);

    // allocate class for io
    try { readwrite = new reg::ReadWriteReg(regopt); }
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // read reference image
    fn = regopt->GetScaFieldFN(2);
    ierr = readwrite->Read(&mR, fn); CHKERRQ(ierr);
    ierr = reg::Assert(mR != NULL, "null pointer"); CHKERRQ(ierr);

    if ( !regopt->SetupDone() ) { ierr = regopt->DoSetup(); CHKERRQ(ierr); }

    // read template image
    fn = regopt->GetScaFieldFN(3);
    ierr = readwrite->Read(&mT, fn); CHKERRQ(ierr);
    ierr = reg::Assert(mT != NULL, "null pointer"); CHKERRQ(ierr);

    ierr = reg::Rescale(mR, 0, 1); CHKERRQ(ierr);
    ierr = reg::Rescale(mT, 0, 1); CHKERRQ(ierr);

    nl = regopt->GetDomainPara().nlocal;

    ierr = VecGetArray(mR, &p_mr); CHKERRQ(ierr);
    ierr = VecGetArray(mT, &p_mt); CHKERRQ(ierr);
#pragma omp parallel
{
#pragma omp for
    for (IntType i = 0; i < nl; ++i) {  // for all grid points
        p_mt[i] = 1.0 - PetscAbs(p_mt[i] - p_mr[i]);
    }
}  // pragma omp
    ierr = VecRestoreArray(mR, &p_mr); CHKERRQ(ierr);
    ierr = VecRestoreArray(mT, &p_mt); CHKERRQ(ierr);

    // write resampled scalar field to file
    fn = regopt->GetScaFieldFN(1);
    ierr = readwrite->Write(mT, fn); CHKERRQ(ierr);

    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (mT != NULL) {ierr = VecDestroy(&mT); CHKERRQ(ierr); mT = NULL;}
    if (mR != NULL) {ierr = VecDestroy(&mR); CHKERRQ(ierr); mR = NULL;}

    regopt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute synthetic velocity field
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeSynVel"
PetscErrorCode ComputeSynVel(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    int problem = 3;
    reg::VecField* v = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    ScalarType *p_vx1 = NULL, *p_vx2 = NULL, *p_vx3 = NULL;
    ScalarType hx[3];
    std::string fnx1, fnx2, fnx3;
    PetscFunctionBegin;

    regopt->Enter(__FUNCT__);

    ierr = regopt->DoSetup(); CHKERRQ(ierr);

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // allocate class for io
    try {v = new reg::VecField(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    for (int i = 0; i < 3; ++i) {
        hx[i] = regopt->GetDomainPara().hx[i];
    }

    ierr = v->GetArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);
#pragma omp parallel
{
    ScalarType x1, x2, x3;
    IntType i1, i2, i3, i;
#pragma omp for
    for (i1 = 0; i1 < regopt->GetDomainPara().isize[0]; ++i1) {  // x1
        for (i2 = 0; i2 < regopt->GetDomainPara().isize[1]; ++i2) {  // x2
            for (i3 = 0; i3 < regopt->GetDomainPara().isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + regopt->GetDomainPara().istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + regopt->GetDomainPara().istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + regopt->GetDomainPara().istart[2]);

                // compute linear / flat index
                i = reg::GetLinearIndex(i1, i2, i3, regopt->GetDomainPara().isize);
                // compute the velocity field
                if (problem == 0) {
                    ScalarType v0 = 0.5;
                    p_vx1[i] = v0*sin(x3)*cos(x2)*sin(x2);
                    p_vx2[i] = v0*sin(x1)*cos(x3)*sin(x3);
                    p_vx3[i] = v0*sin(x2)*cos(x1)*sin(x1);
                } else if (problem == 1) {
                    p_vx1[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                    p_vx2[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                    p_vx3[i] = sin(2.0*x1)*cos(2.0*x2)*sin(2.0*x3);
                } else if (problem == 3) {
                    p_vx1[i] = cos(x1)*sin(x2);
                    p_vx2[i] = cos(x2)*sin(x1);
                    p_vx3[i] = cos(x1)*sin(x3);
                } else if (problem == 4) {
                    p_vx1[i] = cos(x2)*cos(x3);
                    p_vx2[i] = sin(x3)*sin(x1);
                    p_vx3[i] = cos(x1)*cos(x2);
                }
            }  // i1
        }  // i2
    }  // i3
}  // pragma omp parallel
    ierr = v->RestoreArrays(p_vx1, p_vx2, p_vx3); CHKERRQ(ierr);

    // write computed vectorfield to file
    fnx1 = "velocity-field-x1" + regopt->GetReadWriteFlags().extension;
    fnx2 = "velocity-field-x2" + regopt->GetReadWriteFlags().extension;
    fnx3 = "velocity-field-x3" + regopt->GetReadWriteFlags().extension;
    ierr = readwrite->Write(v, fnx1, fnx2, fnx3); CHKERRQ(ierr);

    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (v != NULL) {delete v; v = NULL;}

    regopt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief check the forward solver
 * @param[in] regopt container for user defined options
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "CheckForwardSolve"
PetscErrorCode CheckForwardSolve(reg::RegToolsOpt* regopt) {
    PetscErrorCode ierr = 0;
    IntType nc, nl, ng;
    Vec m0 = NULL, m1 = NULL;
    reg::VecField *v = NULL;
    reg::RegistrationInterface* registration = NULL;
    reg::SynProbRegistration* synprob = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    PetscFunctionBegin;

    regopt->Enter(__FUNCT__);

    nc = 2;

    ierr = regopt->DoSetup(); CHKERRQ(ierr);

    regopt->SetNumImageComponents(nc);

    nc = regopt->GetDomainPara().nc;
    nl = regopt->GetDomainPara().nlocal;
    ng = regopt->GetDomainPara().nglobal;

    // allocation
    try {v = new reg::VecField(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {readwrite = new reg::ReadWriteReg(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {registration = new reg::RegistrationInterface(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    try {synprob = new reg::SynProbRegistration(regopt);}
    catch (std::bad_alloc&) {
        ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    ierr = reg::VecCreate(m0, nc*nl, nc*ng); CHKERRQ(ierr);
    ierr = reg::VecCreate(m1, nc*nl, nc*ng); CHKERRQ(ierr);
    ierr = synprob->ComputeSmoothScalarField(m0, 0); CHKERRQ(ierr);
    ierr = registration->SetInitialGuess(v); CHKERRQ(ierr);
    ierr = registration->SolveForwardProblem(m1, m0); CHKERRQ(ierr);

    ierr = readwrite->WriteMC(m0, "initial-condition.nc"); CHKERRQ(ierr);
    ierr = readwrite->WriteMC(m1, "final-condition.nc"); CHKERRQ(ierr);

    regopt->Exit(__FUNCT__);

    if (v != NULL) {delete v; v = NULL;}
    if (m0 != NULL) {ierr = VecDestroy(&m0); CHKERRQ(ierr); m0 = NULL;}
    if (m1 != NULL) {ierr = VecDestroy(&m1); CHKERRQ(ierr); m1 = NULL;}
    if (synprob != NULL) {delete synprob; synprob = NULL;}
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (registration != NULL) {delete registration; registration = NULL;}

    PetscFunctionReturn(ierr);
}



