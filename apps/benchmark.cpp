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

#include "CLAIREUtils.hpp"
#include "BenchmarkOpt.hpp"
#include "VecField.hpp"
#include "CLAIRE.hpp"
#include <sys/stat.h>

PetscErrorCode RunForwardSolverBenchmark(reg::BenchmarkOpt*);
PetscErrorCode RunGradientBenchmark(reg::BenchmarkOpt*);
PetscErrorCode RunHessianMatvecBenchmark(reg::BenchmarkOpt*);

PetscErrorCode ComputeErrorForwardSolver(reg::BenchmarkOpt*);

PetscErrorCode ComputeSyntheticData(Vec&, reg::BenchmarkOpt*);
PetscErrorCode ComputeSyntheticData(reg::VecField*&, reg::BenchmarkOpt*);




/********************************************************************
 * @brief main function to run registration
 *******************************************************************/
int main(int argc, char **argv) {
    PetscErrorCode ierr = 0;
    reg::BenchmarkOpt* opt = NULL;
    double runtime, value;
    int rval;
    std::stringstream ss;

    // initialize petsc (user is not allowed to set petsc options)
    ierr = PetscInitialize(0, reinterpret_cast<char***>(NULL),
                              reinterpret_cast<char*>(NULL),
                              reinterpret_cast<char*>(NULL)); CHKERRQ(ierr);
    PetscFunctionBegin;

    // allocate class for controlling everything
    try {opt = new reg::BenchmarkOpt(argc, argv);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }
    ierr = opt->DoSetup(); CHKERRQ(ierr);

    IntType n = opt->NumRepeats();
    ss << "number of repetitions " << n;
    ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();

    switch (opt->Benchmark()) {
        case 0:
            ierr = RunForwardSolverBenchmark(opt); CHKERRQ(ierr);
            break;
        case 1:
            ierr = RunGradientBenchmark(opt); CHKERRQ(ierr);
            break;
        case 2:
            ierr = RunHessianMatvecBenchmark(opt); CHKERRQ(ierr);
            break;
        case 3:
            ierr = ComputeErrorForwardSolver(opt); CHKERRQ(ierr);
            break;
        default:
            ierr = reg::ThrowError("benchmark not defined"); CHKERRQ(ierr);
            break;
    }

    // wrap up timers
    ierr = opt->ProcessTimers(); CHKERRQ(ierr);

    // get local runtime
    value = opt->GetRunTime();
    rval = MPI_Reduce(&value, &runtime, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
    ierr = reg::Assert(rval == MPI_SUCCESS, "mpi reduce returned error"); CHKERRQ(ierr);

    // write logfile and display time to solution
    ss << "total runtime (in seconds)   " << std::scientific << runtime;
    ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();
    ss << "average runtime (in seconds) " << std::scientific << runtime/static_cast<double>(n);
    ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
    ss.str(std::string()); ss.clear();


    ierr = opt->WriteLogFile(); CHKERRQ(ierr);
    ierr = opt->DisplayTimeToSolution(); CHKERRQ(ierr);

    if (opt != NULL) {delete opt; opt = NULL;}

    ierr = reg::Finalize(); CHKERRQ(ierr);

    return 0;
}




/********************************************************************
 * @brief perform benchmark for forward solver
 *******************************************************************/
PetscErrorCode RunForwardSolverBenchmark(reg::BenchmarkOpt *opt) {
    PetscErrorCode ierr = 0;
    Vec m = 0;
    reg::VecField* v = NULL;
    reg::CLAIRE* registration = NULL;
    std::stringstream ss;
    double runtime;
    PetscFunctionBegin;

    // make sure we do not store time history
    opt->m_RegFlags.runinversion = false;


    ierr = ComputeSyntheticData(m, opt); CHKERRQ(ierr);
    ierr = ComputeSyntheticData(v, opt); CHKERRQ(ierr);

    try {registration = new reg::CLAIRE(opt);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }
    DBGCHK();

    ierr = registration->SetControlVariable(v); CHKERRQ(ierr);

    DBGCHK();


    ierr = reg::DbgMsg("run forward solver benchmarck"); CHKERRQ(ierr);

    // warm start
    ierr = registration->SolveForwardProblem(NULL, m); CHKERRQ(ierr);

    // reset all timers and counters
    ierr = opt->ResetTimers(); CHKERRQ(ierr);
    ierr = opt->ResetCounters(); CHKERRQ(ierr);

    IntType n = opt->NumRepeats();
    ierr = opt->StartTimer(reg::T2SEXEC); CHKERRQ(ierr);
    runtime = -MPI_Wtime(); // start timer
    for (IntType i = 0; i < n; ++i) {
        if (opt->m_Verbosity > 1) {
            ss  << "forward solve "<< std::setw(3)
                << i << " of " << std::setw(3) << n;
            ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        ierr = registration->SolveForwardProblem(NULL, m); CHKERRQ(ierr);
    }
    runtime += MPI_Wtime(); // stop timer
    ierr = opt->StopTimer(reg::T2SEXEC); CHKERRQ(ierr);
    opt->SetRunTime(runtime);

    if (registration != NULL) {delete registration; registration = NULL;}
    if (m != NULL) {ierr = VecDestroy(&m); CHKERRQ(ierr); m = NULL;}
    if (v != NULL) {delete v; v = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief perform benchmark for forward solver
 *******************************************************************/
PetscErrorCode RunGradientBenchmark(reg::BenchmarkOpt *opt) {
    PetscErrorCode ierr = 0;
    Vec m = 0;
    reg::VecField* v = NULL;
    reg::CLAIRE* registration = NULL;
    std::stringstream ss;
    double runtime;
    PetscFunctionBegin;

    ierr = ComputeSyntheticData(m, opt); CHKERRQ(ierr);
    ierr = ComputeSyntheticData(v, opt); CHKERRQ(ierr);

    try {registration = new reg::CLAIRE(opt);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }
    ierr = registration->SetControlVariable(v); CHKERRQ(ierr);
    ierr = registration->InitializeSolver(); CHKERRQ(ierr);

    ierr = reg::DbgMsg("run gradient benchmarck"); CHKERRQ(ierr);

    // warm start
    ierr = registration->SetReferenceImage(m); CHKERRQ(ierr);
    ierr = registration->SolveForwardProblem(NULL, m); CHKERRQ(ierr);
    ierr = registration->EvaluateGradient(NULL, NULL); CHKERRQ(ierr);

    // reset all timers and counters
    ierr = opt->ResetTimers(); CHKERRQ(ierr);
    ierr = opt->ResetCounters(); CHKERRQ(ierr);

    IntType n = opt->NumRepeats();
    ierr = opt->StartTimer(reg::T2SEXEC); CHKERRQ(ierr);
    runtime = -MPI_Wtime(); // start timer
    for (IntType i = 0; i < n; ++i) {
        if (opt->m_Verbosity > 1) {
            ss  << "gradient evaluation "<< std::setw(3)
                << i << " of " << std::setw(3) << n;
            ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        ierr = registration->EvaluateGradient(NULL, NULL); CHKERRQ(ierr);
    }
    runtime += MPI_Wtime(); // stop timer
    ierr = opt->StopTimer(reg::T2SEXEC); CHKERRQ(ierr);

    opt->SetRunTime(runtime);

    if (registration != NULL) {delete registration; registration = NULL;}
    if (m != NULL) {ierr = VecDestroy(&m); CHKERRQ(ierr); m = NULL;}
    if (v != NULL) {delete v; v = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief perform benchmark for hessian matvec
 *******************************************************************/
PetscErrorCode RunHessianMatvecBenchmark(reg::BenchmarkOpt *opt) {
    PetscErrorCode ierr = 0;
    Vec m = 0;
    reg::VecField *v = NULL, *vtilde = NULL;
    reg::CLAIRE* registration = NULL;
    std::stringstream ss;
    double runtime;
    PetscFunctionBegin;

    ierr = ComputeSyntheticData(m, opt); CHKERRQ(ierr);
    ierr = ComputeSyntheticData(v, opt); CHKERRQ(ierr);
    ierr = ComputeSyntheticData(vtilde, opt); CHKERRQ(ierr);

    try {registration = new reg::CLAIRE(opt);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }
    ierr = registration->SetControlVariable(v); CHKERRQ(ierr);
    ierr = registration->SetIncControlVariable(v); CHKERRQ(ierr);
    ierr = registration->InitializeSolver(); CHKERRQ(ierr);

    ierr = reg::DbgMsg("run hessian matvec benchmarck"); CHKERRQ(ierr);

    // warm start
    ierr = registration->SetReferenceImage(m); CHKERRQ(ierr);
    ierr = registration->SolveForwardProblem(NULL, m); CHKERRQ(ierr);
    ierr = registration->EvaluateGradient(NULL, NULL); CHKERRQ(ierr);
    ierr = registration->HessianMatVec(NULL, NULL); CHKERRQ(ierr);

    // reset all timers and counters
    ierr = opt->ResetTimers(); CHKERRQ(ierr);
    ierr = opt->ResetCounters(); CHKERRQ(ierr);

    IntType n = opt->NumRepeats();
    ierr = opt->StartTimer(reg::T2SEXEC); CHKERRQ(ierr);
    runtime = -MPI_Wtime();
    for (IntType i = 0; i < n; ++i) {
        if (opt->m_Verbosity > 1) {
            ss  << "hessian matvec evaluation "<< std::setw(3)
                << i << " of " << std::setw(3) << n;
            ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();
        }
        ierr = registration->HessianMatVec(NULL, NULL); CHKERRQ(ierr);
    }
    runtime += MPI_Wtime();
    ierr = opt->StopTimer(reg::T2SEXEC); CHKERRQ(ierr);

    opt->SetRunTime(runtime);

    if (v != NULL) {delete v; v = NULL;}
    if (vtilde != NULL) {delete vtilde; vtilde = NULL;}
    if (registration != NULL) {delete registration; registration = NULL;}
    if (m != NULL) {ierr = VecDestroy(&m); CHKERRQ(ierr); m = NULL;}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief perform benchmark for forward solver
 *******************************************************************/
PetscErrorCode ComputeErrorForwardSolver(reg::BenchmarkOpt *opt) {
    PetscErrorCode ierr = 0;
    Vec m0 = NULL, m1 = NULL, m0true = NULL;
    reg::VecField* v = NULL;
    reg::ReadWriteReg* readwrite = NULL;
    reg::CLAIRE* registration = NULL;
    ScalarType val, val0, relval;
    std::stringstream ss;
    double runtime;
    PetscFunctionBegin;

    // make sure we do not store time history
    opt->m_RegFlags.runinversion = false;

    ierr = ComputeSyntheticData(m0true, opt); CHKERRQ(ierr);
    ierr = VecDuplicate(m0true, &m0); CHKERRQ(ierr);
    ierr = VecDuplicate(m0true, &m1); CHKERRQ(ierr);

    ierr = ComputeSyntheticData(v, opt); CHKERRQ(ierr);

    // allocate class for io
    try {readwrite = new reg::ReadWriteReg(opt);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }

    try {registration = new reg::CLAIRE(opt);}
    catch (std::bad_alloc& err) {
        ierr = reg::ThrowError(err); CHKERRQ(ierr);
    }
    ierr = registration->SetControlVariable(v); CHKERRQ(ierr);

    ierr = reg::DbgMsg("computing error for forward solver"); CHKERRQ(ierr);

    // warm start
    ierr = opt->StartTimer(reg::T2SEXEC); CHKERRQ(ierr);
    runtime = -MPI_Wtime();
    ierr = registration->SolveForwardProblem(m1, m0true); CHKERRQ(ierr);
    ierr = v->Scale(-1.0); CHKERRQ(ierr); CHKERRQ(ierr);
    ierr = registration->SetControlVariable(v); CHKERRQ(ierr);
    ierr = registration->SolveForwardProblem(m0, m1); CHKERRQ(ierr);
    runtime += MPI_Wtime();
    ierr = opt->StopTimer(reg::T2SEXEC); CHKERRQ(ierr);
    opt->SetRunTime(runtime);

    ierr = VecAXPY(m0, -1.0, m0true); CHKERRQ(ierr);
    ierr = VecNorm(m0, NORM_2, &val); CHKERRQ(ierr);
    ierr = VecNorm(m0true, NORM_2, &val0); CHKERRQ(ierr);

    /*
    ierr = readwrite->Write(m0, "m0.nc"); CHKERRQ(ierr);
    ierr = readwrite->Write(m1, "m1.nc"); CHKERRQ(ierr);
    ierr = readwrite->Write(m0true, "m0true.nc"); CHKERRQ(ierr);
    */
    relval = val;
    relval /= val0 > 0.0 ? val0 : 1.0;


//    if (opt->m_Verbosity > 0) {
        ss << "numerical error: "<< std::scientific << relval
           << " (absolute " << val << ")";
        ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.clear(); ss.str(std::string());
//    }

    ss << "GPU compute time:"<< std::scientific << opt->m_GPUtime;
    ierr = reg::DbgMsg(ss.str()); CHKERRQ(ierr);
    ss.clear(); ss.str(std::string());

    if (registration != NULL) {delete registration; registration = NULL;}
    if (readwrite != NULL) {delete readwrite; readwrite = NULL;}
    if (m0 != NULL) {ierr = VecDestroy(&m0); CHKERRQ(ierr); m0 = NULL;}
    if (m1 != NULL) {ierr = VecDestroy(&m1); CHKERRQ(ierr); m1 = NULL;}
    if (m0true != NULL) {ierr = VecDestroy(&m0true); CHKERRQ(ierr); m0true = NULL;}
    if (v != NULL) {delete v; v = NULL;}


    PetscFunctionReturn(ierr);
}





/********************************************************************
 * @brief compute synthetic image
 *******************************************************************/
PetscErrorCode ComputeSyntheticData(Vec& m, reg::BenchmarkOpt* opt) {
    PetscErrorCode ierr = 0;
    ScalarType *p_m = NULL;
    ScalarType hx[3], x1, x2, x3;
    IntType i, nl, ng;
    PetscFunctionBegin;

    opt->Enter(__func__);

    // get local and global size
    nl = opt->m_Domain.nl;
    ng = opt->m_Domain.ng;

    // get grid size
    hx[0] = opt->m_Domain.hx[0];
    hx[1] = opt->m_Domain.hx[1];
    hx[2] = opt->m_Domain.hx[2];

    // allocate data
    if (m == NULL) {
        ierr = reg::VecCreate(m, nl, ng); CHKERRQ(ierr);
    }

    ierr = VecGetArray(m, &p_m); CHKERRQ(ierr);
    for (IntType i1 = 0; i1 < opt->m_Domain.isize[0]; ++i1) {  // x1
        for (IntType i2 = 0; i2 < opt->m_Domain.isize[1]; ++i2) {  // x2
            for (IntType i3 = 0; i3 < opt->m_Domain.isize[2]; ++i3) {  // x3
                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + opt->m_Domain.istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + opt->m_Domain.istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + opt->m_Domain.istart[2]);

                // compute linear / flat index
                i = reg::GetLinearIndex(i1, i2, i3, opt->m_Domain.isize);

                p_m[i] =  (PetscSinReal(x1)*PetscSinReal(x1)
                          + PetscSinReal(x2)*PetscSinReal(x2)
                          + PetscSinReal(x3)*PetscSinReal(x3))/3.0;


            }  // i1
        }  // i2
    }  // i3
    ierr = VecRestoreArray(m, &p_m); CHKERRQ(ierr);

    opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief compute synthetic velocity field
 *******************************************************************/
PetscErrorCode ComputeSyntheticData(reg::VecField*& v, reg::BenchmarkOpt* opt) {
    PetscErrorCode ierr = 0;
    ScalarType *p_v1 = NULL, *p_v2 = NULL, *p_v3 = NULL;
    ScalarType hx[3], x1, x2, x3;
    IntType i, vcase = 1;
    PetscFunctionBegin;

    opt->Enter(__func__);

    // get grid size
    hx[0] = opt->m_Domain.hx[0];
    hx[1] = opt->m_Domain.hx[1];
    hx[2] = opt->m_Domain.hx[2];

    // allocate velocity field
    if (v == NULL) {
        try {v = new reg::VecField(opt);}
        catch (std::bad_alloc& err) {
            ierr = reg::ThrowError(err); CHKERRQ(ierr);
        }
    }

    // ierr = v->GetArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = VecGetArray(v->m_X1, &p_v1); CHKERRQ(ierr);
    ierr = VecGetArray(v->m_X2, &p_v2); CHKERRQ(ierr);
    ierr = VecGetArray(v->m_X3, &p_v3); CHKERRQ(ierr);

    for (IntType i1 = 0; i1 < opt->m_Domain.isize[0]; ++i1) {  // x1
        for (IntType i2 = 0; i2 < opt->m_Domain.isize[1]; ++i2) {  // x2
            for (IntType i3 = 0; i3 < opt->m_Domain.isize[2]; ++i3) {  // x3

                // compute coordinates (nodal grid)
                x1 = hx[0]*static_cast<ScalarType>(i1 + opt->m_Domain.istart[0]);
                x2 = hx[1]*static_cast<ScalarType>(i2 + opt->m_Domain.istart[1]);
                x3 = hx[2]*static_cast<ScalarType>(i3 + opt->m_Domain.istart[2]);

                // compute linear / flat index
                i = reg::GetLinearIndex(i1, i2, i3, opt->m_Domain.isize);

                if (vcase == 0) {
                    // compute the velocity field
                    p_v1[i] = 0.5*PetscSinReal(x3)*PetscCosReal(x2)*PetscSinReal(x2);
                    p_v2[i] = 0.5*PetscSinReal(x1)*PetscCosReal(x3)*PetscSinReal(x3);
                    p_v3[i] = 0.5*PetscSinReal(x2)*PetscCosReal(x1)*PetscSinReal(x1);
                } else if (vcase == 1) {
                    // compute the velocity field
                    p_v1[i] = PetscSinReal(x3)*PetscCosReal(x2)*PetscSinReal(x2);
                    p_v2[i] = PetscSinReal(x1)*PetscCosReal(x3)*PetscSinReal(x3);
                    p_v3[i] = PetscSinReal(x2)*PetscCosReal(x1)*PetscSinReal(x1);
                } else if (vcase == 2) {
                    // compute divergence freee velocity field
                    p_v1[i] = PetscCosReal(x2)*PetscCosReal(x3);
                    p_v2[i] = PetscSinReal(x3)*PetscSinReal(x1);
                    p_v3[i] = PetscCosReal(x1)*PetscCosReal(x2);
                } else if (vcase == 3) {
                    p_v1[i] = 0.5;
                    p_v2[i] = 1;
                    p_v3[i] = 0.5;
                } else if (vcase == 4) {
                    p_v1[i] = 0.0;
                    p_v2[i] = 0.0;
                    p_v3[i] = 0.0;
                }
            }  // i1
        }  // i2
    }  // i3

    // ierr = v->RestoreArrays(p_v1, p_v2, p_v3); CHKERRQ(ierr);
    ierr = VecRestoreArray(v->m_X1, &p_v1); CHKERRQ(ierr);
    ierr = VecRestoreArray(v->m_X2, &p_v2); CHKERRQ(ierr);
    ierr = VecRestoreArray(v->m_X3, &p_v3); CHKERRQ(ierr);


    opt->Exit(__func__);

    PetscFunctionReturn(ierr);
}

