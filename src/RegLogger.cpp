/** *  Copyright (c) 2015-2016.
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
 *  along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#ifndef _REGLOGGER_CPP_
#define _REGLOGGER_CPP_

// global includes
#include <string>
#include <vector>

// local includes
#include "RegLogger.hpp"




namespace reg {




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegLogger"
RegLogger::RegLogger() {
    this->Initialize();
}




/********************************************************************
 * @brief standard destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegLogger"
RegLogger::~RegLogger() {
    std::vector<ScalarType>().swap(this->m_KSPResidual);
    std::vector<int>().swap(this->m_KSPIterations);
}




/********************************************************************
 * @brief function to initialize class variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode RegLogger::Initialize() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    for (int i = 0; i < NUMLOGFLAGS; ++i) {
        this->m_Enabled[i] = false;
    }

    this->m_Residual[0] = 0.0;
    this->m_Residual[1] = 0.0;
    this->m_Residual[2] = 0.0;
    this->m_Residual[3] = 0.0;

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief write logs to file
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Write"
PetscErrorCode RegLogger::Write(std::string path) {
    PetscErrorCode ierr = 0;
    int n, rank, nnum;
    std::string fn;
    std::ofstream logwriter;
    std::stringstream ss;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    nnum = 20;

    // log residual
    if (this->m_Enabled[LOGKSPRES]) {
        if (rank == 0) {
            // create output file
            fn = path + "cold-krylov-method-residual.log";
            logwriter.open(fn.c_str());
            ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

            n = static_cast<int>(this->m_KSPResidual.size());
            for (int i = 0; i < n; ++i) {
                ss << std::scientific << std::right
                   << std::setw(2) << this->m_KSPIterations[i]
                   << std::setw(20) << this->m_KSPResidual[i];
                logwriter << ss.str() << std::endl;
                ss.str(std::string()); ss.clear();
            }
            logwriter.close();  // close logger
        }
    }

    // log residual
    if (this->m_Enabled[LOGRES]) {
        if (rank == 0) {
            // create output file
            fn = path + "cold-residual.log";
            logwriter.open(fn.c_str());
            ierr = Assert(logwriter.is_open(), "could not open file for writing"); CHKERRQ(ierr);

            ss  << std::scientific << std::left
                << std::setw(20) << "||mR-mT||_2" << std::right
                << std::setw(nnum) << this->m_Residual[0];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(20) << "||mR-mT||_infty" << std::right
                << std::setw(nnum) << this->m_Residual[1];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(20) << "||mR-m1||_2" << std::right
                << std::setw(nnum) << this->m_Residual[2];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(20) << "||mR-m1||_infty" << std::right
                << std::setw(nnum) << this->m_Residual[3];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            this->m_Residual[0] = (this->m_Residual[0] > 0.0) ? this->m_Residual[0] : 1.0;
            this->m_Residual[1] = (this->m_Residual[1] > 0.0) ? this->m_Residual[1] : 1.0;

            ss  << std::scientific << std::left
                << std::setw(20) << "||mR-m1||_2,rel" << std::right
                << std::setw(nnum) <<  this->m_Residual[2]/this->m_Residual[0];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            ss  << std::scientific << std::left
                << std::setw(20) << "||mR-m1||_infty,rel" << std::right
                << std::setw(nnum) <<  this->m_Residual[3]/this->m_Residual[1];
            logwriter << ss.str() << std::endl;
            ss.clear(); ss.str(std::string());

            // close logger
            logwriter.close();
        }
    }

    PetscFunctionReturn(ierr);
}




}  // namespace reg




#endif  // _REGLOGGER_CPP_
