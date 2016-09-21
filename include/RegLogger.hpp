/**
 *  Copyright (c) 2015-2016.
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


#ifndef _REGLOGGER_H_
#define _REGLOGGER_H_

#include <string>
#include <vector>
#include <fstream>

// local includes
#include "RegUtils.hpp"




namespace reg {



enum LogFlag{
    LOGRES,
    LOGKSPRES,
    LOGLOAD,
    NUMLOGFLAGS
};


class RegLogger {
 public:
    typedef RegLogger Self;
    RegLogger();
    ~RegLogger();

    inline void Enable(LogFlag i){this->m_Enabled[i] = true;}
    inline bool IsEnabled(LogFlag i){return this->m_Enabled[i];}
    inline void Disable(LogFlag i){this->m_Enabled[i] = false;}

    inline void SetKSPResidual(int i, ScalarType value) {
        this->m_KSPResidual.push_back(value);
        this->m_KSPIterations.push_back(i);
    }

    inline void SetResidual(int i, ScalarType value){this->m_Residual[i] = value;}

    PetscErrorCode Write(std::string);


 private:
    PetscErrorCode Initialize();

    bool m_Enabled[NUMLOGFLAGS];

    std::vector<ScalarType> m_KSPResidual;   ///< residual of krylov method
    std::vector<int> m_KSPIterations;        ///< iterations of krylov method
    ScalarType m_Residual[4];
};




}  // namespace reg




#endif  // _REGLOGGER_H_
