/**
 *  Copyright (c) 2015-2016.
 *  All rights reserved.
 *  This file is part of PGLISTR library.
 *
 *  PGLISTR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  PGLISTR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PGLISTR.  If not, see <http://www.gnu.org/licenses/>.
 *
*/



#ifndef _DATAREADWRITEREGISTRATION_H_
#define _DATAREADWRITEREGISTRATION_H_


#include "RegOpt.h"
#include "VecField.h"




namespace reg
{




class DataReadWriteRegistration
{

public:

    typedef DataReadWriteRegistration Self;

    DataReadWriteRegistration(void);
    DataReadWriteRegistration(RegOpt*);
    ~DataReadWriteRegistration(void);

    PetscErrorCode Read(Vec,std::string);
    PetscErrorCode Read(VecField*,std::string,std::string,std::string);
    PetscErrorCode ReadBlock(Vec,int*,std::string);
    PetscErrorCode ReadTimeSeries(Vec,std::string);

    PetscErrorCode Write(Vec,std::string);
    PetscErrorCode Write(VecField*,std::string,std::string,std::string);
    PetscErrorCode WriteBlock(Vec,int*,std::string);
    PetscErrorCode WriteTimeSeries(Vec,std::string);


private:

    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();

    PetscErrorCode ReadNetCDF(Vec,std::string);
    PetscErrorCode ReadTimeSeriesNetCDF(Vec,std::string);
    PetscErrorCode ReadBlockNetCDF(Vec,int*,std::string);

    PetscErrorCode WriteNetCDF(Vec,std::string);
    PetscErrorCode WriteTimeSeriesNetCDF(Vec,std::string);
    PetscErrorCode WriteBlockNetCDF(Vec,int*,std::string);


    RegOpt* m_Opt;


};


} // end of name space

#endif // _DATAREADWRITEREGISTRATION_H_
