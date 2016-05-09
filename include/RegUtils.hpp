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


#ifndef _REGUTILS_H_
#define _REGUTILS_H_

// global includes
#include <iomanip>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <omp.h>

// local includes
#include "petsc.h"
#include "accfft.h"

#define IntType PetscInt
//#define ScalarType PetscScalar

//#define IntType int
#define ScalarType double


//#include "VecField.h"


namespace reg
{


/*! assert (PETSc interface) */
PetscErrorCode Assert(bool, std::string);

/*! throw error (PETSc interface) */
PetscErrorCode ThrowError(std::string);

/*! check if file exists */
bool FileExists(const std::string&);

/*! display message (PETSc interface) */
PetscErrorCode Msg(std::string);

/*! display warning message (PETSc interface) */
PetscErrorCode WrngMsg(std::string);

/*! display dgb message (PETSc interface) */
PetscErrorCode DbgMsg(std::string);

/*! display scalar field */
PetscErrorCode VecView(Vec);

PetscErrorCode Rescale(Vec, ScalarType, ScalarType);
PetscErrorCode GetFileName(std::string&,std::string);
std::vector<unsigned int> String2Vec( const std::string & );

/*! display vector field */
//PetscErrorCode VecView(VecField*);



/*! write data to file */
//PetscErrorCode Write2File(Vec,std::string);


/*! write vector field to file */
//PetscErrorCode Write2File(VecField*,std::string);




/********************************************************************
 * Name: GetLinearIndex
 * Description: map 3d index to linear index (accfft style)
 *******************************************************************/
inline IntType GetLinearIndex(IntType i[3], IntType isize[3])
{
    // row major order (ACCFFT)$
    return i[0]*isize[1]*isize[2]+i[1]*isize[2]+i[2];
}



/********************************************************************
 * Name: GetLinearIndex
 * Description: map 3d index to linear index (accfft style)
 *******************************************************************/
inline IntType GetLinearIndex(IntType i, IntType j, IntType k, IntType isize[3])
{
    // row major order (ACCFFT)$
    return i*isize[1]*isize[2]+j*isize[2]+k;
};


/********************************************************************
 * Name: GetLinearIndex
 * Description: map 3d index to linear index (accfft style)
 *******************************************************************/


/********************************************************************
 * Name: CheckWaveNumbersInv
 * Description: check wave numbers
 *******************************************************************/
inline void CheckWaveNumbersInv(long int w[3],unsigned int n[3])
{

    if     (w[0] >  n[0]/2) w[0]-=n[0];
    if     (w[1] >  n[1]/2) w[1]-=n[1];
    if     (w[2] >  n[2]/2) w[2]-=n[2];

};




/********************************************************************
 * Name: CheckWaveNumbers
 * Description: check wave numbers
 *******************************************************************/
inline void CheckWaveNumbers(long int w[3],unsigned int n[3])
{

    if     (w[0] >  n[0]/2) w[0]-=n[0];
    else if(w[0] == n[0]/2) w[0] = 0;

    if     (w[1] >  n[1]/2) w[1]-=n[1];
    else if(w[1] == n[1]/2) w[1] = 0;

    if     (w[2] >  n[2]/2) w[2]-=n[2];
    else if(w[2] == n[2]/2) w[2] = 0;

};







} // end of namespace

#endif
