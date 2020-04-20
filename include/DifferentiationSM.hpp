/*************************************************************************
 *  Copyright (c) 2018.
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

#ifndef _DIFFERENTIATIONSM_HPP_
#define _DIFFERENTIATIONSM_HPP_

#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "Differentiation.hpp"
#include "DifferentiationKernel.hpp"

//#ifdef REG_HAS_CUDA
//#include <cufft.h>
//#endif

namespace reg {




class DifferentiationSM : public Differentiation {
 public:
    typedef Differentiation SuperClass;
    typedef DifferentiationSM Self;
    
    using SuperClass::Gradient;
    using SuperClass::Divergence;
    using SuperClass::Laplacian;

    DifferentiationSM();
    DifferentiationSM(RegOpt*);
    virtual ~DifferentiationSM();

    virtual PetscErrorCode Gradient(ScalarType*, ScalarType*, ScalarType*, const ScalarType*);
    
    virtual PetscErrorCode Laplacian(ScalarType*, const ScalarType*);
    virtual PetscErrorCode Laplacian(ScalarType*, ScalarType*, ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    
    virtual PetscErrorCode Divergence(ScalarType*, const ScalarType*, const ScalarType*, const ScalarType*);
    
    virtual PetscErrorCode RegLapModOp(VecField*, VecField*, ScalarType);
    virtual PetscErrorCode RegLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode RegBiLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode RegTriLapOp(VecField*, VecField*, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode InvRegLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode InvRegBiLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode InvRegTriLapOp(VecField*, VecField*, bool, ScalarType, ScalarType=0.0);
    virtual PetscErrorCode RegTriLapFunc(VecField*, VecField*, ScalarType, ScalarType=0.0);
    
    virtual PetscErrorCode LerayOperator(VecField*, VecField*, ScalarType, ScalarType);
    virtual PetscErrorCode InvRegLerayOp(VecField*, VecField*, ScalarType, ScalarType, ScalarType);
    
    virtual PetscErrorCode GaussianFilter(ScalarType*, const ScalarType*, const ScalarType*);
    
    virtual PetscErrorCode SetFFT(FourierTransform* fft);
    
    virtual PetscErrorCode Restrict(ScalarType*, const ScalarType*, FourierTransform* coarse);
    virtual PetscErrorCode Prolong(ScalarType*, const ScalarType*, FourierTransform* coarse);
    virtual PetscErrorCode RestrictH0(VecField*, VecField*, FourierTransform* coarse, ScalarType beta);
    virtual PetscErrorCode ProlongH0(VecField*, VecField*, FourierTransform* coarse);
    
 protected:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();
    
    PetscErrorCode SetupData(ComplexType* =nullptr, ComplexType* =nullptr, ComplexType* =nullptr);
    
    PetscErrorCode ComputeForwardFFT(VecField*);
    PetscErrorCode ComputeForwardFFT(const ScalarType*,const ScalarType*,const ScalarType*);
    PetscErrorCode ComputeForwardFFT(const ScalarType*);
    PetscErrorCode ComputeInverseFFT(VecField*);
    PetscErrorCode ComputeInverseFFT(ScalarType*,ScalarType*,ScalarType*);
    PetscErrorCode ComputeInverseFFT(ScalarType*);
        
    DifferentiationKernel m_SpectralKernel;
    std::bitset<3> xyz;
//    double timer[NFFTTIMERS];
    int c_grad;
    int c_div;
    
//#ifdef REG_HAS_CUDA
//    cufftHandle *m_planR2C;
//    cufftHandle *m_planC2R;
//#endif
    
    ComplexType *m_XHat[3];
    
    FourierTransform *m_FFT;
};




}  // end of namespace

#endif  // _DIFFERENTIATION_HPP_
