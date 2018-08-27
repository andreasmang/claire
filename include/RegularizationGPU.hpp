#ifndef __REGULARIZATIONGPU_HPP__
#define __REGULARIZATIONGPU_HPP__

#include "typedef.hpp"

void EvaluateGradientH1SN_GPU(ComplexType *v1hat, ComplexType *v2hat, ComplexType *v3hat, IntType ostart[3], IntType nl[3], IntType nx[3], ScalarType opfact);

#endif
