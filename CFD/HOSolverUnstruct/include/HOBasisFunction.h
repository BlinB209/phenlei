//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      HOBasisFunction.h
//! @brief     basis function, such as monomial, orthogonalized etc.
//! @author    Li Ming, Gong Xiaoquan, Wan Yunbo, Xu Gang.
//! @date      2019-10-17
//! @version   A001

#pragma once

#include "Precision.h"

namespace HOUnstruct
{

PHSPACE::RDouble MonomialBasis(const PHSPACE::RDouble * refCoord, unsigned int i);
PHSPACE::RDouble MonomialBasisGrad(const PHSPACE::RDouble * refCoord, unsigned int i, unsigned int j);

PHSPACE::RDouble MonomialBasis2d(const PHSPACE::RDouble * refCoord, unsigned int i);

void TestBasisFunction();

}