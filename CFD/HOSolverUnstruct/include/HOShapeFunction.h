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
//! @file      HOShapeFunction.h
//! @brief     shape function of all element types supported (such as tri trip2, quad etc.).
//! @author    Li Ming, Gong Xiaoquan, Wan Yunbo, Xu Gang.
//! @date      2019-10-15
//! @version   A001

#pragma once

#include "Precision.h"

namespace HOUnstruct
{

void ShapeFunctionOfTri(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);
void ShapeFunctionOfTriGrad(PHSPACE::RDouble * ret);

void ShapeFunctionOfQuad(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);
void ShapeFunctionOfQuadGrad(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);

void ShapeFunctionOfTetr(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);
void ShapeFunctionOfTetrGrad(PHSPACE::RDouble * ret);

void ShapeFunctionOfHex(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);
void ShapeFunctionOfHexGrad(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);

void ShapeFunctionOfPrism(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);
void ShapeFunctionOfPrismGrad(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);

void ShapeFunctionOfPyra(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);
void ShapeFunctionOfPyraGrad(const PHSPACE::RDouble * refCoord, PHSPACE::RDouble * ret);

}