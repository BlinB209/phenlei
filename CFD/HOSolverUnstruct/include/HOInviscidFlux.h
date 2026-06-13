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
//! @file      HOInviscidFlux.h
//! @brief     basis function, such as monomial, orthogonalized etc.
//! @author    gongxq, liming, Wan yunbo, Xu gang.
//! @date      2019-11-15
//! @version   A001

#pragma once

#include "Precision.h"

using PHSPACE::RDouble;
namespace HOUnstruct
{
    void Vanleer_Scheme(int totalGaussPoint, RDouble ** ql, RDouble ** qr, RDouble *areax, RDouble *areay, RDouble *areaz, RDouble gama, RDouble ** flux);    
}