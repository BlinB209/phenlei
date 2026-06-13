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
//! @file      Glb_Dimension.h
//! @brief     It defines the global grid dimension.
//! @author    Bell.

#pragma once
#include "LIB_Macro.h"

namespace PHSPACE
{
const int ONE_D   = 1;
const int TWO_D   = 2;
const int THREE_D = 3;

//! Return the global grid dimension.
//! - # 2, two dimensional.
//! - # 3, three dimensional.
LIB_EXPORT int GetDim();

//! Assign the global grid dimension.
LIB_EXPORT void SetDim(int dim);
}
