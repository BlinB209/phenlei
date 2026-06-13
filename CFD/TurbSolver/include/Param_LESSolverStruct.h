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
//! @file      Param_LESSolverStruct.h
//! @brief     Record paramters of LESSolverStruct.
//! @author    Zhang Zipei.
#pragma once
#include "Param_LESSolver.h"

namespace PHSPACE
{

class Param_LESSolverStruct : public Param_LESSolver
{
public:
    LIB_EXPORT Param_LESSolverStruct();

    LIB_EXPORT ~Param_LESSolverStruct();

public:

    LIB_EXPORT void Init();

private:

};

}