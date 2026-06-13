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
//! @file      Param_TurbSolverUnstruct.h
//! @brief     Record paramters of TurbSolverstruct.
//! @author    Bell, Zhang Jian, Wan Yunbo, Meng Liyuan.

#pragma once
#include "Param_TurbSolver.h"

namespace PHSPACE
{

class Param_TurbSolverUnstruct : public Param_TurbSolver
{
public:
    LIB_EXPORT Param_TurbSolverUnstruct();

    LIB_EXPORT ~Param_TurbSolverUnstruct();

public:

    LIB_EXPORT void Init();

private:

};

}

