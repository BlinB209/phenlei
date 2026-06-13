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
//! @file      Param_INCompSolverUnstruct.h
//! @brief     parameter files.
//! @author    WanYunbo, Bell.

#pragma once
#include "Constants.h"
#include "Param_CFDSolver.h"

namespace PHSPACE
{

class Param_INCompSolverUnstruct : public Param_CFDSolver
{
public:
    LIB_EXPORT Param_INCompSolverUnstruct();

    LIB_EXPORT ~Param_INCompSolverUnstruct();

public:
    LIB_EXPORT void Init();
    
    void setInComEquaMap();

private:

};

#include"Param_INCompSolverUnstruct.hxx"

}
