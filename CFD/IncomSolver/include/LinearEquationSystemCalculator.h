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
//! @file      LinearEquationSystemCalculator.h
//! @brief     solve Linear Equation System.
//! @author    WanYunbo, Bell.

#pragma once

#include "Geo_UnstructGrid.h"

namespace PHSPACE
{
class LinearEquationSystemCalculator
{

public:
    // constructor & destructor
    LinearEquationSystemCalculator() {};

    virtual ~LinearEquationSystemCalculator() {};

public:
    virtual void explicitSolve(Grid *gridIn, RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace) {};
    virtual void mathLibSolveInit(Grid *gridIn){};
    virtual void mathLibSolve(Grid *gridIn, int iEquation, RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace) {};
    virtual void mathLibSolveFinalize(){};
};

}

