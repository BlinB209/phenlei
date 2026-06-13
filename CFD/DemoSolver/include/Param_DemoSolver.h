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
//! @file      Param_DemoSolver.h
//! @brief     Record parameters of demo solver.
//! @author    Meng Liyuan.

#pragma once
#include "LIB_Macro.h"
#include "Param_CFDSolver.h"

namespace PHSPACE
{

class Param_DemoSolver : public Param_CFDSolver
{
public:
    LIB_EXPORT Param_DemoSolver();

    LIB_EXPORT ~Param_DemoSolver();

public:

    //! Init all parameters.
    LIB_EXPORT void Init();

    //! Get wall temperature.
    RDouble GetWallTemperature() const;
    //! Get angle of attacked.
    RDouble GetAoA() const;
    //! Get angle of slide.
    RDouble GetAngleOfSlide () const;

private:
    //! Wall temperature.
    RDouble wallTemperature;

    //! Angle of attacked.
    double AoA;

    //! Angle of slide.
    RDouble angleSlide;
};

#include "Param_DemoSolver.hxx"
}


