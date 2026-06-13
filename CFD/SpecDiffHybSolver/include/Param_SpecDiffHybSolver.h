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
//! @file      Param_SpecDiffHybSolver.h
//! @brief     Record paramters of SpecDiffHybSolver.
//! @author    zipzhang.
#pragma once
#include "LIB_Macro.h"
#include "Param_SpecSolver.h"

namespace PHSPACE
{

class Param_SpecDiffHybSolver : public Param_SpecSolver
{
public:
    LIB_EXPORT Param_SpecDiffHybSolver();

    LIB_EXPORT ~Param_SpecDiffHybSolver();

public:
    //! Init all parameters.
    LIB_EXPORT void Init();

    RDouble GetPressureGradient() const;

private:
    RDouble pressureGradient;
};

#include "Param_SpecDiffHybSolver.hxx"
}