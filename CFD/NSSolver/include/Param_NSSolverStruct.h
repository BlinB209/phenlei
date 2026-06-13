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
//! @file      Param_NSSolverStruct.h
//! @brief     Record parameters of NS Solver for structured grid.
//! @author    Bell, Zhang Jian, Wan Yunbo, Meng Liyuan.

#pragma once
#include "Param_NSSolver.h"

namespace PHSPACE
{

class Param_NSSolverStruct : public Param_NSSolver
{
public:
    LIB_EXPORT Param_NSSolverStruct();

    LIB_EXPORT ~Param_NSSolverStruct();

public:

    //! Init all parameters.
    LIB_EXPORT void Init();

    //! Get parameter of MUSCL interpolations.
    RDouble GetMUSCLCoefXk() const;

    //! Get the limiter parameter of MUSCL.
    RDouble GetMUSCLCoefXb() const;

    //! Get catalyticCoef.
    RDouble GetCatalyticCoef() const;

    //! Get ChemicalSpectrumRadiusCoef.
    RDouble GetChemicalSpectrumRadiusCoef() const;

    //! Get ChemicalSpectrumRadiusCoef.
    RDouble GetViscousSpectrumRadiusCoef() const;

    //! Get GetInviscidSpectrumRadiusCoef.
    RDouble GetInviscidSpectrumRadiusCoef() const;

    //! Get chemicalRelaxCorf.
    RDouble GetChemicalRelaxCorf() const;

    //! Get staticPressureRelaxCorf.
    RDouble GetStaticPressureRelaxCorf() const;

    //! Get inviscidSchemeName.
    string GetInviscidSchemeName() const;

    //! Get iniSpeedCoef.
    RDouble GetIniSpeedCoef() const;

    //! Get iniSpeedMode.
    int GetIniSpeedMode() const;

private:
    RDouble MUSCLCoefXk;
    RDouble MUSCLCoefXb;

    RDouble catalyticCoef;

    RDouble chemicalSpectrumRadiusCoef;

    RDouble viscousSpectrumRadiusCoef;

    RDouble inviscidSpectrumRadiusCoef;

    RDouble chemicalRelaxCorf;

    RDouble staticPressureRelaxCorf;

    string inviscidSchemeName;

    RDouble iniSpeedCoef;

    int iniSpeedMode;
};

#include "Param_NSSolverStruct.hxx"
}
