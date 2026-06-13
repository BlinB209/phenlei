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
//! @file      Param_LESSolver.h
//! @brief     Record paramters of LESSolver.
//! @author    Zhang Zipei.
#pragma once
#include "LIB_Macro.h"
#include "Param_CFDSolver.h"

namespace PHSPACE
{

class Param_LESSolver : public Param_CFDSolver
{
public:
    LIB_EXPORT Param_LESSolver();

    LIB_EXPORT ~Param_LESSolver();

public:

    LIB_EXPORT void Init();

    RDouble GetEddyViscosityLimit() const;

    int GetWallDampingFunctionType() const;

    int GetSubgridScaleModel() const;

    RDouble * GetTestFilterCoef() const;

    int GetTestFilterWidth() const;

    RDouble GetSmagConstant() const;

    RDouble GetIsotropicConstant() const;

    RDouble GetYploge() const;

    int GetDeltaFunctionType() const;

    RDouble GetRgam() const;

    int GetTurbViscousCutType() const;

    RDouble GetTestFilterScale() const;

    RDouble GetWaleConstant() const;

    RDouble GetSigmaConstant() const;

    int * GetFilterDirection() const;

    int * GetAverageDirection() const;

    int GetAverageWidth() const;

    RDouble GetoPrandtlTurbulence() const;

    RDouble GetPrandtlTurbulence() const;

    RDouble GetFreeStreamViscosity() const;

    string GetGradientName() const;

private:

    RDouble eddyViscosityLimit;

    int wallDampingFunctionType;

    int subgridScaleModel;

    RDouble * testFilterCoef;

    int testFilterWidth;

    RDouble smagConstant;

    RDouble isotropicConstant;

    RDouble yploge;

    int deltaFunctionType;

    RDouble rgam;

    int turbViscousCutType;

    RDouble testFilterScale;

    RDouble waleConstant;

    RDouble sigmaConstant;

    int * filterDirection;

    int * averageDirection;

    int averageWidth;

    RDouble oPrandtlTurbulence;

    RDouble prandtlTurbulence;

    RDouble freeStreamViscosity;

    string gradientName;
};

#include "Param_LESSolver.hxx"
}