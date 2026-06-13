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
//! @file      Param_SpecSolver.h
//! @brief     Record paramters of SpecSolver.
//! @author    zipzhang.
#pragma once
#include "LIB_Macro.h"
#include "Precision.h"

namespace PHSPACE
{

class Param_SpecSolver
{
public:
    LIB_EXPORT Param_SpecSolver();

    LIB_EXPORT virtual ~Param_SpecSolver();

public:
    //! Init all parameters.
    LIB_EXPORT virtual void Init();

    int GetStatisticMethod() const;

    int GetStartStatisticStep() const;

    int GetMaxSimuStep() const;

    int GetIntervalStepRes() const;

    int GetIntervalStepStatistic() const;

    int GetIntervalStepFlow() const;

    RDouble GetRefReNumber() const;

    const string & GetViscousName() const;

    int GetViscousType() const;

    const string & GetFilterName() const;

    RDouble GetSmagConstant() const;

    RDouble GetTestFilterScale() const;

    RDouble GetTimeStep() const;

    //RDouble GetPressureGradient() const;

    int GetSpectrumType() const;

protected:
    int statisticMethod;

    int startStatisticStep;

    int maxSimuStep;

    int intervalStepRes;

    int intervalStepStatistic;

    int intervalStepFlow;

    RDouble refReNumber;

    string viscousName;

    int viscousType;

    string filterName;

    RDouble smagConstant;

    RDouble testFilterScale;

    RDouble timeStep;

    //RDouble pressureGradient;

private:
    int spectrumType;
};

#include "Param_SpecSolver.hxx"
}