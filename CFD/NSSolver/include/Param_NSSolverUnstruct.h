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
//! @file      Param_NSSolverUnstruct.h
//! @brief     Record parameters of NS Solver for unstructured grid.
//! @author    Bell, Zhang Jian, Wan Yunbo, Meng Liyuan.

#pragma once
#include "Param_NSSolver.h"

namespace PHSPACE
{

class Param_NSSolverUnstruct : public Param_NSSolver
{
public:
    LIB_EXPORT Param_NSSolverUnstruct();

    LIB_EXPORT ~Param_NSSolverUnstruct();

public:

    //! Init all parameters.
    LIB_EXPORT void Init();

    //! Get Spatial discretization scheme of Unstructured grid solver.
    string GetSpaceSecheme () const;

    string GetGradientName() const;

    int GetLimitVector() const;

    int GetLimitVariables() const;

    int GetSliceAxis() const;

    RDouble GetSlicePostion() const;

    RDouble GetMgCFLScale() const;

    RDouble GetMgCorrectionLimit() const;

    int GetLimiterType() const;

    void SetLimiterType(int type){this->limiterType = type;}

    //! Get catalyticCoef.
    RDouble GetCatalyticCoef() const;

    //! Get ChemicalSpectrumRadiusCoef.
    RDouble GetChemicalSpectrumRadiusCoef() const;

    //! Get chemicalRelaxCorf.
    RDouble GetChemicalRelaxCorf() const;

    //! Get staticPressureRelaxCorf.
    RDouble GetStaticPressureRelaxCorf() const;
     
private:
    string uns_scheme_name;

    string gradientName;

    int limitVector;

    int limitVariables;

    int sliceAxis;

    RDouble slicePostion;

    RDouble mgCFLScale;

    RDouble mgCorrectionLimit;

    int limiterType;

    RDouble catalyticCoef;

    RDouble chemicalSpectrumRadiusCoef;

    RDouble chemicalRelaxCorf;

    RDouble staticPressureRelaxCorf;
};

#include "Param_NSSolverUnstruct.hxx"
}
