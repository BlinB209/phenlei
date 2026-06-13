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
//! @file      HOSolverUnstructParam.h
//! @brief     Record paramters of NS Solver.
//! @author    Zhang Jian, Wan yunbo, Xu gang.

#pragma once
#include "LIB_Macro.h"
#include "Param_CFDSolver.h"
#include "Precision.h"
namespace PHSPACE
{

class HOSolverUnstruct;

class HOSolverUnstructParam : public Param_CFDSolver
{
public:
    LIB_EXPORT HOSolverUnstructParam();

    LIB_EXPORT ~HOSolverUnstructParam();

    friend class HOSolverUnstruct;

public:
    //! Init all parameters
    LIB_EXPORT void Init();

    //! Get restart flow file name
    const string & GetRestartNSFile() const;

    //! Get res file name
    const string & GetResSaveFile() const;

    //! Get flow field name
    const string & GetFlowFieldTecFile() const;

    //! Get interval steps of dumping force
    int GetIntervalStepForce() const;

    //! Get wall temperature
    RDouble GetWallTemperature() const;

    //! Get number of temperature model
    int GetNTemperatureModel () const;

    //! Get if carry out precondition or not
    int GetIfLowSpeedPrecon() const;

    //! Get attacked of angle
    double GetAoA() const;

    //! Get the entropy fix coefficient 1 of Roe method
    RDouble GetRoeEntropyFixCoef1() const;
    
    //! Get the entropy fix coefficient 2 of Roe method
    RDouble GetRoeEntropyFixCoef2() const;

    //! Get angle of slide
    RDouble GetAngleOfSlide () const;
    
    //! Get inviscos flux type
    string GetSpaceSecheme () const;

    //! Get p_multigrid type
    int GetPMultiGrid () const;

    RDouble GetR() const;

private:
    string restartNSFile;
    string resSaveFile;
    string flowFieldTecFile;
    int intervalStepForce;
    RDouble wallTemperature;
    int nTemperatureModel;
    int ifLowSpeedPrecon;
    RDouble RoeEntropyFixCoef1;
    double AoA;
    RDouble RoeEntropyFixCoef2;
    RDouble angleSlide;

    int dgSolOrder, fvSolOrder;
    int dgSolDof, fvSolDof;
    int pMultiGrid;
    string uns_scheme_name;

    RDouble R; // p = rho R T
};

inline string HOSolverUnstructParam::GetSpaceSecheme() const
{
    return uns_scheme_name;
}

inline int HOSolverUnstructParam::GetPMultiGrid() const
{
    return pMultiGrid;
}

#include "HOSolverUnstructParam.hxx"
}