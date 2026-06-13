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
//! @file      Param_TransitionSolver.h
//! @brief     Record paramters of TransitionSolver.
//! @author    He Kun

#pragma once
#include "Param_CFDSolver.h"

namespace PHSPACE
{

class Param_TransitionSolver : public Param_CFDSolver
{
public:

    LIB_EXPORT ~Param_TransitionSolver();

public:

    LIB_EXPORT void Init();

    RDouble GetTransitionCFLScale() const;

    int GetTransitionType() const;

    int GetNTransitionEquation() const;

    RDouble Getca1() const;

    RDouble Getca2() const;

    RDouble Getce1() const;

    RDouble Getce2() const;

    RDouble Getdct() const;

    RDouble Getdf() const;

    RDouble Getcct() const;

    RDouble Gets1() const;

    RDouble Getccf() const;

    RDouble *GetFreeStreamTransitionVar() const;

private:

    RDouble transitionCFLScale;

    int transitionType;

    int nTransitionEquation;

    RDouble ca1;
    RDouble ce1;
    RDouble ca2;
    RDouble ce2;
    RDouble dct;
    RDouble df;
    RDouble cct;
    RDouble s1;
    RDouble ccf;

    RDouble turbulenceIntensity;

    RDouble freeStreamViscosity;

    RDouble *freeStreamTransitionVar;
};

#include "Param_TransitionSolver.hxx"
}