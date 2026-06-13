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
//! @file      Param_TurbSolver.h
//! @brief     Record paramters of TurbSolver.
//! @author    Bell, Zhang Jian, Wan Yunbo, Meng Liyuan.

#pragma once
#include "Param_CFDSolver.h"

namespace PHSPACE
{

class Param_TurbSolver : public Param_CFDSolver
{
public:
    LIB_EXPORT Param_TurbSolver();

    LIB_EXPORT ~Param_TurbSolver();

public:

    LIB_EXPORT void Init();

    RDouble GetTurbCFLScale() const;

    int GetDESType() const;

    int GetTransitionType() const;

    int GetSAProductType() const;

    int GetSSTProductType() const;

    int GetNTurbulenceEquation() const;

    RDouble GetEddyViscosityLimit() const;

    int GetWallFunctionType() const;

    RDouble GetSST_a1() const;

    RDouble GetSST_alphaw1() const;

    RDouble GetSST_alphaw2() const;

    RDouble GetSST_beta1() const;

    RDouble GetSST_beta2() const;

    RDouble GetSST_betaStar() const;

    RDouble GetSA_cv1_cube() const;

    RDouble GetSSTProductLimit() const;

    RDouble GetKW_sigma() const;

    RDouble GetKW_sigmaW() const;

    RDouble GetKW_sigmaW1() const;

    RDouble GetKW_sigmaW2() const;

    RDouble GetKW_sigmaK() const;

    RDouble GetKW_sigmaK1() const;

    RDouble GetKW_sigmaK2() const;

    RDouble GetFreeStreamViscosity() const;

    RDouble * GetFreeStreamTurbVar() const;

    int GetKindOfTurbSource() const;

    RDouble GetInviscidSpectrumRadiusCoef() const;

    int GetWennSchemeFlag() const;

    //! Return using original ( = 0 ) or conservative variables ( = 1 )\n.
    //! 0: original form.\n
    //! 1: conservative form.\n
    int UsingConservativeForm() const;
private:

    RDouble turbCFLScale;

    int DESType;

    int transitionType;
    int SAProductType;

    int SSTProductType;

    int nTurbulenceEquation;

    RDouble eddyViscosityLimit;

    int wallFunctionType;

    int isTurbulenceModelCorrect;

    RDouble SST_a1;

    RDouble SST_alphaw1;

    RDouble SST_alphaw2;

    RDouble SST_beta1;

    RDouble SST_beta2;

    RDouble SST_betaStar;

    RDouble SA_cv1_cube;

    RDouble SSTProductLimit;

    RDouble KW_sigma;

    RDouble KW_sigmaW;

    RDouble KW_sigmaW1;

    RDouble KW_sigmaW2;

    RDouble KW_sigmaK;

    RDouble KW_sigmaK1;

    RDouble KW_sigmaK2;

    RDouble turbulenceIntensity;

    RDouble freeStreamViscosity;

    RDouble * freeStreamTurbVar;

    int kindOfTurbSource;
    int isWennScheme;

    //! Used to compatible both original ( = 0 ) and conservative variables ( = 1 )\n.
    //! 0: original form.\n
    //! 1: conservative form.\n
    int rhoFlag;

    RDouble inviscidSpectrumRadiusCoef;
};

#include "Param_TurbSolver.hxx"
}