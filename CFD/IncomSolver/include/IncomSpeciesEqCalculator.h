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
//! @file      IncomSpeciesEqCalculator.h
//! @brief     Species solver.
//! @author    WanYunbo, Bell, XuGang.

#pragma once

#include "IncomScalarEqCalculator.h"
#include "IncomGas.h"

namespace PHSPACE
{
class IncomSpeciesEqCalculator : public IncomScalarEqCalculator
{
public:
    IncomSpeciesEqCalculator();
    ~IncomSpeciesEqCalculator();

public:
    void InitFlowAsRestart(Grid *gridIn);
    void IncompressibleInitial(Grid *grid);

    virtual void UpdateProperties(Grid *grid);
    void UpdateAfterIterloop(Grid *grid);
    void UpdateSpBoundary(Grid *grid);

    void GetResidual(Grid *grid, vector<RDouble> &res);
    void InitialUnsteadyVar(Grid *grid);
    void UpdateUnsteadyVariable(Grid *grid);
    void UpdateBCValue(Grid *grid);
    void SetDiffusionCoeff(Grid *gridIn);
};

}