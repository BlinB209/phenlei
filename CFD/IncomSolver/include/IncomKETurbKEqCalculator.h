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
//! @file      IncomKETurbKEqCalculator.h
//! @brief     Kinetic equation of KE solver.
//! @author    WanYunbo, Bell, XuGang.


#pragma once

#include <vector>
#include "IncomScalarEqCalculator.h"
#include "CFDSolver.h"
#include "TK_Exit.h"

using namespace std;
namespace PHSPACE
{
class IncomKETurbKEqCalculator : public IncomScalarEqCalculator
{
public:
    IncomKETurbKEqCalculator();
    virtual ~IncomKETurbKEqCalculator();

    void AllocateGlobalVar(Grid *grid);
    void InitFlowAsRestart(Grid *grid);
    void IncompressibleInitial(Grid *grid);
    void SetDiffusionCoeff(Grid *grid);
    void UpdateProperties(Grid *grid);

    void CalcOtherMatrixACoeff(Grid *grid);
    void CalcOtherbCoeff(Grid *grid, int iEquation);
    
    void GetResidual(Grid *grid, vector<RDouble> &res);
    void InitialUnsteadyVar(Grid *grid);
    void UpdateUnsteadyVariable(Grid *grid);
    void UpdateBCValue(Grid *grid);

};

}
