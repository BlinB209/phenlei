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
//! @file      IncomKETurbEpsilonEqCalculator.h
//! @brief     Epsilon equation of KE solver.
//! @author    WanYunbo, Bell, XuGang.

#pragma once

#include <vector>
#include "IncomScalarEqCalculator.h"

using namespace std;

namespace PHSPACE
{
class IncomKETurbEpsilonEqCalculator : public IncomScalarEqCalculator
{
public:
    IncomKETurbEpsilonEqCalculator();
    virtual ~IncomKETurbEpsilonEqCalculator();

    void SetDiffusionCoeff(Grid *gridIn);
    void InitFlowAsRestart(Grid *gridIn);
    
    void solveScalarEquation(Grid *grid, int iEquation);
    
    void UpdateProperties(Grid *grid);
    
    void CalcOtherMatrixACoeff(Grid *grid);
    void CalcOtherbCoeff(Grid *grid, int iEquation);
    
    void GetResidual(Grid *grid, vector<RDouble> &res);
    void InitialUnsteadyVar(Grid *grid);
    void UpdateUnsteadyVariable(Grid *grid);
    void UpdateBCValue(Grid *grid);

public:
    void calWallBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void WallModification(Grid *grid);
    
    void CalcGenTerm(Grid *grid);
    void UpdateWallGen(Grid *grid);

public:
};

}
