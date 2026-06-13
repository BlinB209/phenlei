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
//! @file      IncomSATurbKEqCalculator.h
//! @brief     Kinetic equation of SA solver.
//! @author    WanYunbo, Bell, XuGang.
#pragma once

#include <vector>
#include "IncomScalarEqCalculator.h"
#include "TK_Exit.h"

using namespace std;

namespace PHSPACE
{
class IncomSATurbKEqCalculator : public IncomScalarEqCalculator
{
public:
    IncomSATurbKEqCalculator();
    virtual ~IncomSATurbKEqCalculator();

    void InitFlowAsRestart(Grid *gridIn);
    void AllocateGlobalVar(Grid *gridIn);
    void IncompressibleInitial(Grid *grid);
    void SetDiffusionCoeff(Grid *gridIn);
    void UpdateAfterIterloop(Grid *grid);    //! Updates the viscosity of the momentum equation.

    void calWallBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void calVelocityInletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void calPressureOutletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void CalcOtherMatrixACoeff(Grid *grid);
    void CalcOtherbCoeff(Grid *grid, int iEquation);
    
    void GetResidual(Grid *grid, vector<RDouble> &res);

    void InitialUnsteadyVar(Grid *grid);
    void UpdateUnsteadyVariable(Grid *grid);
    void UpdateBCValue(Grid *grid);

    void CalcGrad(Grid *grid);
};

}

