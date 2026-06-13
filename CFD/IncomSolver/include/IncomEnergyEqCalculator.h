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
//! @file      IncomEnergyEqCalculator.h
//! @brief     energy equation solver.
//! @author    WanYunbo, Bell, XuGang.

#pragma once

#include "IncomScalarEqCalculator.h"
#include "IncomGas.h"

namespace PHSPACE
{
class IncomEnergyEqCalculator : public IncomScalarEqCalculator
{
public:
    IncomEnergyEqCalculator();
    ~IncomEnergyEqCalculator();

public:

    void InitFlowAsRestart(Grid *gridIn);
    void AllocateGlobalVar(Grid *gridIn);
    void IncompressibleInitial(Grid *grid);
    void UpdateProperties(Grid *grid);

    void CalcGrad(Grid *grid);
    void CalcOtherbCoeff(Grid *grid, int iEquation);

    void calWallBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void calVelocityInletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void calPressureInletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);
    
    void calPressureOutletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);
    
    void calFarfieldBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void calInterfaceBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void GetResidual(Grid *grid, vector<RDouble> &res);
    void InitialUnsteadyVar(Grid *grid);
    void UpdateUnsteadyVariable(Grid *grid);
    void UpdateBCValue(Grid *grid);
};

}