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
//! @file      IncomScalarEqCalculator.h
//! @brief     scalar solver.
//! @author    WanYunbo, Bell, XuGang.
#pragma once

#include "IO_FileName.h"
#include "IncomCalculator.h"

#include "TK_Exit.h"


namespace PHSPACE
{
class IncomScalarEqCalculator : public IncomCalculator
{
public:

public:
    IncomScalarEqCalculator();
    ~IncomScalarEqCalculator();

public:
    virtual void GetResidual(Grid *grid, vector<RDouble>& res);

    virtual void InitFlowAsRestart(Grid *gridIn);
    virtual void AllocateGlobalVar(Grid *gridIn);
    virtual void InitialUnsteadyVar(Grid *grid) {};
    virtual void IncompressibleInitial(Grid *grid);
    virtual void UpdateUnsteadyVariable(Grid *grid);

public:
    virtual void solveScalarEquation(Grid *grid, int iEquation);
    virtual void constructMatrixACoeff(Grid *grid, int iEquation);
    virtual void constructBCoeff(Grid *grid, int iEquation);
    
    virtual void UpdateBCValue(Grid *grid);
    virtual void UpdateProperties(Grid *grid) {};
    virtual void UpdateAfterIterloop(Grid *grid) {};
    virtual void UpdateUnsteadyProperties(Grid *grid) {};
    virtual void UpdateUnsteadyFlux(Grid *grid) {};

    void Post() {};

public:
    virtual void AllocateUnsteadyVar(Grid *grid);
    virtual void CalcGrad(Grid *grid);
    virtual void SetDiffusionCoeff(Grid *grid);
    virtual void CalcOtherMatrixACoeff(Grid *grid) {};
    virtual void CalcOtherbCoeff(Grid *grid, int iEquation) {};
    
public:

    void treatBC(Grid *grid, int iVariable, RDouble matrixCoeff);

    virtual void calWallBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);
    virtual void calSymmetryBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    virtual void calVelocityInletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    virtual void calPressureOutletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    virtual void calInterfaceBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);
    
    virtual void calPressureInletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    virtual void calFarfieldBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff);

    void calcInvisSourceTerm(Grid *grid, int solverIndex);

    void FirstUpwindInvisTerm(Grid *grid);

    void SecondUpWindSourceTerm(Grid *grid, int solverIndex);

    void QUICKSourceTerm(Grid *grid, int solverIndex);

    void SecondCentralSourceTerm(Grid *grid, int solverIndex);

    void VisMatrixTerm(Grid *grid, int solverIndex);

    void VisSourceTerm(Grid *grid, int solverIndex);

    void calcTransMatrixTerm(Grid *grid);

    void TransMatrixTerm_1st(Grid *grid);

    void TransMatrixTerm_2nd(Grid *grid);

    void calcTransSourceTerm(Grid *grid, int iVariable);

    void TransSourceTerm_1st(Grid *grid, int iVariable);

    void TransSourceTerm_2nd(Grid *grid, int iVariable);

    void relaxMatrixTerm(Grid *grid, int iVariable);

    void relaxSourceTerm(Grid *grid, int iVariable);

    void ImplicitEuler_AllocateMemory(Grid *grid, vector<string> &phiName);
    void ImplicitEuler_ReInitTimeVar(Grid *grid, vector<string> &phiName);
    void ImplicitEuler_SaveOldTimeValue(Grid *grid, vector<string> &phiName);
    void ImplicitEuler_InitRestartTimeVar(Grid *grid, vector<string> &nameLists);
    void ImplicitEuler_ReadRestartFile(Grid *grid, vector<string> &nameLists, hid_t grploc);

    void Implicit2ndOrder_AllocateMemory(Grid *grid, vector<string> &phiName);
    void Implicit2ndOrder_ReInitTimeVar(Grid *grid, vector<string> &phiName);
    void Implicit2ndOrder_SaveOldTimeValue(Grid *grid, vector<string> &phiName);
    void Implicit2ndOrder_InitRestartTimeVar(Grid *grid, vector<string> &nameLists);
    void Implicit2ndOrder_ReadRestartFile(Grid *grid, vector<string> &nameLists, hid_t grploc);
};

}
