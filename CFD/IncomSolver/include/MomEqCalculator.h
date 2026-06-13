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
//! @file      MomEqCalculator.h
//! @brief     mom equation solver.
//! @author    WanYunbo, Bell.
#pragma once

#include <map>
#include <iostream>
#include "IncomGas.h"
#include "IncomScalarEqCalculator.h"

namespace PHSPACE
{
    class UnstructGrid;
    class MomEqCalculator : public IncomScalarEqCalculator
    {

    public:
        MomEqCalculator();
        ~MomEqCalculator() {};

    public:

        void CalcGrad(Grid *grid);
        void solveMomentumEquations(Grid *grid);
        void solveScalarEquation(Grid *grid, int iEquation);
        
        void CalcVelocityRes(UnstructGrid *grid);

        void CoefPrepForPres(Grid *grid);
        void CorrectFaceFlux(Grid *grid);
        void calcFaceFlux(Grid *grid);
        void updateBCValue(Grid *grid);
        void pressureSourceTerm(Grid *grid, int iVariable);
        void BulkVisSourceTerm(Grid *grid, int iVariable);
        void GravitySourceTerm(Grid *grid, int iVariable);

        virtual void CalcOtherbCoeff(Grid *gridIn, int iVariable);
        
        void SolveExplicitVelocityEquations(Grid *grid);
        void UpdateUnsteadyFlux(Grid *grid);
        void IncompressibleInitial(Grid *grid);
        void CalcWalldistOnWall(Grid *grid);
        void UpdateBoundaryPressureAndVelocity(Grid *grid);
        void InitFlowAsRestart(Grid *grid);
        void AllocateGlobalVar(Grid *grid);
        
        void AllocateUnsteadyVar(Grid *grid);
        void InitialUnsteadyVar(Grid *grid);
        void UpdateUnsteadyProperties(Grid *grid);
        void UpdateProperties(Grid *grid);
        
        void GetResidual(Grid *grid, vector<RDouble> &res);
        void UpdateUnsteadyVariable(Grid *grid);
    };

    void FlowFarField_GetBoundrayData(Grid *grid, Data_Param *bcData);

}
