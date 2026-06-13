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
//! @file      IncomCalculator.h
//! @brief     SIMPLE class solver.
//! @author    WanYunbo, Bell, XuGang.

#pragma once

#include <map>
#include <iostream>
#include "IncomCalculator.h"
#include "CFDSolver.h"
#include "hypreCalculator.h"


namespace PHSPACE
{
    class UnstructGrid;

    class IncomCalculator
    {

    public:
        IncomCalculator() {};
        ~IncomCalculator() {};

    public:

        virtual void SetSolverIndex(int solverIndex) { this->solverIndex = solverIndex; }
        virtual int GetSolverIndex() { return this->solverIndex; }
        virtual int GetSolverIndex(int iEquation);

        virtual void calculateLinearEquation(Grid *grid, int iEquation,
            RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace);

        virtual void SetAxEqualbCalculator(LinearEquationSystemCalculator *hypreAxEqbSolverIn);

        virtual void CalcRes(Grid *grid, int iEquation, RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace);

        void SkewnessCorrection(Grid *grid, RDouble *phi, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz);

        void FaceGradCorrection(Grid *grid, RDouble *dphidxf, RDouble *dphidyf, RDouble *dphidzf, RDouble *phi, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz);

        void FaceValueCorretion(Grid *grid, RDouble *dphidxf, RDouble *dphidyf, RDouble *dphidzf, RDouble *phif, RDouble *phi);

        void GradCorrection(Grid *grid, RDouble *phif, RDouble *phi, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz);

        virtual void GradientCalculation(Grid *grid, string phiName, string dphidxName, string dphidyName, string dphidzName, bool isComm = false);
        
        void calcFaceWeight(Grid *grid);
    public:

        void CompressAnInterfaceVar(RDouble *field);
        void CommunicateAnInterfaceVar();
        void DecompressAnInterfaceVar(RDouble *field);
        void CommunicateAnInterfaceVar(RDouble *field);

        void CompressAnInterfaceVar(DataContainer *&dataContainer, RDouble *field, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex);
        void DecompressAnInterfaceVar(DataContainer *&dataContainer, RDouble *field, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex);

        void InitializeMatrixACoeff(Grid *gridIn);
        void InitializeBCoeff(Grid *gridIn);

    public:
        LinearEquationSystemCalculator *AxEqualbCalculator;

        int solverIndex;

        vector <PH_Request > requestContainer;
        vector <vector <DataContainer *> > receivedDataBuffer;
        vector <vector <DataContainer *> > sendDataBuffer;

    };
}
