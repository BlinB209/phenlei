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
//! @file      IncomUnstructSolver.h
//! @brief     Species solver.
//! @author    Wan Yunbo, Bell, Xu Gang.

#pragma once

#include <map>
#include <iostream>

#include "IncomSpeciesEqCalculator.h"
#include "IncomEnergyEqCalculator.h"
#include "IncomKETurbKEqCalculator.h"
#include "IncomKETurbEpsilonEqCalculator.h"
#include "IncomSATurbKEqCalculator.h"
#include "MomEqCalculator.h"
#include "PPEqCalculator.h"
#include "hypreCalculator.h"
#include "Param_INCompSolverUnstruct.h"

namespace PHSPACE
{

    class IncomUnstructSolver : public CFDSolver
    {
    public:
        IncomUnstructSolver();
        ~IncomUnstructSolver();

    public:
        void DumpResultFile(Grid *grid, int level = 0);
        void DumpResidual(Grid *grid);
        
        void InitFlowAsReadingRestart(ActionKey *actkey);
        void ReadRestartH5(ActionKey *actkey);
        void DumpRestartH5(ActionKey *actkey);
        void CreateH5RestartFile(ActionKey *actkey);

        void AirForceCoef(ActionKey *actkey);
        void AirForceCoefParts(ActionKey *actkey, int partID, int Coordinate = GlobalCoordinate);
        void CpDistriCoef(ActionKey *actkey);

        void IncompressibleInitial(Grid *grid);
        void UpdateUnsteadyProperties(Grid *grid);
        void UpdateUnsteadyFlux(Grid *grid);
        void UpdateUnsteadyVariable(Grid *grid);
        void GetResidual(Grid *grid, vector<RDouble> &res);
        void GetPhiName(vector<string> &phiName);
        void SolveIncomSteadyField();
        void SolveTurbulenceEquation(Grid *gridIn);
        void SolveEnergyEquation(Grid *gridIn);
        void SolveSpeciesEquation(Grid *gridIn);

        void SolutionIsConverged(Grid *grid, bool *flag);
        virtual void AllocateGlobalVar(Grid *gridIn);
        void InitControlParameters();
        LIB_EXPORT Param_INCompSolverUnstruct * GetControlParameters() const;
        virtual void ComputePostVisualVariables(Post_Visual *postVisualization);
        
        void CommunicationInterfaceData() {};

        void InitFlowAsRestart();
        bool JudgeIfRestart();
        void InitDependentVariables();
        void ConstructLocalCell2GlobalCellMap(Grid *grid);
        void SetAxEqbSolver (Grid *grid);

        void UpdateUnsteadyFlow(Grid *grid);

    public:

        MomEqCalculator *momEqCalculator;
        PPEqCalculator *PPEquaCalculator;
        IncomScalarEqCalculator *IncomTurbKSolver;
        IncomScalarEqCalculator *IncomTurbEplisionSolver;
        IncomScalarEqCalculator *IncomTurbSASolver;
        IncomScalarEqCalculator *IncomEnergySolver;
        IncomScalarEqCalculator **IncomSpeciesSolver;

        LinearEquationSystemCalculator *AxEqualbCalculator;

    public:
        virtual void Post() {}
    };

}
