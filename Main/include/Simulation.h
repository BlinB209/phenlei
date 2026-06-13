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
//! @file      Simulation.h
//! @brief     Explain this file briefly.
//! @author    PHengLEI team.

#pragma once
#include "Precision.h"

namespace PHSPACE
{
class Zone;
class TimeSpan;
class Region;

class Simulation
{
private:
    Region *region;
public:
    Simulation(int dimension);
    ~Simulation();
public:
    //! Select task according to taskType.
    void Start();

    //! Main part of simulation.
    void SolveSimulation();

    //! Post process of simulation.
    void PostSimulation();

private:
    //! Main work flow of CFD simulation process.
    void RunCFD();

    //! Compute wall distance for turbulent flow while there was no .wdt file.
    void ComputeWallDistance();

    //! Need convert process or not.
    void IsNeedConvert();

    //! The work flow to convert, partition and compute.
    void RunIntegrativeSolver();

    //! Run SpecDiffHyb solver.
#ifdef USE_SpecDiffHybSolver
    void RunSpecDiffHybSolver();
    void RunSpecSolver();
#endif

    //! Identify the motive derivatives from the unsteady result file "./results/aircoef.dat".

    //! Run post process without flow field simulation.
    void RunPostProcessing();

    //! Read the grid such as grid coordinate, boundary information, wall distance for turbulent flow and interpoint information.
    void ReadGeometry();

    //! Update all zones parameter.
    void UpdateAllZonesPara();

    //! Initialize the geometry.
    void InitGeometry();

    //! Communicate the grid metrics on the interfaces.
    void SwapGeometry();

    //! Initialize overset grid.
    void InitOversetGrid();

    //! Run CFD using the self-adaption method.
    void SelfAdaptionRunCFD();
    void FirstStepRunCFD(bool isStart = true);
    void SecondStepRunCFD(bool isStart = false);
    void ThirdStepRunCFD(bool isStart = false);
    void FinalStepRunCFD();

    //! The numerical iteration process similar to the function named SolveSimulation.
    void SelfAdaptionSolve();

    //! Release the dynamic memories related to the solvers.
    void RefreshSolvers();

    //! Check whether the current computation is restart process.
    bool IsRestartProcess();

    //! Define solver names for all zones and create the newly defined solvers.
    //! IMPORTANT for solver developers:
    //! Define your solver name which is same with the solver class name here!
    void DefineAndCreateSolvers();

    //! Initialize global boundary condition.
    void InitGlobalBoundaryCondition();

    void InitGlobalVolumeCondition();

    //! Initialize each solver on each zone.
    void InitializeSolvers();

    //! Initionalize flow field in the first flowInitStep steps with coarse grid when multigrid simulation.
    void MultiGridInitFlow();

    //!
    void OversetGridConfig();

    //!
    void UpdateOuterIterationStep(int &outerIterStep);

    //! 
    void SolveOneOuterStep();

    //!
    void PreSolve();

    //!
    void MainSolve();

    //!
    void PostSolve();

    //!
    void SolveSteady();

    //!
    void SolveSteadyByGMRES();

    //!
    void SolveUnsteady();

    //!
    bool IsSubIterationConvergence(const int &innerIterStep, const int &minSubIteration, const int &maxSubIteration, const RDouble &tolalSubIteration, const RDouble &averageQuotientMaxInSolvers);

};
}

