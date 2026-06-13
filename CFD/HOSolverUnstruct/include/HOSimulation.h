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
//! @file      HOSimulation.h
//! @brief     High Order Solver Simulation.
//! @author    xxx

#pragma once
#include "Region.h"
using namespace PHSPACE;
class Zone;
class TimeSpan;
class Region;

#include "HODefine.h"

namespace HOUnstruct
{
    
class HOSimulation
{
private:
    PHSPACE::Region * region;
    PHSPACE::TimeSpan * simu_ts;

    int isUnsteady;

public:
    HOSimulation(int dimension=3);
    ~HOSimulation();
public:
    void Start();
    void SolveSimu();
    void PostSimu();
    void MultiGridInitFlow();
    void InitializeProtoTable();
    void Run();

private:
    void InitGeometry();
    void ReadGeometry();
    void SwapGeometry();
    void InitOversetGrid() {}
    void OversetGridConfig();
    void ReadGrid();
    void ReadSimulationInfo();

    //! Define solver names for all zones and create the newly defined solvers.
    //! IMPORTANT for solver developers:
    //! Define your solver name which is same with the solver class name here!
    void DefineAndCreateSolvers();

    //! Initialize each solver on each zone.
    void InitializeSolvers();

    void ComputeWallDistance();

    void UpdateAllZonesPara();
    void RunPostProcessing();

    void SolveOneOuterStep();
    void UpdateOuterIterationStep(int & outnstep);
    void PreSolve();
    void MainSolve();
    void PostSolve();

    void InitGlobalBoundaryCondition();

    void SolveUnsteady();
    void SolveSteady();
    void SolveSteadyField();

    bool IsSubIterationConvergence(const int & innstep, const int & min_sub_iter, const int & max_sub_iter, const RDouble & tol_sub_iter, const RDouble & averageQuotientMaxInSolvers);
};

}