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
//! @file      DemoSolver.h
//! @brief     Demo solver.
//! @author    Meng Liyuan.
#pragma once
#include "CFDSolver.h"
namespace PHSPACE
{
class Grid;
class Param_DemoSolver;
class DemoSolver : public CFDSolver
{
public:
    DemoSolver();
    ~DemoSolver();
public:
    LIB_EXPORT void InitControlParameters();
    LIB_EXPORT Param_DemoSolver *GetControlParameters() const;
    void AllocateGlobalVar(Grid * grid);
    void DeAllocateGlobalVar(Grid *grid);
    void ZeroResiduals(Grid * grid);
    void TimeStep(Grid * grid);
    void ComputeResidual(Grid *grid, FieldProxy *rhs_proxy);
    void SourceFlux(Grid *grid);
    void InviscidFlux(Grid *grid);
    void ViscousFlux(Grid *grid);
    void RungeKutta(Grid *grid, FieldProxy *rhs_proxy);
    void Post();
    void PostSolve(Grid * grid, int stage, int level);

private:
};

}