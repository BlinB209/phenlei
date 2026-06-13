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
//! @file      LESSolver.h
//! @brief     Large-eddy simulation solver.
//! @author    Zhang Zipei.

#pragma once
#include "CFDSolver.h"
#include "Param_LESSolver.h"

namespace PHSPACE
{
class ActionKey;
class ActionTag;
class Gradient;
class Grid;

class LESSolver : public CFDSolver
{
private:
    static bool init_param;
    static bool free_param;

    static int nTurbEquations;

public:
    LESSolver();
    ~LESSolver();

public:
    void Action(ActionKey *actkey);
    void FillActionKey(ActionKey *actkey, int action, int level);
    void Post();
    void PostSolve(Grid *grid_in, int stage, int level = 0);
    void DumpResultFile(Grid *grid, int level = 0);

    static void InitParameter();
    static void FreeParameter();

    virtual void ComputeViscousCoeff(Grid *grid) {};

    virtual void ObtainViscosity(Grid *grid) {};

    //int GetNumberOfEquations();

    //! Get the file name for Restart info dumping.
    virtual const string GetRestartFileName();

    //! Compare the outnstep in the flowfield file between the NS and LES.
    void CompareOutStepOfFlowfieldFile(int outnstepofNS, int outnstepofTurb) const;

    LIB_EXPORT Param_LESSolver * GetControlParameters();

private:
    void RegisterCFDSolverInterfaceField();

    void InitDependentVariables();
    
};

}