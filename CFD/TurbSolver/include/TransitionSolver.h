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
//! @file      TransitionSolver.h
//! @brief     transition solver.
//! @author    He Kun

#pragma once
#include "CFDSolver.h"
#include "Param_TransitionSolver.h"

namespace PHSPACE
{
class ActionKey;

void Transition_MxDQ(RDouble *mat,RDouble *dq,int n_transition,RDouble *df);

class TransitionSolver : public CFDSolver
{
private:

public:
    TransitionSolver();
    ~TransitionSolver();

public:
    void Solve();
    void Action(ActionKey *actkey);
    void FillActionKey(ActionKey *actkey, int action, int level);
    void Post();
    void PostSolve(Grid *gridIn, int stage, int level = 0);
    void DumpResultFile(Grid *grid, int level = 0);

    void UpdateResiduals(Grid *gridIn);
    virtual void ResetWallScalar(Grid *gridIn) {};

    virtual void InviscidFlux  (Grid *gridIn){};
    virtual void ViscousFlux   (Grid *gridIn){};
    virtual void DualTimeSource(Grid *gridIn){};

    virtual void ZeroResidualOfSpecialCells(Grid *gridIn){};

    void SourceFlux(Grid *gridIn);
    virtual void FreeGradientProxy(Grid *gridIn){};
    virtual void Diagonal(Grid *gridIn){};

    virtual void InitSpectrum(Grid *gridIn){};

    virtual void SourceFluxTwoEquation(Grid *gridIn){};

    virtual void InitCGrid(Grid *fineGrid, Grid *coarseGrid){};

    string CastName(const string &name);

    //! Get the file name for Residual dumping.
    virtual const string GetResidualFileName();
    //! Get the file name for Restart info dumping.
    virtual const string GetRestartFileName();

    int GetNumberOfEquations();

    //! Compare the outnstep in the flowfield file between the NS and turbulence.
    void CompareOutStepOfFlowfieldFile(int outnstepofNS, int outnstepofTransition) const;

    //! Get control paramters.
    LIB_EXPORT Param_TransitionSolver *GetControlParameters();
    
private:
    void RegisterCFDSolverInterfaceField();

    //! Register the interpoint field in the turbulent solver.
    void RegisterCFDSolverInterpointField();
    void RegisterOversetField();
    void DeAllocateOversetInterfaceVar(Data_ParamFieldSuite *dataStore);
    void InitCoarseGridsFlow ();

    //! Prepare overset interface data for communication.
    void PrepareOversetInterfaceData(Data_ParamFieldSuite *datastore, InterfaceDataProxy *interfaceDataProxy);

    void InitDependentVariables();
};

}