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
//! @file      HOSolverUnstruct.h
//! @brief     High Order Solver.
//! @author    Zhang Jian, Wan yunbo, Xu gang.

#pragma once
#include <vector>

#include "CFDSolver.h"
#include "HODefine.h"
#include "HOGeometryStructure.h"
#include "HOStandardElement.h"

using namespace HOUnstruct;
//using namespace std;

namespace PHSPACE
{
class ActionKey;
class ActionTag;
class Grid;
class FaceProxy;

//! HighOrderSolverUnsruct solver
class HOSolverUnstruct : public CFDSolver
{
private:
    //[hlevel]
    std::vector< HighOrderGrid > highOrderGrids;

public:
    HOSolverUnstruct();
    ~HOSolverUnstruct();
    
public:
    //! These are cfd solver interface:

    //! Create and init ControlParameters
    LIB_EXPORT void InitControlParameters();

    //! Update some specific parameters, not in control parameter files.
    void ReadParameter();

    //! Allocate varibales and update to grid:
    void AllocateGlobalVar(Grid * gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);

    //! Register variables that need to MPI communication.
    void RegisterCFDSolverInterfaceField();

    //! Initialize flow field.
    //! For NS Solver as example, Will call some other virtual funcions:
    //! 1. InitFlowAsRestart
    //! 2. InitFlowAsReadingRestart
    //! 3. InitDependentVariables
    //void InitFlow();
    bool JudgeIfRestart();
    bool JudgeIfReadAverage();
    void InitFlowAsRestart();
    //void InitFlowAsRestart();

    void InitMemory();
    void ReleaseMemory();
    void AllocateGlobalVariables();
    void InitFlow();

    virtual void InitDependentVariables();
    virtual void Boundary(Grid *gridIn);
    void CompViscousCoef(Grid *gridIn);
    void InitViscousCoefTurb(Grid *gridIn);

    void CompGamaAndTField(Grid *gridIn);
    void ComputeGaussPointQ(Grid *gridIn);
    void ComputeQAverage(Grid *gridIn);

    //! Init coarse grids flow
    void InitCoarseGridsFlow();
    virtual void RestrictAllQ(Grid *fgridIn, Grid *cgridIn);

    void LocalTimeStep(Grid *gridIn);
    void GlobalTimeStep(Grid *gridIn);
    void LocalGlobalTimeStep(Grid *gridIn);
    RDouble ComputeCFL(int hLevel);
    void InvSpectrumRadius(Grid *grid, RDouble * invSpectrum);
    void VisSpectrumRadius(Grid *grid, RDouble * visSpectrum);


    //! Register action for communicate interface data.
    //! Usually, using ActionReflect function.
    //void CommunicationInterfaceData();
    void UploadInterfaceData(ActionKey *actkey);
    void DownloadInterfaceData(ActionKey *actkey);

    virtual void ActionReflect(ActionTag *acttag);

    //! Set resiudals to zero
    //! "res" is a pointer store in grid.
    void ZeroResiduals(Grid *gridIn);
    FieldProxy * GetResidualProxy(Grid *gridIn);

    virtual FieldProxy * CreateFieldProxy(Grid *gridIn);
    virtual void StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy);

    void LoadQ(Grid *gridIn, FieldProxy *qProxy);
    void Lhs(Grid *gridIn, FieldProxy *dqProxy, double coef);

    //! Compute local time step "dt".
    //! "dt" is a pointer stored in grid.
    void TimeStep(Grid *gridIn);
    void ComputeMinTimeStep(Grid *grid, RDouble & minDt, RDouble & maxDt);
    virtual void ReduceMaxTimeStep(Grid *gridIn, RDouble globalMinDt);

    //! Integrate one time step
    //! Including ComputeResidual and time integration(e.g. ForwardLUSGS).
    //void ForwardStep(Grid *grid, FieldProxy *rhsProxy);

    //! Compute residuals.
    //! For NS Solver as example, will call some other virutal functions:
    //! LoadResidual, UpdateResidual
    //! UpdateResidual includes Boundary, RightHandSide
    //! RightHandSide includes InviscidFlux, ViscousFlux, SourceFlux
    //void ComputeResidual(Grid *grid, FieldProxy *rhsProxy);

    virtual void LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy);
    virtual void UpdateResiduals(Grid *gridIn);
    virtual void RightHandSide(Grid *gridIn);

    //! Inviscid flux.
    virtual void InviscidFlux(Grid *gridIn);

    //! Viscous flux.
    virtual void ViscousFlux (Grid *gridIn);

    //! Source flux, such as chemical source term.
    virtual void SourceFlux   (Grid *gridIn);

    //! res = res/vol;
    void ResDivideVol   (Grid *gridIn);

    virtual void RungeKutta(Grid *grid, FieldProxy *rhsProxy);
    void InviscidFluxForFace(Grid *grid_in);
    //! Inviscid flux of cell
    void InviscidFluxForCell(Grid *grid_in);

    //! Time integration using lusgs scheme.
    void ForwardLUSGS(Grid *gridIn);

    //! sweeps == 1.
    virtual void FillField(Grid *gridIn, FieldProxy *field1_proxy, FieldProxy *field2_proxy);
    virtual void SolveLUSGS(Grid *gridIn, FieldProxy *dq_proxy);

    //! sweeps > 1.
    virtual void FillField(Grid *gridIn, FieldProxy *field_proxy, RDouble value);
    virtual void SolveLUSGS(Grid *gridIn, FieldProxy *dq_proxy, FieldProxy *rhsProxy, int sweeps, double epsilon);

    virtual void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dq_proxy);

    virtual void RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy);

    //! After one time step inner integration, do some post solve work
    //! e.g CommunicationInterfaceData, dump restart file, dump aerodynamic coefficient et, al.
    virtual int GetNPostSolve();
    void PostSolve(Grid *gridIn, int stage, int level);
    int GetSchemeID(const string &scheme_name);

    virtual void GetResidual(ActionKey *actkey);

    void Post() {};
    //! Prepare post visualization variables, store it in to  postVisualization database.
    LIB_EXPORT void ComputePostVisualVariables(Post_Visual * postVisualization);

public:
    //! These are HOSolverUnstruct defined interface.

};

}