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
//! @file      NSSolver.h
//! @brief     The base NS solver, which will derive both structured and unstructured NS solver,
//!            these two solvers are actually used instead of this NSSolver.
//! @author    He Xin, Bell, He Kun, He Xianyao, Wan Yunbo, Liu Jian, Xu Qingxin, Zhang Zipei, Li peng, Ma Yankai, 
//!            Min Yaobing, Meng Liyuan, Xu Gang, Guo Yongyan, Zhang Yang.

#pragma once
#include "CFDSolver.h"
#include "Flux_Inviscid.h"

namespace PHSPACE
{
class ActionKey;
class Grid;
class FaceProxy;
class Param_NSSolver;

//! @brief The base class of all other NS solvers.
//!        The inherit order is: NSSolver->NSSolverStructFV/NSSolverStructFD/NSSolverUnstruct.\n
//!        Note that, the finally actually used solver are NSSolverStructFV/NSSolverStructFD/NSSolverUnstruct,
//!        rather than the NSSolver. This solver contains the common components and functions which are 
//!        independent on the grid type is structured or unstructured.
class NSSolver : public CFDSolver
{
public:
    NSSolver();
    ~NSSolver();

public:
    //! Return the control parameters of the NS solver.
    //! It is important that, each level of inherit class has their own control parameter.
    LIB_EXPORT Param_NSSolver *GetControlParameters() const;

    //! To get the number of equations to be solve in the solver.
    //! which is the instantiation of the one in Solver.
    int GetNumberOfEquations();    
    
    //! Compute viscous coefficient only for the perfect gas.
    //! Here it is a pure virtual function, it will be re-defined in NSSolverStructFV/NSSolverStructFD/NSSolverUnstruct.
    virtual void ComputeViscousCoefficientWithoutChemical(Grid *grid)  = 0;
    virtual void ComputeViscousCoefficient(Grid *grid) {};
    //! Compute specific heat ratio 'gama' and temperature.
    //! Here it is a pure virtual function, it will be re-defined in NSSolverStructFV/NSSolverStructFD/NSSolverUnstruct.
    virtual void ComputeGamaAndTemperature(Grid *grid) = 0;

    //! Compute time step of steady simulation. Using Local Time Step method default.
    void TimeStep(Grid *grid);

    //! Compute time step using local method, meaning that the time step of each cell element
    //! depends on it's volume.
    //! Here it is a pure virtual function, it will be re-defined in NSSolverStructFV/NSSolverStructFD/NSSolverUnstruct.
    virtual void LocalTimeStep(Grid *grid_in) = 0;

    //! Compute time step using global method, meaning that using the minimum time step of cells.
    //! Here it is a pure virtual function, it will be re-defined in NSSolverStructFV/NSSolverStructFD/NSSolverUnstruct.
    virtual void GlobalTimeStep(Grid *grid_in) = 0;

    virtual void LocalGlobalTimeStep(Grid *grid_in) = 0;
    
    //! Compute the residual of the Right-Hand-Side of the equations.
    //! The residual is the flux actually, including inviscid flux, viscous flux,
    //! source flux, etc. These fluxes are defined in NSSolverStructFV/NSSolverStructFD/NSSolverUnstruct.
    //! This function is an very important main function of the NS solver.
    void UpdateResiduals(Grid *grid);

    //! =======Some special post-process: END==============

    //! Return the Residual file name for dumping.
    virtual const string GetResidualFileName();

    //! Do some post-process work after simulation over.
    void Post();

    //! Do some post-process work after each iteration over.
    void PostSolve(Grid *grid, int stage, int level = 0);

    void DumpResultFile(Grid *grid, int level = 0);

    //! switch the compute method used for self-adaption calculation.
    void CheckResult(Grid *grid, int level = 0);

    void CheckResiduals(Grid *grid, int level = 0);

    //! dump CFL number.
    void DumpCFLNumber(Grid *grid, int level = 0);

    //! exam the surface heating change.
    void CheckSurfaceHeatingChange(Grid *grid, int level = 0);

    //! Fill the action.
    //! @param[in ]: action     The task type.
    //! @param[in ]: level      The multi-grid level.
    void FillActionKey(ActionKey *actkey, int action, int level);

    //! The following three functions are used in LUSGS, which will be consider more latter.
    virtual void GetResidual(ActionKey *actkey);
    
    //! The following three functions are used in overset grid, consider more latter.
    void DeAllocateOversetInterfaceVar(Data_ParamFieldSuite *dataStore);    
    void PrepareOversetInterfaceData(Data_ParamFieldSuite *dataStore, InterfaceDataProxy *interfaceDataProxy);

    virtual void ZeroResidualOfSpecialCells(Grid * gridIn) {};

protected:
    //! Compute CFL number at each iteration.
    RDouble ComputeCFL(RDouble partialCFL = -1.0);

    //! Source flux, such as chemical source, unsteady source, volume force, etc.
    void SourceFlux(Grid *grid);

    //! Source flux: Compute the chemical source terms.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    virtual void ChemicalSource(Grid *grid) = 0;

    //! Source flux: compute the rotating source term.
    virtual void RotatingSource(Grid *grid) = 0;

    //! Source flux: Compute the gravity source terms.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    virtual void gravitySource(Grid *grid) = 0;

    //! Source flux: Compute the porous medium source terms.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    virtual void PorousMediumSource(Grid *grid) = 0;

    //! Source flux: compute the unsteady dual time source.
    virtual void DualTimeSource(Grid *grid) = 0;

    virtual void Turb_Sengy(Grid *grid) {};

    //! Return if need to dump out the field data.
    bool WantVisualField(Grid *grid);

    //! Register interface variables into field of NS Solver.
    //! The common part of structured and unstructured solver.
    virtual void RegisterCFDSolverInterfaceField();

    virtual void ReleaseCFDSolverInterfaceField();

    //! Register interface variables into field of NS Solver.
    //! The unstructured part.
    virtual void RegisterCFDSolverInterfaceFieldUnstruct();

    void RegisterOversetField();

    //! Register the interpoint field in the NS solver.
    virtual void RegisterCFDSolverInterpointField();

    //! Dump out typical standard case: turbulent flat plate.
    virtual void Turbulence_Flat_Plate_Output_NASA(Grid *gridIn) {};

private:
    //! Initialize the flow on the coarse grid if using multigrid.
    void InitCoarseGridsFlow();

    //! Initialize some dependent variables, such as gama, viscous coefficients, temperature,
    //! before each iteration.
    void InitDependentVariables();
};

void SetInviscidSchemeParameters(InviscidSchemeParameter *invSchemePara, FaceProxy *face_proxy, Param_NSSolver *parameters);

}
