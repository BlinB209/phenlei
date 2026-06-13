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
//! @file      Controller.h
//! @brief     the class 'Controller' is used to assemble grid zones and \n
//!            solvers together to simulate. All simulations are controlled \n
//!            by the class 'Controller'.
//! @author    Bell.

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
#include <vector>
namespace PHSPACE
{
class Zone;
class Region;
class FieldProxy;
class CFDSolver;
class ZoneFieldProxy;

//! @brief  It defines the class 'Controller', which is used to assemble grid zones and
//!         solvers together to simulate. All simulations are controlled by this class.
//!         All solvers are implemented in the way of multi-grid cycles, which is defined
//!         in Controller.
class Controller
{
private:
    //! Computational region, in which the local computational grids are defined.
    //! this is the very important base of the CFD solver.
    Region *region;

    //! Flow field variables proxy, actually, it stores the fields on different grid levels in MG.
    ZoneFieldProxy *qProxy;

    //! Right Hand Side (RHS) proxy, actually, it stores the RHS on different grid levels in MG.
    ZoneFieldProxy *rhsProxy;

    //! Solver index.
    int solverID;

public:
    Controller(int iSolver, Region *region);
    ~Controller();

public:
    LIB_EXPORT void Initialize(int iSolver);
    LIB_EXPORT void InitCoarseGridsFlow(int iSolver, int level);

    LIB_EXPORT void CleanUp(int iSolver);
    
    //! Steady simulation by implementing multi-grid cycle.
    //! Each solver simulation will be controlled by this multi-grid cycle, and this is the right
    //! reason why this class is named as 'Controller'.
    //! The standard simulation order is:
    //! Zero residual --> go into Multi-Grid cycle, if single level grid used, go to Relaxation directly.
    LIB_EXPORT void SolveSteadyField();

    //! Steady simulation by implementing multi-grid cycle.
    //! Each solver simulation will be controlled by this multi-grid cycle, and this is the right
    //! reason why this class is named as 'Controller'.
    //! The standard simulation order is:
    //! Zero residual --> go into Multi-Grid cycle, if single level grid used, go to Relaxation directly.
    LIB_EXPORT void SolveIncomSteadyField(int iSolver, int level = 0);

    //! Steady simulation by implementing multi-grid cycle.
    //! Each solver simulation will be controlled by this multi-grid cycle, and this is the right
    //! reason why this class is named as 'Controller'.
    //! The standard simulation order is:
    //! Zero residual --> go into Multi-Grid cycle, if single level grid used, go to Relaxation directly.
    LIB_EXPORT void SolveHyrbridDBPBSteadyField(int iSolver, int level = 0);

    //! Unsteady simulation using Dual Time Step method.
    //! In which, the steady simulation process is used as the virtual dual time step (or sub-iteration),
    //! and the steady relaxation sub-iterations are nested inside the real physical time step.
    LIB_EXPORT void SolveUnsteadyField();

    //! Update unsteady flow.
    //! The the flow field of n+1 physical time step is updated by using the last sub-iteration results.
    void UpdateUnsteadyFlow(int iSolver, int level = 0);
    //! Multi-Grid (MG) acceleration cycle control.
    //! All MG stages are defined in it, such as relaxation, restriction, prolongation, correction.
    void MultiGrid(int iSolver, int level = 0);

    //! Implicit and explicit time integration.
    //! Time step is first computed, and then using Runge-Kutta or LUSGS or other methods to advancing one step.
    void Relaxation(int iSolver, int level);

    //! Computing time step of each iteration, in which CFL number is used.
    //! Both local and global methods are defined. Computing process:\n
    //!   1. Compute time step of each zone.\n
    //!   2. Compute the global minimum time step.\n
    //!   3. Balance time step by reducing time step to no more than k-times of global minimum time step.
    void TimeStep(int iSolver, int level);

    //! Compute the inviscous/viscous/chemical spectrum radius, \n
    //! the results stored in array invSpectralRadius and visSpectralRadius and srs.
    void SpectrumRadius(int iSolver, int level);

    //! Compute the diagonal 'matrix' for LU-SGS. Although it is named as 'matrix',
    //! only scale variable is used since the Jacobian matrix is simplified to be spectrum radius.
    void Diagonal(int iSolver, int level);

    //! Compute and return the global minimum time step.
    //! For parallel computing, all_reduce is used to communicate the global minimum value.
    RDouble ComputeMinTimeStep(int iSolver, int level);

    //! Balance time step by reducing time step to no more than k-times of global minimum time step.
    void ReduceMaxTimeStep(int iSolver, int level);

    //! Set residual to be zero.
    void ZeroResiduals(int iSolver, int level = 0);

    //! Compute solver residual, by computing inviscid and viscous flux and source terms.
    void LoadResiduals(int iSolver, int level);

    //! Compute inviscid and viscous flux and source terms.
    void UpdateResiduals(int iSolver, int level);

    void UpdateResidualsOnlyforMixGrid(int iSolver, int level);

    //! Compute pressure factor for entropy fix (Method 2).
    void PressureFactor(int iSolver, int level);

    //! Compute the gama and temperature of flowfield on the level-th multigrid.
    void ComputeGamaAndTField(int iSolver, int level);

    //! Compute the independent flowfield value (such as gama/temperature/laminar viscous/gradient) on the level-th multigrid.
    void UpdateIndependentField(int iSolver, int level);

    //! Set the boundary condition of flowfield on the level-th multigrid.
    void SetBoundaryValueOfField(int iSolver, int level);

    void MixingPlaneBoundary(int iSolver, int level);

    void PreSolveforStructHighOrder(int iSolver, int level);

    void GetIndependentVariablesforStructHighOrder(int iSolver, int level);

    void GetDependentVariablesforStructHighOrder(int iSolver, int level);

    //! Compute the gama and temperature of flowfield on the level-th multigrid.
    void ComputeViscousCoeff(int iSolver, int level);

    //! Compute the cross and blending function F1 and F2 term of turbulence flow.
    void ComputeBlendFunctionSST(int iSolver, int level);

    //! Compute the precondition coefficient of low mach flow.
    void ComputePreconditionCoefficient(int iSolver, int level);

    //! Compute gradient for viscous flux computation.
    void ComputeGradient(int iSolver, int level);

    void ComputeLimiter(int iSolver, int level);

    //! Compute node bctype.
    void ComputeNodeBCType();
    void CommunicationNodeBCType();

    //! Relaxation by explicit Runge-Kutta method.
    void RungeKutta(int iSolver, int level);

    //! Define the post-processing for each solver.
    void PostSolve(int level = 0);

    //! Flow field interpolation to fine grid from computing results on coarse grid.
    LIB_EXPORT void InterpolatFineGrid(int iSolver, int level);

    //! Flow initialization using coarse level if Multi-Grid is used.
    LIB_EXPORT void MultiGridInitFlow(int level);

    //! Compute the sub-iteration norm.
    RDouble ComputeSubIterationNorm(int iSolver, int level);

private:
    //! Get the current being used solver index.
    int  GetSolverIndex();

    //! Get the number of post-processing which is need for each solver.
    int  GetNPostSolve(int iSolver);

    //! Get flow field on multi-grid of 'level' of iSolver-th solver.
    FieldProxy **GetQProxy  (int level);

    //! Get Right Hand Side on multi-grid of 'level' of iSolver-th solver.
    FieldProxy **GetRHSProxy(int level);

    //!
    Zone * GetZone(int iZone);

    CFDSolver * GetCFDSolver(int iZone, int iSolver);

#ifdef USE_GMRESSOLVER
public:
    //! GMRES solver -- a linear system solver
    //! This method can be referenced "Sen Zhang, Boqian Wang, et al. Implementation of a Newton-Krylov Algorithm in the Open-source Solver PHengLEI[C]."
    //! GPPS Hong Kong , October 17-19, 2023.
    void GMRESSolver(int iSolver, int level);
    void GMRESSolver_Coupled(int level);  // GMRESCoupled
    void GMRESSolver_backupversion(int iSolver, int level);
    void GMRESForZone(RDouble **dRdq, std::vector<int>& AI, std::vector<int>& AJ, RDouble **res, FieldProxy *dqProxy, const int& nTotalcells, const int& nEquations);  //GMRESSparse GMRESCSR
#endif

private:
    //! For New LUSGS with one sweep, Sub-routines of LUSGS method, which is divided into 4 steps.
    void LUSGS(int iSolver, int level);
    void DualTimeStepsLUSGS(int iSolver, int level);
    void LUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ);
    void BLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ);
    void MatrixLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ);
    void ForwardSweep(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    void BackwardSweep(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    void UpdateFlowField(int iSolver, int level, FieldProxy **dqProxy);
    void DetermineCFLNumber(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **dqProxy1);

    void CommunicationInterfaceData(int iSolver);
    void CommunicationInterfaceDataOneLevel(int iSolver, int level);
    void CommunicationDQ(int iSolver, int level, FieldProxy **dqProxy);
    void CommunicateGenericArray(int iSolver, int level);
    void CommunicateGradientArray(int iSolver, int level);
    void CommunicateGradientandLimitOnlyForMixGrid();
    void CommunicateQlQrOnlyforMixGrid();
    void CommunicateFacefluxOnlyforMixGrid();

    void LHS(int iSolver, int level, FieldProxy **resProxy, const RDouble &coef);
    void FillField(int iSolver, int level, FieldProxy **targetProxy, FieldProxy **sourceProxy);

    //! To start the block LU-SGS Method.
    //! @param[in ]: iSolver is the serial number of the current solver.
    //! @param[in ]: level denotes the serial number of the layer where the variables are compute in the multigrid method.
    //void BLUSGS(int iSolver, int level);

    //! Block Lower-Upper Symmetric Gauss-Seidel Scheme for unstruct grid.
    void BLUSGS(int iSolver, int level);

    //! Forward sweep in BLUSGS methods.
    void ForwardSweepBLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal , int &iSweep);
    
    //! Backward sweep in BLUSGS methods.
    void BackwardSweepBLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal);

    //! Matix Lower-Upper Symmetric Gauss-Seidel Scheme for struct grid.
    void MatrixLUSGS(int iSolver, int level);

    //! Forward sweep in MatrixLUSGS methods.
    void ForwardSweepMatrixLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal , int &iSweep);
    
    //! Backward sweep in MatrixLUSGS methods.
    void BackwardSweepMatrixLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal);

    //! Line implicit Lower-Upper Symmetric Gauss-Seidel Scheme for unstruct grid.
    void LineLUSGS(int iSolver, int level);

    void LineLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ);

    //! Forward sweep in line BLUSGS methods.
    void ForwardSweepLineLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal, int &iSweep);

    //! Backward sweep in line LUSGS methods.
    void BackwardSweepLineLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal);


    //! To execute the forward sweep of the block LU-SGS method.
    //! @param[in ]: iSolver is the serial number of the current solver.
    //! @param[in ]: level denotes the serial number of the layer where the variables are compute in the multigrid method.
    //! @param[out]: dqProxy is an array that stores the difference values between the conservative variables at n and n+1 time step.
    //void ForwardSweepBLUSGS(int iSolver, int level, FieldProxy **dqProxy);

    //! To execute the backward sweep of the block LU-SGS method.
    //! @param[in ]: iSolver is the serial number of the current solver.
    //! @param[in ]: level denotes the serial number of the layer where the variables are compute in the multigrid method.
    //! @param[out]: dqProxy is an array that stores the difference values between the conservative variables at n and n+1 time step.
    //void BackwardSweepBLUSGS(int iSolver, int level, FieldProxy **dqProxy);

    void RungeKuttaThirdTVD(int iSolver, int level);
    void CreatqProxyfor3rdRK(int iSolver, int level, FieldProxy **q0Proxy);
    void Stagefor3rdRK(int iSolver, int level, FieldProxy **q0Proxy, int nstage);
    
    //! MG stages.
    void InitResidual(int iSolver, int level);
    void LoadQ(int iSolver, int level);
    void RestrictAllQ(int iSolver, int level);
    void RestrictDefect(int iSolver, int level);
    void PutCorrection(int iSolver, int level);
    void PutCorrectionBack(int iSolver, int level);
    void CorrectFineGrid(int iSolver, int level);

    void CreateQProxy   (int iSolver, int level);
    void DestroyQProxy  (int level);
    void CreateRHSProxy (int iSolver, int level);
    void DestroyRHSProxy(int level);

    void StoreRhsByResidual(int iSolver, int level);
    void RecoverResidual   (int iSolver, int level);
    bool CVGNorm(int iSolver, int level = 0);
    bool ISCoarsestGrid(int level);

    bool JudgeIfCalTheZoneForDensityBasedMethod(int zoneID, int iSolver);
    bool JudgeIfCalTheZoneForPressureBasedMethod(int zoneID, int iSolver);
};

class MultiGridFieldProxy
{
private:
    int iSolver, level, nzone;
    FieldProxy **field_proxy;
public:
    MultiGridFieldProxy(int iSolver, int level, int nZone);
    ~MultiGridFieldProxy();
public:
    FieldProxy ** GetFieldProxy() { return field_proxy; }
};

class ZoneFieldProxy
{
private:
    int iSolver, nzone;
    uint_t nlevel;
    MultiGridFieldProxy **field_proxy;
public:
    ZoneFieldProxy(int iSolver, uint_t nLevel, int nZone);
    ~ZoneFieldProxy();
public:
    FieldProxy ** GetFieldProxy(int level);
};

#include "Controller.hxx"

}
