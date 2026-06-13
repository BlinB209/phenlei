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
//! @file      TurbSolverUnstr.h
//! @brief     turbulence solver for unstruct grid.
//! @author    Bell, He Xin, Zhang Laiping, Dr. Wang, Zhang Yang, He Kun, 
//!            Wan Yunbo, Xu Gang, Zhang Yaobing.

#pragma once
#include "TransitionSolver.h"
#include "Limiter.h"
#include "FaceProxy.h"

namespace PHSPACE
{
class Param_TransitionSolver;

class TransitionSolverUnstr: public TransitionSolver
{
private:
    //! Gradient of variables on turbulent equation.
    Gradient *gradientTransitionField;

    //! Gradient of velocity.
    Gradient *gradientVelocity;

public:
    TransitionSolverUnstr();
    ~TransitionSolverUnstr();

    //! Allocate and delete global variables.
    void AllocateGlobalVar(Grid *gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);

    //! Load rhs to residual which stored in grid
    void LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy);

    //! Variables such as specTurb are set to zero.
    void InitSpectrum(Grid *gridIn);

    //! Boundary condition.
    void Boundary(Grid *gridIn);

    //! Outflow Boundary condition, suitable for mass out, pressure out, and velocity out in turbulence flows.
    void OutflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    //! Inflow Boundary condition, suitable for mass in, pressure in, and velocity in for turbulence flows.
    void InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Wall Boundary condition.
    void WallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Farfield Boundary condition.
    void FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Residuals are set to zero
    void ZeroResiduals(Grid *gridIn);

    //! Judge if the flow field variables file already exists in Turb solver.
    bool JudgeIfRestart();
    bool JudgeIfReadAverage();

    //! Initialize the flow field as simulation restarting.
    void InitFlowAsRestart();

    //! Write and read restart file.
    void DumpRestartH5(ActionKey *actkey);
    void ReadRestartH5(ActionKey *actkey);

    //! Compute the residual of turbulent equation.
    void GetResidual(ActionKey *actkey);

    void RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation);

    void InitMixingPlane(RDouble ***MixingPlaneVar, int Dim1, int Dim2, int Dim3, RDouble Value);

    void AverageMixingPlane(Grid *grid);
    void MixingPlaneDataTransfer(Grid *grid, Grid *NeighborGrid);
    void NonReflective(Grid *grid, Grid *NeighborGrid);
    void SetMixingPlaneData(Grid *grid);

    //! Restrict defect from fine to coarse grid.
    void RestrictDefect(Grid *fineGridIn, Grid *coarseGridIn);

    //! Correct variable in fine grid using correction in coarse grid
    void CorrectFineGrid(Grid *fineGridIn, Grid *coarseGridIn);

    //! Interpolate variables of turbulence model in coarse grid to that in fine grid.
    void InterpolatFineGrid(Grid *fineGridIn, Grid *coarseGridIn);

    //! Put correction back on the coarse grid.
    void PutCorrectionBack(Grid *gridIn, FieldProxy *qProxy);

    //! Load NS variables stored in grid to q.
    void LoadQ(Grid *gridIn, FieldProxy *qProxy);

    //! Compute the residuals with RungeKutta method.
    void RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coefficient);

    //! Update flow field variables of turbulene equation.
    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);

    //! LU-SGS forward sweep.
    void SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    
    //! LU-SGS backward sweep.
    void SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    
    //! Set the values at ghost cells.
    void SetGhostDQLUSGS(Grid *gridIn, RDouble **dq);    // Bell 20130401 add
    
    void SetWallBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq);

    void SetFarfieldBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq);
    
    //! Smooth variables of transition equation.
    void SMoothTransition(Grid *gridIn);
    void SMoothTransitionPoint(Grid *gridIn, RDouble *primitiveVariable, int i);
    
    //! Compute spectrumRadius of transition equation.
    void Diagonal(Grid *grid);
    
    //! Compute inviscid flux of transition equation.
    void InviscidFlux(Grid *gridIn);
    
    //! Compute viscous flux of transition equation.
    void ViscousFlux(Grid *gridIn);
    
    //! Compute source term of dual time step.
    void DualTimeSource(Grid *gridIn);
    
    //! Update terms used for unsteady flow simulation.
    void UpdateUnsteadyFlow(Grid *gridIn);
    
    //! Convergence norm for the unsteady flow simulation.
    RDouble UnsteadyConvergence(Grid *gridIn);

    //! Restrict variables of turbulent equation on fine and coarse grid.
    void RestrictAllQ(Grid *fineGridIn, Grid *coarseGridIn);
    
    //! Compute flux of the turbulent equation.
    void LoadFlux(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Set value to flowfiled variables.
    void FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value);
    void FillField(Grid *gridIn, FieldProxy *field1Proxy, FieldProxy *field2Proxy);
    
    //! Create proxy for variables of turbulent equation located in cell center.
    FieldProxy *CreateFieldProxy(Grid *gridIn);
    
    //! Get the corresopnding proxy.
    FieldProxy *GetFieldProxy(Grid *gridIn, const string &fieldName);
    
    //! Get the proxy of residuals.
    FieldProxy *GetResidualProxy(Grid *gridIn);
    
    //! Load the residuals of turbulent equation.
    void RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy);
    
    //! Initialize the residuals of turbulent equation.
    void InitResidual(Grid *gridIn, FieldProxy *rhsProxy);
    
    //! Rhs is set equal to residuals of turbulent equation. 
    void StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy);
    
    //! Create proxy for variables of turbulent and NS equation that located in face.
    FaceProxy *CreateFaceProxy(Grid *gridIn);
    
    //! Variables of transition equation that located in face.
    TransitionFaceValue *CreateTransitionFaceValue(Grid *gridIn);

    Limiter * CreateLimiter(Grid *grid_in);//!modified by zhangjian

    //! Compute the node value in the transition flow.
    //! @param[in] grid         the mesh in the computing.
    void ComputeQTransitionNodeValue(Grid *grid_in);

    //! modify the node value for the interpoint in the turbulent flow.
    //! @param[in] grid         the mesh in the computing.
    void ModifyQTransitionNodeValue(Grid *grid_in);
public:
    void SourceFluxTwoEquation(Grid *gridIn);

    //! Get face value for viscous flux of computation.
    void GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);

    //! Compute viscous flux of turbulent equations.
    void ComputeVisflux (Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);

    //! Initialize wall distance for coarse grid.
    void InitCGrid(Grid *fineGridIn, Grid *coarseGridIn);
    
    //! Upload interface value for variables in turbulence model. 
    void UploadInterfaceValue(ActionKey *actkey);
    void UploadInterfaceData(ActionKey *actkey);    // Bell 20120910 add
    
    //! Download interface value for variables in turbulence model. 
    void DownloadInterfaceValue(ActionKey *actkey);
    void DownloadInterfaceData(ActionKey *actkey);    // Bell 20120910 add

    //! Upload the interpoint variables for send.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void UploadInterpointData(ActionKey *actkey);

    //! Download the interpoint variables from the receive data.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void DownloadInterpointData(ActionKey *actkey);

    void CommunicationInterpointWeight(ActionKey *actkey);
    void DownloadInterpointWeight(ActionKey *actkey);

    //! Upload overset data at interface for variables in turbulence model. 
    void UploadOversetData(ActionKey *actkey);

    //! Download overset data at interface for variables in turbulence model.
    void DownloadOversetData(ActionKey *actkey);

    //! Compute inviscid flux of turbulence model equation.
    void ComputeInviscidFlux(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Inviscid flux scheme for turbulence model equation.
    void NNDFlux(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Get left and right variables of turbulence model and NS equation.
    void GetQlQr(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);

    //! Obtain the gradient of variables in turbulence model.
    void GetGradientField(Grid *gridIn);
    
    //! Compute gradient of variables in turbulence model.
    void ComputeGradient(Grid *grid);

    //! Create and initialize controlparameters
    LIB_EXPORT void InitControlParameters();

    //! Get control paramters.
    LIB_EXPORT Param_TransitionSolver *GetControlParameters();

    ////! Implement: to compress transition dq data from grid to date dataContainer.
    void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);

private:

    //! Get variables at face interpolated from left and right side for turbulence model equation.
    void GetQlQrTransition(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Get variables at face interpolated from left and right side for NS equation.
    void GetQlQrNS(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Compute the weights used to average dqdx, dqdy, etc.
    void ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);

    //! Read the averaged variables of statistical flow field.
    void ReadStatisticalFlow(ActionKey *actkey);

    //! Create proxy for variables of NS equations that located in face.
    FaceProxy * CreateFaceProxyNS(Grid *gridIn);
    
    //! Create proxy for variables of transition equations that located in face.
    FaceProxy * CreateFaceProxyTransition(Grid *gridIn);

private:

    int machZeroExp;
};

}
