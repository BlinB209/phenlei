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
//! @file      TurbSolverStruct.h
//! @brief     turbulence solver for struct grid.
//! @author    He Xin, He Kun, He Xianyao, Liu Jian, Zhang Zipei, Li peng, Ma Yankai.

#pragma once
#include "TransitionSolver.h"

namespace PHSPACE
{
class StructBC;
class Param_TransitionSolver;

class TransitionSolverStr : public TransitionSolver
{
public:
    TransitionSolverStr();
    ~TransitionSolverStr();

    //! To obtain the turbulent residual according to the flux between the cell face.
    //! @param[in out ]: gridIn denotes the current computational regions of the grids, in which turbulent residual are stored.
    //! @param[in ]: fluxProxy denotes the flux proxy of the cell face.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void LoadFlux(Grid *gridIn,  FieldProxy *fluxProxy, int iSurface);

    //! To obtain the turbulent viscous flux of the cell face.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out ]: fluxProxy denotes the flux proxy of the cell face.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CompVisFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);

    //! To obtain the two equation turbulent viscid flux of the cell face.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out ]: fluxProxy denotes the flux proxy of the cell face.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CompVisFluxTwoEquation(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);

    //! To set the turbulent residual to zero of the Special cell, which is inherited from class NSSolver.
    //! @param[in ]: gridIn denotes the current computational regions of the grids, in which turbulent residual are stored.
    void ZeroResidualOfSpecialCells(Grid * gridIn);

    //! To set the transition variables from primitive form to conservative form in two equation model.
    //! @param[in ]: gridIn denotes the current computational regions of the grids, in which turbulent variables are stored.
    void ChangeTransitionQ (Grid *gridIn);

    //! To set the transition variables from conservative form to primitive form in two equation model.
    //! @param[in ]: gridIn denotes the current computational regions of the grids, in which turbulent variables are stored.
    void RecoverTransitionQ(Grid *gridIn);

    //! To obtain the two equation turbulent inviscid flux of the cell face.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: qLeftProxy denotes the turbulent q proxy to the left of the cell face.
    //! @param[in ]: qRightProxy denotes the turbulent q proxy to the right of the cell face.
    //! @param[in ]: qpmvl_proxy denotes the laminar q proxy to the left of the cell face.
    //! @param[in ]: qpmvr_proxy denotes the laminar q proxy to the right of the cell face.
    //! @param[out ]: fluxProxy denotes the flux proxy of the cell face.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void ComputeInviscidFlux(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qpmvl_proxy, FieldProxy *qpmvr_proxy, FieldProxy *fluxProxy, int iSurface);

    //! To obtain the two equation turbulent inviscid flux of the cell face by NND method, called by function ComputeInviscidFlux.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: qLeftProxy denotes the turbulent q proxy to the left of the cell face.
    //! @param[in ]: qRightProxy denotes the turbulent q proxy to the right of the cell face.
    //! @param[in ]: qLaminarLeftProxy denotes the laminar q proxy to the left of the cell face.
    //! @param[in ]: qpmvr_proxy denotes the laminar q proxy to the right of the cell face.
    //! @param[out ]: fluxProxy denotes the flux proxy of the cell face.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void NNDFlux(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, FieldProxy *fluxProxy, int iSurface);

    //! To obtain the turbulent MUSCL interpolate value of the cell face in iSurface direction.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out ]: qLeftProxy denotes the turbulent q proxy to the left of the cell face.
    //! @param[out ]: qRightProxy denotes the turbulent q proxy to the right of the cell face.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void GetInviscidFaceValue(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, int iSurface);

    //! To obtain the laminar velocity of the cell face  for turbulent equation in iSurface direction.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out ]: qLaminarLeftProxy denotes the velocity proxy to the left of the cell face.
    //! @param[out ]: qLaminarRightProxy denotes the velocity proxy to the right of the cell face.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void GetInvFaceQValueforTransition(Grid *gridIn, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, int iSurface);

    //! To modify the interpolate turbulent vaule on the boundary, especially on the SOLID SURFACE.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in out ]: qLeftProxy denotes the turbulent q proxy to the left of the cell face.
    //! @param[in out ]: qRightProxy denotes the turbulent q proxy to the right of the cell face.
    //! @param[in out ]: qLaminarLeftProxy denotes the velocity proxy to the left of the cell face.
    //! @param[in out ]: qLaminarRightProxy denotes the velocity proxy to the right of the cell face.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CorrectFaceVar(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, int iSurface); 

    //! To get the gradient of turbulent variables at the cell center.
    //! @param[in ]: gridIn denotes the current computational regions of the grids, on which dependent variables are stored.
    void ComputeGradientCellCenter(Grid *gridIn);

    //! To get the gradient of turbulent variables at the cell center.
    //! @param[in ]: gridIn denotes the current computational regions of the grids, on which dependent variables are stored.
    //void ReconGrad(Grid *gridIn, int iSurface);

    //! To modified the boundary value of the turbulent variables's gradient at the cell center.
    //! @param[in ]: gridIn denotes the current computational regions of the grids, on which dependent variables are stored.
    void SetCellCenterGradientBoundary(Grid *gridIn);

    //! To get the gradient of turbulent variables at the cell face from the cell center gradient.
    //! @param[in ]: gridIn denotes the current computational regions of the grids, on which dependent variables are stored.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void ReconGradVolWeightedwithCorrection(Grid *gridIn, int iSurface);

    //! To correct the gradient of turbulent variables at the special physical boundaries, such as viscous wall, symmetry.
    //! @param[in ]: gridIn denotes the current computational regions of the grids, on which dependent variables are stored.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void GetGradientAtFace_CorrectionAtPhysicalBoundary(Grid *gridIn, int iSurface);

    void RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation);//zhangyong

    void AverageMixingPlane(Grid *grid) {};
    void MixingPlaneDataTransfer(Grid *grid, Grid *NeighborGrid) {};
    void NonReflective(Grid *grid, Grid *NeighborGrid) {};
    void SetMixingPlaneData(Grid *grid) {};

    //! Interface: to put the boundary data to DataContainer for parallel communication.
    //! @param[in ]: actkey contains the control information for parallel communication.
    void UploadInterfaceData(ActionKey *actkey);

    //! implement: to put the boundary data to DataContainer for parallel communication, called by function UploadInterfaceData.
    //! @param[in ]: actkey contains the control information for parallel communication.
    void UploadInterfaceValue(ActionKey *actkey){};

    //! Interface: to get the boundary data from DataContainer for parallel communication.
    //! @param[in ]: actkey contains the control information for parallel communication.
    void DownloadInterfaceData(ActionKey *actkey);

    //! implement: to get the boundary data from DataContainer for parallel communication, called by function DownloadInterfaceData.
    //! @param[in ]: actkey contains the control information for parallel communication.
    void DownloadInterfaceValue(ActionKey *actkey){};

    //! Implement: to compress turbulent dq data from grid to date dataContainer.
    void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);

    //! Todo: overset information communication hasn't carried out yet.
    //! @param[in ]: actkey contains the control information for parallel communication.
    void UploadOversetData(ActionKey *actkey);

    void DownloadOversetData(ActionKey *actkey);

    //! To set all of the value of the turbulent field proxy equal to the input value
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: fieldProxy denotes the target turbulent field proxy which will be covered by input value .
    //! @param[in ]: value denotes the source value.
    void FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value);

    //! To set all of the value of the fieldProxy equal to the input value
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: targetFieldProxy denotes the target turbulent field proxy which will be covered by the source field proxy.
    //! @param[in ]: sourceFieldProxy denotes the source turbulent field proxy.
    void FillField(Grid *gridIn, FieldProxy *targetFieldProxy, FieldProxy *sourceFieldProxy);

    //! To creat a turbulent fieldProxy and get data by field name.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in]: fieldName denotes the field name to be set into fieldProxy.
    //! @return: the turbulent field proxy created.
    FieldProxy *GetFieldProxy(Grid *gridIn, const string &fieldName);

    //! To creat a turbulent fieldProxy without initionalization.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @return: the turbulent field proxy created.
    FieldProxy *CreateFieldProxy(Grid *gridIn);

    //! To creat a laminar fieldProxy without initionalization.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @return: the laminar field proxy created.
    FieldProxy *CreateFlowFieldProxy(Grid *gridIn);

    //! To creat a turbulent residual proxy and set by "res_turb" data in gridIn.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @return: the turbulent residual proxy created.
    FieldProxy *GetResidualProxy(Grid *gridIn);

    //! To set the turbulent right hand side proxy by "res_turb" in gridIn.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: rightHandSideProxy: the turbulent right hand side proxy.
    void StoreRhsByResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);

    //! To initionalization the "res_turb" data in gridIn by turbulent right hand side proxy.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: rightHandSideProxy: the turbulent right hand side proxy.
    void InitResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);

    //! To set the "res_turb" data in gridIn by turbulent right hand side proxy.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: rightHandSideProxy: the turbulent right hand side proxy.
    void RecoverResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);

    //! Create and initialize control parameters.
    LIB_EXPORT void InitControlParameters();

    //! Get control parameters.
    LIB_EXPORT Param_TransitionSolver *GetControlParameters();

    //! Judge if the flow field variables file already exists in Turb solver.
    bool JudgeIfRestart();
    bool JudgeIfReadAverage();

    bool JudgeIfProtectedRestart();

    void InitFlowAsReadingProtectedRestart(ActionKey *actkey);

    //! To write the flow data to restart file by HDF5 "fts" form for continue to simulation.
    //! @param[in ]: actkey contains the control information.
    void DumpRestartH5(ActionKey *actkey);

    //! To read the flow data from restart file by HDF5 "fts" form  for continue to simulation.
    //! @param[in ]: actkey contains the control information.
    void ReadRestartH5(ActionKey *actkey);

private:

    //! To allocate the internal storage of field variables/residual/spectum array and other global data array.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void AllocateGlobalVar(Grid *gridIn);

    //! To deallocate the internal storage of array.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void DeAllocateGlobalVar(Grid *gridIn);

    //! To initionalize the flow flield according to the inflow condition
    void InitFlowAsRestart();

    //! To set the flow field value of the external boundary condition of the computational region. \n
    //! It looped all the internal boundary and was realized by call the specific boundary functions.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void Boundary(Grid *gridIn);

    //! To set the inflow boundary condition by setting the flow field value equal to the far field value.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: structBC denotes the current boundary condition region which need to be set.
    void InFlowBC  (Grid *gridIn, StructBC *structBC);

    //! To set the outflow boundary condition by 2nd order interpolation.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: structBC denotes the current boundary condition region which need to be set.
    void OutFlowBC (Grid *gridIn, StructBC *structBC);

    //! To set the symmetry boundary condition by 1st order interpolation.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: structBC denotes the current boundary condition region which need to be set.
    void SymmetryBC(Grid *gridIn, StructBC *structBC);

    //! To set the viscid wall boundary condition by setting the turbulent variables equal to zero or other specific value.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: structBC denotes the current boundary condition region which need to be set.
    void VisWall(Grid *gridIn, StructBC *structBC);

    void VisWallWithWallFunctionStandard(Grid *gridIn, StructBC *structBC);
    void VisWallWithWallFunctionPAB3D(Grid *gridIn, StructBC *structBC);

    //! To set the farfield boundary condition according to riemann method.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: structBC denotes the current boundary condition region which need to be set.
    void FarFieldBC(Grid *gridIn, StructBC *structBC);

    //! To implement the standard wall function while y+ isn't meet the need.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: structBC denotes the current boundary condition region which need to be set.
    void WallFunctionStandard(Grid *gridIn, StructBC *structBC, int is, int js, int ks, int it, int jt, int kt, int i, int j, int k);

    //! To implement the Pab3D wall function while y+ isn't meet the need.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: structBC denotes the current boundary condition region which need to be set.
    void WallFunctionPab3D(Grid *gridIn, StructBC *structBC, int is, int js, int ks, int it, int jt, int kt, int i, int j, int k);

    //! To obtain the turbulent residual according to the right hand side.
    //! @param[in out ]: gridIn denotes the current computational regions of the grids, in which turbulent residual are stored.
    //! @param[in ]: rightHandSideProxy denotes the right hand side proxy of the cell face.
    void LoadResiduals(Grid *gridIn, FieldProxy *rightHandSideProxy);

    //! To set the turbulent residual zero.
    //! @param[in out ]: gridIn denotes the current computational regions of the grids, in which turbulent residual are stored.
    void ZeroResiduals(Grid *gridIn);

    //! To set the turbulent spectum zero.
    //! @param[in out ]: gridIn denotes the current computational regions of the grids, in which turbulent spectum are stored.
    void InitSpectrum(Grid *gridIn);

    //! Interface: to get the turbulent source flux of one equation model.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void SourceFluxOneEquation(Grid *gridIn);

    //! Implement: to get the turbulent source flux of one equation model by original method.
    //! Don't use yet.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void SourceFluxOneEquationOriginal(Grid *gridIn);

    //! Implement: to get the turbulent source flux of one equation model by new method.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void SourceFluxOneEquationNew(Grid *gridIn);

    //! To get the turbulent source flux of two one equation model.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void SourceFluxTwoEquation(Grid *gridIn);

    //! To get the turbulent inviscid flux.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void InviscidFlux(Grid *gridIn);

    //! To get the turbulent viscid flux.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void ViscousFlux(Grid *gridIn);

    //! To get the turbulent dual time source residual.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void DualTimeSource(Grid *gridIn);

    //! To update the unsteady turbulent flow field variables such as qTurbulence, residual, volume.
    //! q(n+1) = q(n); q(n) = q(n-1)..
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void UpdateUnsteadyFlow(Grid *gridIn);

    RDouble UnsteadyConvergence(Grid *gridIn);

    //! Interface: to get the turbulent spectum.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void Diagonal(Grid *gridIn);

    //! Implement: to get the turbulent spectum of one equation model, called by function Spectrum.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void SpectrumRadiusOfOneEquation(Grid *gridIn);

    //! Implement: to get the turbulent spectum of two equation model, called by function Spectrum.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void SpectrumRadiusOfTwoEquation(Grid *gridIn);

    //! To load flow variables stored in grid to qProxy.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: qProxy denotes the turbulent field proxy.
    void LoadQ(Grid *gridIn, FieldProxy *qProxy);

    //! To get the turbulent residual for Runge-Kutta method.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coef);

    //! Implicit LU-SGS time integration method.
    void SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    void SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);

    //! To update the turbulent flow field variables for steady flow.
    //! q(n+1) = q(n)+dq.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: qProxy denotes the turbulent field proxy.
    //! @param[in ]: dqProxy denotes the of turbulent field proxy increment of one iteration step.
    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);

    //! To set turbulent flow field variables's value at the corner point.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void CornerPoint(Grid *gridIn);

    //! To get the distance between the two point(i, j, k) and (i+1, j+1, k+1).
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: i, j, k denotes grid point index.
    //! @return: the distance between the two point(i, j, k) and (i+1, j+1, k+1).
    RDouble GetDistance(Grid *gridIn, int i, int j, int k);

    //! To get maximum residual for output.
    //! @param[in ]: actkey contains the control parameter.
    void GetResidual(ActionKey *actkey);
};

}
