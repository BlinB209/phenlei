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
#include "TurbSolver.h"
#include "Limiter.h"
#include "FaceProxy.h"

namespace PHSPACE
{
class Param_TurbSolverUnstruct;

class TurbSolverUnstr: public TurbSolver
{
private:
    static const int FORMATION_ORGSA_SST = 0;
    static const int FORMATION_SA_NSMB   = 1;
    static const int FORMATION_SA_CFL3D  = 2;
private:
    //! Length scale of XXDES, such as DES/DDES/IDDES.
    RDouble  *DESLength;

    //! Gradient of variables on turbulent equation.
    Gradient *gradientTurbField;

    //! Gradient of velocity.
    Gradient *gradientVelocity;

public:
    TurbSolverUnstr();
    ~TurbSolverUnstr();

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

    //! Overlapping boundary : The flow field variables are interpolated from the interpolating cell to the invalid cell 
    void ReSetOversetBoundary(Grid* gridIn);

    //! Inflow Boundary condition, suitable for mass in, pressure in, and velocity in for turbulence flows.
    void InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Wall Boundary condition.
    void WallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    //! Wall Boundary condition.
    void VisWallWithWallFunctionStandard(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    //! Wall Boundary condition.
    void VisWallWithWallFunctionPAB3D(Grid *gridIn, UnstructBC *bcRegionUnstruct);

#ifdef USE_GMRESSOLVER
    //! GMRESturb Wall Boundary condition.
    void GMRES_WallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    //! GMRESturb Farfield Boundary condition.
    void GMRES_FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
#endif

    //! Farfield Boundary condition.
    void FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    //! Residuals are set to zero
    void ZeroResiduals(Grid *gridIn);

    void RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation);

    void InitMixingPlane(RDouble ***MixingPlaneVar, int Dim1, int Dim2, int Dim3, RDouble Value);

    void AverageMixingPlane(Grid *grid);
    void MixingPlaneDataTransfer(Grid *grid, Grid *NeighborGrid);
    void NonReflective(Grid *grid, Grid *NeighborGrid);
    void SetMixingPlaneData(Grid *grid);

    //! Judge if the flow field variables file already exists in Turb solver.
    bool JudgeIfRestart();
    bool JudgeIfReadAverage();

    //! Initialize the flow field as simulation restarting.
    void InitFlowAsRestart();

    //! Write the file for flowfield restarting.
    void DumpRestart(ActionKey *actkey);

    //! Read the file for flowfield restarting.
    void ReadRestart(ActionKey *actkey);

    //! Write and read restart file.
    void DumpRestartH5(ActionKey *actkey);
    void ReadRestartH5(ActionKey *actkey);

    //! Compute the residual of turbulent equation.
    void GetResidual(ActionKey *actkey);

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
    
    //! Smooth variables of turbulent equation.
    void SMoothTurbulence(Grid *gridIn);
    void SMoothTurbulencePoint(Grid *gridIn, RDouble *primitiveVariable, int i);
    
    //! Compute spectrumRadius of turbulent equation.
    void Diagonal(Grid *grid);
    
    //! Compute inviscid flux of turbulent equation.
    void InviscidFlux(Grid *gridIn);
    
    //! Compute viscous flux of turbulent equation.
    void ViscousFlux(Grid *gridIn);
    
    //! Compute source term of dual time step.
    void DualTimeSource(Grid *gridIn);
    
    //! Update terms used for unsteady flow simulation.
    void UpdateUnsteadyFlow(Grid *gridIn);
    
    //! Convergence norm for the unsteady flow simulation.
    RDouble UnsteadyConvergence(Grid *gridIn);
    
    //! Modify residual values on wall boundary.
    void ModifyResiduals(Grid *grid);
    void ModifyResidualOnWallBoundary(Grid *gridIn);
    
    //! Reset variables of turbulent equation on wall boundary.
    void ResetWallScalar(Grid *gridIn);
    
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
    
    //! Variables of turbulent equation that located in face.
    TurbFaceValue *CreateTurbFaceValue(Grid *gridIn);

    Limiter * CreateLimiter(Grid *grid_in);//!modified by zhangjian

    //! Compute the node value in the turbulent flow.
    //! @param[in] grid         the mesh in the computing.
    void ComputeQTurbNodeValue(Grid *grid_in);

    //! modify the node value for the interpoint in the turbulent flow.
    //! @param[in] grid         the mesh in the computing.
    void ModifyQTurbNodeValue(Grid *grid_in);
private:
    FaceProxy * CreateFaceProxyTurb(Grid *grid_in);
public:
    void SourceFluxOneEquation(Grid *gridIn);
    void SourceFluxTwoEquation(Grid *gridIn);

    //! Source term of one equation turbulence model with many variants.
    void SourceFluxOneEquationOriginal(Grid *gridIn);
 
 #ifdef USE_GMRESSOLVER
    //! GMRESturb
    void GMRES_SourceFluxOneEquationOriginal(Grid *gridIn);
#endif

    void SourceFluxOneEquationOriginalCFL3D(Grid *gridIn);
    //void SourceFluxOneEquationEdwards(Grid *gridIn);
    //void SourceFluxOneEquationNew1(Grid *gridIn);
    //void SourceFluxOneEquationNew(Grid *gridIn);

    //! Get face value for viscous flux of computation.
    void GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
  
#ifdef USE_GMRESSOLVER
    //! GMRESCoupled
    template <typename T>
    void GMRES_GetVisFaceValue(UnstructBCSet *unstructBCSet, int iFace, int nBoundFace, int nTurbulenceEquation, RDouble KW_sigma,  T rhol, T rhor, T pressurel, T pressurer, T *turbQl, T *turbQr, T& mul, T& mlt, T& viscousLaminarl, T& viscousLaminarr);
   
    //! GMRESturb
    void GMRES_ComputeVisflux (Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
 #endif

    //! Compute viscous flux of turbulent equations.
    void ComputeVisflux (Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Compute the spectrum radius of one equation turbulence model.
    void SpectrumRadiusOfOneEquation(Grid *gridIn);
    
    //! Compute the spectrum radius of two equation turbulence model.
    void SpectrumRadiusOfTwoEquation(Grid *gridIn);
    //! Compute the blending function for EASM model.
    //void EASMBlending(Grid *gridIn, RDouble *cross, RDouble *blend);
    
    //! Compute the crossing term in two equation turbulence model.
    void Crossing(Grid *gridIn);
    
    //! Compute the blending function in two equation turbulence model.
    void Blending(Grid *gridIn);

    //void Crossing(Grid *gridIn, RDouble *cross);
    //void Blending(Grid *gridIn, RDouble *cross, RDouble *blend);

    //void SourceFluxEASM(Grid *gridIn);
    //void Gatski_Speziale_Model(Grid *gridIn);
    //void EASMko2003_Model(Grid *gridIn);
    //void EASMko2005_Model(Grid *gridIn);

    //! Compute turbulent viscous coefficient. 
    void ComputeViscousCoeff(Grid *grid);
    void Viscosity(Grid *grid);
    
    //! Initialize wall distance for coarse grid.
    void InitCGrid(Grid *fineGridIn, Grid *coarseGridIn);
    
    //! Upload interface value for variables in turbulence model. 
    void UploadInterfaceValue(ActionKey *actkey);
    void UploadInterfaceData(ActionKey *actkey);    //! Bell 20120910 add
    
    //! Download interface value for variables in turbulence model. 
    void DownloadInterfaceValue(ActionKey *actkey);
    void DownloadInterfaceData(ActionKey *actkey);    //! Bell 20120910 add

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

#ifdef USE_GMRESSOLVER
    //! GMRESturb
    void GMRES_NNDFlux(Grid *gridIn, FaceProxy *FaceProxy, int nst, int ned);

    //! GMRESCoupled fix left and right variables of turbulence model and NS equations, required for GMRES by AD
    template <typename T>
    void FixBoundaryQlQr(UnstructBCSet *unstructBCSet, int iFace, int nTurbulenceEquation, int nEquation, T *primQl, T *primQr, T *turbQl, T *turbQr);

    //! GMRES solver -- a linear system solver
    //! This method can be referenced "Sen Zhang, Boqian Wang, et al. Implementation of a Newton-Krylov Algorithm in the Open-source Solver PHengLEI[C]."
    //! GPPS Hong Kong , October 17-19, 2023.
    void GMRESSolver(Grid *gridIn, FieldProxy *dqProxy);
#endif

    //! Get left and right variables of turbulence model and NS equation.
    void GetQlQr(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);

    //! Obtain the gradient of variables in turbulence model.
    void GetGradientField(Grid *gridIn);
    
    //! Compute gradient of variables in turbulence model.
    void ComputeGradient(Grid *grid);

    //! Create and initialize controlparameters
    LIB_EXPORT void InitControlParameters();

    //! Get control paramters.
    LIB_EXPORT Param_TurbSolverUnstruct *GetControlParameters();
    
    //! Write data needed communication at interface to datacontainer.
    void CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn);
    
    //! Read data from datacontainer after communication.
    void DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn);

    ////! Implement: to compress turbulenct dq data from grid to date dataContainer.
    //void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    //void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);

private:

    //! Get variables at face interpolated from left and right side for turbulence model equation.
    void GetQlQrTurbulence(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Get variables at face interpolated from left and right side for NS equation.
    void GetQlQrNS(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Compute the weights used to average dqdx, dqdy, etc.
    void ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int nst, int ned);
    
    //! Write the interface information. 
    void DumpInterfaceInfoBell();         // Bell 20120910 add
    
    //!  Decide whether modify the computation of term with coefficient cb2 in original SA turbulence model.
    //!  If modify the SA turbulence model: For the original SA turbulence model,
    //!  transform the gradient square term into the viscous flux, using CFL3D or NSMB method.\n
    //!    -# 0: SST or original SA.
    //!    -# 1: NSMB formation SA.
    //!    -# 2: CFL3D formation SA.
    int TurbEquationFormation() {return this->Turb_Equation_Formation;}
    
    //! Compute length scale of one equation turbulence model.
    void ComputeLengthScaleofOneEquationModel(Grid *gridIn);
    
    //! Compute length scale with original Deatched Eddy Simulation method.
    void ComputeDESLength(Grid *gridIn);
    
    //! Compute length scale with Delayed Deatched Eddy Simulation method.
    void ComputeDDESLength(Grid *gridIn);
    
    //! Compute length scale with improved Delayed Deatched Eddy Simulation method.
    void ComputeIDDESLength(Grid *gridIn);
    
    //! Obtaion the length scale in turbulence simulation.
    RDouble *GetLengthScale(Grid *gridIn);
    
    //! Obtaion the length scale of DES.
    RDouble *GetDESLength() { return this->DESLength; }
    
    //! Set the length scale of DES.
    void SetDESLength(RDouble *length) { this->DESLength = length; }
    
    //! Read the averaged variables of statistical flow field.
    void ReadStatisticalFlow(ActionKey *actkey);

    //! Create proxy for variables of NS equations that located in face.
    FaceProxy * CreateFaceProxyNS(Grid *gridIn);
    
    //! Create proxy for variables of turbulent equations that located in face.
    FaceProxy * CreateFaceProxyTurbulence(Grid *gridIn);

private:
    //! SA equation formation/variety.
    //!    -# 0: SST
    //!    -# 1: NSMB formation
    //!    -# 2: CFL3D formation
    int Turb_Equation_Formation;

    int machZeroExp;
};

//! Compute length scale of LES for DES based on SA turbulence model.
void ComputeLESLengthofSA(Grid *gridIn, RDouble *lengthOfLES, RDouble *lowReynoldsNumberCorrection, int DESType);

//! Compute low Reynolds number correction term for SA turbulence model when implementing DES.
void ComputeLowReynoldsNumberCorrection(Grid *gridIn, RDouble *lowReynoldsNumberCorrection);
}
