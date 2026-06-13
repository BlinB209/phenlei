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
//! @file      NSSolverStruct.h
//! @brief     NS solver for struct grid.
//! @author    He Xin, Bell, He Kun, He Xianyao, Wan Yunbo, Liu Jian, Xu Qingxin, Zhang Zipei, Li peng, Ma Yankai, 
//!            Min Yaobing, Meng Liyuan, Xu Gang, Guo Yongyan, Zhang Yang.

#pragma once
#include "NSSolver.h"
using namespace std;

namespace PHSPACE
{
class GeomProxy;
class StructBC;
class Param_NSSolverStruct;
class GradientCellCenter;
class Gradient;

class StructIter
{
public:
    int i0,j0,k0;
    int i1,j1,k1;
    int il1,jl1,kl1;
    int ist, ied, jst, jed, kst, ked;
    int nsurf;
    int nsize;
  //  int surfaceSize; //  Size of  surface
};

class NSSolverStruct : public NSSolver
{
private:
    GradientCellCenter *gradientCellCenterUVWT[5];
    int systemgridtype;
public:
    int nWallBC; //The number of surface boundary conditions.
    int nMaxSurfaceCell; //The number of surface cells.
    int nIonizationFlag;
    NSSolverStruct();
    ~NSSolverStruct();
private:
    void ReadParameter();
    void AllocateGlobalVar  (Grid *gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);
    void GetSurfaceCellNumber(Grid *gridIn);
    void InitializeSurfaceParameters(Grid *gridIn);
    void InitializeSurfaceParametersFromFile(Grid *gridIn);
    void ChangeSurfaceTemperature(Grid *gridIn);
    void ChangeSurfaceMassFractions(Grid *gridIn);

    //! Judge if the flow field variables file already exists in NS solver.
    bool JudgeIfRestart();
    bool JudgeIfReadAverage();

    bool JudgeIfProtectedRestart();

    void InitFlowAsRestart();
    void InitialSetPartialParameter(StructGrid *grid, bool flowTag = true);

    void InitGrad(Grid *gridIn);
    void InitCGrid(Grid *fineGridIn, Grid *coarseGridIn);
    void ZeroResiduals(Grid *grid);

    void CreateStandardADT(Grid *grid, DataStruct_AdtTree<int,RDouble> *&adtTree, RDouble &tol);

    void CompNodeVar(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy);
    void InterfaceValueFix(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy);
    void WallValueFix(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy);
    
    void CompNodeVar(Grid *gridIn, RDouble4D &qn, int m, RDouble4D &q, int n);
    void CompNodeVar(Grid *gridIn, RDouble4D &qn, int m, RDouble3D &q);
    
    void LoadQ(Grid *grid, FieldProxy *qProxy);
    void LoadResiduals(Grid *gridIn, FieldProxy *rightHandSideProxy);
    void ComputeMinTimeStep(Grid *grid, RDouble &minDt, RDouble &maxDt);
    void ReduceMaxTimeStep(Grid *grid, RDouble globalMinDt);
    void LocalTimeStep(Grid *gridIn);
    void GlobalTimeStep(Grid *gridIn);
    void LocalGlobalTimeStep(Grid *gridIn);
    void SpectrumRadius(Grid *grid);
    void SpectrumRadiusInviscid(Grid *gridIn);
    void SpectrumRadiusViscous(Grid *gridIn);
    void ConservativePreconMatrix(Grid *gridIn, int sign);

    void Diagonal(Grid *gridIn);    //! Source linear.
    //Compute the local time scale.
    void ComputeTimeScale(Grid *gridIn);
    //Compute the local time scale of flow.
    void ComputeFlowTimeScale(Grid *gridIn);

    //! Set the boundary conditions.
    void Boundary(Grid *grid);
    void CalculateBoundaryData();
    void CalPressureFactor(Grid *gridIn);
    void GetHeatTransferCoeff(Grid *gridIn);
    void InviscidFlux(Grid *gridIn);
    void ZeroResidualOfSpecialCells(Grid *gridIn);    //! Added By Guo Yongheng 20160825.
    void ViscousFlux(Grid *gridIn);
    void DualTimeSource(Grid *gridIn);

    //! Compute the chemical source terms.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void ChemicalSource(Grid *gridIn);  //! Source.
    
    void gravitySource(Grid* gridIn) {};  //! Source.
    void PorousMediumSource(Grid *gridIn) {};//! Source

    void RotatingSource(Grid* gridIn) {};  //! Source.

    void AverageMixingPlane(Grid *grid) {};
    void MixingPlaneDataTransfer(Grid *grid, Grid *Neighborgrid) {};
    void NonReflective(Grid *grid, Grid *NeighborGrid) {};
    void SetMixingPlaneData(Grid *grid) {};

    //! Compute the distance between the center of the grid cell on the first layer and the center of the surface cell.
    void ComputeFirstLayerGridHeight(Grid *gridIn);
    void UpdateUnsteadyFlow(Grid *gridIn);
    RDouble UnsteadyConvergence(Grid *gridIn);

    void ComputeViscousCoeff(Grid *gridIn);
    void ComputeViscousCoefficientWithoutChemical(Grid *gridIn);
    void ComputePreconditionCoefficient(Grid *gridIn);
    void ComputePreconditionCoefficientUnsteady(Grid *gridIn);

    void ComputeGamaAndTemperature(Grid *gridIn);                 //Gama, Temperature
    void RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coef);
    
    void SolutionFix(Grid *gridIn, RDouble *prim, int i, int j, int k);
    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);    //! Primitive vs Conservative.

    void UpdateQlQrOnlyforMixGrid(Grid *grid);
    void UpdateInviscidfluxOnlyforMixGrid(Grid *grid);
    void UpdateViscousfluxOnlyforMixGrid(Grid *grid);
    void LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn);

    void LoadFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    void CompVisFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    //! To correct the viscous fluxes of solid boundaries under the non-equilibrium condition.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: fluxProxy is an array that stores the inviscid fluxes of grid faces.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CorrectViscousFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);

    void CompVisFluxnew(Grid *gridIn, Gradient *gradientUVWT, FieldProxy *fluxProxy, int iSurface);
    void CompVisFluxnew(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    void CompVisFluxLES(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    void ComputeInviscidFlux(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);
    void InviscidFluxWrap(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, INVScheme inviscidScheme, int iSurface);
    //! To correct the inviscid fluxes of solid boundaries under the non-equilibrium condition.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: fluxProxy is an array that stores the inviscid fluxes of grid faces.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CorrectInviscidFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);

    //! To correct the primitive variables of faces nearby the wall.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: inviscidFaceProxy is an array that stores the primitive variables of grid faces.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CorrectFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);
    void CorrectFaceVarR(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);
    void CorrectChemicalFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);

    void GetInvFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);
    
    //! Compute gradient of primitive variables & temperature for viscous flux computation.
    void ComputeGradient(Grid *grid);
    void ComputeGradientCellCenter(Grid *grid);
    void ComputeGradientCellCenterForLES(Grid *grid);
    void ComputeGradientCellCenterOfruvwptOnlyForMixGrid(Grid *grid);
    void ReconGrad(Grid *gridIn, int iSurface);

    void RegisterCFDSolverInterfaceField();
    void RegisterOversetField();
    void ReleaseCFDSolverInterfaceField();
    void RotateVectorFromInterface(Grid *grid, const int &neighborZoneIndex, const int &neqn);

    void CompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);

    void RestrictAllQ(Grid *fineGridIn, Grid *coarseGridIn);
    void RestrictDefect(Grid *fineGridIn, Grid *coarseGridIn);
    void PutCorrectionBack(Grid *grid, FieldProxy *qProxy);
    void CorrectFineGrid(Grid *fineGridIn, Grid *coarseGridIn);
    void InterpolatFineGrid(Grid *fineGridIn, Grid *coarseGridIn);
    void prol   (Grid *fineGridIn, Grid *coarseGridIn, FieldProxy *workProxy);
    void prolhic_eric(Grid *fineGridIn, Grid *coarseGridIn, FieldProxy *workProxy);
    void PositiveUpdate(Grid *grid, FieldProxy *qProxy, FieldProxy *dqProxy);

    void SolveLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    void SolveLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);

    //! To determine the CFL number of the next iteration step.
    void DetermineCFLNumber(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *dqProxy1);


    //! To execute the forward sweep of the Matrix LU-SGS method.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: dqProxy is an array that stores the difference values between the conservative variables at n and n+1 time step.
    void SolveMatrixLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
    //! To execute the backward sweep of the Matix LU-SGS method.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: dqProxy is an array that stores the difference values between the conservative variables at n and n+1 time step.
    void SolveMatrixLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

    void SetGhostDQLUSGS(Grid *gridIn, RDouble4D &dq);

    //calculate the total area of all BC regions with the same BC name.
    void calculateBCArea(Grid *gridIn, StructBC *structBC);

    //calculate the mass flux ratio of all mass outlet BC regions with the same BC name.
    void CalculateMassFluxRatio();

    void ViscousAdiabaticWall(Grid *gridIn, StructBC *structBC);
    void ViscousIsotropicWall(Grid *gridIn, StructBC *structBC, int iWallBC);
    void SymmetryBC3D(Grid *gridIn, StructBC *structBC);
    void SingularAxisBC3D(Grid *gridIn, StructBC *structBC);
    void SingularAxisFullPlaneBC3D(Grid *gridIn, StructBC *structBC);
    void FarFieldRiemannInvariants(Grid *gridIn, StructBC *structBC, RDouble refGama);
    void InFlowBC3D(Grid *gridIn, StructBC *structBC);
    void OutFlowBC3D(Grid *gridIn, StructBC *structBC);
    void PressureInletBCRiemann(Grid *gridIn, StructBC *structBC);
    void PressureOutletBC(Grid *gridIn, StructBC *structBC);

    //! Get the primitiveVar of Inflow BC by Interrpolate from trajectory data.
    void SetBCInitialValuesByInterpolation(StructBC *structBC, RDouble *primitiveVar);

    //! Modify primitiveVar of Inflow BC by SpeedCoef.
    void SetBCInitialValuesBySpeedCoef(int nNumberOfSpeedStep, RDouble *primitiveVar);

    //! Mass flow in boundary condition.
    void MassFlowInletBC(Grid *gridIn, StructBC *structBC);

    //! Mass flow out boundary condition.
    void MassFlowOutletBC(Grid *gridIn, StructBC *structBC);

    void NozzleInletFlow(Grid *gridIn, StructBC *structBC);
    void NozzleOutFlow(Grid *gridIn, StructBC *structBC);

    void CornerPoint(Grid *gridIn);
    
    void FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value);
    void FillField(Grid *gridIn, FieldProxy *field1Proxy, FieldProxy *field2Proxy);
    
    FieldProxy * GetFieldProxy(Grid *gridIn, const string &field_name);
    FieldProxy * CreateFieldProxy(Grid *gridIn);
    FieldProxy * GetResidualProxy(Grid *gridIn);

    void ReconGradVolWeightedwithCorrection(Grid *gridIn, int nSurface);
    void GetGradientAtFace_CorrectionAtPhysicalBoundary(Grid *gridin, int iSurface);

    void StoreRhsByResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);
    void InitResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);
    void RecoverResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);
    
    FaceProxy * CreateFaceProxy(Grid *gridIn);
    GeomProxy * CreateGeomProxy(Grid *gridIn);
    
    StructIter * CreateInvStructIter(Grid *gridIn, int iSurface);
    
    bool GetInviscidFaceValue(Grid *gridIn, FaceProxy *faceProxy, StructIter *structIter, FieldProxy *qlProxy, FieldProxy *qrProxy, int nlen);
    bool LoadFaceFlux(Grid *gridIn, FaceProxy *faceProxy, StructIter *structIter, FieldProxy *fluxProxy, int nlen);

    //! Restart from perfect gas flow field for chemical reaction simulation.
    void ReadPrefectGasRestartH5(ActionKey *actkey);

    void DumpRestart(ActionKey *actkey);
    void ReadRestart(ActionKey *actkey);

    void DumpRestartH5(ActionKey *actkey);
    void ReadRestartH5(ActionKey *actkey);
    void ReadRestartFlowH5(ActionKey *actkey);

    void AirForceCoef(ActionKey *actkey);

    //! Set the temperatures using the variables reading from the previous file.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: nChemical denotes the flow type of the previous computation.
    //! @param[in ]: nTemperatureModel denotes the temperature type of the previous computation.
    //! @param[in ]: primTemperature denotes the temperatures of the previous computation.
    void SetFlowfieldTemperature(Grid *grid, int nChemical, int nTemperatureModel, RDouble4D *primTemperature);

    //! Set the primitive variables using the variables reading from the previous file.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: nChemical denotes the flow type of the previous computation.
    //! @param[in ]: nTemperatureModel denotes the temperature type of the previous computation.
    //! @param[in ]: numberOfSpecies denotes the species number of the previous computation.
    //! @param[in ]: speciesOrder denotes the species order of the previous computation.
    //! @param[in ]: primVariables denotes the primitive variables of the previous computation.
    void SetPrimitiveVariables(Grid *grid, int nChemical, int nTemperatureModel, int numberOfSpecies, Int1D *speciesOrder, RDouble4D *primVariables);

    //! Set the surface mass fractions using the variables reading from the previous file.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: nChemical denotes the flow type of the previous computation.
    //! @param[in ]: numberOfSpecies denotes the species number of the previous computation.
    //! @param[in ]: speciesOrder denotes the species order of the previous computation.
    //! @param[in ]: primMassFractions denotes the surface mass fractions of the previous computation.
    void SetSurfaceMassFractions(Grid *grid, int nChemical, int numberOfSpecies, Int1D *speciesOrder, RDouble3D *primMassFractions);

    //! Set the surface slip variables using the variables reading from the previous file.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[in ]: nChemical denotes the flow type of the previous computation.
    //! @param[in ]: numberOfSpecies denotes the species number of the previous computation.
    //! @param[in ]: speciesOrder denotes the species order of the previous computation.
    //! @param[in ]: primSlipVar denotes the surface slip variables of the previous computation.
    void SetSurfaceSlipVariables(Grid *grid, int nChemical, int numberOfSpecies, Int1D *speciesOrder, RDouble3D *primSlipVar);

private:
    void AirForceCoefParts(ActionKey *actkey, int partID, int Coordinate = GlobalCoordinate);
    void AirForceCoefBodies(ActionKey *actkey);
    
    //! Write the variables of Cp, Cf, Qw, y+ etc. on surface to the file.
    void CpDistriCoef(ActionKey *actkey);
    
    //! To write variables on surface such as p, rho, T etc. to the file. added by Li Peng on Mar.5, 2019.
    //! Add this function in order to verify the chemical codes and the related algorithm of computing thermal distribution.
    void ExportSurfaceVariables(ActionKey *actkey);

    //! Compute the surface heating change and obtain the maximum change of the surface heating ratio.
    void GetSurfaceHeatingChange(ActionKey *actkey, RDouble &localMaxHeatChange);

    //! Write the surface information to the file.
    void WriteSurfaceInfo(ActionKey *actkey);

    //! To obtain the derivatives of the designate primary variable whose sequence number in the array is marked with the valued of nIndex.
    void GetDerivatives(int nIndex, int iSurface, int ndim, int nlr, int i, int j, int k, int ni,int nj, int nk, RDouble4D &q,
                        RDouble4D &xvn, RDouble4D &yvn, RDouble4D &zvn, RDouble3D &vol, RDouble &dvdx, RDouble &dvdy, RDouble &dvdz);
    
    void GetResidual(ActionKey *actkey);
    void GetResidual(ActionKey *actkey, RDouble &localMaxRes);
    void GetResidual(Grid *gridIn, RDouble &localMaxRes);

    //! Obtain the CFL number.
    void ObtainCFLNumber(Grid *gridIn, RDouble &globalMinCFL, RDouble &globalMaxCFL);

    void VisualizationAverageFlow(ActionKey *actkey);

    void ComputeHeatFlux(ActionKey *actkey, RDouble *HeatFlux);

public:
    RDouble4D * CompMachNumber(Grid *gridIn);
    RDouble4D * CompMachNodeNumber(Grid* gridIn, Post_Visual *postVisualization);
    RDouble4D * ComputeMachNumberField(Grid *gridIn);

    void ComputeVorticitybyQCriteria(Grid *gridIn, RDouble3D *&vorticity_x, RDouble3D *&vorticity_y, RDouble3D *&vorticity_z, RDouble3D *&vorticityMagnitude,  RDouble3D *&strain_rate, RDouble3D *&Q_criteria);
    void ComputeCp(Grid *gridIn, RDouble4D *&cp);

    //! Compute the dimensional variables for density, pressure, temperature and the number density of electron.
    RDouble4D * ComputeDimensionalElectron(Grid *gridIn);

    RDouble4D * ComputeDimensionalVariables(Grid *gridIn);

    RDouble4D * ComputePrimitiveVariablesWithMoleFraction(Grid *gridIn);

    //! Compute the non-equilibrium number, such as the Damkohler number, Vibrational non-equilibrium number, etc.
    RDouble4D * ComputeNonequiliriumNumber(Grid *gridIn);

public:
    void VisualizationVariables(Grid *gridIn, FieldProxy *workProxy, const string &filename);

    void UploadInterfaceData  (ActionKey *actkey);

    void DownloadInterfaceData(ActionKey *actkey);

    void UploadOversetData  (ActionKey *actkey);

    void DownloadOversetData(ActionKey *actkey);

    GradientCellCenter *GetGradientCellCenterUVWT(Grid *gridIn, FieldProxy *q);
    GradientCellCenter *GetGradientCellCenterUVWT(Grid *gridIn);

    //! Prepare post visualization variables, store it in to grid's postVisualization database.
    LIB_EXPORT void ComputePostVisualVariables(Post_Visual * postVisualization);
    //! Prepare monitored probes variables for post-processing , store it into postProbesVar database.
    LIB_EXPORT void ComputePostProbesVariables(Post_Probes * postProbesVar);
    //! Create and init control parameters.
    LIB_EXPORT void InitControlParameters();

    //! Get controlParamters(a Para_NSSolver type pointer).
    LIB_EXPORT Param_NSSolverStruct *GetControlParameters() const;

    void Monitor(Grid *gridIn, RDouble4D *dq=NULL);
};

    bool iterijk(int &i, int &j, int &k, int ist, int ied, int jst, int jed, int kst, int ked);
}

