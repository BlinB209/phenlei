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
//! @file      NSSolverUnstruct.h
//! @brief     NS solver for unstructured grid.
//! @author    Bell, He Xin, Zhang Laiping, Dr. Wang, Zhang Yang, He Kun, 
//!            Wan Yunbo, Xu Gang, Meng Liyuan, Zhang Yaobing.

#pragma once
#include "NSSolver.h"
#include "PHMatrix.h"
#include "Limiter.h"

namespace PHSPACE
{
class GeomProxy;
class Gradient;
class Param_NSSolverUnstruct;

class NSSolverUnstruct : public NSSolver
{
private:
    Gradient * gradientNSField;
    Gradient * gradientTemperature;
    int systemgridtype;
    Limiter *limiterOnlyForMixGrid;
    FaceProxy *InviscidfaceProxyOnlyForMixGrid;
    FaceProxy *ViscousfaceProxyOnlyForMixGrid;
public:
    NSSolverUnstruct();
    ~NSSolverUnstruct();
public:
    void ReadParameter();

    void AllocateGlobalVar(Grid *gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);

#ifdef USE_GMRESSOLVER
    // GMRESCSR
    void FindNeighborOfRoot(int root, int& neighbor, int nTotalFace, int nBoundFace, int* leftCellofFace, int* rightCellofFace, vector<int>& AI, vector<int>& AJ, int counter);
    void AllocateJacobianMatrix4GMRES(Grid *gridIn);
    void AllocateJacobianMatrix4GMRES_Neighbor(Grid *gridIn);
#endif

    //! Judge if the flow field variables file already exists in NS solver.
    bool JudgeIfRestart();
    bool JudgeIfInterpolate();
    bool JudgeIfReadAverage();

    void DumpRestart(ActionKey *actkey);
    void ReadRestart(ActionKey *actkey);

    void DumpRestartH5(ActionKey *actkey);
    void ReadRestartH5(ActionKey *actkey);
    void InterpolateFromRestartH5(ActionKey *actkey);

    void InitFlowAsRestart();
    void LoadQ(Grid *gridIn, FieldProxy *qProxy);
    void RestrictDefect(Grid *fgrid_in, Grid *cgrid_in);
    void PutCorrectionBack(Grid *gridIn, FieldProxy *qProxy);

    //! Prolongation operator, from coarse mesh to fine mesh.
    void CorrectFineGrid(Grid *fgrid_in, Grid *cgrid_in);

    //! Prolongation by the current coarse cell directly.
    void CorrectFineGridZeroOrder(Grid *fgrid_in, Grid *cgrid_in);

    //! Prolongation by the current and the neighbor coarse cells.
    void CorrectFineGridHighOrder(Grid *fgrid_in, Grid *cgrid_in);

    void InterpolatFineGrid(Grid *fgrid_in, Grid *cgrid_in);


    void LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy);   

    void Boundary(Grid *gridIn);
    void CornerPoint(Grid *gridIn);
    void CalculateBoundaryData();
    void ComputeMassFlow(Grid *gridIn);
    void ComputeMassInFlow(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    void ComputeMassOutFlow(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    void ComputeGamaAndTemperature(Grid *gridIn);

    //! Compute the node value in the flow.
    //! @param[in] gridIn         the mesh in the computing.
    void ComputeNodeValue(Grid *gridIn);

    //! Modify the node value for the interpoint.
    //! @param[in] gridIn         the mesh in the computing.
    void ModifyNodeValue(Grid *gridIn);

    void CompNodeVarForGGNodeLaplacian(Grid *gridIn);

    void GetQlQr  (Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void GetGamaLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void GetPressureFactorLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void GetTimeCoefficientLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void GetPreconCoefficientLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void ReConstructFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd);
    void ZeroResiduals(Grid *gridIn);
    void SpectrumRadiusInviscid(Grid *grid);
    void SpectrumRadiusViscous(Grid *grid);
    void SpectrumRadiusChemical(Grid *grid);
    void ConservativePreconMatrix(Grid *gridIn, int sign);
    void Diagonal(Grid *gridIn);
    void ComputeMinTimeStep(Grid *grid, RDouble & minDt, RDouble & maxDt);
    void ReduceMaxTimeStep(Grid *grid, RDouble globalMinDt);
    void LocalTimeStep(Grid *grid);
    void LimitCFL(Grid *gridIn, RDouble *CFLCell);

    void GlobalTimeStep(Grid *grid);
    void LocalGlobalTimeStep(Grid *grid);
    void RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coef);
    void LoadFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);

    void ZeroResidualOfSpecialCells(Grid *gridIn);

    void UpdateQlQrOnlyforMixGrid(Grid *grid);
    void UpdateInviscidfluxOnlyforMixGrid(Grid *grid);
    void UpdateViscousfluxOnlyforMixGrid(Grid *grid);
    void LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn);
    void RestrictAllQ(Grid *fgrid, Grid *cgrid);
    void ComputeViscousCoeff(Grid *gridIn);
    void ComputeViscousCoefficient(Grid *gridIn);
    void ComputeViscousCoefficientWithoutChemical(Grid *gridIn);

    void BoundaryQlQrFix(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void InitCGrid(Grid *fgrid_in, Grid *cgrid_in);
    RDouble * CompSoundSpeed(Grid *gridIn);
    RDouble * CompMachNumber(Grid *gridIn);
    RDouble ** ComputeDimensionalVariables(Grid *gridIn);
    RDouble ** ComputePrimitiveVariablesWithMoleFraction(Grid *gridIn);

    RDouble * ComputeVorticitybyQCriteria(Grid *gridIn);
    void ComputeVorticitybyQCriteria(Grid *gridIn, RDouble * vorticity_x, RDouble * vorticity_y, RDouble * vorticity_z, RDouble * vorticityMagnitude,  RDouble * strain_rate,  RDouble * Q_criteria);
    RDouble * ComputeCp(Grid *gridIn);
    RDouble * ComputeCpAverage(Grid *gridIn);

    void ComputeGradient(Grid *grid);
    void ComputeLimiter(Grid *grid);
    void ComputelimiterOnlyForMixGrid(Grid *grid);
    void InitiallimiterOnlyForMixGrid(Grid *grid);
    void InitialInviscidfaceproxyOnlyForMixGrid(Grid *grid);
    void InitialViscousfaceproxyOnlyForMixGrid(Grid *grid);

    void FillField(Grid *gridIn, FieldProxy *field_proxy, RDouble value);
    void FillField(Grid *gridIn, FieldProxy *fieldTarget, FieldProxy *fieldSource);
    FieldProxy * CreateFieldProxy(Grid *gridIn);
    FieldProxy * GetFieldProxy(Grid *gridIn, const string &field_name);
    FieldProxy * GetResidualProxy(Grid *gridIn);

    void RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy);
    void InitResidual(Grid *gridIn, FieldProxy *rhsProxy);
    void StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy);

    FaceProxy     * CreateFaceProxy    (Grid *gridIn);
    GeomProxy     * CreateGeomProxy    (Grid *gridIn);
    NSFaceValue   * CreateNSFaceValue  (Grid *gridIn);

    void FillGeomProxy(Grid *gridIn, GeomProxy *geom_proxy, int localStart, int localEnd);

    void CalPressureFactor(Grid *gridIn);
    void GetHeatTransferCoeff(Grid *gridIn);

    void FreeGradientProxy(Grid *grid);

    void DumpLeakageInformation(Grid *grid);

    void InviscidFlux  (Grid *gridIn);
    void ViscousFlux   (Grid *gridIn);
    void ChemicalSource(Grid *gridIn);
    void gravitySource(Grid *gridIn);
    void PorousMediumSource(Grid *gridIn);
    void RotatingSource(Grid* gridIn);
    void DualTimeSource(Grid *gridIn);
    void StoreHistoricalResidual(Grid *gridIn);
    void UpdateUnsteadyFlow(Grid *gridIn);
    void ComputePreconditionCoefficient(Grid *grid);
    void ComputePreconditionCoefficientUnsteady(Grid *gridIn);
    RDouble UnsteadyConvergence(Grid *gridIn);

    //Compute the distance between the center of the grid cell on the first layer and the center of the surface cell.
    void ComputeFirstLayerGridHeight(Grid *gridIn);

    void Turb_Sengy(Grid *gridIn);
    void GetInviscidFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd);
    void GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
 
    void ComputeVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void ModifyVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);

    #ifdef USE_GMRESSOLVER
    // GMRESVis
    void ComputeVisflux_FD(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void ComputeVisflux_GMRES(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    // GMRESVis GMRESVisFD
    void Cal_GMRES_Visflux(RDouble* fL, RDouble* fR, RDouble* primitiveVariableFace, RDouble* dqdxL,
                           RDouble* dqdyL, RDouble* dqdzL, RDouble* dqdxR, RDouble* dqdyR, RDouble* dqdzR, int iFace,
                           int jFace, RDouble* flux, Grid *gridIn, FaceProxy *faceProxy);
    // GMRESAD
    template <typename T>
    void Cal_GMRES_Visflux_AD(T* fL, T* fR, T* dqdxL,
                           T* dqdyL, T* dqdzL, T* dqdxR, T* dqdyR, T* dqdzR, int iFace,
                           int jFace, T* flux, Grid *gridIn, FaceProxy *faceProxy, Param_NSSolverUnstruct *parameters);

    // GMRESCoupled
    template <typename T>
    void Cal_GMRES_Visflux_AD_Coupled(T* fL, T* fR, T fLturb, T fRturb, T* dqdxL,
                                    T* dqdyL, T* dqdzL, T* dqdxR, T* dqdyR, T* dqdzR, int iFace,
                                    int jFace, T* flux, Grid *gridIn, FaceProxy *faceProxy, Param_NSSolverUnstruct *parameters);


    //! GMRES solver -- a linear system solver
    //! This method can be referenced "Sen Zhang, Boqian Wang, et al. Implementation of a Newton-Krylov Algorithm in the Open-source Solver PHengLEI[C]."
    //! GPPS Hong Kong , October 17-19, 2023.
    void GMRESSolver(Grid *gridIn, FieldProxy *dqProxy);
    void GMRESSolverRow(Grid *gridIn, FieldProxy *dqProxy);
    void GMRESSolverSingle(Grid *gridIn, FieldProxy *dqProxy);
#endif

    void ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void ComputeInviscidFlux(Grid *grid, FaceProxy *faceProxy, int localStart, int localEnd);
    void SolutionFix(Grid *gridIn, RDouble *prim, int iCell);

    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);
    void SetGhostDQLUSGS(Grid *gridIn, RDouble **dq, RDouble **q);
    void SetSymmBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q);
    void SetWallBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q);
    void SetExtraBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q);
    
    void SwapNeighborsDQNonblocking(Grid *gridIn, FieldProxy *dqProxy);
    void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    RDouble GetRad(RDouble *prim, RDouble *nxyz, RDouble gama);
    void SolveLUSGSForward(Grid *grid, FieldProxy * dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    void SolveLUSGSBackward(Grid *grid, FieldProxy * dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);


    void AirForceCoef(ActionKey *actkey);
    void AirForceCoefParts(ActionKey *actkey, int partID, int Coordinate = GlobalCoordinate);
    void AirForceCoefBodies(ActionKey *actkey);
    void CpDistriCoef(ActionKey *actkey);
    void GetResidual(ActionKey *actkey);

    //! Dump out typical standard case: turbulent flat plate.
    void Turbulence_Flat_Plate_Output_NASA(Grid *gridIn);
    void CreateNodeVariable(Grid *gridIn, RDouble **nodeVariable);

    void ComputeResidual(ActionKey *actkey);
    void VisualizationAverageFlow(ActionKey *actkey);

    void ComputeAverageVisualNodeField(Grid *gridIn, RDouble **qn);
    void BoundaryVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot);

    void BoundaryVisualization(Grid *gridIn, std::ostringstream &oss, vector<string> &title, RDouble **qn, int nvarplot);
    void SliceVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot);

    void SpectrumRadius(Grid *gridIn);

    void RestrictDQ(RDouble *dqtry, int nm, RDouble lim_max);
    void ModifyDensity(RDouble *dq, int nNSEquation,int nLaminar,int nEquation,RDouble &densityRelax);

    void UploadInterfaceData(ActionKey *actkey);
    void DownloadInterfaceData(ActionKey *actkey);

    //! Upload the interpoint variables for send.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void UploadInterpointData(ActionKey *actkey);

    //! Download the interpoint variables from the receive data.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void DownloadInterpointData(ActionKey *actkey);

    void CommunicationInterpointWeight(ActionKey *actkey);
    void DownloadInterpointWeight(ActionKey *actkey);

    void UploadOversetData  (ActionKey *actkey);
    void DownloadOversetData(ActionKey *actkey);

    //Gradient * GetGradientNSField(Grid *gridIn);
    //Gradient * GetGradientTemperature(Grid *gridIn);

    void GetGradientField(Grid *gridIn);

    //! Prepare post visualization variables, store it in to postVisualization database.
    LIB_EXPORT void ComputePostVisualVariables(Post_Visual * postVisualization);
    //! Prepare monitored probes variables for post-processing , store it into postProbesVar database.
    LIB_EXPORT void ComputePostProbesVariables(Post_Probes * postProbesVar);
    //! Create and init ControlParameters
    LIB_EXPORT void InitControlParameters();

    //! Get controlParamters(a Para_NSSolver type pointer)
    LIB_EXPORT Param_NSSolverUnstruct * GetControlParameters() const;

    void CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn);
    void DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn);

    void CompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    //! Out flow boundary condition for one Region of a bounday.
    //! A BC Region Contains lots of faces with the same boundary name.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void OutflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    //! Overlapping boundary : The flow field variables are interpolated from the interpolating cell to the invalid cell 
    void ReSetOversetBoundary(Grid* gridIn);

    //! Symmetry boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void SymmetryBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

#ifdef USE_GMRESSOLVER
    // GMRESBoundary
    //! GMRES Symmetry boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    // GMRESAD GMRESBoundary perturbation about the primitive variables with AD
    void GMRES_SymmetryBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstr);

    // GMRESPV GMRESBoundary perturbation about the primitive variables with FD
    void GMRES_SymmetryBCRegion_FD(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    // GMRESBoundary perturbation about the conservative variables
    void GMRES_SymmetryBCRegion_CV(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    void Cal_GMRES_SymmetryBCRegion(Grid *gridIn, RDouble *PrimL, RDouble *PrimR, const int &iFace);  

     // GMRESVis GMRESBoundary GMRESAD
    //! Viscous Adiabatic wall boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    // GMRESBoundary perturbation about the primitive variables with AD
    void GMRES_ViscousAdiabaticWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    // GMRESBoundary perturbation about the primitive variables with FD
    void GMRES_ViscousAdiabaticWallBCRegion_FD(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    void Cal_GMRES_ViscousAdiabaticWallBCRegion(Grid *gridIn, RDouble *PrimL, RDouble *PrimR, int iFace, RDouble uwall, RDouble vwall, RDouble wwall);

    //! Viscous Isotropic wall boundary condition for one Region of a bounday. 
    //! for GMRES  GMRES3D
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void GMRES_ViscousIsotropicWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    // GMRESBoundary
    //! GMRES Farfield boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    // GMRESBoundary GMRESAD
    void GMRES_FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstrct);

    // GMRESPV GMRESBoundary
    void GMRES_FarfieldBCRegion_FD(Grid *gridIn, UnstructBC *bcRegionUnstruct); 
    
    // GMRESBoundary using the conservative variables
    void GMRES_FarfieldBCRegion_CV(Grid *gridIn, UnstructBC *bcRegionUnstruct); 

    // GMRESBoundary
    //! Calculate GMRES Farfield boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    template <typename T>
    void Cal_GMRES_FarfieldBCRegion(T* qL,T* qR,int le, int iFace, Grid *gridIn, UnstructBC *bcRegionUnstruct);
#endif 

    //! Viscous Adiabatic wall boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void ViscousAdiabaticWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Viscous Isotropic wall boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void ViscousIsotropicWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Farfield boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Inflow boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);   
    
    //! PressureOutlet boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void PressureOutletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);   
    
    //! PressureInlet boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void PressureInletBCRiemannRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct); 
    
    //! calculate the total area of all BC regions with the same BC name.
    void calculateBCArea(Grid *gridIn, UnstructBC *bcRegionUnstruct); 

    //! calculate the mass flux ratio of all mass outlet BC regions with the same BC name.
    void CalculateMassFluxRatio();

    //! mass flow in boundary condition.
    void MassFlowInletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct); 
    
    //! mass flow out boundary condition.
    void MassFlowOutletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    //! mixing plane for multi row turbomachinery.
    void MixingPlaneInterface(Grid *gridIn);

private:
    void InviscidFluxWrap(Grid *gridIn, FaceProxy *faceProxy, INVScheme inv_scheme, int localStart, int localEnd);
    void RegisterCFDSolverInterfaceField();
    //! Register the interpoint field in the NS unstructured solver.
    void RegisterCFDSolverInterpointField();
    void RegisterOversetField();

    void RotateVelocity(Grid *grid, RDouble **q);
    void RotateVectorFromInterface(Grid *grid, const int &neighborZoneIndex, const int &neqn);

    void InitMixingPlane(RDouble ***MixingPlaneVar,int Dim1,int Dim2, int Dim3,RDouble Value);

    void AverageMixingPlane(Grid *grid);
    void MixingPlaneDataTransfer(Grid *grid, Grid *NeighborGrid);
    void NonReflective(Grid *grid, Grid *NeighborGrid);
    void SetMixingPlaneData(Grid *grid);

    RDouble P0Wall(RDouble p1, RDouble gama, RDouble mach);
    RDouble R0Wall(RDouble r1, RDouble gama, RDouble mach);
    Limiter * CreateLimiter(Grid *gridIn);
};

}

