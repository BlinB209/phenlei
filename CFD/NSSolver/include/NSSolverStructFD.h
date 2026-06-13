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
//! @file      NSSolverStructFD.h
//! @brief     NS solver for struct grid with finite difference method.
//! @author    Min Yaobing, Ma Yankai.

//! CHEN J Q, MA Y K, MIN Y B, et al. Design and development of homogeneous hybrid 
//! solvers on National Numerical Windtunnel(NNW) PHengLEI[J]. Acta Aerodynamica Sinica,
//! 2020, 38(6): 1103-1110.

#pragma once
#include "NSSolverStruct.h"

namespace PHSPACE
{
class NSSolverStructFD : public NSSolverStruct
{
public:
    NSSolverStructFD();
    ~NSSolverStructFD();
private:
    void AllocateGlobalVar  (Grid *gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);
    void GetSurfaceCellNumber(Grid *gridIn, int &nWall, int &nCell);
    void ComputeFirstLayerGridHeight(Grid *gridIn, RDouble2D *firstLayerHeight);
    
    void Primitive2Conservative(RDouble prim[5], RDouble q[5]);
    void Conservative2Primitive(RDouble q[5], RDouble prim[5]);
    
    void InFlowBC3D(Grid *gridIn, StructBC *structBC);
    void OutFlowBC3D(Grid *gridIn, StructBC *structBC);
    void SymmetryBC3D(Grid *gridIn, StructBC *structBC);
    void ViscousAdiabaticWall(Grid *gridIn, StructBC *structBC);
    void ViscousIsotropicWall(Grid *gridIn, StructBC *structBC);
    void FarFieldRiemannInvariants(Grid *gridIn, StructBC *structBC);
    void CornerPoint(Grid *gridIn);
    
    void GetIndependentVariablesforStructHighOrder(Grid *gridIn);
    void GetDependentVariablesforStructHighOrder(Grid *gridIn);
    void ObtainBoundaryValue(Grid *gridIn);
    void ObtainGamaAndTemperature(Grid *gridIn);
    void ObtainPressureFactor(Grid *gridIn);
    void ObtainViscousCoef(Grid *gridIn);
    void ObtainGradientCellCenter(Grid *gridIn);
    
    void Boundary(Grid *gridIn);
    void ComputeGamaAndTemperature(Grid *gridIn);
    void CalPressureFactor(Grid *gridIn);
    void GetHeatTransferCoeff(Grid *gridIn);
    void ComputeViscousCoeff(Grid *gridIn);
    void ComputeViscousCoefficientWithoutChemical(Grid *gridIn);
    
    void ZeroResiduals(Grid *grid);
    void InitFlowAsRestart();
    void InitLocalMemory(Grid *gridIn);
    void SetWeightedInterpolationCoef(Grid *gridIn);
    
    int total_negative;
    string solvername;//!WCNS,HDCS
    double EPSILON;//! epsilon within weighted interpolation
    string fluxname;
    FieldProxy *CellfluxProxy;//! only used for HDCS

    RDouble6D *WeightedInterpolationCoef;
    
    string gradient_method;
    FieldProxy *UVWTproxy;
    FieldProxy *UVWTfaceIproxy;
    FieldProxy *UVWTfaceJproxy;
    FieldProxy *UVWTfaceKproxy;
    FieldProxy *dqdxFaceProxy;
    FieldProxy *dqdyFaceProxy;
    FieldProxy *dqdzFaceProxy;
    
    void InviscidFlux(Grid *grid);
    
    void GetFaceVariable_WCNS(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf);
    void GetFaceVariable_HDCS_firstlayerExplicit(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf);
    void GetFaceVariable_HDCS_otherlayerImplicit(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf);
    void GetFaceVariable_HDCS_otherlayerImplicit_nsurf1(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy);
    void GetFaceVariable_HDCS_otherlayerImplicit_nsurf2(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy);
    void GetFaceVariable_HDCS_otherlayerImplicit_nsurf3(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy);
    void SolveTridiagonalEquations(int Nod, RDouble BL, RDouble BR, RDouble *F);

    RDouble LinearExplicitInterpolation5th(RDouble prim[5]);
    RDouble LinearImplicitInterpolation7th(RDouble prim[5]);
    void GetFaceVariableCorrectionAtPhysicalBoundary(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf); 
    
    void InterpolationRUVWP5th(RDouble Coef[30], RDouble PrimPlate[5][5], RDouble Primhalf[5]);
    void InterpolationRUVWP3rd(RDouble Prim[5][5], RDouble Primhalf[5]);

    RDouble BasicWeightedInterpolation5th(RDouble Coef[30], RDouble QmPlate[5]);
    RDouble BasicInterpolation3rd(RDouble Qp[5], int m);
    RDouble Limiter_ThirdSmooth(const RDouble &x, const RDouble &y);
    
    void GetFaceInviscidFlux(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, int nsurf);
    void FluxSteger(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, int nsurf);
    void FluxRoe(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, int nsurf);
    void FluxDifference     (Grid *gridIn, FieldProxy *fluxProxy, int nsurf);
    void FluxDifference_WCNS(Grid *gridIn, FieldProxy *fluxProxy, int nsurf);
    void FluxDifference_HDCS(Grid *gridIn, FieldProxy *fluxProxy, int nsurf);
    void GetCellInviscidFlux(Grid *gridIn, int nsurf);
    void GetCellViscousFlux(Grid *gridIn, int nsurf);
    
    void GetUVWTproxy(Grid *gridIn);
    void GetUVWTfaceproxy(Grid *gridIn, int nsurf);
    void GetUVWTfaceproxy_HDCS(Grid *gridIn, int nsurf);
    void GetdqdkxietactaCellCenter(Grid *gridIn, int nsurf);
    void GetdqdkxietactaCellCenter_HDCS(Grid *gridIn, int nsurf);
    void ComputeGradientCellCenter(Grid *gridIn);
    
    void ViscousFlux(Grid *gridIn);
    void GetFaceViscousFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    void GetGradientAtFace(Grid *gridIn, int nSurface);
    
    void GetGradientAtFaceconservation(Grid *gridIn, int nSurface);
    void GetGradientAtFaceconservation_HDCS(Grid *gridIn, int nSurface);
    void GetGradientAtFacechainrule(Grid *gridIn, int nSurface);
    void GetGradientAtFaceCorrectionAtPhysicalBoundary(Grid *gridin, int iSurface);

    void SpectrumRadius(Grid *grid);
    void SpectrumRadiusInviscid(Grid *gridIn);
    void SpectrumRadiusViscous(Grid *gridIn);
    void TimeStep(Grid *gridIn);
    void ReduceMaxTimeStep(Grid *gridIn, RDouble globalMinDt);
    void Diagonal(Grid *gridIn);
    void SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    void SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    void SolveLUSGSbyJacobiIteration(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *dqtmpProxy, RDouble &sweepNormal);
    void LUSGSInitializationStructHighOrder(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ);
    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);

    void GetMXDQ(RDouble *prim, RDouble KX, RDouble KY, RDouble KZ, RDouble *dq, RDouble *f, RDouble radius, int ipn);
    void SolutionFix(FieldProxy *qProxy, RDouble prim[5], int i, int j, int k);

    void UpdateQlQrOnlyforMixGrid(Grid *gridIn);
    void UpdateInviscidfluxOnlyforMixGrid(Grid *gridIn);
    void UpdateViscousfluxOnlyforMixGrid(Grid *gridIn);
    void LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn);

    void GetFaceViscousFluxLES(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
};

}
