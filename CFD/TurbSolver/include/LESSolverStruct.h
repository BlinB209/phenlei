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
//! @file      LESSolverStruct.h
//! @brief     Large-eddy simulation solver for struct grid.
//! @author    Zhang Zipei.

#pragma once
#include "LESSolver.h"
#include "AleModel.h"

namespace PHSPACE
{
class Grid;
class StructBC;
class Param_LESSolverStruct;

class LESSolverStruct : public LESSolver
{
private:

public:
    LESSolverStruct();
    ~LESSolverStruct();

public:
    //! Create and initialize controlparameters.
    LIB_EXPORT void InitControlParameters();

    //! Get control paramters.
    LIB_EXPORT Param_LESSolverStruct * GetControlParameters();

    //! Judge if the flow field variables file already exists in LES solver.
    bool JudgeIfRestart();
    bool JudgeIfReadAverage();

    //! To write the flow data to restart file by HDF5 "fts" form for continue to simulation.
    //! @param[in ]: actkey contains the control information.
    void DumpRestartH5(ActionKey *actkey);

    //! To read the flow data from restart file by HDF5 "fts" form  for continue to simulation.
    //! @param[in ]: actkey contains the control information.
    void ReadRestartH5(ActionKey *actkey);

private:
    void AllocateGlobalVar(Grid *gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);

    void Init(Grid *gridIn);

    void InitFlowAsRestart();

    void ComputeCellLengthScale(Grid *gridIn);

    void ComputeWallFunction(Grid *gridIn);

    void WallFunctionofVanDriest(Grid *gridIn);
    void WallFunctionofDXB(Grid *gridIn);
    void WallFunctionofPiomelli(Grid *gridIn);

    void ComputeViscousCoeff(Grid *gridIn);

    void GetDependentVariablesforStructHighOrder(Grid *gridIn);

    void ComputeGradientCellCenter(Grid *gridIn);
    void ObtainGradientCellCenterHighOrder(Grid *gridIn);
    void GetUVWTproxy(Grid *gridIn);
    void GetUVWTfaceproxy(Grid *gridIn, int nsurf);
    void GetUVWTfaceproxy_HDCS(Grid *gridIn, int nsurf);
    void GetdqdkxietactaCellCenter(Grid *gridIn, int nsurf);
    void GetdqdkxietactaCellCenter_HDCS(Grid *gridIn, int nsurf);

    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);

    //! Interface: to get the turbulent viscosity for struct high-order solver.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void ObtainViscosity(Grid *gridIn);

    //! To set turbulent flow field variables's value at the corner point.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void CornerPoint(Grid *gridIn);

    //! To set the flow field value of the external boundary condition of the computational region. \n
    //! It looped all the internal boundary and was realized by call the specific boundary functions.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void Boundary(Grid *gridIn);
    void ObtainBoundaryValue(Grid *gridIn);
    void VisWall(Grid *gridIn, StructBC *structBC);
    void FarFieldBC(Grid *gridIn, StructBC *structBC);
    void SymmetryBC(Grid *gridIn, StructBC *structBC);
    void InFlowBC(Grid *gridIn, StructBC *structBC);
    void OutFlowBC(Grid *gridIn, StructBC *structBC);

    void Smagorinsky(Grid *gridIn);

    void DynamicSmagViscosity(Grid *gridIn);
    void DynamicSmagViscosityCompressible(Grid *gridIn);
    //! CHANG-YUE XU, LI-WEI CHEN AND XI-YUN LU, 2010, Large-eddy simulation of the compressible flow past a wavy cylinder.

    void DynamicSmagViscosityCompressibleFD(Grid *gridIn);

    void DynamicSmagViscosityCompressibleOld(Grid *gridIn);
    void DynamicSmagViscosityIncompressible(Grid *gridIn);
    void DynamicSmagViscosityIncompressibleOld(Grid *gridIn);
    void ComputeTestFilteredQ(Grid *gridIn);
    void ComputeTestFilteredQFD(Grid *gridIn);
    //void ComputeTestFilteredDensity(Grid *gridIn);
    void ComputeTestFilteredTemperature(Grid *gridIn);
    void ComputeTestFilteredTemperatureFD(Grid *gridIn);

    void ComputeStrainRateTensor(Grid *gridIn);

    void ComputeGradT(Grid *gridIn);
    void ComputeGradTFD(Grid *gridIn);

    void ComputeTestFilteredStrainRateTensor(Grid *gridIn);
    void ComputeTestFilteredStrainRateTensorFD(Grid *gridIn);

    void ComputeTestFilteredRhoUiUj(Grid *gridIn);
    void ComputeTestFilteredRhoUiUjFD(Grid *gridIn);

    void ComputeLeonardStress(Grid *gridIn);
    void ComputeModeledStress(Grid *gridIn);
    void ComputeTestFilteredAlphaIJ(Grid *gridIn);
    void ComputeTestFilteredAlphaIJFD(Grid *gridIn);
    void ComputeBetaIJ(Grid *gridIn);
    void ComputeAnisotropicConstant(Grid *gridIn);
    void ComputeIsotropicConstant(Grid *gridIn);
    void ComputePrandtlNumber(Grid *gridIn);
    void ComputeLeonardTemperatureTerm(Grid *gridIn);
    void ComputeModeledTemperatureTerm(Grid *gridIn);

    void ComputeTestFilteredRhoUT(Grid *gridIn);
    void ComputeTestFilteredRhoUTFD(Grid *gridIn);
    
    void ComputeTestFilteredGradT(Grid *gridIn);

    void ComputeTestFilteredRhoStrainRateMagnitudeGradT(Grid *gridIn);
    void ComputeTestFilteredRhoStrainRateMagnitudeGradTFD(Grid *gridIn);

    void WALE(Grid *grid);
    //! F. Nicoud and F. Ducros, 1999, Subgrid-Scale Stress Modelling Based on the Square of the Velocity Gradient Tensor.
    //! Christian Wollblad et.al., 2006, Large Eddy Simulation of Transonic Flow with Shock Wave/Turbulence Boundary Layer Interaction.

    void Sigma(Grid *gridIn);

    void ComputeDifferentialSigma(RDouble ** velocityGradTensor, RDouble &differentialSigma);

    void filter(const int & nn, RFloat1D * fu1, RFloat1D * fu2, RFloat1D * fu3
                ,RFloat1D * fuifuj11, RFloat1D * fuifuj22, RFloat1D * fuifuj33, RFloat1D * fuifuj12, RFloat1D * fuifuj13, RFloat1D * fuifuj23
                ,RFloat1D * fs11    , RFloat1D * fs22    , RFloat1D * fs33    , RFloat1D * fs12    , RFloat1D * fs13    , RFloat1D * fs23    
                ,RFloat1D * fsfs11  , RFloat1D * fsfs22  , RFloat1D * fsfs33  , RFloat1D * fsfs12  , RFloat1D * fsfs13  , RFloat1D * fsfs23);
    void TestFilter(const int &nn,RFloat1D *fu);
    void TestFilter(const int &n0, const int &nn, RDouble1D *&tf, RDouble1D *&f);
    void TestFilter(const int &n0, const int &nn, RDouble *tf, RDouble *f);

    //! To obtain the grid scale using in sub-grid model.
    void GetDeltaBar(const int iDeltaBar, const RDouble deltai, const RDouble deltaj, const RDouble deltak, RDouble &deltaBar);

    //! To calculate the non-isotriopic fix of filter width of Scotti&Meneveau (scotti1993, scotti1997, physics of fluid).
    void GetFilterWidthOfScotti(const RDouble deltai, const RDouble deltaj, const RDouble deltak, RDouble &deltaBar);

public:
    void UploadInterfaceValue(ActionKey *actkey);
    void UploadInterfaceData(ActionKey *actkey);

    void DownloadInterfaceValue(ActionKey *actkey);
    void DownloadInterfaceData(ActionKey *actkey);

    void RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation);

private:
    FieldProxy *UVWTproxy;
    FieldProxy *UVWTfaceIproxy;
    FieldProxy *UVWTfaceJproxy;
    FieldProxy *UVWTfaceKproxy;
};

}