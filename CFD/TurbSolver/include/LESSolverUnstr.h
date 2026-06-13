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
//! @file      LESSolverUnstr.h
//! @brief     Large-eddy simulation solver for unstruct grid.
//! @author    Zhang Zipei.

#pragma once
#include "LESSolver.h"
#include "AleModel.h"

namespace PHSPACE
{
class Param_LESSolverUnstruct;

class LESSolverUnstr : public LESSolver
{

public:
    LESSolverUnstr();
    ~LESSolverUnstr();

    void AllocateGlobalVar(Grid *gridIn);

    void DeAllocateGlobalVar(Grid *gridIn);

    LIB_EXPORT void InitControlParameters();

    LIB_EXPORT Param_LESSolverUnstruct *GetControlParameters();

private:
    void Init(Grid *gridIn);

    void InitFlowAsRestart();

    void ComputeCellLengthScale(Grid *gridIn);

    void ComputeGradient(Grid *gridIn);

    void GetGradientField(Grid *gridIn);

    void ComputeViscousCoeff(Grid *gridIn);

    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);

    void ObtainViscosity(Grid *gridIn);

    void ObtainBoundaryValue(Grid *gridIn);

    void ComputeStrainRateTensor(Grid *gridIn);

    void ComputeWallFunction(Grid *gridIn);

    void WallFunctionofVanDriest(Grid *gridIn);

    void Boundary(Grid *gridIn);

    void Smagorinsky(Grid *gridIn);

    void WALE(Grid *gridIn);

    void Sigma(Grid *gridIn);

    void ComputeDifferentialSigma(RDouble ** velocityGradTensor, RDouble &differentialSigma);

public:
    void UploadInterfaceValue(ActionKey *actkey);
    void UploadInterfaceData(ActionKey *actkey);

    void DownloadInterfaceValue(ActionKey *actkey);
    void DownloadInterfaceData(ActionKey *actkey);

    FieldProxy * GetFieldProxy(Grid *gridIn, const string &field_name);
};

}