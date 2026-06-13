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
//! @file      LES.h
//! @brief     LES Solver.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "Complex.h"
#include "SpecDiffHybGrid.h"
//#include "SpecGrid.h"
#include "CompactDifferenceFirstDerivative.h"
//#include "Param_SpecSolver.h"

using namespace std;

namespace PHSPACE
{
class LESHIT
{
public:
    LESHIT();
    ~LESHIT();

public:
    Complex4D *complexDiffTau;

protected:
    RDouble1D *gkx, *gky, *gkz;
    Complex4D *complexVarFilter, *testFilterGradUij, *complexTauij;
    RDouble4D *realVarFilter, *realTauij;

private:
    RDouble delta,dynamicCoef;

public:
    void InitLESData(SpecGrid *GridData);

protected:
    void ComputeFilterFunction(const RDouble, RDouble, RDouble &);

private:
    void AllocLESData(SpecGrid *GridData);

    void GetDelta(SpecGrid *GridData);

    void SetTestFilter(SpecGrid *GridData);  
};

class LES : public LESHIT
{
public:
    LES();
    ~LES();

private:
    RDouble1D *delta;
    RDouble1D *wallDistance, *wallDistancePlus;
    RDouble1D *vanDriestFunc;
    RDouble1D *dynamicCoef;

public:
    void InitLESData(SpecDiffHybGrid *GridData);

    void Smagorinsky(SpecDiffHybGrid *, RDouble1D *, RDouble4D *);

    void DynamicSmagorinsky(SpecDiffHybGrid *, CompactDifferenceFirstDerivative *, RDouble1D *, RDouble3D *, RDouble4D *, Complex3D *, Complex3D *, Complex3D *);

    void Sigma(SpecDiffHybGrid *, RDouble1D *, RDouble4D *);

    //LIB_EXPORT Param_SpecSolver * GetControlParameters() const;

    //LIB_EXPORT void InitControlParameters();

private:
    void AllocLESData(SpecDiffHybGrid *);

    void GetWallDistance(SpecDiffHybGrid *);

    void GetDelta(SpecDiffHybGrid *);

    void SetTestFilter(SpecDiffHybGrid *);  

    void FilterXYInSpectralSpace(PHComplex *, const int, SpecDiffHybGrid *);

    void GetDiffTau(SpecDiffHybGrid *, CompactDifferenceFirstDerivative *);

    void ComputeDifferentialSigma(RDouble ** velocityGradTensor, RDouble &differentialSigma);

    //void ComputeFilterFunction(const RDouble, RDouble, RDouble &);

};

}



