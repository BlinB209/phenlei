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
//! @file      IncomGas.h
//! @brief     Explain this file briefly.
//! @author    Ming Pingjian, Sun Huawen, Huang Jin.

#pragma once
#include "Precision.h"
#include "PerfectGas.h"
#include "Geo_UnstructGrid.h"
#include <string>
using namespace std;

namespace PHSPACE
{
namespace GAS_SPACE
{
typedef struct {
    RDouble molarMass;
    RDouble rho;
    RDouble mu;
    RDouble cp;
    RDouble k;
} propertySp;

typedef enum
{
    rhoConstant = 0,
    rhoIncompressibleIdealGas,
    rhoIdealGas,
    rhoVolumnWeightedMixingLaw
}RhoProperty;

typedef enum
{
    muConstant = 0,
    muSutherLand,
    muIdealGasMixingLaw,
    muMassWeightedMixingLaw,
    muPiecewiseLinear,
    muPiecewisePolynomial,
    muPolynomial,
    muPowerLaw
}MuProperty;

typedef enum
{
    cpConstant = 0,
    cpMixing,
    cpPiecewiseLinear,
    cpPiecewisePolynomial,
    cpPolynomial
}CpProperty;

typedef enum
{
    kConstant = 0,
    kMassWeightedMixingLaw,
    kIdealGasMixingLaw,
    kPiecewiseLinear,
    kPiecewisePolynomial,
    kPolynomial
}KProperty;

typedef enum
{
    massDiffConstantDiluteApprox = 0,
    massDiffDiluteApprox,
    massDiffMulticomponent,
    massDiffUnityLewisNumber
}MassDiffProperty;

class IncomGas :public PerfectGas
{
public:
    IncomGas();
    ~IncomGas();

public:
    void InitCommonParameter();

    void InitParameterForMu();
    void InitParameterForK();
    void InitParameterForMassdiff();

    void SetGasName(string speciesNameIncom);
    void SetGasMolarMass(string gasMolarMass);
    void SetGasRho(string gasRho);
    void SetGasMu(string gasMu);
    void SetGasCp(string gasCp);
    void SetGasK(string gasK);

    string * GetGasName() { return gasNameList; }
    RDouble * GetGasMolarMass() { return gasMolarMassList; }
    RDouble * GetGasRho() { return gasRhoList; }
    RDouble * GetGasMu() { return gasMuList; }
    RDouble * GetGasCp() { return gasCpList; }
    RDouble * GetGasK() { return gasKList; }

public:
    void UpdateRg(UnstructGrid *grid, RDouble *Rg);
    void UpdateRho(UnstructGrid *grid, RDouble *rho);
    void UpdateMu(UnstructGrid *grid, RDouble *mu);
    void UpdateCp(UnstructGrid *grid, RDouble *cp);
    void UpdateK(UnstructGrid *grid, RDouble *k);
    void UpdateMassDiff(UnstructGrid *grid, int iGas, RDouble *massdiff);

protected:
    void ComputePhiIJ();
    void ComputeAmountOfSubstance(UnstructGrid *grid, RDouble **amountOfSubstance, RDouble **amountFracOfSubstance);
    void SetPropertyLib();

private:
    vector<string> componentLib;
    map<string, propertySp> propertyLib;
};

}
}