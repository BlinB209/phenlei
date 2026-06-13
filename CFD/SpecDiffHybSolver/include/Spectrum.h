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
//! @file      Spectrum.h
//! @brief     Spectrum data for SpecDiffHybSolver.
//! @author    ZhangZipei.

#pragma once
//#include "TypeDefine.h"
//#include "Complex.h"
#include "LIB_Macro.h"
#include "Precision.h"
//#include "Param_Spectrum.h"
//#include "SpecSolver.h"

using namespace std;

namespace PHSPACE
{

class Spectrum
{
public:
    Spectrum();
    ~Spectrum();

private:
    int spectrumType;

public:
    RDouble *k_CBC, *eK_CBC, *logK_CBC, *logEK_CBC, *d2logEKdK2_CBC;

    int nPoints;

private:
    void AllocSpectrumData(int nP);

    void InitCBCSpectrum();

    void ReadCBCSpectrum();

    void WriteCBCSpectrum();

    void Spline(RDouble *x, RDouble *y, int n, RDouble yp1, RDouble ypn, RDouble *y2);

public:
    void InitSpectrumData();

    //void GetAssumedSpectrum(RDouble ek, RDouble modk);

    //LIB_EXPORT void InitControlParameters();

    //LIB_EXPORT Param_Spectrum *GetControlParameters() const;
};

}

