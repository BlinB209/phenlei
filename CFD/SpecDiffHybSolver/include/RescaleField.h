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
//! @file      RescaleField.h
//! @brief     Rescale velocity every time step if necessary.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "SpecGrid.h"
#include "SpecDiffHybGrid.h"
//#include "Param_SpecSolver.h"

using namespace std;

namespace PHSPACE
{
class RescaleFieldHIT
{
public:
    RescaleFieldHIT();
    ~RescaleFieldHIT();

public:
    void InitRescaleData();

private:
    void AllocRescaleData();

private:
    RDouble uFluctuation, vFluctuation, wFluctuation;
    RDouble1D *rescaleStatistic, *rescaleStatisticLocal;
};

class RescaleField
{
public:
    RescaleField();
    ~RescaleField();

private: 
    RDouble uavgDivideByUtau;
    RDouble1D *uAverage;
    RDouble1D *uAverageOfReference, *uFluctuationOfReference, *vFluctuationOfReference, *wFluctuationOfReference;
    RDouble1D *uFluctuation, *vFluctuation, *wFluctuation;
    RDouble2D *rescaleStatistic, *rescaleStatisticLocal;

public:
    void InitRescaleData(SpecDiffHybGrid *GridData);

    void GetUVWMagnitude(SpecDiffHybGrid *, Complex3D *, Complex3D *, Complex3D *);
    
    void RescaleUAverage(SpecDiffHybGrid *, Complex3D *);
    
    void RescaleFluctuation(SpecDiffHybGrid *, Complex3D *, Complex3D *, Complex3D *);
    
    void FixUAverage(SpecDiffHybGrid *, Complex3D *);

    //LIB_EXPORT Param_SpecSolver * GetControlParameters() const;

    //LIB_EXPORT void InitControlParameters();

private:
    void AllocRescaleData(SpecDiffHybGrid *GridData);

    void SetUavg();

    void GetUVWMagnitudeOfReference(SpecDiffHybGrid *);

    void InputUVWMagnitudeOfReference( SpecDiffHybGrid * );

    void GetUVWLocal(SpecDiffHybGrid *, Complex3D *, Complex3D *, Complex3D *);

    void Interpolation( RDouble *, RDouble *, int, RDouble *, RDouble *, int );
};

}



