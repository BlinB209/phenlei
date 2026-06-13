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
//! @file      CorrectWallDivergence.h
//! @brief     Ensure divergence on wall equals zero for SpecDiffHybSolver.
//! @author    ZhangZipei.

#pragma once

#include "Complex.h"
#include "SpecDiffHybGrid.h"
#include "CompactDifferenceFirstDerivative.h"
#include "CompactDifferenceSecondDerivative.h"
#include "PoissonSolver.h"
#include "SpecDiffHybSolver.h"
#include "Param_SpecDiffHybSolver.h"

using namespace std;

namespace PHSPACE
{

class CorrectWallDivergence
{
public:
    CorrectWallDivergence();
    ~CorrectWallDivergence();

public:
    void InitCorrectWallDivergenceData( SpecDiffHybGrid *, CompactDifferenceFirstDerivative *, CompactDifferenceSecondDerivative *, PoissonSolver *, RDouble1D * );

    void rotProjection();

private:   

public:
    Complex3D *complexDeltaU1, *complexDeltaV1, *complexDeltaW1, *complexDeltaP1;
    Complex3D *complexDeltaU2, *complexDeltaV2, *complexDeltaW2, *complexDeltaP2;
    Complex2D *complexDivDeltaU1LowerBound, *complexDivDeltaU1UpperBound, *complexDivDeltaU2LowerBound, *complexDivDeltaU2UpperBound;

private:
    void SetCorrectWallDivergenceData( SpecDiffHybGrid *, CompactDifferenceFirstDerivative *, CompactDifferenceSecondDerivative *, PoissonSolver *, RDouble1D * );
    void AllocCorrectWallDivergenceData( SpecDiffHybGrid * );
};

}



