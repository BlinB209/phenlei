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
//! @file      unapCalculator.h
//! @brief     call unap lib.
//! @author    WanYunbo, Bell.

#pragma once

#include "Geo_UnstructGrid.h"
#include <sstream>
#include "LinearEquationSystemCalculator.h"

#ifdef USE_UnapLib
#include "PBiCGStab.hpp"
#include "PCG.hpp"
#include "GMRES.hpp"
#include "chebySmoother.hpp"
#include "lduAgglomeration.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"
#include "lduDiagPrecond.hpp"
#include "lduGaussSeidelSmoother.hpp"
#include "matrixConversion.hpp"
#include "printUNAP.hpp"
#include "rcmf.h"
#include "readFromOpenFOAM.hpp"
#include "unapMultigrid.hpp"
#include "directionalAgglomeration.hpp"
#include "lduMGPrecond.hpp"
#include "unapLduMatrix.hpp"

using namespace UNAP;
#endif

namespace PHSPACE
{
class unapCalculator : public LinearEquationSystemCalculator
{
#ifdef USE_UnapLib
private:
    LduMatrix* lduA;

    labelVector* upperAddr;
    labelVector* lowerAddr;
    scalarVector* upper;
    scalarVector* lower;
    scalarVector* diag;
    scalarVector* b;
    scalarVector* x;
    SolverPerformance solverPerf;

    UNAP::Solver* unapSolver;
    UNAP::Preconditioner* unapPrecond;

    label nNeiProcs;
    label* destRank;
    label* offDiagRows;
    label* offDiagCols;
    label* offDiagStarts;
    scalar* offDiagCoeffs;
#endif
public:

    unapCalculator();

    virtual ~unapCalculator();

public:
    void mathLibSolveInit(Grid *gridIn);
    void mathLibSolve(Grid *gridIn, int iEquation, RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace);
    void mathLibSolveFinalize();
};


}

