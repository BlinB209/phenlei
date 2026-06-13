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
//! @file      PoissonSolver.h
//! @brief     Poisson solver for SpecDiffHybSolver.
//! @author    ZhangZipei.

#pragma once
#include "SpecDiffHybGrid.h"
#include "ExplicitDifferenceBoundary.h"
#include "CompactDifferenceSecondDerivative.h"
#include "Complex.h"
//#include "SpecDiffHybSolver.h"
#include "Param_SpecDiffHybSolver.h"

using namespace std;

namespace PHSPACE
{

class PoissonSolver
{
public:
    PoissonSolver();
    ~PoissonSolver();

public:
    RDouble4D *coefMatrixVelocityDirichlet, *coefMatrixPressureNeumann, *boundaryMatrixPressureNeumann;

public:
    void InitPoissonSolverData(SpecDiffHybGrid *, ExplicitDifferenceBoundary *, CompactDifferenceSecondDerivative*, RDouble1D *, RDouble1D *);
    
    void PoissonSolverNeumannBC( const int, RDouble *, PHComplex *, int *, int *, const int, const PHComplex, const PHComplex, PHComplex *, RDouble *, RDouble * );
    
    static void PoissonSolverDirichletBC( const int , RDouble *, RDouble *, PHComplex *, const RDouble, int *, int *, const PHComplex, const PHComplex, RDouble *, PHComplex * );
    
    static void PoissonSolverVelocity( const int, RDouble *, PHComplex *, int *, int *, const PHComplex, const PHComplex, RDouble *, RDouble *, PHComplex * );

    static void PoissonSolverVelocity(PHComplex &uLocal, const RDouble timeStep, const RDouble k2, const RDouble phi);
    
    static void LUDecomposeNeumannBC( int, RDouble *, int, RDouble * );
    
    static void LUDecompose( const int, RDouble * );
    
    static void LUSolve( const int, RDouble *, RDouble * );
    
    static void LUSolve( const int, RDouble *, PHComplex * );
    
    static void LUSolveNeumannBC( const int, RDouble *, const int, RDouble *, PHComplex * );
    
    void GetCoefMatrixVelocityDirichlet( SpecDiffHybGrid *, ExplicitDifferenceBoundary *, CompactDifferenceSecondDerivative *, RDouble1D *, RDouble1D *);

private:
    void AllocPoissonSolverData(SpecDiffHybGrid *, ExplicitDifferenceBoundary *);
    
    void GetCoefMatrixPressureNeumann( SpecDiffHybGrid *, ExplicitDifferenceBoundary *, CompactDifferenceSecondDerivative *);
    
    void DecomCoefMatrixPressureNeumannBC( const int, RDouble *, RDouble *, const double, int *, int *, const int, RDouble *, RDouble *, RDouble *, RDouble * );
    
    void DecomCoefMatrixVelocityDirichletBC( const int, RDouble *, RDouble *, const double, int *, int *, const double, RDouble *, RDouble * );
};

}



