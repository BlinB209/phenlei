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
//! @file      Param_INCompSolverUnstruct.h
//! @brief     pressure poission solver.
//! @author    WanYunbo, Bell.

#pragma once

#include <map>
#include <iostream>
#include "PPEqCalculator.h"
#include "IncomCalculator.h"


namespace PHSPACE
{
    class UnstructGrid;

    class PPEqCalculator : public IncomCalculator
    {

    public:
        PPEqCalculator();
        ~PPEqCalculator();

    public:
        void correctPressureAndVelocity(Grid *grid);
        void correctBoundaryVelocity(Grid *grid);
        void UpdateBCValue(Grid *grid);
        
        void solvePPEquation(Grid *grid);
        void Post() {};
        
        void treatBC(Grid *grid);
        void constructMatrixACoeff(Grid *grid);
        void constructBCoeff(Grid *grid);
        void calcTransMatrixTerm(Grid *grid);
        void TransMatrixTerm_1st(Grid *grid);
        void TransMatrixTerm_2nd(Grid *grid);

    };
}
