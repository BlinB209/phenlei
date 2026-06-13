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
//! @brief     call yhamg lib.
//! @author    WanYunbo, Bell.

#pragma once

#include "Geo_UnstructGrid.h"
#include "LinearEquationSystemCalculator.h"
#ifdef USE_YhamgLib
#include "yhamg.h"
using namespace YHAMG;
#endif
namespace PHSPACE
{
    
    class yhamgCalculator : public LinearEquationSystemCalculator
    {
#ifdef USE_YhamgLib
    private:
        Grid *grid;
        ParOperator *A;
        ParOperator **PrecondList;
        vector<int> *ColVectors;
        int *exter_colind2Face;

        ParVector *b;
        ParVector *x;
#endif

    public:
        // constructor & destructor
        yhamgCalculator();

        virtual ~yhamgCalculator();

    public:
        void mathLibSolveInit(Grid *gridIn);
        void mathLibSolve(Grid *gridIn, int iEquation, RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace);
        void mathLibSolveFinalize();

    };

}

