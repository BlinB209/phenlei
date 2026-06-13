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
//! @file      hypreCalculator.h
//! @brief     call hypre lib.
//! @author    WanYunbo, Bell, XuGang.

#pragma once

#include "Geo_UnstructGrid.h"
#include "LinearEquationSystemCalculator.h"

#ifdef USE_HypreLib
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "PHMpi.h"
#include "_hypre_parcsr_ls.h"
#endif

namespace PHSPACE
{
class hypreCalculator : public LinearEquationSystemCalculator
{

#ifdef USE_HypreLib
private:
    Grid *grid;
    HYPRE_IJMatrix A_hypre;
    HYPRE_ParCSRMatrix parcsr_A_hypre;
    HYPRE_IJVector b_hypre;
    HYPRE_ParVector par_b_hypre;
    HYPRE_IJVector x_hypre;
    HYPRE_ParVector par_x_hypre;
    HYPRE_Solver **solverList;
    HYPRE_Solver **precondList;
    vector<int> *ColVectors;
    vector<double> *ValuesVector;

    int ilower;
    int iupper;
    int localSize;
    int *rows;
    RDouble *colValueOfOneRow;
    int *colNumberOfOneRow;
#endif
public:
    // constructor & destructor
    hypreCalculator();

    virtual ~hypreCalculator();

public:
    void mathLibSolveInit(Grid *gridIn);
    void mathLibSolve(Grid *gridIn, int iEquation, RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace);
    void mathLibSolveFinalize();

#ifdef USE_HypreLib
    void SetAMGDDPrecondParameters(HYPRE_Solver &amg_precond);
    void destoryTempDataInAMGDD(HYPRE_Solver &amgddsolver);
#endif
};

}

