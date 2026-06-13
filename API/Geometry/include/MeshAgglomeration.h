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
//! @file      Mesh_Agglomeration.h
//! @brief     Mesh coarsen using agglomeration method for unstructured grid.
//! @author    Bell.

#pragma once
#include "Geo_UnstructGrid.h"
using namespace std;

namespace PHSPACE
{

//! @brief Mesh coarsen using agglomeration method for unstructured grid.
class Mesh_Agglomeration
{
public:
    //! @param[in] fineGridIn     The pre-fine grid, which is need to be coarsen.\n
    //! @param[in] fineGridLevel  The level of the fineGrid.   
    LIB_EXPORT Mesh_Agglomeration(UnstructGrid *fineGridIn, int fineGridLevelIn);

    LIB_EXPORT ~Mesh_Agglomeration();

    //! Coarse the fine grid, this function is recommended,
    //! since the coarse grid quality is guaranteed.
    LIB_EXPORT void CoarsenGridOnce();

    LIB_EXPORT void CoarsenGridOncebyMGridgen();

    //! Coarse the fine grid, fast but not robust, this function is NOT recommended.
    LIB_EXPORT void CoarsenGridOnceWithoutSkewnessLimit();

    //! Return the agglomerated coarse grid.
    LIB_EXPORT UnstructGrid * GetCoarseGrid() { return this->coarseGrid; }

private:
    //! The pre-fine grid, which is need to be coarsen.
    UnstructGrid *fineGrid;

    //! The coarse grid, which is agglomerated from fineGrid.
    UnstructGrid *coarseGrid;

    //! The level of the fineGrid.
    int fineGridLevel;

private:
    void AspectRatioModifyBySkewness(RDouble *ratio, RDouble *xcc, RDouble *ycc, RDouble *zcc, const RDouble &limitAngle);
};

}