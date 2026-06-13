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
//! @file      HOBoundary.h
//! @brief     basis function, such as monomial, orthogonalized etc.
//! @author    gongxq, liming, Wan yunbo, Xu gang.
//! @date      2019-11-09
//! @version   A001

#pragma once

#include "Precision.h"

namespace HOUnstruct
{
    void SymmetryBoundary(int iGausspoint, int numberGaussPoint, RDouble xNormal, RDouble yNormal, RDouble zNormal, vector < vector < RDouble > > & conservationQ);
    void SymmetryBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem);
    void ViscousAdiabaticWallBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xfv, RDouble yfv, RDouble zfv, RDouble refMachNumber, int nl, int nchem);
    void ViscousIsotropicWallBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xfv, RDouble yfv, RDouble zfv, RDouble refMachNumber, RDouble tw, int nl, int nchem);
    void FarfieldBC(RDouble *prims, RDouble *primInflow, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem);
    void FarFieldRiemannInvariants(RDouble *prims, RDouble *primInflow, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem);
    void InflowBC(RDouble *prims, RDouble *primt, int nl, int nchem);
    void OutflowBC(RDouble *prims, RDouble *primt, int nl, int nchem);
}