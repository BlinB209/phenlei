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
//! @file      GridType.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "LIB_Macro.h"
#include "TypeDefine.h"
#include "Geo_GridIndex.h"

using namespace std;

namespace PHSPACE
{

const int UNSTRUCTGRID = 0;
const int STRUCTGRID   = 1;
const int MIXGRID      = 2;

Grid * CreateGridGeneral(int gridtype, GridID *index, int level, int nDim);
int GetSystemGridType(set<int> &grid_type);
}
