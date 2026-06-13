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
//! @file      Geometry.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "LIB_Macro.h"
#include "TypeDefine.h"
#include "Geo_Grid.h"

using namespace std;

namespace PHSPACE
{

class PHGeometry
{
private:
    vector < Grid * > *grids;
public:
    PHGeometry(vector < Grid * > *grids);
    ~PHGeometry();
public:
    Grid * GetGrid(int level = 0);
    int GetNumberOfMultiGrid() const { return static_cast<int>(grids->size()); }
    void AddGrid(Grid *grid);
    vector < Grid * > * GetVectorGrid();
};

}

