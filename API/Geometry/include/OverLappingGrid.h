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
//! @file      OverLappingGrid.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "Geo_Grid.h"
#include "GridType.h"

namespace PHSPACE
{
class Grid;
class StructGrid;
class UnstructGrid;

class Geo_OversetGrid;
void DumpFlags(Geo_OversetGrid **oversetGrids, int nBlocks, const string &filename);
void DumpOversetGrid(Geo_OversetGrid **oversetGrids, int nBlocks, const string &filename);

}
