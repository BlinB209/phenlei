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
//! @file      GridFactory.h
//! @brief     Some functions about grid.
//! @author    He Xin, Bell, Baka.

#pragma once
#include "MixGrid.h"

namespace PHSPACE
{
//! 
void CreateGrid();

void ConvertGrid();

void RefineGrid();

void CombineGrid();

void RepairGrid();

void MirrorGrid();

void StructToUnstructGrid();

//!
void OversetGridView();

void Fantasy2Fluent();
}
