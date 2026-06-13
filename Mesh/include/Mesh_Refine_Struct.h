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
//! @file      Mesh_Refine_Struct.h
//! @brief     Structured grid refine.
//! @author    Baka, He Xianyao.

#pragma once
#include "Region.h"
#include "Geometry.h"
#include "Mesh_Refine.h"

using namespace std;

namespace PHSPACE
{

class Mesh_Refine_Struct : public Mesh_Refine
{
public:
    LIB_EXPORT Mesh_Refine_Struct(const string &from_gfile);
    LIB_EXPORT ~Mesh_Refine_Struct();

private:
    //! Initialize: allocate memory.
    void AllocateMemory();

    //! Grid refine.
    void RefineGrid();
    void RefineGrid(int iZone, StructGrid *grid_in);

    //! Generate refined grid, and then dump grid.
    void GenerateAndDumpComputationalGrid();
    void GenerateBCRegion();

private:
    vector<Grid *> OrdinaryGrid;
};

}