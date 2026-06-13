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
//! @file      Pre_StructToUnstruct.h
//! @brief     Grid conversion from struct grid to unstruct grid.
//! @author    Baka.

#pragma once

using namespace std;

namespace PHSPACE
{
class Grid;
class StructGrid;
class Pre_StructToUnstruct
{
private:
    int numberOfZones;
    Grid **structGrid;
    Grid **structToUnstruct;

public:
    Pre_StructToUnstruct(Grid **structGridIn, int nZonesIn);
    ~Pre_StructToUnstruct();

public:
    void Run();
    Grid ** GetUnstructGrids();

private:
    void InitMemory();
    void CopyCoordinate();
    void CopyIBlock();
    void ConstructFaceTopo();
    void ConstructBCType();
    void CopyInterface();
    int GetCellIndex(StructGrid *strGrid, int iCell, int jCell, int kCell);
    int GetNodeIndex(StructGrid *strGrid, int iNode, int jNode, int kNode);
    void ConstructCellMapStrToUns();

};

}