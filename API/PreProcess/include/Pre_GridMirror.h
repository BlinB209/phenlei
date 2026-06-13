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
//! @file      Pre_GridMirror.h
//! @brief     It defines grid mirror function for pre-processing.
//! @author    Baka.

#pragma once
#include "Region.h"
#include "Geometry.h"

using namespace std;

namespace PHSPACE
{
class Pre_GridMirror
{
public:
    string gridFileName;

public:
    LIB_EXPORT Pre_GridMirror(const string &from_gfile);
    LIB_EXPORT ~Pre_GridMirror();

public:
    LIB_EXPORT void Run();

private:
    //! Revert boundary condition information.
    void RevertBCInfo();

    void ReadLnkFile(const string &fileName_in);

    //! Calculate the symmetry face.
    void CalSymmetryFace();

    //! Calculate the blocks number of the whole mirror grid.
    void CalBlocksNumberOfMirrorGrid();

    //! Mirror the grid.
    void GridMirror();

    //! Mirror the structured grid.
    void MirrorGridStr(StructGrid *grid_in);

    //! Mirror the boundary condition region.
    void MirrorBCRegion();

    //! Mirror the unstructured grid.
    void MirrorGridUnstr(UnstructGrid *grid_in);

    //! Mirror the coordinates of unstructured grid.
    void MirrorCoordinate(UnstructGrid *grid_in);

    //! Mirror the node number of each face of unstructured grid.
    void MirrorNodeNumberOfEachFace(UnstructGrid *grid_in);

    void MirrorFace2Node(UnstructGrid *grid_in);

    //! Mirror the cell of each face of unstructured grid.
    void MirrorCellOfFace(UnstructGrid *grid_in);

    void MirrorCell2Node(UnstructGrid *grid_in);

    void TradeGridInfo();
    void DumpStrGrid(string out_gfile);
    void DumpUnstrGrid(string out_gfile);

private:
    int numberOfZones;
    int mirrorindex;
    int Dimension;
    int numberOfOriginalBlocks;
    int SymmetryFaceVector;
    int gridtype;
    int BlockofEachProc;
    bool LnkFileExist;

    vector<Grid *> OrdinaryGrid;
    //Grid ** OrdinaryGrid;
    vector<DataContainer *> Gridlist;
    DataContainer *gridinfo;
    /*vector<Grid *> MirrorGrids;
    Grid ** OutputGrid;*/
};

const int X_axis = 0;
const int Y_axis = 1;
const int Z_axis = 2;

}