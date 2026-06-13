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
//! @file      Geo_CellTopo_Unstruct.h
//! @brief     It defines the cell topology of the unstructured grid,
//!            such as the cell-face and cell-node connected relationship.
//! @author    Bell, He Xin, Tang Jin.
#pragma once
#include <vector>
#include "LIB_Macro.h"

namespace PHSPACE
{
//! @brief Geo_CellTopo_Unstruct class defines the cell-face, cell-node connected relationship of unstructured grid.\n
//! 1: faces of each cell.\n
//! 2: nodes of each cell.\n
//! 3: neighbor cells of each cell.
class Geo_CellTopo_Unstruct
{
private:
    //! Number of faces per cell.
    int *face_number_of_each_cell;

    //! Cell to face: face index of each cell.
    int **cell2face;

    //! Number of cells per cell.
    int *node_number_of_each_cell;

    //! Cell to node: node index of each cell.
    int *cell2node;

    //! Cell to node: node index of each cell.
    int **cell2nodeArray;

    //! Number of neighbor cells per cell;
    int *number_of_neighbor_cell;

    //! Cell to cell: neighbor cell index of each cell.
    vector<int> *cell2cell;
    // GMRESVis
    vector<int> * neighborCells;
    vector<int> * neighborFaces;
    vector<int> * neighborLR;

public:
    LIB_EXPORT Geo_CellTopo_Unstruct();
    LIB_EXPORT ~Geo_CellTopo_Unstruct();

public:
    //! Set the number of faces per cell.
    void SetFaceNumberOfEachCell(int *face_number_of_each_cell);

    //! Set the cell to face.
    void SetCell2Face(int **cell2faceIn);

    //! Set the number of cells per cell.
    void SetNodeNumberOfEachCell(int *node_number_of_each_cell);

    //! Set the cell to node.
    void SetCell2Node(int *cell2node);

    //! Set the cell to node.
    void SetCell2NodeArray(int **cell2nodeArrayIn);

    //! Set number of neighbor cells per cell.
    void SetNumberOfNeighborCell(int *number_of_neighbor_cell);

    //! Set the cell to cell: neighbor cell index of each cell.
    void SetCell2Cell(vector<int> *cell2cell);
    
    // GMRESVis
    void GMRESSetNeighborCells(vector<int> * neighborCells);
    void GMRESSetNeighborFaces(vector<int> * neighborFaces);
    void GMRESSetNeighborLR(vector<int> * neighborLR);

    //! Get the number of faces per cell.
    int * GetFaceNumberOfEachCell();

    //! Get the cell to face: face index of each cell.
    int ** GetCell2Face();

    //! Get number of nodes per cell.
    int * GetNodeNumberOfEachCell();

    //! Get the cell to node: node index of each cell.
    int * GetCell2Node();

    //! Get the cell to node: node index of each cell.
    int ** GetCell2NodeArray();

    //! Get number of neighbor cells per cell.
    int * GetNumberOfNeighborCell();

    //! Get the cell to cell: neighbor cell index of each cell.
    vector<int> * GetCell2Cell();
    
    // GMRESVis
    vector<int> * GMRESGetNeighborCells();
    vector<int> * GMRESGetNeighborFaces();
    vector<int> * GMRESGetNeighborLR();
};

#include "Geo_CellTopo_Unstruct.hxx"
}
