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
//! @file      Mesh_CellTopo.h
//! @brief     Store cell information of grid, used in unstructured grid refining.
//! @author    Baka, Bell.

#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{
//! @brief Mesh_CellTopo store cell information of grid, used in unstructured grid refining.
class Mesh_CellTopo
{
private:
    //! Cell type, such as TETRA_4, PENTA_6.
    PHVectorInt1D cellType;

    //! Cell level, parent cell or child cell.
    PHVectorInt1D cellLevel;

    //! Cell refine type, ISOTROPICREFINE or ANISOTROPICREFINE.
    PHVectorInt1D refineType;

    //! Parent cell index of each cell.
    PHVectorInt1D parentCellIndex;

    //! Children cell index of each cell.
    PHVectorInt2D childrenCellIndexes;

    //! Cell index map between refined cell and global cell.
    PHVectorInt1D computationalToGlobalCellIndex;
    PHVectorInt1D globalToComputationalCellIndex;

    //! Cell modified status.
    PHVectorInt1D computationalCellModifiedStatus;
    PHVectorInt1D cellModifiedStatus;
    PHVectorInt1D cellComputationalStatus;

    //! Cell topo information.
    PHVectorInt2D cell2Node;
    PHVectorInt2D cell2Face;
    PHVector1D< set< int > > node2Cell;
    PHVector1D< set< int > > cell2Cell;

public:
    Mesh_CellTopo();
    ~Mesh_CellTopo();

public:
    PHVectorInt1D & GetCellType() { return cellType; }
    PHVectorInt1D & GetCellLevel(){ return cellLevel; }
    PHVectorInt1D & GetRefineType() { return this->refineType; }
    PHVectorInt1D & GetComputationalToGlobalCellIndex() { return computationalToGlobalCellIndex; }
    PHVectorInt1D & GetGlobalToComputationalCellIndex() { return globalToComputationalCellIndex; }
    PHVectorInt1D & GetParentCellIndex() { return parentCellIndex; }
    PHVectorInt1D & GetComputationalCellModifiedStatus(){ return computationalCellModifiedStatus; }
    PHVectorInt1D & GetCellModifiedStatus() { return cellModifiedStatus; }
    PHVectorInt1D & GetCellComputationalStatus() { return cellComputationalStatus; }

    PHVectorInt2D & GetCell2Node() { return cell2Node; }
    PHVectorInt2D & GetCell2Face() { return cell2Face; }
    PHVectorInt1D & GetCell2Node(int cellIndex) { return cell2Node[cellIndex]; }
    PHVectorInt1D & GetChildrenCellIndexe(int cellIndex) { return childrenCellIndexes[cellIndex]; }

    PHVector1D< set< int > > & GetNode2Cell() { return node2Cell; }
    PHVector1D< set< int > > & GetCell2Cell() { return cell2Cell; }

    int GetNTotalCell();
    int GetCellType(int cellIndex) { return cellType[cellIndex]; }
    int GetCellLevel(int cellIndex) { return cellLevel[cellIndex]; }
    int GetRefineType(int cellIndex) { return refineType[cellIndex]; }
    int GetCellComputationalStatus(int cellIndex) { return cellComputationalStatus[cellIndex]; }
    int GetCellModifiedStatus(int cellIndex) { return cellModifiedStatus[cellIndex]; }

    void SetAllCellLevelToZero();
    void SetAllParentIndexToRootState();
    void SetCellType(int cellIndex, int cellTypeIn);
    void SetCellModifiedStatus(int cellIndex, int cellModifiedStatusIn);
    void SetComputationalStatusOfCell(int cellIndex, int computationalStatusOfCell);
    void SetCell2Node(int cellIndex, PHVectorInt1D &cell2NodeIn);
    void SetParentCellIndex(int cellIndex, int parentCellIndexIn);
    void SetCellLevel(int cellIndex, int cellLevelIn);

    //! Add new cell.
    void AddCell(PHVectorInt1D &cell2Node, int cellType, int cellLevel, int parentCellIndex);
    void ResizeCellNumber(int nTotalCellIn);

    //! Return if need split.
    bool ElementNeedSplit(int cellIndex);
};

const int SPLIT = 0;
const int NOCHANGE = 1;
const int MERGE_ELEMENT = 2;
const int DELETE = 3;
const int HIDDEN = 4;
const int ON = 5;
const int BOUNDARYFACE = 6;
const int HANGINGFACE = 7;
const int GENERALFACE = 8;

const int ISOTROPICREFINE = 0;
const int ANISOTROPICREFINE = 1;

}