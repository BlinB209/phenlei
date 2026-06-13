#include "Mesh_CellTopo.h"
#include "PHHeader.h"
#include "Geo_Element.h"

namespace PHSPACE
{

Mesh_CellTopo::Mesh_CellTopo()
{
    cellType.resize(0);
    cellLevel.resize(0);
    parentCellIndex.resize(0);

    computationalToGlobalCellIndex.resize(0);
    globalToComputationalCellIndex.resize(0);
}

Mesh_CellTopo::~Mesh_CellTopo()
{

}

int Mesh_CellTopo::GetNTotalCell()
{
    return static_cast<int>(cellType.size());
}

void Mesh_CellTopo::SetAllCellLevelToZero()
{
    int nTotalCell = GetNTotalCell();

    PHVectorInt1D &cellLevel = GetCellLevel();
    cellLevel.resize(nTotalCell);

    SetField(cellLevel, 0, nTotalCell);
}

void Mesh_CellTopo::SetAllParentIndexToRootState()
{
    int nTotalCell = GetNTotalCell();

    PHVectorInt1D &parentCellIndex = GetParentCellIndex();
    parentCellIndex.resize(nTotalCell);

    SetField(parentCellIndex, PHSPACE::INVALID_INDEX, nTotalCell);
}

void Mesh_CellTopo::SetCellType(int cellIndex, int cellTypeIn)
{
    this->cellType[cellIndex] = cellTypeIn;
}

void Mesh_CellTopo::SetCellModifiedStatus(int cellIndex, int cellModifiedStatusIn)
{
    this->cellModifiedStatus[cellIndex] = cellModifiedStatusIn;
}

void Mesh_CellTopo::SetComputationalStatusOfCell(int cellIndex, int computationalStatusOfCell)
{
    this->cellComputationalStatus[cellIndex] = computationalStatusOfCell;
}

void Mesh_CellTopo::SetCell2Node(int cellIndex, PHVectorInt1D &cell2NodeIn)
{
    this->cell2Node[cellIndex] = cell2NodeIn;
}

void Mesh_CellTopo::SetParentCellIndex(int cellIndex, int parentCellIndexIn)
{
    this->parentCellIndex[cellIndex] = parentCellIndexIn;
}

void Mesh_CellTopo::SetCellLevel(int cellIndex, int cellLevelIn)
{
    this->cellLevel[cellIndex] = cellLevelIn;
}

bool Mesh_CellTopo::ElementNeedSplit(int cellIndex)
{
    return cellModifiedStatus[cellIndex] == SPLIT;
}

void Mesh_CellTopo::AddCell(PHVectorInt1D &cell2Node, int cellType, int cellLevel, int parentCellIndex)
{
    int cellIndex = GetNTotalCell();
    ResizeCellNumber(cellIndex + 1);

    SetCell2Node(cellIndex, cell2Node);
    SetCellType(cellIndex, cellType);
    SetParentCellIndex(cellIndex, parentCellIndex);
    SetCellLevel(cellIndex, cellLevel);

    SetCellModifiedStatus(cellIndex, NOCHANGE);
    SetComputationalStatusOfCell(cellIndex, ON);

    PHVectorInt1D &childrenCellIndexe = GetChildrenCellIndexe(parentCellIndex);
    childrenCellIndexe.push_back(cellIndex);

    globalToComputationalCellIndex[cellIndex] = INVALID_INDEX;
}

void Mesh_CellTopo::ResizeCellNumber(int nTotalCellIn)
{
    this->cell2Node.resize(nTotalCellIn);
    this->cellType.resize(nTotalCellIn);
    this->cellLevel.resize(nTotalCellIn);
    this->cellModifiedStatus.resize(nTotalCellIn);
    this->cellComputationalStatus.resize(nTotalCellIn);
    this->parentCellIndex.resize(nTotalCellIn);
    this->childrenCellIndexes.resize(nTotalCellIn);
    this->globalToComputationalCellIndex.resize(nTotalCellIn);
}

}