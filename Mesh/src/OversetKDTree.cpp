#include "OversetKDTree.h"
#include "Geo_Grid.h"
#include "Pointer.h"

namespace PHSPACE
{

BasicKDTree::BasicKDTree(RDouble **dataIn,
    int         numberOfSearchElementsIn,
    int         numberOfSearchDirectionIn) :
    data(dataIn),
    numberOfSearchElements(numberOfSearchElementsIn),
    numberOfSearchDirection(numberOfSearchDirectionIn)
{
    leftSonIndex = new int[numberOfSearchElements];
    rightSonIndex = new int[numberOfSearchElements];
    SetField(leftSonIndex, - 1, numberOfSearchElements);
    SetField(rightSonIndex, - 1, numberOfSearchElements);

    if (numberOfSearchElements == 1)
    {
        rootIndex = 0;
        return;
    }

    BuildBasicKDTree();
}

BasicKDTree::~BasicKDTree()
{
    data = NULL;

    delete[] leftSonIndex;
    delete[] rightSonIndex;

    leftSonIndex = NULL;
    rightSonIndex = NULL;
}

void BasicKDTree::BuildBasicKDTree()
{
    partitionMethod = & BasicKDTree::PartitionSinglePivotPreConditionDoubleEndScanQuickSort;

    elementOrder = new int[numberOfSearchElements];
    for (int iElement = 0; iElement < numberOfSearchElements; ++ iElement)
    {
        elementOrder[iElement] = iElement;
    }

    int iLayer = 0;
    int iStart = 0;
    int iFinal = numberOfSearchElements - 1;

    TreeSortSinglePivotQuickSort(iStart, iFinal, iLayer);

    rootIndex = elementOrder[iFinal / 2];

    delete[] elementOrder;
    elementOrder = NULL;
}

void BasicKDTree::SwapElementIndex(int indexA, int indexB)
{
    int      tempOrder = elementOrder[indexA];
    elementOrder[indexA] = elementOrder[indexB];
    elementOrder[indexB] = tempOrder;
}

int  BasicKDTree::PartitionSinglePivotDoubleEndScanQuickSort(int iStart, int iFinal, int iLayer)
{
    int     iPivot = iFinal --;

    while (iStart <= iFinal)
    {
        while (iStart <= iFinal && data[iLayer][elementOrder[iStart]] < data[iLayer][elementOrder[iPivot]])
        {
            iStart ++;
        }

        while (iStart <= iFinal && data[iLayer][elementOrder[iFinal]] >= data[iLayer][elementOrder[iPivot]])
        {
            iFinal --;
        }

        if (iStart < iFinal)
        {
            SwapElementIndex(iStart, iFinal);
        }
    }

    SwapElementIndex(iStart, iPivot);

    return iStart;
}

int  BasicKDTree::PartitionSinglePivotPreConditionDoubleEndScanQuickSort(int iStart, int iFinal, int iLayer)
{
    SwapElementIndex(iFinal, (iStart + iFinal) / 2);

    int     iPivot = iFinal --;

    while (iStart <= iFinal)
    {
        while (iStart <= iFinal && data[iLayer][elementOrder[iStart]] < data[iLayer][elementOrder[iPivot]])
        {
            iStart ++;
        }

        while (iStart <= iFinal && data[iLayer][elementOrder[iFinal]] >= data[iLayer][elementOrder[iPivot]])
        {
            iFinal --;
        }

        if (iStart < iFinal)
        {
            SwapElementIndex(iStart, iFinal);
        }
    }

    SwapElementIndex(iStart, iPivot);

    return iStart;
}

int  BasicKDTree::PartitionSinglePivotForwardScanQuickSort(int iPivot, int iFinal, int iLayer)
{
    int       indexA = iPivot;

    for (int indexB = iPivot + 1; indexB <= iFinal; ++ indexB)
    {
        if (data[iLayer][elementOrder[indexB]] < data[iLayer][elementOrder[iPivot]])
        {
            SwapElementIndex(++ indexA, indexB);
        }
    }

    SwapElementIndex(iPivot, indexA);

    return indexA;
}

int  BasicKDTree::PartitionSinglePivotPreConditionForwardScanQuickSort(int iPivot, int iFinal, int iLayer)
{
    SwapElementIndex(iPivot, (iPivot + iFinal) / 2);

    int       indexA = iPivot;

    for (int indexB = iPivot + 1; indexB <= iFinal; ++ indexB)
    {
        if (data[iLayer][elementOrder[indexB]] < data[iLayer][elementOrder[iPivot]])
        {
            SwapElementIndex(++ indexA, indexB);
        }
    }

    SwapElementIndex(iPivot, indexA);

    return indexA;
}

int  BasicKDTree::FindHalfSinglePivotQuickSort(int iStart, int iFinal, int iLayer)
{
    int middle = (iStart + iFinal) / 2;

    int iPoint = - 1;

    while (iPoint != middle)
    {
        iPoint = (this->*partitionMethod)(iStart, iFinal, iLayer);

        if (iPoint < middle)
        {
            iStart = iPoint + 1;
        }
        else if (iPoint > middle)
        {
            iFinal = iPoint - 1;
        }
        else
        {
            break;
        }
    }

    return iPoint;
}

void BasicKDTree::TreeSortSinglePivotQuickSort(int iStart, int iFinal, int iLayer)
{
    iLayer %= numberOfSearchDirection;

    int middle = FindHalfSinglePivotQuickSort(iStart, iFinal, iLayer);

    if (middle == iStart)
    {
        leftSonIndex[elementOrder[middle]] = - 1;
    }
    else if (middle == iStart + 1)
    {
        leftSonIndex[elementOrder[middle]] = elementOrder[iStart];
    }
    else
    {
        TreeSortSinglePivotQuickSort(iStart, middle - 1, iLayer + 1);
        leftSonIndex[elementOrder[middle]] = elementOrder[(iStart + middle - 1) / 2];
    }

    if (middle == iFinal - 1)
    {
        rightSonIndex[elementOrder[middle]] = elementOrder[iFinal];
    }
    else
    {
        TreeSortSinglePivotQuickSort(middle + 1, iFinal, iLayer + 1);
        rightSonIndex[elementOrder[middle]] = elementOrder[(middle + 1 + iFinal) / 2];
    }
}

GridKDTree::GridKDTree(UnstructGrid *gridIn, int *keySearchOfCellsIn, RDouble *minBoxOfZoneIn, RDouble *maxBoxOfZoneIn) : grid(gridIn)
{
    geometricDimension = GetDim();
    rootIndex = 0;
    leftSonIndex = NULL;
    rightSonIndex = NULL;

    ifAuxiliaryGrid = false;
    if (NULL == keySearchOfCellsIn)
    {
        ifAuxiliaryGrid = true;
    }

    int numberOfTotalCells = grid->GetNTotalCell();

    numberOfSearchCells = 0;
    if (ifAuxiliaryGrid)
    {
        for (int iCell = 0; iCell < numberOfTotalCells; ++ iCell)
        {
            searchCellIndex.push_back(iCell);
        }
        numberOfSearchCells = numberOfTotalCells;
    }
    else
    {
        for (int iCell = 0; iCell < numberOfTotalCells; ++ iCell)
        {
            if (keySearchOfCellsIn[iCell])
            {
                searchCellIndex.push_back(iCell);
                ++ numberOfSearchCells;
            }
        }
    }

    if (!numberOfSearchCells)
    {
        return;
    }

    switch (geometricDimension)
    {
    case TWO_D:
        inCellCheckMethod = & GridKDTree::InCellCheck2D;
        break;
    case THREE_D:
        inCellCheckMethod = & GridKDTree::InCellCheck3D;
        break;
    default:
        exit(- 1);
    }

    numberOfSearchDirection = 2 * geometricDimension;

    RDouble tolerance = GlobalDataBase::GetDoubleParaFromDB("toleranceForOversetBox");
    if (ifAuxiliaryGrid)
    {
        ComputeGlobalBox(0.0);
    }
    else
    {
        minBoxOfZone = minBoxOfZoneIn;
        maxBoxOfZone = maxBoxOfZoneIn;
    }

    ComputeCellBox(tolerance);

    basicKDTree = new BasicKDTree(cellBox, numberOfSearchCells,  2 * geometricDimension);

    rootIndex = basicKDTree->GetRootIndex();
    leftSonIndex = basicKDTree->GetLeftSonIndex();
    rightSonIndex = basicKDTree->GetRightSonIndex();
}

GridKDTree::~GridKDTree()
{
    if (!numberOfSearchCells)
    {
        return;
    }

    delete basicKDTree;

    leftSonIndex = NULL;
    rightSonIndex = NULL;

    if (ifAuxiliaryGrid)
    {
        delete[] minBoxOfZone;
        delete[] maxBoxOfZone;
        minBoxOfZone = NULL;
        maxBoxOfZone = NULL;
    }
    else
    {
        minBoxOfZone = NULL;
        maxBoxOfZone = NULL;
    }

    DelPointer2(cellBox);
    cellBox = NULL;
}

void GridKDTree::ComputeGlobalBox(RDouble tolerance)
{
    minBoxOfZone = new RDouble[geometricDimension];
    maxBoxOfZone = new RDouble[geometricDimension];
    SetField(minBoxOfZone, LARGE, geometricDimension);
    SetField(maxBoxOfZone, - LARGE, geometricDimension);

    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    int numberOfNodes = grid->GetNTotalNode();
    for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
    {
        RDouble &minInThisDimension = minBoxOfZone[iDimension];
        RDouble &maxInThisDimension = maxBoxOfZone[iDimension];
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            RDouble &nodeCoordinateInThisDimension = nodeCoordinate[iDimension][iNode];
            minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
            maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
        }
        minInThisDimension -= tolerance;
        maxInThisDimension += tolerance;
    }
}

void GridKDTree::ComputeCellBox(RDouble tolerance)
{
    cellBox = NewPointer2< RDouble >(numberOfSearchDirection, numberOfSearchCells);

    int **cellNodeIndex = grid->GetCell2NodeArray();
    int *cellNodeNumber = grid->GetNodeNumberOfEachCell();
    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
    {
        SetField(cellBox[iDimension], LARGE, numberOfSearchCells);
        SetField(cellBox[iDimension + geometricDimension], - LARGE, numberOfSearchCells);
        for (int iCell = 0; iCell < numberOfSearchCells; ++ iCell)
        {
            RDouble &minInThisDimension = cellBox[iDimension][iCell];
            RDouble &maxInThisDimension = cellBox[iDimension + geometricDimension][iCell];
            int &cellIndex = searchCellIndex[iCell];
            int &numberOfNodesInThisCell = cellNodeNumber[cellIndex];
            for (int iNode = 0; iNode < numberOfNodesInThisCell; ++ iNode)
            {
                RDouble &nodeCoordinateInThisDimension = nodeCoordinate[iDimension][cellNodeIndex[cellIndex][iNode]];
                minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
                maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
            }
            minInThisDimension -= tolerance;
            maxInThisDimension += tolerance;
        }
    }
}
void GridKDTree::EnlargeCellBox(const int &cellIndex, const RDouble &eps)
{
    for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
    {
        RDouble &minInThisDimension = cellBox[iDimension][cellIndex];
        RDouble &maxInThisDimension = cellBox[iDimension + geometricDimension][cellIndex];
        minInThisDimension -= eps;
        maxInThisDimension += eps;
    }
}

int  GridKDTree::SearchNodeInKDTree(RDouble *nodeCoordinateIn)
{
    int donorCellIndex = -1;

    if (numberOfSearchCells)
    {
        if (IfNodeIsInBox(nodeCoordinateIn, minBoxOfZone, maxBoxOfZone, geometricDimension))
        {
            SearchNodeInKDTree(nodeCoordinateIn, rootIndex, donorCellIndex, 0);
        }
    }

    return donorCellIndex;
}

void GridKDTree::SearchNodeInKDTree(RDouble *nodeCoordinateIn, int rootIndexIn, int &donorCellIndexIn, int iLayerIn)
{
    if (rootIndexIn == - 1)
    {
        return;
    }

    iLayerIn %= numberOfSearchDirection;

    if (IfNodeIsInCellBox(nodeCoordinateIn, rootIndexIn))
    {
        if ((this->*inCellCheckMethod)(nodeCoordinateIn, rootIndexIn, donorCellIndexIn)) return;

        SearchNodeInKDTree(nodeCoordinateIn, leftSonIndex[rootIndexIn], donorCellIndexIn, iLayerIn + 1);
        if (donorCellIndexIn != - 1) return;
        SearchNodeInKDTree(nodeCoordinateIn, rightSonIndex[rootIndexIn], donorCellIndexIn, iLayerIn + 1);
    }
    else
    {
        if (IfInLeftHalfTree(nodeCoordinateIn, rootIndexIn, iLayerIn))
        {
            SearchNodeInKDTree(nodeCoordinateIn, leftSonIndex[rootIndexIn], donorCellIndexIn, iLayerIn + 1);
        }
        else if (IfInRightHalfTree(nodeCoordinateIn, rootIndexIn, iLayerIn))
        {
            SearchNodeInKDTree(nodeCoordinateIn, rightSonIndex[rootIndexIn], donorCellIndexIn, iLayerIn + 1);
        }
        else
        {
            SearchNodeInKDTree(nodeCoordinateIn, leftSonIndex[rootIndexIn], donorCellIndexIn, iLayerIn + 1);
            if (donorCellIndexIn != - 1) return;
            SearchNodeInKDTree(nodeCoordinateIn, rightSonIndex[rootIndexIn], donorCellIndexIn, iLayerIn + 1);
        }
    }
}

bool GridKDTree::InCellCheck2D(RDouble *nodeCoordinateIn, int rootIndexIn, int &donorCellIndexIn)
{
    int *leftCellIndex = grid->GetLeftCellOfFace();

    int *cellFaceNumber = grid->GetFaceNumberOfEachCell();
    int **cellFaceIndex = grid->GetCell2Face();

    RDouble eps = GetToleranceForOversetSearch();

    RDouble *faceNormalX = grid->GetFaceNormalX();
    RDouble *faceNormalY = grid->GetFaceNormalY();

    RDouble *faceCenterX = grid->GetFaceCenterX();
    RDouble *faceCenterY = grid->GetFaceCenterY();

    int &cellIndex = searchCellIndex[rootIndexIn];
    int &numberOfFacesInCell = cellFaceNumber[cellIndex];

    for (int iFace = 0; iFace < numberOfFacesInCell; ++ iFace)
    {
        int &faceIndex = cellFaceIndex[cellIndex][iFace];

        RDouble &thisFaceNormalX = faceNormalX[faceIndex];
        RDouble &thisFaceNormalY = faceNormalY[faceIndex];

        RDouble &thisFaceCenterX = faceCenterX[faceIndex];
        RDouble &thisFaceCenterY = faceCenterY[faceIndex];

        RDouble nominalAngle = (nodeCoordinateIn[0] - thisFaceCenterX) * thisFaceNormalX
            + (nodeCoordinateIn[1] - thisFaceCenterY) * thisFaceNormalY;

        if (leftCellIndex[faceIndex] == cellIndex) nominalAngle = - nominalAngle;

        if (nominalAngle < - eps) return false;

    }

    donorCellIndexIn = cellIndex;

    return true;
}

bool GridKDTree::InCellCheck3D(RDouble *nodeCoordinateIn, int rootIndexIn, int &donorCellIndexIn)
{
    int *cellFaceNumber = grid->GetFaceNumberOfEachCell();
    int **cellFaceIndex = grid->GetCell2Face();
    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();

    RDouble eps = GetToleranceForOversetSearch();

    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    RDouble *faceCenterX = grid->GetFaceCenterX();
    RDouble *faceCenterY = grid->GetFaceCenterY();
    RDouble *faceCenterZ = grid->GetFaceCenterZ();

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    int &cellIndex = searchCellIndex[rootIndexIn];
    int &numberOfFacesInCell = cellFaceNumber[cellIndex];

    RDouble &thisCellCenterX = cellCenterX[cellIndex];
    RDouble &thisCellCenterY = cellCenterY[cellIndex];
    RDouble &thisCellCenterZ = cellCenterZ[cellIndex];

    for (int iFace = 0; iFace < numberOfFacesInCell; ++ iFace)
    {
        int &faceIndex = cellFaceIndex[cellIndex][iFace];
        int &faceNodesNumber = faceNodeNumber[faceIndex];

        RDouble &thisFaceCenterX = faceCenterX[faceIndex];
        RDouble &thisFaceCenterY = faceCenterY[faceIndex];
        RDouble &thisFaceCenterZ = faceCenterZ[faceIndex];

        for (int iNode = 0; iNode < faceNodesNumber; ++ iNode)
        {
            int nodeIndex1 = faceNodeIndex[faceIndex][iNode];
            int nodeIndex2 = faceNodeIndex[faceIndex][(iNode + 1) % faceNodesNumber];

            RDouble &nodeCoordinateX1 = nodeCoordinateX[nodeIndex1];
            RDouble &nodeCoordinateY1 = nodeCoordinateY[nodeIndex1];
            RDouble &nodeCoordinateZ1 = nodeCoordinateZ[nodeIndex1];
            RDouble &nodeCoordinateX2 = nodeCoordinateX[nodeIndex2];
            RDouble &nodeCoordinateY2 = nodeCoordinateY[nodeIndex2];
            RDouble &nodeCoordinateZ2 = nodeCoordinateZ[nodeIndex2];

            if (In4PointsCell(nodeCoordinateX1, nodeCoordinateY1, nodeCoordinateZ1,
                nodeCoordinateX2, nodeCoordinateY2, nodeCoordinateZ2,
                thisFaceCenterX, thisFaceCenterY, thisFaceCenterZ,
                thisCellCenterX, thisCellCenterY, thisCellCenterZ,
                nodeCoordinateIn, eps))
            {
                donorCellIndexIn = cellIndex;
                return true;
            }
        }
    }
    return false;
}

bool GridKDTree::In4PointsCell(const RDouble &x1, const RDouble &y1, const RDouble &z1,
    const RDouble &x2, const RDouble &y2, const RDouble &z2,
    const RDouble &x3, const RDouble &y3, const RDouble &z3,
    const RDouble &x4, const RDouble &y4, const RDouble &z4,
    RDouble *nodeCoordinateIn, const RDouble &eps)
{
    RDouble volume = ABS(ComputeVolumeBy4Ppoints(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4));

    RDouble volume1 = ABS(ComputeVolumeBy4Ppoints(x1, y1, z1, x2, y2, z2, x3, y3, z3,
        nodeCoordinateIn[0], nodeCoordinateIn[1], nodeCoordinateIn[2]));
    RDouble volume2 = ABS(ComputeVolumeBy4Ppoints(x2, y2, z2, x3, y3, z3, x4, y4, z4,
        nodeCoordinateIn[0], nodeCoordinateIn[1], nodeCoordinateIn[2]));
    RDouble volume3 = ABS(ComputeVolumeBy4Ppoints(x3, y3, z3, x4, y4, z4, x1, y1, z1,
        nodeCoordinateIn[0], nodeCoordinateIn[1], nodeCoordinateIn[2]));
    RDouble volume4 = ABS(ComputeVolumeBy4Ppoints(x4, y4, z4, x1, y1, z1, x2, y2, z2,
        nodeCoordinateIn[0], nodeCoordinateIn[1], nodeCoordinateIn[2]));

    if(ABS((volume1 + volume2 + volume3 + volume4) / volume - 1.0) > eps)
    {
        return false;
    }
    return true;
}

RDouble GridKDTree::ComputeVolumeBy4Ppoints(const RDouble &x1, const RDouble &y1, const RDouble &z1,
    const RDouble &x2, const RDouble &y2, const RDouble &z2,
    const RDouble &x3, const RDouble &y3, const RDouble &z3,
    const RDouble &x4, const RDouble &y4, const RDouble &z4)
{
    RDouble xx1 = x2 - x1;
    RDouble yy1 = y2 - y1;
    RDouble zz1 = z2 - z1;
    RDouble xx2 = x3 - x1;
    RDouble yy2 = y3 - y1;
    RDouble zz2 = z3 - z1;
    RDouble xx3 = x4 - x1;
    RDouble yy3 = y4 - y1;
    RDouble zz3 = z4 - z1;

    RDouble volume = xx1 * yy2 * zz3+ xx2 * yy3 * zz1+ xx3 * yy1 * zz2
        - xx3 * yy2 * zz1- xx2 * yy1 * zz3- xx1 * yy3 * zz2;
    return volume + SMALL;
}

RDouble GetToleranceForOversetSearch()
{
    return GlobalDataBase::GetDoubleParaFromDB("toleranceForOversetSearch");
}

}