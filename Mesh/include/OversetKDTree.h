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
//! @file      OversetKDTree.h
//! @brief     Unstructured grid overset KDT tree.
//! @author    Chang Xinhua, He Kun, Zhang Laiping.

#pragma once
#include "TypeDefine.h"
#include "DataContainer.h"
#include "Geo_Grid.h"
#include "Glb_Dimension.h"
#include "Region.h"
namespace PHSPACE
{
class BasicKDTree;
class GridKDTree;

class BasicKDTree
{
    typedef int (BasicKDTree:: *PartitionMethod)(int, int, int);
public:
    BasicKDTree(RDouble **dataIn,
        int         numberOfSearchElementsIn,
        int         numberOfSearchDirectionIn);
    ~BasicKDTree();
private:
    PartitionMethod partitionMethod;

    int numberOfSearchDirection;
    int numberOfSearchElements;
    RDouble **data;

    int *elementOrder;

    int *leftSonIndex;
    int *rightSonIndex;
    int rootIndex;
public:
    void BuildBasicKDTree();
private:
    void SwapElementIndex(int indexA, int indexB);
    int  PartitionSinglePivotDoubleEndScanQuickSort(int iStart, int iFinal, int iLayer);
    int  PartitionSinglePivotPreConditionDoubleEndScanQuickSort(int iStart, int iFinal, int iLayer);
    int  PartitionSinglePivotForwardScanQuickSort(int iPivot, int iFinal, int iLayer);
    int  PartitionSinglePivotPreConditionForwardScanQuickSort(int iPivot, int iFinal, int iLayer);
    int  FindHalfSinglePivotQuickSort(int iStart, int iFinal, int iLayer);
    void TreeSortSinglePivotQuickSort(int iStart, int iFinal, int iLayer);
public:
    int   GetRootIndex() { return this->rootIndex; }
    int *GetLeftSonIndex() { return this->leftSonIndex; }
    int *GetRightSonIndex() { return this->rightSonIndex; }
};

class GridKDTree
{
    typedef bool (GridKDTree:: *InCellCheckMethod)(RDouble *, int, int &);
public:
    GridKDTree(UnstructGrid *GridIn,
        int *keySearchOfCellsIn = NULL,
        RDouble *minBoxOfZoneIn = NULL,
        RDouble *maxBoxOfZoneIn = NULL);
    ~GridKDTree();
private:
    int geometricDimension;
    UnstructGrid *grid;
    BasicKDTree *basicKDTree;
    InCellCheckMethod           inCellCheckMethod;
    int *leftSonIndex;
    int *rightSonIndex;
    int                         rootIndex;

    int                         numberOfSearchDirection;
    int                         numberOfSearchCells;
    vector< int >               searchCellIndex;

    RDouble *minBoxOfZone;
    RDouble *maxBoxOfZone;
    RDouble **cellBox;

    bool                        ifAuxiliaryGrid;
private:
    void ComputeGlobalBox(RDouble tolerance);
    void ComputeCellBox(RDouble tolerance);
public:
    void EnlargeCellBox(const int &cellIndex, const RDouble &eps);
    int  SearchNodeInKDTree(RDouble *nodeCoordinateIn);
private:
    void SearchNodeInKDTree(RDouble *nodeCoordinateIn, int rootIndexIn, int &donorCellIndexIn, int iLayerIn);
private:
    bool InCellCheck2D(RDouble *nodeCoordinateIn, int rootIndexIn, int &donorCellIndexIn);
    bool InCellCheck3D(RDouble *nodeCoordinateIn, int rootIndexIn, int &donorCellIndexIn);
private:
    bool     In4PointsCell(const RDouble &x1, const RDouble &y1, const RDouble &z1,
        const RDouble &x2, const RDouble &y2, const RDouble &z2,
        const RDouble &x3, const RDouble &y3, const RDouble &z3,
        const RDouble &x4, const RDouble &y4, const RDouble &z4,
        RDouble *nodeCoordinateIn, const RDouble &eps);
    RDouble ComputeVolumeBy4Ppoints(const RDouble &x1, const RDouble &y1, const RDouble &z1,
        const RDouble &x2, const RDouble &y2, const RDouble &z2,
        const RDouble &x3, const RDouble &y3, const RDouble &z3,
        const RDouble &x4, const RDouble &y4, const RDouble &z4);
private:
    inline bool IfNodeIsInCellBox(RDouble *nodeCoordinateIn, int rootIndexIn)
    {
        for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
        {
            RDouble &nodeCoordinateInThisDimension = nodeCoordinateIn[iDimension];
            if (nodeCoordinateInThisDimension < cellBox[iDimension][rootIndexIn]) return false;
            if (nodeCoordinateInThisDimension > cellBox[iDimension + geometricDimension][rootIndexIn]) return false;
        }
        return true;
    }

    inline bool IfInLeftHalfTree(RDouble *nodeCoordinateIn, int rootIndexIn, int iDirectionIn)
    {
        if (iDirectionIn < geometricDimension)
        {
            if (cellBox[iDirectionIn][rootIndexIn] > nodeCoordinateIn[iDirectionIn]) return true;
        }
        return false;
    }

    inline bool IfInRightHalfTree(RDouble *nodeCoordinateIn, int rootIndexIn, int iDirectionIn)
    {
        if (iDirectionIn >= geometricDimension)
        {
            if (cellBox[iDirectionIn][rootIndexIn] < nodeCoordinateIn[iDirectionIn - geometricDimension]) return true;
        }
        return false;
    }
};

inline bool IfNodeIsInBox(RDouble *nodeCoordinateIn, RDouble *minBoxIn, RDouble *maxBoxIn, int geometricDimensionIn)
{
    for (int iDimension = 0; iDimension < geometricDimensionIn; ++ iDimension)
    {
        RDouble &nodeCoordinateInThisDimension = nodeCoordinateIn[iDimension];
        if (nodeCoordinateInThisDimension < minBoxIn[iDimension]) return false;
        if (nodeCoordinateInThisDimension > maxBoxIn[iDimension]) return false;
    }
    return true;
};

inline bool IfNodeIsInBox(vector< RDouble * > &nodeCoordinateIn, int nodeIndexIn, RDouble *minBoxIn, RDouble *maxBoxIn, int geometricDimensionIn)
{
    for (int iDimension = 0; iDimension < geometricDimensionIn; ++ iDimension)
    {
        RDouble &nodeCoordinateInThisDimension = nodeCoordinateIn[iDimension][nodeIndexIn];
        if (nodeCoordinateInThisDimension < minBoxIn[iDimension]) return false;
        if (nodeCoordinateInThisDimension > maxBoxIn[iDimension]) return false;
    }
    return true;
};

RDouble GetToleranceForOversetSearch();

}