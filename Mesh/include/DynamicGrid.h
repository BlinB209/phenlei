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
//! @file      DynamicGrid.h
//! @brief     Dynamic grid structures, for arbitrary mixed cells.
//! @author    Bell, Baka.

#pragma once
#include "Region.h"
#include "Geo_Grid.h"

using namespace std;
namespace PHSPACE
{
class DYCell;
class DYFace;
class DYNode;

class DynamicGrid
{
private:
    //! Original grid.
    Grid *stationalGrid;

    //! Number of total nodes.
    int nTotalNode;

    //! The number of Total Faces.
    int nTotalFace;

    //! The number of Total Cells.
    int nTotalCell;

    //! The number of boundary faces which include interface faces.
    int nBoundFace;

    //! Node information.
    DYNode *nodeArray;

    //! Face information.
    DYFace *faceArray;

    //! Cell information.
    DYCell *cellArray;

    //! Location of each node. When node on second segment, value = 1.
    int *secondSegment;

    //! The node index which is rotate center of second segment.
    int rotateNodeIndex;

public:
    //! @param[in] stationalGridIn  Original grid.
    DynamicGrid(Grid *stationalGridIn);
    ~DynamicGrid();

public:
    int GetNTotalNode() { return this->nTotalNode; }
    int GetNTotalFace() { return this->nTotalFace; }
    int GetNTotalCell() { return this->nTotalCell; }
    int GetNBoundFace() { return this->nBoundFace; }
    int GetRotateNodeIndex() { return this->rotateNodeIndex; }

    DYNode * GetNodeArray() { return this->nodeArray; }
    DYFace * GetFaceArray() { return this->faceArray; }
    DYCell * GetCellArray() { return this->cellArray; }

    int * GetSecondSegment() { return this->secondSegment; }
    void SetSecondSegment(int *secondSegmentIn) { this->secondSegment = secondSegmentIn; }

public:
    //! Computer spring coefficient.
    void CalSpringK();

    //! Construct node topo.
    void FindNode2Node();

    //! Construct cell topo.
    void ReconstrcutImplicitGeomeInfor();

    //! Distill grid information from original grid.
    void StationalData2DynamicData();

    void SetSymmetryToZero();
};

const int STATIC_POINT        = 0;
const int DYNAMIC_POINT_WALL  = 1;
const int DYNAMIC_POINT_FIELD = 2;

}
