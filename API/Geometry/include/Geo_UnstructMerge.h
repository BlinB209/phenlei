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
//! @file      Geo_UnstructMerge.h
//! @brief     It is the class 'OriginalUnstructGridMerge', which is type of Unstructured geometry grid operation for grid merging.
//! @author    Xu Fengding, Xu Gang .

#pragma once
#include <iostream>
#include "Mesh_Deformation.h"
#include "Mesh_DeformationSPRING.h"
#include "Mesh_DeformationRBF.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "GeometryUnit.h"
#include "Geo_Grid.h"
#include "Geo_UnstructBC.h"
#include <vector>
using namespace std;
namespace PHSPACE
{
class OriginalUnstructGridMerge
{
private:
    int                                 nTotalBlock;
    int                                 nTotalZone;
    int                                 nDimensions;
    int                                 *nTotalNode;
    int                                 *nTotalFace;
    int                                 *nTotalCell;
    RDouble                             ***coordinates;
    int                                 **nodeNumberOfEachFace;
    int                                 **face2Node;
    int                                 **leftCellOfFace;
    int                                 **rightCellOfFace;
    int                                 **nodeNumberOfEachCell;
    int                                 **cell2Node;
    int                                 *nIFace;
    InterfaceInfo                       **interfaceInfo;
    int                                 *nBoundFace;
    vector < vector < UnstructBC* > >   dumpOriginalBCRegion;
    UnstructBCSet                       ***bcInfo;
    int                                 **keyActiveOfCells;
    vector < vector < UnstructGrid* > > originalgridofBlock;
    UnstructGrid                        **mergegrid;
    vector < vector < int > >           oriNodeIndex;
    vector < vector < int > >           oriFaceIndex;
    vector < vector < int > >           oriCellIndex;
    int                                 *ordinaryGridIndexOfEachZoneGlobal;
    vector < vector < int > >           oriGridIDOfEachZoneGlobal;
public:
    OriginalUnstructGridMerge();

    ~OriginalUnstructGridMerge();

    //! Pass in the mesh that needs to be merged
    void SetOriginalGridOfBlock(vector < vector < Grid* > > originalgridofBlockIn, int nTotalblock, int nTotalZone, int nDim);

    void SetOrdinaryGridIndexOfEachZoneGlobal(int *ordinaryGridIndexOfEachZoneGlobal);

    void SetOriGridIDOfEachZoneGlobal(int blockID, int zoneID);

    //! Initial all information need to merge.
    void Initial();

    void MergeBoundary();

    void MergeCell();

    void MergeFace();

    void MergeNode();

    //! Compute the mergegrid.
    void Run();

    //! Build a grid based on grid information.
    void ConstructGrid();

    //! Return the mergegrid.
    void GetMergeGrid(UnstructGrid **mergeGridOut);

    //! Return nTotalblock.
    int GetnTotalblock() { return nTotalBlock; };
};
}