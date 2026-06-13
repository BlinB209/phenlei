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
//! @file      Pre_GridCombine.h
//! @brief     Explain this file briefly.
//! @author    Xu Gang.

#pragma once
#include "PHMpi.h"
#include "Region.h"

using namespace std;

namespace PHSPACE
{

class CombinGrid
{
public:
    CombinGrid(Region *region, Grid **grid_in, int nZones_in);
    ~CombinGrid();
public:
    void RunCombinGrid();
    Grid ** GetGridAll() const;
    void DumpMergedGrid(const string &gridFileName, Grid **grids, int nBlocks);
    void DumpMergedGridP2P(const string &gridFileName, Grid **grids, int nBlocks);

private:
    void CheckifIndependentZoneExist();
    Grid * CreateCombinationGrid();
    void ResetCoordinate();
    void ResetNodeNumberOfEachFace();
    void ResetFace2Node();
    void ResetLeftRightCell();
    void MergeCell2Node();
    void MergeBCType();
    void DumpCombinationGrid();

    void SwapInterfaceData();
    void CompressInterfaceData(int zoneID, DataContainer *sendBuffer, int neighborZone);
    void DecompressInterfaceData(int zoneID, DataContainer *recieveBuffer, int neighborZone);
    
    void Init();
    void InitLocal2All();
    void InitNode2All();
    void InitFace2All();
    void InitCell2All();
    void InitIsBCFaceMerged();

    //! Establish the nodes' connection between zone and whole.
    void MergeNode2All();
    //void MergeNode2AllNew();

    //! Establish the cells' connection between zone and whole.
    void MergeCell2All();

    //! Establish the faces' connection between zone and whole.
    void MergeFace2All();

    //! Establish the interfaces' connection between zone and whole.
    void MergeInterface2All();

    //! Whether this node is boundary node.
    bool IsBcNode(int *p, int node_in);

    int * GetNode2All(int iZone) { return node2all[iZone]; }
    int * GetFace2All(int iZone) { return face2all[iZone]; }
    int * GetCell2All(int iZone) { return cell2all[iZone]; }
    //bool * GetIsBCFaceMerged(int iZone) { return isBCFaceMerged[iZone]; }

    int GetNTotalCell() const { return nTotalCell; }
    int GetNTotalNode() const { return nTotalNode; }
    int GetNTotalFace() const { return nTotalFace; }
    int GetNIFace() const { return nIFace; }

    vector<RDouble> GetGlobalx() { return Globalx; }
    vector<RDouble> GetGlobaly() { return Globaly; }
    vector<RDouble> GetGlobalz() { return Globalz; }

    int * GetNodeNumberOfAllFace() { return node_number_of_all_face; }
    vector <int> GetAllFace2Node() { return allface2node; }
    int * GetLeftCellOfAllFace() { return left_cell_of_all_face; }
    int * GetRightCellOfAllFace() { return right_cell_of_all_face; }

private:
    Region *region;

    //! Original multi-block grids.
    Grid **grid_in;

    //! Combined single block grid.
    UnstructGrid *grid_all;

    //! The number of Local zones.
    int numberOfLocalZones;

    //! The number of total cells.
    int nTotalCell;

    //! The number of total nodes.
    int nTotalNode;

    //! The number of total faces.
    int nTotalFace;

    //! The number of total interfaces.
    int nIFace;

    //! The nodes' connection between zone and whole.
    int **node2all;

    //! The faces' connection between zone and whole.
    int **face2all;

    //! The cells' connection between zone and whole.
    int **cell2all;

    //! The grid's information in each zone.
    int *nTotalNodeEachLocalZone, *nTotalFaceEachLocalZone, *nTotalCellEachLocalZone, *nIFaces;

    //! Global coordinates.
    vector<RDouble> Globalx, Globaly, Globalz;

    //! Left cell index of each face after combination.
    int *left_cell_of_all_face;

    //! Right cell index of each face after combination.
    int *right_cell_of_all_face;

    //! Number of nodes per face after combination.
    int *node_number_of_all_face;

    bool *isFaceMerged;
    int  *face2NewFaceafterReorder;
    int  *face2Interface;
    int  *face2InterfaceafterReorder;
    int  *interface2BoundaryFaceafterReorder;
    bool *isGlobalInterFace;

    //! Face to node: node index of each face after combination.
    vector <int> allface2node;

    set <int> zoneIDinCurrentProc;
    map <int, int> zoneIndexMap;
};

}