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
//! @file      GridView.h
//! @brief     Show relationship between different overset grids.
//! @author    Guo Yongheng.

#pragma once
#include <fstream>
#include <vector>
#include "Precision.h"
using namespace std;

namespace PHSPACE
{
const int WHITE = 0, BLUE = 1, RED = 2, PURPLE = 3;
const int ASPECTX = 0, ASPECTY = 1, ASPECTZ = 2;
const int EXTERNAL_BC = 9;

class PHSimpleFace
{
public:
    PHSimpleFace();
    ~PHSimpleFace();
public:
    void SetKeyParameter(int *keyNumber);
    PHSimpleFace * GetNextPointer();
    void CreateNewSimpleFace();
    void ShowSimpleFace();
    int GetKeyParameter();
    int * GetLocalParameter();
protected:
    int keyNumber[7];
    PHSimpleFace *simpleFaceNext;
};

class BasalGrid
{
public:
    BasalGrid();
    ~BasalGrid();
public:
    void GenerateBasalSpace(int iDimension, int jDimension, int kDimension);
    void SetCoordinate(int aspect, int i, int j, int k, RDouble coordinate);
    int GetIDimension();
    int GetJDimension();
    int GetKDimension();
    RDouble GetCoordinate(int aspect, int i, int j, int k);
    RDouble * GetX() { return xCoordinate; }
    RDouble * GetY() { return yCoordinate; }
    RDouble * GetZ() { return zCoordinate; }
protected:
    int iDimension, jDimension, kDimension;
    int pointNumber;
    RDouble *xCoordinate;
    RDouble *yCoordinate;
    RDouble *zCoordinate;
};

class BoundaryCluster
{
public:
    BoundaryCluster();
    ~BoundaryCluster();
public:
    void SetBlockLabel(int blockLabel);
    void SetBlockDimension(int iDimension, int jDimension, int kDimension);
    void SetProcessLabel(int processLabel);
    void SetOverallFaceNumber(int overallFaceNumber);
    void SetExternalFaceNumber(int externalFaceNumber);
    int GetOverallFaceNumber();
    int GetProcessLabel();
    int GetExternalFaceNumber();
    vector< vector< int > > ExtractPhysicalBoundaryIndex();
public:
    void SetKeyNumberPointer(int *keyNumber);
    void CreateNewSimpleFace(int iFace);
    void DeleteLinkedTable();
public:
    PHSimpleFace * GetSimpleFaceHead();
    void RecordSolidFaceInformation();
    void OutputHoleSurfaceDimension(ofstream &outfile);
    void OutputHoleSurfaceGrid(ofstream &outfile, BasalGrid *basalGrid);
    int GetNumberOfSolidFaces() { return numberOfSolidFaces; }
    void ShowBoundaryCluster();
    void SetSlaveBlockList();
    int * GetSlaveBlockList();
protected:
    int iDimension, jDimension, kDimension;
    int overallFaceNumber, externalFaceNumber;
    int blockLabel, processLabel;
    PHSimpleFace *simpleFaceHead;
    PHSimpleFace *simpleFaceGuard;
    int *keyNumber;
    int *slaveBlockList;

    int numberOfSolidFaces;
    vector< int > ist, ied, jst, jed, kst, ked;
};

class LinkedBlocks
{
public:
    LinkedBlocks();
    ~LinkedBlocks();
public:
    void SetBlockLabel(int blockLabel);
    int GetBlockLabel();
    LinkedBlocks * GetNextPointer();
    void CreateNewLinkedBlock();
protected:
    int blockLabel;
    LinkedBlocks *linkedBlocksNext;
};

class BlockCluster
{
public:
    BlockCluster();
    ~BlockCluster();
public:
    void CreatNewLinkedBlock();
    void DeleteLinkedTable();
    void SetBlockLabel(int blockLabel);
    void ShowLinkedBlocks();
    void SetListColor(int listColor);
    void SetGroupLabel(int groupLabel);
    void SetNewBlockLabel(int newBlockLabel);
    int GetListColor();
    int GetGroupLabel();
    int GetNewBlockLabel();
    LinkedBlocks * GetLinkedBlocksHead();
protected:
    LinkedBlocks *linkedBlocksHead;
    LinkedBlocks *linkedBlocksGuard;
    int listColor, groupLabel, newBlockLabel;
};

class BlockGroupManager
{
public:
    BlockGroupManager();
    ~BlockGroupManager();
public:
    void Run();
protected:
    void InputOriginalGrid();
    void InputOriginalBoundary();
    void ShowOriginalBoundary();
    void GatherLinkedBlocks();
    void ShowLinkedBlocks();
    void MarkerBlocksRelation();
    void PermuteBlockLabel();
    void SeparateGrid();
    void OutputNewBlockGroup(int iGroup, ofstream &myFile);
    int SearchColorfulList(int spectrum);
    void OutputNewSurfaceGrid();
    void OutputNewGrid();
    void OutputNewBoundary();
    void OutputStandardBoundary();
    void OutputZoneInverseMapping();
    void OutputHoleGroup(int iGroup);
protected:
    int globalBlockNumber, processNumber, groupNumber, fileNumber;
    int originalExternalLabel;
    int keyNumber[7];

    bool readProc;

    BoundaryCluster *boundaryCluster;
    BlockCluster *blockCluster;
    BasalGrid *basalGrid;
    int *InverseMapping;
    int *blockStart;
    int *blockEnd;
};

}

