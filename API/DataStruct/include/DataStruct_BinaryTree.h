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
//! @file      DataStruct_BinaryTree.h
//! @brief     Plant binary tree for overset grid.
//! @author    Guo Yongheng.

#include <vector>
using namespace std;
namespace PHSPACE
{
struct DataStuct_BinaryNode
{
    int inf, sup, zoneLabel;
    DataStuct_BinaryNode *leftChild, *rightChild;
};

class DataStruct_BinaryTree
{
private:
    vector<int> startPointLabel;
    int numberOfZones, numberOfPoints;
    DataStuct_BinaryNode *root;

public:
    //!
    DataStruct_BinaryTree(int *startPointLabel, int nZones, int numberOfPoints);

    //!
    ~DataStruct_BinaryTree();

public:
    //!
    int ComputeZoneLabel(int globalPointIndex);

private:
    int ComputeZoneLabel(DataStuct_BinaryNode *binaryNode, int globalPointIndex);
    void CopyOriginalVector(int *startPointLabel, int nZones);
    DataStuct_BinaryNode * Create(DataStuct_BinaryNode *binaryNode, int leftLabel, int rightLabel);
    void Release(DataStuct_BinaryNode *binaryNode);
};

}