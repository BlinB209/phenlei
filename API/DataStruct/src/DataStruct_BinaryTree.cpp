#include "DataStruct_BinaryTree.h"
using namespace std;

namespace PHSPACE
{
DataStruct_BinaryTree::DataStruct_BinaryTree(int *startPointLabel, int nZones, int numberOfPoints)
{
    CopyOriginalVector(startPointLabel, nZones);
    this->numberOfPoints = numberOfPoints;
    root = Create(root, 0, nZones-1);
    return;
}

DataStruct_BinaryTree::~DataStruct_BinaryTree()
{
    Release(root);
    return;
}

int DataStruct_BinaryTree::ComputeZoneLabel(int globalPointIndex)
{
    return ComputeZoneLabel(root, globalPointIndex);
}

int DataStruct_BinaryTree::ComputeZoneLabel(DataStuct_BinaryNode *binaryNode, int globalPointIndex)
{
    if (globalPointIndex < binaryNode->inf)
    {
        return ComputeZoneLabel(binaryNode->leftChild, globalPointIndex);
    }
    else if (globalPointIndex > binaryNode->sup)
    {
        return ComputeZoneLabel(binaryNode->rightChild, globalPointIndex);
    }
    else
    {
        return binaryNode->zoneLabel;
    }
}

void DataStruct_BinaryTree::CopyOriginalVector(int *startPointLabel, int nZones)
{
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        this->startPointLabel.push_back(startPointLabel[iZone]);
    }
    this->numberOfZones = nZones;
    return;
}

DataStuct_BinaryNode * DataStruct_BinaryTree::Create(DataStuct_BinaryNode *binaryNode, int leftLabel, int rightLabel)
{
    if (leftLabel > rightLabel) return 0;

    int middleLabel = (leftLabel + rightLabel) / 2;
    binaryNode = new DataStuct_BinaryNode();
    binaryNode->inf = startPointLabel[middleLabel];
    if (middleLabel < numberOfZones-1)
    {
        binaryNode->sup = startPointLabel[middleLabel+1] - 1;
    }
    else
    {
        binaryNode->sup = numberOfPoints - 1;
    }
    binaryNode->zoneLabel = middleLabel;

    binaryNode->leftChild  = Create(binaryNode->leftChild, leftLabel, middleLabel-1);
    binaryNode->rightChild = Create(binaryNode->rightChild, middleLabel+1, rightLabel);

    return binaryNode;
}

void DataStruct_BinaryTree::Release(DataStuct_BinaryNode *binaryNode)
{
    if (!binaryNode) return;

    Release(binaryNode->leftChild);
    Release(binaryNode->rightChild);
    delete binaryNode;
    return;
}

}