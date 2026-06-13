#include "DynamicGrid.h"
#include "GeometryUnit.h"
#include "HyList.h"

namespace PHSPACE
{
DynamicGrid::DynamicGrid(Grid *stationalGridIn)
{
    stationalGrid = stationalGridIn;

    nodeArray = 0;
    faceArray = 0;
    cellArray = 0;

    nTotalNode = 0;
    nTotalCell = 0;
    nBoundFace = 0;
    nTotalFace = 0;

    secondSegment = 0;
    rotateNodeIndex = 0;
}

DynamicGrid::~DynamicGrid()
{
    delete [] nodeArray;
    delete [] faceArray;
    delete [] cellArray;
    delete [] secondSegment;
}

void DynamicGrid::CalSpringK()
{
    RDouble dx, dy ,dz, dd; 
    double kk;

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        int nodeIndex1 = iNode;
        int nNodes = nodeArray[iNode].nodeNumberAround;
        nodeArray[iNode].node2nodeK = new double [nNodes];

        for (int jNode = 0; jNode < nNodes; ++ jNode)
        {
            int nodeIndex2 = nodeArray[iNode].node2node[jNode];

            dx = nodeArray[nodeIndex1].x - nodeArray[nodeIndex2].x;
            dy = nodeArray[nodeIndex1].y - nodeArray[nodeIndex2].y;
            dz = nodeArray[nodeIndex1].z - nodeArray[nodeIndex2].z;
            dd = sqrt(dx * dx + dy * dy + dz * dz);

            kk = static_cast <double> (1.0 / dd);
            nodeArray[iNode].node2nodeK[jNode] = kk;
        }
    }
}

void DynamicGrid::FindNode2Node()
{
    if (nodeArray[0].node2node != 0) return;

    HyList <int> *node2node = new HyList <int> [nTotalNode];
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        node2node[iNode].SetAverNnode(16);
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nNodes = faceArray[iFace].nodeNumber;
        int *face2node = faceArray[iFace].face2node;

        for (int iNode = 0; iNode < nNodes; ++ iNode)
        {
            int nodeIndex1 = face2node[iNode];
            int nodeIndex2 = face2node[(iNode + 1) % nNodes];

            node2node[nodeIndex1].insert(nodeIndex2);
            node2node[nodeIndex2].insert(nodeIndex1);
        }
    }

    int minsize = 10000, maxsize = 0, aversize = 0;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        int nNodes = node2node[iNode].size();
        nodeArray[iNode].nodeNumberAround = nNodes;

        nodeArray[iNode].node2node = new int[nNodes];
        for (int jNode = 0; jNode < nNodes; ++ jNode)
        {
            nodeArray[iNode].node2node[jNode] = node2node[iNode].GetData(jNode);
        }

        minsize = MIN(minsize, nNodes);
        maxsize = MAX(maxsize, nNodes);
        aversize += nNodes;
    }

    delete [] node2node;
    CalSpringK();
}

void DynamicGrid::ReconstrcutImplicitGeomeInfor()
{
    cellArray = new DYCell[nTotalCell];

    HyList <int> *cell2face = new HyList <int> [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2face[iCell].SetAverNnode(6);
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int leftCell = faceArray[iFace].leftCell;
        cell2face[leftCell].insert(iFace);
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int leftCell  = faceArray[iFace].leftCell;
        int rightCell = faceArray[iFace].rightCell;

        cell2face[leftCell ].insert(iFace);
        cell2face[rightCell].insert(iFace);
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nFaces = cell2face[iCell].size();
        cellArray[iCell].faceNumber = nFaces;

        cellArray[iCell].cell2face = new int[nFaces];
        for (int iFace = 0; iFace < nFaces; ++ iFace)
        {
            cellArray[iCell].cell2face[iFace] = cell2face[iCell].GetData(iFace);
        }
    }
    delete [] cell2face;

    FindNode2Node();
}

void DynamicGrid::SetSymmetryToZero()
{
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = faceArray[iFace].bcType;
        if (bcType == PHENGLEI::SYMMETRY)
        {
            int *face2node = faceArray[iFace].face2node;
            for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
            {
                int nodeIndex = face2node[iNode];
                nodeArray[nodeIndex].z    = 0.0;
                nodeArray[nodeIndex].znew = 0.0;
            }
        }
    }
}

void DynamicGrid::StationalData2DynamicData()
{
    UnstructGrid *grid = UnstructGridCast(stationalGrid);

    nTotalNode = grid->GetNTotalNode();
    nTotalFace = grid->GetNTotalFace();
    nTotalCell = grid->GetNTotalCell(); 
    nBoundFace = grid->GetNBoundFace();

    //! Node info.
    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    DYNode *nodeList = new DYNode[nTotalNode];
    this->nodeArray = nodeList;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        nodeList[iNode].x = x[iNode];
        nodeList[iNode].y = y[iNode];
        nodeList[iNode].z = z[iNode];

        nodeList[iNode].xnew = x[iNode];
        nodeList[iNode].ynew = y[iNode];
        nodeList[iNode].znew = z[iNode];
    }

    //! Face info.
    int *face2node           = grid->GetFace2Node();
    int *nodeNumOfEaceFace   = grid->GetNodeNumberOfEachFace();
    int *leftCellOfEachFace  = grid->GetLeftCellOfFace();
    int *rightCellOfEachFace = grid->GetRightCellOfFace();
    UnstructBCSet **bcr = grid->GetBCRecord();

    int count = 0;
    DYFace *faceList = new DYFace[nTotalFace];
    this->faceArray = faceList;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nNodes = nodeNumOfEaceFace[iFace];
        faceList[iFace].nodeNumber = nNodes;
        faceList[iFace].face2node  = new int[nNodes];

        for (int iNode = 0; iNode < nNodes; ++ iNode)
        {
            faceList[iFace].face2node[iNode] = face2node[count ++];
        }

        faceList[iFace].leftCell  = leftCellOfEachFace [iFace];
        faceList[iFace].rightCell = rightCellOfEachFace[iFace];

        if (iFace < nBoundFace)
        {
            faceList[iFace].bcType = bcr[iFace]->GetKey(); 
        }
        else
        {
            faceList[iFace].bcType = PHENGLEI::INTERIOR;
        }
    }

    //! Find the index of the rotate point of second segment.
    double secondSegmentRotatePostionX = 0.0, secondSegmentRotatePostionY = 0.0, secondSegmentRotatePostionZ = 0.0;
    int *secondSegment = new int[nTotalNode];
    SetField(secondSegment, 0, nTotalNode);
    this->secondSegment = secondSegment;

    SimpleBC *boundaryCondition = NULL;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bctype = faceArray[iFace].bcType;
        if (bctype != PHENGLEI::SOLID_SURFACE)
        {
            continue;
        }

        boundaryCondition = bcr[iFace]->GetBoundaryCondition();
        int isSecondSegment = 0;
        if (boundaryCondition->CheckParamData("isSecondSegment"))
        {
            boundaryCondition->GetParamData("isSecondSegment", &isSecondSegment, PHINT, 1);
        }

        if (!isSecondSegment)
        {
            continue;
        }

        int *face2nodeFaceA = faceArray[iFace].face2node;
        for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
        {
            int nodeIndex = face2nodeFaceA[iNode];
            secondSegment[nodeIndex] = 1;
        }

        boundaryCondition->GetParamData("secondSegmentRotatePostionX", &secondSegmentRotatePostionX, PHINT, 1);
        boundaryCondition->GetParamData("secondSegmentRotatePostionY", &secondSegmentRotatePostionY, PHINT, 1);
        boundaryCondition->GetParamData("secondSegmentRotatePostionZ", &secondSegmentRotatePostionZ, PHINT, 1);
    }

    rotateNodeIndex = -1;
    double minDist = 1000.0;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (secondSegment[iNode] != 1)
        {
            continue;
        }

        double dx, dy, dz;
        dx = nodeList[iNode].x - secondSegmentRotatePostionX;
        dy = nodeList[iNode].y - secondSegmentRotatePostionY;
        dz = nodeList[iNode].z - secondSegmentRotatePostionZ;

        double dd = sqrt(dx * dx + dy * dy + dz * dz);
        minDist = MIN(minDist, dd);
        if (dd <= 5.0e-07)
        {
            rotateNodeIndex = iNode;
            break;
        }

        if (rotateNodeIndex >= 0)
        {
            break;
        }
    }
}

}