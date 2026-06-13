#include <iostream>
#include <Geo_UnstructMerge.h>
#include <iostream>
#include "Mesh_Deformation.h"
#include "Mesh_DeformationSPRING.h"
#include "Mesh_DeformationRBF.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "GeometryUnit.h"
#include "Geo_Grid.h"
#include "PHMpi.h"
#include "Glb_Dimension.h"
#include "Pre_HDF5File.h"
#include <algorithm> 
#include <string.h>
#include <Geo_SimpleBC.h>
#include "GridType.h"
#include <Region.h>
#include <Post_WriteTecplot.h>
#include <TK_Log.h>
using namespace std;

namespace PHSPACE
{
using namespace PHMPI;

OriginalUnstructGridMerge::OriginalUnstructGridMerge()
{
    nTotalZone                        = 0;
    nTotalBlock                       = 0;
    nDimensions                       = 0;
    nTotalNode                        = 0;
    nTotalFace                        = 0;
    nTotalCell                        = 0;
    coordinates                       = 0;
    nodeNumberOfEachFace              = 0;
    face2Node                         = 0;
    leftCellOfFace                    = 0;
    rightCellOfFace                   = 0;
    nodeNumberOfEachCell              = 0;
    cell2Node                         = 0;
    nIFace                            = 0;
    interfaceInfo                     = 0;
    nBoundFace                        = 0;
    bcInfo                            = 0;
    keyActiveOfCells                  = 0;
    mergegrid                         = 0;
    ordinaryGridIndexOfEachZoneGlobal = 0;
}

OriginalUnstructGridMerge::~OriginalUnstructGridMerge()
{
    nTotalZone = 0;
    if (mergegrid)
    {
        DelPointer(mergegrid);
        DelPointer(nTotalNode);
        DelPointer(nTotalFace);
        DelPointer(nTotalCell);
        for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
        {
            DelPointer(coordinates[iBlock]);
        }

        DelPointer(coordinates);

        DelPointer(nodeNumberOfEachFace);
        DelPointer(face2Node);
        DelPointer(leftCellOfFace);
        DelPointer(rightCellOfFace);
        DelPointer(nodeNumberOfEachCell);
        DelPointer(keyActiveOfCells);
        DelPointer(cell2Node);
        DelPointer(nIFace);
        DelPointer(interfaceInfo);
        DelPointer(nBoundFace);
        DelPointer(bcInfo);
    }
}

void OriginalUnstructGridMerge::SetOriginalGridOfBlock(vector < vector < Grid* > > originalgridofBlockIn, int nTotalBlockIn, int nTotalZone, int nDim)
{
    this->nTotalBlock = nTotalBlockIn;
    this->nTotalZone  = nTotalZone;
    originalgridofBlock.resize(nTotalBlock);
    nDimensions = nDim;
    for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
    {
        if (!originalgridofBlockIn[iBlock].size())
        {
            continue;
        }

        for (int iZone = 0; iZone < originalgridofBlockIn[iBlock].size(); ++ iZone)
        {
            originalgridofBlock[iBlock].push_back(UnstructGridCast(originalgridofBlockIn[iBlock][iZone]));
        }
    }
}

void OriginalUnstructGridMerge::SetOrdinaryGridIndexOfEachZoneGlobal(int *ordinaryGridIndexOfEachZoneGlobal)
{
    this->ordinaryGridIndexOfEachZoneGlobal = ordinaryGridIndexOfEachZoneGlobal;
}

void OriginalUnstructGridMerge::SetOriGridIDOfEachZoneGlobal(int blockID, int zoneID)
{
    this->oriGridIDOfEachZoneGlobal[blockID].push_back(zoneID);
}

void OriginalUnstructGridMerge::Initial()
{
    oriNodeIndex.resize(nTotalBlock);
    oriFaceIndex.resize(nTotalBlock);
    oriCellIndex.resize(nTotalBlock);

    nTotalNode           = new int [nTotalBlock];
    nTotalFace           = new int [nTotalBlock];
    nTotalCell           = new int [nTotalBlock];
    coordinates          = new RDouble ** [nTotalBlock];
    nodeNumberOfEachFace = new int * [nTotalBlock];
    face2Node            = new int * [nTotalBlock];
    leftCellOfFace       = new int * [nTotalBlock];
    rightCellOfFace      = new int * [nTotalBlock];
    nodeNumberOfEachCell = new int * [nTotalBlock];
    cell2Node            = new int * [nTotalBlock];
    nIFace               = new int [nTotalBlock];
    interfaceInfo        = new InterfaceInfo * [nTotalBlock];
    nBoundFace           = new int [nTotalBlock];
    dumpOriginalBCRegion.resize(nTotalBlock);
    bcInfo           = new UnstructBCSet ** [nTotalBlock];
    keyActiveOfCells = new int* [nTotalBlock];
    originalgridofBlock.resize(nTotalBlock);
    mergegrid = new UnstructGrid * [nTotalBlock];
    for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
    {
        coordinates[iBlock] = new RDouble * [3];
    }

    for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
    {
        if (!originalgridofBlock[iBlock].size())
        {
            continue;
        }

        nTotalNode[iBlock] = 0;
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            for (int iNode = 0; iNode < (originalgridofBlock[iBlock][iZone]->GetNTotalNode()); ++ iNode)
            {
                nTotalNode[iBlock] = max(nTotalNode[iBlock], (originalgridofBlock[iBlock][iZone]->GetOrdinaryNodeIndex())[iNode]);
            }
        }
        ++ (nTotalNode[iBlock]);

        nTotalFace[iBlock] = 0;
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            for (int iFace = 0; iFace < (originalgridofBlock[iBlock][iZone]->GetNTotalFace()); ++ iFace)
            {
                nTotalFace[iBlock] = max(nTotalFace[iBlock], (originalgridofBlock[iBlock][iZone]->GetOrdinaryFaceIndex())[iFace]);
            }
        }
        ++ (nTotalFace[iBlock]);

        nTotalCell[iBlock] = 0;
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            for (int iCell = 0; iCell < (originalgridofBlock[iBlock][iZone]->GetNTotalCell()); ++ iCell)
            {
                nTotalCell[iBlock] = max(nTotalCell[iBlock], (originalgridofBlock[iBlock][iZone]->GetOrdinaryCellIndex())[iCell]);
            }
        }
        ++ (nTotalCell[iBlock]);
    }

    oriGridIDOfEachZoneGlobal.resize(nTotalBlock);
    for (int iZone = 0; iZone < nTotalZone; ++iZone)
    {
        SetOriGridIDOfEachZoneGlobal(ordinaryGridIndexOfEachZoneGlobal[iZone], iZone);
    }
}

void OriginalUnstructGridMerge::MergeNode()
{
    for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
    {
        if (!originalgridofBlock[iBlock].size())
        {
            continue;
        }

        for (int iDim = 0; iDim < THREE_D; ++ iDim)
        {
            coordinates[iBlock][iDim] = new RDouble[nTotalNode[iBlock]];
        }

        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            for (int iNode = 0; iNode < (originalgridofBlock[iBlock][iZone]->GetNTotalNode()); ++ iNode)
            {
                int nodeID = (originalgridofBlock[iBlock][iZone]->GetOrdinaryNodeIndex())[iNode];
                RDouble connX, connY, connZ;
                connX = (originalgridofBlock[iBlock][iZone]->GetX())[iNode];
                connY = (originalgridofBlock[iBlock][iZone]->GetY())[iNode];
                connZ = (originalgridofBlock[iBlock][iZone]->GetZ())[iNode];
                coordinates[iBlock][ONE_D - 1]  [nodeID] = connX;
                coordinates[iBlock][TWO_D - 1]  [nodeID] = connY;
                coordinates[iBlock][THREE_D - 1][nodeID] = connZ;
            }
        }
    }

}

void OriginalUnstructGridMerge::MergeFace()
{
    for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
    {
        if (!originalgridofBlock[iBlock].size())
        {
            continue;
        }

        nodeNumberOfEachFace[iBlock] = new int [nTotalFace[iBlock]];
        vector<vector<int>>dumpFace2Node;
        dumpFace2Node.resize(nTotalFace[iBlock]);
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            for (int iFace = 0; iFace < (originalgridofBlock[iBlock][iZone]->GetNTotalFace()); ++ iFace)
            {
                int faceID = (originalgridofBlock[iBlock][iZone]->GetOrdinaryFaceIndex())[iFace];
                nodeNumberOfEachFace[iBlock][faceID] = (originalgridofBlock[iBlock][iZone]->GetNodeNumberOfEachFace())[iFace];
                dumpFace2Node[faceID].resize(nodeNumberOfEachFace[iBlock][faceID]);
                for (int iNode = 0; iNode < nodeNumberOfEachFace[iBlock][faceID]; ++ iNode)
                {
                    int oriFace2NodeNodeID = (originalgridofBlock[iBlock][iZone]->GetFace2NodeArray())[iFace][iNode];
                    dumpFace2Node[faceID][iNode] = (originalgridofBlock[iBlock][iZone]->GetOrdinaryNodeIndex())[oriFace2NodeNodeID];
                }
            }
        }
        int dataSize = 0;
        for (int iFace = 0; iFace < nTotalFace[iBlock]; ++ iFace)
        {
            dataSize = dataSize + nodeNumberOfEachFace[iBlock][iFace];
        }
        face2Node[iBlock] = new int [dataSize];
        int nodeCount = 0;
        for (int iFace = 0; iFace < nTotalFace[iBlock]; ++ iFace)
        {
            for (int iNode = 0; iNode < nodeNumberOfEachFace[iBlock][iFace]; ++ iNode)
            {
                face2Node[iBlock][nodeCount] = dumpFace2Node[iFace][iNode];
                ++ nodeCount;
            }
        }

        leftCellOfFace [iBlock] = new int [nTotalFace[iBlock]];
        rightCellOfFace[iBlock] = new int [nTotalFace[iBlock]];
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            for (int iFace = 0; iFace < (originalgridofBlock[iBlock][iZone]->GetNTotalFace()); ++ iFace)
            {
                int faceID         = (originalgridofBlock[iBlock][iZone]->GetOrdinaryFaceIndex())[iFace];
                int oriLeftCellID  = (originalgridofBlock[iBlock][iZone]->GetLeftCellOfFace())[iFace];
                int oriRightCellID = (originalgridofBlock[iBlock][iZone]->GetRightCellOfFace())[iFace];
                leftCellOfFace[iBlock][faceID] = (originalgridofBlock[iBlock][iZone]->GetOrdinaryCellIndex())[oriLeftCellID];

                if (oriRightCellID >= 0 && oriRightCellID < (originalgridofBlock[iBlock][iZone]->GetNTotalCell()))
                {
                    rightCellOfFace[iBlock][faceID] = (originalgridofBlock[iBlock][iZone]->GetOrdinaryCellIndex())[oriRightCellID];
                }
                else if (oriRightCellID < 0)
                {
                    rightCellOfFace[iBlock][faceID] = -1;
                }
                else if (oriRightCellID >= (originalgridofBlock[iBlock][iZone]->GetNTotalCell()))
                {
                    int isInterFace = 0;
                    int neighborZoneID, iFaceID, neighborZoneIFaceID, neighborZoneFaceID;
                    for (int iIFace = 0; iIFace < (originalgridofBlock[iBlock][iZone]->GetNIFace()); ++ iIFace)
                    {
                        if (iFace == (originalgridofBlock[iBlock][iZone]->GetInterfaceInfo()->GetInterFace2BoundaryFace())[iIFace])
                        {
                            isInterFace = 1;
                            neighborZoneID = (originalgridofBlock[iBlock][iZone]->GetInterfaceInfo()->GetInterFace2ZoneID())[iIFace];
                            iFaceID = iIFace;
                            if (ordinaryGridIndexOfEachZoneGlobal[neighborZoneID] != iBlock)
                            {
                                rightCellOfFace[iBlock][faceID] = nTotalCell[iBlock] + faceID;
                                continue;
                            }
                            else
                            {
                                int zoneID;
                                for (int iZoneFind = 0; iZoneFind < originalgridofBlock[iBlock].size(); ++iZoneFind)
                                {
                                    if (neighborZoneID == oriGridIDOfEachZoneGlobal[iBlock][iZoneFind])
                                    {
                                        zoneID = iZoneFind;
                                        break;
                                    }
                                }
                                neighborZoneIFaceID = (originalgridofBlock[iBlock][iZone]->GetInterfaceInfo()->GetInterFace2InterFaceID())[iFaceID];
                                neighborZoneFaceID = (originalgridofBlock[iBlock][zoneID]->GetInterfaceInfo()->GetInterFace2BoundaryFace())[neighborZoneIFaceID];
                                int oriNeighborZoneCellID = (originalgridofBlock[iBlock][zoneID]->GetLeftCellOfFace())[neighborZoneFaceID];
                                rightCellOfFace[iBlock][faceID] = (originalgridofBlock[iBlock][zoneID]->GetOrdinaryCellIndex())[oriNeighborZoneCellID];
                            }
                            continue;
                        }
                    }
                    if (!isInterFace)
                    {
                        rightCellOfFace[iBlock][faceID] = nTotalCell[iBlock] + faceID;
                    }
                }
            }
        }
    }
}

void OriginalUnstructGridMerge::MergeCell()
{
    for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
    {
        if (!originalgridofBlock[iBlock].size())
        {
            continue;
        }

        nodeNumberOfEachCell[iBlock] = new int [nTotalCell[iBlock]];
        vector<vector<int>>dumpCell2Node;
        dumpCell2Node.resize(nTotalCell[iBlock]);
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            for (int iCell = 0; iCell < (originalgridofBlock[iBlock][iZone]->GetNTotalCell()); ++ iCell)
            {
                int cellID = (originalgridofBlock[iBlock][iZone]->GetOrdinaryCellIndex())[iCell];
                nodeNumberOfEachCell[iBlock][cellID] = (originalgridofBlock[iBlock][iZone]->GetNodeNumberOfEachCell())[iCell];
                dumpCell2Node[cellID].resize(nodeNumberOfEachCell[iBlock][cellID]);
                for (int iNode = 0; iNode < nodeNumberOfEachCell[iBlock][cellID]; ++ iNode)
                {
                    int oriCell2NodeNodeID = (originalgridofBlock[iBlock][iZone]->GetCell2NodeArray())[iCell][iNode];
                    dumpCell2Node[cellID][iNode] = (originalgridofBlock[iBlock][iZone]->GetOrdinaryNodeIndex())[oriCell2NodeNodeID];
                }
            }
        }
        int dataSize = 0;
        for (int iCell = 0; iCell < nTotalCell[iBlock]; ++ iCell)
        {
            dataSize = dataSize + nodeNumberOfEachCell[iBlock][iCell];
        }
        cell2Node[iBlock] = new int [dataSize];
        int nodeCount = 0;
        for (int iCell = 0; iCell < nTotalCell[iBlock]; ++ iCell)
        {
            for (int iNode = 0; iNode < nodeNumberOfEachCell[iBlock][iCell]; ++ iNode)
            {
                cell2Node[iBlock][nodeCount] = dumpCell2Node[iCell][iNode];
                ++nodeCount;
            }
        }
        keyActiveOfCells[iBlock] = new int [nTotalCell[iBlock]];
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++iZone)
        {
            for (int iCell = 0; iCell < (originalgridofBlock[iBlock][iZone]->GetNTotalCell()); ++iCell)
            {
                int cellID = (originalgridofBlock[iBlock][iZone]->GetOrdinaryCellIndex())[iCell];
                keyActiveOfCells[iBlock][cellID] = (originalgridofBlock[iBlock][iZone]->GetBlankIndex())[iCell];
            }
        }
    }
}

void OriginalUnstructGridMerge::MergeBoundary()
{
    for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
    {
        if (!originalgridofBlock[iBlock].size())
        {
            continue;
        }

        nIFace[iBlock] = 0;
        int nRedundantIFace = 0;
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            InterfaceInfo *dumpInterfaceInfo = originalgridofBlock[iBlock][iZone]->GetInterfaceInfo();
            int nIFaceOfEachZone = originalgridofBlock[iBlock][iZone]->GetNIFace();
            for (int iIFace = 0; iIFace < nIFaceOfEachZone; ++ iIFace)
            {
                int neightborZoneID = (dumpInterfaceInfo->GetInterFace2ZoneID())[iIFace];
                if (ordinaryGridIndexOfEachZoneGlobal[neightborZoneID] == iBlock)
                {
                    ++ nRedundantIFace;
                }
                else
                {
                    ++ nIFace[iBlock];
                }
            }
        }

        nBoundFace[iBlock] = 0;
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            nBoundFace[iBlock] = nBoundFace[iBlock] + (originalgridofBlock[iBlock][iZone]->GetNBoundFace());
        }
        nBoundFace[iBlock] = nBoundFace[iBlock] - nRedundantIFace;
        bcInfo[iBlock] = new UnstructBCSet * [nBoundFace[iBlock]];
        for (int iBCFace = 0; iBCFace < nBoundFace[iBlock]; ++ iBCFace)
        {
            bcInfo[iBlock][iBCFace] = new UnstructBCSet();
        }
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++iZone)
        {
            for (int iFace = 0; iFace < (originalgridofBlock[iBlock][iZone]->GetNBoundFace()); ++ iFace)
            {
                int faceID    = (originalgridofBlock[iBlock][iZone]->GetOrdinaryFaceIndex())[iFace];
                if (faceID < nBoundFace[iBlock])
                {
                    int    key    = (originalgridofBlock[iBlock][iZone]->GetBCRecord())[iFace]->GetKey();
                    string bcName = (originalgridofBlock[iBlock][iZone]->GetBCRecord())[iFace]->GetBCName();
                    bcInfo[iBlock][faceID]->SetKey(key);
                    bcInfo[iBlock][faceID]->SetBCName(bcName);
                }
            }
        }
    }
}

void OriginalUnstructGridMerge::ConstructGrid()
{
    for (int iBlock = 0; iBlock < nTotalBlock; ++iBlock)
    {
        if (!originalgridofBlock[iBlock].size())
        {
            continue;
        }

        GridID *index = new GridID(iBlock);
        mergegrid[iBlock] = UnstructGridCast(CreateGridGeneral(UNSTRUCTGRID, index, 0, 3));
        mergegrid[iBlock]->SetDim(nDimensions);
        mergegrid[iBlock]->SetNTotalNode(nTotalNode[iBlock]);
        mergegrid[iBlock]->SetNTotalFace(nTotalFace[iBlock]);
        mergegrid[iBlock]->SetNTotalCell(nTotalCell[iBlock]);
        mergegrid[iBlock]->SetX(coordinates[iBlock][ONE_D - 1]);
        mergegrid[iBlock]->SetY(coordinates[iBlock][TWO_D - 1]);
        mergegrid[iBlock]->SetZ(coordinates[iBlock][THREE_D - 1]);
        mergegrid[iBlock]->RotateAboutAxis();
        mergegrid[iBlock]->ComputeMinMaxBox();
        mergegrid[iBlock]->SetNodeNumberOfEachFace(nodeNumberOfEachFace[iBlock]);
        mergegrid[iBlock]->SetFace2Node(face2Node[iBlock]);
        mergegrid[iBlock]->SetLeftCellOfFace(leftCellOfFace[iBlock]);
        mergegrid[iBlock]->SetRightCellOfFace(rightCellOfFace[iBlock]);
        mergegrid[iBlock]->SetNodeNumberOfEachCell(nodeNumberOfEachCell[iBlock]);
        mergegrid[iBlock]->SetCell2Node(cell2Node[iBlock]);
        mergegrid[iBlock]->SetNIFace(nIFace[iBlock]);
        mergegrid[iBlock]->SetNBoundFace(nBoundFace[iBlock]);
        mergegrid[iBlock]->SetBCRecord(bcInfo[iBlock]);
        mergegrid[iBlock]->ConstructBCRegion();
        mergegrid[iBlock]->SpecifyRightCellofBC();
        mergegrid[iBlock]->SetBlankIndex(keyActiveOfCells[iBlock]);
    }
}

void OriginalUnstructGridMerge::Run()
{
    Initial();
    MergeNode();
    MergeFace();
    MergeCell();
    MergeBoundary();
    ConstructGrid();
}

void OriginalUnstructGridMerge::GetMergeGrid(UnstructGrid **mergeGridOut)
{
    for (int iBlock = 0; iBlock < nTotalBlock; ++ iBlock)
    {
        mergeGridOut[iBlock] = mergegrid[iBlock];
    }
}

}