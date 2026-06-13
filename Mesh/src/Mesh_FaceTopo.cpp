#include "Mesh_FaceTopo.h"
#include "Geo_Element.h"
#include "TK_Exit.h"

namespace PHSPACE
{
Mesh_FaceTopo::Mesh_FaceTopo()
{

}

Mesh_FaceTopo::~Mesh_FaceTopo()
{

}

int Mesh_FaceTopo::GetNTotalFaces()
{
    return static_cast<int>(faceType.size());
}

int Mesh_FaceTopo::GetNumberOfTotalEdgesWithMiddleNode()
{
    return static_cast<int>(middleNodeOfEdge.size());
}

int Mesh_FaceTopo::GetNBoundaryFaces()
{
    int nTotalFaces = this->GetNTotalFaces();
    int nBoundaryFaces = 0;

    for (int iFace = 0; iFace < nTotalFaces; ++ iFace)
    {
        int rightCell = rightCellIndex[iFace];
        bool isBoundaryFace = (childrenFaceIndex[iFace].size() == 0);
        isBoundaryFace = isBoundaryFace && PHSPACE::IsBoundaryFace(rightCell);

        if (isBoundaryFace)
        {
            ++ nBoundaryFaces;
        }
    }

    return nBoundaryFaces;
}

int Mesh_FaceTopo::FindFace(MultiIndex face)
{
    set< MultiIndex >::iterator iter = faceList.find(face);
    if (iter == faceList.end())
    {
        return PHSPACE::INVALID_INDEX;
    }
    return iter->id;
}

int Mesh_FaceTopo::FindFace(PHVectorInt1D face2Node)
{
    uint_t nodeNumber = face2Node.size();

    PHVectorInt1D faceOfSort = face2Node;
    sort(faceOfSort.begin(), faceOfSort.end());

    MultiIndex faceOfMultiIndex(static_cast<int>(nodeNumber), GetNTotalFaces());
    faceOfMultiIndex.SetData(&faceOfSort[0]);

    return FindFace(faceOfMultiIndex);
}

int Mesh_FaceTopo::FindCompositeFace(MultiIndex compositeFace)
{
    set< MultiIndex >::iterator iter = compositeFacesList.find(compositeFace);
    if (iter == compositeFacesList.end())
    {
        return PHSPACE::INVALID_INDEX;
    }
    return iter->id;
}

int Mesh_FaceTopo::FindEdgeWithMiddleNode(MultiIndex edgeWithMiddleNode)
{
    set< MultiIndex >::iterator iter = edgesWithMiddleNode.find(edgeWithMiddleNode);
    if (iter == edgesWithMiddleNode.end())
    {
        return INVALID_INDEX;
    }
    return iter->id;
}

void Mesh_FaceTopo::InsertFace(MultiIndex face)
{
    faceList.insert(face);
}

void Mesh_FaceTopo::InsertCompositeFace(MultiIndex compositeface)
{
    compositeFacesList.insert(compositeface);
}

void Mesh_FaceTopo::InsertEdgeWithMiddleNode(MultiIndex edgeWithMiddleNode)
{
    edgesWithMiddleNode.insert(edgeWithMiddleNode);
}

int Mesh_FaceTopo::AddFace(PHVectorInt1D &face2NodeIn, int faceTypeIn, int parentFaceIndex)
{
    uint_t nodeNumber = face2NodeIn.size();

    PHVectorInt1D sortedFaceNodeIndex = face2NodeIn;
    sort(sortedFaceNodeIndex.begin(), sortedFaceNodeIndex.end());

    MultiIndex faceOfMultiIndex(static_cast<int>(nodeNumber), GetNTotalFaces());
    faceOfMultiIndex.SetData(&sortedFaceNodeIndex[0]);

    int faceIndex = FindFace(faceOfMultiIndex);
    if (faceIndex == PHSPACE::INVALID_INDEX)
    {
        InsertFace(faceOfMultiIndex);
        faceIndex = this->GetNTotalFaces();

        ResizeFaceNumber(faceIndex + 1);
        SetDefaultFaceInformation(faceIndex, faceTypeIn, face2NodeIn, parentFaceIndex);
        AddFaceHierachical(faceIndex, parentFaceIndex);
    }
    return faceIndex;
}

void Mesh_FaceTopo::ResizeFaceNumber(int faceNumber)
{
    faceType             .resize(faceNumber);
    face2Node            .resize(faceNumber);
    leftCellIndex        .resize(faceNumber);
    rightCellIndex       .resize(faceNumber);
    parentFaceIndex      .resize(faceNumber);
    childrenFaceIndex    .resize(faceNumber);
    faceBoundaryCondition.resize(faceNumber);
}

void Mesh_FaceTopo::SetDefaultFaceInformation(int faceIndex, int faceTypeIn, PHVectorInt1D &face2NodeIn, int parentFaceIndex)
{
    this->leftCellIndex[faceIndex] = PHSPACE::INVALID_INDEX;
    this->rightCellIndex[faceIndex] = PHSPACE::INVALID_INDEX;

    SimpleBC *parentBC = this->faceBoundaryCondition[parentFaceIndex];
    this->faceBoundaryCondition[faceIndex] = parentBC;

    this->faceType[faceIndex] = faceTypeIn;
    this->face2Node[faceIndex] = face2NodeIn;
}

void Mesh_FaceTopo::AddFaceHierachical(int faceIndex, int parentFaceIndexIn)
{
    this->parentFaceIndex[faceIndex] = parentFaceIndexIn;
    this->childrenFaceIndex[parentFaceIndexIn].push_back(faceIndex);
}

void Mesh_FaceTopo::ScanChildFace(PHVectorInt1D cell2Node, int cellType, int cellIndex)
{
    BasicElement *basicElement = GetBasicElement(cellType);
    uint_t faceNumber = basicElement->GetElementFaceNumber();
    PHVectorInt1D absoluteFace2Node;

    for (int iFace = 0; iFace < faceNumber; ++ iFace)
    {
        vector< int > &relativeFace2Node = basicElement->GetElementFace(iFace);
        int iFaceType = basicElement->GetFaceElementType(iFace);

        uint_t nodeNumber = relativeFace2Node.size();
        absoluteFace2Node.resize(nodeNumber);

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            int relativeIndex = relativeFace2Node[iNode];
            int nodeIndex = cell2Node[relativeIndex];
            absoluteFace2Node[iNode] = nodeIndex;
        }

        PHVectorInt1D sortedFace2Node = absoluteFace2Node;
        sort(sortedFace2Node.begin(), sortedFace2Node.end());

        MultiIndex faceOfMultiIndex(static_cast<int>(nodeNumber), GetNTotalFaces());
        faceOfMultiIndex.SetData(&sortedFace2Node[0]);

        int globalFaceIndex = FindFace(faceOfMultiIndex);
        if (globalFaceIndex == INVALID_INDEX)
        {
            uint_t globalFaceNumber = faceList.size();
            uint_t faceNumber1 = leftCellIndex.size();
            if (globalFaceNumber != faceNumber1)
            {
                cout << "globalFaceNumber != faceNumber " << globalFaceNumber << " " << faceNumber1 << "\n";
                TK_Exit::PrintDebugInfoExit("");
            }

            InsertFace(faceOfMultiIndex);

            leftCellIndex.push_back(cellIndex);
            rightCellIndex.push_back(INVALID_INDEX);
            faceBoundaryCondition.push_back(0);
            faceType.push_back(iFaceType);

            face2Node.push_back(absoluteFace2Node);
            parentFaceIndex.push_back(INVALID_INDEX);
            childrenFaceIndex.resize(childrenFaceIndex.size() + 1);
        }
        else
        {
            if (leftCellIndex[globalFaceIndex] == INVALID_INDEX)
            {
                leftCellIndex[globalFaceIndex] = cellIndex;
            }
            else
            {
                int rightCell = rightCellIndex[globalFaceIndex];
                bool isBoundaryFace = PHSPACE::IsBoundaryFace(rightCell);

                if (isBoundaryFace)
                {
                    if (leftCellIndex[globalFaceIndex] != cellIndex)
                    {
                        rightCellIndex[globalFaceIndex] = cellIndex;
                    }
                }
            }
        }
    }
}

}