#include "Mesh_RefineController.h"
#include "Geo_Element.h"
#include "TK_Log.h"

namespace PHSPACE
{

Mesh_RefineController::Mesh_RefineController()
{
    zoneIndex = 0;
    nTotalNode = 0;
    nBoundaryFace = 0;
    oldNumberOfCell = 0;
    hasBeenRefined = false;

    cellTopo = new Mesh_CellTopo();
    nodeTopo = new Geo_PointFactory();
    faceTopo = new Mesh_FaceTopo();
    computationalFaceTopo = new Mesh_FaceTopo();
    interfaceInformation = 0;
}

Mesh_RefineController::~Mesh_RefineController()
{
    delete cellTopo;
    delete nodeTopo;
    delete faceTopo;
    delete computationalFaceTopo;
}

void Mesh_RefineController::Initialize()
{
    this->SetNBoundaryFace(faceTopo->GetNBoundaryFaces());
    this->SetOldNumberOfCell(cellTopo->GetNTotalCell());
}

void Mesh_RefineController::ConstructNode2Cell()
{
    int nTotalNode = static_cast<int>(nodeTopo->GetNumberOfPoints());
    this->SetNTotalNode(nTotalNode);

    PHVector1D < set< int > > &node2Cell = cellTopo->GetNode2Cell();

    PHVectorInt1D &leftCellIndex  = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex = faceTopo->GetRightCellIndex();
    PHVectorInt2D &face2NodeArray = faceTopo->GetFace2Node();

    node2Cell.resize(nTotalNode);

    int nTotalFace = faceTopo->GetNTotalFaces();
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int leftCell  = leftCellIndex [iFace];
        int rightCell = rightCellIndex[iFace];

        uint_t nodeNumber = face2NodeArray[iFace].size();
        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            int nodeIndex = face2NodeArray[iFace][iNode];

            node2Cell[nodeIndex].insert(leftCell);
            if (!PHSPACE::IsBoundaryFace(rightCell))
            {
                node2Cell[nodeIndex].insert(rightCell);
            }
        }
    }
}

void Mesh_RefineController::ConstructCell2Cell()
{
    PHVector1D < set< int > > &cell2Cell = cellTopo->GetCell2Cell();
    PHVector1D < set< int > > &node2Cell = cellTopo->GetNode2Cell();
    PHVectorInt2D &cell2Node = cellTopo->GetCell2Node();

    int nTotalCell = cellTopo->GetNTotalCell();
    cell2Cell.resize(nTotalCell);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        PHVectorInt1D &iCell2Node = cell2Node[iCell];

        int nodeNumber = static_cast<int>(iCell2Node.size());
        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            int nodeIndex = iCell2Node[iNode];;
            set< int > &iNode2Cell = node2Cell[nodeIndex];

            for (set< int >::iterator iterN = iNode2Cell.begin(); iterN != iNode2Cell.end(); ++ iterN)
            {
                int cellIndex = *iterN;
                if (cellIndex != iCell)
                {
                    cell2Cell[iCell].insert(cellIndex);
                }
            }
        }
    }
}

void Mesh_RefineController::InitAnisoRefineType()
{
    int nTotalCell = cellTopo->GetNTotalCell();
    int nBoundFace = faceTopo->GetNBoundaryFaces();
    int nTotal = nTotalCell + nBoundFace;

    PHVectorInt1D &cellType = cellTopo->GetCellType();
    PHVectorInt1D &refineType = cellTopo->GetRefineType();
    refineType.resize(nTotal);

    int anisoRefine = GlobalDataBase::GetIntParaFromDB("anisoRefine");
    if (anisoRefine == ISOTROPICREFINE)
    {
        //! Set all cells to be isotropic refine.
        SetField(refineType, ISOTROPICREFINE, nTotal);
    }
    else
    {
        //! Set all cells to be anisotropic refine.
        SetField(refineType, ANISOTROPICREFINE, nTotal);

        //! Set all the non-prism cells to be isotropic refine.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            if (cellType[iCell] != PENTA_6)
            {
                refineType[iCell] = ISOTROPICREFINE;
            }
        }
    }
}

void Mesh_RefineController::FixAnisoRefineType()
{
    int nTotalCell = cellTopo->GetNTotalCell();
    int nBoundFace = faceTopo->GetNBoundaryFaces();
    int nTotalFace = faceTopo->GetNTotalFaces();

    PHVectorInt2D &cell2Face = cellTopo->GetCell2Face();
    PHVectorInt1D &refineType = cellTopo->GetRefineType();

    PHVectorInt1D &faceType = faceTopo->GetFaceType();
    PHVectorInt1D &leftCellIndex = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex = faceTopo->GetRightCellIndex();
    vector<SimpleBC *> &faceBoundaryCondition = faceTopo->GetFaceBoundaryCondition();

    deque<int> cellQueue;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (faceType[iFace] != QUAD_4) continue;
        if (faceBoundaryCondition[iFace]->GetBCType() != PHENGLEI::INTERFACE) continue;

        int leftCell  = leftCellIndex[iFace];
        int rightCell = rightCellIndex[iFace];

        if (leftCell >= nTotalCell)
        {
            SWAP(leftCell, rightCell);
        }
        if (rightCell < 0)
        {
            rightCell = -rightCell;
        }

        if (refineType[leftCell] != ISOTROPICREFINE && refineType[rightCell] == ISOTROPICREFINE)
        {
            cellQueue.push_back(leftCell);
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if (faceType[iFace] != QUAD_4) continue;

        int leftCell  = leftCellIndex[iFace];
        int rightCell = rightCellIndex[iFace];

        if (refineType[leftCell] == ISOTROPICREFINE && refineType[rightCell] != ISOTROPICREFINE)
        {
            cellQueue.push_back(rightCell);
        }
        if (refineType[leftCell] != ISOTROPICREFINE && refineType[rightCell] == ISOTROPICREFINE)
        {
            cellQueue.push_back(leftCell);
        }
    }

    while (!cellQueue.empty())
    {
        int currentCellID = cellQueue.front();
        cellQueue.pop_front();

        if (refineType[currentCellID] == ISOTROPICREFINE)
        {
            continue;
        }

        refineType[currentCellID] = ISOTROPICREFINE;

        int nFace = static_cast<int>(cell2Face[currentCellID].size());
        for (int iFace = 0; iFace < nFace; ++ iFace)
        {
            int faceID = cell2Face[currentCellID][iFace];

            if (faceID < nBoundFace) continue;
            if (faceType[faceID] != QUAD_4) continue;

            int leftCell  = leftCellIndex[faceID];
            int rightCell = rightCellIndex[faceID];

            int nextCell = rightCell;
            if (rightCell == currentCellID)
            {
                nextCell = leftCell;
            }

            cellQueue.push_back(nextCell);
        }
    }
}

void Mesh_RefineController::RefineGrid()
{
    WriteLogFile("Start Splitting Grid!");

    int numberOfGridSelfRefine = 1;
    for (int iGridSelfRefine = 0; iGridSelfRefine < numberOfGridSelfRefine; ++ iGridSelfRefine)
    {
        this->SplitAllElements();
    }

    this->SetRefineStatus(true);

    WriteLogFile("End Splitting Grid!");

    ConstructBasicData();
}

void Mesh_RefineController::ConstructBasicData()
{
    WriteLogFile("Start Constructing Basic Data!");

    ScanAllElementsToReconstructFaceDataStructure();

    GenerateCellIndex();

    GenerateGlobalToComputationalCellFaceIndex();

    GenerateComputationalCellFaceInformation();

    GenerateComputationalFaceNodeIndex();
}

void Mesh_RefineController::ScanAllElementsToReconstructFaceDataStructure()
{
    int nTotalCell = cellTopo->GetNTotalCell();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int iCellType = cellTopo->GetCellType(iCell);
        if (!IsBasicVolumeElementType(iCellType))
        {
            continue;
        }

        faceTopo->ScanChildFace(cellTopo->GetCell2Node(iCell), iCellType, iCell);
    }

    PrintToWindow("  Total Element Number = ", nTotalCell, "\n");
}

void Mesh_RefineController::GenerateCellIndex()
{
    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();
    PHVectorInt1D &globalToComputationalCellIndex = cellTopo->GetGlobalToComputationalCellIndex();

    int nTotalCell = cellTopo->GetNTotalCell();
    globalToComputationalCellIndex.resize(nTotalCell);

    int computationalCell = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int computationalStatusOfCell = cellTopo->GetCellComputationalStatus(iCell);
        if (computationalStatusOfCell == ON)
        {
            ++computationalCell;
        }
    }

    int numberOfComputationalCell = computationalCell;
    computationalToGlobalCellIndex.resize(numberOfComputationalCell);

    computationalCell = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int computationalStatusOfCell = cellTopo->GetCellComputationalStatus(iCell);
        if (computationalStatusOfCell == ON)
        {
            computationalToGlobalCellIndex[computationalCell] = iCell;
            globalToComputationalCellIndex[iCell] = computationalCell;

            ++ computationalCell;
        }
    }
}

void Mesh_RefineController::GenerateComputationalFaceNodeIndex()
{
    PHVectorInt1D &computationalToGlobalNodeIndex = nodeTopo->GetComputationalToGlobalNodeIndexMapping();
    PHVectorInt1D &globalToComputationalNodeIndex = nodeTopo->GetGlobalToComputationalNodeIndexMapping();

    PHVectorInt2D &intermediateFace2Node  = computationalFaceTopo->GetIntermediateFace2Node();
    PHVectorInt2D &computationalFace2Node = computationalFaceTopo->GetFace2Node();

    int numberOfFacesForComputation = computationalFaceTopo->GetNTotalFaces();

    set< int > computationalNodeIndexSet;
    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        PHVectorInt1D faceNodesIndexes = intermediateFace2Node[iFaceForComputation];
        uint_t faceNodesIndexesSize = faceNodesIndexes.size();
        for (int iNode = 0; iNode < faceNodesIndexesSize; ++ iNode)
        {
            int globalNodeIndex = faceNodesIndexes[iNode];
            computationalNodeIndexSet.insert(globalNodeIndex);
        }
    }

    int nTotalFaces = faceTopo->GetNTotalFaces();
    PHVectorInt2D &face2Node = faceTopo->GetFace2Node();

    set< int > globalNodeIndexSet;
    for (int iFace = 0; iFace < nTotalFaces; ++ iFace)
    {
        PHVectorInt1D faceNodesIndexes = face2Node[iFace];
        uint_t faceNodesIndexesSize = faceNodesIndexes.size();
        for (int iNode = 0; iNode < faceNodesIndexesSize; ++ iNode)
        {
            int globalNodeIndex = faceNodesIndexes[iNode];
            globalNodeIndexSet.insert(globalNodeIndex);
        }
    }

    uint_t numberOfNodesForComputation = computationalNodeIndexSet.size();
    uint_t numberOfNodes = globalNodeIndexSet.size();

    computationalToGlobalNodeIndex.resize(numberOfNodesForComputation);
    globalToComputationalNodeIndex.resize(numberOfNodes);
    SetField(globalToComputationalNodeIndex, INVALID_INDEX, static_cast<int>(numberOfNodes));

    int iNodeForComputation = 0;
    for (set< int >::iterator iter = computationalNodeIndexSet.begin(); iter != computationalNodeIndexSet.end(); ++ iter)
    {
        int nodeIndex = *iter;

        computationalToGlobalNodeIndex[iNodeForComputation] = nodeIndex;
        globalToComputationalNodeIndex[nodeIndex] = iNodeForComputation;
        ++ iNodeForComputation;
    }

    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        PHVectorInt1D &faceNodeIndexArray = intermediateFace2Node[iFaceForComputation];
        PHVectorInt1D faceNodeIndexArrayForComputation;
        uint_t numberOfNodesOfFace = faceNodeIndexArray.size();
        for (int iNode = 0; iNode < numberOfNodesOfFace; ++ iNode)
        {
            int nodeIndex = faceNodeIndexArray[iNode];
            int nodeIndexForComputation = globalToComputationalNodeIndex[nodeIndex];
            faceNodeIndexArrayForComputation.push_back(nodeIndexForComputation);
        }
        computationalFace2Node.push_back(faceNodeIndexArrayForComputation);
    }
}

void Mesh_RefineController::GenerateComputationalCellFaceInformation()
{
    GenerateComputationalElementFaceInformationFirstPart();

    GenerateComputationalElementFaceInformationSecondPart();
}

void Mesh_RefineController::GenerateComputationalElementFaceInformationSecondPart()
{
    PHVectorInt1D &parentFaceIndex              = faceTopo->GetParentFaceIndex();
    PHVectorInt1D &childrenToCompositeFaceIndex = faceTopo->GetChildrenToCompositeFaceIndex();

    PHVectorInt1D &faceTypeForComputation         = computationalFaceTopo->GetFaceType();
    PHVectorInt1D &computationalToGlobalFaceIndex = computationalFaceTopo->GetComputationalToGlobalFaceIndex();

    int numberOfFacesForComputation = computationalFaceTopo->GetNTotalFaces();
    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        if (faceTypeForComputation[iFaceForComputation] == HANGINGFACE)
        {
            int faceIndex = computationalToGlobalFaceIndex[iFaceForComputation];
            int faceParentID = parentFaceIndex[faceIndex];

            if (faceParentID == INVALID_INDEX)
            {
                cout << "error Generate EdgesWithMiddlePoints \n";
            }

            int iCompositeFace = childrenToCompositeFaceIndex[faceParentID];
            GenerateEdgesWithMiddlePoints(iCompositeFace);
        }
    }

    PHVectorInt2D &intermediateFace2Node = computationalFaceTopo->GetIntermediateFace2Node();

    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        PHVectorInt1D &middleNodeOfEdge = faceTopo->GetMiddleNodeOfEdge();

        PHVectorInt1D &iFace2Node = intermediateFace2Node[iFaceForComputation];
        PHVectorInt1D iFace2NodeForComputation = iFace2Node;

        uint_t nodeNumber = iFace2Node.size();
        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            int iNode0 = iNode;
            int iNode1 = (iNode + 1) % nodeNumber;

            int nodeIndex0 = iFace2Node[iNode0];
            int nodeIndex1 = iFace2Node[iNode1];

            PHVectorInt1D edgeNodeIndexes;
            edgeNodeIndexes.push_back(nodeIndex0);
            edgeNodeIndexes.push_back(nodeIndex1);

            PHVectorInt1D sortedEdgeNodeIndexes = edgeNodeIndexes;
            sort(sortedEdgeNodeIndexes.begin(), sortedEdgeNodeIndexes.end());

            int edgeWithMiddlePointCandidateIndex = faceTopo->GetNumberOfTotalEdgesWithMiddleNode();
            int edgeNodeNumber = 2;

            MultiIndex edgeOfMultiIndex(edgeNodeNumber, edgeWithMiddlePointCandidateIndex);

            edgeOfMultiIndex.SetData(&sortedEdgeNodeIndexes[0]);

            int edgeWithMiddlePointIndex = faceTopo->FindEdgeWithMiddleNode(edgeOfMultiIndex);

            if (edgeWithMiddlePointIndex == INVALID_INDEX) continue;

            int middlePointIndex = middleNodeOfEdge[edgeWithMiddlePointIndex];

            PHVectorInt1D::iterator iter = find(iFace2NodeForComputation.begin(), iFace2NodeForComputation.end(), nodeIndex1);

            if (iter == iFace2NodeForComputation.begin())
            {
                iFace2NodeForComputation.push_back(middlePointIndex);
            }
            else
            {
                iFace2NodeForComputation.insert(iter, middlePointIndex);
            }
        }
        intermediateFace2Node[iFaceForComputation] = iFace2NodeForComputation;
    }
}

void Mesh_RefineController::GenerateEdgesWithMiddlePoints(int iCompositeFace)
{
    PHVectorInt1D &middleNodeOfEdge   = faceTopo->GetMiddleNodeOfEdge();
    PHVectorInt1D &compositeFaceType  = faceTopo->GetCompositeFaceType();
    PHVectorInt2D &compositeFace2Node = faceTopo->GetCompositeFace2Node();

    int iCompositeFaceType = compositeFaceType[iCompositeFace];
    int simpleFaceType = GetSimpleElementType(iCompositeFaceType);
    int simpleFaceNodeNumber = GetElementNodeNumbers(simpleFaceType);

    int edgeNodeNumber = 2;
    for (int iNode = 0; iNode < simpleFaceNodeNumber; ++ iNode)
    {
        int iMiddleNode = simpleFaceNodeNumber + iNode;

        int iNode0 = iNode;
        int iNode1 = (iNode + 1) % simpleFaceNodeNumber;

        PHVectorInt1D edgeNodeIndexes;
        edgeNodeIndexes.push_back(compositeFace2Node[iCompositeFace][iNode0]);
        edgeNodeIndexes.push_back(compositeFace2Node[iCompositeFace][iNode1]);

        PHVectorInt1D sortedEdgeNodeIndexes = edgeNodeIndexes;
        sort(sortedEdgeNodeIndexes.begin(), sortedEdgeNodeIndexes.end());

        int edgeWithMiddleNodeIndexCandidate = faceTopo->GetNumberOfTotalEdgesWithMiddleNode();

        MultiIndex edgeOfMultiIndex(edgeNodeNumber, edgeWithMiddleNodeIndexCandidate);
        edgeOfMultiIndex.SetData(&sortedEdgeNodeIndexes[0]);

        int edgeWithMiddleNodeIndex = faceTopo->FindEdgeWithMiddleNode(edgeOfMultiIndex);
        if (edgeWithMiddleNodeIndex == INVALID_INDEX)
        {
            faceTopo->InsertEdgeWithMiddleNode(edgeOfMultiIndex);
            middleNodeOfEdge.push_back(compositeFace2Node[iCompositeFace][iMiddleNode]);
        }
    }
}

void Mesh_RefineController::GenerateComputationalElementFaceInformationFirstPart()
{
    PHVectorInt1D &computationalToGlobalFaceIndex = computationalFaceTopo->GetComputationalToGlobalFaceIndex();
    PHVectorInt1D &faceTypeForComputation         = computationalFaceTopo->GetFaceType();

    PHVectorInt2D &face2Node = faceTopo->GetFace2Node();
    PHVectorInt2D &intermediateFace2NodeForComputation = computationalFaceTopo->GetIntermediateFace2Node();

    int numberOfFacesForComputation = computationalFaceTopo->GetNTotalFaces();
    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        int faceIndex = computationalToGlobalFaceIndex[iFaceForComputation];
        int faceType = faceTypeForComputation[iFaceForComputation];

        if (faceType == BOUNDARYFACE)
        {
            ProcessBoundaryFace(iFaceForComputation);
        }
        else if (faceType == GENERALFACE)
        {
            ProcessGereralFace(iFaceForComputation);
        }
        else if (faceType == HANGINGFACE)
        {
            ProcessHangingFace(iFaceForComputation);
        }

        PHVectorInt1D &iFace2Node = face2Node[faceIndex];
        intermediateFace2NodeForComputation.push_back(iFace2Node);
    }
}

void Mesh_RefineController::ProcessBoundaryFace(int iFaceForComputation)
{
    PHVectorInt1D &cellComputationalStatus        = cellTopo->GetCellComputationalStatus();
    PHVectorInt1D &globalToComputationalCellIndex = cellTopo->GetGlobalToComputationalCellIndex();

    PHVectorInt1D &leftCellIndex  = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex = faceTopo->GetRightCellIndex();

    PHVectorInt1D &leftCellIndexForComputation    = computationalFaceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndexForComputation   = computationalFaceTopo->GetRightCellIndex();
    PHVectorInt1D &computationalToGlobalFaceIndex = computationalFaceTopo->GetComputationalToGlobalFaceIndex();

    int faceIndex = computationalToGlobalFaceIndex[iFaceForComputation];
    int leftCell = leftCellIndex[faceIndex];
    int rightCell = rightCellIndex[faceIndex];

    bool isRightElementExist = !IsElementNotExist(rightCell);
    if (cellComputationalStatus[leftCell] == HIDDEN)
    {
        string err = " Error: ProcessBoundaryFace error computationalStatus of left Element ";
        TK_Exit::ExceptionExit(err);
    }
    else if (isRightElementExist)
    {
        string err = " Error: ProcessBoundaryFace error computationalStatus of right Element ";
        TK_Exit::ExceptionExit(err);
    }
    else
    {
        int leftCellForComputation = globalToComputationalCellIndex[leftCell];
        leftCellIndexForComputation.push_back(leftCellForComputation);
        rightCellIndexForComputation.push_back(INVALID_INDEX);
    }
}

void Mesh_RefineController::ProcessGereralFace(int iFaceForComputation)
{
    PHVectorInt1D &globalToComputationalCellIndex = cellTopo->GetGlobalToComputationalCellIndex();

    PHVectorInt1D &leftCellIndex  = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex = faceTopo->GetRightCellIndex();

    PHVectorInt1D &leftCellIndexForComputation    = computationalFaceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndexForComputation   = computationalFaceTopo->GetRightCellIndex();
    PHVectorInt1D &computationalToGlobalFaceIndex = computationalFaceTopo->GetComputationalToGlobalFaceIndex();

    int faceIndex = computationalToGlobalFaceIndex[iFaceForComputation];
    int leftCell  = leftCellIndex[faceIndex];
    int rightCell = rightCellIndex[faceIndex];

    int leftCellForComputation  = globalToComputationalCellIndex[leftCell];
    int rightCellForComputation = globalToComputationalCellIndex[rightCell];

    leftCellIndexForComputation.push_back(leftCellForComputation);
    rightCellIndexForComputation.push_back(rightCellForComputation);
}

void Mesh_RefineController::ProcessHangingFace(int iFaceForComputation)
{
    PHVectorInt1D &cellComputationalStatus        = cellTopo->GetCellComputationalStatus();
    PHVectorInt1D &globalToComputationalCellIndex = cellTopo->GetGlobalToComputationalCellIndex();

    PHVectorInt1D &leftCellIndex   = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex  = faceTopo->GetRightCellIndex();
    PHVectorInt1D &parentFaceIndex = faceTopo->GetParentFaceIndex();

    PHVectorInt1D &leftCellIndexForComputation    = computationalFaceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndexForComputation   = computationalFaceTopo->GetRightCellIndex();
    PHVectorInt1D &computationalToGlobalFaceIndex = computationalFaceTopo->GetComputationalToGlobalFaceIndex();

    int faceIndex = computationalToGlobalFaceIndex[iFaceForComputation];
    int leftCell = leftCellIndex[faceIndex];
    int rightCell = rightCellIndex[faceIndex];

    int parentFace = parentFaceIndex[faceIndex];
    if (parentFace == INVALID_INDEX)
    {
        string err = " Error: ProcessHangingFace error parentFaceIndex ";
        TK_Exit::ExceptionExit(err);
    }

    if (cellComputationalStatus[leftCell] == HIDDEN)
    {
        if (rightCell == INVALID_INDEX)
        {
            string err = " Error: ProcessHangingFace error rightElementIndex ";
            TK_Exit::ExceptionExit(err);
        }
        else
        {
            if (cellComputationalStatus[rightCell] == HIDDEN)
            {
                string err = " Error: ProcessHangingFace error rightElementStatus ";
                TK_Exit::ExceptionExit(err);
            }

            int leftCellParent  = leftCellIndex[parentFace];
            int rightCellParent = rightCellIndex[parentFace];

            if (cellComputationalStatus[leftCellParent] == ON)
            {
                int leftCellForComputation = globalToComputationalCellIndex[leftCellParent];
                leftCellIndexForComputation.push_back(leftCellForComputation);
            }
            else if (cellComputationalStatus[rightCellParent] == ON)
            {
                int leftCellForComputation = globalToComputationalCellIndex[rightCellParent];
                leftCellIndexForComputation.push_back(leftCellForComputation);
            }
            else
            {
                string err = " Error: ProcessHangingFace error leftElementIndexForComputation ";
                TK_Exit::ExceptionExit(err);
            }

            int rightCellForComputation = globalToComputationalCellIndex[rightCell];
            rightCellIndexForComputation.push_back(rightCellForComputation);
        }
    }
    else if (cellComputationalStatus[leftCell] == ON)
    {
        int leftCellForComputation = globalToComputationalCellIndex[leftCell];
        leftCellIndexForComputation.push_back(leftCellForComputation);

        if (rightCell == INVALID_INDEX)
        {
            int leftCellParent  = leftCellIndex[parentFace];
            int rightCellParent = rightCellIndex[parentFace];

            if (cellComputationalStatus[leftCellParent] == ON)
            {
                int rightCellForComputation = globalToComputationalCellIndex[leftCellParent];
                rightCellIndexForComputation.push_back(rightCellForComputation);
            }
            else if (cellComputationalStatus[rightCellParent] == ON)
            {
                int rightCellForComputation = globalToComputationalCellIndex[rightCellParent];
                rightCellIndexForComputation.push_back(rightCellForComputation);
            }
            else
            {
                string err = " Error: ProcessHangingFace error rightElementIndexForComputation 0 ";
                TK_Exit::ExceptionExit(err);
            }
        }
        else
        {
            if (cellComputationalStatus[rightCell] == HIDDEN)
            {
                int leftCellParent  = leftCellIndex[parentFace];
                int rightCellParent = rightCellIndex[parentFace];

                if (cellComputationalStatus[leftCellParent] == ON)
                {
                    int rightCellForComputation = globalToComputationalCellIndex[leftCellParent];
                    rightCellIndexForComputation.push_back(rightCellForComputation);
                }
                else if (cellComputationalStatus[rightCellParent] == ON)
                {
                    int rightCellForComputation = globalToComputationalCellIndex[rightCellParent];
                    rightCellIndexForComputation.push_back(rightCellForComputation);
                }
                else
                {
                    string err = " ProcessHangingFace error rightElementIndexForComputation 1 ";
                    TK_Exit::ExceptionExit(err);
                }
            }
            else
            {
                string err = " ProcessHangingFace error rightElementIndexForComputation 2 ";
                TK_Exit::ExceptionExit(err);
            }
        }
    }
    else
    {
        string err = " ProcessHangingFace error ";
        TK_Exit::ExceptionExit(err);
    }
}

void Mesh_RefineController::GenerateGlobalToComputationalCellFaceIndex()
{
    PHVectorInt1D &globalToComputationalFaceIndex = faceTopo->GetGlobalToComputationalFaceIndex();

    int numberOfTotalFaces = faceTopo->GetNTotalFaces();
    globalToComputationalFaceIndex.resize(numberOfTotalFaces);

    int iComputationalFaceCount = 0;
    int iBoundaryFace = 0;

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        GenerateGlobalToComputationalCellFaceIndex(iFace, iComputationalFaceCount, iBoundaryFace);
    }
}

void Mesh_RefineController::GenerateGlobalToComputationalCellFaceIndex(int iFace, int &iComputationalFaceCount, int &iBoundaryFace)
{
    PHVectorInt1D &parentFaceIndex = faceTopo->GetParentFaceIndex();
    int parentFace = parentFaceIndex[iFace];

    if (parentFace == INVALID_INDEX)
    {
        GenerateGlobalToComputationalCellFaceIndexNoParentCase(iFace, iComputationalFaceCount, iBoundaryFace);
    }
    else
    {
        GenerateGlobalToComputationalCellFaceIndexWithParentCase(iFace, iComputationalFaceCount, iBoundaryFace);
    }
}

void Mesh_RefineController::GenerateGlobalToComputationalCellFaceIndexWithParentCase(int iFace, int &iComputationalFaceCount, int &iBoundaryFace)
{
    PHVectorInt1D &cellComputationalStatus = cellTopo->GetCellComputationalStatus();

    PHVectorInt1D &parentFaceIndex = faceTopo->GetParentFaceIndex();
    PHVectorInt1D &leftCellIndex   = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex  = faceTopo->GetRightCellIndex();

    PHVectorInt1D      &faceTypeForComputation              = computationalFaceTopo->GetFaceType();
    vector<SimpleBC *> &faceBoundaryCondition               = faceTopo->GetFaceBoundaryCondition();
    vector<SimpleBC *> &faceBoundaryConditionForComputation = computationalFaceTopo->GetFaceBoundaryCondition();

    PHVectorInt1D &computationalToGlobalFaceIndex = computationalFaceTopo->GetComputationalToGlobalFaceIndex();
    PHVectorInt1D &globalToComputationalFaceIndex = faceTopo->GetGlobalToComputationalFaceIndex();

    int leftCell   = leftCellIndex[iFace];
    int rightCell  = rightCellIndex[iFace];
    int parentFace = parentFaceIndex[iFace];

    bool isLeftCellNotExist = IsElementNotExist(leftCell);
    if (isLeftCellNotExist)
    {
        TK_Exit::PrintDebugInfoExit(" GenerateGlobalToComputationalElementFaceIndexMappingWithParentCase error\n");
    }

    bool isRightCellNotExist = IsElementNotExist(rightCell);
    if (isRightCellNotExist)
    {
        int rightCellIndexOfParent = rightCellIndex[parentFace];

        bool isBoundaryFaceOfParent = PHSPACE::IsBoundaryFace(rightCellIndexOfParent);
        if (isBoundaryFaceOfParent)
        {
            if (cellComputationalStatus[leftCell] == HIDDEN)
            {
                return;
            }
            else
            {
                computationalToGlobalFaceIndex.push_back(iFace);
                globalToComputationalFaceIndex[iFace] = iComputationalFaceCount;

                faceBoundaryConditionForComputation.push_back(faceBoundaryCondition[iFace]);
                faceTypeForComputation.push_back(BOUNDARYFACE);

                ++ iComputationalFaceCount;
                ++ iBoundaryFace;
                return;
            }
        }
        else
        {
            if (cellComputationalStatus[leftCell] == HIDDEN)
            {
                return;
            }
            else if (cellComputationalStatus[leftCell] == ON)
            {
                computationalToGlobalFaceIndex.push_back(iFace);
                globalToComputationalFaceIndex[iFace] = iComputationalFaceCount;

                faceTypeForComputation.push_back(PHSPACE::HANGINGFACE);
                ++ iComputationalFaceCount;

                return;
            }
        }
    }
    else
    {
        if ((cellComputationalStatus[leftCell]  == HIDDEN) &&
            (cellComputationalStatus[rightCell] == ON))
        {
            int leftParent  = leftCellIndex[parentFace];
            int rightParent = rightCellIndex[parentFace];

            if (cellComputationalStatus[leftParent]  == HIDDEN &&
                cellComputationalStatus[rightParent] == HIDDEN)
            {
                return;
            }
            else
            {
                computationalToGlobalFaceIndex.push_back(iFace);
                globalToComputationalFaceIndex[iFace] = iComputationalFaceCount;

                faceTypeForComputation.push_back(HANGINGFACE);

                ++ iComputationalFaceCount;
                return;
            }
        }
        else if ((cellComputationalStatus[leftCell]  == ON) &&
                 (cellComputationalStatus[rightCell] == HIDDEN))
        {
            int leftParent  = leftCellIndex[parentFace];
            int rightParent = rightCellIndex[parentFace];
            if (cellComputationalStatus[leftParent]  == HIDDEN &&
                cellComputationalStatus[rightParent] == HIDDEN)
            {
                return;
            }
            else
            {
                computationalToGlobalFaceIndex.push_back(iFace);
                globalToComputationalFaceIndex[iFace] = iComputationalFaceCount;

                faceTypeForComputation.push_back(HANGINGFACE);

                ++ iComputationalFaceCount;

                return;
            }
        }
        else if ((cellComputationalStatus[leftCell]  == HIDDEN) &&
                 (cellComputationalStatus[rightCell] == HIDDEN))
        {
            return;
        }
        else if ((cellComputationalStatus[leftCell]  == ON) &&
                 (cellComputationalStatus[rightCell] == ON))
        {
            computationalToGlobalFaceIndex.push_back(iFace);
            globalToComputationalFaceIndex[iFace] = iComputationalFaceCount;

            faceTypeForComputation.push_back(GENERALFACE);

            ++ iComputationalFaceCount;

            return;
        }
    }
}

void Mesh_RefineController::GenerateGlobalToComputationalCellFaceIndexNoParentCase(int iFace, int &iComputationalFaceCount, int &iBoundaryFace)
{
    PHVectorInt1D &cellComputationalStatus = cellTopo->GetCellComputationalStatus();
                  
    PHVectorInt1D &leftCellIndex  = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex = faceTopo->GetRightCellIndex();

    PHVectorInt1D      &faceTypeForComputation              = computationalFaceTopo->GetFaceType();
    vector<SimpleBC *> &faceBoundaryCondition               = faceTopo->GetFaceBoundaryCondition();
    vector<SimpleBC *> &faceBoundaryConditionForComputation = computationalFaceTopo->GetFaceBoundaryCondition();

    PHVectorInt1D &computationalToGlobalFaceIndex = computationalFaceTopo->GetComputationalToGlobalFaceIndex();
    PHVectorInt1D &globalToComputationalFaceIndex = faceTopo->GetGlobalToComputationalFaceIndex();

    int leftCell  = leftCellIndex [iFace];
    int rightCell = rightCellIndex[iFace];

    bool isLeftCellNotExist = IsElementNotExist(leftCell);
    if (isLeftCellNotExist)
    {
        TK_Exit::PrintDebugInfoExit("error\n");
    }
    else
    {
        bool isRightCellNotExist = IsElementNotExist(rightCell);
        if (isRightCellNotExist)
        {
            if (cellComputationalStatus[leftCell] == HIDDEN)
            {
                return;
            }
            else
            {
                computationalToGlobalFaceIndex.push_back(iFace);
                globalToComputationalFaceIndex[iFace] = iComputationalFaceCount;

                faceBoundaryConditionForComputation.push_back(faceBoundaryCondition[iFace]);//*/
                faceTypeForComputation.push_back(BOUNDARYFACE);

                ++ iComputationalFaceCount;
                ++ iBoundaryFace;

                return;
            }
        }
        else
        {
            if (cellComputationalStatus[leftCell]  == HIDDEN ||
                cellComputationalStatus[rightCell] == HIDDEN)
            {
                return;
            }
            else
            {
                computationalToGlobalFaceIndex.push_back(iFace);
                globalToComputationalFaceIndex[iFace] = iComputationalFaceCount;

                faceTypeForComputation.push_back(GENERALFACE);
                ++ iComputationalFaceCount;

                return;
            }
        }
    }
}

void Mesh_RefineController::SplitAllElements()
{
    int nTotalCell = cellTopo->GetNTotalCell();
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        SplitElement(iCell);
    }
}

void Mesh_RefineController::SplitElement(int cellIndex)
{
    if (!cellTopo->ElementNeedSplit(cellIndex)) return;

    int cellType = cellTopo->GetCellType(cellIndex);
    if (!IsBasicVolumeElementType(cellType))
    {
        ReactivateChildElement(cellIndex);
        return;
    }

    int refineType                = cellTopo->GetRefineType(cellIndex);
    int typeOfCompositeCell       = GetRefineElemType(cellType, refineType);
    int NodeNumberOfCell          = GetElementNodeNumbers(cellType);
    int NodeNumberOfCompositeCell = GetElementNodeNumbers(typeOfCompositeCell);

    cellTopo->SetCellType(cellIndex, typeOfCompositeCell);

    PHVectorInt2D &cell2Node = cellTopo->GetCell2Node();
    cell2Node[cellIndex].resize(NodeNumberOfCompositeCell);

    for (int iNode = NodeNumberOfCell; iNode < NodeNumberOfCompositeCell; ++ iNode)
    {
        vector< int > &relatedPointListForMiddlePointComputation = GetRelatedPointListForMiddlePointComputation(typeOfCompositeCell, iNode);

        Geo_PointFactory::point_type *centerPoint = nodeTopo->GetCenterPoint(cell2Node[cellIndex], relatedPointListForMiddlePointComputation);
        bool isNodeExist;
        int index = nodeTopo->AddPoint(centerPoint, isNodeExist);
        cell2Node[cellIndex][iNode] = index;
    }

    AddChildElement(cellIndex);
    SplitElementFace(cellIndex);
}

void Mesh_RefineController::ReactivateChildElement(int cellIndex)
{
    PHVectorInt1D &childrenCellIndexe = cellTopo->GetChildrenCellIndexe(cellIndex);
    uint_t childrenNumber = childrenCellIndexe.size();
    for (int iChildCell = 0; iChildCell < childrenNumber; ++ iChildCell)
    {
        int iChildCellIndex = childrenCellIndexe[iChildCell];

        int childCellModifiedStatus = cellTopo->GetCellModifiedStatus(iChildCellIndex);
        if (childCellModifiedStatus != NOCHANGE)
        {
            TK_Exit::PrintDebugInfoExit("ReactivateChildElement error\n");
        }

        cellTopo->SetComputationalStatusOfCell(iChildCellIndex, ON);
    }
}

void Mesh_RefineController::AddChildElement(int parentElementIndex)
{
    int parentCellType  = cellTopo->GetCellType(parentElementIndex);
    int parentCellLevel = cellTopo->GetCellLevel(parentElementIndex);

    BasicElement *basicElement = GetBasicElement(parentCellType);

    //! Judge the dividing method.
    bool usingDividingCase2 = false;
    if (parentCellType == TETRA_10)
    {
        PHVectorInt1D &parentCell2Node = cellTopo->GetCell2Node(parentElementIndex);

        int node7 = parentCell2Node[7 - 1];
        int node9 = parentCell2Node[9 - 1];
        int node5 = parentCell2Node[5 - 1];
        int node10 = parentCell2Node[10 - 1];

        RDouble dx, dy, dz;
        dx = nodeTopo->GetPoint(node7).X() - nodeTopo->GetPoint(node9).X();
        dy = nodeTopo->GetPoint(node7).Y() - nodeTopo->GetPoint(node9).Y();
        dz = nodeTopo->GetPoint(node7).Z() - nodeTopo->GetPoint(node9).Z();
        RDouble dist1 = dx * dx + dy * dy + dz * dz;

        dx = nodeTopo->GetPoint(node5).X() - nodeTopo->GetPoint(node10).X();
        dy = nodeTopo->GetPoint(node5).Y() - nodeTopo->GetPoint(node10).Y();
        dz = nodeTopo->GetPoint(node5).Z() - nodeTopo->GetPoint(node10).Z();
        RDouble dist2 = dx * dx + dy * dy + dz * dz;

        if (dist1 > dist2)
        {
            usingDividingCase2 = true;
        }
    }

    uint_t numberOfChildCell = GetChildElementNumber(parentCellType);

    PHVectorInt1D childCell2Node;
    for (int iChildCell = 0; iChildCell < numberOfChildCell; ++ iChildCell)
    {
        int childCellType = basicElement->GetChildElementType(iChildCell);
        PHVectorInt1D &parentCell2Node = cellTopo->GetCell2Node(parentElementIndex);

        vector< int > childElementRelativeNodeIndex = basicElement->GetChildElementRelativeNodeIndex(iChildCell);
        if (usingDividingCase2)
        {
            childElementRelativeNodeIndex = basicElement->GetChildElementRelativeNodeIndexCase2(iChildCell);
        }

        int nodeNumberOfChildCell = GetElementNodeNumbers(childCellType);
        childCell2Node.resize(nodeNumberOfChildCell);
        for (int iNode = 0; iNode < nodeNumberOfChildCell; ++ iNode)
        {
            childCell2Node[iNode] = parentCell2Node[childElementRelativeNodeIndex[iNode]];
        }

        cellTopo->AddCell(childCell2Node, childCellType, parentCellLevel + 1, parentElementIndex);
    }
}

void Mesh_RefineController::SplitElementFace(int cellIndex)
{
    int cellType = cellTopo->GetCellType(cellIndex);

    BasicElement *basicElement = GetBasicElement(cellType);

    uint_t numberOfFaces = basicElement->GetElementPhysicsFaceNumber();
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int iFaceType = basicElement->GetElementPhysicsFaceType(iFace);

        vector< int > &relativeFace2Node = basicElement->GetElementPhysicsFace(iFace);
        PHVectorInt1D &parentCell2Node = cellTopo->GetCell2Node(cellIndex);

        PHVectorInt1D iFace2Node;
        uint_t nodeNumber = relativeFace2Node.size();
        iFace2Node.resize(nodeNumber);
        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            iFace2Node[iNode] = parentCell2Node[relativeFace2Node[iNode]];
        }

        SplitFace(iFace2Node, iFaceType);
    }
}

void Mesh_RefineController::SplitFace(PHVectorInt1D &iCompositeFace2Node, int iCompositeFaceTypeIn)
{
    BasicElement *basicElement = GetBasicElement(iCompositeFaceTypeIn);

    uint_t numberOfChildFaces  = GetChildElementNumbers(iCompositeFaceTypeIn);
    int chindrenFaceType       = GetSimpleElementType(iCompositeFaceTypeIn);
    int chindrenFaceNodeNumber = GetElementNodeNumbers(chindrenFaceType);

    PHVectorInt1D chindrenFace2Node;
    chindrenFace2Node.resize(chindrenFaceNodeNumber);
    for (int iNode = 0; iNode < chindrenFaceNodeNumber; ++ iNode)
    {
        chindrenFace2Node[iNode] = iCompositeFace2Node[iNode];
    }

    PHVectorInt2D &compositeFace2Node = faceTopo->GetCompositeFace2Node();

    PHVectorInt1D sortedCompositeFace2Node = iCompositeFace2Node;
    sort(sortedCompositeFace2Node.begin(), sortedCompositeFace2Node.end());

    uint_t nodeNumberOfCompositeFace = iCompositeFace2Node.size();

    int candidateCompositeFaceIndex = faceTopo->GetNTotalCompositeFaces();
    MultiIndex compositeFaceOfMultiIndex(static_cast<int>(nodeNumberOfCompositeFace), candidateCompositeFaceIndex);
    compositeFaceOfMultiIndex.SetData(&sortedCompositeFace2Node[0]);

    PHVectorInt1D &compositeFaceType            = faceTopo->GetCompositeFaceType();
    PHVectorInt1D &childrenToCompositeFaceIndex = faceTopo->GetChildrenToCompositeFaceIndex();
    PHVectorInt1D &compositeToChildrenFaceIndex = faceTopo->GetCompositeToChildrenFaceIndex();

    compositeToChildrenFaceIndex.resize(candidateCompositeFaceIndex);

    int compositeFaceIndex = faceTopo->FindCompositeFace(compositeFaceOfMultiIndex);
    int parentFaceIndex = faceTopo->FindFace(chindrenFace2Node);

    if (compositeFaceIndex == INVALID_INDEX)
    {
        faceTopo->InsertCompositeFace(compositeFaceOfMultiIndex);
        compositeFace2Node.push_back(iCompositeFace2Node);

        compositeToChildrenFaceIndex.resize(candidateCompositeFaceIndex + 1);

        int newSize = MAX(parentFaceIndex + 1, static_cast<int>(childrenToCompositeFaceIndex.size()));
        childrenToCompositeFaceIndex.resize(newSize);

        childrenToCompositeFaceIndex[parentFaceIndex] = candidateCompositeFaceIndex;
        compositeToChildrenFaceIndex[candidateCompositeFaceIndex] = parentFaceIndex;

        compositeFaceType.push_back(iCompositeFaceTypeIn);
    }

    for (int iChildFace = 0; iChildFace < numberOfChildFaces; ++ iChildFace)
    {
        int childFaceType = basicElement->GetChildElementType(iChildFace);
        int childFaceNodeNumber = GetElementNodeNumbers(childFaceType);

        vector< int > &childFaceElementRelativeNodeIndex = basicElement->GetChildElementRelativeNodeIndex(iChildFace);

        PHVectorInt1D childFace2Node;
        childFace2Node.resize(childFaceNodeNumber);
        for (int iNode = 0; iNode < childFaceNodeNumber; ++ iNode)
        {
            childFace2Node[iNode] = iCompositeFace2Node[childFaceElementRelativeNodeIndex[iNode]];
        }

        faceTopo->AddFace(childFace2Node, childFaceType, parentFaceIndex);
    }
}

void Mesh_RefineController::GenerateComputationalGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int globalZoneIndex = PHMPI::GetLocalZoneIDToGlobalZoneID(this->zoneIndex);
    grid->SetZoneID(globalZoneIndex);

    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();
    int numberOfComputationalElements = static_cast<int>(computationalToGlobalCellIndex.size());
    int numberOfNodes = static_cast<int>(nodeTopo->GetNumberOfComputationalNodes());

    PrintToWindow("  nTCell of Computational: ", numberOfComputationalElements, "\n");
    grid->SetNTotalCell(numberOfComputationalElements);
    grid->SetNTotalNode(numberOfNodes);

    RDouble *x = new RDouble[numberOfNodes];
    RDouble *y = new RDouble[numberOfNodes];
    RDouble *z = new RDouble[numberOfNodes];

    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    PHVectorInt1D &computationalToGlobalNodeIndexMapping = nodeTopo->GetComputationalToGlobalNodeIndexMapping();
    uint_t numberOfNodesForComputation = computationalToGlobalNodeIndexMapping.size();
    for (int iNode = 0; iNode < numberOfNodesForComputation; ++ iNode)
    {
        int nodeIndex = computationalToGlobalNodeIndexMapping[iNode];

        x[iNode] = nodeTopo->GetPoint(nodeIndex).X();
        y[iNode] = nodeTopo->GetPoint(nodeIndex).Y();
        z[iNode] = nodeTopo->GetPoint(nodeIndex).Z();
    }

    this->GenerateComputationalUnstructuredGrid(grid);
    this->GenerateCellToNode(grid);
}

void Mesh_RefineController::GenerateCellToNode(UnstructGrid *grid)
{
    PHVectorInt2D &cell2Node                      = cellTopo->GetCell2Node();
    PHVectorInt1D &cellComputationalStatus        = cellTopo->GetCellComputationalStatus();
    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();

    int nTotalCell = cellTopo->GetNTotalCell();
    int numberOfComputationalCell = nTotalCell;
    uint_t isComputationalCellExist = computationalToGlobalCellIndex.size();
    if (isComputationalCellExist != 0)
    {
        numberOfComputationalCell = static_cast<int>(computationalToGlobalCellIndex.size());
    }

    vector< vector< int > > cellToNodeArray;
    cellToNodeArray.resize(numberOfComputationalCell);

    int countCell = 0;
    int countNode = 0;
    for (int elementIndex = 0; elementIndex < nTotalCell; ++ elementIndex)
    {
        if (isComputationalCellExist != 0)
        {
            if (cellComputationalStatus[elementIndex] == HIDDEN)
            {
                continue;
            }
        }

        PHVectorInt1D &elementNodeIndexes = cell2Node[elementIndex];
        int numberOfNodes = static_cast <int> (elementNodeIndexes.size());

        ASSERT(countCell < numberOfComputationalCell);

        cellToNodeArray[countCell].resize(numberOfNodes);
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            cellToNodeArray[countCell][iNode] = elementNodeIndexes[iNode];
        }

        countNode += numberOfNodes;
        ++ countCell;
    }

    int *cellNodeNumberContainer = new int[numberOfComputationalCell];
    int *cellNodeIndexContainer  = new int[countNode];

    countNode = 0;
    for (int iCell = 0; iCell < numberOfComputationalCell; ++ iCell)
    {
        int nodeNumber = static_cast<int>(cellToNodeArray[iCell].size());
        cellNodeNumberContainer[iCell] = nodeNumber;

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            cellNodeIndexContainer[countNode ++] = cellToNodeArray[iCell][iNode];
        }
    }

    grid->SetNodeNumberOfEachCell(cellNodeNumberContainer);
    grid->SetCell2Node(cellNodeIndexContainer);
}

void Mesh_RefineController::GenerateComputationalUnstructuredGrid(UnstructGrid *grid)
{
    vector<SimpleBC *> &faceBoundaryCondition = computationalFaceTopo->GetFaceBoundaryCondition();

    int numberOfFaces = computationalFaceTopo->GetNTotalFaces();
    int numberOfBoundaryFaces = static_cast<int>(faceBoundaryCondition.size());
    grid->SetNBoundFace(numberOfBoundaryFaces);
    grid->SetNTotalFace(numberOfFaces);

    UnstructBCSet **bcr = new UnstructBCSet *[numberOfBoundaryFaces];
    grid->SetBCRecord(bcr);

    set<string> bcNameSet;
    uint_long iCount = 0;

    for (int iFace = 0; iFace < numberOfBoundaryFaces; iFace++)
    {  
        SimpleBC *bc = faceBoundaryCondition[iFace];
        if (bc == 0)
        {
            TK_Exit::ExceptionExit("boundary Condition == 0 in generate unstruct grid!");
        }
        if (bc->GetBCName() == "")
        {
            TK_Exit::ExceptionExit("bcName = null in generate unstruct grid!");
        }

        bcr[iCount] = new UnstructBCSet();
        bcr[iCount]->SetKey(bc->GetBCType());
        bcr[iCount]->SetBCName(bc->GetBCName());
        bcr[iCount]->SetBoundaryCondition(bc);

        bcNameSet.insert(bc->GetBCName());

        ++ iCount;
    }
    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());

    grid->CreateUnstructBCSet(nBCRegionUnstruct);

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = new int[numberOfBoundaryFaces];
    set<string>::iterator iter;
    int count = 0;

    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); iter++)
    {
        UnstructBC *unstructBC = new UnstructBC(count);
        unstructBCSet->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < numberOfBoundaryFaces; iFace++)
        {
            SimpleBC *bc = faceBoundaryCondition[iFace];
            if (bc->GetBCName() == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(bc->GetBCType());
                bcRegionIDofBCFace[iFace] = count;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }

        }
        count++;
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);

    ReorderFaceToNodeIndexes(grid);
    ReorderFaceToCellIndexes(grid);
}

void Mesh_RefineController::ReorderFaceToNodeIndexes(UnstructGrid *grid)
{
    int numberOfFaces = grid->GetNTotalFace();

    PHVectorInt1D &rightCellIndex = computationalFaceTopo->GetRightCellIndex();
    PHVectorInt2D &face2Node      = computationalFaceTopo->GetFace2Node();

    int *nodeNumberOfEachFace = new int[numberOfFaces];

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCell = rightCellIndex[iFace];
        if (rightCell == INVALID_INDEX)
        {
            nodeNumberOfEachFace[iCount] = static_cast<int>(face2Node[iFace].size());
            ++ iCount;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCell = rightCellIndex[iFace];
        if (rightCell != INVALID_INDEX)
        {
            nodeNumberOfEachFace[iCount] = static_cast<int>(face2Node[iFace].size());
            ++ iCount;
        }
    }
    grid->SetNodeNumberOfEachFace(nodeNumberOfEachFace);

    int face2NodeSize = 0;
    iCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        face2NodeSize += nodeNumberOfEachFace[iCount];
        ++ iCount;
    }

    int *face2NodeArray = new int[face2NodeSize];
    grid->SetFace2Node(face2NodeArray);

    int iFace2NodeCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCell = rightCellIndex[iFace];
        if (rightCell == INVALID_INDEX)
        {
            PHVectorInt1D &faceNodeIndexes = face2Node[iFace];
            uint_t nodeNumber = faceNodeIndexes.size();
            for (int iNode = 0; iNode < nodeNumber; ++ iNode)
            {
                face2NodeArray[iFace2NodeCount] = faceNodeIndexes[iNode];
                ++ iFace2NodeCount;
            }
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCell = rightCellIndex[iFace];
        if (rightCell != INVALID_INDEX)
        {
            PHVectorInt1D &faceNodeIndexes = face2Node[iFace];
            uint_t nodeNumber = faceNodeIndexes.size();
            for (int iNode = 0; iNode < nodeNumber; ++ iNode)
            {
                face2NodeArray[iFace2NodeCount] = faceNodeIndexes[iNode];
                ++ iFace2NodeCount;
            }
        }
    }
}

void Mesh_RefineController::ReorderFaceToCellIndexes(UnstructGrid *grid)
{
    PHVectorInt1D &leftCellIndex  = computationalFaceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex = computationalFaceTopo->GetRightCellIndex();

    int numberOfFaces = grid->GetNTotalFace();
    int numberOfBoundaryFaces = grid->GetNBoundFace();

    int *leftCellOfFace  = new int[numberOfFaces];
    int *rightCellOfFace = new int[numberOfFaces];

    grid->SetLeftCellOfFace(leftCellOfFace);
    grid->SetRightCellOfFace(rightCellOfFace);

    int iFaceCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCell = rightCellIndex[iFace];
        if (rightCell == INVALID_INDEX)
        {
            leftCellOfFace [iFaceCount] = leftCellIndex [iFace];
            rightCellOfFace[iFaceCount] = rightCellIndex[iFace];
            ++ iFaceCount;
        }
    }

    if (iFaceCount != numberOfBoundaryFaces)
    {
        TK_Exit::PrintDebugInfoExit(" FATAL ERROR : iFaceCount != numberOfBoundaryFaces \n");
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCell = rightCellIndex[iFace];
        if (rightCell != INVALID_INDEX)
        {
            leftCellOfFace [iFaceCount] = leftCellIndex [iFace];
            rightCellOfFace[iFaceCount] = rightCellIndex[iFace];
            ++ iFaceCount;
        }
    }
}

void Mesh_RefineController::CompressRefineType(DataContainer *dataContainer, int neighborZoneIndex)
{
    if (!interfaceInformation) return;

    PHVectorInt1D &refineType = cellTopo->GetRefineType();

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->MoveToBegin();

    int *dataTemp = new int[interfaceNumberBetweenTwoNeighboringZone];
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        this->GetSourceCellIndex(iFace, sourceCell);

        dataTemp[iLocalFace] = refineType[sourceCell];
    }
    PHWrite(dataContainer, dataTemp, interfaceNumberBetweenTwoNeighboringZone);
    delete [] dataTemp;
}

void Mesh_RefineController::DecompressRefineType(DataContainer *dataContainer, int neighborZoneIndex)
{
    if (!interfaceInformation) return;

    PHVectorInt1D &refineType = cellTopo->GetRefineType();

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    dataContainer->MoveToBegin();

    int *dataTemp = new int[interfaceNumberBetweenTwoNeighboringZone];
    PHRead(dataContainer, dataTemp, interfaceNumberBetweenTwoNeighboringZone);

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        int targetCell;
        this->GetTargetIndex(iFace, targetCell);

        refineType[targetCell] = dataTemp[iLocalFace];
    }

    delete [] dataTemp;
}

bool Mesh_RefineController::IsNeedInfectInterfaceIsotropicRefineType()
{
    int nTotalCell            = cellTopo->GetNTotalCell();
    PHVectorInt1D &refineType = cellTopo->GetRefineType();

    int nBoundFace                            = faceTopo->GetNBoundaryFaces();
    PHVectorInt1D &faceType                   = faceTopo->GetFaceType();
    PHVectorInt1D &leftCellIndex              = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex             = faceTopo->GetRightCellIndex();
    vector<SimpleBC *> &faceBoundaryCondition = faceTopo->GetFaceBoundaryCondition();

    bool doseInterfaceIsotropicRefineTypeExist = false;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (faceBoundaryCondition[iFace]->GetBCType() != PHENGLEI::INTERFACE) continue;
        if (faceType[iFace] != QUAD_4) continue;

        int leftCell  = leftCellIndex[iFace];
        int rightCell = rightCellIndex[iFace];

        if (leftCell >= nTotalCell)
        {
            SWAP(leftCell, rightCell);
        }
        if (rightCell < 0)
        {
            rightCell = -rightCell;
        }

        if (refineType[leftCell] != ISOTROPICREFINE && refineType[rightCell] == ISOTROPICREFINE)
        {
            doseInterfaceIsotropicRefineTypeExist = true;
            break;
        }
    }

    return doseInterfaceIsotropicRefineTypeExist;
}

void Mesh_RefineController::GetSourceCellIndex(int iFace, int &sourceIndex)
{
    int *localInterfaceToBoundaryFaceIndexesMapping = this->interfaceInformation->GetInterFace2BoundaryFace();

    PHVectorInt1D &leftCellIndex = faceTopo->GetLeftCellIndex();
    int boundaryFaceIndex = localInterfaceToBoundaryFaceIndexesMapping[iFace];
    sourceIndex = leftCellIndex[boundaryFaceIndex];
}

void Mesh_RefineController::GetTargetIndex(int iFace, int &targetIndex)
{
    int *localInterfaceToBoundaryFaceIndexesMapping = interfaceInformation->GetInterFace2BoundaryFace();

    PHVectorInt1D &rightCellIndex = faceTopo->GetRightCellIndex();
    int boundaryFaceIndex = localInterfaceToBoundaryFaceIndexesMapping[iFace];
    targetIndex = rightCellIndex[boundaryFaceIndex];
    targetIndex = -targetIndex;
}

}