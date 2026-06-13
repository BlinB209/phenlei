#include "Mesh_Refine_Unstruct.h"
#include "GlobalDataBase.h"
#include "TK_Log.h"
#include "TK_Time.h"
#include "Glb_Dimension.h"
#include "Geo_Element.h"
#include "GridType.h"
#include "PHHeader.h"
#include "ComputationalGrid.h"

using namespace std;
namespace PHSPACE
{
Mesh_Refine_Unstruct::Mesh_Refine_Unstruct(const string &gridFileNameIn) :
    Mesh_Refine(gridFileNameIn)
{
    gridRefineControllers = 0;
}

Mesh_Refine_Unstruct::~Mesh_Refine_Unstruct()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete gridRefineControllers[iZone];
    }
    delete [] gridRefineControllers;

}

void Mesh_Refine_Unstruct::GenerateAndDumpComputationalGrid()
{
    WriteLogFile("Start Generating Grid!");

    AllocateComputationalGrid();

    GenerateComputationalGrid();

    GenerateGlobalInterfaceLink();

    DumpComputationalGrid();
}

void Mesh_Refine_Unstruct::GenerateGlobalInterfaceLink()
{
    ComputationalGridProcess *computationalGridProcess = new ComputationalGridProcess(region);
    computationalGridProcess->Initialize(refinedGrids, numberOfZones);
    computationalGridProcess->PostProcessMultiZoneComputationalGrids();

    delete computationalGridProcess;
}

void Mesh_Refine_Unstruct::GenerateComputationalGrid()
{
    int nTCellLocal = 0;
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        gridRefineControllers[iZone]->GenerateComputationalGrid(refinedGrids[iZone]);

        nTCellLocal += refinedGrids[iZone]->GetNTotalCell();
    }

    RDouble nTCellGlobal = nTCellLocal * 1.0;
    RDouble ntCellLocalTemp = nTCellLocal * 1.0;
    PH_Reduce(&ntCellLocalTemp, &nTCellGlobal, 1, PH_SUM);

    int maxNTCellLocal, minNTCellLocal;
    PH_Reduce(&nTCellLocal, &maxNTCellLocal, 1, PH_MAX);
    PH_Reduce(&nTCellLocal, &minNTCellLocal, 1, PH_MIN);
    PrintToWindow("  Number of Global cell: ", nTCellGlobal, "\n");
    PrintToWindow("  Min & Max number of cell after adapt: ", minNTCellLocal, maxNTCellLocal, "\n");
}

void Mesh_Refine_Unstruct::AllocateComputationalGrid()
{
    refinedGrids = new Grid *[numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, new GridID(iZone), 0, GetDim());
        refinedGrids[iZone] = grid;
    }
}

void Mesh_Refine_Unstruct::AllocateMemory()
{
    numberOfZones = PHMPI::GetNumberofLocalZones();
    gridRefineControllers = new Mesh_RefineController * [numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Mesh_RefineController *gridRefineController = new Mesh_RefineController();
        gridRefineController->SetZoneIndex(iZone);
        gridRefineControllers[iZone] = gridRefineController;
    }
}

void Mesh_Refine_Unstruct::ConstructGridTopo()
{
    WriteLogFile("Start Constructing Grid Info!");

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int globalZoneIndex = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(globalZoneIndex, 0);
        grid->ComputeMetrics();
        grid->SkewnessSummary();

        ConstructGridTopo(grid);
    }

    WriteLogFile("End Constructing Grid Info!");
}

void Mesh_Refine_Unstruct::ConstructGridTopo(Grid *grid)
{
    ConstructCellTopo(grid);

    ConstructNodeTopo(grid);

    ConstructFaceTopo(grid);

    ConstructMetrics(grid);

    ConstructInterFaceTopo(grid);
}

void Mesh_Refine_Unstruct::ConstructMetrics(Grid *grid)
{
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    RDouble *faceAreaContainer    = unstructuredGrid->GetFaceArea();
    RDouble *faceNormalXContainer = unstructuredGrid->GetFaceNormalX();
    RDouble *faceNormalYContainer = unstructuredGrid->GetFaceNormalY();
    RDouble *faceNormalZContainer = unstructuredGrid->GetFaceNormalZ();

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];
    Mesh_FaceTopo *faceTopo = currentController->GetFaceFactory();

    PHVectorRDouble1D &faceArea    = faceTopo->GetFaceArea();
    PHVectorRDouble1D &faceNormalX = faceTopo->GetFaceNormalX();
    PHVectorRDouble1D &faceNormalY = faceTopo->GetFaceNormalY();
    PHVectorRDouble1D &faceNormalZ = faceTopo->GetFaceNormalZ();

    int nTotalFace = unstructuredGrid->GetNTotalFace();
    faceArea.resize(nTotalFace);
    faceNormalX.resize(nTotalFace);
    faceNormalY.resize(nTotalFace);
    faceNormalZ.resize(nTotalFace);

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        faceArea   [iFace] = faceAreaContainer   [iFace];
        faceNormalX[iFace] = faceNormalXContainer[iFace];
        faceNormalY[iFace] = faceNormalYContainer[iFace];
        faceNormalZ[iFace] = faceNormalZContainer[iFace];
    }
}

void Mesh_Refine_Unstruct::ConstructCellTopo(Grid *grid)
{
    int dimension = GetDim();
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    int nTotalCell = unstructuredGrid->GetNTotalCell();
    int *nodeNumberOfEachCell = unstructuredGrid->GetNodeNumberOfEachCell();
    int *cell2Node  = unstructuredGrid->GetCell2Node();

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();
    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();
    PHVectorInt1D &globalToComputationalCellIndex = cellTopo->GetGlobalToComputationalCellIndex();

    computationalToGlobalCellIndex.resize(nTotalCell);
    globalToComputationalCellIndex.resize(nTotalCell);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        computationalToGlobalCellIndex[iCell] = iCell;
        globalToComputationalCellIndex[iCell] = iCell;
    }

    PHVectorInt1D &cellType = cellTopo->GetCellType();
    cellType.resize(nTotalCell);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nodeNumber = nodeNumberOfEachCell[iCell];
        int iCellType = GetCellType(dimension, nodeNumber);
        cellType[iCell] = iCellType;
    }

    cellTopo->SetAllCellLevelToZero();
    cellTopo->SetAllParentIndexToRootState();

    PHVectorInt2D &cell2NodeArray = cellTopo->GetCell2Node();
    cell2NodeArray.resize(nTotalCell);
    int count = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nodeNumber = nodeNumberOfEachCell[iCell];
        cell2NodeArray[iCell].resize(nodeNumber);
        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            cell2NodeArray[iCell][iNode] = cell2Node[count ++];
        }
    }

    int anisoRefine = this->refineParameter->GetAnisoRefineType();
    if (anisoRefine)
    {
        int *faceNumberOfEachCell = unstructuredGrid->GetFaceNumberOfEachCell();
        int **cell2Face = unstructuredGrid->GetCell2Face();

        PHVectorInt2D &cell2FaceArray = cellTopo->GetCell2Face();
        cell2FaceArray.resize(nTotalCell);

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int faceNumber = faceNumberOfEachCell[iCell];
            cell2FaceArray[iCell].resize(faceNumber);
            for (int iFace = 0; iFace < faceNumber; ++ iFace)
            {
                cell2FaceArray[iCell][iFace] = cell2Face[iCell][iFace];
            }
        }
    }
}

void Mesh_Refine_Unstruct::ConstructNodeTopo(Grid *grid)
{
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    RDouble *x = unstructuredGrid->GetX();
    RDouble *y = unstructuredGrid->GetY();
    RDouble *z = unstructuredGrid->GetZ();

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];
    Geo_PointFactory *nodeTopo = currentController->GetNodeFactory();

    RDouble minDistance, maxDistance;
    unstructuredGrid->GetMinMaxDS(minDistance, maxDistance);
    nodeTopo->SetMinDistance(minDistance);

    int nTotalNode = unstructuredGrid->GetNTotalNode();
    for (int iPoint = 0; iPoint < nTotalNode; ++ iPoint)
    {
        Point3D *point = new Point3D(x[iPoint], y[iPoint], z[iPoint], iPoint);

        bool isNodeExist;
        nodeTopo->AddPoint(point, isNodeExist);
    }

    int oldNumberOfComputationalNodes = nTotalNode;
    nodeTopo->SetOldNumberOfComputationalNodes(oldNumberOfComputationalNodes);
}

void Mesh_Refine_Unstruct::ConstructFaceTopo(Grid *grid)
{
    int dimension = GetDim();
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    int nTotalFace            = unstructuredGrid->GetNTotalFace();
    int nBoundFace            = unstructuredGrid->GetNBoundFace();
    int nTotalCell            = unstructuredGrid->GetNTotalCell();
    int *face2Node            = unstructuredGrid->GetFace2Node();
    int *leftCellOfFace       = unstructuredGrid->GetLeftCellOfFace();
    int *rightCellOfFace      = unstructuredGrid->GetRightCellOfFace();
    int *nodeNumberOfEachFace = unstructuredGrid->GetNodeNumberOfEachFace();

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];
    Mesh_FaceTopo *faceTopo = currentController->GetFaceFactory();

    PHVectorInt1D &computationalToGlobalFaceIndex = faceTopo->GetComputationalToGlobalFaceIndex();
    computationalToGlobalFaceIndex.resize(nTotalFace);

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        computationalToGlobalFaceIndex[iFace] = iFace;
    }

    PHVectorInt1D &leftCellIndex  = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex = faceTopo->GetRightCellIndex();
    leftCellIndex.resize(nTotalFace);
    rightCellIndex.resize(nTotalFace);

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int leftCell  = leftCellOfFace [iFace];
        int rightCell = rightCellOfFace[iFace];
        if (rightCell >= nTotalCell)
        {
            rightCell = - rightCell;
        }

        leftCellIndex [iFace] = leftCell;
        rightCellIndex[iFace] = rightCell;
    }

    int nodeCount = 0;
    for (int iFace = 0; iFace < nTotalFace; ++iFace)
    {
        int nodeNumber = nodeNumberOfEachFace[iFace];
        PHVectorInt1D iFace2Node;
        iFace2Node.resize(nodeNumber);
        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            iFace2Node[iNode] = face2Node[nodeCount ++];
        }

        sort(iFace2Node.begin(), iFace2Node.end());

        MultiIndex faceOfMultiIndex(nodeNumber, iFace);
        faceOfMultiIndex.SetData(&iFace2Node[0]);
        int id = faceTopo->FindFace(faceOfMultiIndex);
        if (id == PHSPACE::INVALID_INDEX)
        {
            faceTopo->InsertFace(faceOfMultiIndex);
        }
    }

    PHVectorInt2D &face2NodeArray = faceTopo->GetFace2Node();
    face2NodeArray.resize(nTotalFace);
    nodeCount = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nodeNumber = nodeNumberOfEachFace[iFace];
        face2NodeArray[iFace].resize(nodeNumber);

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            face2NodeArray[iFace][iNode] = face2Node[nodeCount ++];
        }
    }

    vector<SimpleBC * > &faceBoundaryCondition = faceTopo->GetFaceBoundaryCondition();
    faceBoundaryCondition.resize(nTotalFace);
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        faceBoundaryCondition[iFace] = 0;
    }
    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iFace]);

        SimpleBC *iFaceBoundaryCondtion = new SimpleBC();
        int boundaryType = bcRegion->GetBCType();
        iFaceBoundaryCondtion->SetBCType(boundaryType);
        if (boundaryType != PHENGLEI::INTERFACE)
        {
            string bcName = bcRegion->GetBCName();
            iFaceBoundaryCondtion->SetBCName(bcName);
        }
        else
        {
            string bcName = "Interface";
            iFaceBoundaryCondtion->SetBCName(bcName);
        }
        faceBoundaryCondition[iFace] = iFaceBoundaryCondtion;
    }

    PHVectorInt1D &faceType = faceTopo->GetFaceType();
    faceType.resize(nTotalFace);
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nodeNumber  = nodeNumberOfEachFace[iFace];
        int iFaceType   = GetFaceType(dimension, nodeNumber);
        faceType[iFace] = iFaceType;
    }

    PHVectorInt1D &parentFaceIndex = faceTopo->GetParentFaceIndex();
    parentFaceIndex.resize(nTotalFace);
    SetField(parentFaceIndex, PHSPACE::INVALID_INDEX, nTotalFace);

    PHVectorInt2D &childrenFaceIndex = faceTopo->GetChildrenFaceIndex();
    childrenFaceIndex.resize(nTotalFace);

    PHVectorInt1D &childrenToCompositeFaceIndex = faceTopo->GetChildrenToCompositeFaceIndex();
    childrenToCompositeFaceIndex.resize(nTotalFace);
}

void Mesh_Refine_Unstruct::ConstructInterFaceTopo(Grid *grid)
{
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];

    currentController->SetInterfaceInformation(unstructuredGrid->GetInterfaceInfo());
}

void Mesh_Refine_Unstruct::RefineGrid()
{
    WriteLogFile("Start Grid Adaptive!");
    PrintToWindow("Start Grid Adaptive!\n");

    PH_Barrier();
    TimeSpan timeSpan;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        currentController->RefineGrid();
    }

    PH_Barrier();
    WriteLogFile("Time consuming of grid refinement: ", timeSpan.GetTimeSpan());

    WriteLogFile("End Grid Adaptive!");
    PrintToWindow("End Grid Adaptive!\n");
}

void Mesh_Refine_Unstruct::FixAnisoRefineType()
{
    InitAnisoRefineType();

    if (GetDim() == PHSPACE::TWO_D) return;

    int anisoRefine = this->refineParameter->GetAnisoRefineType();
    if (anisoRefine != ANISOTROPICREFINE) return;

    WriteLogFile("Start to fix cell type!");
    TimeSpan fixTime;

    int infectTime = 0;
    int isContinueInfect = 1;
    while (isContinueInfect)
    {
        isContinueInfect = 0;

        if (infectTime > 0)
        {
            PrintToWindow("\nInfect isotropic refine type time: ", infectTime, "......");
            WriteLogFile("\nInfect isotropic refine type time: ", infectTime, "......");
        }

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Mesh_RefineController *currentController = gridRefineControllers[iZone];
            currentController->FixAnisoRefineType();
        }

        SwapInterfaceRefineType();

        isContinueInfect = IsNeedInfectInterfaceIsotropicRefineType();

        PH_CompareMaxMin(isContinueInfect, 1);

        ++ infectTime;
    }

    PrintToWindow("\n");

    fixTime.ShowTimeSpanToLogFile("Fix cell refine type");
    WriteLogFile("End fix cell type!");
}

bool Mesh_Refine_Unstruct::IsNeedInfectInterfaceIsotropicRefineType()
{
    bool isContinue = false;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];

        if (currentController->IsNeedInfectInterfaceIsotropicRefineType())
        {
            isContinue = true;
            break;
        }
    }

    return isContinue;
}

void Mesh_Refine_Unstruct::SwapInterfaceRefineType()
{
    WriteLogFile("Start Swapping refine type!");

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer * > > receivedDataBuffer;
    vector< vector< DataContainer * > > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Compress the send information into the actkey.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
                Mesh_RefineController *currentController = gridRefineControllers[localZoneIndex];
                currentController->CompressRefineType(sendBuffer, neighborZone);
            }
        }
    }

    //! Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->Size();
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                std::size_t numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (std::size_t iNeighbor1 = 0; iNeighbor1 < numberOfNeighborTemp; ++ iNeighbor1)
                {
                    uint_t neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(iNeighbor1);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = static_cast<int>(iNeighbor1);
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->Size();
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    uint_t lastID = receivedDataBuffer[iZone].size();
                    SWAP(receivedDataBuffer[iZone][lastID - 1], sendDataBuffer[neighborZone][neighborOrer]);
                }
            }
        }
    }

    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Decompress the interface data from datacontainer.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(neighborZone);
                Mesh_RefineController *currentController = gridRefineControllers[localZoneIndex];
                currentController->DecompressRefineType(receiveData, iZone);

                ++ count;
            }
        }
    }

    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();

    WriteLogFile("End Swapping refine type!");
}

void Mesh_Refine_Unstruct::InitAnisoRefineType()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        currentController->InitAnisoRefineType();
    }
}

void Mesh_Refine_Unstruct::BuildRefineProperty()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        SetDefaultRefineProperty(iZone);

        ComputeRefineProperty(iZone);

        ResetCellProperty(iZone);
    }
}

void Mesh_Refine_Unstruct::ResetCellProperty(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    PHVectorInt1D &cellComputationalStatus = cellTopo->GetCellComputationalStatus();
    PHVectorInt1D &cellModifiedStatus = cellTopo->GetCellModifiedStatus();

    SetField(cellModifiedStatus, SPLIT, static_cast<int>(cellModifiedStatus.size()));
    SetField(cellComputationalStatus, HIDDEN, static_cast<int>(cellComputationalStatus.size()));
}

void Mesh_Refine_Unstruct::ComputeRefineProperty(int iZone)
{
    PrepareDataStructure(iZone);

    SetCellModifiedStatusToNoChange(iZone);

    SetCellComputationalStatusToHiddenState(iZone);

    SetCurrentCellModifiedStatusByOldCellModifiedStatus(iZone);

    SetCurrentCellComputationalStatus(iZone);
}

void Mesh_Refine_Unstruct::PrepareDataStructure(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];

    currentController->Initialize();
    currentController->ConstructNode2Cell();
    currentController->ConstructCell2Cell();
}

void Mesh_Refine_Unstruct::SetCellModifiedStatusToNoChange(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    int numberOfTotal = cellTopo->GetNTotalCell() + currentController->GetNBoundaryFace();

    PHVectorInt1D &cellModifiedStatus = cellTopo->GetCellModifiedStatus();
    cellModifiedStatus.resize(numberOfTotal);
    SetField(cellModifiedStatus, NOCHANGE, numberOfTotal);
}

void Mesh_Refine_Unstruct::SetCellComputationalStatusToHiddenState(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    int numberOfTotal = cellTopo->GetNTotalCell() + currentController->GetNBoundaryFace();

    PHVectorInt1D &cellComputationalStatus = cellTopo->GetCellComputationalStatus();
    cellComputationalStatus.resize(numberOfTotal);
    SetField(cellComputationalStatus, HIDDEN, numberOfTotal);
}

void Mesh_Refine_Unstruct::SetCurrentCellModifiedStatusByOldCellModifiedStatus(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    PHVectorInt1D &computationalCellModifiedStatus = cellTopo->GetComputationalCellModifiedStatus();
    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();

    int oldNumberOfCell = currentController->GetOldNumberOfCell();
    for (int iCell = 0; iCell < oldNumberOfCell; ++ iCell)
    {
        int cellIndex = computationalToGlobalCellIndex[iCell];
        int cellModifiedStatus = computationalCellModifiedStatus[iCell];
        cellTopo->SetCellModifiedStatus(cellIndex, cellModifiedStatus);
    }
}

void Mesh_Refine_Unstruct::SetCurrentCellComputationalStatus(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    PHVectorInt1D &computationalCellModifiedStatus = cellTopo->GetComputationalCellModifiedStatus();
    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();

    int oldNumberOfCell = currentController->GetOldNumberOfCell();
    for (int iCell = 0; iCell < oldNumberOfCell; ++ iCell)
    {
        int cellIndex = computationalToGlobalCellIndex[iCell];
        int cellModifiedStatus = computationalCellModifiedStatus[iCell];

        if (cellModifiedStatus == NOCHANGE)
        {
            cellTopo->SetComputationalStatusOfCell(cellIndex, ON);
        }
    }
}

void Mesh_Refine_Unstruct::SetDefaultRefineProperty(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    int nTotalCell = cellTopo->GetNTotalCell();
    PHVectorInt1D &computationalCellModifiedStatus = cellTopo->GetComputationalCellModifiedStatus();
    computationalCellModifiedStatus.resize(nTotalCell);
    SetField(computationalCellModifiedStatus, NOCHANGE, nTotalCell);
}

}


