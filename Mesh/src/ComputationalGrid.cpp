#include "ComputationalGrid.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"
#include "InterfaceLinkStructure.h"
using namespace std;

namespace PHSPACE
{
ZoneInterface::ZoneInterface()
{
    pmin = 0;
    pmax = 0;
    x    = 0;
    y    = 0;
    z    = 0;
    zoneIndex = 0;
    numberOfNodes = 0;
    disMin = 0.0;
    disMax = 0.0;
}

ZoneInterface::~ZoneInterface()
{
    delete [] pmin;
    delete [] pmax;
    delete [] x;
    delete [] y;
    delete [] z;
}

void ZoneInterface::Encode(DataContainer *&dataContainer)
{
    dataContainer->MoveToBegin();

    int zoneIndex = this->GetZoneIndex();
    PHWrite(dataContainer, zoneIndex);

    int dimension = 3;
    RDouble *minBox = this->GetMinBox();
    RDouble *maxBox = this->GetMaxBox();
    PHWrite(dataContainer, minBox, dimension);
    PHWrite(dataContainer, maxBox, dimension);
    PHWrite(dataContainer, disMin);
    PHWrite(dataContainer, disMax);

    int numberOfNodes = this->GetNumberOfNodes();
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    PHWrite(dataContainer, numberOfNodes);
    PHWrite(dataContainer, x, numberOfNodes);
    PHWrite(dataContainer, y, numberOfNodes);
    PHWrite(dataContainer, z, numberOfNodes);

    uint_t numberOfFaces = this->GetNumberOfFaces();
    PHVectorInt2D &faceNodeIndexContainer = this->GetFaceNodeIndexContainer();
    PHWrite(dataContainer, numberOfFaces);
    PHWrite(dataContainer, faceNodeIndexContainer, numberOfFaces);
}

void ZoneInterface::Decode(DataContainer *&dataContainer)
{
    dataContainer->MoveToBegin();

    int zoneIndex;
    PHRead(dataContainer, zoneIndex);
    this->SetZoneIndex(zoneIndex);

    int dimension = 3;
    RDouble *minBox = new RDouble [dimension];
    RDouble *maxBox = new RDouble [dimension];
    PHRead(dataContainer, minBox, dimension);
    PHRead(dataContainer, maxBox, dimension);
    PHRead(dataContainer, disMin);
    PHRead(dataContainer, disMax);
    this->SetMinBox(minBox);
    this->SetMaxBox(maxBox);

    int numberOfNodes;
    PHRead(dataContainer, numberOfNodes);

    RDouble *x = new RDouble [numberOfNodes];
    RDouble *y = new RDouble [numberOfNodes];
    RDouble *z = new RDouble [numberOfNodes];
    PHRead(dataContainer, x, numberOfNodes);
    PHRead(dataContainer, y, numberOfNodes);
    PHRead(dataContainer, z, numberOfNodes);
    this->SetNumberOfNodes(numberOfNodes);
    this->SetX(x);
    this->SetY(y);
    this->SetZ(z);

    uint_t numberOfFaces;
    PHRead(dataContainer, numberOfFaces);

    PHVectorInt2D &faceNodeIndexContainer = this->GetFaceNodeIndexContainer();
    faceNodeIndexContainer.resize(numberOfFaces);
    PHRead(dataContainer, faceNodeIndexContainer, numberOfFaces);
}

ComputationalGridProcess::ComputationalGridProcess(Region *region)
{
    this->region = region;

    interfaceLinkStructure = 0;
}

ComputationalGridProcess::~ComputationalGridProcess()
{
    delete interfaceLinkStructure;
}

void ComputationalGridProcess::Initialize(Grid **gridContainer, int nZones)
{
    this->numberOfZones = nZones;
    this->gridContainer = gridContainer;

    boundaryTypeVectorList.resize(nZones);
    boundaryNameVectorList.resize(nZones);
    newLocalInterfaceIndexToGlobalInterfaceIndexMappingList.resize(nZones);
    faceNodeIndexVectorList.resize(nZones);
    newBoundaryLeftCellIndexList.resize(nZones);
    newLeftCellIndexVectorList.resize(nZones);
    newRightCellIndexVectorList.resize(nZones);
    zoneInterfacesTopologyList.resize(nZones);
}

void ComputationalGridProcess::PostProcessMultiZoneComputationalGrids()
{
    WriteLogFile("Start Generating Global Interface Link!");

    CreateNeighborZoneInterfacesTopology();

    CreateInterfaceLinkStructure();

    PreprocessBoundaryCondition();

    GenerateLocalGlobalInterfaceMappingByNeighborZones();

    ReconstructInterfaceTopology();

    ReGenerateLocalGlobalInterfaceMapping();

    MatchInterfaceTopology();

    PHSPACE::WriteLogFile("End Generating Global Interface Link!");
}

void ComputationalGridProcess::CreateInterfaceLinkStructure()
{
    RDouble mindis, maxdis, pmin[3], pmax[3];

    GetNeighborZoneBoundingBox(pmin, pmax);
    GetNeighborZoneMinMaxDistance(mindis, maxdis);

    PrintToWindow(" mindis = ", mindis, "\n");

    if (GetDim() == TWO_D)
    {
        pmin[ 2 ] -= mindis;
        pmax[ 2 ] += maxdis;
    }

    vector< ZoneInterface * > &neighborZoneInterfaces = this->GetNeighborZoneInterfaces();
    int maxNumberOfNeighborZones = 0;
    for (std::size_t iZone = 0; iZone < neighborZoneInterfaces.size(); ++ iZone)
    {
        int zoneNumber = neighborZoneInterfaces[iZone]->GetZoneIndex();
        maxNumberOfNeighborZones = MAX(maxNumberOfNeighborZones, zoneNumber);
    }

    interfaceLinkStructure = new InterfaceLinkStructure(pmin, pmax, mindis / 5, maxNumberOfNeighborZones + 1);
}

void ComputationalGridProcess::PreprocessBoundaryCondition()
{
    ChangeBoundaryConditionTypeIfNeeded();
}

void ComputationalGridProcess::ChangeBoundaryConditionTypeIfNeeded()
{
    int ignoreNoBoundaryCondition = GlobalDataBase::GetIntParaFromDB("omit_no_bound_bc");

    if (ignoreNoBoundaryCondition == 1)
    {
        return;
    }

    Grid **gridContainer = this->GetGridContainer();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Grid *grid = gridContainer[iZone];
        grid->ChangeBCType(PHENGLEI::NO_BOUNDARY_CONDITION, PHENGLEI::INTERFACE);
    }
}

void ComputationalGridProcess::GenerateLocalGlobalInterfaceMappingByNeighborZones()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        UnstructGrid *grid = UnstructGridCast(gridContainer[iZone]);

        currentZoneID = iZone;
        currentGrid   = grid;

        this->CreateInterfaceInformation();
    }

    uint_t numberOfNeighborZone = this->neighborZoneInterfaces.size();
    for (int iZone = 0; iZone < numberOfNeighborZone; ++ iZone)
    {
        this->GenerateLocalGlobalInterfaceMappingByNeighborZones(iZone);
    }
}

void ComputationalGridProcess::CreateInterfaceInformation()
{
    UnstructGrid *unstructuredGrid = currentGrid;

    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();
   
    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    InterfaceInfo *interfaceInformation = 0;
    int numberOfInterfaces = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iFace]);

        int boundaryConditionType = bcRegion->GetBCType();
        if (boundaryConditionType == PHENGLEI::INTERFACE || boundaryConditionType < 0)
        {
            ++ numberOfInterfaces;
        }
    }
    if (numberOfInterfaces == 0)
    {
        unstructuredGrid->SetInterfaceInfo(interfaceInformation);
        return;
    }

    interfaceInformation = new InterfaceInfo(numberOfInterfaces);
    unstructuredGrid->SetInterfaceInfo(interfaceInformation);

    CreateLocalInterfaceToBoundaryFaceIndexesMapping();
}

void ComputationalGridProcess::CreateLocalInterfaceToBoundaryFaceIndexesMapping()
{
    UnstructGrid *unstructuredGrid = currentGrid;

    InterfaceInfo *interfaceInformation = unstructuredGrid->GetInterfaceInfo();
    if (! interfaceInformation) return;

    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int *localInterfaceToBoundaryFaceIndexesMapping = interfaceInformation->GetInterFace2BoundaryFace();
    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();

    int iLocalInterfaceCount  = 0;
    for (int iBoundaryFace = 0; iBoundaryFace < numberOfBoundaryFaces; ++ iBoundaryFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iBoundaryFace]);
        int boundaryConditionType = bcRegion->GetBCType();
        if (!IsInterface(boundaryConditionType))
        {
            continue;
        }

        localInterfaceToBoundaryFaceIndexesMapping[iLocalInterfaceCount] = iBoundaryFace;
        ++ iLocalInterfaceCount;
    }
}

void ComputationalGridProcess::GenerateLocalGlobalInterfaceMappingByNeighborZones(int iZone)
{
    InterfaceLinkStructure *interfaceLinkStructure = this->GetInterfaceLinkStructure();

    vector< ZoneInterface * > &neighborZoneInterfaces = this->GetNeighborZoneInterfaces();
    ZoneInterface *zoneInterfacesTopology = neighborZoneInterfaces[iZone];
    int zoneIndex = zoneInterfacesTopology->GetZoneIndex();

    interfaceLinkStructure->InitializeInterface(zoneInterfacesTopology);

    RDouble *x = zoneInterfacesTopology->GetX();
    RDouble *y = zoneInterfacesTopology->GetY();
    RDouble *z = zoneInterfacesTopology->GetZ();

    PHVectorInt2D &faceNodeIndexContainer = zoneInterfacesTopology->GetFaceNodeIndexContainer();
    
    uint_t numberOfInterfaces = zoneInterfacesTopology->GetNumberOfFaces();

    if (numberOfInterfaces == 0) 
    {
        return;
    }

    int maxNumberOfFaceNodes = -1;
    for (int iFace = 0; iFace < numberOfInterfaces; ++ iFace)
    {
        maxNumberOfFaceNodes = MAX(maxNumberOfFaceNodes, (int)faceNodeIndexContainer[iFace].size());
    }

    RDouble *xList = new RDouble[maxNumberOfFaceNodes];
    RDouble *yList = new RDouble[maxNumberOfFaceNodes];
    RDouble *zList = new RDouble[maxNumberOfFaceNodes];

    PHVectorInt1D facePointIndexes;
    PHVectorInt1D faceNodeIndexes(maxNumberOfFaceNodes);

    int iLocalInterfaceCount = 0;
    int pointCount = interfaceLinkStructure->GetNumberOfInterfacePoints();

    for (int iFace = 0; iFace < numberOfInterfaces; ++ iFace)
    {
        int faceNodeNumber = static_cast<int>(faceNodeIndexContainer[iFace].size());

        for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
        {
            faceNodeIndexes[iNode] = faceNodeIndexContainer[iFace][iNode];
        }
        facePointIndexes.resize(faceNodeNumber);

        GetFaceCoordinateList(faceNodeIndexes, faceNodeNumber, xList, yList, zList, x, y, z);
        GetCoordinateIndexList(interfaceLinkStructure, pointCount, xList, yList, zList, faceNodeNumber, facePointIndexes);
        CreateLinkInformation(facePointIndexes, zoneIndex, iLocalInterfaceCount, interfaceLinkStructure);

        ++ iLocalInterfaceCount;
    }

    delete [] xList;
    delete [] yList;
    delete [] zList;
}

void ComputationalGridProcess::ReGenerateLocalGlobalInterfaceMapping()
{
    WriteLogFile("Start ReGenerateLocalGlobalInterfaceMapping");

    InterfaceLinkStructure *interfaceLinkStructure = this->GetInterfaceLinkStructure();
    Grid **gridContainer = this->GetGridContainer();
    InitializeNewLocalGlobalInterfaceMapping();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        UnstructGrid *grid = UnstructGridCast(gridContainer[iZone]);

        currentZoneID = iZone;
        currentGrid   = grid;

        ReGenerateLocalGlobalInterfaceMapping(interfaceLinkStructure);
    }

    UpdateLocalGlobalInterfaceMapping();
    UpdateOtherTopologyTerm();
}

void ComputationalGridProcess::ReGenerateLocalGlobalInterfaceMapping(InterfaceLinkStructure *interfaceLinkStructure)
{
    UnstructGrid *unstructuredGrid = currentGrid;

    if (!ExistInterfaces(unstructuredGrid))
    {
        return;
    }

    ModifyFaceNodesIndexes(interfaceLinkStructure);
    ModifyBoundaryInformation(interfaceLinkStructure);
}

void ComputationalGridProcess::ModifyBoundaryInformation(InterfaceLinkStructure *interfaceLinkStructure)
{
    UnstructGrid *unstructuredGrid = currentGrid;

    FaceSearchingManager *faceSearchingManager = interfaceLinkStructure->GetFaceSearchingManager();
    
    PHVectorInt2D &globalChildFaceIndexes = faceSearchingManager->GetGlobalChildFaceIndexes();
    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = interfaceLinkStructure->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();
    int localInterfaceIndex = 0;

    PHVectorInt1D localInterfaceIndexToBoundaryFaceIndexMapping;

    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iFace]);

        int boundaryConditionType = bcRegion->GetBCType();
        if (IsInterface(boundaryConditionType))
        {
            localInterfaceIndexToBoundaryFaceIndexMapping.push_back(iFace);
            ++ localInterfaceIndex;
        }
    }

    int zoneIndex = unstructuredGrid->GetZoneID();
    uint_t numberOfInterfaces = localInterfaceIndexToBoundaryFaceIndexMapping.size();
    uint_t newNumberOfInterfaces = numberOfInterfaces;

    PHVectorInt1D &newLocalInterfaceIndexToGlobalInterfaceIndexMapping = this->newLocalInterfaceIndexToGlobalInterfaceIndexMappingList[currentZoneID];

    PHVectorInt1D numberOfChildFacesOfInterface;
    numberOfChildFacesOfInterface.resize(numberOfInterfaces);
    ZeroField(numberOfChildFacesOfInterface, numberOfInterfaces);

    for (int interfaceIndex = 0; interfaceIndex < numberOfInterfaces; ++ interfaceIndex)
    {
        int globalInterfaceIndex  = localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex][interfaceIndex];
        uint_t numberOfChildFaces = globalChildFaceIndexes[globalInterfaceIndex].size();

        if (numberOfChildFaces > 0)
        {
            numberOfChildFacesOfInterface[interfaceIndex] = static_cast<int>(numberOfChildFaces);
            for (int iChildFace = 0; iChildFace < numberOfChildFaces; ++ iChildFace)
            {
                int childFaceIndex = globalChildFaceIndexes[globalInterfaceIndex][iChildFace];
                newLocalInterfaceIndexToGlobalInterfaceIndexMapping.push_back(childFaceIndex);
                numberOfChildFacesOfInterface.push_back(0);
            }
            ++ newNumberOfInterfaces;
        }
        else
        {
            newLocalInterfaceIndexToGlobalInterfaceIndexMapping.push_back(globalInterfaceIndex);
        }
    }

    ResetNumberOfBoundaryCondition(interfaceLinkStructure, numberOfChildFacesOfInterface);
}

void ComputationalGridProcess::ResetNumberOfBoundaryCondition(InterfaceLinkStructure *interfaceLinkStructure, PHVectorInt1D &numberOfChildFacesOfInterface)
{
    UnstructGrid *unstructuredGrid = currentGrid;

    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();

    PHVectorInt1D &boundaryTypeVector = boundaryTypeVectorList[currentZoneID];
    vector< string > &boundaryNameVector = boundaryNameVectorList[currentZoneID];

    boundaryTypeVector.resize(0);
    boundaryNameVector.resize(0);

    int localInterfaceIndex = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iFace]);

        int boundaryConditionType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();
        if (IsInterface(boundaryConditionType))
        {
            int numberOfChildFaces = numberOfChildFacesOfInterface[localInterfaceIndex];

            if (numberOfChildFaces > 0)
            {
                for (int iChildFace = 0; iChildFace < numberOfChildFaces; ++ iChildFace)
                {
                    boundaryTypeVector.push_back(boundaryConditionType);
                    boundaryNameVector.push_back(bcName);
                }
            }
            else
            {
                boundaryTypeVector.push_back(boundaryConditionType);
                boundaryNameVector.push_back(bcName);
            }

            ++ localInterfaceIndex;
        }
        else
        {
            boundaryTypeVector.push_back(boundaryConditionType);
            boundaryNameVector.push_back(bcName);
        }
    }

    DoFurtherWorkUsingBoundaryTypeVector(interfaceLinkStructure);
}

void ComputationalGridProcess::DoFurtherWorkUsingBoundaryTypeVector(InterfaceLinkStructure *interfaceLinkStructure)
{
    UnstructGrid *unstructuredGrid = currentGrid;

    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = interfaceLinkStructure->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();
    PHVectorInt1D &newLocalInterfaceIndexToGlobalInterfaceIndexMapping = this->newLocalInterfaceIndexToGlobalInterfaceIndexMappingList[currentZoneID];

    int zoneIndex = unstructuredGrid->GetZoneID();
    localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex] = newLocalInterfaceIndexToGlobalInterfaceIndexMapping;

    ProcessGlobalToLocalInterfaceTerm(interfaceLinkStructure);
}

void ComputationalGridProcess::ProcessGlobalToLocalInterfaceTerm(InterfaceLinkStructure *interfaceLinkStructure)
{
    UnstructGrid *unstructuredGrid = currentGrid;

    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = interfaceLinkStructure->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();
    PHVectorInt1D &newLocalInterfaceIndexToGlobalInterfaceIndexMapping = this->newLocalInterfaceIndexToGlobalInterfaceIndexMappingList[currentZoneID];

    PHVectorInt1D tmp = newLocalInterfaceIndexToGlobalInterfaceIndexMapping;
    sort(tmp.begin(), tmp.end());

    int zoneIndex = unstructuredGrid->GetZoneID();
    localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex] = newLocalInterfaceIndexToGlobalInterfaceIndexMapping;

    uint_t numberOfNewInterfaces = newLocalInterfaceIndexToGlobalInterfaceIndexMapping.size();

    PHVectorInt2D &newGlobalInterfaceIndexToZoneIndexesMapping = interfaceLinkStructure->GetNewGlobalInterfaceIndexToZoneIndexesMapping();
    PHVectorInt2D &newGlobalInterfaceIndexToLocalInterfaceIndexesMapping = interfaceLinkStructure->GetNewGlobalInterfaceIndexToLocalInterfaceIndexesMapping();

    for (int localInterfaceIndex = 0; localInterfaceIndex < numberOfNewInterfaces; ++ localInterfaceIndex)
    {
        int globalInterfaceIndex = newLocalInterfaceIndexToGlobalInterfaceIndexMapping[localInterfaceIndex];

        int currentSize = static_cast<int>(newGlobalInterfaceIndexToZoneIndexesMapping.size());
        currentSize = MAX(globalInterfaceIndex + 1, currentSize);

        newGlobalInterfaceIndexToZoneIndexesMapping.resize(currentSize);
        newGlobalInterfaceIndexToLocalInterfaceIndexesMapping.resize(currentSize);

        uint_t numberOfFaceZone = newGlobalInterfaceIndexToZoneIndexesMapping[globalInterfaceIndex].size();
        if (numberOfFaceZone < 2)
        {
            newGlobalInterfaceIndexToZoneIndexesMapping[globalInterfaceIndex].push_back(zoneIndex);
            newGlobalInterfaceIndexToLocalInterfaceIndexesMapping[globalInterfaceIndex].push_back(localInterfaceIndex);
        }
    }
}

void ComputationalGridProcess::ModifyFaceNodesIndexes(InterfaceLinkStructure *interfaceLinkStructure)
{
    UnstructGrid *unstructuredGrid = currentGrid;

    int numberOfFaces = unstructuredGrid->GetNTotalFace();

    int *faceNodeNumberContainer = unstructuredGrid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = unstructuredGrid->GetFace2Node();

    PHVectorInt2D &faceNodeIndexVector = this->faceNodeIndexVectorList[currentZoneID];

    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        int numberOfFaceNodes = faceNodeNumberContainer[iFace];
        iCount += numberOfFaceNodes;
    }

    int innerFaceStartPosition = iCount;

    SetBoundaryFaceNodeIndexVector(interfaceLinkStructure);
    SetNewFaceToCellIndexesContainer(interfaceLinkStructure);

    for (int iFace = numberOfBoundaryFaces; iFace < numberOfFaces; ++ iFace)
    {
        int numberOfFaceNodes = faceNodeNumberContainer[iFace];
        PHVectorInt1D tmpVector;
        for (int iNode = 0; iNode < numberOfFaceNodes; ++ iNode)
        {
            int nodeIndex = faceNodeIndexContainer[innerFaceStartPosition ++];
            tmpVector.push_back(nodeIndex);
        }
        faceNodeIndexVector.push_back(tmpVector);
    }
}

void ComputationalGridProcess::SetNewFaceToCellIndexesContainer(InterfaceLinkStructure *interfaceLinkStructure)
{
    UnstructGrid *unstructuredGrid = currentGrid;

    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();

    PHVectorInt1D &newLeftCellIndexVector  = this->newLeftCellIndexVectorList[currentZoneID];
    PHVectorInt1D &newRightCellIndexVector = this->newRightCellIndexVectorList[currentZoneID];

    int *leftCellIndexContainer  = unstructuredGrid->GetLeftCellOfFace();
    int *rightCellIndexContainer = unstructuredGrid->GetRightCellOfFace();

    FaceSearchingManager *faceSearchingManager = interfaceLinkStructure->GetFaceSearchingManager();
    PHVectorInt2D &globalChildFaceIndexes = faceSearchingManager->GetGlobalChildFaceIndexes();
    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = interfaceLinkStructure->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

    int zoneIndex = unstructuredGrid->GetZoneID();

    int localInterfaceIndex = 0;
  
    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iFace]);

        int boundaryConditionType = bcRegion->GetBCType();

        if (IsInterface(boundaryConditionType))
        {
            int globalInterfaceIndex  = localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex][localInterfaceIndex];
            uint_t numberOfChildFaces = globalChildFaceIndexes[globalInterfaceIndex].size();

            if (numberOfChildFaces > 0)
            {
                for (int iChildFace = 0; iChildFace < numberOfChildFaces; ++ iChildFace)
                {
                    newLeftCellIndexVector.push_back(leftCellIndexContainer[iFace]);
                }
            }
            else
            {
                newLeftCellIndexVector.push_back(leftCellIndexContainer[iFace]);
            }
            ++ localInterfaceIndex;
        }
        else
        {
            newLeftCellIndexVector.push_back(leftCellIndexContainer[iFace]);
        }
    }

    uint_t newNumberOfBoundaryFaces = newLeftCellIndexVector.size();
    int numberOfCells = unstructuredGrid->GetNTotalCell();

    for (int iFace = 0; iFace < newNumberOfBoundaryFaces; ++ iFace)
    {
        newRightCellIndexVector.push_back(iFace + numberOfCells);
    }
    newRightCellIndexVector.resize(newNumberOfBoundaryFaces);

    int numberOfFaces = unstructuredGrid->GetNTotalFace();
    for (int iFace = numberOfBoundaryFaces; iFace < numberOfFaces; ++ iFace)
    {
        newLeftCellIndexVector .push_back(leftCellIndexContainer [iFace]);
        newRightCellIndexVector.push_back(rightCellIndexContainer[iFace]);
    }
}

void ComputationalGridProcess::SetBoundaryFaceNodeIndexVector(InterfaceLinkStructure *interfaceLinkStructure)
{
    UnstructGrid *unstructuredGrid = currentGrid;

    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();

    int *faceNodeNumberContainer = unstructuredGrid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = unstructuredGrid->GetFace2Node();

    PHVectorInt1D &newBoundaryLeftCellIndex = this->newBoundaryLeftCellIndexList[currentZoneID];

    int *leftCellIndexContainer = unstructuredGrid->GetLeftCellOfFace();

    FaceSearchingManager *faceSearchingManager  = interfaceLinkStructure->GetFaceSearchingManager();
    PHVectorInt2D &relativeChildFaceNodeIndexes = faceSearchingManager->GetRelativeChildFaceNodeIndexes();
    PHVectorInt2D &globalChildFaceIndexes       = faceSearchingManager->GetGlobalChildFaceIndexes();

    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = interfaceLinkStructure->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

    int zoneIndex = unstructuredGrid->GetZoneID();

    PHVectorInt2D &faceNodeIndexVector = this->faceNodeIndexVectorList[currentZoneID];
    faceNodeIndexVector.resize(0);

    int localInterfaceIndex = 0;
    int facePointStartPosition = 0;

    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iFace]);

        int boundaryConditionType = bcRegion->GetBCType();
        int numberOfFaceNodes = faceNodeNumberContainer[iFace];
        PHVectorInt1D tmpVector;

        if (IsInterface(boundaryConditionType))
        {
            int globalInterfaceIndex  = localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex][localInterfaceIndex];
            uint_t numberOfChildFaces = globalChildFaceIndexes[globalInterfaceIndex].size();

            if (numberOfChildFaces > 0)
            {
                for (int iChildFace = 0; iChildFace < numberOfChildFaces; ++ iChildFace)
                {
                    newBoundaryLeftCellIndex.push_back(leftCellIndexContainer[iFace]);

                    int childFaceIndex = globalChildFaceIndexes[globalInterfaceIndex][iChildFace];
                    uint_t numberOfChildFaceNodes = relativeChildFaceNodeIndexes[childFaceIndex].size();

                    tmpVector.resize(0);
                    for (int iNode = 0; iNode < numberOfChildFaceNodes; ++ iNode)
                    {
                        int relativeNodeIndex = facePointStartPosition + relativeChildFaceNodeIndexes[childFaceIndex][iNode];
                        int nodeIndex = faceNodeIndexContainer[relativeNodeIndex];
                        tmpVector.push_back(nodeIndex);
                    }
                    faceNodeIndexVector.push_back(tmpVector);
                }
            }
            else
            {
                newBoundaryLeftCellIndex.push_back(leftCellIndexContainer[iFace]);

                for (int iNode = 0; iNode < numberOfFaceNodes; ++ iNode)
                {
                    int nodeIndex = faceNodeIndexContainer[facePointStartPosition + iNode];
                    tmpVector.push_back(nodeIndex);
                }
                faceNodeIndexVector.push_back(tmpVector);
            }
            ++ localInterfaceIndex;
        }
        else
        {
            newBoundaryLeftCellIndex.push_back(leftCellIndexContainer[iFace]);

            for (int iNode = 0; iNode < numberOfFaceNodes; ++ iNode)
            {
                int nodeIndex = faceNodeIndexContainer[facePointStartPosition + iNode];
                tmpVector.push_back(nodeIndex);
            }
            faceNodeIndexVector.push_back(tmpVector);
        }

        facePointStartPosition += numberOfFaceNodes;
    }
}

void ComputationalGridProcess::InitializeNewLocalGlobalInterfaceMapping()
{
    InterfaceLinkStructure *interfaceLinkStructure = this->GetInterfaceLinkStructure();
    interfaceLinkStructure->InitializeNewLocalGlobalInterfaceMapping();
}

void ComputationalGridProcess::UpdateLocalGlobalInterfaceMapping()
{
    InterfaceLinkStructure *interfaceLinkStructure = this->GetInterfaceLinkStructure();
    interfaceLinkStructure->UpdateLocalGlobalInterfaceMapping();
}

void ComputationalGridProcess::UpdateOtherTopologyTerm()
{
    InterfaceLinkStructure *interfaceLinkStructure = this->GetInterfaceLinkStructure();
    Grid **gridContainer = this->GetGridContainer();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Grid *grid = gridContainer[iZone];

        currentZoneID = iZone;
        currentGrid   = UnstructGridCast(grid);

        UpdateOtherTopologyTerm(interfaceLinkStructure);
    }
}

void ComputationalGridProcess::UpdateOtherTopologyTerm(InterfaceLinkStructure *interfaceLinkStructure)
{
    UnstructGrid  *unstructuredGrid = currentGrid;
    InterfaceInfo *interfaceInformation = unstructuredGrid->GetInterfaceInfo();

    if (!interfaceInformation) return;

    UpdateOtherBoundaryTopologyTerm();
    UpdateOtherGridFaceTopologyTerm();

    int currentZoneIndex = unstructuredGrid->GetZoneID();

    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = interfaceLinkStructure->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

    uint_t numberOfNewInterfaces = localInterfaceIndexToGlobalInterfaceIndexMapping[currentZoneIndex].size();
    interfaceInformation->ReSize(static_cast<int>(numberOfNewInterfaces));

    CreateLocalInterfaceToBoundaryFaceIndexesMapping();
}

void ComputationalGridProcess::UpdateOtherGridFaceTopologyTerm()
{
    UnstructGrid *unstructuredGrid = currentGrid;
    int numberOfFaces = unstructuredGrid->GetNTotalFace();

    int *leftCellIndexContainer  = unstructuredGrid->GetLeftCellOfFace();
    int *rightCellIndexContainer = unstructuredGrid->GetRightCellOfFace();
    delete [] leftCellIndexContainer;
    delete [] rightCellIndexContainer;

    PHVectorInt1D &newLeftCellIndexVector  = this->newLeftCellIndexVectorList [currentZoneID];
    PHVectorInt1D &newRightCellIndexVector = this->newRightCellIndexVectorList[currentZoneID];

    int newNumberOfFaces = static_cast<int>(newLeftCellIndexVector.size());
    unstructuredGrid->SetNTotalFace(newNumberOfFaces);

    leftCellIndexContainer  = new int[newNumberOfFaces];
    rightCellIndexContainer = new int[newNumberOfFaces];
    unstructuredGrid->SetLeftCellOfFace(leftCellIndexContainer);
    unstructuredGrid->SetRightCellOfFace(rightCellIndexContainer);

    for (int iFace = 0; iFace < newNumberOfFaces; ++ iFace)
    {
        leftCellIndexContainer [iFace] = newLeftCellIndexVector [iFace];
        rightCellIndexContainer[iFace] = newRightCellIndexVector[iFace];
    }

    int *faceNodeNumberContainer = unstructuredGrid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = unstructuredGrid->GetFace2Node();
    delete [] faceNodeNumberContainer;    faceNodeNumberContainer = NULL;
    delete [] faceNodeIndexContainer;     faceNodeIndexContainer  = NULL;

    PHVectorInt2D &faceNodeIndexVector = faceNodeIndexVectorList[currentZoneID];
    uint_t sizeOfFaceNodeIndexContainer = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        sizeOfFaceNodeIndexContainer += faceNodeIndexVector[iFace].size();
    }
    faceNodeNumberContainer = new int[numberOfFaces];
    faceNodeIndexContainer  = new int[sizeOfFaceNodeIndexContainer];
    unstructuredGrid->SetFace2Node(faceNodeIndexContainer );
    unstructuredGrid->SetNodeNumberOfEachFace(faceNodeNumberContainer);

    int faceNodeStartPosition = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int numberOfFaceNodes = static_cast<int>(faceNodeIndexVector[iFace].size());
        faceNodeNumberContainer[iFace] = numberOfFaceNodes;

        for (int iNode = 0; iNode < numberOfFaceNodes; ++ iNode)
        {
            faceNodeIndexContainer[faceNodeStartPosition + iNode] = faceNodeIndexVector[iFace][iNode];
        }
        faceNodeStartPosition += numberOfFaceNodes;
    }
}

void ComputationalGridProcess::UpdateOtherBoundaryTopologyTerm()
{
    UnstructGrid *unstructuredGrid = currentGrid;

    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();
   
    UnstructBCSet **bcr = unstructuredGrid->GetBCRecord();
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        delete bcr[iFace];
    }
    delete [] bcr;

    PHVectorInt1D &boundaryTypeVector    = boundaryTypeVectorList[currentZoneID];
    vector< string > &boundaryNameVector = boundaryNameVectorList[currentZoneID];

    int newNumberOfBoundaryFaces = static_cast<int>(boundaryTypeVector.size());
    unstructuredGrid->SetNBoundFace(newNumberOfBoundaryFaces);

    UnstructBCSet **bcrNew = new UnstructBCSet *[newNumberOfBoundaryFaces];
 
    //! Create UnstructbcRegion.
    set<string> bcNameSet;
    for (int iFace = 0; iFace < newNumberOfBoundaryFaces; iFace++)
    {   
        bcrNew[iFace] = new UnstructBCSet;
        bcrNew[iFace]->SetKey(boundaryTypeVector[iFace]);
        bcrNew[iFace]->SetBCName(boundaryNameVector[iFace]);

        bcNameSet.insert(boundaryNameVector[iFace]);
    }

    unstructuredGrid->SetBCRecord(bcrNew);

    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());

    unstructuredGrid->CreateUnstructBCSet(nBCRegionUnstruct);

    UnstructBCSet *unstructBCSetNew = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = new int[newNumberOfBoundaryFaces];
    set<string>::iterator iter;
    int count = 0;
    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); iter++)
    {
        UnstructBC *unstructBC = new UnstructBC(count);
        unstructBCSetNew->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < newNumberOfBoundaryFaces; iFace++)
        {
            if (boundaryNameVector[iFace] == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(boundaryTypeVector[iFace]);
                bcRegionIDofBCFace[iFace] = count;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }

        }
        count++;
    }
    unstructBCSetNew->SetBCFaceInfo(bcRegionIDofBCFace);
}

void ComputationalGridProcess::MatchInterfaceTopology()
{
    WriteLogFile("Start MatchInterfaceTopology");

    InterfaceLinkStructure *interfaceLinkStructure = this->GetInterfaceLinkStructure();

    Grid **gridContainer = this->GetGridContainer();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Grid *grid = gridContainer[iZone];

        currentZoneID = iZone;
        currentGrid   = UnstructGridCast(grid);

        interfaceLinkStructure->MatchInterfaceTopology(grid);
    }
}

void ComputationalGridProcess::ReconstructInterfaceTopology()
{
    InterfaceLinkStructure *interfaceLinkStructure = this->GetInterfaceLinkStructure();
    interfaceLinkStructure->ReconstructInterfaceTopology();
}

void ComputationalGridProcess::GetNeighborZoneBoundingBox(RDouble *pmin, RDouble *pmax)
{
    pmin[0] = LARGE;
    pmin[1] = LARGE;
    pmin[2] = LARGE;

    pmax[0] = - LARGE;
    pmax[1] = - LARGE;
    pmax[2] = - LARGE;

    vector< ZoneInterface * > &neighborZoneInterfaces = this->GetNeighborZoneInterfaces();
    uint_t numberOfNeighborZones = neighborZoneInterfaces.size();
    for (int iZone = 0; iZone < numberOfNeighborZones; ++ iZone)
    {
        ZoneInterface *zoneInterface = neighborZoneInterfaces[iZone];
        RDouble *localPmin = zoneInterface->GetMinBox();
        RDouble *localPmax = zoneInterface->GetMaxBox();

        pmin[0] = MIN(pmin[0], localPmin[0]);
        pmin[1] = MIN(pmin[1], localPmin[1]);
        pmin[2] = MIN(pmin[2], localPmin[2]);

        pmax[0] = MAX(pmax[0], localPmax[0]);
        pmax[1] = MAX(pmax[1], localPmax[1]);
        pmax[2] = MAX(pmax[2], localPmax[2]);
    }
}

void ComputationalGridProcess::GetNeighborZoneMinMaxDistance(RDouble &mindis, RDouble &maxdis)
{
    mindis =   LARGE;
    maxdis = - LARGE;

    vector< ZoneInterface * > &neighborZoneInterfaces = this->GetNeighborZoneInterfaces();
    uint_t numberOfNeighborZones = neighborZoneInterfaces.size();
    for (int iZone = 0; iZone < numberOfNeighborZones; ++ iZone)
    {
        RDouble dismin, dismax;
        ZoneInterface *zoneInterface = neighborZoneInterfaces[iZone];

        zoneInterface->GetMinMaxDistance(dismin, dismax);

        mindis = MIN(mindis, dismin);
        maxdis = MAX(maxdis, dismax);
    }
}

void ComputationalGridProcess::CreateNeighborZoneInterfacesTopology()
{
    PHSPACE::WriteLogFile("Start Creating Neighbor Zone Interfaces Topology!");

    BuildZoneInterfacesOfThisProcessor();

    InitZoneInterfacesOfThisProcessor();

    SwapZoneInterfaceTopologyToNeighborZones();

    PHSPACE::WriteLogFile("End Creating Neighbor Zone Interfaces Topology!");
}

void ComputationalGridProcess::BuildZoneInterfacesOfThisProcessor()
{
    Grid **gridContainer = this->GetGridContainer();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        UnstructGrid *grid = UnstructGridCast(gridContainer[iZone]);
        BuildZoneInterfacesTopology(grid, iZone);
    }
}

void ComputationalGridProcess::BuildZoneInterfacesTopology(UnstructGrid *unstructuredGrid, int iZone)
{
    ZoneInterface *zoneInterfacesTopology = new ZoneInterface();
    zoneInterfacesTopologyList[iZone] = zoneInterfacesTopology;

    zoneInterfacesTopology->SetZoneIndex( unstructuredGrid->GetZoneID());

    RDouble *minBox = unstructuredGrid->GetMinBox();
    RDouble *maxBox = unstructuredGrid->GetMaxBox();
    RDouble *pMin = new RDouble[3];
    RDouble *pMax = new RDouble[3];
    SetField(pMin, minBox, 3);
    SetField(pMax, maxBox, 3);
    zoneInterfacesTopology->SetMinBox(pMin);
    zoneInterfacesTopology->SetMaxBox(pMax);

    RDouble disMin, disMax;
    unstructuredGrid->GetMinMaxDS(disMin, disMax);
    zoneInterfacesTopology->SetMinMaxDistance(disMin, disMax);

    int numberOfNodes = unstructuredGrid->GetNTotalNode();
    vector< bool > isInterfaceNode(numberOfNodes);

    int *faceNodeNumberContainer = unstructuredGrid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = unstructuredGrid->GetFace2Node();
    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();

    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int interfaceIndex = 0;
    int nodePosition = 0;
    for (int iBoundaryFace = 0; iBoundaryFace < numberOfBoundaryFaces; ++ iBoundaryFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iBoundaryFace]);

        int boundaryConditionType = bcRegion->GetBCType();
        int faceNodeNumber = faceNodeNumberContainer[iBoundaryFace];
        if (!IsInterface(boundaryConditionType))
        {
            nodePosition += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
        {
            int nodeIndex = faceNodeIndexContainer[nodePosition + iNode];
            isInterfaceNode[nodeIndex] = true;
        }
        ++ interfaceIndex;
        nodePosition += faceNodeNumber;
    }

    int numberOfNodesOnInterfaces = 0;
    PHVectorInt1D globalNodeIndexToLocalNodeIndex(numberOfNodes);
    SetField(globalNodeIndexToLocalNodeIndex, -1, numberOfNodes);
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (isInterfaceNode[iNode])
        {
            globalNodeIndexToLocalNodeIndex[iNode] = numberOfNodesOnInterfaces;
            ++ numberOfNodesOnInterfaces;
        }
    }

    RDouble *x    = unstructuredGrid->GetX();
    RDouble *y    = unstructuredGrid->GetY();
    RDouble *z    = unstructuredGrid->GetZ();
    RDouble *xNew = new RDouble[numberOfNodesOnInterfaces];
    RDouble *yNew = new RDouble[numberOfNodesOnInterfaces];
    RDouble *zNew = new RDouble[numberOfNodesOnInterfaces];
    zoneInterfacesTopology->SetX(xNew);
    zoneInterfacesTopology->SetY(yNew);
    zoneInterfacesTopology->SetZ(zNew);

    numberOfNodesOnInterfaces = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (isInterfaceNode[iNode])
        {
            xNew[numberOfNodesOnInterfaces] = x[iNode];
            yNew[numberOfNodesOnInterfaces] = y[iNode];
            zNew[numberOfNodesOnInterfaces] = z[iNode];

            ++ numberOfNodesOnInterfaces;
        }
    }
    zoneInterfacesTopology->SetNumberOfNodes(numberOfNodesOnInterfaces);

    PHVectorInt2D &faceNodeIndexContainerOfZoneInterfaces = zoneInterfacesTopology->GetFaceNodeIndexContainer();
    faceNodeIndexContainerOfZoneInterfaces.resize(interfaceIndex);

    interfaceIndex = 0;
    nodePosition = 0;
    for (int iBoundaryFace = 0; iBoundaryFace < numberOfBoundaryFaces; ++ iBoundaryFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iBoundaryFace]);

        int boundaryConditionType = bcRegion->GetBCType();
        int faceNodeNumber = faceNodeNumberContainer[iBoundaryFace];
        if (!IsInterface(boundaryConditionType))
        {
            nodePosition += faceNodeNumber;
            continue;
        }
        faceNodeIndexContainerOfZoneInterfaces[interfaceIndex].resize(faceNodeNumber);

        for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
        {
            int globalNodeIndex = faceNodeIndexContainer[nodePosition + iNode];
            int localNodeIndex  = globalNodeIndexToLocalNodeIndex[globalNodeIndex];
            faceNodeIndexContainerOfZoneInterfaces[interfaceIndex][iNode] = localNodeIndex;
        }

        ++ interfaceIndex;
        nodePosition += faceNodeNumber;
    }
}

void ComputationalGridProcess::InitZoneInterfacesOfThisProcessor()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneInterface *zoneInterfaceTopology = this->zoneInterfacesTopologyList[iZone];
        this->AddZoneInterfacesTopology(zoneInterfaceTopology);
    }
}

void ComputationalGridProcess::SwapZoneInterfaceLengthToNeighborZones(vector<CharVecSizeType> &neighborLength)
{
    PHSPACE::WriteLogFile("Start Swapping Length of Neighbor Zones ...");

    int nZones = PHMPI::GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = this->GetGrid(0);
        if (region->ZoneBelongToCurrentProcessor(iZone))
        {
            int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
            grid = this->GetGrid(localZoneIndex);
        }

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessorIndex    = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessorIndex = PHMPI::GetZoneProcessorID(neighborZone);

            DataContainer *dataContainer = new DataContainer();
            if (PHMPI::GetCurrentProcessorID() == sendProcessorIndex)
            {
                DataContainer *tempData = new DataContainer();
                this->CompressZoneInterfaceIntoDataContainer(tempData, iZone);

                dataContainer->MoveToBegin();
                CharVecSizeType len = tempData->Size();
                PHWrite(dataContainer, len);
                delete tempData;
            }

            PHMPI::PH_Trade(dataContainer, sendProcessorIndex, receiveProcessorIndex);

            if (PHMPI::GetCurrentProcessorID() == receiveProcessorIndex)
            {
                dataContainer->MoveToBegin();
                CharVecSizeType len;
                PHRead(dataContainer, len);
                neighborLength.push_back(len);
            }

            delete dataContainer;
        }
    }

    PHSPACE::WriteLogFile("End Swapping Length of Neighbor Zones !");
}

void ComputationalGridProcess::SwapZoneInterfaceTopologyToNeighborZones()
{
    PHSPACE::WriteLogFile("Start Swapping Zone Interface Topology To Neighbor Zones!");

    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        SwapZoneInterfaceTopologyToNeighborZonesNonBlocking();
    }
    else
    {
        SwapZoneInterfaceTopologyToNeighborZonesBlocking();
    }

    PHSPACE::WriteLogFile("End Swapping Zone Interface Topology To Neighbor Zones!");
}

void ComputationalGridProcess::SwapZoneInterfaceTopologyToNeighborZonesBlocking()
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = this->GetGrid(0);
        if (region->ZoneBelongToCurrentProcessor(iZone))
        {
            int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
            grid = this->GetGrid(localZoneIndex);
        }

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessorIndex    = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessorIndex = PHMPI::GetZoneProcessorID(neighborZone);

            this->SwapInterfaceTopologyBetweenTwoNeighboringZone(iZone, sendProcessorIndex, receiveProcessorIndex);
        }
    }
}

void ComputationalGridProcess::SwapZoneInterfaceTopologyToNeighborZonesNonBlocking()
{
    vector<CharVecSizeType> neighborLength;
    SwapZoneInterfaceLengthToNeighborZones(neighborLength);

    PHSPACE::WriteLogFile("Start Swapping Zone Interface Topology To Neighbor Zones!");
    
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer * > > receivedDataBuffer;
    vector< vector< DataContainer * > > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    int neighborCount = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = this->GetGrid(0);
        if (region->ZoneBelongToCurrentProcessor(iZone))
        {
            int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
            grid = this->GetGrid(localZoneIndex);
        }

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor    = PHMPI::GetZoneProcessorID(iZone);
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
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone].push_back(sendBuffer);
            }

            int tag = iZone;
            if (currentProcessor == sendProcessor)
            {
                this->CompressZoneInterfaceIntoDataContainer(sendBuffer, iZone);
            }
            if (sendProcessor == receiveProcessor) continue;
            
            if (currentProcessor == sendProcessor)
            {
                streamsize nlen = sendBuffer->Size();
                send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
            }
            else if (currentProcessor == receiveProcessor)
            {
                CharVecSizeType nlen = neighborLength[neighborCount++];
                receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
            }
        }
    }

    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = this->GetGrid(0);
        if (region->ZoneBelongToCurrentProcessor(iZone))
        {
            int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
            grid = this->GetGrid(localZoneIndex);
        }

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
                this->AnasysZoneInterfaceFromDataContainer(receiveData);

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
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }

    PH_Barrier();
}

void ComputationalGridProcess::SwapInterfaceTopologyBetweenTwoNeighboringZone(int zoneIndex, int sendProcessorIndex, int receiveProcessorIndex)
{
    DataContainer *dataContainer = new DataContainer();

    if (PHMPI::GetCurrentProcessorID() == sendProcessorIndex)
    {
        this->CompressZoneInterfaceIntoDataContainer(dataContainer, zoneIndex);
    }

    PHMPI::PH_Trade(dataContainer, sendProcessorIndex, receiveProcessorIndex);

    if (PHMPI::GetCurrentProcessorID() == receiveProcessorIndex)
    {
        this->AnasysZoneInterfaceFromDataContainer(dataContainer);
    }

    delete dataContainer;
}

void ComputationalGridProcess::CompressZoneInterfaceIntoDataContainer(DataContainer *&dataContainer, int globalZoneIndex)
{
    int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(globalZoneIndex);
    ZoneInterface *zoneInterface = zoneInterfacesTopologyList[localZoneIndex];
    zoneInterface->Encode(dataContainer);
}

void ComputationalGridProcess::AnasysZoneInterfaceFromDataContainer(DataContainer *&dataContainer)
{
    this->DecodeZoneInterfacesTopology(dataContainer);
}

void ComputationalGridProcess::DecodeZoneInterfacesTopology(DataContainer *&dataContainer)
{
    ZoneInterface *zoneInterfacesTopology = new ZoneInterface();
    zoneInterfacesTopology->Decode(dataContainer);

    this->AddZoneInterfacesTopology(zoneInterfacesTopology);
}

void ComputationalGridProcess::AddZoneInterfacesTopology(ZoneInterface *zoneInterfaceIn)
{
    int zoneIndex = zoneInterfaceIn->GetZoneIndex();

    vector< ZoneInterface * > &neighborZoneInterfaces = this->GetNeighborZoneInterfaces();
    for (std::size_t iZone = 0; iZone < neighborZoneInterfaces.size(); ++ iZone)
    {
        if (neighborZoneInterfaces[iZone]->GetZoneIndex() == zoneIndex)
        {
            delete zoneInterfaceIn;
            return;
        }
    }
    
    neighborZoneInterfaces.push_back(zoneInterfaceIn);
}

bool ExistInterfaces(UnstructGrid *unstructuredGrid)
{
    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();
    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();


    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iFace]);

        int boundaryConditionType = bcRegion->GetBCType();
        if (IsInterface(boundaryConditionType))
        {
            return true;
        }
    }

    return false;
}

void GetFaceCoordinateList(PHVectorInt1D &faceNodeIndexes, int numberOfPoints, RDouble *xList, RDouble *yList, RDouble *zList, RDouble *x, RDouble *y, RDouble *z)
{
    for (int iNode = 0; iNode < numberOfPoints; ++ iNode)
    {
        xList[iNode] = x[faceNodeIndexes[iNode]];
        yList[iNode] = y[faceNodeIndexes[iNode]];
        zList[iNode] = z[faceNodeIndexes[iNode]];
    }
}

void GetCoordinateIndexList(InterfaceLinkStructure *interfaceLinkStructure, int &pointCount, RDouble *xList, RDouble *yList, RDouble *zList, int numberOfPoints, PHVectorInt1D &pointIndexes)
{
    DataStruct_AdtTree< int, RDouble > *coordinateTree = &interfaceLinkStructure->GetCoordinateTree();
    RDouble diffTolerance = interfaceLinkStructure->GetTolerance();

    RDouble coordinate[3], minWindow[3], maxWindow[3];

    for (int iNode = 0; iNode < numberOfPoints; ++ iNode)
    {
        coordinate[0] = xList[iNode];
        coordinate[1] = yList[iNode];
        coordinate[2] = zList[iNode];

        minWindow[0] = coordinate[0] - diffTolerance;
        minWindow[1] = coordinate[1] - diffTolerance;
        minWindow[2] = coordinate[2] - diffTolerance;

        maxWindow[0] = coordinate[0] + diffTolerance;
        maxWindow[1] = coordinate[1] + diffTolerance;
        maxWindow[2] = coordinate[2] + diffTolerance;

        GetCoordinateIndex(coordinateTree, coordinate, minWindow, maxWindow, pointCount, pointIndexes[iNode]);
    }
}

int GetCoordinateIndex(DataStruct_AdtTree< int, RDouble > *adtTree, RDouble *coordinate, RDouble *minWindow, RDouble *maxWindow, int &pointCount, int &pointIndex)
{
    typedef DataStruct_AdtTree< int, RDouble > AdtTree;
    typedef DataStruct_AdtNode< int, RDouble > AdtNode;

    AdtNode *node;
    AdtTree::AdtNodeList nodeList;

    nodeList.resize(0);
    adtTree->FindNodesInRegion(minWindow, maxWindow, nodeList);

    if (nodeList.size() == 0)
    {
        node = new AdtNode(3, coordinate, pointCount);
        adtTree->AddNode(node);
        pointIndex = pointCount;
        ++ pointCount;
        return 0;
    }
    else
    {
        if (nodeList.size() > 1)
        {
            cout << " impossible nodeList.size() = " << nodeList.size() << "\n";
            TK_Exit::PrintDebugInfoExit("");
        }
        node = nodeList[0];
        pointIndex = node->GetData();
        return 1;
    }
}

}