#include "InterfaceLinkStructure.h"
#include "Geo_Element.h"
#include "ComputationalGrid.h"

namespace PHSPACE
{

FaceSearching::FaceSearching(const PHVectorInt1D &nodeIndexes, int faceIndex)
{
    this->faceIndex         = faceIndex;
    this->nodeIndexes       = nodeIndexes;
    this->sortedNodeIndexes = nodeIndexes;
    std::sort(sortedNodeIndexes.begin(), sortedNodeIndexes.end());
}

FaceSearching::~FaceSearching()
{

}

FaceSearchingManager::FaceSearchingManager()
{

}

FaceSearchingManager::~FaceSearchingManager()
{
    referenceFacesSet.clear();
    for (std::size_t i = 0; i < referenceFacesArray.size(); ++ i)
    {
        delete referenceFacesArray[i];
    }
}

void FaceSearchingManager::AddFace(const PHVectorInt1D &faceNodeIndexes)
{
    int faceIndex = static_cast<int>(referenceFacesArray.size());
    FaceSearching *faceSearching = new FaceSearching(faceNodeIndexes, faceIndex);

    set<FaceSearching *, CompareFaceByMethod1>::iterator iter = referenceFacesSet.find(faceSearching);

    if (iter == referenceFacesSet.end())
    {
        referenceFacesSet.insert(faceSearching);
        referenceFacesArray.push_back(faceSearching);
    }
    else
    {
        delete faceSearching;
    }
}

void FaceSearchingManager::ComputeNewFaceIndex()
{
    int numberOfFaces = static_cast<int>(referenceFacesArray.size());
    if (numberOfFaces == 0)
    {
        return;
    }

    statusContainerOfFaces.resize(numberOfFaces);
    SetField(statusContainerOfFaces, INVALID_INDEX, numberOfFaces);

    childFaceIndexes.resize(numberOfFaces);
    relativeChildFaceNodeIndexes.resize(numberOfFaces);

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        FaceSearching *faceSearching = referenceFacesArray[iFace];
        SplitQuadrilateralToTriangle(faceSearching);
    }
}

set<FaceSearching *, CompareFaceByMethod1>::iterator FaceSearchingManager::FindFace(FaceSearching *faceSearching)
{
    return referenceFacesSet.find(faceSearching);
}

void FaceSearchingManager::GetLocalTriangleIndexes(PHVectorInt2D &localTriangleIndexesContainer)
{
    PHVectorInt1D localIndexes1(3), localIndexes2(3);
    localIndexes1[0] = 0;
    localIndexes1[1] = 1;
    localIndexes1[2] = 2;

    localIndexes2[0] = 2;
    localIndexes2[1] = 3;
    localIndexes2[2] = 0;

    localTriangleIndexesContainer.push_back(localIndexes1);
    localTriangleIndexesContainer.push_back(localIndexes2);

    localIndexes1[0] = 3;
    localIndexes1[1] = 0;
    localIndexes1[2] = 1;

    localIndexes2[0] = 1;
    localIndexes2[1] = 2;
    localIndexes2[2] = 3;

    localTriangleIndexesContainer.push_back(localIndexes1);
    localTriangleIndexesContainer.push_back(localIndexes2);
}

void FaceSearchingManager::GetTriangleIndexes(FaceSearching *parentFaceSearching, PHVectorInt2D &localTriangleIndexesContainer, PHVectorInt2D &triangleIndexesContainer)
{
    const PHVectorInt1D &originalFaceNodeIndexes = parentFaceSearching->GetNodeIndexes();

    triangleIndexesContainer.resize(localTriangleIndexesContainer.size());
    uint_t localTriangleIndexesContainerSize = localTriangleIndexesContainer.size();
    for (int i = 0; i < localTriangleIndexesContainerSize; ++ i)
    {
        PHVectorInt1D &triangleIndexes = triangleIndexesContainer[i];
        for (int iNode = 0; iNode < 3; ++ iNode)
        {
            triangleIndexes.push_back(originalFaceNodeIndexes[localTriangleIndexesContainer[i][iNode]]);
        }
    }
}

void FaceSearchingManager::SplitQuadrilateralToTriangle(FaceSearching *parentFaceSearching)
{
    uint_t numberOfNodes = parentFaceSearching->GetNumberOfNodes();
    if (numberOfNodes <= 3) return;

    PHVectorInt2D localTriangleIndexesContainer;
    GetLocalTriangleIndexes(localTriangleIndexesContainer);

    PHVectorInt2D triangleIndexesContainer;
    GetTriangleIndexes(parentFaceSearching, localTriangleIndexesContainer, triangleIndexesContainer);
    uint_t triangleIndexesContainerSize = triangleIndexesContainer.size();
    for (int i = 0; i < triangleIndexesContainerSize; ++ i)
    {
        PHVectorInt1D &triangleIndexes = triangleIndexesContainer[i];

        FaceSearching *faceSearching = new FaceSearching(triangleIndexes);
        set<FaceSearching *, CompareFaceByMethod1>::iterator iter = FindFace(faceSearching);

        if (iter != referenceFacesSet.end())
        {
            int faceIndex = (*iter)->GetFaceIndex();

            int parentFaceIndex = parentFaceSearching->GetFaceIndex();
            childFaceIndexes[parentFaceIndex].push_back(faceIndex);

            relativeChildFaceNodeIndexes[faceIndex] = localTriangleIndexesContainer[i];
        }

        delete faceSearching;
    }
}

InterfaceLinkStructure::InterfaceLinkStructure(RDouble *pmin, RDouble *pmax, RDouble tolerance, int nZones) : coordinateTree(3, pmin, pmax)
{
    this->tolerance = tolerance;
    localInterfaceIndexToGlobalInterfaceIndexMapping.resize(nZones);

    faceSearchingManager = new FaceSearchingManager();
}

InterfaceLinkStructure::~InterfaceLinkStructure()
{
    delete faceSearchingManager;
}

void InterfaceLinkStructure::CreateLinkInformation(PHVectorInt1D &facePointIndexes, int zoneIndex, int iLocalInterfaceCount)
{
    set< DataStruct_Sort< PHVectorInt1D > > &referenceInterfaceListForSearching = this->GetReferenceInterfaceListForSearching();
    PHVectorInt2D &globalInterfaceIndexToZoneIndexesMapping           = this->GetGlobalInterfaceIndexToZoneIndexesMapping();
    PHVectorInt2D &globalInterfaceIndexToLocalInterfaceIndexesMapping = this->GetGlobalInterfaceIndexToLocalInterfaceIndexesMapping();

    int numberOfGlobalInterfaces = static_cast<int>(globalInterfaceIndexToZoneIndexesMapping.size());

    AddFace(facePointIndexes);

    sort(facePointIndexes.begin(), facePointIndexes.end());
    DataStruct_Sort< PHVectorInt1D > face(facePointIndexes, numberOfGlobalInterfaces);
    set< DataStruct_Sort< PHVectorInt1D > >::iterator iter = referenceInterfaceListForSearching.find(face);

    if (iter == referenceInterfaceListForSearching.end())
    {
        PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = this->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

        referenceInterfaceListForSearching.insert(face);
        localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex][iLocalInterfaceCount] = numberOfGlobalInterfaces;
        PHVectorInt1D zoneIndexVector;
        PHVectorInt1D localInterfaceIndexVector;
        zoneIndexVector.push_back(zoneIndex);
        localInterfaceIndexVector.push_back(iLocalInterfaceCount);

        globalInterfaceIndexToZoneIndexesMapping.push_back(zoneIndexVector);
        globalInterfaceIndexToLocalInterfaceIndexesMapping.push_back(localInterfaceIndexVector);
    }
    else
    {
        PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = this->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

        int globalInterfaceIndex = iter->index;
        localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex][iLocalInterfaceCount] = globalInterfaceIndex;
        globalInterfaceIndexToZoneIndexesMapping[globalInterfaceIndex].push_back(zoneIndex);
        globalInterfaceIndexToLocalInterfaceIndexesMapping[globalInterfaceIndex].push_back(iLocalInterfaceCount);
        referenceInterfaceListForSearching.erase(iter);
    }
}

void InterfaceLinkStructure::CreateLinkInformationBySet(PHVectorInt1D &facePointIndexes, int zoneIndex, int iLocalInterfaceCount)
{
    set < DataStruct_Sort< set< int > > > &referenceInterfaceSetForSearching = this->GetReferenceInterfaceSetForSearching();
    PHVectorInt2D &globalInterfaceIndexToZoneIndexesMapping = this->GetGlobalInterfaceIndexToZoneIndexesMapping();
    PHVectorInt2D &globalInterfaceIndexToLocalInterfaceIndexesMapping = this->GetGlobalInterfaceIndexToLocalInterfaceIndexesMapping();

    int numberOfGlobalInterfaces = static_cast<int>(globalInterfaceIndexToZoneIndexesMapping.size());

    set< int > faceSet;
    uint_t facePointIndexesSize = facePointIndexes.size();
    for (int iNode = 0; iNode < facePointIndexesSize; ++ iNode)
    {
        faceSet.insert(facePointIndexes[iNode]);
    }

    DataStruct_Sort< set<int> > sortFaceSet(faceSet, numberOfGlobalInterfaces);

    set < DataStruct_Sort< set< int > > > ::iterator iter = referenceInterfaceSetForSearching.find(sortFaceSet);
    if (iter == referenceInterfaceSetForSearching.end())
    {
        PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = this->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

        referenceInterfaceSetForSearching.insert(sortFaceSet);
        localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex][iLocalInterfaceCount] = numberOfGlobalInterfaces;

        PHVectorInt1D zoneIndexVector;
        PHVectorInt1D localInterfaceIndexVector;
        zoneIndexVector.push_back(zoneIndex);
        localInterfaceIndexVector.push_back(iLocalInterfaceCount);

        globalInterfaceIndexToZoneIndexesMapping.push_back(zoneIndexVector);
        globalInterfaceIndexToLocalInterfaceIndexesMapping.push_back(localInterfaceIndexVector);
    }
    else
    {
        PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = this->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

        int globalInterfaceIndex = iter->index;
        localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex][iLocalInterfaceCount] = globalInterfaceIndex;
        globalInterfaceIndexToZoneIndexesMapping [globalInterfaceIndex].push_back(zoneIndex);
        globalInterfaceIndexToLocalInterfaceIndexesMapping[globalInterfaceIndex].push_back(iLocalInterfaceCount);
        referenceInterfaceSetForSearching.erase(iter);
    }
}

void InterfaceLinkStructure::MatchInterfaceTopology(Grid *grid)
{
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;
    
    int numberOfInterfaces = interfaceInformation->GetNIFace();
    int currentZoneGlobalIndex = grid->GetZoneID();

    int nPeoridic = 0;
    
    PHVectorInt2D &globalInterfaceIndexToZoneIndexesMapping           = this->GetGlobalInterfaceIndexToZoneIndexesMapping();
    PHVectorInt2D &globalInterfaceIndexToLocalInterfaceIndexesMapping = this->GetGlobalInterfaceIndexToLocalInterfaceIndexesMapping();
    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping   = this->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();

    uint_t numberOfNewInterfaces = localInterfaceIndexToGlobalInterfaceIndexMapping[currentZoneGlobalIndex].size();
    
    int *neighborZoneIndexContainer               = interfaceInformation->GetInterFace2ZoneID();
    int *neighborZoneLocalInterfaceIndexContainer = interfaceInformation->GetInterFace2InterFaceID();

    UnstructGrid *unstrGrid = UnstructGridCast(grid);
    unstrGrid->ComputeMetrics();
    RDouble *xfc = unstrGrid->GetFaceCenterX();
    RDouble *yfc = unstrGrid->GetFaceCenterY();
    RDouble *zfc = unstrGrid->GetFaceCenterZ();

    for (int iLocalInterface = 0; iLocalInterface < numberOfNewInterfaces; ++ iLocalInterface)
    {
        int  globalInterfaceIndex   = localInterfaceIndexToGlobalInterfaceIndexMapping[currentZoneGlobalIndex][iLocalInterface];
        bool findRequestedInterface = false;
        uint_t numberOfZonesBelongToCurrentInterface = globalInterfaceIndexToZoneIndexesMapping[globalInterfaceIndex].size();

        if (numberOfZonesBelongToCurrentInterface != 2)
        {
            if (numberOfZonesBelongToCurrentInterface > 2)
            {
                cout << " More than two faces coincide\n";
            }
            else
            {
                cout << " Less than two faces coincide\n";
            }
            cout << " Current ZoneIndex  = " << currentZoneGlobalIndex << "\n";
            cout << " numberOfZonesBelongToCurrentInterface = " << numberOfZonesBelongToCurrentInterface << "\n";
            cout << " LocalInterface Index = " << iLocalInterface << " numberOfInterfaces = " << numberOfInterfaces << "\n";
        }
        

        for (int iZoneBelongToCurrentInterface = 0; iZoneBelongToCurrentInterface < numberOfZonesBelongToCurrentInterface; ++ iZoneBelongToCurrentInterface)
        {
            int candidateNeighborZoneIndex               = globalInterfaceIndexToZoneIndexesMapping[globalInterfaceIndex][iZoneBelongToCurrentInterface];
            int candidateNeighborZonelocalInterfaceIndex = globalInterfaceIndexToLocalInterfaceIndexesMapping[globalInterfaceIndex][iZoneBelongToCurrentInterface];

            if ((candidateNeighborZoneIndex != currentZoneGlobalIndex) || (candidateNeighborZonelocalInterfaceIndex != iLocalInterface))
            {
                neighborZoneIndexContainer              [iLocalInterface] = candidateNeighborZoneIndex;
                neighborZoneLocalInterfaceIndexContainer[iLocalInterface] = candidateNeighborZonelocalInterfaceIndex;
                findRequestedInterface = true;
                break;
            }
        }

        if (!findRequestedInterface)
        {
            int bcFaceID = interfaceInformation->GetInterFace2BoundaryFace()[globalInterfaceIndex];
            cout << "   Face center: " << xfc[bcFaceID] << " " << yfc[bcFaceID] << " " << zfc[bcFaceID] << "  ";
            cout << "iLocalInterface = " << iLocalInterface << ", neighbor zone can not been found!!!\n";
            ++ nPeoridic;
            TK_Exit::ExceptionExit("Wrong neighbor zone interface!");
        }
    }
}

void InterfaceLinkStructure::ReconstructInterfaceTopology()
{
    FaceSearchingManager *faceSearchingManager = GetFaceSearchingManager();
    faceSearchingManager->ComputeNewFaceIndex();
}

void InterfaceLinkStructure::InitializeInterface(Grid *grid)
{
    int zoneIndex = grid->GetZoneID();
    int numberOfInterfaces = grid->GetNIFace();

    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = this->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();
    localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex].resize(numberOfInterfaces);
}

void InterfaceLinkStructure::InitializeInterface(ZoneInterface *zoneInterfaceTopology)
{
    int zoneIndex = zoneInterfaceTopology->GetZoneIndex();
    uint_t numberOfInterfaces = zoneInterfaceTopology->GetNumberOfFaces();

    PHVectorInt2D &localInterfaceIndexToGlobalInterfaceIndexMapping = this->GetLocalInterfaceIndexToGlobalInterfaceIndexMapping();
    localInterfaceIndexToGlobalInterfaceIndexMapping[zoneIndex].resize(numberOfInterfaces);
}

void InterfaceLinkStructure::InitializeNewLocalGlobalInterfaceMapping()
{
    PHVectorInt2D &globalInterfaceIndexToZoneIndexesMapping           = this->GetGlobalInterfaceIndexToZoneIndexesMapping();
    PHVectorInt2D &globalInterfaceIndexToLocalInterfaceIndexesMapping = this->GetGlobalInterfaceIndexToLocalInterfaceIndexesMapping();

    PHVectorInt2D &newGlobalInterfaceIndexToZoneIndexesMapping           = this->GetNewGlobalInterfaceIndexToZoneIndexesMapping();
    PHVectorInt2D &newGlobalInterfaceIndexToLocalInterfaceIndexesMapping = this->GetNewGlobalInterfaceIndexToLocalInterfaceIndexesMapping ();

    newGlobalInterfaceIndexToZoneIndexesMapping           = globalInterfaceIndexToZoneIndexesMapping;
    newGlobalInterfaceIndexToLocalInterfaceIndexesMapping = globalInterfaceIndexToLocalInterfaceIndexesMapping;
}

void InterfaceLinkStructure::UpdateLocalGlobalInterfaceMapping()
{
    PHVectorInt2D &globalInterfaceIndexToZoneIndexesMapping           = this->GetGlobalInterfaceIndexToZoneIndexesMapping();
    PHVectorInt2D &globalInterfaceIndexToLocalInterfaceIndexesMapping = this->GetGlobalInterfaceIndexToLocalInterfaceIndexesMapping();

    PHVectorInt2D &newGlobalInterfaceIndexToZoneIndexesMapping           = this->GetNewGlobalInterfaceIndexToZoneIndexesMapping();
    PHVectorInt2D &newGlobalInterfaceIndexToLocalInterfaceIndexesMapping = this->GetNewGlobalInterfaceIndexToLocalInterfaceIndexesMapping ();

    globalInterfaceIndexToZoneIndexesMapping           = newGlobalInterfaceIndexToZoneIndexesMapping;
    globalInterfaceIndexToLocalInterfaceIndexesMapping = newGlobalInterfaceIndexToLocalInterfaceIndexesMapping;
}

void InterfaceLinkStructure::InitializeStatusContainerOfGlobalFaces()
{
    PHVectorInt1D &statusContainerOfGlobalFaces = GetStatusContainerOfGlobalFaces();

    PHVectorInt2D &globalInterfaceIndexToLocalInterfaceIndexesMapping = this->GetGlobalInterfaceIndexToLocalInterfaceIndexesMapping();
    int numberOfTotalInterfaces = static_cast<int>(globalInterfaceIndexToLocalInterfaceIndexesMapping.size());

    statusContainerOfGlobalFaces.resize(numberOfTotalInterfaces);
    SetField(statusContainerOfGlobalFaces, INVALID_INDEX, numberOfTotalInterfaces);
}

void InterfaceLinkStructure::AddFace(const PHVectorInt1D &facePointIndexes)
{
    FaceSearchingManager *faceSearchingManager = GetFaceSearchingManager();
    faceSearchingManager->AddFace(facePointIndexes);
}

void CreateLinkInformation(PHVectorInt1D &facePointIndexes, int zoneIndex, int iLocalInterfaceCount, InterfaceLinkStructure *interfaceLinkStructure)
{
    interfaceLinkStructure->CreateLinkInformation(facePointIndexes, zoneIndex, iLocalInterfaceCount);
}

}
