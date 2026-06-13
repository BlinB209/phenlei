#include "Mesh_ConstructInterfaceLink.h"
#include "PHMpi.h"
#include "Geo_SimpleBC.h"
#include "PHHeader.h"
#include "GridType.h"
#include "Glb_Dimension.h"
#include "Pre_GridConversion.h"
#include "TK_Log.h"

namespace PHSPACE
{

Mesh_ConstructInterfaceLink::Mesh_ConstructInterfaceLink(Grid **gridIn, int nZonesInCurrentProcessIn)
{
    this->nZonesInCurrentProcess = nZonesInCurrentProcessIn;
    this->gridContainer          = gridIn;
    this->nZonesGlobal           = 0;
    this->zoneProc               = 0;
    this->localZoneIndex         = 0;

    this->zoneInterfacesTopologyList.resize(0);
}

Mesh_ConstructInterfaceLink::~Mesh_ConstructInterfaceLink()
{
    if (zoneProc)
    {
        delete [] zoneProc; zoneProc = nullptr;
    }

    if (localZoneIndex)
    {
        delete [] localZoneIndex; localZoneIndex = nullptr;
    }

    for (int iZone = 0; iZone < zoneInterfacesTopologyList.size(); ++ iZone)
    {
        delete zoneInterfacesTopologyList[iZone];
    }
}

void Mesh_ConstructInterfaceLink::Run()
{
    ReorderGridIndex();

    ChangeNoBoundaryCondition();

    BuildZoneInterfacesInfo();

    ConstructInterfaceLink();
}

void Mesh_ConstructInterfaceLink::CollectUnsGridInterfaceInfo(Grid *grid, DataContainer *cdata)
{
    UnstructGrid *unsGrid = UnstructGridCast(grid);

    int gridIndex = unsGrid->GetZoneID();
    PHWrite(cdata, gridIndex);

    int nIFace = unsGrid->CompNIFace();
    PHWrite(cdata, nIFace);

    InterfaceInfo *interfaceInfo = 0;
    if (0 == nIFace)
    {
        unsGrid->SetInterfaceInfo(interfaceInfo);
        return;
    }

    interfaceInfo = new InterfaceInfo(nIFace);
    unsGrid->SetInterfaceInfo(interfaceInfo);
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

    RDouble *minBox = unsGrid->GetMinBox();
    RDouble *maxBox = unsGrid->GetMaxBox();
    PHWrite(cdata, minBox, 3);
    PHWrite(cdata, maxBox, 3);

    RDouble disMin, disMax;
    unsGrid->GetMinMaxDS(disMin, disMax);
    PHWrite(cdata, disMin);
    PHWrite(cdata, disMax);

    int numberOfNodes = unsGrid->GetNTotalNode();
    vector< bool > isInterfaceNode(numberOfNodes);
    SetField(isInterfaceNode, false, numberOfNodes);

    int *faceNodeNumberContainer = unsGrid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = unsGrid->GetFace2Node();
    int  numberOfBoundaryFaces   = unsGrid->GetNBoundFace();

    UnstructBCSet **bcRegion = unsGrid->GetBCRecord();

    int *nodeNumOfFaceLocal = new int[nIFace];
    int faceNodeLocalSize = 0;
    int interfaceIndex = 0;
    int nodePosition = 0;
    for (int iBFace = 0; iBFace < numberOfBoundaryFaces; ++ iBFace)
    {
        int bcType = bcRegion[iBFace]->GetKey();
        int faceNodeNumber = faceNodeNumberContainer[iBFace];
        if (!IsInterface(bcType))
        {
            nodePosition += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
        {
            int nodeIndex = faceNodeIndexContainer[nodePosition + iNode];
            isInterfaceNode[nodeIndex] = true;
        }

        interFace2BoundaryFace[interfaceIndex] = iBFace;
        nodeNumOfFaceLocal[interfaceIndex] = faceNodeNumber;
        faceNodeLocalSize += faceNodeNumber;
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

    RDouble *x    = unsGrid->GetX();
    RDouble *y    = unsGrid->GetY();
    RDouble *z    = unsGrid->GetZ();
    RDouble *xNew = new RDouble[numberOfNodesOnInterfaces];
    RDouble *yNew = new RDouble[numberOfNodesOnInterfaces];
    RDouble *zNew = new RDouble[numberOfNodesOnInterfaces];

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

    PHWrite(cdata, numberOfNodesOnInterfaces);
    PHWrite(cdata, xNew, numberOfNodesOnInterfaces);
    PHWrite(cdata, yNew, numberOfNodesOnInterfaces);
    PHWrite(cdata, zNew, numberOfNodesOnInterfaces);
    delete [] xNew; xNew = nullptr;
    delete [] yNew; yNew = nullptr;
    delete [] zNew; zNew = nullptr;

    int *face2NodeLocal = new int[faceNodeLocalSize];
    int nodePos = 0;
    nodePosition = 0;
    for (int iBFace = 0; iBFace < numberOfBoundaryFaces; ++ iBFace)
    {
        int bcType = bcRegion[iBFace]->GetKey();
        int faceNodeNumber = faceNodeNumberContainer[iBFace];
        if (!IsInterface(bcType))
        {
            nodePosition += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
        {
            int globalNodeIndex = faceNodeIndexContainer[nodePosition + iNode];
            int localNodeIndex  = globalNodeIndexToLocalNodeIndex[globalNodeIndex];
            face2NodeLocal[nodePos] = localNodeIndex;
            nodePos ++;
        }

        nodePosition += faceNodeNumber;
    }

    PHWrite(cdata, nodeNumOfFaceLocal, nIFace);
    PHWrite(cdata, face2NodeLocal, faceNodeLocalSize);
    delete [] nodeNumOfFaceLocal;    nodeNumOfFaceLocal = nullptr;
    delete [] face2NodeLocal;    face2NodeLocal = nullptr;
}

void Mesh_ConstructInterfaceLink::DecodeInterfaceInfo(ZoneInterface *zoneInterface, DataContainer *cdata)
{
    cdata->MoveToBegin();

    int zoneIndex = -1;
    PHRead(cdata, zoneIndex);
    zoneInterface->SetZoneIndex(zoneIndex);

    int nIFace = 0;
    PHRead(cdata, nIFace);

    RDouble *minBox = new RDouble[3];
    RDouble *maxBox = new RDouble[3];
    PHRead(cdata, minBox, 3);
    PHRead(cdata, maxBox, 3);
    zoneInterface->SetMinBox(minBox);
    zoneInterface->SetMaxBox(maxBox);

    RDouble disMin, disMax;
    PHRead(cdata, disMin);
    PHRead(cdata, disMax);
    zoneInterface->SetMinMaxDistance(disMin, disMax);

    int numberOfNodesOnInterfaces = 0;
    PHRead(cdata, numberOfNodesOnInterfaces);
    zoneInterface->SetNumberOfNodes(numberOfNodesOnInterfaces);

    RDouble *xCoor = new RDouble[numberOfNodesOnInterfaces];
    RDouble *yCoor = new RDouble[numberOfNodesOnInterfaces];
    RDouble *zCoor = new RDouble[numberOfNodesOnInterfaces];
    PHRead(cdata, xCoor, numberOfNodesOnInterfaces);
    PHRead(cdata, yCoor, numberOfNodesOnInterfaces);
    PHRead(cdata, zCoor, numberOfNodesOnInterfaces);
    zoneInterface->SetX(xCoor);
    zoneInterface->SetY(yCoor);
    zoneInterface->SetZ(zCoor);

    int *nodeNumOfFace = new int[nIFace];
    PHRead(cdata, nodeNumOfFace, nIFace);

    int faceNodeSize = 0;
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        faceNodeSize += nodeNumOfFace[iFace];
    }

    int *face2Node = new int[faceNodeSize];
    PHRead(cdata, face2Node, faceNodeSize);

    PHVectorInt2D &faceNodeIndexContainerOfZoneInterfaces = zoneInterface->GetFaceNodeIndexContainer();
    faceNodeIndexContainerOfZoneInterfaces.resize(nIFace);

    int nodepos = 0;
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        int nodeNumber = nodeNumOfFace[iFace];
        faceNodeIndexContainerOfZoneInterfaces[iFace].resize(nodeNumber);

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            faceNodeIndexContainerOfZoneInterfaces[iFace][iNode] = face2Node[nodepos];
            nodepos ++;
        }
    }

    delete [] nodeNumOfFace; nodeNumOfFace = nullptr;
    delete [] face2Node; face2Node = nullptr;
}

void Mesh_ConstructInterfaceLink::BuildZoneInterfacesInfo()
{
    int myid  = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < nZonesGlobal; ++ iZone)
    {
        DataContainer *cdata = new DataContainer();
        cdata->MoveToBegin();

        int sendProc = zoneProc[iZone];
        int tag = iZone;

        if (sendProc == myid)
        {
            int localGridIndex = this->localZoneIndex[iZone];
            int gridType = gridContainer[localGridIndex]->Type();
            if (gridType == PHSPACE::UNSTRUCTGRID)
            {
                CollectUnsGridInterfaceInfo(gridContainer[localGridIndex], cdata);
            }
            else
            {

            }

        }

        PHMPI::PH_Bcast(cdata, sendProc, tag);
        cdata->MoveToBegin();

        int neighborZoneIndex = -1;
        PHRead(cdata, neighborZoneIndex);

        int nIFace = 0;
        PHRead(cdata, nIFace);

        if (0 == nIFace)
        {
            delete cdata;    cdata = nullptr;
            continue;
        }

        RDouble *neighborMinBox = new RDouble[3];
        RDouble *neighborMaxBox = new RDouble[3];
        PHRead(cdata, neighborMinBox, 3);
        PHRead(cdata, neighborMaxBox, 3);

        RDouble disMin, disMax;
        PHRead(cdata, disMin);
        PHRead(cdata, disMax);

        RDouble dist = half * disMin;

        bool currentZoneIsNeighbor = false;
        for (int iZone = 0; iZone < nZonesInCurrentProcess; ++ iZone)
        {
            RDouble *minBox = gridContainer[iZone]->GetMinBox();
            RDouble *maxBox = gridContainer[iZone]->GetMaxBox();

            if (!(neighborMinBox[0] > (maxBox[0] + dist) || neighborMaxBox[0] < (minBox[0] - dist) ||
                  neighborMinBox[1] > (maxBox[1] + dist) || neighborMaxBox[1] < (minBox[1] - dist) ||
                  neighborMinBox[2] > (maxBox[2] + dist) || neighborMaxBox[2] < (minBox[2] - dist)))
            {
                currentZoneIsNeighbor = true;
                break;
            }
        }
        delete [] neighborMinBox; neighborMinBox = nullptr;
        delete [] neighborMaxBox; neighborMaxBox = nullptr;

        if (!currentZoneIsNeighbor)
        {
            delete cdata;    cdata = nullptr;
            continue;
        }

        ZoneInterface *zoneInterfacesTopology = new ZoneInterface();
        DecodeInterfaceInfo(zoneInterfacesTopology, cdata);

        zoneInterfacesTopologyList.push_back(zoneInterfacesTopology);
        delete cdata;    cdata = nullptr;
    }

}

void Mesh_ConstructInterfaceLink::BuildGridLink(int iZone, LinkStruct *link)
{
    ZoneInterface *currentZoneInterfaceInfo = zoneInterfacesTopologyList[iZone];

    int zoneIndex = currentZoneInterfaceInfo->GetZoneIndex();
    int nIFace    = currentZoneInterfaceInfo->GetNumberOfFaces();

    VVInt &facemap = link->GetFaceMap();
    facemap[zoneIndex].resize(nIFace);

    RDouble *x = currentZoneInterfaceInfo->GetX();
    RDouble *y = currentZoneInterfaceInfo->GetY();
    RDouble *z = currentZoneInterfaceInfo->GetZ();

    PHVectorInt2D &faceNodeIndexContainerOfZoneInterfaces = currentZoneInterfaceInfo->GetFaceNodeIndexContainer();

    int maxlist = 0;
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        int nodeNum = faceNodeIndexContainerOfZoneInterfaces[iFace].size();
        maxlist = MAX(maxlist, nodeNum);
    }

    RDouble *xlist = new RDouble[maxlist];
    RDouble *ylist = new RDouble[maxlist];
    RDouble *zlist = new RDouble[maxlist];

    vector<int> index(maxlist), pindex(maxlist);

    RDouble                          tol       = link->GetTolerance();
    VVInt                           &zoneid    = link->GetZoneID();
    VVInt                           &faceid    = link->GetFaceID();
    LinkStruct::AdtTree             &coor_tree = link->GetCoordinateTree();
    set < DataStruct_Sort< VInt > > &facelist  = link->GetFaceList();

    uint_t fcount = zoneid.size();
    int    pcount = coor_tree.GetNodeNum();
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        int nodeNum = faceNodeIndexContainerOfZoneInterfaces[iFace].size();
        for (int iNode = 0; iNode < nodeNum; ++ iNode)
        {
            index[iNode] = faceNodeIndexContainerOfZoneInterfaces[iFace][iNode];
        }
        pindex.resize(nodeNum);

        GetFaceCoorList(index, nodeNum, xlist, ylist, zlist, x, y, z);
        GetCoorIndexList(&coor_tree, tol, pcount, xlist, ylist, zlist, nodeNum, pindex);
        Create_Link_Info(pindex, zoneIndex, iFace, fcount, facelist, zoneid, faceid, link);
    }

    delete [] xlist; xlist = nullptr;
    delete [] ylist; ylist = nullptr;
    delete [] zlist; zlist = nullptr;
}

void Mesh_ConstructInterfaceLink::ConstructInterfaceLink()
{
    RDouble mindis, maxdis, pmin[3], pmax[3];

    GetBoundingBox(pmin, pmax);
    GetMinMaxDS(mindis, maxdis);

    if (GetDim() == TWO_D)
    {
        pmin[2] -= mindis;
        pmax[2] += maxdis;
    }

    LinkStruct link(pmin, pmax, mindis / 2.0, this->nZonesGlobal);
    for (int iZone = 0; iZone < zoneInterfacesTopologyList.size(); ++ iZone)
    {
        BuildGridLink(iZone, &link);
    }

    for (int iZone = 0; iZone < nZonesInCurrentProcess; ++ iZone)
    {
        MatchInterface(gridContainer[iZone], &link);
    }
}

void Mesh_ConstructInterfaceLink::GetBoundingBox(RDouble *pmin, RDouble *pmax)
{
    pmin[0] = LARGE;
    pmin[1] = LARGE;
    pmin[2] = LARGE;

    pmax[0] = - LARGE;
    pmax[1] = - LARGE;
    pmax[2] = - LARGE;

    for (int iZone = 0; iZone < zoneInterfacesTopologyList.size(); ++ iZone)
    {
        RDouble *local_pmin = zoneInterfacesTopologyList[iZone]->GetMinBox();
        RDouble *local_pmax = zoneInterfacesTopologyList[iZone]->GetMaxBox();

        pmin[0] = MIN(pmin[0], local_pmin[0]);
        pmin[1] = MIN(pmin[1], local_pmin[1]);
        pmin[2] = MIN(pmin[2], local_pmin[2]);

        pmax[0] = MAX(pmax[0], local_pmax[0]);
        pmax[1] = MAX(pmax[1], local_pmax[1]);
        pmax[2] = MAX(pmax[2], local_pmax[2]);
    }
}

void Mesh_ConstructInterfaceLink::GetMinMaxDS(RDouble &mindis, RDouble &maxdis)
{
    mindis =   LARGE;
    maxdis = - LARGE;

    for (int iZone = 0; iZone < zoneInterfacesTopologyList.size(); ++ iZone)
    {
        RDouble dismin, dismax;
        zoneInterfacesTopologyList[iZone]->GetMinMaxDistance(dismin, dismax);

        mindis = MIN(mindis, dismin);
        maxdis = MAX(maxdis, dismax);
    }
}

void Mesh_ConstructInterfaceLink::ReorderGridIndex()
{
    int nProc = PHMPI::GetNumberOfProcessor();
    int myid  = PHMPI::GetCurrentProcessorID();

    int *nZonesOfEachProcessor = new int[nProc]();
    PH_AllGather(&nZonesInCurrentProcess, 1, nZonesOfEachProcessor, 1);

    int zoneInsexStart = 0;
    for (int iProc = 0; iProc < myid; ++ iProc)
    {
        zoneInsexStart += nZonesOfEachProcessor[iProc];
    }

    PH_AllReduce(&nZonesInCurrentProcess, &nZonesGlobal, 1, PH_SUM);
    this->localZoneIndex = new int[nZonesGlobal];
    this->zoneProc       = new int[nZonesGlobal];

    int *tempZoneProc = new int[nZonesGlobal];
    SetField(tempZoneProc, 0, nZonesGlobal);

    for (int iZone = 0; iZone < nZonesInCurrentProcess; ++ iZone)
    {
        int gridIndex = iZone + zoneInsexStart;
        gridContainer[iZone]->SetZoneID(gridIndex);
        gridContainer[iZone]->SetZoneLocalID(iZone);

        tempZoneProc[gridIndex]   = myid;
        localZoneIndex[gridIndex] = iZone;
    }

    PH_AllReduce(tempZoneProc, zoneProc, nZonesGlobal, MPI_SUM);
    delete [] tempZoneProc; tempZoneProc = nullptr;
}

void Mesh_ConstructInterfaceLink::ChangeNoBoundaryCondition()
{
    int ignoreNoBoundaryCondition = GlobalDataBase::GetIntParaFromDB("omit_no_bound_bc");
    if (ignoreNoBoundaryCondition == 1)
    {
        return;
    }

    for (int iZone = 0; iZone < nZonesInCurrentProcess; ++ iZone)
    {
        Grid *grid = gridContainer[iZone];
        grid->ChangeBCType(PHENGLEI::NO_BOUNDARY_CONDITION, PHENGLEI::INTERFACE);
    }
}




}