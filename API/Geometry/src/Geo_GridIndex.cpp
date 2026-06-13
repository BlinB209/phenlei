#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_UnstructGrid.h"
#include "PHMpi.h"
#pragma warning(disable:4100)

using namespace std;

namespace PHSPACE
{
vector< vector< Grid * > * > ggrids;
ZoneConnectivity *zoneConnectivity;
ZoneConnectivityForPoint *zoneConnectivityForPoint;
RDouble globalGridMinBox[3], globalGridMaxBox[3];

Grid * GetGrid(int iZone, int level)
{
    if (!ggrids[iZone])
    {
        return 0;
    }

    return (*ggrids[iZone])[level];
}

void DeleteGGrids()
{
    ggrids.clear();
}

InterfaceInfo * GetInterfaceInfo(int iZone, int level)
{
    Grid *grid = GetGrid(iZone, level);
    if (!grid) return 0;
    return grid->GetInterfaceInfo();
}

InterpointInformation * GetInterpointInfo(int iZone, int level)
{
    Grid *grid = GetGrid(iZone, level);
    if (! grid)
    {
        return 0;
    }
    //! Temporarily debug to obtain the interpoint information of unstructured grid.
    return (static_cast<UnstructGrid *>(grid))->GetInterpointInfo();
}

OversetInfoProxy * GetOversetInfoProxy(int iZone, int level)
{
    Grid *grid = GetGrid(iZone, level);
    if (!grid) return 0;
    return grid->GetOversetInfoProxy();
}

uint_t GetNumberOfMultiGrid(int iZone)
{
    return ggrids[iZone]->size();
}

uint_t GetNumberOfMultiGrid()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (ggrids[iZone])
        {
            return ggrids[iZone]->size();
        }
    }
    return 0;
}

void AddGridToGlobal(int iZone, vector <Grid *> *grids)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    if (ggrids.size() == 0)
    {
        ggrids.resize(nZones, 0);
    }
    ggrids[iZone] = grids;
}

void GetGlobalMinMaxBox(RDouble *globalMinBox, RDouble *globalMaxBox)
{
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        globalMinBox[iDim] = globalGridMinBox[iDim];
        globalMaxBox[iDim] = globalGridMaxBox[iDim];
    }
}

void SetGlobalMinMaxBox(RDouble * globalMinBox, RDouble * globalMaxBox)
{
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        globalGridMinBox[iDim] = globalMinBox[iDim];
        globalGridMaxBox[iDim] = globalMaxBox[iDim];
    }
}

ZoneConnectivityForPoint::ZoneConnectivityForPoint()
{
    numberOfZonesForPoint = 0;
    zoneNeighborForPoint  = 0;
};

ZoneConnectivityForPoint::~ZoneConnectivityForPoint()
{
    for (int iZone = 0; iZone < numberOfZonesForPoint; ++ iZone)
    {
        delete zoneNeighborForPoint[iZone];
    }
    delete [] zoneNeighborForPoint;
}

int  ZoneConnectivityForPoint::GetNumberOfZonesForPoint() const
{ 
    return numberOfZonesForPoint;
};

void ZoneConnectivityForPoint::SetNumberOfZonesForPoint(int numberOfZonesForPoint)
{
    this->numberOfZonesForPoint = numberOfZonesForPoint;
};

void ZoneConnectivityForPoint::CreateAllZones(int numberOfZonesForPoint)
{
    SetNumberOfZonesForPoint(numberOfZonesForPoint);
    zoneNeighborForPoint = new ZoneNeighbor *[numberOfZonesForPoint];
    for (int iZone = 0; iZone < numberOfZonesForPoint; ++ iZone)
    {
        zoneNeighborForPoint[iZone] = new ZoneNeighbor();
    }
}

ZoneNeighbor * ZoneConnectivityForPoint::GetZoneNeighborForPoint(int iZone) const 
{
    return zoneNeighborForPoint[iZone];
}

GridID::GridID(int index, int local_index)
{
    this->index       = index;
    this->local_index = local_index;
}

GridID::~GridID()
{
}

OversetInfo::OversetInfo()
{
    cellIndexForTransfer = 0;
    dataStorageForTransfer = new Data_ParamFieldSuite();
}

OversetInfo::~OversetInfo()
{
    delete dataStorageForTransfer;
    delete [] cellIndexForTransfer;
}

void OversetInfo::AllocateData()
{
    int numberOfCellForTransfer = this->GetNumberOfCellForTransfer();
    int * cellIndexForTransfer = new int [ numberOfCellForTransfer ];
    this->SetCellIndexForTransfer(cellIndexForTransfer);
}

void OversetInfo::Encode()
{
    ;
}

void OversetInfo::Decode(DataContainer *cdata)
{
    int zoneIndexOfNeighbor = 0;
    cdata->Read(reinterpret_cast<char *>(&zoneIndexOfNeighbor), sizeof(int));
    SetZoneIndexOfNeighbor(zoneIndexOfNeighbor);

    int numberOfCellForTransfer;
    cdata->Read(reinterpret_cast<char *>(&numberOfCellForTransfer), sizeof(int));

    int * cellIndexForTransfer = new int [numberOfCellForTransfer];
    cdata->Read(reinterpret_cast<char *>(cellIndexForTransfer), numberOfCellForTransfer * sizeof(int));
    SetCellIndexForTransfer(cellIndexForTransfer);
}

CompositeOversetInfo::CompositeOversetInfo()
{
}

CompositeOversetInfo::~CompositeOversetInfo()
{
}


void CompositeOversetInfo::Encode(DataContainer *cdata)
{
    int numberOfZone = this->GetNumberOfZone();
    cdata->Write(reinterpret_cast<char *>(&numberOfZone), sizeof(int));

    int numberOfCell = this->GetNumberOfCell();
    cdata->Write(reinterpret_cast<char *>(&numberOfCell), sizeof(int));

    for (int iZone = 0; iZone < numberOfZone; ++ iZone)
    {
        OversetInfo *oversetInfo = GetOversetInfo(iZone);
        oversetInfo->Encode();
    }
}

void CompositeOversetInfo::Decode(DataContainer *cdata)
{
    int numberOfZone;
    cdata->Read(reinterpret_cast<char *>(&numberOfZone), sizeof(int));

    SetNumberOfZone(numberOfZone);

    int numberOfCell;
    cdata->Read(reinterpret_cast<char *>(&numberOfCell), sizeof(int));
    SetNumberOfCell(numberOfCell);

    for (int iZone = 0; iZone < numberOfZone; ++ iZone)
    {
        OversetInfo *oversetInfo = new OversetInfo();
        AddOversetInfo(oversetInfo);
        oversetInfo->Decode(cdata);
    }
}

OversetInfoProxy::OversetInfoProxy()
{
    sendIndexProxy = 0;
    recvIndexProxy = 0;
}

OversetInfoProxy::~OversetInfoProxy()
{
}

void OversetInfoProxy::AllocateNeighborInfo()
{
}

void OversetInfoProxy::Encode(DataContainer *cdata)
{
    sendInfo->Encode(cdata);
    recvInfo->Encode(cdata);
}

void OversetInfoProxy::Decode(DataContainer *cdata)
{
    sendInfo->Decode(cdata);
    recvInfo->Decode(cdata);
}

int OversetInfoProxy::GetNIFaceOfNeighbor(int ineighbor) const
{
    return -1;
}

void OversetInfoProxy::InitReceivingInformation(vector<int> &cellIndexRef1, vector<int> &zoneIndexRef2, vector<int> &cellIndexRef2)
{
    this->cellIndexRef1 = cellIndexRef1;
    this->zoneIndexRef2 = zoneIndexRef2;
    this->cellIndexRef2 = cellIndexRef2;

    set<int> neighborZoneIndexSet;
    GetNeighborZoneIndexSet(neighborZoneIndexSet);

    IndexProxy * recvIndexProxy      = GetRecvIndexProxy();
    IndexProxy * temporaryIndexProxy = GetTemporaryIndexProxy();

    for (set<int>::iterator iter = neighborZoneIndexSet.begin(); iter != neighborZoneIndexSet.end(); ++ iter)
    {
        //! jZone is the neighbor, which is the interpolated zone overlapped with iZone.
        int jZone = *iter;

        vector<int> cellIndexI, cellIndexJ;

        GetPartOfArrayByCondition(zoneIndexRef2, jZone, cellIndexRef1, cellIndexI);
        GetPartOfArrayByCondition(zoneIndexRef2, jZone, cellIndexRef2, cellIndexJ);

        IndexSlice *indexSlice = new IndexSlice();

        indexSlice->SetTag(jZone);
        indexSlice->SetIndex(cellIndexI);
        recvIndexProxy->AddIndexSlice(indexSlice);

        IndexSlice *temporaryIndexSlice = new IndexSlice();

        temporaryIndexSlice->SetTag(jZone);
        temporaryIndexSlice->SetIndex(cellIndexJ);
        temporaryIndexProxy->AddIndexSlice(temporaryIndexSlice);
    }
}

void OversetInfoProxy::PrepareSendingInformation(Geo_OversetGrid **oversetGrids)
{
    IndexProxy * temporaryIndexProxy = GetTemporaryIndexProxy();

    uint_t numberOfNeighbor = temporaryIndexProxy->size();
    int zoneIndex = this->GetGrid()->GetZoneID();
    
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        IndexSlice * temporaryIndexSlice = temporaryIndexProxy->GetIndexSlice(iNeighbor);
        int jZone = temporaryIndexSlice->GetTag();
        OversetInfoProxy * jOversetInfoProxy = GetOversetInfoProxy(jZone);
        IndexProxy *sendIndexProxy = jOversetInfoProxy->GetSendIndexProxy();

        IndexSlice * sendIndexSlice = new IndexSlice();
        sendIndexSlice->SetTag(zoneIndex);
        sendIndexSlice->SetIndex(temporaryIndexSlice->GetIndex());
        sendIndexProxy->AddIndexSlice(sendIndexSlice);
    }
}

void OversetInfoProxy::PostProcessRecv()
{
    IndexProxy * recvIndexProxy = GetRecvIndexProxy();
    int numberOfCellForRecv = 0;

    uint_t numberOfNeighbor = recvIndexProxy->size();
    
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        int iLength = recvIndexProxy->GetLength(iNeighbor);
        recvIndexProxy->SetIndexMap(iNeighbor, iLength);
        numberOfCellForRecv += iLength;
    }

    this->SetNumberOfCellForRecv(numberOfCellForRecv);
}

void OversetInfoProxy::PostProcessSend()
{
    IndexProxy * sendIndexProxy = GetSendIndexProxy();

    int numberOfCellForSend = 0;
    uint_t numberOfNeighbor = sendIndexProxy->size();
    
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        int iLength = sendIndexProxy->GetLength(iNeighbor);
        sendIndexProxy->SetIndexMap(iNeighbor, iLength);
        numberOfCellForSend += iLength;
    }

    this->SetNumberOfCellForSend(numberOfCellForSend);
}

void OversetInfoProxy::PostProcess(Geo_OversetGrid **oversetGrids)
{
    PostProcessRecv();
    PostProcessSend();
}

void OversetInfoProxy::GetNeighborZoneIndexSet(set<int> &neighborZoneIndexSet)
{
    Array2Set(zoneIndexRef2, neighborZoneIndexSet);
}

void OversetInfoProxy::Process()
{
    set<int> neighborZoneIndexSet;
    GetNeighborZoneIndexSet(neighborZoneIndexSet);

    IndexProxy * recvIndexProxy      = GetRecvIndexProxy();
    IndexProxy * temporaryIndexProxy = GetTemporaryIndexProxy();

    vector<int> & cellIndexRef1 = GetCellIndexRef1();
    vector<int> & cellIndexRef2 = GetCellIndexRef2();
    vector<int> & zoneIndexRef2 = GetZoneIndexRef2();

    for (set<int>::iterator iter = neighborZoneIndexSet.begin(); iter != neighborZoneIndexSet.end(); ++ iter)
    {
        //! jZone is the neighbor,which is the interpolated zone overlapped with iZone.
        int jZone = *iter;

        vector<int> cellIndexI, cellIndexJ;

        GetPartOfArrayByCondition(zoneIndexRef2, jZone, cellIndexRef1, cellIndexI);
        GetPartOfArrayByCondition(zoneIndexRef2, jZone, cellIndexRef2, cellIndexJ);

        IndexSlice * indexSlice = new IndexSlice();

        indexSlice->SetTag(jZone);
        indexSlice->SetIndex(cellIndexI);
        recvIndexProxy->AddIndexSlice(indexSlice);

        IndexSlice * temporaryIndexSlice = new IndexSlice();

        temporaryIndexSlice->SetTag(jZone);
        temporaryIndexSlice->SetIndex(cellIndexJ);
        temporaryIndexProxy->AddIndexSlice(temporaryIndexSlice);
    }
}

}
