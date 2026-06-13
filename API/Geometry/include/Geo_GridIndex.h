//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Geo_GridIndex.h
//! @brief     This file is part of PHengLEI Technology software library.
//!            It defines the class 'Geo_GridIndex', which is the index information used in Grid.
//! @author    He Xin, Bell.

#pragma once
#include "Data_ParamFieldSuite.h"
#include "Geo_SimpleGrid.h"
#include "GlobalDataBase.h"
#include "DataContainer.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"

namespace PHSPACE
{
class Grid;
class ActionKey;
class ActionTag;
class LinkStruct;
class Geo_OversetGrid;

class IndexSlice
{
private:
    int tag;
    int nlen;
    vector<int> index;
public:
    int  GetTag() { return tag; }
    void SetTag(int tag) { this->tag = tag; }

    int  GetLength() { return nlen; }
    void SetLength(int nlen) { this->nlen = nlen; }

    uint_t size() { return index.size(); }

    int * GetIndexPointer() { return &index[0]; }

    vector<int> & GetIndex() { return index; }
    void SetIndex(vector<int> &index) { this->index = index; }
};

class IndexProxy
{
private:
    vector <IndexSlice *> indexArray;
    vector <IndexSlice *> indexMapArray;
    vector <int> zoneIndexArray;
    int init;
private:
    void InitZoneIndexArray()
    {
        if (init) return;
        for (std::size_t i = 0; i < indexArray.size(); ++ i)
        {
            zoneIndexArray.push_back(indexArray[i]->GetTag());
        }
        init = 1;
    }
public:
    IndexProxy()
    {
        init = 0;
    }

    ~IndexProxy()
    {
        ;
    }
public:
    IndexSlice * GetIndexSlice(int m) { return indexArray[m]; }
    void  AddIndexSlice(IndexSlice *indexSlice) { indexArray.push_back(indexSlice); }
    int   GetLength(int m) { return indexArray[m]->GetLength(); }
    int * GetIndexPointer(int m) { return indexArray[m]->GetIndexPointer(); }
    int * GetIndexMapPointer(int m) { return indexMapArray[m]->GetIndexPointer(); }

    void SetIndexMap(int m, int start)
    {
        int *indexMap = GetIndexMapPointer(m);
        int length = GetLength(m);
        for (int i = 0; i < length; ++ i)
        {
            indexMap[i] = start + i;
        }
    }

    uint_t size() const { return indexArray.size(); }

    vector < int > & GetZoneIndexArray()
    {
        InitZoneIndexArray();
        return zoneIndexArray; 
    }
};

class OversetInfo
{
private:
    int zoneid, zoneIndexOfNeighbor;
    int numberOfCellForTransfer;
    int *cellIndexForTransfer;
private:
    Data_ParamFieldSuite *dataStorageForTransfer;
public:
    OversetInfo();
    ~OversetInfo();
public:
    int * GetCellIndexForTransfer() const { return cellIndexForTransfer; }
    void SetCellIndexForTransfer(int *cellIndexForTransfer) { this->cellIndexForTransfer = cellIndexForTransfer; }

    int  GetNumberOfCellForTransfer() const { return numberOfCellForTransfer; }
    void SetNumberOfCellForTransfer(int numberOfCellForTransfer) { this->numberOfCellForTransfer = numberOfCellForTransfer; }

    int  GetZoneIndexOfNeighbor() const { return zoneIndexOfNeighbor; }
    void SetZoneIndexOfNeighbor(int zoneIndexOfNeighbor) { this->zoneIndexOfNeighbor = zoneIndexOfNeighbor; }
public:
    void AllocateData();
public:
    void Encode();
    void Decode(DataContainer *cdata);
};

class CompositeOversetInfo
{
private:
    int numberOfZone;
    int numberOfCell;
    int *cellIndex;
    vector <OversetInfo *> oversetInfoArray;
    Data_ParamFieldSuite *dataStorage;
public:
    CompositeOversetInfo();
    ~CompositeOversetInfo();
public:
    int  GetNumberOfZone() const { return numberOfZone; }
    void SetNumberOfZone(int numberOfZone) { this->numberOfZone = numberOfZone; }

    int  GetNumberOfCell() const { return numberOfCell; }
    void SetNumberOfCell(int numberOfCell) { this->numberOfCell = numberOfCell; }

    OversetInfo * GetOversetInfo(int ineighbor) { return oversetInfoArray[ineighbor]; }

    int GetNumberOfCell(int ineighbor)
    { 
        OversetInfo *oversetInfo = GetOversetInfo(ineighbor);
        return oversetInfo->GetNumberOfCellForTransfer();
    };

    int * GetCellIndex(int ineighbor)
    {
        OversetInfo *oversetInfo = GetOversetInfo(ineighbor);
        return oversetInfo->GetCellIndexForTransfer();
    };

    int FindIthNeighbor() { return -1; }

    Data_ParamFieldSuite * GetDataStorage() const { return dataStorage; }

    void AddOversetInfo(OversetInfo *oversetInfo) { oversetInfoArray.push_back(oversetInfo); }
public:
    void Encode(DataContainer *cdata);
    void Decode(DataContainer *cdata);
};

class OversetInfoProxy
{
private:
    CompositeOversetInfo *sendInfo;
    CompositeOversetInfo *recvInfo;
    Grid *grid;
private:
    int zoneIndex;

    IndexProxy *sendIndexProxy;
    IndexProxy *recvIndexProxy;

    IndexProxy *temporaryIndexProxy;

    vector<int> cellIndexRef1, cellIndexRef2, zoneIndexRef2;

    int numberOfCellForRecv, numberOfCellForSend;
public:
    vector<int> & GetCellIndexRef1() { return cellIndexRef1; }
    vector<int> & GetCellIndexRef2() { return cellIndexRef2; }
    vector<int> & GetZoneIndexRef2() { return zoneIndexRef2; }

    void Process();
    void GetNeighborZoneIndexSet(set<int> &neighborZoneIndexSet);
public:
    IndexProxy * GetSendIndexProxy() { return sendIndexProxy; }
    IndexProxy * GetRecvIndexProxy() { return recvIndexProxy; }

    IndexProxy * GetTemporaryIndexProxy() { return temporaryIndexProxy; }
public:
    vector<int> & GetNeighborZoneIndexArray() { return sendIndexProxy->GetZoneIndexArray(); }    // for send
public:
    OversetInfoProxy();
    ~OversetInfoProxy();
public:
    uint_t GetNumberOfNeighborForSend() { return sendIndexProxy->size(); }
    uint_t GetNumberOfNeighborForRecv() { return recvIndexProxy->size(); }

    int GetNumberOfCellForRecv() { return numberOfCellForRecv; }
    int GetNumberOfCellForSend() { return numberOfCellForSend; }

    void SetNumberOfCellForRecv(int numberOfCellForRecv) { this->numberOfCellForRecv = numberOfCellForRecv; }
    void SetNumberOfCellForSend(int numberOfCellForSend) { this->numberOfCellForSend = numberOfCellForSend; }

    int GetNumberOfCellForRecv(int iNeighbor) { return recvIndexProxy->GetLength(iNeighbor); }
    int GetNumberOfCellForSend(int iNeighbor) { return sendIndexProxy->GetLength(iNeighbor); }

    int * GetCellIndexForRecv(int iNeighbor) { return recvIndexProxy->GetIndexPointer(iNeighbor); }
    int * GetCellIndexForSend(int iNeighbor) { return sendIndexProxy->GetIndexPointer(iNeighbor); }

    int * GetCellIndexMapForRecv(int iNeighbor) { return recvIndexProxy->GetIndexMapPointer(iNeighbor); }
    int * GetCellIndexMapForSend(int iNeighbor) { return sendIndexProxy->GetIndexMapPointer(iNeighbor); }

    Data_ParamFieldSuite * GetRecvDataStorage() const { return recvInfo->GetDataStorage(); }
    Data_ParamFieldSuite * GetSendDataStorage() const { return sendInfo->GetDataStorage(); }

    void PostProcessRecv();
    void PostProcessSend();

    int FindIthNeighbor(vector<int> &zoneIndexArray, int zone)
    {
        for (std::size_t i = 0; i < zoneIndexArray.size(); ++ i)
        {
            if (zoneIndexArray[i] == zone)
            {
                return int(i);
            }
        }
        return -1; 
    };

    int FindIthNeighborOfRecv(int zone)
    {
        vector<int> &zoneIndexArray = recvIndexProxy->GetZoneIndexArray();
        return FindIthNeighbor(zoneIndexArray, zone);
    };

    int FindIthNeighborOfSend(int zone)
    {
        vector<int> &zoneIndexArray = sendIndexProxy->GetZoneIndexArray();
        return FindIthNeighbor(zoneIndexArray, zone);

    };
public:
    void SetGrid(Grid *grid) { this->grid = grid; }
    Grid * GetGrid() { return grid; }
    void PrepareSendingInformation(Geo_OversetGrid **oversetGrids);
    void PostProcess(Geo_OversetGrid **oversetGrids);
    void InitReceivingInformation(vector<int> &cellIndexRef1, vector<int> &zoneIndexRef2, vector<int> &cellIndexRef2);
    CompositeOversetInfo * GetRecvInfo() const { return recvInfo; }
    CompositeOversetInfo * GetSendInfo() const { return sendInfo; }

    int  GetNIFaceOfNeighbor(int ineighbor) const;

    void AllocateNeighborInfo();
public:
    void Encode(DataContainer *cdata);
    void Decode(DataContainer *cdata);
};

class ZoneNeighbor
{
public:
    ZoneNeighbor()
    {
    };

    ~ZoneNeighbor()
    {
    };
public:
    void SetNumberOfNeighbor(int neighbor) { data.resize(neighbor); }
    int  GetNumberOfNeighbor() const { return static_cast<int>(data.size()); }
    void AddZoneIndexOfNeighbor(int value) { data.push_back(value); }
    int * GetZoneIndexOfNeighbor() { return &data[0]; }
    int  GetZoneIndexOfNeighbor(const uint_t &ineighbor) const { return data[ineighbor]; }
    bool hasNeighborZone(int zoneIndex)
    {
        return find(data.begin(), data.end(), zoneIndex) == data.end() ? false : true;
    }
private:
    vector<int> data;
};

class ZoneConnectivity
{
public:
    ZoneConnectivity()
    { 
        numberOfZones = 0;
        zoneNeighbor  = 0;
    };
    ~ZoneConnectivity()
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            delete zoneNeighbor[iZone];
        }
        delete [] zoneNeighbor;
    }
    int  GetNumberOfZones() const { return numberOfZones; }
    void SetNumberOfZones(int nZones) { this->numberOfZones = nZones; }
    void CreateAllZones(int nZones)
    {
        SetNumberOfZones(nZones);
        zoneNeighbor = new ZoneNeighbor *[nZones];
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            zoneNeighbor[iZone] = new ZoneNeighbor();
        } 
    }
    ZoneNeighbor * GetZoneNeighbor(int iZone) const { return zoneNeighbor[iZone]; }
private:
    int numberOfZones;
    ZoneNeighbor **zoneNeighbor;
};

//! @brief the zone connectivity for point.\n
class ZoneConnectivityForPoint
{
public:
    ZoneConnectivityForPoint();
    ~ZoneConnectivityForPoint();
    int  GetNumberOfZonesForPoint() const;
    void SetNumberOfZonesForPoint(int numberOfZonesForPoint);
    void CreateAllZones(int numberOfZonesForPoint);
    ZoneNeighbor * GetZoneNeighborForPoint(int iZone) const;
private:
    int numberOfZonesForPoint;
    ZoneNeighbor **zoneNeighborForPoint;
};

class GridID
{
public:
    GridID(int index = 0, int local_index = 0);
    ~GridID();
public:
    void SetIndex(int index) { this->index = index; }
    void SetLocalIndex(int local_index) { this->local_index = local_index; }
    int GetIndex() const { return index; }
    int GetLocalIndex() const { return local_index; }
private:
    int index;
    int local_index;
};

//! The level  number of multigrid.
uint_t GetNumberOfMultiGrid();
uint_t GetNumberOfMultiGrid(int iZone);

Grid * GetGrid(int iZone, int level = 0);
void DeleteGGrids();
void AddGridToGlobal(int iZone, vector <Grid *> *grids);
void GetGlobalMinMaxBox(RDouble *globalMinBox, RDouble *globalMaxBox);
void SetGlobalMinMaxBox(RDouble *globalMinBox, RDouble *globalMaxBox);

InterfaceInfo * GetInterfaceInfo(int iZone, int level = 0);
//! Return the interpoint information of the iZone.
//! @param[in] iZone : the zone number.
//! @param[in] level : the number of ghost cell layers.
//! @param[out] the object of the class InterpointInformation.
InterpointInformation * GetInterpointInfo(int iZone, int level);
OversetInfoProxy * GetOversetInfoProxy(int iZone, int level = 0);

}
