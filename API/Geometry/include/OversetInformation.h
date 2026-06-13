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
//! @file      OversetInformation.h
//! @brief     It define overset information reconstruction(unstructured grid).
//! @author    He Xin, Chang Xinghua .

#pragma once

#include "Geo_SimpleGrid.h"
namespace PHSPACE
{

class Grid;
class DataContainer;
class Geo_OversetGrid;
class OversetCell;
class OversetInterplateCell;
class OversetDonorCell;
class OversetInterplateCellManager;
class OversetDonorCellManager;
class Data_ParamFieldSuite;

class OversetSwapManager
{
public:
    OversetSwapManager ();
    ~OversetSwapManager();
protected:
    OversetInterplateCellManager *oversetInterplateCellManager;
    OversetDonorCellManager      *oversetDonorCellManager;
public:
    OversetInterplateCellManager * GetOversetInterplateCellManager() { return oversetInterplateCellManager; }
    OversetDonorCellManager * GetOversetDonorCellManager() { return oversetDonorCellManager; }

    void AddOversetInterplateCell(OversetInterplateCell *oversetInterplateCell);
    void AddOversetDonorCell(OversetDonorCell *oversetDonorCell);
public:
    void PostProcessReceivingOversetInformation();
    void PostProcessSendingOversetInformation();

    void EncodeOversetInformation(DataContainer *dataContainer);
    void DecodeOversetInformation(DataContainer *dataContainer);
protected:
    void EncodeOversetInterplateCellInformation(DataContainer *dataContainer);
    void EncodeOversetDonorCellInformation(DataContainer *dataContainer);

    void DecodeOversetInterplateCellInformation(DataContainer *dataContainer);
    void DecodeOversetDonorCellInformation(DataContainer *dataContainer);
};

class OversetDataProxy
{
public:
    OversetDataProxy();
    ~OversetDataProxy();
protected:
    int sendReceiveFlag; //! 0: send ; 1: receive;
    OversetSwapManager *oversetSwapManager;
public:
    void  Initialize(OversetSwapManager *oversetSwapManager, int sendReceiveFlag);
    int   GetNumberOfNeighbors();
    int   GetNumberOfCells();
    int   GetNumberOfCells(int neighborIndex);
    int * GetCellIndexContainer(int neighborIndex);
    int * GetCellStorageIndexMapping(int neighborIndex);

    RDouble * GetInterCellCenterX(int neighborIndex);
    RDouble * GetInterCellCenterY(int neighborIndex);
    RDouble * GetInterCellCenterZ(int neighborIndex);
    RDouble * GetInterCellVolume (int neighborIndex);
    Data_ParamFieldSuite * GetDataStorage();
public:
    int FindNeighborIndex(int zoneIndex);
    PHVectorInt1D GetNeighborZoneIndexArray();
protected:
    int FindNeighborIndex(PHVectorInt1D &neighborZoneIndexArray, int neighborZoneIndexForSearching);
};

class OversetInformationProxy
{
public:
    OversetInformationProxy ();
    ~OversetInformationProxy();
protected:
    Grid *grid;

    OversetSwapManager *oversetSwapManager;

    PHVector1D < OversetCell * > *oversetCellContainer;

    OversetDataProxy *oversetDataProxyForSend;
    OversetDataProxy *oversetDataProxyForReceive;
public:
    int GetZoneIndex();
    OversetSwapManager * GetOversetSwapManager() { return oversetSwapManager; }
    PHVector1D < OversetCell * > * GetOversetCellContainer() { return oversetCellContainer; }
    void SetOversetCellContainer(PHVector1D < OversetCell * > *oversetCellContainer) { this->oversetCellContainer = oversetCellContainer; }

    PHVectorInt1D GetNeighborZoneIndexArray();
public:
    void GenerateNeighborZoneIndexSet(set < int > &neighborZoneIndexSet);
    void GetNeighborZoneIndexContainer(PHVectorInt1D &neighborZoneIndexContainer);
    void FreeOversetCellContainer(PHVector1D < OversetCell * > *oversetCellContainer);
public:
    OversetDonorCellManager * GetOversetDonorCellManager() { return oversetSwapManager->GetOversetDonorCellManager(); }
    OversetInterplateCellManager * GetOversetInterplateCellManager() { return oversetSwapManager->GetOversetInterplateCellManager(); }

    OversetDataProxy * GetOversetDataProxyForSend() { return oversetDataProxyForSend; }
    OversetDataProxy * GetOversetDataProxyForReceive() { return oversetDataProxyForReceive; }
public:
    void PostProcessReceivingOversetInformation();
    void PostProcessSendingOversetInformation();
public:

    void SetGrid(Grid *grid) { this->grid = grid; }
    Grid * GetGrid() { return grid; }

    void PostProcessOversetInformation();
    void InitializeReceivingInformation(PHVector1D < OversetCell * > *oversetCellContainer);

public:
    void GetCellIndexContainer(int donorZoneIndexIn, PHVectorInt1D &cellIndexContainer);
    void GetDonorCellIndexContainer(int donorZoneIndexIn, PHVectorInt1D &donorCellIndexContainer,
                                                          PHVectorRDouble1D &interCellCenterX,
                                                          PHVectorRDouble1D &interCellCenterY,
                                                          PHVectorRDouble1D &interCellCenterZ,
                                                          PHVectorRDouble1D &interCellVolume);
public:
    void EncodeOversetInformation(DataContainer *dataContainer);
    void DecodeOversetInformation(DataContainer *dataContainer);
};

class ActionKey;
class OversetInformationProxy;

OversetInformationProxy * GetOversetInformationProxy(int iZone, int gridLevelIndex = 0);

class OversetTopologyManager
{
public:
    OversetTopologyManager();
    ~OversetTopologyManager();
protected:
    PHVectorInt2D oversetTopology;
    PHVectorInt2D oversetTopologyForSend;
    PHVectorInt2D oversetTopologyForReceive;
public:
    PHVectorInt2D & GetOversetTopology();
    PHVectorInt2D & GetOversetTopologyForSend();
    PHVectorInt2D & GetOversetTopologyForReceive();

    uint_t GetNumberOfOversetNeighbors(int zoneIndex);
    int GetNeighborOversetZoneIndex(int zoneIndex, int ithOversetNeighbor);
public:
    void InitializeOversetInterfaceTopology();
    void InitializeOversetInterfaceTopologyServerCollectAndBcast();
};

OversetTopologyManager * GetOversetTopologyManager();
void InitializeOversetInterfaceTopology();

class OversetCell
{
public:
    OversetCell();
    ~OversetCell();
protected:
    int zoneIndex;    //!zone index to be interpolated
    int cellIndex;    //!cell index to be interpolated
    int donorZoneIndex;
    int donorCellIndex;
public:
    void SetZoneIndex     (int zoneIndex     ) { this->zoneIndex      = zoneIndex     ; }
    void SetCellIndex     (int cellIndex     ) { this->cellIndex      = cellIndex     ; }
    void SetDonorZoneIndex(int donorZoneIndex) { this->donorZoneIndex = donorZoneIndex; }
    void SetDonorCellIndex(int donorCellIndex) { this->donorCellIndex = donorCellIndex; }

    int GetZoneIndex() { return zoneIndex; }
    int GetCellIndex() { return cellIndex; }
    int GetDonorZoneIndex() { return donorZoneIndex; }
    int GetDonorCellIndex() { return donorCellIndex; }
};

class OversetInterplateCell
{
public:
    OversetInterplateCell();
    ~OversetInterplateCell();
protected:
    int zoneIndex;
    int donorZoneIndex;
    PHVectorInt1D cellIndexContainer;    //!the set of elements to be interpolated
    PHVectorInt1D storageIndexContainer;
    Data_ParamFieldSuite *dataStorage;
public:
    Data_ParamFieldSuite * GetDataStorage() { return dataStorage; }

    void SetZoneIndex     (int zoneIndex     ) { this->zoneIndex      = zoneIndex; }
    void SetDonorZoneIndex(int donorZoneIndex) { this->donorZoneIndex = donorZoneIndex; }

    int GetZoneIndex     () { return zoneIndex; }
    int GetDonorZoneIndex() { return donorZoneIndex; }

    PHVectorInt1D & GetCellIndexContainer() { return cellIndexContainer; }

    int GetNumberOfOversetInterplateCells() { return static_cast<int>(cellIndexContainer.size()); }

    PHVectorInt1D & GetStorageIndexContainer() { return storageIndexContainer; }
public:
    void EncodeOversetInterplateCellInformation(DataContainer *dataContainer);
    void DecodeOversetInterplateCellInformation(DataContainer *dataContainer);
};

class OversetDonorCell
{
public:
    OversetDonorCell ();
    ~OversetDonorCell();
protected:
    int zoneIndex;
    int interplateZoneIndex;
    PHVectorInt1D donorCellIndexContainer;
    PHVectorInt1D storageIndexContainer;
    Data_ParamFieldSuite *dataStorage;
    PHVectorRDouble1D interCellCenterX;
    PHVectorRDouble1D interCellCenterY;
    PHVectorRDouble1D interCellCenterZ;
    PHVectorRDouble1D interCellVolume;
public:
    Data_ParamFieldSuite * GetDataStorage() { return dataStorage; }

    void SetZoneIndex(int zoneIndex) { this->zoneIndex = zoneIndex; }
    int  GetZoneIndex() { return zoneIndex; }

    void SetInterplateZoneIndex(int interplateZoneIndex) { this->interplateZoneIndex = interplateZoneIndex; }
    int  GetInterplateZoneIndex() { return interplateZoneIndex; }

    PHVectorInt1D     & GetDonorCellIndexContainer() { return donorCellIndexContainer; }
    PHVectorRDouble1D & GetInterCellCenterX() { return interCellCenterX; }
    PHVectorRDouble1D & GetInterCellCenterY() { return interCellCenterY; }
    PHVectorRDouble1D & GetInterCellCenterZ() { return interCellCenterZ; }
    PHVectorRDouble1D & GetInterCellVolume()  { return interCellVolume; }

    void SetDonorCellIndexContainer(PHVectorInt1D     & donorCellIndexContainerIn) { this->donorCellIndexContainer = donorCellIndexContainerIn; }
    void SetInterCellCenterX       (PHVectorRDouble1D & interCellCenterXIn       ) { this->interCellCenterX        = interCellCenterXIn; }
    void SetInterCellCenterY       (PHVectorRDouble1D & interCellCenterYIn       ) { this->interCellCenterY        = interCellCenterYIn; }
    void SetInterCellCenterZ       (PHVectorRDouble1D & interCellCenterZIn       ) { this->interCellCenterZ        = interCellCenterZIn; }
    void SetInterCellVolume        (PHVectorRDouble1D & interCellVolumeIn       ) { this->interCellVolume        = interCellVolumeIn; }

    int GetNumberOfOversetDonorCells() { return static_cast<int>(donorCellIndexContainer.size()); }

    PHVectorInt1D & GetStorageIndexContainer() { return storageIndexContainer; }
public:
    void EncodeOversetDonorCellInformation(DataContainer *dataContainer);
    void DecodeOversetDonorCellInformation(DataContainer *dataContainer);
};

class OversetInterplateCellManager
{
public:
    OversetInterplateCellManager ();
    ~OversetInterplateCellManager();
protected:
    PHVector1D < OversetInterplateCell * > oversetInterplateCellContainer;
    int totalNumberOfOversetInterplateCells;
    Data_ParamFieldSuite *dataStorage;
public:
    PHVector1D < OversetInterplateCell * > & GetOversetInterplateCellContainer() { return oversetInterplateCellContainer; }

    int GetNumberOfNeighbors() { return static_cast<int>(oversetInterplateCellContainer.size()); }

    int  GetTotalNumberOfOversetInterplateCells() { return totalNumberOfOversetInterplateCells; }
    void SetTotalNumberOfOversetInterplateCells(int totalNumberOfOversetInterplateCells) { this->totalNumberOfOversetInterplateCells = totalNumberOfOversetInterplateCells; }

    int GetNumberOfOversetInterplateCells(int iNeighbor);

    int * GetOversetInterplateCell(int iNeighbor);
    int * GetInterplateCellStorageIndex(int iNeighbor);

    Data_ParamFieldSuite * GetDataStorage() { return dataStorage; }

    PHVectorInt1D GetNeighborDonorZoneIndexArray();
public:
    void PostProcessReceivingOversetInformation();
    void EncodeOversetInterplateCellInformation(DataContainer *dataContainer);
    void DecodeOversetInterplateCellInformation(DataContainer *dataContainer);
public:
    void AddOversetInterplateCell(OversetInterplateCell *oversetInterplateCell);
};

class OversetDonorCellManager
{
public:
    OversetDonorCellManager();
    ~OversetDonorCellManager();
protected:
    PHVector1D < OversetDonorCell * > oversetDonorCellContainer;
    int totalNumberOfOversetDonorCells;
    Data_ParamFieldSuite *dataStorage;
public:
    PHVector1D < OversetDonorCell * > & GetOversetDonorCellContainer() { return oversetDonorCellContainer; }

    int GetNumberOfNeighbors() { return static_cast<int>(oversetDonorCellContainer.size()); }

    int  GetTotalNumberOfOversetDonorCells() { return totalNumberOfOversetDonorCells; }
    void SetTotalNumberOfOversetDonorCells(int totalNumberOfOversetDonorCells) { this->totalNumberOfOversetDonorCells = totalNumberOfOversetDonorCells; }

    int GetNumberOfOversetDonorCells(int iNeighbor);

    int * GetOversetDonorCell(int iNeighbor);
    int * GetDonorCellStorageIndex(int iNeighbor);

    RDouble * GetInterCellCenterX(int iNeighbor);
    RDouble * GetInterCellCenterY(int iNeighbor);
    RDouble * GetInterCellCenterZ(int iNeighbor);
    RDouble * GetInterCellVolume(int iNeighbor);

    Data_ParamFieldSuite * GetDataStorage() { return dataStorage; }

    PHVectorInt1D GetNeighborInterplateZoneIndexArray();
public:
    void PostProcessSendingOversetInformation();
    void EncodeOversetDonorCellInformation(DataContainer *dataContainer);
    void DecodeOversetDonorCellInformation(DataContainer *dataContainer);
public:
    void AddOversetDonorCell(OversetDonorCell *oversetDonorCell);
};

uint_t GetNumberOfOversetNeighbors(int zoneIndex);
int GetNeighborOversetZoneIndex(int zoneIndex, int iNeighbor);

}
