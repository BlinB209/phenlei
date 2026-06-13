#include "PHMpi.h"
#include "OversetInformation.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "Geo_UnstructBC.h"
#include "Geo_UnstructGrid.h"
#include "Geo_StructGrid.h"
#include "Geo_FaceTopo_Unstruct.h"

using namespace std;

namespace PHSPACE
{

OversetInformationProxy *GetOversetInformationProxy(int iZone, int gridLevelIndex)
{
    Grid *grid = GetGrid(iZone, gridLevelIndex);
    if (! grid) return 0;
    UnstructGrid *gridUnstr = UnstructGridCast(grid);
    return gridUnstr->GetOversetInformationProxy();
}

OversetSwapManager::OversetSwapManager()
{
    oversetInterplateCellManager = new OversetInterplateCellManager();
    oversetDonorCellManager = new OversetDonorCellManager();
}

OversetSwapManager::~OversetSwapManager()
{
    delete oversetInterplateCellManager;
    delete oversetDonorCellManager;
};

void OversetSwapManager::AddOversetInterplateCell(OversetInterplateCell *oversetInterplateCell)
{
    OversetInterplateCellManager *oversetInterplateCellManager = this->GetOversetInterplateCellManager();
    oversetInterplateCellManager->AddOversetInterplateCell(oversetInterplateCell);
}

void OversetSwapManager::AddOversetDonorCell(OversetDonorCell *oversetDonorCell)
{
    OversetDonorCellManager *oversetDonorCellManager = this->GetOversetDonorCellManager();
    oversetDonorCellManager->AddOversetDonorCell(oversetDonorCell);
}

void OversetSwapManager::PostProcessReceivingOversetInformation()
{
    OversetInterplateCellManager *oversetInterplateCellManager = this->GetOversetInterplateCellManager();
    oversetInterplateCellManager->PostProcessReceivingOversetInformation();
}

void OversetSwapManager::PostProcessSendingOversetInformation()
{
    OversetDonorCellManager *oversetDonorCellManager = this->GetOversetDonorCellManager();
    oversetDonorCellManager->PostProcessSendingOversetInformation();
}

void OversetSwapManager::EncodeOversetInterplateCellInformation(DataContainer *dataContainer)
{
    OversetInterplateCellManager *oversetInterplateCellManager = this->GetOversetInterplateCellManager();
    oversetInterplateCellManager->EncodeOversetInterplateCellInformation(dataContainer);
}

void OversetSwapManager::EncodeOversetDonorCellInformation(DataContainer *dataContainer)
{
    OversetDonorCellManager *oversetDonorCellManager = this->GetOversetDonorCellManager();
    oversetDonorCellManager->EncodeOversetDonorCellInformation(dataContainer);
}

void OversetSwapManager::DecodeOversetInterplateCellInformation(DataContainer *dataContainer)
{
    OversetInterplateCellManager *oversetInterplateCellManager = this->GetOversetInterplateCellManager();
    oversetInterplateCellManager->DecodeOversetInterplateCellInformation(dataContainer);
}

void OversetSwapManager::DecodeOversetDonorCellInformation(DataContainer *dataContainer)
{
    OversetDonorCellManager *oversetDonorCellManager = this->GetOversetDonorCellManager();
    oversetDonorCellManager->DecodeOversetDonorCellInformation(dataContainer);
}

void OversetSwapManager::EncodeOversetInformation(DataContainer *dataContainer)
{
    this->EncodeOversetInterplateCellInformation(dataContainer);
    this->EncodeOversetDonorCellInformation(dataContainer);
}

void OversetSwapManager::DecodeOversetInformation(DataContainer *dataContainer)
{
    this->DecodeOversetInterplateCellInformation(dataContainer);
    this->DecodeOversetDonorCellInformation(dataContainer);
}

OversetDataProxy::OversetDataProxy()
{
    ;
}

OversetDataProxy::~OversetDataProxy()
{
    ;
}

void OversetDataProxy::Initialize(OversetSwapManager *oversetSwapManager, int sendReceiveFlag)
{
    this->oversetSwapManager = oversetSwapManager;
    this->sendReceiveFlag = sendReceiveFlag;
}

int OversetDataProxy::GetNumberOfNeighbors()
{
    if (sendReceiveFlag == 0)     //! Send.
    {
        return oversetSwapManager->GetOversetDonorCellManager()->GetNumberOfNeighbors();
    }
    else     //! Receive.
    {
        return oversetSwapManager->GetOversetInterplateCellManager()->GetNumberOfNeighbors();
    }
};

int OversetDataProxy::GetNumberOfCells()
{
    if (sendReceiveFlag == 0)     //! Send.
    {
        return oversetSwapManager->GetOversetDonorCellManager()->GetTotalNumberOfOversetDonorCells();
    }
    else     //! Receive.
    {
        return oversetSwapManager->GetOversetInterplateCellManager()->GetTotalNumberOfOversetInterplateCells();
    }
}

int OversetDataProxy::GetNumberOfCells(int neighborIndex)
{
    if (sendReceiveFlag == 0)     //! Send.
    {
        return oversetSwapManager->GetOversetDonorCellManager()->GetNumberOfOversetDonorCells(neighborIndex);
    }
    else      //! Receive.
    {
        return oversetSwapManager->GetOversetInterplateCellManager()->GetNumberOfOversetInterplateCells(neighborIndex);
    }
};

RDouble *OversetDataProxy::GetInterCellCenterX(int neighborIndex)
{
    return oversetSwapManager->GetOversetDonorCellManager()->GetInterCellCenterX(neighborIndex);
}

RDouble *OversetDataProxy::GetInterCellCenterY(int neighborIndex)
{
    return oversetSwapManager->GetOversetDonorCellManager()->GetInterCellCenterY(neighborIndex);
}

RDouble *OversetDataProxy::GetInterCellCenterZ(int neighborIndex)
{
    return oversetSwapManager->GetOversetDonorCellManager()->GetInterCellCenterZ(neighborIndex);
}
RDouble *OversetDataProxy:: GetInterCellVolume(int neighborIndex)
{
    return oversetSwapManager->GetOversetDonorCellManager()->GetInterCellVolume(neighborIndex);
}
int *OversetDataProxy::GetCellIndexContainer(int neighborIndex)
{
    if (sendReceiveFlag == 0)     //! Send.
    {
        return oversetSwapManager->GetOversetDonorCellManager()->GetOversetDonorCell(neighborIndex);
    }
    else      //! Receive.
    {
        return oversetSwapManager->GetOversetInterplateCellManager()->GetOversetInterplateCell(neighborIndex);
    }
};

int *OversetDataProxy::GetCellStorageIndexMapping(int neighborIndex)
{
    if (sendReceiveFlag == 0)      //! Send.
    {
        return oversetSwapManager->GetOversetDonorCellManager()->GetDonorCellStorageIndex(neighborIndex);
    }
    else      //! Receive.
    {
        return oversetSwapManager->GetOversetInterplateCellManager()->GetInterplateCellStorageIndex(neighborIndex);
    }
};

Data_ParamFieldSuite *OversetDataProxy::GetDataStorage()
{
    if (sendReceiveFlag == 0)     //! Send.
    {
        return oversetSwapManager->GetOversetDonorCellManager()->GetDataStorage();
    }
    else      //! Receive.
    {
        return oversetSwapManager->GetOversetInterplateCellManager()->GetDataStorage();
    }
};

PHVectorInt1D OversetDataProxy::GetNeighborZoneIndexArray()
{
    if (sendReceiveFlag == 0)
    {
        PHVectorInt1D zoneIndexArray = this->oversetSwapManager->GetOversetDonorCellManager()->GetNeighborInterplateZoneIndexArray();
        return zoneIndexArray;
    }
    else
    {
        PHVectorInt1D zoneIndexArray = this->oversetSwapManager->GetOversetInterplateCellManager()->GetNeighborDonorZoneIndexArray();
        return zoneIndexArray;
    }
}

int OversetDataProxy::FindNeighborIndex(int zoneIndex)
{
    PHVectorInt1D zoneIndexArray = GetNeighborZoneIndexArray();

    return this->FindNeighborIndex(zoneIndexArray, zoneIndex);
}

int OversetDataProxy::FindNeighborIndex(PHVectorInt1D &neighborZoneIndexArray, int neighborZoneIndexForSearching)
{
    uint_t numberOfNeighborZones = neighborZoneIndexArray.size();
    for (int iNeighbor = 0; iNeighbor < numberOfNeighborZones; ++ iNeighbor)
    {
        int currentNeighborZoneIndex = neighborZoneIndexArray[iNeighbor];
        if (currentNeighborZoneIndex == neighborZoneIndexForSearching)
        {
            return iNeighbor;
        }
    }
    return - 1;
};

OversetInformationProxy::OversetInformationProxy()
{
    grid = 0;
    oversetCellContainer = 0;
    oversetSwapManager = new OversetSwapManager();

    oversetDataProxyForSend = new OversetDataProxy();
    oversetDataProxyForReceive = new OversetDataProxy();

    oversetDataProxyForSend->Initialize(oversetSwapManager, 0);
    oversetDataProxyForReceive->Initialize(oversetSwapManager, 1);
}

OversetInformationProxy::~OversetInformationProxy()
{
    FreeOversetCellContainer(oversetCellContainer);

    FreePointer(oversetSwapManager);
    FreePointer(oversetDataProxyForSend);
    FreePointer(oversetDataProxyForReceive);
}

int OversetInformationProxy::GetZoneIndex()
{
    return grid->GetZoneID();
}

PHVectorInt1D OversetInformationProxy::GetNeighborZoneIndexArray()
{
    return oversetSwapManager->GetOversetDonorCellManager()->GetNeighborInterplateZoneIndexArray();
}

void OversetInformationProxy::FreeOversetCellContainer(PHVector1D < OversetCell * > *oversetCellContainer)
{
    if (! oversetCellContainer) return;

    uint_t numberOfElements = oversetCellContainer->size();
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        delete (* oversetCellContainer)[iElement];

    }
    delete oversetCellContainer;
}

void OversetInformationProxy::EncodeOversetInformation(DataContainer *dataContainer)
{
    OversetSwapManager *oversetSwapManager = this->GetOversetSwapManager();
    oversetSwapManager->EncodeOversetInformation(dataContainer);
}

void OversetInformationProxy::DecodeOversetInformation(DataContainer *dataContainer)
{
    OversetSwapManager *oversetSwapManager = this->GetOversetSwapManager();
    oversetSwapManager->DecodeOversetInformation(dataContainer);
}

void OversetInformationProxy::GetCellIndexContainer(int donorZoneIndexIn, PHVectorInt1D &cellIndexContainer)
{
    PHVector1D < OversetCell * > *oversetCellContainer = this->GetOversetCellContainer();
    cellIndexContainer.resize(0);

    uint_t numberOfElements = oversetCellContainer->size();
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        OversetCell *oversetCell = (* oversetCellContainer)[iElement];
        int cellIndex = oversetCell->GetCellIndex();
        int donorZoneIndex = oversetCell->GetDonorZoneIndex();

        if (donorZoneIndex == donorZoneIndexIn)
        {
            cellIndexContainer.push_back(cellIndex);
        }
    }
}

void OversetInformationProxy::GetDonorCellIndexContainer(int donorZoneIndexIn, PHVectorInt1D &donorCellIndexContainer,
    PHVectorRDouble1D &interCellCenterX,
    PHVectorRDouble1D &interCellCenterY,
    PHVectorRDouble1D &interCellCenterZ,
    PHVectorRDouble1D &interCellVolume)
{
    //! Need further revise for structure mesh.
    UnstructGrid *unstructuredGrid = PHSPACE::UnstructGridCast(grid);
    RDouble *xc = unstructuredGrid->GetCellCenterX();
    RDouble *yc = unstructuredGrid->GetCellCenterY();
    RDouble *zc = unstructuredGrid->GetCellCenterZ();
    RDouble *vol  = unstructuredGrid->GetCellVolume();

    int nBoundaryFace = unstructuredGrid->GetNBoundFace();
    int *leftCellofFace  = unstructuredGrid->GetLeftCellOfFace();
    int *rightCellofFace = unstructuredGrid->GetRightCellOfFace();
    UnstructBCSet **unstructBCSet = unstructuredGrid->GetBCRecord();

    PHVector1D < OversetCell * > *oversetCellContainer = this->GetOversetCellContainer();
    donorCellIndexContainer.resize(0);

    uint_t numberOfElements = oversetCellContainer->size();

    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    int count = 0;
    if (isOversetSlip)
    {
        for (int iFace = 0; iFace < nBoundaryFace; iFace++)
        {
            int bcType = unstructBCSet[iFace]->GetKey();
            if (bcType == PHENGLEI::OVERSET)
            {
                int le = leftCellofFace[iFace];
                OversetCell *oversetCell = (* oversetCellContainer)[count];
                int donorCellIndex = oversetCell->GetDonorCellIndex();
                int donorZoneIndex = oversetCell->GetDonorZoneIndex();

                if (donorZoneIndex == donorZoneIndexIn)
                {
                    donorCellIndexContainer.push_back(donorCellIndex);
                    interCellCenterX.push_back(xc[le]);
                    interCellCenterY.push_back(yc[le]);
                    interCellCenterZ.push_back(zc[le]);
                    interCellVolume.push_back(vol[le]);
                }
                count++;
            }
        }
    }
    else 
    {
        for (int iElement = 0; iElement < numberOfElements; ++ iElement)
        {
            OversetCell *oversetCell = (* oversetCellContainer)[iElement];
            int cellIndex = oversetCell->GetCellIndex();
            int donorCellIndex = oversetCell->GetDonorCellIndex();
            int donorZoneIndex = oversetCell->GetDonorZoneIndex();

            if (donorZoneIndex == donorZoneIndexIn)
            {
                donorCellIndexContainer.push_back(donorCellIndex);
                interCellCenterX.push_back(xc[cellIndex]);
                interCellCenterY.push_back(yc[cellIndex]);
                interCellCenterZ.push_back(zc[cellIndex]);
                interCellVolume.push_back(vol[cellIndex]);
            }
        }
    }
}

void OversetInformationProxy::InitializeReceivingInformation(PHVector1D < OversetCell * > *oversetCellContainer)
{
    this->SetOversetCellContainer(oversetCellContainer);
}

void OversetInformationProxy::PostProcessReceivingOversetInformation()
{
    OversetSwapManager *oversetSwapManager = this->GetOversetSwapManager();
    oversetSwapManager->PostProcessReceivingOversetInformation();
}

void OversetInformationProxy::PostProcessSendingOversetInformation()
{
    OversetSwapManager *oversetSwapManager = this->GetOversetSwapManager();
    oversetSwapManager->PostProcessSendingOversetInformation();
}

void OversetInformationProxy::PostProcessOversetInformation()
{
    this->PostProcessReceivingOversetInformation();
    this->PostProcessSendingOversetInformation();
}

void OversetInformationProxy::GenerateNeighborZoneIndexSet(set< int > &neighborZoneIndexSet)
{
    PHVector1D < OversetCell * > *oversetCellContainer = this->GetOversetCellContainer();

    uint_t numberOfElements = oversetCellContainer->size();
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        OversetCell *oversetCell = (* oversetCellContainer)[iElement];
        int donorZoneIndex = oversetCell->GetDonorZoneIndex();
        neighborZoneIndexSet.insert(donorZoneIndex);
    }
}

void OversetInformationProxy::GetNeighborZoneIndexContainer(PHVectorInt1D &neighborZoneIndexContainer)
{
    set< int > neighborZoneIndexSet;
    this->GenerateNeighborZoneIndexSet(neighborZoneIndexSet);

    for (set< int >::iterator iter = neighborZoneIndexSet.begin(); iter != neighborZoneIndexSet.end(); ++ iter)
    {
        neighborZoneIndexContainer.push_back(* iter);
    }
}

OversetTopologyManager oversetTopologyManager;

OversetTopologyManager *GetOversetTopologyManager()
{
    return & oversetTopologyManager;
}

void InitializeOversetInterfaceTopology()
{
    OversetTopologyManager *oversetTopologyManagerTmp = PHSPACE::GetOversetTopologyManager();
    oversetTopologyManagerTmp->InitializeOversetInterfaceTopology();
}

OversetTopologyManager::OversetTopologyManager()
{
    ;
}

OversetTopologyManager::~OversetTopologyManager()
{
    ;
}

PHVectorInt2D &OversetTopologyManager::GetOversetTopology()
{
    return oversetTopology;
};

PHVectorInt2D &OversetTopologyManager::GetOversetTopologyForSend()
{
    return oversetTopologyForSend;
};

PHVectorInt2D &OversetTopologyManager::GetOversetTopologyForReceive()
{
    return oversetTopologyForReceive;
};

uint_t OversetTopologyManager::GetNumberOfOversetNeighbors(int zoneIndex)
{
    PHVectorInt2D &oversetTopology = this->GetOversetTopology();
    return oversetTopology[zoneIndex].size();
}

int OversetTopologyManager::GetNeighborOversetZoneIndex(int zoneIndex, int ithOversetNeighbor)
{
    PHVectorInt2D &oversetTopology = this->GetOversetTopology();
    return oversetTopology[zoneIndex][ithOversetNeighbor];
}

void OversetTopologyManager::InitializeOversetInterfaceTopology()
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    PHVectorInt2D &oversetTopology = this->GetOversetTopology();
    PHVectorInt2D &oversetTopologyForSend = this->GetOversetTopologyForSend();
    PHVectorInt2D &oversetTopologyForReceive = this->GetOversetTopologyForReceive();

    oversetTopologyForSend.resize(nZones);
    oversetTopologyForReceive.resize(nZones);

    oversetTopology.resize(nZones);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        int currentProcessor = PHMPI::GetCurrentProcessorID();

        if (processorIndex != currentProcessor) continue;

        OversetInformationProxy *oversetInformationProxy = GetOversetInformationProxy(iZone, 0);
        if (oversetInformationProxy == 0) continue;

        OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();

        PHVectorInt1D neighborZoneIndexArray = oversetDataProxy->GetNeighborZoneIndexArray();

        oversetTopology[iZone] = neighborZoneIndexArray;
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        int numberOfNeighbors = static_cast<int>(oversetTopology[iZone].size());

        PHMPI::PH_Bcast(&numberOfNeighbors, sizeof(int), processorIndex);

        oversetTopology[iZone].resize(numberOfNeighbors);

        if (numberOfNeighbors <= 0) continue;

        PHMPI::PH_Bcast(&oversetTopology[iZone][0], numberOfNeighbors * sizeof(int), processorIndex);
    }
}

void OversetTopologyManager::InitializeOversetInterfaceTopologyServerCollectAndBcast()
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int serverProcessor = PHMPI::GetServerSepMode();

    PHVectorInt2D &oversetTopology = this->GetOversetTopology();
    PHVectorInt2D &oversetTopologyForSend = this->GetOversetTopologyForSend();
    PHVectorInt2D &oversetTopologyForReceive = this->GetOversetTopologyForReceive();

    oversetTopologyForSend.resize(nZones);
    oversetTopologyForReceive.resize(nZones);
    oversetTopology.resize(nZones);

    DataContainer *dataContainerLocal = new DataContainer();
    DataContainer *dataContainerGlobal = new DataContainer();

    int nZonesInThisProc = 0;

    //! Get the number of blocks in each process nZonesInThisProc.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessor == processorIndex)
        {
            ++ nZonesInThisProc;
        }
    }

    PHSPACE::PHWrite(dataContainerLocal, nZonesInThisProc);

    //! For each process,get the interCell information of all zones and write into dataContainer.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessor == processorIndex)
        {
            OversetInformationProxy *oversetInformationProxy = GetOversetInformationProxy(iZone, 0);

            OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();

            PHVectorInt1D neighborZoneIndexArray = oversetDataProxy->GetNeighborZoneIndexArray();

            int numberOfNeigibors = static_cast<int>( neighborZoneIndexArray.size() );

            PHSPACE::PHWrite(dataContainerLocal, iZone);

            PHSPACE::PHWrite(dataContainerLocal, numberOfNeigibors);

            if (numberOfNeigibors != 0)
            {
                PHSPACE::PHWrite(dataContainerLocal, & neighborZoneIndexArray[0], numberOfNeigibors);
            }
        }
    }

    //! The main process collects and read information.
    int numberOfprocessors = PHMPI::GetNumberOfProcessor();

    if (currentProcessor != serverProcessor)
    {
        send(dataContainerLocal, serverProcessor, currentProcessor);
    }

    if (currentProcessor == serverProcessor)
    {
        dataContainerLocal->MoveToBegin();

        dataContainerGlobal->Write(dataContainerLocal);

        for (int iProc = 0; iProc < numberOfprocessors; ++ iProc)
        {
            if (iProc == serverProcessor)
            {
                continue;
            }

            DataContainer *dataContainerRevc = new DataContainer();

            receive(dataContainerRevc, iProc, iProc);

            dataContainerRevc->MoveToBegin();

            dataContainerGlobal->Write(dataContainerRevc);

            FreePointer(dataContainerRevc);
        }
    }

    FreePointer(dataContainerLocal);

    PHMPI::PH_BcastByServerSepMode(dataContainerGlobal);
    dataContainerGlobal->MoveToBegin();

    for (int iProc = 0; iProc < numberOfprocessors; ++ iProc)
    {
        PHSPACE::PHRead(dataContainerGlobal, nZonesInThisProc);

        for (int iZone = 0; iZone < nZonesInThisProc; ++ iZone)
        {
            int zoneIndex;
            int numberOfNeigibors;

            PHSPACE::PHRead(dataContainerGlobal, zoneIndex);
            PHSPACE::PHRead(dataContainerGlobal, numberOfNeigibors);

            if (numberOfNeigibors != 0)
            {
                PHVectorInt1D neighborList;
                neighborList.resize(numberOfNeigibors);
                PHSPACE::PHRead(dataContainerGlobal, &neighborList[0], numberOfNeigibors);
                oversetTopology[zoneIndex] = neighborList;
            }
        }
    }

    FreePointer(dataContainerGlobal);
}

OversetCell::OversetCell()
{
    ;
}

OversetCell::~OversetCell()
{
    ;
}

OversetInterplateCell::OversetInterplateCell()
{
    dataStorage = new Data_ParamFieldSuite();
}

OversetInterplateCell::~OversetInterplateCell()
{
    FreePointer(dataStorage);
}

void OversetInterplateCell::EncodeOversetInterplateCellInformation(DataContainer *dataContainer)
{
    int donorZoneIndex = this->GetDonorZoneIndex();

    PHSPACE::PHWrite(dataContainer, donorZoneIndex);

    int numberOfOversetInterplateCells = this->GetNumberOfOversetInterplateCells();
    PHSPACE::PHWrite(dataContainer, numberOfOversetInterplateCells);

    PHVectorInt1D &cellIndexContainer = this->GetCellIndexContainer();

    PHSPACE::PHWrite(dataContainer, cellIndexContainer, numberOfOversetInterplateCells);

    PHVectorInt1D &storageIndexContainer = this->GetStorageIndexContainer();

    PHSPACE::PHWrite(dataContainer, storageIndexContainer, numberOfOversetInterplateCells);
}

void OversetInterplateCell::DecodeOversetInterplateCellInformation(DataContainer *dataContainer)
{
    int donorZoneIndex = - 1;

    PHSPACE::PHRead(dataContainer, donorZoneIndex);

    this->SetDonorZoneIndex(donorZoneIndex);

    int numberOfOversetInterplateCells;
    PHSPACE::PHRead(dataContainer, numberOfOversetInterplateCells);

    PHVectorInt1D &cellIndexContainer = this->GetCellIndexContainer();
    cellIndexContainer.resize(numberOfOversetInterplateCells);

    PHSPACE::PHRead(dataContainer, cellIndexContainer, numberOfOversetInterplateCells);

    PHVectorInt1D &storageIndexContainer = this->GetStorageIndexContainer();
    storageIndexContainer.resize(numberOfOversetInterplateCells);

    PHSPACE::PHRead(dataContainer, storageIndexContainer, numberOfOversetInterplateCells);
}

OversetDonorCell::OversetDonorCell()
{
    dataStorage = new Data_ParamFieldSuite();
}

OversetDonorCell::~OversetDonorCell()
{
    FreePointer(dataStorage);
}

void OversetDonorCell::EncodeOversetDonorCellInformation(DataContainer *dataContainer)
{
    int interplateZoneIndex = this->GetInterplateZoneIndex();

    PHSPACE::PHWrite(dataContainer, interplateZoneIndex);

    int numberOfOversetDonorCells = this->GetNumberOfOversetDonorCells();
    PHSPACE::PHWrite(dataContainer, numberOfOversetDonorCells);

    PHVectorInt1D &donorCellIndexContainer = this->GetDonorCellIndexContainer();

    PHSPACE::PHWrite(dataContainer, donorCellIndexContainer, numberOfOversetDonorCells);

    PHSPACE::PHWrite(dataContainer, interCellCenterX, numberOfOversetDonorCells);
    PHSPACE::PHWrite(dataContainer, interCellCenterY, numberOfOversetDonorCells);
    PHSPACE::PHWrite(dataContainer, interCellCenterZ, numberOfOversetDonorCells);

    PHVectorInt1D &storageIndexContainer = this->GetStorageIndexContainer();

    PHSPACE::PHWrite(dataContainer, storageIndexContainer, numberOfOversetDonorCells);
}

void OversetDonorCell::DecodeOversetDonorCellInformation(DataContainer *dataContainer)
{
    int interplateZoneIndex = - 1;

    PHSPACE::PHRead(dataContainer, interplateZoneIndex);

    this->SetInterplateZoneIndex(interplateZoneIndex);

    int numberOfOversetDonorCells;
    PHSPACE::PHRead(dataContainer, numberOfOversetDonorCells);

    PHVectorInt1D &donorCellIndexContainer = this->GetDonorCellIndexContainer();

    donorCellIndexContainer.resize(numberOfOversetDonorCells);

    interCellCenterX.resize(numberOfOversetDonorCells);
    interCellCenterY.resize(numberOfOversetDonorCells);
    interCellCenterZ.resize(numberOfOversetDonorCells);

    PHSPACE::PHRead(dataContainer, donorCellIndexContainer, numberOfOversetDonorCells);

    PHSPACE::PHRead(dataContainer, interCellCenterX, numberOfOversetDonorCells);
    PHSPACE::PHRead(dataContainer, interCellCenterY, numberOfOversetDonorCells);
    PHSPACE::PHRead(dataContainer, interCellCenterZ, numberOfOversetDonorCells);

    PHVectorInt1D &storageIndexContainer = this->GetStorageIndexContainer();
    storageIndexContainer.resize(numberOfOversetDonorCells);

    PHSPACE::PHRead(dataContainer, storageIndexContainer, numberOfOversetDonorCells);
}

OversetInterplateCellManager::OversetInterplateCellManager()
{
    dataStorage = new Data_ParamFieldSuite();
    totalNumberOfOversetInterplateCells = 0;
}

OversetInterplateCellManager::~OversetInterplateCellManager()
{
    for (unsigned int i = 0; i < oversetInterplateCellContainer.size(); ++ i)
    {
        delete oversetInterplateCellContainer[i];
    }
    FreePointer(dataStorage);
}

int OversetInterplateCellManager::GetNumberOfOversetInterplateCells(int iNeighbor)
{
    return oversetInterplateCellContainer[iNeighbor]->GetNumberOfOversetInterplateCells();
}

int *OversetInterplateCellManager::GetOversetInterplateCell(int iNeighbor)
{
    PHVectorInt1D &cellIndexContainer = oversetInterplateCellContainer[iNeighbor]->GetCellIndexContainer();
    return & cellIndexContainer[0];
}

int *OversetInterplateCellManager::GetInterplateCellStorageIndex(int iNeighbor)
{
    PHVectorInt1D &storageIndexContainer = oversetInterplateCellContainer[iNeighbor]->GetStorageIndexContainer();
    return & storageIndexContainer[0];
}

PHVectorInt1D OversetInterplateCellManager::GetNeighborDonorZoneIndexArray()
{
    PHVectorInt1D neighborDonorZoneIndexArray;

    uint_t numberOfNeighbors = this->GetNumberOfNeighbors();

    for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
    {
        int donorZoneIndex = oversetInterplateCellContainer[iNeighbor]->GetDonorZoneIndex();
        neighborDonorZoneIndexArray.push_back(donorZoneIndex);
    }

    return neighborDonorZoneIndexArray;
}

void OversetInterplateCellManager::PostProcessReceivingOversetInformation()
{
    PHVector1D < OversetInterplateCell * > &oversetInterplateCellContainer = this->GetOversetInterplateCellContainer();
    uint_t numberOfNeighbors = oversetInterplateCellContainer.size();

    int totalNumberOfOversetInterplateCells = 0;

    for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
    {
        OversetInterplateCell *oversetInterplateCell = oversetInterplateCellContainer[iNeighbor];

        PHVectorInt1D &storageIndexContainer = oversetInterplateCell->GetStorageIndexContainer();
        uint_t numberOfOversetInterplateCells = oversetInterplateCell->GetNumberOfOversetInterplateCells();

        for (int iOversetInterplateCells = 0; iOversetInterplateCells < numberOfOversetInterplateCells; ++ iOversetInterplateCells)
        {
            storageIndexContainer.push_back(totalNumberOfOversetInterplateCells ++);
        }
    }

    this->SetTotalNumberOfOversetInterplateCells(totalNumberOfOversetInterplateCells);
}

void OversetInterplateCellManager::EncodeOversetInterplateCellInformation(DataContainer *dataContainer)
{
    PHVector1D < OversetInterplateCell * > &oversetInterplateCellContainer = this->GetOversetInterplateCellContainer();

    int numberOfOversetNeighbors = static_cast<int>( oversetInterplateCellContainer.size() );
    PHSPACE::PHWrite(dataContainer, numberOfOversetNeighbors);

    int totalNumberOfOversetInterplateCells = this->GetTotalNumberOfOversetInterplateCells();
    PHSPACE::PHWrite(dataContainer, totalNumberOfOversetInterplateCells);

    for (int iOversetNeighbor = 0; iOversetNeighbor < numberOfOversetNeighbors; ++ iOversetNeighbor)
    {
        OversetInterplateCell *oversetInterplateCell = oversetInterplateCellContainer[iOversetNeighbor];
        oversetInterplateCell->EncodeOversetInterplateCellInformation(dataContainer);
    }
}

void OversetInterplateCellManager::DecodeOversetInterplateCellInformation(DataContainer *dataContainer)
{
    PHVector1D < OversetInterplateCell * > &oversetInterplateCellContainer = this->GetOversetInterplateCellContainer();

    int numberOfOversetNeighbors;
    PHSPACE::PHRead(dataContainer, numberOfOversetNeighbors);

    oversetInterplateCellContainer.resize(numberOfOversetNeighbors);

    int totalNumberOfOversetInterplateCells;
    PHSPACE::PHRead(dataContainer, totalNumberOfOversetInterplateCells);

    this->SetTotalNumberOfOversetInterplateCells(totalNumberOfOversetInterplateCells);

    for (int iOversetNeighbor = 0; iOversetNeighbor < numberOfOversetNeighbors; ++ iOversetNeighbor)
    {
        OversetInterplateCell *oversetInterplateCell = new OversetInterplateCell();

        oversetInterplateCell->DecodeOversetInterplateCellInformation(dataContainer);

        oversetInterplateCellContainer[iOversetNeighbor] = oversetInterplateCell;
    }
}

void OversetInterplateCellManager::AddOversetInterplateCell(OversetInterplateCell *oversetInterplateCell)
{
    PHVector1D < OversetInterplateCell * > &oversetInterplateCellContainer = this->GetOversetInterplateCellContainer();
    oversetInterplateCellContainer.push_back(oversetInterplateCell);
}

OversetDonorCellManager::OversetDonorCellManager()
{
    dataStorage = new Data_ParamFieldSuite();
    totalNumberOfOversetDonorCells = 0;
}

OversetDonorCellManager::~OversetDonorCellManager()
{
    for (unsigned int i = 0; i < oversetDonorCellContainer.size(); ++ i)
    {
        delete oversetDonorCellContainer[i];
    }
    delete dataStorage;
}

int OversetDonorCellManager::GetNumberOfOversetDonorCells(int iNeighbor)
{
    return oversetDonorCellContainer[iNeighbor]->GetNumberOfOversetDonorCells();
}

int *OversetDonorCellManager::GetOversetDonorCell(int iNeighbor)
{
    PHVectorInt1D &donorCellIndexContainer = oversetDonorCellContainer[iNeighbor]->GetDonorCellIndexContainer();
    return & donorCellIndexContainer[0];
}

RDouble *OversetDonorCellManager::GetInterCellCenterX(int iNeighbor)
{
    PHVectorRDouble1D &interCellCenterX = oversetDonorCellContainer[iNeighbor]->GetInterCellCenterX();
    return & interCellCenterX[0];
}

RDouble *OversetDonorCellManager::GetInterCellCenterY(int iNeighbor)
{
    PHVectorRDouble1D &interCellCenterY = oversetDonorCellContainer[iNeighbor]->GetInterCellCenterY();
    return & interCellCenterY[0];
}

RDouble *OversetDonorCellManager::GetInterCellCenterZ(int iNeighbor)
{
    PHVectorRDouble1D &interCellCenterZ = oversetDonorCellContainer[iNeighbor]->GetInterCellCenterZ();
    return & interCellCenterZ[0];
}
RDouble *OversetDonorCellManager::GetInterCellVolume(int iNeighbor)
{
    PHVectorRDouble1D &interCellVolume = oversetDonorCellContainer[iNeighbor]->GetInterCellVolume();
    return & interCellVolume[0];
}
int *OversetDonorCellManager::GetDonorCellStorageIndex(int iNeighbor)
{
    PHVectorInt1D &storageIndexContainer = oversetDonorCellContainer[iNeighbor]->GetStorageIndexContainer();
    return & storageIndexContainer[0];
}

PHVectorInt1D OversetDonorCellManager::GetNeighborInterplateZoneIndexArray()
{
    PHVectorInt1D neighborInterplateZoneIndexArray;

    uint_t numberOfNeighbors = this->GetNumberOfNeighbors();

    for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
    {
        int interplateZoneIndex = oversetDonorCellContainer[iNeighbor]->GetInterplateZoneIndex();
        neighborInterplateZoneIndexArray.push_back(interplateZoneIndex);
    }

    return neighborInterplateZoneIndexArray;
}

void OversetDonorCellManager::PostProcessSendingOversetInformation()
{
    PHVector1D < OversetDonorCell * > &oversetDonorCellContainer = this->GetOversetDonorCellContainer();
    uint_t numberOfNeighbors = oversetDonorCellContainer.size();

    int totalNumberOfOversetDonorCells = 0;
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
    {
        OversetDonorCell *oversetDonorCell = oversetDonorCellContainer[iNeighbor];

        PHVectorInt1D &storageIndexContainer = oversetDonorCell->GetStorageIndexContainer();
        int numberOfOversetDonorCells = oversetDonorCell->GetNumberOfOversetDonorCells();

        for (int iOversetDonorCell = 0; iOversetDonorCell < numberOfOversetDonorCells; ++ iOversetDonorCell)
        {
            storageIndexContainer.push_back(totalNumberOfOversetDonorCells ++);
        }
    }

    this->SetTotalNumberOfOversetDonorCells(totalNumberOfOversetDonorCells);
}

void OversetDonorCellManager::EncodeOversetDonorCellInformation(DataContainer *dataContainer)
{
    PHVector1D < OversetDonorCell * > &oversetDonorCellContainer = this->GetOversetDonorCellContainer();

    int numberOfOversetNeighbors = static_cast<int>( oversetDonorCellContainer.size() );
    PHSPACE::PHWrite(dataContainer, numberOfOversetNeighbors);

    int totalNumberOfOversetDonorCells = this->GetTotalNumberOfOversetDonorCells();
    PHSPACE::PHWrite(dataContainer, totalNumberOfOversetDonorCells);

    for (int iOversetNeighbor = 0; iOversetNeighbor < numberOfOversetNeighbors; ++ iOversetNeighbor)
    {
        OversetDonorCell *oversetDonorCell = oversetDonorCellContainer[iOversetNeighbor];
        oversetDonorCell->EncodeOversetDonorCellInformation(dataContainer);
    }
}

void OversetDonorCellManager::DecodeOversetDonorCellInformation(DataContainer *dataContainer)
{
    PHVector1D < OversetDonorCell * > &oversetDonorCellContainer = this->GetOversetDonorCellContainer();

    int numberOfOversetNeighbors;
    PHSPACE::PHRead(dataContainer, numberOfOversetNeighbors);

    oversetDonorCellContainer.resize(numberOfOversetNeighbors);

    int totalNumberOfOversetDonorCells;
    PHSPACE::PHRead(dataContainer, totalNumberOfOversetDonorCells);

    this->SetTotalNumberOfOversetDonorCells(totalNumberOfOversetDonorCells);

    for (int iOversetNeighbor = 0; iOversetNeighbor < numberOfOversetNeighbors; ++ iOversetNeighbor)
    {
        OversetDonorCell *oversetDonorCell = new OversetDonorCell();

        oversetDonorCell->DecodeOversetDonorCellInformation(dataContainer);

        oversetDonorCellContainer[iOversetNeighbor] = oversetDonorCell;
    }
}

void OversetDonorCellManager::AddOversetDonorCell(OversetDonorCell *oversetDonorCell)
{
    PHVector1D < OversetDonorCell * > &oversetDonorCellContainer = GetOversetDonorCellContainer();
    oversetDonorCellContainer.push_back(oversetDonorCell);
}

uint_t GetNumberOfOversetNeighbors(int zoneIndex)
{
    OversetTopologyManager *oversetTopologyManagerTmp = PHSPACE::GetOversetTopologyManager();
    uint_t numberOfOversetNeighbors = oversetTopologyManagerTmp->GetNumberOfOversetNeighbors(zoneIndex);
    return numberOfOversetNeighbors;
}

int GetNeighborOversetZoneIndex(int zoneIndex, int iNeighbor)
{
    OversetTopologyManager *oversetTopologyManagerTmp = PHSPACE::GetOversetTopologyManager();
    int neighborOversetZoneIndex = oversetTopologyManagerTmp->GetNeighborOversetZoneIndex(zoneIndex, iNeighbor);
    return neighborOversetZoneIndex;
}

}
