#include "Task_ServerUpdateInterface.h"
#include "Task.h"
#include "Geo_Grid.h"
#include "TK_Log.h"

namespace PHSPACE
{

void MainTaskBlocking(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);
            PH_Interface(actkey, iZone, neighborZoneIndex);
        }
    }
}

void MainTaskNonBlocking(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    vector <vector <DataContainer *> > receivedDataBuffer;
    vector <vector <DataContainer *> > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    vector <PH_Request> requestContainer;

    //! Each processor runs along the same order, which is from 0 to nZones.
    //! Each simulation zone is looped over to judge if the communication is necessary.
    //! Step 1: Compressing data into actkey and Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessID = GetZoneProcessorIDSepMode(iZone);
            int receiveProcessID = GetZoneProcessorIDSepMode(neighborZoneIndex);

            if (currentProcessorID != sendProcessID && currentProcessorID != receiveProcessID)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessorID == receiveProcessID)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessorID == sendProcessID)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone].push_back(sendBuffer);
            }

            //! Swapping.
            PH_InterfaceNonBlocking(actkey, sendBuffer, receivedBuffer, iZone, neighborZoneIndex, requestContainer);
        }
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translate the data container.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);
            int receiveProcessID = GetZoneProcessorIDSepMode(neighborZoneIndex);
            if (currentProcessorID == receiveProcessID)
            {
                DataContainer *originalData = actkey->GetData();

                actkey->SetData(receivedDataBuffer[iZone][count]);
                actkey->ipos = iZone;
                GeneralTranslateAction(actkey, neighborZoneIndex);

                actkey->SetData(originalData);
                ++ count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
        receivedDataBuffer[iDim].clear();
    }
    receivedDataBuffer.clear();

    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
        sendDataBuffer[iDim].clear();
    }

    sendDataBuffer.clear();
    requestContainer.clear();
}

}