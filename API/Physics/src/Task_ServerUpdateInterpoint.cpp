#include "Task_ServerUpdateInterpoint.h"
#include "Geo_Grid.h"
#pragma warning(disable:4100)
namespace PHSPACE
{
TaskServerUpdateInterpoint::TaskServerUpdateInterpoint() {}

TaskServerUpdateInterpoint::~TaskServerUpdateInterpoint() {}

void TaskServerUpdateInterpoint::PreTask(ActionKey *actkey) {}

void TaskServerUpdateInterpoint::PostTask(ActionKey *actkey) {}

void TaskServerUpdateInterpoint::MainTask(ActionKey *actkey)
{
    this->MainTaskNonBlocking(actkey);
}

void TaskServerUpdateInterpoint::MainTaskBlocking(ActionKey *actkey)
{
    //using namespace PHMPI;
    //int nZones = GetNumberofGlobalZones();

    //for (int iZone = 0; iZone < nZones; ++ iZone)
    //{
    //    ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);

    //    std::size_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

    //    for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    //    {
    //        int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);

    //        PH_Interpoint(actkey, iZone, neighborZoneIndex);
    //    }
    //}
}

void TaskServerUpdateInterpoint::MainTaskNonBlocking(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int currentProcessor = GetCurrentProcessorID();

    vector<vector<DataContainer *> > receivedDataBuffer;
    vector<vector<DataContainer *> > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    vector<PH_Request> requestContainer;

    //! Each processor runs along the same order, which is from 0 to nZones.
    //! Each simulation zone is looped over to judge if the communication is necessary.
    //! Step 1: Compressing data into actkey and Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZoneIndexForPoint = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor = GetZoneProcessorIDSepMode(iZone);
            int receiveProcessor = GetZoneProcessorIDSepMode(neighborZoneIndexForPoint);

            if (currentProcessor != sendProcessor &&currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
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

            //! Swapping.
            InterpointNonBlocking(actkey, sendBuffer, receivedBuffer, iZone, neighborZoneIndexForPoint, requestContainer);
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
        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t ineighbor = 0; ineighbor < numberOfNeighbor; ++ ineighbor)
        {
            int neighborZoneIndexForPoint = zoneNeighborForPoint->GetZoneIndexOfNeighbor(ineighbor);

            int receiveProcessor = GetZoneProcessorIDSepMode(neighborZoneIndexForPoint);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *originalData = actkey->GetData();

                actkey->SetData(receivedDataBuffer[iZone][count]);
                actkey->ipos = iZone;
                GeneralTranslateAction(actkey,neighborZoneIndexForPoint);

                actkey->SetData(originalData);
                ++ count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (std::size_t iDimension = 0; iDimension < receivedDataBuffer.size(); ++ iDimension)
    {
        for (std::size_t jDimension = 0; jDimension < receivedDataBuffer[iDimension].size(); ++ jDimension)
        {
            delete receivedDataBuffer[iDimension][jDimension];
        }
        receivedDataBuffer[iDimension].clear();
    }
    receivedDataBuffer.clear();

    for (std::size_t iDimension = 0; iDimension < sendDataBuffer.size(); ++ iDimension)
    {
        for (std::size_t jDimension = 0; jDimension < sendDataBuffer[iDimension].size(); ++ jDimension)
        {
            delete sendDataBuffer[iDimension][jDimension];
        }
        sendDataBuffer[iDimension].clear();
    }
    sendDataBuffer.clear();

    requestContainer.clear();
}

void TaskServerUpdateInterpoint::InterpointNonBlocking(ActionKey *actkey, DataContainer *sendBuffer, DataContainer *receivedBuffer, int iZone, int jZone, vector<PH_Request> &requestContainer)
{
    using namespace PHMPI;

    int currentProcessor = GetCurrentProcessorID();
    int sendProcessor = GetZoneProcessorIDSepMode(iZone);
    int receiveProcessor = GetZoneProcessorIDSepMode(jZone);

    int tag = GetSendRecvTag(actkey, iZone);

    //! Compress the send information into the actkey.
    if (currentProcessor == sendProcessor)
    {    
        DataContainer *originalData = actkey->GetData();

        DataContainer *compressData = sendBuffer;
        if (sendProcessor == receiveProcessor) compressData = receivedBuffer;

        //! Fill the send buffer by actkey.
        actkey->SetData(compressData);
        actkey->ipos = jZone;
        GeneralAction(actkey, iZone);

        actkey->SetData(originalData);
    }


    if (sendProcessor == receiveProcessor)
    {
        return;
    }

    //! Communication.
    if (currentProcessor == sendProcessor)
    {
        streamsize nlen = sendBuffer->Size();

        //! Send the data to neighbors.
        send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
    }
    else if (currentProcessor == receiveProcessor)
    {
        //! Pre-calculating the data length of the received data.
        CharVecSizeType nlen = InterpointGetLength(actkey, iZone, jZone);

        //! Receive data from neighbors.
        receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
    }
}

CharVecSizeType TaskServerUpdateInterpoint::InterpointGetLength(ActionKey *actkey, int iZone, int jZone)
{
    using namespace PHMPI;
    CharVecSizeType nlen = 0;

    int sendProcessor = GetZoneProcessorIDSepMode(iZone);
    int receiveProcessor = GetZoneProcessorIDSepMode(jZone);

    int myID = GetCurrentProcessorID();    
    if (myID == sendProcessor)
    {
        DataContainer *cdata = actkey->GetData();
        nlen = cdata->Size();
    }
    if (myID == receiveProcessor)
    {
        actkey->ipos = iZone;
        nlen = static_cast<CharVecSizeType>(GeneralTranslateActionLength(actkey, jZone));
    }
    return nlen;
}

}