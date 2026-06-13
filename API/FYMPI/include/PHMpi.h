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
//! @file      PHMPI.h
//! @brief     It defines many basic global variables and interfaces which are related to MPI.
//! @author    Bell, He Xin, Xu Qingxin.
#pragma once

#ifdef PH_PARALLEL
    #include "mpi.h"
#endif

#include "LIB_Macro.h"
#include "DataStruct_BinaryTree.h"
#include "ActionKey.h"
#include "Parallel.h"

using namespace std;

namespace PHSPACE
{

//! Parallel or file mode:
//!    CS_MODE(0) : Client-Server mode.
//!    P2P_MODE(1): Peer-To-Peer mode.
const int CS_MODE  = 0;
const int P2P_MODE = 1;

//! Zone block distribution method.
//! -# DETERMIN_BY_PARTITION(0): each zone is assigned to the one that defined in grid partition procedure.
//! -# REDISTRIBUTION(1)       : random assigned for each zone or by some else ways.
const int DETERMIN_BY_PARTITION = 0;
const int REDISTRIBUTION        = 1;

namespace PHMPI
{
extern int csmode;
extern int server;
extern VVInt zone_over_link;

typedef void (*DATA_COMPRESS)(DataContainer *&cdata);
typedef void (*DATA_DECOMPRESS)(DataContainer *cdata);

//! PHengLEI solver initialization.
void Initialization();

void SetServerProcessorID(int serverId);

void SetGridServerProcessorID(int serverId);

//! Set the number of global zones, which make up the whole computation domain,
//! no matter the domain is divided into many partition or not. The number of global zone
//! is equal to the sum of number of local zone on each processor.
LIB_EXPORT void SetNumberOfGlobalZones(int nblock);

//! Set the number of processors used in this simulation process.
void SetNumberOfProcessor(int nproc);

void SetNumberOfGridProcessor(int nproc);

//! Set parallel Strategy:
//! -# 0 : each zone is assigned to the one that defined in grid partition procedure.
//! -# 1 : random assigned for each zone or by some else ways.
void SetParallelStrategy(int parallelStrategy);
int  GetParallelStrategy();

//! Set the number of local zones on the current processor.
//! The sum of local zones of each processor is equal to the number of global zones.
LIB_EXPORT void SetNumberofLocalZones(int nblock);

//! Set the number of files, that are the groups to manage grid/visualization/flow files, etc.
void SetNumberofFiles(int fileNumber);

//! Set the corresponding processor ID to be 'id' of iZone-th zone.
LIB_EXPORT void SetZoneProcessorID(int iZone, int id);

void SetZoneProcessID_INP(int *data);

void SetIsFileServer(int isServer);

//! Set zone block distribution method.
void SetZoneDistributionMethod(int method);

void InsertGlobalZones(int globalZoneID);

void BuildProcessorTree();

//! Special functions for HYBRID solver.
void SetNumberOfProcStructUsed(int number);

int GetServerProcessorID();

//! Return the processor index of the current processor, which is start from 0 and 
//! end with (nProcessor-1).
LIB_EXPORT int GetCurrentProcessorID();

int GetGridServerProcessorID();

//! Return the number of global zones, which make up the whole computation domain,
//! no matter the domain is divided into many partition or not. The number of global zone
//! is equal to the sum of number of local zone on each processor.
LIB_EXPORT int GetNumberofGlobalZones();

//! Special functions for HYBRID solver.
int * GetNZonesInProcessor();

//! Return the number of processors used in this simulation process.
LIB_EXPORT int GetNumberOfProcessor();

LIB_EXPORT int GetNumberOfGridProcessor();

//! Return the number of local zones on the current processor.
//! The sum of local zones of each processor is equal to the number of global zones.
LIB_EXPORT int GetNumberofLocalZones();

//! Return the number of files, that are the groups to manage grid/visualization/flow files, etc.
LIB_EXPORT int GetNumberofFiles();

//! Return the zone processor ID array pointer.
LIB_EXPORT int * GetZoneProcessorID();

//! Return the corresponding processor ID of iZone-th zone.
LIB_EXPORT int GetZoneProcessorID(int iZone);

LIB_EXPORT int * GetZoneProcessorIDForGrid();

LIB_EXPORT int GetZoneProcessorIDForGrid(int iZone);

//! Return the zone grid type array pointer.
LIB_EXPORT int * GetZoneGridType();

//! Return the grid type of iZone-th zone.
LIB_EXPORT int GetZoneGridType(int iZone);

//! Return the zone block ID array pointer.
LIB_EXPORT int * GetZoneGridID();

//! Return the zone file ID array pointer.
LIB_EXPORT int * GetZoneFileID();

//! Return the corresponding file manager processor ID of iZone-th zone.
//! File manage processor: which is used to read or write information of iZone-th zone.
//! The zone ID is not necessary equal to the file manager processor ID.
LIB_EXPORT int GetZoneFileID(int iZone);

//! This function is special used in HYBRID solver.
LIB_EXPORT int * GetZoneProcessorID_INP();

//! This function is special used in HYBRID solver.
LIB_EXPORT int * GetZoneProcessorIDDump();

LIB_EXPORT int * GetZoneGroupID();

//! Return zone block distribution method.
//! -# 0 : each zone is assigned to the one that defined in grid partition procedure.
//! -# 1 : random assigned for each zone or by some else ways.
LIB_EXPORT int GetZoneDistributionMethod();

//! Return the global zone index of local zone index = localZoneID.
LIB_EXPORT int GetLocalZoneIDToGlobalZoneID(int localZoneID);

//! Return the local zone index of global zone index = localZoneID.
LIB_EXPORT int GetGlobalZoneIDToLocalZoneID(int globalZoneID);

//! Special functions for HYBRID solver.
int GetNumberOfProcStructUsed();

void CreateZoneProcessorID(int nZonesTmp);

void CreatetZoneProcessorID_INP(int nZonesTmp);

void CreateZoneGridID(int nZonesTmp);

void CreateZoneGridType(int nZonesTmp);

void CreateZoneFileID(int nZonesTmp);

//! Special functions for HYBRID solver.
void CreateZoneProcessorIDDump(int nZonesTmp);

void CreateZoneGroupID(int nZonesTmp);

void FreeBlockData();

void CreateZoneStartPointContainer(int nZonesTmp);

int * GetZoneStartPointContainer();

void DeleteZoneStartPointContainer();

void CreateZoneStartCenterContainer(int nZonesTmp);

int * GetZoneStartCenterContainer();

void DeleteZoneStartCenterContainer();

void GenerateGlobalBinaryTree(int * zoneStartPointContainerTmp, int nZonesTmp, int numberOfPoints);

//! The following 8 functions should be moved to overset file.
int GetZoneIndexAccordingToGlobalPointLabel(int globalPointLabel);

void DeleteGlobalBinaryTree();

//! Return true if number of processors > 1.
LIB_EXPORT bool IsParallelRun();

LIB_EXPORT bool CurrentProcessorIsGridProcessor();

//! Return the interval of processors that manage files.
LIB_EXPORT int GetIntervalofProcessorFiles();

//! Return the file index that manage current processor data.
LIB_EXPORT int GetFileIndexofCurrentProcessor();

LIB_EXPORT int GetZoneProcessorIDSepMode(int iZone);

//! Return the send or receive processor ID of iZone-th zone.
LIB_EXPORT void GetSendRecvProc(int iZone, int &send_proc, int &recv_proc);

int GetServerSepMode();

LIB_EXPORT int * GetZoneProcessorIDSepMode();

//! Return the file manager processor ID, two cases of file mode:
//! Case CS_MODE (0): the processor ID is the server ID (0).
//! Case P2P_MODE(1): the processor ID is the current processor.
LIB_EXPORT void GetFileProc(int &file_proc);

//! Return if the current processor if file server.
LIB_EXPORT int IsCurrentProcessorFileServer();

void PH_Read(fstream &file, void *data, int size, int proc);

void PH_Read(fstream &file, void *data, int size);

void PH_ReadBcastInteger(fstream &file, void *data, int size, int proc);

void PH_ReadBcastDouble(fstream &file, void *data, int size, int proc);

void PH_Read_Bcast(fstream &file, void *data, int size, int proc);

void PH_BcastToClinet(void *data, int size, int proc, int key, int tag = 0);

void PH_BcastSepMode(void *data, int size, int proc, int tag = 0);

void PH_Bcast(void *data, int size, int proc, int tag = 0);

void PH_Bcast(DataContainer *cdata, int proc, int tag = 0);

void PH_BcastByServer(DataContainer *cdata, int tag = 0);

void PH_GatherByServer(DataContainer *cdata, int tag = 0);

void PH_BcastByServerSepMode(DataContainer *cdata, int tag = 0);

void PH_BcastByServerToClinet(DataContainer *cdata, int key, int tag = 0);

void PH_BcastNonBlocking(DataContainer *&cdata, int bcastProcessor, int tag);

void PH_BcastNlenData(int &nlen, char *&data, int proc, int tag = 0);

void PH_BcastComposite(DATA_COMPRESS dt_c, DATA_DECOMPRESS dt_d, int proc, int tag = 0);

void PH_BcastString(string &cs, int proc, int tag);

void PH_CollectString(string &cs, int proc, int tag);

void ServerRead(fstream &file, void *data, int size);

void PH_Trade(void *data, int size, int send_proc, int recv_proc, int tag = 0);

void PH_Trade(DataContainer *cdata, int send_proc, int recv_proc, int tag = 0);

void PH_Trade(DataContainer *cdata, int send_proc, int recv_proc, streamsize nlen, int tag = 0);

void PH_Trade(ActionKey *actkey, int send_proc, int recv_proc, int tag = 0);

void PH_Trade(ActionKey *actkey, int send_proc, int recv_proc, streamsize nlen, int tag);

template < typename T >
void PH_SmartBcastNewVersion(T *field, int numberOfElements, int processorIndex)
{
#ifdef PH_PARALLEL
    //! In order to avoid special situation.
    if (numberOfElements <= 0) return;
    int bufferSize = numberOfElements * sizeof(T);
    MPI_Bcast(field, bufferSize, MPI_CHAR, processorIndex, MPI_COMM_WORLD);
#endif

//! MPI_BCAST(buffer, count, datatype, root, comm)
//! IN/OUT buffer    the start address of communication buffer(optional data types).
//! IN count         number of data to be broadcast or received (int type).
//! IN datatype      types of data to be broadcast or received (handle).
//! IN root          the ID of broadcast data root process(int type).
//! IN comm          communication field (handle).
//int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
}

template < typename T >
void PH_SmartSend(T *field, int numberOfElements, int processorIndex, int tag)
{
#ifdef PH_PARALLEL
    //! In order to avoid special situation.
    if (numberOfElements <= 0) return;
    MPI_Send(field, numberOfElements * sizeof(T), MPI_CHAR, processorIndex, tag, MPI_COMM_WORLD);
#endif
}

template < typename T >
void PH_SmartRecv(T *field, int numberOfElements, int processorIndex, int tag)
{
#ifdef PH_PARALLEL
    //! In order to avoid special situation.
    if (numberOfElements <= 0) return;

    MPI_Status status;
    MPI_Recv(field, numberOfElements * sizeof(T), MPI_CHAR, processorIndex, tag, MPI_COMM_WORLD, &status);
#endif
}

template < typename T >
void FillGobalData(T localData, T *&globalData)
{
    int nTProcessor = PHMPI::GetNumberOfProcessor();
    int myid        = PHMPI::GetCurrentProcessorID();

    for (int iProc = 0; iProc < nTProcessor; ++ iProc)
    {
        int procID = iProc;
        T tempData;

        if (myid == procID)
        {
            tempData = localData;
        }

        PHMPI::PH_Bcast(&tempData, sizeof(T), procID);

        globalData[iProc] = tempData;
    }
}

template < typename T >
void FillGobalData(T *localData, int *localSize, T **&globalData)
{
    int nTProcessor = PHMPI::GetNumberOfProcessor();
    int myid        = PHMPI::GetCurrentProcessorID();

    for (int iProc = 0; iProc < nTProcessor; ++ iProc)
    {
        int procID = iProc;

        int nData = localSize[iProc];
        T *tempData = new T [nData];

        if (myid == procID)
        {
            for (int iData = 0; iData < nData; ++ iData)
            {
                tempData[iData] = localData[iData];
            }
        }

        PHMPI::PH_Bcast(tempData, nData * sizeof(T), procID);

        globalData[iProc] = tempData;
    }
}

template < typename T >
void FillGobalDataCollectBcast(T *localData, int *localSize, T **&globalData)
{
    int nTProcessor = PHMPI::GetNumberOfProcessor() + PHMPI::GetNumberOfGridProcessor();
    int myid        = PHMPI::GetCurrentProcessorID();
    int serverTmp      = PHMPI::GetServerProcessorID();
    int tag = 0;

    //! Step1: collection from each other processors by server.
    int totalDataSize = 0;
    int *stride = new int [nTProcessor];
    stride[0] = 0;
    for (int iProc = 0; iProc < nTProcessor; ++ iProc)
    {
        totalDataSize += localSize[iProc];

        if (iProc > 0) stride[iProc] = stride[iProc - 1] + localSize[iProc - 1];
    }

    T * globalDataTemp = new int [totalDataSize];

    PH_Gatherv(localData, localSize[myid], globalDataTemp, localSize, stride, serverTmp);

    delete [] stride;

    if (myid == serverTmp)
    {
        int dataCount = 0;
        for (int iProc = 0; iProc < nTProcessor; ++ iProc)
        {
            int nData = localSize[iProc];
            T *tempData = new T [nData];

            for (int jData = 0; jData < nData; ++ jData)
            {
                tempData[jData] = globalDataTemp[dataCount ++];
            }

            globalData[iProc] = tempData;
        }
    }

    delete [] globalDataTemp;

    //! Step2: broad casting to each other processors by server.
    DataContainer *cdata = new DataContainer();

    if (myid == serverTmp)
    {
        cdata->MoveToBegin();
        for (int iProc = 0; iProc < nTProcessor; ++ iProc)
        {
            T *temp = globalData[iProc];
            PHWrite(cdata, temp, localSize[iProc]);
        }
    }

    PHSPACE::PHMPI::PH_BcastByServerToClinet(cdata, 2, tag);

    if (myid != serverTmp)
    {
        cdata->MoveToBegin();
        for (int iProc = 0; iProc < nTProcessor; ++ iProc)
        {
            globalData[iProc] = new T [localSize[iProc]];
            PHRead(cdata, globalData[iProc], localSize[iProc]);
        }
    }

    delete cdata;    cdata = nullptr;
}

//! Broad cast a data with size of nData to every other processor.
template < typename T >
void PH_BcastNonBlocking(T *data, int nData, int bcastProcessor, int tag)
{
    using namespace PHMPI;

    vector< T * > sendDataBuffer;
    vector< T * > receivedDataBuffer;
    vector< PH_Request > requestContainer;

    int sendProcessor = bcastProcessor;
    int myid = GetCurrentProcessorID();
    int number_of_processor = GetNumberOfProcessor();

    if (myid == bcastProcessor)
    {
        sendDataBuffer.reserve(number_of_processor);
        requestContainer.reserve(number_of_processor);

        //! Step 1: Server sending.
        for (int iProc = 0; iProc < number_of_processor; ++ iProc)
        {
            if (iProc == bcastProcessor) continue;

            int receiveProcessor = iProc;

            //! Allocate the buffers for sending and receiving.
            T *sendBuffer = new T [nData];
            sendDataBuffer.push_back(sendBuffer);

            for (int iData = 0; iData < nData; ++ iData)
            {
                sendBuffer[iData] = data[iData];
            }

            //! Send the data to others.
            requestContainer.push_back(MPI_REQUEST_NULL);
            PHSPACE::PH_Send(sendBuffer, sizeof(T) * nData, receiveProcessor, tag, &requestContainer.back());
        }
    }
    else
    {
        T *receivedBuffer = new T[nData];
        receivedDataBuffer.push_back(receivedBuffer);

        //! Receive the data from server.
        requestContainer.push_back(MPI_REQUEST_NULL);
        PHSPACE::PH_Receive(receivedBuffer, sizeof(T) * nData, sendProcessor, tag, &requestContainer.back());
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translating.
    if (myid != bcastProcessor)
    {
        T *receivedBuffer = receivedDataBuffer[0];
        for (int iData = 0; iData < nData; ++ iData)
        {
            data[iData] = receivedBuffer[iData];
        }
    }

    //! Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        delete receivedDataBuffer[iDim];
    }
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        delete sendDataBuffer[iDim];
    }

    requestContainer.clear();
    receivedDataBuffer.clear();
    sendDataBuffer.clear();
}

}

bool IsNeedNonBlockingCommunication(const ActionKey *actkey);

template < typename T >
void PHRead(DataContainer *dataContainer, T *field, int numberOfElements)
{
    dataContainer->Read(field, numberOfElements * sizeof(T));
}

template < typename T >
void PHRead(DataContainer *dataContainer, T *field, uint_t numberOfElements)
{
    dataContainer->Read(field, numberOfElements * sizeof(T));
}

template < typename T >
void PHRead(DataContainer *dataContainer, vector<T> &field, int numberOfElements)
{
    dataContainer->Read(&field[0], numberOfElements * sizeof(T));
}

template < typename T >
void PHRead(DataContainer *dataContainer, T &value)
{
    dataContainer->Read(&value, sizeof(T));
}

template < typename T >
void PHWrite(DataContainer *dataContainer, T &value)
{
    dataContainer->Write(&value, sizeof(T));
}

template < typename T >
void PHWrite(DataContainer *dataContainer, T *field, int numberOfElements)
{
    dataContainer->Write(field, numberOfElements * sizeof(T));
}

template < typename T >
void PHWrite(DataContainer *dataContainer, T *field, uint_t numberOfElements)
{
    dataContainer->Write(field, numberOfElements * sizeof(T));
}

template < typename T >
void PHWrite(DataContainer *dataContainer, vector<T> &field, int numberOfElements)
{
    PHWrite(dataContainer, &field[0], numberOfElements);
}

template < typename T >
void PHRead(DataContainer *dataContainer, PHVector1D <T> &field, int numberOfElements)
{
    dataContainer->Read(&field[0], numberOfElements * sizeof(T));
}

template < typename T >
void PHRead(DataContainer *dataContainer, PHVector2D <T> &field2D, int numberOfElements)
{
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        PHVector1D <T> &field = field2D[iElement];

        int numberOfSubElements = 0;
        PHRead(dataContainer, numberOfSubElements);

        field.resize(numberOfSubElements);
        PHRead(dataContainer, field, numberOfSubElements);
    }
}

template < typename T >
void PHRead(DataContainer *dataContainer, PHVector2D <T> &field2D, uint_t numberOfElements)
{
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        PHVector1D <T> &field = field2D[iElement];

        int numberOfSubElements = 0;
        PHRead(dataContainer, numberOfSubElements);

        field.resize(numberOfSubElements);
        PHRead(dataContainer, field, numberOfSubElements);
    }
}

template < typename T >
void PHWrite(DataContainer *dataContainer, PHVector1D <T> &field, int numberOfElements)
{
    PHWrite(dataContainer, &field[0], numberOfElements);
}

template < typename T >
void PHWrite(DataContainer *dataContainer, PHVector2D <T> &field2D, int numberOfElements)
{
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        PHVector1D <T> &field = field2D[iElement];

        int numberOfSubElements = field.size();
        PHWrite(dataContainer, numberOfSubElements);
        PHWrite(dataContainer, field, numberOfSubElements);
    }
}

template < typename T >
void PHWrite(DataContainer *dataContainer, PHVector2D <T> &field2D, uint_t numberOfElements)
{
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        PHVector1D <T> &field = field2D[iElement];

        int numberOfSubElements = static_cast<int>(field.size());
        PHWrite(dataContainer, numberOfSubElements);
        PHWrite(dataContainer, field, numberOfSubElements);
    }
}

template < typename T >
void ZeroField(T *field, uint_t numberOfElements)
{
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        field[iElement] = T(0);
    }
}

template < typename T >
void ZeroField(T **field2D, uint_t numberOfEquations, uint_t numberOfElements)
{
    for (uint_t iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
    {
        ZeroField(field2D[iEquation], numberOfElements);
    }
}

template < typename T >
void ZeroField(PHVector1D <T> &field, uint_t numberOfElements)
{
    ZeroField(&field[0], numberOfElements);
}

void StreamToActionKey(ActionKey *actkey, ostringstream &oss);
void DataCharToString(DataContainer *cdata, string &str);

int GetSendRecvTag(ActionKey *actkey, int iZone);
int GetSendRecvTag(int tagKind, int index, int iZone);

void PH_AllreduceSepMode(int *sendingBuffer, int *receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op);
void PH_AllreduceSepMode(float *sendingBuffer, float *receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op);
void PH_AllreduceSepMode(double *sendingBuffer, double *receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op);

int PH_Reduce(int *sendbuf, int *recvbuf, int count, PH_Op op);
int PH_Reduce(float *sendbuf, float *recvbuf, int count, PH_Op op);
int PH_Reduce(double *sendbuf, double *recvbuf, int count, PH_Op op);

void PH_CompareMaxMin(int &local_data, int flag);
void PH_CompareMaxMin(RDouble &local_data, int flag);

int PH_BarrierSepMode();
}
