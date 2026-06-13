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
//! @file      Parallel.h
//! @brief     Explain this file briefly.
//! @author    Bell, He Xin, Xu Qingxin.

#pragma once
#include "DataContainer.h"

#ifdef PH_PARALLEL
    #include "mpi.h"

    typedef MPI_Request PH_Request;
    typedef MPI_Op      PH_Op;

    #define PH_REQUEST_NULL MPI_REQUEST_NULL
    #define PH_MAX          MPI_MAX
    #define PH_MIN          MPI_MIN
    #define PH_SUM          MPI_SUM
#else
    typedef int PH_Request;
    typedef int PH_Op;

    #define PH_REQUEST_NULL 0
    #define PH_MAX          0
    #define PH_MIN          0
    #define PH_SUM          0
#endif

using namespace std;

namespace PHSPACE
{
int  PH_Init();
int  PH_Init(int &argc, char ***argv);
void PH_Rank(int &rank);
void PH_Size(int &size);
void PH_Finalize();

void PH_GetProcessorName(string &proc_name);

//! Non-blocking standard send process of MPICH.
//! @param[in] data    address of the send buffer.
//! @param[in] size    size of the send buffer.
//! @param[in] to      the send target process number.
//! @param[in] tag     information tag.
void PH_Send(void *data, int size, int to, int tag);

//! Non-blocking standard receive process of MPICH.
//! @param[in] data    address of the receive buffer.
//! @param[in] size    size of the receive buffer.
//! @param[in] to      the receive source process number.
//! @param[in] tag     information tag.
void PH_Receive(void *data, int size, int from, int tag);
int  PH_Receive(void *data, int size, int tag);
//!
void PH_SendString(string &cs, int proc, int tag);
void PH_ReceiveString(string &cs, int proc, int tag);

//! Non-blocking standard send process of MPICH.
//! @param[in] data    address of the send buffer.
//! @param[in] size    size of the send buffer.
//! @param[in] to      the send target process number.
//! @param[in] tag     information tag.
//! @param[in] request PH_Request pointer, communication object handle.
int PH_Send(void *data, int size, int to, int tag, PH_Request *request);

//! Non-blocking standard receive process of MPICH.
//! @param[in] data    address of the receive buffer.
//! @param[in] size    size of the receive buffer.
//! @param[in] from    the receive source process number.
//! @param[in] tag     information tag.
//! @param[in] request PH_Request pointer, communication object handle.
int PH_Receive(void *data, int size, int from, int tag, PH_Request *request);

//! Non-blocking completed communication process of MPICH.
//! @param[in] request    PH_Request pointer, communication object handle.
//! @param[in] count      number of communication objects to be completed.
int PH_Wait(PH_Request *request);
int PH_Wait(int count, PH_Request *array_of_requests);

//! Global reduction operation of MPICH.
//! @param[in] sendbuf    address of the send buffer.
//! @param[in] recvbuf    address of the receive buffer.
//! @param[in] count      number of communication objects to be reduced.
//! @param[in] op         reduction operation handle.
int PH_AllReduce(void   *sendbuf, void   *recvbuf, int count, PH_Op op);
int PH_AllReduce(double *sendbuf, double *recvbuf, int count, PH_Op op);
int PH_AllReduce(float  *sendbuf, float  *recvbuf, int count, PH_Op op);
int PH_AllReduce(int    *sendbuf, int    *recvbuf, int count, PH_Op op);

int PH_AllGather(int *sendbuf, int sendcount, int *recvbuf, int recvcount);

//int PH_Reduce(void *sendbuf, void *recvbuf, int count, PH_Op op, int root);
int PH_Barrier();

//!
int PH_Gatherv(void *sendBuffer, int sendCount, void *receiveBuffer, int *receiveCount, int *displs, int root);

//!
void send(DataChar *dataChar, int to, int tag);
void send(DataChar *dataChar, int to, vector<PH_Request> &requestContainer, int tag);
void receive(DataChar *dataChar, int from, int tag);
void receive(DataChar *dataChar, int from, vector<PH_Request> &requestContainer, int tag);

//!
void send(DataContainer *dataContainer, int to, int tag);
void send(DataContainer *dataContainer, int to, streamsize nlen, int tag);
void send(DataContainer *dataContainer, int to, streamsize nlen, vector<PH_Request> &requestContainer, int tag);

//!
void receive(DataContainer *dataContainer, int from, int tag);
void receive(DataContainer *dataContainer, int from, CharVecSizeType nlen, int tag);
void receive(DataContainer *dataContainer, int from, CharVecSizeType nlen, vector<PH_Request> &requestContainer, int tag);
}
