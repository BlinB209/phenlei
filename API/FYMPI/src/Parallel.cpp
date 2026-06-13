#include "Parallel.h"
#include "PHMpi.h"

namespace PHSPACE
{

int PH_Init()
{
    int err      = 0;
    int argc     = 0;
    char ***argv = 0;
#ifdef PH_PARALLEL
    err = MPI_Init(&argc, argv);
#endif
    return err;
}

int PH_Init(int &argc, char ***argv)
{
    int err = 0;
#ifdef PH_PARALLEL
    err = MPI_Init(&argc, argv);
#endif
    return err;
}

void PH_Rank(int &rank)
{
#ifdef PH_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
}

void PH_Size(int &size)
{
#ifdef PH_PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
}

void PH_Finalize()
{
#ifdef PH_PARALLEL
    MPI_Finalize();
#endif
}

void PH_GetProcessorName(string &proc_name)
{
#ifdef PH_PARALLEL
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int nlen = 0;

    MPI_Get_processor_name(processor_name, &nlen);
    proc_name = processor_name;
#endif
}

void PH_Send(void *data, int size, int to, int tag)
{
#ifdef PH_PARALLEL
    //! In order to avoid special situation.
    if (size <= 0) return;
    MPI_Send(data, size, MPI_CHAR, to, tag, MPI_COMM_WORLD);
#endif
}

void PH_Receive(void *data, int size, int from, int tag)
{
#ifdef PH_PARALLEL
    //! In order to avoid special situation.
    if (size <= 0) return;

    MPI_Status status;
    MPI_Recv(data, size, MPI_CHAR, from, tag, MPI_COMM_WORLD, &status);
#endif
}

int PH_Receive(void *data, int size, int tag)
{
    int err_code = 0;

#ifdef PH_PARALLEL
    //! In order to avoid special situation.
    if (size <= 0) return -1;

    MPI_Status status;
    MPI_Recv(data, size, MPI_CHAR, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
    err_code = status.MPI_SOURCE;
#endif

    return err_code;
}

void PH_SendString(string &cs, int proc, int tag)
{
    int nlen = static_cast<int>(cs.length());

    PH_Send(&nlen, sizeof(int), proc, tag);

    char *data = new char[nlen+1];

    cs.copy(data, nlen);
    data[nlen] = '\0';

    PH_Send(data, nlen+1, proc, tag);

    delete [] data;
}

void PH_ReceiveString(string &cs, int proc, int tag)
{
    int nlen = 0;
    PH_Receive(&nlen, sizeof(int), proc, tag);

    char *data = new char[nlen+1];

    PH_Receive(data, nlen+1, proc, tag);

    cs = data;

    delete [] data;
}

//! Non-blocking send/receive function:
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function :   int PH_Send(void *data, int size, int to, int tag, PH_Request *request)                                  +
// Author   :   Xu QingXin,  2013.02.16,   Modified 2013.02.16.                                                          +
// Purpose  :   Non-blocking standard send function,package Non-blocking standard send process function of MPICH.        +
// Description: data, input type parameter, optional type, address of the send buffer.                                   +
//              size, input type parameter, int type,size of the send buffer in bytes.                                   +
//              to  , input type parameter, int type,the send target process number.                                     +
//              tag , input type parameter, int type,information tag.                                                    +
//              comm, input optional parameter, MPI_Comm type, communication field handle, default is MPI_COMM_WORLD.    +
//              request, output type parameter, PH_Request type pointer, communication object handle.                    +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int PH_Send(void *data, int size, int to, int tag, PH_Request *request)
{
    int errcode = 0;

    if (size <= 0)
    {
        errcode = -1;
        return errcode;
    }

#ifdef PH_PARALLEL
    errcode = MPI_Isend(data, size, MPI_CHAR, to, tag, MPI_COMM_WORLD, request);
#endif

    return errcode;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function    :   int PH_Receive(void *data, int size, int from, int tag, PH_Request *request)                               +
// Author      :   Xu QingXin,  2013.02.07, Modified 2013.02.07.                                                              +
// Purpose     :   Non-blocking standard receive function, package Non-blocking standard receive process function of MPICH.   +
// Description :   data, input type parameter, optional type, address of the receive buffer.                                  +
//                 size, input type parameter, int type, size of the receive buffer in bytes.                                 +
//                 from, input type parameter, int type, the receive source process number.                                   +
//                 tag , input type parameter, int type,information tag.                                                      +
//                 request, output type parameter, PH_Request type pointer, communication object handle.                      +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int PH_Receive(void *data, int size, int from, int tag, PH_Request *request)
{
    int errcode = 0;

    if (size <= 0)
    {
        errcode = -1;
        return errcode;
    }

#ifdef PH_PARALLEL
    errcode = MPI_Irecv(data, size, MPI_CHAR, from, tag, MPI_COMM_WORLD, request);
#endif
    
    return errcode;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function    :   int PH_Wait(PH_Request *request)                                                                                                   +
// Author      :   Xu QingXin, 2013.02.16, Modified 2013.02.16.                                                                                       +
// Purpose     :   Non-blocking completed communication operation function, package Non-blocking completed communication process function of MPICH.   +
// Description :request, input/output type parameter, PH_Request type pointer, communication object handle.                                           +
//               status, output type parameter, PH_Status type pointer, state object,default is MPI_STATUS_IGNORE.                                    +
//                count, input type parameter, int type, number of communication objects to be completed.                                             +
//    array_of_requests, input/output type parameter, PH_Request type one dimension array, store count communication objects handle.                  +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int PH_Wait(PH_Request *request)
{
    int errcode = 0;
    
#ifdef PH_PARALLEL
    errcode = MPI_Wait(request, MPI_STATUS_IGNORE);
#endif
    
    return errcode;
}

int PH_Wait(int count, PH_Request *array_of_requests)
{
    int errcode = 0;
 
#ifdef PH_PARALLEL
    errcode = MPI_Waitall(count, array_of_requests, MPI_STATUSES_IGNORE);
#endif
    
    return errcode;
}

//! Global reduction function:
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function    :   int PH_Reduce(void *sendbuf, void *recvbuf, int count, PH_Op op)                    +
//                 int PH_Reduce(void *sendbuf, void *recvbuf, int count, PH_Op op, int root)          +
// Author      :   Xu QingXin, 2013.04.02, Modified 2013.04.02.                                        +
// Purpose     :   Global reduction function, package Global reduction operation function of MPICH.    +
// Description : sendbuf, input type parameter, optional type, address of the send buffer.             +
//               recvbuf, input type parameter, optional type, address of the receive buffer.          +
//               count  , input type parameter, int type,number of communication objects to be reduced.+
//               op     , input type parameter, MPI_Op type, reduction operation handle.               +
//               root   , input type parameter, int type,root process number.                          +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int PH_AllReduce(void *sendbuf, void *recvbuf, int count, PH_Op op)
{
    int errcode = 0;

#ifdef PH_PARALLEL
    errcode = MPI_Allreduce(sendbuf, recvbuf, count, MPI_CHAR, op, MPI_COMM_WORLD);
#endif

    return errcode;
}

int PH_AllReduce(double *sendbuf, double *recvbuf, int count, PH_Op op)
{
    int errcode = 0;
    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        for (int iData = 0; iData < count; ++ iData)
        {
            recvbuf[iData] = sendbuf[iData];
        }
        return 1;
    }

#ifdef PH_PARALLEL
    errcode = MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, MPI_COMM_WORLD);
#endif

    return errcode;
}

int PH_AllReduce(float *sendbuf, float *recvbuf, int count, PH_Op op)
{
    int errcode = 0;

#ifdef PH_PARALLEL
    errcode = MPI_Allreduce(sendbuf, recvbuf, count, MPI_FLOAT, op, MPI_COMM_WORLD);
#endif

    return errcode;
}

int PH_AllReduce(int *sendbuf, int *recvbuf, int count, PH_Op op)
{
    int errcode = 0;
    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        for (int iData = 0; iData < count; ++ iData)
        {
            recvbuf[iData] = sendbuf[iData];
        }
        return 1;
    }

#ifdef PH_PARALLEL
    errcode = MPI_Allreduce(sendbuf, recvbuf, count, MPI_INT, op, MPI_COMM_WORLD);
#endif

    return errcode;
}

int PH_AllGather(int *sendbuf, int sendcount, int *recvbuf, int recvcount)
{
    int errcode = 0;
    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        for (int iData = 0; iData < sendcount; ++ iData)
        {
            recvbuf[iData] = sendbuf[iData];
        }
        return 1;
    }

#ifdef PH_PARALLEL
    errcode = MPI_Allgather(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, MPI_INT, MPI_COMM_WORLD);
#endif

    return errcode;
}

// int PH_Reduce(void *sendbuf, void *recvbuf, int count, PH_Op op, int root)
// {
//      int errcode = 0;
// 
// #ifdef PH_PARALLEL
//      errcode = MPI_Reduce(sendbuf, recvbuf, count, MPI_CHAR, op, root, MPI_COMM_WORLD);
// #endif
// 
//      return errcode;
// }

int PH_Barrier()
{
    int errorCode = 0;
#ifdef PH_PARALLEL
    errorCode = MPI_Barrier(MPI_COMM_WORLD);
#endif
    return errorCode;
}

int PH_Gatherv(void *sendBuffer, int sendCount, void *receiveBuffer, int *receiveCount, int *displs, int root)
{
    int errcode = 0;

#ifdef PH_PARALLEL
    errcode = MPI_Gatherv(sendBuffer, sendCount, MPI_INT, receiveBuffer, receiveCount, displs, MPI_INT, root, MPI_COMM_WORLD);
#endif

    return errcode;
}

void send(DataChar *dataChar, int to, int tag)
{
    int nlen = static_cast<int>(dataChar->Size());

    //! It is necessary to judge whether the data length is zero.
    if (nlen <= 0) return;

    PH_Send(dataChar->Begin(), nlen, to, tag);
}

void send(DataChar *dataChar, int to, vector<PH_Request> &requestContainer, int tag)
{
    int nlen = static_cast<int>(dataChar->Size());

    //! It is necessary to judge whether the data length is zero.
    if (nlen <= 0) return;

    requestContainer.push_back(MPI_REQUEST_NULL);

    PH_Send(dataChar->Begin(), nlen, to, tag, &requestContainer.back());
}

void receive(DataChar *dataChar, int from, int tag)
{
    int nlen = static_cast<int>(dataChar->Size());

    //! It is necessary to judge whether the data length is zero.
    if (nlen <= 0) return;

    PH_Receive(dataChar->Begin(), nlen, from, tag);
}

void receive(DataChar *dataChar, int from, vector<PH_Request> &requestContainer, int tag)
{
    int nlen = static_cast<int>(dataChar->Size());

    //! It is necessary to judge whether the data length is zero.
    if (nlen <= 0) return;

    requestContainer.push_back(MPI_REQUEST_NULL);

    PH_Receive(dataChar->Begin(), nlen, from, tag, &requestContainer.back());
}

void send(DataContainer *dataContainer, int to, int tag)
{
    CharVecSizeType nlen = dataContainer->Size();

    PH_Send(&nlen, sizeof(streamsize), to, tag);

    //! It is necessary to judge whether the data length is zero.
    if (nlen <= 0) return;

    for (ContainerSizeType i = 0; i < dataContainer->ElementSize(); ++ i)
    {
        send(dataContainer->GetIter(i), to, tag);
    }
}

void send(DataContainer *dataContainer, int to, streamsize nlen, int tag)
{
    //! It is necessary to judge whether the data length is zero.
    if (nlen <= 0) return;

    for (ContainerSizeType i = 0; i < dataContainer->ElementSize(); ++ i)
    {
        send(dataContainer->GetIter(i), to, tag);
    }
}

void send(DataContainer *dataContainer, int to, streamsize nlen, vector<PH_Request> &requestContainer, int tag)
{
    //! It is necessary to judge whether the data length is zero.
    if (nlen <= 0) return;

    for (ContainerSizeType i = 0; i < dataContainer->ElementSize(); ++ i)
    {
        send(dataContainer->GetIter(i), to, requestContainer, tag);
    }
}

void receive(DataContainer *dataContainer, int from, int tag)
{
    streamsize nlen = 0;
    PH_Receive(&nlen, sizeof(streamsize), from, tag);

    if (nlen <= 0) return;

    dataContainer->SecureAbsoluteSpace(static_cast< CharVecSizeType >(nlen));

    for (ContainerSizeType i = 0; i < dataContainer->ElementSize(); ++ i)
    {
        receive(dataContainer->GetIter(i), from, tag);
    }
}

void receive(DataContainer *dataContainer, int from, CharVecSizeType nlen, int tag)
{
    if (nlen <= 0) return;

    dataContainer->SecureAbsoluteSpace(nlen);

    for (ContainerSizeType i = 0; i < dataContainer->ElementSize(); ++ i)
    {
        receive(dataContainer->GetIter(i), from, tag);
    }
}

void receive(DataContainer *dataContainer, int from, CharVecSizeType nlen, vector<PH_Request> &requestContainer, int tag)
{
    if (nlen <= 0) return;

    dataContainer->SecureAbsoluteSpace(nlen);

    for (ContainerSizeType i = 0; i < dataContainer->ElementSize(); ++ i)
    {
        receive(dataContainer->GetIter(i), from, requestContainer, tag);
    }
}

}