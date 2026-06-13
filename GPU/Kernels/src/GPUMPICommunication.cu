#include "GPUDeviceControl.h"
#include "GPUKernelTestPart2.h"
#include "GPUMPICommunication.h"
#include "OutputDebug.h"
#include "cudaErrorHandle.h"
#include "mpi.h"

using namespace std;
namespace GPUKernels
{
    void CompressInterfaceTurbDataSmallGPUBufferOverlap(const int nameID, const int nTotal, const int neqn,
                                                        const int nIFace, const int nTotalNgbZones,
                                                        InterfaceInfo *iinfo, RFloat **GPUDataSendSmallBuffer)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        for (int ngbID = nTotalNgbZones - 1; ngbID >= 0; ngbID--)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = offsetNgbLBTurb[ngbID * nameTurbBuffer->size() + nameID];

            int gridSize  = 1;
            int blockSize = 1;
            int loopLen   = nIFaceOfNeighbor;

            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompressSmallBuf, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompressSmallBuf, 0, gridSize, blockSize);
#endif

            //! get initial address of variable in GPUSendBufferTurbNDLE_API_ERR(
            int offsetVar = offsetBufferTurb[nameID];
            //! get initial address of variable in GPUSendBufferTurb
            int dimData = dimBufferTurb[nameID];
            //! set trafer data size
            size_t sizeBuf = dimData * nIFaceOfNeighbor * sizeof(double);
            //! set HostSendBuffer
            double *HostSendBuffer = HostMPIDataSendTurb + offsetNgb + offsetNgbLengthBuffer;
            GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                offsetFaces, d_faceIndexForSend, GPUSendBufferTurb,
                                                                GPUDataSendSmallBuffer[ngbID]);
            HANDLE_API_ERR(cudaEventRecord(compressDoneTurb, 0));
            //! HANDLE_API_ERR(cudaStreamWaitEvent(dataStreamTurb, compressDoneTurb, 0));
            HANDLE_API_ERR(cudaStreamWaitEvent(dataStreamTurb, compressDoneTurb, 0));
            HANDLE_API_ERR(cudaMemcpyAsync(HostSendBuffer, GPUDataSendSmallBuffer[ngbID], sizeBuf,
                                           cudaMemcpyDeviceToHost, dataStreamTurb));
            //! HANDLE_API_ERR(cudaMemcpy(HostSendBuffer, GPUDataSendSmallBuffer[ngbID],
            //! sizeBuf, cudaMemcpyDeviceToHost));
        }
    }

    void CompressInterfaceNSDataSmallGPUBufferOverlap(const int nameID, const int nTotal, const int neqn,
                                                      const int nIFace, const int nTotalNgbZones, InterfaceInfo *iinfo,
                                                      RFloat **GPUDataSendSmallBuffer)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        for (int ngbID = nTotalNgbZones - 1; ngbID >= 0; ngbID--)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = offsetNgbLBNS[ngbID * nameNSBuffer->size() + nameID];

            int gridSize  = 1;
            int blockSize = 1;
            int loopLen   = nIFaceOfNeighbor;

            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompressSmallBuf, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompressSmallBuf, 0, gridSize, blockSize);
#endif

            //! get initial address of variable in GPUSendBufferNSNDLE_API_ERR(
            int offsetVar = offsetBufferNS[nameID];
            //! get initial address of variable in GPUSendBufferNS
            int dimData = dimBufferNS[nameID];
            //! set trafer data size
            size_t sizeBuf = dimData * nIFaceOfNeighbor * sizeof(double);
            //! set HostSendBuffer
            double *HostSendBuffer = HostMPIDataSendNS + offsetNgb + offsetNgbLengthBuffer;
            GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                offsetFaces, d_faceIndexForSend, GPUSendBufferNS,
                                                                GPUDataSendSmallBuffer[ngbID]);
            HANDLE_API_ERR(cudaEventRecord(compressDoneNS, 0));
            //! HANDLE_API_ERR(cudaStreamWaitEvent(dataStreamNS, compressDoneNS, 0));
            HANDLE_API_ERR(cudaStreamWaitEvent(dataStreamNS, compressDoneNS, 0));
            HANDLE_API_ERR(cudaMemcpyAsync(HostSendBuffer, GPUDataSendSmallBuffer[ngbID], sizeBuf,
                                           cudaMemcpyDeviceToHost, dataStreamNS));
            //! HANDLE_API_ERR(cudaMemcpy(HostSendBuffer, GPUDataSendSmallBuffer[ngbID],
            //! sizeBuf, cudaMemcpyDeviceToHost));
        }
    }

    void CallCompressInterfaceDataSmallGPUBufferOverlap(const int nTotal, const int neqn, const int nIFace,
                                                        const string &name, InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        int     nameNSID, nameTurbID;
        RFloat *deviceFieldVar;
        double *deviceSendBuffer;
        int     offsetBuffer;
        int     nTotalNgbZones = numNgbZones[0];
        if ("q" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_q);
        }
        else if ("t" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_t);
        }
        else if ("turb::q" == name)
        {
            nameTurbID = GetIndexBufferTurb(name);
            CompressInterfaceTurbDataSmallGPUBufferOverlap(nameTurbID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                           GPUDataSendTurb_q);
        }
        else if ("dqdx" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_dqdx);
        }
        else if ("dqdy" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_dqdy);
        }
        else if ("dqdz" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_dqdz);
        }
        else if ("limit" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_limit);
        }
        else if ("dtdx" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_dtdx);
        }
        else if ("dtdy" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_dtdy);
        }
        else if ("dtdz" == name)
        {
            nameNSID = GetIndexBufferNS(name);
            CompressInterfaceNSDataSmallGPUBufferOverlap(nameNSID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                         GPUDataSendNS_dtdz);
        }
        else if ("turb::dqdx" == name)
        {
            nameTurbID = GetIndexBufferTurb(name);
            CompressInterfaceTurbDataSmallGPUBufferOverlap(nameTurbID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                           GPUDataSendTurb_dqdx);
        }
        else if ("turb::dqdy" == name)
        {
            nameTurbID = GetIndexBufferTurb(name);
            CompressInterfaceTurbDataSmallGPUBufferOverlap(nameTurbID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                           GPUDataSendTurb_dqdy);
        }
        else if ("turb::dqdz" == name)
        {
            nameTurbID = GetIndexBufferTurb(name);
            CompressInterfaceTurbDataSmallGPUBufferOverlap(nameTurbID, nTotal, neqn, nIFace, nTotalNgbZones, iinfo,
                                                           GPUDataSendTurb_dqdz);
        }
        else
        {
            cout << "Error:" << name << " is not included in CallGPUUploadInterfaceValue" << endl;
            exit(1);
        }
    }
    void DecompressInterfaceNSDataGPUSmallBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
#ifdef MPIOVERLAP
        if (GPUMPIDataRecvTurb)
        {
            MPI_Waitall(nTotalRequestTurb, requestContainerTurbRemains, MPI_STATUSES_IGNORE);
            int nTotalNgbZonesTurb = numNgbZones[0];
            int nIFaceTurb         = iinfo->GetNIFace();

            for (int ngbID = nTotalNgbZonesTurb - 1; ngbID >= 0; ngbID--)
            {
                //! get neighbor zone's global index (jZone)
                int ngbZoneID = recvGlobalZone[ngbID];
                //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
                //! ngbID=iNeighbor
                int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
                //! get total faces of jZone's interface
                int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                //! get initial address of jZone in d_send
                int offsetFaces = offsetFaceIndex[ngbID];
                //! get inital address of jZone in GPUMPIDataSendTurb
                int offsetNgb = offsetMPIDataTurb[ngbID];
                //! get offset of each variables in buffer of jZone
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompressSmallBuf, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompressSmallBuf, 0, gridSize, blockSize);
#endif

                for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                {
                    //! get name of varialbe
                    const string &name = (*nameTurbBuffer)[nameID];
                    //! get initial address of variable in GPUSendBufferTurb
                    int offsetVar = offsetBufferTurb[nameID];
                    //! get initial address of variable in GPUSendBufferTurb
                    int dimData = dimBufferTurb[nameID];
                    //! set trafer data size
                    size_t sizeBuf = dimData * nIFaceOfNeighbor * sizeof(double);
                    //! set HostRecvBuffer
                    double *HostRecvBuffer = HostMPIDataRecvTurb + offsetNgb + offsetNgbLengthBuffer;
                    if (name == "turb::dqdx")
                    {
#ifdef HDCOVERLAP
                        HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvTurb_dqdx[ngbID], HostRecvBuffer, sizeBuf,
                                                       cudaMemcpyHostToDevice, dataStreamTurb));
                        GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamTurb>>>(
                            offsetVar, dimData, nIFaceTurb, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                            GPUDataRecvTurb_dqdx[ngbID], GPURecvBufferTurb);
#else
                        //! HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvTurb_dqdx[ngbID],
                        //! HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                        HANDLE_API_ERR(
                            cudaMemcpy(GPUDataRecvTurb_dqdx[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                        GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                            offsetVar, dimData, nIFaceTurb, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                            GPUDataRecvTurb_dqdx[ngbID], GPURecvBufferTurb);
#endif
                    }
                    else if (name == "turb::dqdy")
                    {
#ifdef HDCOVERLAP
                        HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvTurb_dqdy[ngbID], HostRecvBuffer, sizeBuf,
                                                       cudaMemcpyHostToDevice, dataStreamTurb));
                        GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamTurb>>>(
                            offsetVar, dimData, nIFaceTurb, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                            GPUDataRecvTurb_dqdy[ngbID], GPURecvBufferTurb);
#else
                        HANDLE_API_ERR(
                            cudaMemcpy(GPUDataRecvTurb_dqdy[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                        GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                            offsetVar, dimData, nIFaceTurb, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                            GPUDataRecvTurb_dqdy[ngbID], GPURecvBufferTurb);
#endif
                    }
                    else if (name == "turb::dqdz")
                    {
#ifdef HDCOVERLAP
                        HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvTurb_dqdz[ngbID], HostRecvBuffer, sizeBuf,
                                                       cudaMemcpyHostToDevice, dataStreamTurb));
                        GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamTurb>>>(
                            offsetVar, dimData, nIFaceTurb, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                            GPUDataRecvTurb_dqdz[ngbID], GPURecvBufferTurb);
#else
                        HANDLE_API_ERR(
                            cudaMemcpy(GPUDataRecvTurb_dqdz[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                        GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                            offsetVar, dimData, nIFaceTurb, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                            GPUDataRecvTurb_dqdz[ngbID], GPURecvBufferTurb);
#endif
                    }
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
            }
        }
#endif //! end ifdef MPIOVERLAP                                                                                        \
    //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! Get total interface size of iZone
        int nIFace = iinfo->GetNIFace();
        //! for (int ngbID = nTotalNgbZones-1; ngbID >= 0; ngbID--){
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataDecompressSmallBuf, loopLen, 0, gridSize, blockSize, blockSizeLimit);

            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameNSBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferNS
                int offsetVar = offsetBufferNS[nameID];
                //! get initial address of variable in GPUSendBufferNS
                int dimData = dimBufferNS[nameID];
                //! set trafer data size
                size_t sizeBuf = dimData * nIFaceOfNeighbor * sizeof(double);
                //! set HostSendBuffer
                double *HostRecvBuffer = HostMPIDataRecvNS + offsetNgb + offsetNgbLengthBuffer;
                if (name == "q")
                {
                    HANDLE_API_ERR(cudaMemcpy(GPUDataRecvNS_q[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_q[ngbID], GPURecvBufferNS);
                }
                else if (name == "t")
                {
                    HANDLE_API_ERR(cudaMemcpy(GPUDataRecvNS_t[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_t[ngbID], GPURecvBufferNS);
                }
#ifndef MPIOVERLAP //! Afterwards operations are not necessary instantly.
                else if (name == "dqdx")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dqdx[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_dqdx[ngbID], GPURecvBufferNS);
                }
                else if (name == "dqdy")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dqdy[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_dqdy[ngbID], GPURecvBufferNS);
                }
                else if (name == "dqdz")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dqdz[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_dqdz[ngbID], GPURecvBufferNS);
                }
                else if (name == "limit")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_limit[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_limit[ngbID], GPURecvBufferNS);
                }
                else if (name == "dtdx")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dtdx[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_dtdx[ngbID], GPURecvBufferNS);
                }
                else if (name == "dtdy")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dtdy[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_dtdy[ngbID], GPURecvBufferNS);
                }
                else if (name == "dtdz")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dtdz[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvNS_dtdz[ngbID], GPURecvBufferNS);
                }
#endif //! end ifdef MPIOVERLAP
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

#ifdef TESTNEWORGMPI
        //! Test HostMPIDataRecvNS and HostRecvLargeBuffer
        double *HostRecvBufferNSAWARE = new double[lengthBufferNS];
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForRecv      = iinfo->GetFaceIndexForRecv(iNeighbor);
            int  offsetNgbLengthBuffer = 0;
            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                int dimData = dimBufferNS[nameID];
                TestGPUMPIDataDecompressNS(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                           offsetNgbLengthBuffer, faceIndexForRecv, HostMPIDataRecvNS,
                                           HostRecvBufferNS);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }
#endif
    }

    __global__ void GPUMPIDataDecompressSmallBuf(const int offsetVar, const int neqn, const int nIFaceTotal,
                                                 const int nIFaceOfNeighbor, const int offsetFaceIndex,
                                                 const int *faceIndexForRecv, const double *GPUMPIDataRecvSmall,
                                                 double *GPURecvBufferLarge)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int faceID;
        int equationID;
        int ngbZoneFaceID;
        int largeStorageID;
        int faceIndexLargeStore;
        for (faceID = threadID; faceID < nIFaceOfNeighbor; faceID += blockDim.x * gridDim.x)
        {
            faceIndexLargeStore = faceIndexForRecv[offsetFaceIndex + faceID];

            for (equationID = 0; equationID < neqn; equationID++)
            {
                ngbZoneFaceID                      = equationID * nIFaceOfNeighbor + faceID;
                largeStorageID                     = offsetVar + equationID * nIFaceTotal + faceIndexLargeStore;
                GPURecvBufferLarge[largeStorageID] = GPUMPIDataRecvSmall[ngbZoneFaceID];
            }
        }
    }

    void DecompressInterfaceTurbDataGPUSmallBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

#ifdef MPIOVERLAP
        MPI_Waitall(nTotalRequestNS, requestContainerNSRemains, MPI_STATUSES_IGNORE);
        int nTotalNgbZonesNS = numNgbZones[0];
        int nIFaceNS         = iinfo->GetNIFace();

        for (int ngbID = 0; ngbID < nTotalNgbZonesNS; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataDecompressSmallBuf, loopLen, 0, gridSize, blockSize, blockSizeLimit);

            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameNSBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferNS
                int offsetVar = offsetBufferNS[nameID];
                //! get initial address of variable in GPUSendBufferNS
                int dimData = dimBufferNS[nameID];
                //! set trafer data size
                size_t sizeBuf = dimData * nIFaceOfNeighbor * sizeof(double);
                //! set HostSendBuffer
                double *HostRecvBuffer = HostMPIDataRecvNS + offsetNgb + offsetNgbLengthBuffer;
                if (name == "dqdx")
                {
#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvNS_dqdx[ngbID], HostRecvBuffer, sizeBuf,
                                                   cudaMemcpyHostToDevice, dataStreamNS));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamNS>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dqdx[ngbID], GPURecvBufferNS);
#else
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dqdx[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dqdx[ngbID], GPURecvBufferNS);
#endif
                }
                else if (name == "dqdy")
                {
#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvNS_dqdy[ngbID], HostRecvBuffer, sizeBuf,
                                                   cudaMemcpyHostToDevice, dataStreamNS));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamNS>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dqdy[ngbID], GPURecvBufferNS);
#else
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dqdy[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dqdy[ngbID], GPURecvBufferNS);
#endif
                }
                else if (name == "dqdz")
                {
#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvNS_dqdz[ngbID], HostRecvBuffer, sizeBuf,
                                                   cudaMemcpyHostToDevice, dataStreamNS));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamNS>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dqdz[ngbID], GPURecvBufferNS);
#else
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dqdz[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dqdz[ngbID], GPURecvBufferNS);

#endif
                }
                else if (name == "limit")
                {
#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvNS_limit[ngbID], HostRecvBuffer, sizeBuf,
                                                   cudaMemcpyHostToDevice, dataStreamNS));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamNS>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_limit[ngbID], GPURecvBufferNS);
#else
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_limit[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_limit[ngbID], GPURecvBufferNS);
#endif
                }
                else if (name == "dtdx")
                {
#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvNS_dtdx[ngbID], HostRecvBuffer, sizeBuf,
                                                   cudaMemcpyHostToDevice, dataStreamNS));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamNS>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dtdx[ngbID], GPURecvBufferNS);
#else
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dtdx[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dtdx[ngbID], GPURecvBufferNS);

#endif
                }
                else if (name == "dtdy")
                {
#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvNS_dtdy[ngbID], HostRecvBuffer, sizeBuf,
                                                   cudaMemcpyHostToDevice, dataStreamNS));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamNS>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dtdy[ngbID], GPURecvBufferNS);
#else
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dtdy[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dtdy[ngbID], GPURecvBufferNS);
#endif
                }
                else if (name == "dtdz")
                {
#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaMemcpyAsync(GPUDataRecvNS_dtdz[ngbID], HostRecvBuffer, sizeBuf,
                                                   cudaMemcpyHostToDevice, dataStreamNS));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize, 0, dataStreamNS>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dtdz[ngbID], GPURecvBufferNS);
#else
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvNS_dtdz[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFaceNS, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvNS_dtdz[ngbID], GPURecvBufferNS);
#endif
                }
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }
#endif //! end ifdef MPIOVERLAP

        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! Get total interface size of iZone
        int nIFace = iinfo->GetNIFace();
        for (int ngbID = nTotalNgbZones - 1; ngbID >= 0; ngbID--)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompressSmallBuf, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompressSmallBuf, 0, gridSize, blockSize);
#endif

            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameTurbBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferTurb
                int offsetVar = offsetBufferTurb[nameID];
                //! get initial address of variable in GPUSendBufferTurb
                int dimData = dimBufferTurb[nameID];
                //! set trafer data size
                size_t sizeBuf = dimData * nIFaceOfNeighbor * sizeof(double);
                //! set HostRecvBuffer
                double *HostRecvBuffer = HostMPIDataRecvTurb + offsetNgb + offsetNgbLengthBuffer;
                if (name == "turb::q")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvTurb_q[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                          offsetFaces, d_faceIndexForRecv,
                                                                          GPUDataRecvTurb_q[ngbID], GPURecvBufferTurb);
                }
#ifndef MPIOVERLAP
                else if (name == "turb::dqdx")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvTurb_dqdx[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFace, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvTurb_dqdx[ngbID], GPURecvBufferTurb);
                }
                else if (name == "turb::dqdy")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvTurb_dqdy[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFace, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvTurb_dqdy[ngbID], GPURecvBufferTurb);
                }
                else if (name == "turb::dqdz")
                {
                    HANDLE_API_ERR(
                        cudaMemcpy(GPUDataRecvTurb_dqdz[ngbID], HostRecvBuffer, sizeBuf, cudaMemcpyHostToDevice));
                    GPUMPIDataDecompressSmallBuf<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFace, nIFaceOfNeighbor, offsetFaces, d_faceIndexForRecv,
                        GPUDataRecvTurb_dqdz[ngbID], GPURecvBufferTurb);
                }
#endif //! end ifdef MPIOVERLAP
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

#ifdef TESTNEWORGMPI
        //! Test HostMPIDataRecvTurb and HostRecvLargeBuffer
        size_t sizeAWARETest = sizeof(double) * lengthBufferTurb;
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForRecv      = iinfo->GetFaceIndexForRecv(iNeighbor);
            int  offsetNgbLengthBuffer = 0;
            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                int dimData = dimBufferTurb[nameID];
                TestGPUMPIDataDecompressTurb(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                             offsetNgbLengthBuffer, faceIndexForRecv, HostMPIDataRecvTurb,
                                             HostRecvBufferTurb);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }
#endif
    }

    void DecompressInterfaceDataSmallGPUBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! printf("Program is running in GPUMPIDataDecompress\n");
        if (solverID == 0)
        { //! NS equation
            DecompressInterfaceNSDataGPUSmallBuffer(solverID, iinfo);
        }
        else if (solverID == 1)
        { //! Turb equation
            DecompressInterfaceTurbDataGPUSmallBuffer(solverID, iinfo);
        }
    }

    void TransferInterfaceDataNSMPI(const int solverID, InterfaceInfo *iinfo)
    {
        //! Add for test
        //! printf("Progam is running in TransferInterfaceDataNSMPI\n");
        //! Add end
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! for string total number of MPI operations
        int countRequest = 0;
#ifdef HDCOVERLAP
        cudaStreamSynchronize(dataStreamNS);
#endif

#ifdef MPIOVERLAP
        int countRequestFirst   = 0;
        int countRequestRemains = 0;
#endif

#ifndef MPIOVERLAP
        //! loop all of neighbor zones for MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthNS[ngbID];
            //! Get initial address of send data in GPUMPIDataSendNS
            double *sendData = GPUMPIDataSendNS + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            //! Add for test
            //! printf("iZone = %d, ngbID = %d, ngbProcess=%d\n", iZone, ngbID,
            //! ngbProcess); Add end
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestNS)
        {
            printf("countRequest is not equal to nTotalRequestNS, countRequest = %d, "
                   "nTotalRequestNS = %d\n",
                   countRequest, nTotalRequestNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerNS, MPI_STATUSES_IGNORE);
#else //! ifdef MPIOVERLAP                                                                          \
      //! loop all of neighbor zones for MPI operations on First NS buffer
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthNSFirst[ngbID];
            //! Get initial address of send data in GPUMPIDataSendNS
            double *sendData = GPUMPIDataSendNS + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerNSFirst[countRequestFirst] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerNSFirst[countRequestFirst]);
            countRequestFirst++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerNSFirst[countRequestFirst] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerNSFirst[countRequestFirst]);
            countRequestFirst++;
#endif
        }
        //! test total number of MPI operations
        if (countRequestFirst != nTotalRequestNS)
        {
            printf("countRequestFirst is not equal to nTotalRequestNS, countRequest = %d, "
                   "nTotalRequestNS = %d\n",
                   countRequestFirst, nTotalRequestNS);
            exit(1);
        }
        MPI_Waitall(countRequestFirst, requestContainerNSFirst, MPI_STATUSES_IGNORE);

        //! loop all of neighbor zones for MPI operations on Remain NS buffer
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get initial address of jZone in d_send
            //! int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthNSRemains[ngbID];
            //! Get initial address of send data in GPUMPIDataSendNS
            double *sendData = GPUMPIDataSendNS + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerNSRemains[countRequestRemains] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendNS + ngbBufferLengthNSFirst[ngbID] + offsetNgb, ngbLength * sizeof(double),
                      MPI_CHAR, ngbProcess, tagIZone, MPI_COMM_WORLD, &requestContainerNSRemains[countRequestRemains]);
            countRequestRemains++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerNSRemains[countRequestRemains] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvNS + ngbBufferLengthNSFirst[ngbID] + offsetNgb, ngbLength * sizeof(double),
                      MPI_CHAR, ngbProcess, tagJZone, MPI_COMM_WORLD, &requestContainerNSRemains[countRequestRemains]);
            countRequestRemains++;
#endif
        }
        //! test total number of MPI operations
        if (countRequestRemains != nTotalRequestNS)
        {
            printf("countRequestRemains is not equal to nTotalRequestNS, countRequest = "
                   "%d, nTotalRequestNS = %d\n",
                   countRequestRemains, nTotalRequestNS);
            exit(1);
        }
        //! It should be placed on Decompress of Turb interface data.
        //! MPI_Waitall(countRequestRemains, requestContainerNSRemains,
        //! MPI_STATUSES_IGNORE);
#endif //! end ifndef MPIOVERLAP

        /*
  //!It is used for finishing turb remain MPI data decompress and cudaMemcpy. It
  should be not here. It should be placed in a position between DecompressNS and
  Turb compute. #ifdef HDCOVERLAP cudaStreamSynchronize(dataStreamTurb); #endif
  */
    }

    void TransferInterfaceDataTurbMPI(const int solverID, InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! for string total number of MPI operations
        int countRequest = 0;

#ifdef HDCOVERLAP
        cudaStreamSynchronize(dataStreamTurb);
#endif

#ifdef MPIOVERLAP
        int countRequestFirst   = 0;
        int countRequestRemains = 0;
#endif

#ifndef MPIOVERLAP
        //! loop all of neighbor zones only for MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthTurb[ngbID];
            //! Get initial address of send data in GPUMPIDataSendTurb
            double *sendData = GPUMPIDataSendTurb + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestTurb)
        {
            printf("countRequest is not equal to nTotalRequestTurb, countRequest = %d, "
                   "nTotalRequestTurb = %d\n",
                   countRequest, nTotalRequestTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerTurb, MPI_STATUSES_IGNORE);
#else //! ifdef MPIOVERLAP

        //! loop all of neighbor zones only for MPI First operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthTurbFirst[ngbID];
            //! Get initial address of send data in GPUMPIDataSendTurb
            double *sendData = GPUMPIDataSendTurb + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurbFirst[countRequestFirst] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerTurbFirst[countRequestFirst]);
            countRequestFirst++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurbFirst[countRequestFirst] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerTurbFirst[countRequestFirst]);
            countRequestFirst++;
#endif
        }
        //! test total number of MPI operations
        if (countRequestFirst != nTotalRequestTurb)
        {
            printf("countRequestFirst is not equal to nTotalRequestTurb, countRequest = "
                   "%d, nTotalRequestTurb = %d\n",
                   countRequestFirst, nTotalRequestTurb);
            exit(1);
        }
        MPI_Waitall(countRequestFirst, requestContainerTurbFirst, MPI_STATUSES_IGNORE);

        //! loop all of neighbor zones only for MPI Remain operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthTurbRemains[ngbID];
            //! Get initial address of send data in GPUMPIDataSendTurb
            double *sendData = GPUMPIDataSendTurb + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurbRemains[countRequestRemains] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendTurb + ngbBufferLengthTurbFirst[ngbID] + offsetNgb, ngbLength * sizeof(double),
                      MPI_CHAR, ngbProcess, tagIZone, MPI_COMM_WORLD,
                      &requestContainerTurbRemains[countRequestRemains]);
            countRequestRemains++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurbRemains[countRequestRemains] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvTurb + ngbBufferLengthTurbFirst[ngbID] + offsetNgb, ngbLength * sizeof(double),
                      MPI_CHAR, ngbProcess, tagJZone, MPI_COMM_WORLD,
                      &requestContainerTurbRemains[countRequestRemains]);
            countRequestRemains++;
#endif
        }
        //! test total number of MPI operations
        if (countRequestRemains != nTotalRequestTurb)
        {
            printf("countRequestRemains is not equal to nTotalRequestTurb, countRequest = "
                   "%d, nTotalRequestTurb = %d\n",
                   countRequestRemains, nTotalRequestTurb);
            exit(1);
        }
        //! It should be placed into Decompress NS Interface part
        //! MPI_Waitall(countRequestRemains, requestContainerTurbRemains,
        //! MPI_STATUSES_IGNORE);
#endif //! ifdef MPIOVERLAP

        /*It is used for finishing NS interface remain MPI data. It should not be
  here. It should be a position between DecompressTurb and start of NS equation.
  #ifdef HDCOVERLAP
          cudaStreamSynchronize(dataStreamNS);
  #endif
  */
    }

    void TransferInterfaceDataMPI(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        if (solverID == 0)
        { //! NS equation
            TransferInterfaceDataNSMPI(solverID, iinfo);
        }
        else if (solverID == 1)
        { //! Turb equation
            TransferInterfaceDataTurbMPI(solverID, iinfo);
        }
    }

    void CompressInterfaceDataNSGPUSmallBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        //! Step1: Compress interface data into small gpu buffer
        //! and Transfer data from small gpu buffer into large host buffer
        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! Get total interface size of iZone
        int nIFace = iinfo->GetNIFace();
        for (int ngbID = nTotalNgbZones - 1; ngbID >= 0; ngbID--)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompressSmallBuf, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompressSmallBuf, 0, gridSize, blockSize);
#endif

            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameNSBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferNS
                int offsetVar = offsetBufferNS[nameID];
                //! get initial address of variable in GPUSendBufferNS
                int dimData = dimBufferNS[nameID];
                //! set trafer data size
                size_t sizeBuf = dimData * nIFaceOfNeighbor * sizeof(double);
                //! set HostSendBuffer
                double *HostSendBuffer = HostMPIDataSendNS + offsetNgb + offsetNgbLengthBuffer;
                if (name == "q")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_q[ngbID]);
                    HANDLE_API_ERR(cudaMemcpy(HostSendBuffer, GPUDataSendNS_q[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "dqdx")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_dqdx[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendNS_dqdx[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "dqdy")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_dqdy[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendNS_dqdy[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "dqdz")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_dqdz[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendNS_dqdz[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "limit")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_limit[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendNS_limit[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "dtdx")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_dtdx[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendNS_dtdx[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "dtdy")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_dtdy[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendNS_dtdy[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "dtdz")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_dtdz[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendNS_dtdz[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "t")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferNS, GPUDataSendNS_t[ngbID]);
                    HANDLE_API_ERR(cudaMemcpy(HostSendBuffer, GPUDataSendNS_t[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

//! test HostMPIDataSendNS and HostSendLargeBufferNS
#ifdef TESTNEWORGMPI
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForSend      = iinfo->GetFaceIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                int dimData = dimBufferNS[nameID];
                //! Get start position for variable in HostSendBufferNSAWARE.
                TestCompressAndTransferInterpolateDataNS(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                         offsetNgbLengthBuffer, faceIndexForSend, HostMPIDataSendNS,
                                                         HostSendBufferNS);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

#endif
    }

    void CompressInterfaceDataTurbGPUSmallBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! Step1: Compress interface data into small gpu buffer
        //! and Transfer data from small gpu buffer into large host buffer
        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! Get total interface size of iZone
        int nIFace = iinfo->GetNIFace();
        for (int ngbID = nTotalNgbZones - 1; ngbID >= 0; ngbID--)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompressSmallBuf, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompressSmallBuf, 0, gridSize, blockSize);
#endif

            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameTurbBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferTurb
                int offsetVar = offsetBufferTurb[nameID];
                //! get initial address of variable in GPUSendBufferTurb
                int dimData = dimBufferTurb[nameID];
                //! set trafer data size
                size_t sizeBuf = dimData * nIFaceOfNeighbor * sizeof(double);
                //! set HostRecvBuffer
                double *HostSendBuffer = HostMPIDataSendTurb + offsetNgb + offsetNgbLengthBuffer;
                if (name == "turb::q")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferTurb, GPUDataSendTurb_q[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendTurb_q[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "turb::dqdx")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferTurb, GPUDataSendTurb_dqdx[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendTurb_dqdx[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "turb::dqdy")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferTurb, GPUDataSendTurb_dqdy[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendTurb_dqdy[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                else if (name == "turb::dqdz")
                {
                    GPUMPIDataCompressSmallBuf<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, nIFaceOfNeighbor,
                                                                        offsetFaces, d_faceIndexForSend,
                                                                        GPUSendBufferTurb, GPUDataSendTurb_dqdz[ngbID]);
                    HANDLE_API_ERR(
                        cudaMemcpy(HostSendBuffer, GPUDataSendTurb_dqdz[ngbID], sizeBuf, cudaMemcpyDeviceToHost));
                }
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }
//! Test HostMPIDataSendTurb and HostSendBufferTurb
#ifdef TESTNEWORGMPI
        size_t sizeAWARETest = sizeof(double) * lengthBufferTurb;
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForSend      = iinfo->GetFaceIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                int dimData = dimBufferTurb[nameID];
                //! Get start position for variable in HostSendBufferTurbAWARE.
                TestCompressAndTransferInterpolateDataTurb(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                           offsetNgbLengthBuffer, faceIndexForSend, HostMPIDataSendTurb,
                                                           HostSendBufferTurb);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

#endif
    }

    __global__ void GPUMPIDataCompressSmallBuf(const int offsetVar, const int neqn, const int nIFaceTotal,
                                               const int nIFaceOfNeighbor, const int offsetFaceIndex,
                                               const int *faceIndexForSend, const double *GPUSendBufferLarge,
                                               double *GPUMPIDataSendSmall)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int faceID;
        int equationID;
        int ngbZoneFaceID;
        int largeStorageID;
        int faceIndexLargeStore;
        for (faceID = threadID; faceID < nIFaceOfNeighbor; faceID += blockDim.x * gridDim.x)
        {
            faceIndexLargeStore = faceIndexForSend[offsetFaceIndex + faceID];

            for (equationID = 0; equationID < neqn; equationID++)
            {
                ngbZoneFaceID                      = equationID * nIFaceOfNeighbor + faceID;
                largeStorageID                     = offsetVar + equationID * nIFaceTotal + faceIndexLargeStore;
                GPUMPIDataSendSmall[ngbZoneFaceID] = GPUSendBufferLarge[largeStorageID];
            }
        }
    }

    void CompressInterfaceDataSmallGPUBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        if (solverID == 0)
        { //! NS equation
            CompressInterfaceDataNSGPUSmallBuffer(solverID, iinfo);
        }
        else if (solverID == 1)
        { //! Turb equation
            CompressInterfaceDataTurbGPUSmallBuffer(solverID, iinfo);
        }
    }

    void SetGlobalMinDt(double &globalMinDt)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
#ifdef TESTNEWORGMPI
#ifdef CUDA_AWARE_MPI
        HANDLE_API_ERR(cudaMemcpy(&h_localMinDt, d_globalMinDt, 1 * sizeof(double), cudaMemcpyDeviceToHost));
#endif
        TestAllreduceMinDt(h_localMinDt, globalMinDt);
#endif

        //! Set globalMinDt by localMinDt
        globalMinDt = h_localMinDt;
    }

    void ComputeMinimumDtAllreduce(double &globalMinDt)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
//! Get minimum dt in localMinDt, noting that first parameter is send buffer and
//! the second parameter is recv buffer
#ifdef CUDA_AWARE_MPI
        //! cudaDeviceSynchronize();
        cudaStreamSynchronize(0);
        MPI_Allreduce(d_minDt, d_globalMinDt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
        MPI_Allreduce(&globalMinDt, &h_localMinDt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

#endif
    }

    void DecompressInterpointData(const int solverID, InterpointInformation *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        if (solverID == 0)
        { //! NS equation
#ifndef CUDA_AWARE_MPI
            size_t sizeAWARE = sizeof(double) * lengthInterpointBufferNS;
            HANDLE_API_ERR(
                cudaMemcpy(GPUMPIDataRecvInterpointNS, HostMPIDataRecvInterpointNS, sizeAWARE, cudaMemcpyHostToDevice));
#endif
            int iZone          = sendGlobalZoneForPoint[0];
            int nTotalNgbZones = numNgbZonesForPoint[0];
            int nIPoint        = iinfo->GetNumberOfInterpoints();

            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int ngbZoneID             = recvGlobalZoneForPoint[ngbID];
                int iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int nIPointOfNeighbor     = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int offsetNgbLengthBuffer = 0;
                int offsetPoints          = offsetPointIndex[ngbID];
                int offsetNgb             = offsetMPIDataInterpointNS[ngbID];

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIPointOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataDecompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataDecompress, 0, gridSize, blockSize);
#endif

                for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSInterpointBuffer)[nameID];
                    int           offsetVar = offsetInterpointBufferNS[nameID];
                    int           dimData   = dimInterpointBufferNS[nameID];
                    GPUMPIDataDecompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor, offsetNgbLengthBuffer, offsetPoints,
                        d_pointIndexForRecv, GPUMPIDataRecvInterpointNS, GPURecvInterpointBufferNS);

#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
                }

            } //! end for ngbID

#ifdef HDCOVERLAP
            cudaStreamSynchronize(dataStreamTurb);
#endif

#ifdef TESTNEWORGMPI
            //! Test GPUMPIDataRecvNS and HostRecvLargeBuffer
            double *HostRecvInterpointBufferNSAWARE = new double[lengthInterpointBufferNS];
            size_t  sizeAWARETest                   = sizeof(double) * lengthInterpointBufferNS;
            HANDLE_API_ERR(cudaMemcpy(HostRecvInterpointBufferNSAWARE, GPUMPIDataRecvInterpointNS, sizeAWARETest,
                                      cudaMemcpyDeviceToHost));
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int  ngbZoneID             = recvGlobalZoneForPoint[ngbID];
                int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int  nIPointOfNeighbor     = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForRecv     = iinfo->GetPointIndexForRecv(iNeighbor);
                int  offsetNgbLengthBuffer = 0;
                for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
                {
                    int dimData = dimInterpointBufferNS[nameID];
                    TestGPUMPIDataDecompressInterpointNS(iZone, ngbZoneID, nameID, ngbID, nIPoint, nIPointOfNeighbor,
                                                         offsetNgbLengthBuffer, pointIndexForRecv,
                                                         HostRecvInterpointBufferNSAWARE, HostRecvInterpointBufferNS);

                    offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
                }
            }

            delete[] HostRecvInterpointBufferNSAWARE;
#endif
        }
        else if (solverID == 1)
        { //! Turbulence equation
#ifndef CUDA_AWARE_MPI
            size_t sizeAWARE = sizeof(double) * lengthInterpointBufferTurb;
            HANDLE_API_ERR(cudaMemcpy(GPUMPIDataRecvInterpointTurb, HostMPIDataRecvInterpointTurb, sizeAWARE,
                                      cudaMemcpyHostToDevice));
#endif
            int iZone          = sendGlobalZoneForPoint[0];
            int nTotalNgbZones = numNgbZonesForPoint[0];
            int nIPoint        = iinfo->GetNumberOfInterpoints();
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int ngbZoneID             = recvGlobalZoneForPoint[ngbID];
                int iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int nIPointOfNeighbor     = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int offsetNgbLengthBuffer = 0;
                int offsetPoints          = offsetPointIndex[ngbID];
                int offsetNgb             = offsetMPIDataInterpointTurb[ngbID];

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIPointOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataDecompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataDecompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
                {
                    const string &name      = (*nameTurbInterpointBuffer)[nameID];
                    int           offsetVar = offsetInterpointBufferTurb[nameID];
                    int           dimData   = dimInterpointBufferTurb[nameID];
                    GPUMPIDataDecompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor, offsetNgbLengthBuffer, offsetPoints,
                        d_pointIndexForRecv, GPUMPIDataRecvInterpointTurb, GPURecvInterpointBufferTurb);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif

                    offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
                }

            } //! end for ngbID

#ifdef HDCOVERLAP
            cudaStreamSynchronize(dataStreamNS);
#endif

#ifdef TESTNEWORGMPI
            //! Test GPUMPIDataRecvTurb and HostRecvLargeBuffer
            double *HostRecvInterpointBufferTurbAWARE = new double[lengthInterpointBufferTurb];
            size_t  sizeAWARETest                     = sizeof(double) * lengthInterpointBufferTurb;
            HANDLE_API_ERR(cudaMemcpy(HostRecvInterpointBufferTurbAWARE, GPUMPIDataRecvInterpointTurb, sizeAWARETest,
                                      cudaMemcpyDeviceToHost));
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int  ngbZoneID             = recvGlobalZoneForPoint[ngbID];
                int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int  nIPointOfNeighbor     = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForRecv     = iinfo->GetPointIndexForRecv(iNeighbor);
                int  offsetNgbLengthBuffer = 0;
                for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
                {
                    int dimData = dimInterpointBufferTurb[nameID];
                    TestGPUMPIDataDecompressInterpointTurb(
                        iZone, ngbZoneID, nameID, ngbID, nIPoint, nIPointOfNeighbor, offsetNgbLengthBuffer,
                        pointIndexForRecv, HostRecvInterpointBufferTurbAWARE, HostRecvInterpointBufferTurb);

                    offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
                }
            }

            delete[] HostRecvInterpointBufferTurbAWARE;
#endif
        }
    }

    void DecompressInterpolateData(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! printf("Program is running in GPUMPIDataDecompress\n");
        if (solverID == 0)
        { //! NS equation
#ifndef CUDA_AWARE_MPI
            size_t sizeAWARE = sizeof(double) * lengthBufferNS;
            HANDLE_API_ERR(cudaMemcpy(GPUMPIDataRecvNS, HostMPIDataRecvNS, sizeAWARE, cudaMemcpyHostToDevice));
#endif
            int iZone          = sendGlobalZone[0];
            int nTotalNgbZones = numNgbZones[0];
            int nIFace         = iinfo->GetNIFace();
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int ngbZoneID             = recvGlobalZone[ngbID];
                int iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int offsetNgbLengthBuffer = 0;
                int offsetFaces           = offsetFaceIndex[ngbID];
                int offsetNgb             = offsetMPIDataNS[ngbID];

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataDecompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataDecompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSBuffer)[nameID];
                    int           offsetVar = offsetBufferNS[nameID];
                    int           dimData   = dimBufferNS[nameID];
                    GPUMPIDataDecompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor, offsetNgbLengthBuffer, offsetFaces,
                        d_faceIndexForRecv, GPUMPIDataRecvNS, GPURecvBufferNS);

#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }

            } //! end for ngbID
#ifdef TESTNEWORGMPI
            //! Test GPUMPIDataRecvNS and HostRecvLargeBuffer
            double *HostRecvBufferNSAWARE = new double[lengthBufferNS];
            size_t  sizeAWARETest         = sizeof(double) * lengthBufferNS;
            HANDLE_API_ERR(cudaMemcpy(HostRecvBufferNSAWARE, GPUMPIDataRecvNS, sizeAWARETest, cudaMemcpyDeviceToHost));
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int  ngbZoneID             = recvGlobalZone[ngbID];
                int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int *faceIndexForRecv      = iinfo->GetFaceIndexForRecv(iNeighbor);
                int  offsetNgbLengthBuffer = 0;
                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    int dimData = dimBufferNS[nameID];
                    TestGPUMPIDataDecompressNS(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                               offsetNgbLengthBuffer, faceIndexForRecv, HostRecvBufferNSAWARE,
                                               HostRecvBufferNS);

                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
            }

            delete[] HostRecvBufferNSAWARE;
#endif
        }
        else if (solverID == 1)
        { //! Turbulence equation
#ifndef CUDA_AWARE_MPI
            size_t sizeAWARE = sizeof(double) * lengthBufferTurb;
            HANDLE_API_ERR(cudaMemcpy(GPUMPIDataRecvTurb, HostMPIDataRecvTurb, sizeAWARE, cudaMemcpyHostToDevice));
#endif
            int iZone          = sendGlobalZone[0];
            int nTotalNgbZones = numNgbZones[0];
            int nIFace         = iinfo->GetNIFace();
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int ngbZoneID             = recvGlobalZone[ngbID];
                int iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int offsetNgbLengthBuffer = 0;
                int offsetFaces           = offsetFaceIndex[ngbID];
                int offsetNgb             = offsetMPIDataTurb[ngbID];

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataDecompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataDecompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                {
                    const string &name      = (*nameTurbBuffer)[nameID];
                    int           offsetVar = offsetBufferTurb[nameID];
                    int           dimData   = dimBufferTurb[nameID];
                    GPUMPIDataDecompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor, offsetNgbLengthBuffer, offsetFaces,
                        d_faceIndexForRecv, GPUMPIDataRecvTurb, GPURecvBufferTurb);

#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }

            } //! end for ngbID
#ifdef TESTNEWORGMPI
            //! Test GPUMPIDataRecvTurb and HostRecvLargeBuffer
            double *HostRecvBufferTurbAWARE = new double[lengthBufferTurb];
            size_t  sizeAWARETest           = sizeof(double) * lengthBufferTurb;
            HANDLE_API_ERR(
                cudaMemcpy(HostRecvBufferTurbAWARE, GPUMPIDataRecvTurb, sizeAWARETest, cudaMemcpyDeviceToHost));
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int  ngbZoneID             = recvGlobalZone[ngbID];
                int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int *faceIndexForRecv      = iinfo->GetFaceIndexForRecv(iNeighbor);
                int  offsetNgbLengthBuffer = 0;
                for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                {
                    int dimData = dimBufferTurb[nameID];
                    TestGPUMPIDataDecompressTurb(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                 offsetNgbLengthBuffer, faceIndexForRecv, HostRecvBufferTurbAWARE,
                                                 HostRecvBufferTurb);

                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
            }

            delete[] HostRecvBufferTurbAWARE;
#endif
        }
    }

    void DecompressInterpolateDataMPIOverlap(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! printf("Program is running in GPUMPIDataDecompress\n");
        if (solverID == 0)
        { //! NS equation
#ifdef MPIOVERLAP
            if (GPUMPIDataRecvTurb)
            {
                MPI_Waitall(nTotalRequestTurb, requestContainerTurbRemains, MPI_STATUSES_IGNORE);
#ifndef CUDA_AWARE_MPI
                size_t sizeRemains = sizeof(double) * lengthBufferTurb;
                HANDLE_API_ERR(
                    cudaMemcpy(GPUMPIDataRecvTurb, HostMPIDataRecvTurb, sizeRemains, cudaMemcpyHostToDevice));
#endif
                int nTotalNgbZonesTurb = numNgbZones[0];
                int nIFaceTurb         = iinfo->GetNIFace();
                for (int ngbID = 0; ngbID < nTotalNgbZonesTurb; ngbID++)
                {
                    int ngbZoneID             = recvGlobalZone[ngbID];
                    int iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                    int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                    int offsetNgbLengthBuffer = 0;
                    int offsetFaces           = offsetFaceIndex[ngbID];
                    int offsetNgb             = offsetMPIDataTurb[ngbID];

                    int gridSize       = 1;
                    int blockSize      = 1;
                    int loopLen        = nIFaceOfNeighbor;
                    int blockSizeLimit = 0;
                    KernelLaunchPara((void *)GPUMPIDataDecompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                    ReportKernelPara((void *)GPUMPIDataDecompress, 0, gridSize, blockSize);
#endif
                    for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                    {
                        const string &name      = (*nameTurbBuffer)[nameID];
                        int           offsetVar = offsetBufferTurb[nameID];
                        int           dimData   = dimBufferTurb[nameID];
                        GPUMPIDataDecompress<<<gridSize, blockSize>>>(
                            offsetVar, dimData, nIFaceTurb, offsetNgb, nIFaceOfNeighbor, offsetNgbLengthBuffer,
                            offsetFaces, d_faceIndexForRecv, GPUMPIDataRecvTurb, GPURecvBufferTurb);

#ifdef KERNELLAUNCHTEST
                        HANDLE_KERNEL_ERR();
#endif
                        offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                    }

                } //! end for ngbID
            }
#endif //! end ifdef MPIOVERLAP
#ifndef CUDA_AWARE_MPI
            size_t sizeAWARE = sizeof(double) * lengthBufferNS;
            HANDLE_API_ERR(cudaMemcpy(GPUMPIDataRecvNS, HostMPIDataRecvNS, sizeAWARE, cudaMemcpyHostToDevice));
#endif
            int iZone          = sendGlobalZone[0];
            int nTotalNgbZones = numNgbZones[0];
            int nIFace         = iinfo->GetNIFace();
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int ngbZoneID             = recvGlobalZone[ngbID];
                int iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int offsetNgbLengthBuffer = 0;
                int offsetFaces           = offsetFaceIndex[ngbID];
                int offsetNgb             = offsetMPIDataNS[ngbID];

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataDecompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataDecompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSBuffer)[nameID];
                    int           offsetVar = offsetBufferNS[nameID];
                    int           dimData   = dimBufferNS[nameID];
                    GPUMPIDataDecompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor, offsetNgbLengthBuffer, offsetFaces,
                        d_faceIndexForRecv, GPUMPIDataRecvNS, GPURecvBufferNS);

#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }

            } //! end for ngbID
#ifdef TESTNEWORGMPI
            //! Test GPUMPIDataRecvNS and HostRecvLargeBuffer
            double *HostRecvBufferNSAWARE = new double[lengthBufferNS];
            size_t  sizeAWARETest         = sizeof(double) * lengthBufferNS;
            HANDLE_API_ERR(cudaMemcpy(HostRecvBufferNSAWARE, GPUMPIDataRecvNS, sizeAWARETest, cudaMemcpyDeviceToHost));
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int  ngbZoneID             = recvGlobalZone[ngbID];
                int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int *faceIndexForRecv      = iinfo->GetFaceIndexForRecv(iNeighbor);
                int  offsetNgbLengthBuffer = 0;
                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    int dimData = dimBufferNS[nameID];
                    TestGPUMPIDataDecompressNS(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                               offsetNgbLengthBuffer, faceIndexForRecv, HostRecvBufferNSAWARE,
                                               HostRecvBufferNS);

                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
            }

            delete[] HostRecvBufferNSAWARE;
#endif
        }
        else if (solverID == 1)
        { //! Turbulence equation
#ifdef MPIOVERLAP
            MPI_Waitall(nTotalRequestNS, requestContainerNSRemains, MPI_STATUSES_IGNORE);
#ifndef CUDA_AWARE_MPI
            size_t sizeRemains = sizeof(double) * lengthBufferNS;
            HANDLE_API_ERR(cudaMemcpy(GPUMPIDataRecvNS, HostMPIDataRecvNS, sizeRemains, cudaMemcpyHostToDevice));
#endif
            int nTotalNgbZonesNS = numNgbZones[0];
            int nIFaceNS         = iinfo->GetNIFace();
            for (int ngbID = 0; ngbID < nTotalNgbZonesNS; ngbID++)
            {
                int ngbZoneID             = recvGlobalZone[ngbID];
                int iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int offsetNgbLengthBuffer = 0;
                int offsetFaces           = offsetFaceIndex[ngbID];
                int offsetNgb             = offsetMPIDataNS[ngbID];

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataDecompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataDecompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSBuffer)[nameID];
                    int           offsetVar = offsetBufferNS[nameID];
                    int           dimData   = dimBufferNS[nameID];
                    GPUMPIDataDecompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFaceNS, offsetNgb, nIFaceOfNeighbor, offsetNgbLengthBuffer, offsetFaces,
                        d_faceIndexForRecv, GPUMPIDataRecvNS, GPURecvBufferNS);

#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }

            } //! end for ngbID
#endif //! end ifdef MPIOVERLAP
#ifndef CUDA_AWARE_MPI
            size_t sizeAWARE = sizeof(double) * lengthBufferTurb;
            HANDLE_API_ERR(cudaMemcpy(GPUMPIDataRecvTurb, HostMPIDataRecvTurb, sizeAWARE, cudaMemcpyHostToDevice));
#endif
            int iZone          = sendGlobalZone[0];
            int nTotalNgbZones = numNgbZones[0];
            int nIFace         = iinfo->GetNIFace();
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int ngbZoneID             = recvGlobalZone[ngbID];
                int iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int offsetNgbLengthBuffer = 0;
                int offsetFaces           = offsetFaceIndex[ngbID];
                int offsetNgb             = offsetMPIDataTurb[ngbID];

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataDecompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataDecompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                {
                    const string &name      = (*nameTurbBuffer)[nameID];
                    int           offsetVar = offsetBufferTurb[nameID];
                    int           dimData   = dimBufferTurb[nameID];
                    GPUMPIDataDecompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor, offsetNgbLengthBuffer, offsetFaces,
                        d_faceIndexForRecv, GPUMPIDataRecvTurb, GPURecvBufferTurb);

#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }

            } //! end for ngbID
#ifdef TESTNEWORGMPI
            //! Test GPUMPIDataRecvTurb and HostRecvLargeBuffer
            double *HostRecvBufferTurbAWARE = new double[lengthBufferTurb];
            size_t  sizeAWARETest           = sizeof(double) * lengthBufferTurb;
            HANDLE_API_ERR(
                cudaMemcpy(HostRecvBufferTurbAWARE, GPUMPIDataRecvTurb, sizeAWARETest, cudaMemcpyDeviceToHost));
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int  ngbZoneID             = recvGlobalZone[ngbID];
                int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int *faceIndexForRecv      = iinfo->GetFaceIndexForRecv(iNeighbor);
                int  offsetNgbLengthBuffer = 0;
                for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                {
                    int dimData = dimBufferTurb[nameID];
                    TestGPUMPIDataDecompressTurb(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                 offsetNgbLengthBuffer, faceIndexForRecv, HostRecvBufferTurbAWARE,
                                                 HostRecvBufferTurb);

                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
            }

            delete[] HostRecvBufferTurbAWARE;
#endif
        }
    }
    void CompressAndTransferInterpolateData(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! printf("Program is running in CompressAndTransferInterpolateData\n");
        //! printf("solverID = %d\n", solverID);
        /*
  //!test zone connectivity
  int iZone = sendGlobalZone[0];
  int nTotalNgbZones = numNgbZones[0];
  printf("Zone %d owns %d neighbor zones\n", iZone, nTotalNgbZones);
  for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++){
          printf("Zone %d owns %d neighbor zone in process %d\n", iZone,
  recvGlobalZone[ngbID], recvProcess[ngbID]);
  }
  */
        if (solverID == 0)
        { //! NS equation
            //! rearange data in GPUSendBufferNS
            int iZone          = sendGlobalZone[0];
            int nTotalNgbZones = numNgbZones[0];
            int nIFace         = iinfo->GetNIFace();
            //! vector <MPI_Request> requestContainer;
            int countRequest = 0;
#ifndef ORGISENDIRECVORDER
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int ngbZoneID = recvGlobalZone[ngbID];
                int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
                //! it may be a bug?? here or creation part
                int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                //! faceIndexForSend should be put in the device
                //! int *faceIndexForSend = iterfaceInfo->GetFaceIndexForSend(ineighbor)
                int offsetNgbLengthBuffer = 0;
                int offsetFaces           = offsetFaceIndex[ngbID];
                int offsetNgb             = offsetMPIDataNS[ngbID];

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSBuffer)[nameID];
                    int           offsetVar = offsetBufferNS[nameID];
                    int           dimData   = dimBufferNS[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                                offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                                GPUSendBufferNS, GPUMPIDataSendNS);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
                //! send data
                double *sendData = GPUMPIDataSendNS + offsetNgb;
                //! double ** sendData = &GPUMPIDataSendNS + offsetNgb;
                int ngbProcess = recvProcess[ngbID];
                //! size_t ngbSize = offsetNgbLengthBuffer * sizeof(double);
                int ngbLength = offsetNgbLengthBuffer;
                //! int offsetTagIJ = 1000;
                int offsetTagIJ = offsetTagInterfaceNS;
                int tagIZone    = iZone + offsetTagIJ;
//! MPI_Isend(sendData, ngbLength, MPI_DOUBLE, ngbProcess, tagIZone,
//! MPI_COMM_WORLD, &MPI_REQUEST_NULL); test printf("iZone = %d, ngbProcess = %d,
//! ngbLength = %d, tagIZone = %d\n", iZone, ngbProcess, ngbLength, tagIZone);
#ifndef CUDA_AWARE_MPI
                size_t sizeAWARE = sizeof(double) * lengthBufferNS;
                HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendNS, GPUMPIDataSendNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif
#ifdef ISENDIRECVMPI
                //! requestContainer.push_back(MPI_REQUEST_NULL);
                //! MPI_Isend(sendData, ngbLength, MPI_DOUBLE, ngbProcess, tagIZone,
                //! MPI_COMM_WORLD, &requestContainer.back());
#ifdef CUDA_AWARE_MPI
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                //! MPI_Isend(sendData, ngbLength, MPI_DOUBLE, ngbProcess, tagIZone,
                //! MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                MPI_Isend(GPUMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#else
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                //! MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength, MPI_DOUBLE,
                //! ngbProcess, tagIZone, MPI_COMM_WORLD,
                //! &requestContainerNS[countRequest]);
                MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#endif
#endif
                //! recv data
                double *recvData = GPUMPIDataRecvNS + offsetNgb;
                //! double ** recvData = &GPUMPIDataRecvNS + offsetNgb;
                int tagJZone = ngbZoneID + offsetTagIJ;
//! MPI_Irecv(recvData, ngbLength, MPI_DOUBLE, ngbProcess, tagJZone,
//! MPI_COMM_WORLD, &MPI_REQUEST_NULL); test printf("iZone = %d, ngbProcess = %d,
//! ngbLength = %d, tagJZone = %d\n", iZone, ngbProcess, ngbLength, tagJZone);
#ifdef ISENDIRECVMPI
#ifdef CUDA_AWARE_MPI
                //! requestContainer.push_back(MPI_REQUEST_NULL);
                //! MPI_Irecv(recvData, ngbLength, MPI_DOUBLE, ngbProcess, tagJZone,
                //! MPI_COMM_WORLD, &requestContainer.back());
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                //! MPI_Irecv(recvData, ngbLength, MPI_DOUBLE, ngbProcess, tagJZone,
                //! MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                MPI_Irecv(GPUMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#else
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                //! MPI_Irecv(recvData, ngbLength, MPI_DOUBLE, ngbProcess, tagJZone,
                //! MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#endif
#endif
#ifndef ISENDIRECVMPI
                MPI_Status status;
//! MPI_Sendrecv(sendData, ngbLength, MPI_DOUBLE, ngbProcess, tagIZone, recvData,
//! ngbLength, MPI_DOUBLE, ngbProcess, tagJZone, MPI_COMM_WORLD, &status);
#ifdef CUDA_AWARE_MPI
                MPI_Sendrecv(GPUMPIDataSendNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess, tagIZone,
                             GPUMPIDataRecvNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess, tagJZone, MPI_COMM_WORLD,
                             &status);
#else
                MPI_Sendrecv(HostMPIDataSendNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess, tagIZone,
                             HostMPIDataRecvNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess, tagJZone, MPI_COMM_WORLD,
                             &status);

#endif
                //! check MPI_Sendrecv results
                int recvElems;
                MPI_Get_count(&status, MPI_DOUBLE, &recvElems);
                if (recvElems != ngbLength)
                {
                    printf("MPI_Sendrecv Error in iZone=%d, ngbZoneID = %d, ngbProcess = %d, "
                           "recvElems = %d, ngbLength=%d\n",
                           iZone, ngbZoneID, ngbProcess, recvElems, ngbLength);
                    exit(1);
                }
/*
//!dead lock
MPI_Send(GPUMPIDataSendNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess,
tagIZone, MPI_COMM_WORLD); MPI_Recv(GPUMPIDataRecvNS + offsetNgb, ngbLength,
MPI_DOUBLE, ngbProcess, tagJZone, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
*/
#endif
            } //! end for ngbID
#ifdef ISENDIRECVMPI
            //! int countRequest = static_cast<int>(requestContainer.size());
            //! MPI_Waitall(countRequest, &(requestContainer[0]), MPI_STATUSES_IGNORE);
            //! noting that when MPI_Waitall is finished, requestContainer can be
            //! cleared. Otherwise, MPI_Waitall cannot wait for request results, leading
            //! to incomplete transfer. Furthermore, std::vector<Object> will not clear
            //! memory. So, it should be used as a global variable.
            //! requestContainer.clear();
            if (countRequest != nTotalRequestNS)
            {
                printf("countRequest is not equal to nTotalRequestNS, countRequest = %d, "
                       "nTotalRequestNS = %d\n",
                       countRequest, nTotalRequestNS);
                exit(1);
            }
            MPI_Waitall(countRequest, requestContainerNS, MPI_STATUSES_IGNORE);
//! just for test
//! MPI_Barrier(MPI_COMM_WORLD);
#endif //! end ifdef ISENDIRECVMPI

#else //! else ORGISENDIRECVORDER
            int nTotalMPIOperations = nTotalRequestNS;
            for (int mpiID = 0; mpiID < nTotalRequestNS; mpiID++)
            {
                //! get neighbor local label zone ID
                int ngbID            = sendRecvNgbID[mpiID]; //! from 0 to nTotalNgbZones-1
                int sendOrRecv       = isSendOrRecv[mpiID];
                int ngbZoneID        = recvGlobalZone[ngbID];
                int iNeighbor        = iinfo->FindIthNeighbor(ngbZoneID);
                int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int offsetFaces      = offsetFaceIndex[ngbID];
                int offsetNgb        = offsetMPIDataNS[ngbID];
                int ngbProcess       = recvProcess[ngbID];
                int ngbLength        = ngbBufferLengthNS[ngbID];
                //! int offsetTagIJ = 1000;
                int offsetTagIJ = offsetTagInterfaceNS;

                //! MPI_Isend operation
                if (sendOrRecv == 1)
                {
                    int offsetNgbLengthBuffer = 0;

                    int gridSize       = 1;
                    int blockSize      = 1;
                    int loopLen        = nIFaceOfNeighbor;
                    int blockSizeLimit = 0;
                    KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                    ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                    for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                    {
                        const string &name      = (*nameNSBuffer)[nameID];
                        int           offsetVar = offsetBufferNS[nameID];
                        int           dimData   = dimBufferNS[nameID];
                        GPUMPIDataCompress<<<gridSize, blockSize>>>(
                            offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor, offsetNgbLengthBuffer, offsetFaces,
                            d_faceIndexForSend, GPUSendBufferNS, GPUMPIDataSendNS);
#ifdef KERNELLAUNCHTEST
                        HANDLE_KERNEL_ERR();
#endif
                        offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                    }

                    double *sendData = GPUMPIDataSendNS + offsetNgb;
                    int     tagIZone = iZone + offsetTagIJ;
#ifndef CUDA_AWARE_MPI
                    size_t sizeAWARE = sizeof(double) * lengthBufferNS;
                    HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendNS, GPUMPIDataSendNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
                    requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                    MPI_Isend(GPUMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                              MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                    countRequest++;
#else
                    requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                    MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                              MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                    countRequest++;
#endif
                }
                else if (sendOrRecv == 0)
                { //! receive data
                    double *recvData = GPUMPIDataRecvNS + offsetNgb;
                    int     tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                    requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                    MPI_Irecv(GPUMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                              MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                    countRequest++;
#else
                    requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                    MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                              MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                    countRequest++;
#endif
                }
            } //! end for mpiID

#ifdef ISENDIRECVMPI
            if (countRequest != nTotalRequestNS)
            {
                printf("countRequest is not equal to nTotalRequestNS, countRequest = %d, "
                       "nTotalRequestNS = %d\n",
                       countRequest, nTotalRequestNS);
                exit(1);
            }
            MPI_Waitall(countRequest, requestContainerNS, MPI_STATUSES_IGNORE);
#endif //! end ifdef ISENDIRECVMPI
#endif //! end ifdef ORGISENDIRECVORDER                                                                                \
       //! just for test                                                                                               \
       //! MPI_Barrier(MPI_COMM_WORLD);

#ifdef TESTNEWORGMPI
            //! test GPUMPIDataSendNS and HostSendLargeBufferNS
            double *HostSendBufferNSAWARE = new double[lengthBufferNS];
            size_t  sizeAWARETest         = sizeof(double) * lengthBufferNS;
            HANDLE_API_ERR(cudaMemcpy(HostSendBufferNSAWARE, GPUMPIDataSendNS, sizeAWARETest, cudaMemcpyDeviceToHost));
            for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
            {
                int  ngbZoneID             = recvGlobalZone[ngbID];
                int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
                int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                int *faceIndexForSend      = iinfo->GetFaceIndexForSend(iNeighbor);
                int  offsetNgbLengthBuffer = 0;

                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    int dimData = dimBufferNS[nameID];
                    //! Get start position for variable in HostSendBufferNSAWARE.
                    TestCompressAndTransferInterpolateDataNS(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                             offsetNgbLengthBuffer, faceIndexForSend,
                                                             HostSendBufferNSAWARE, HostSendBufferNS);

                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
            }

            delete[] HostSendBufferNSAWARE;
#endif
        }
        else if (solverID == 1)
        { //! Turb equation
        }
    }

    void CompressAndTransferInterpointData(const int solverID, InterpointInformation *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        if (solverID == 0)
        {            //! NS equation
#ifdef ISENDIRECVMPI //! MPI_Isend and MPI_Irecv
#ifdef CATOVERLAP
            CATOverlapInterpointIsendIrecvNS(iinfo);
#else
            CATSeparateInterpointIsendIrecvNS(iinfo);
#endif
#else //! MPI_Sendrecv
            CATInterpointSendrecvNS(iinfo);
#endif
        }
        else if (solverID == 1)
        {            //! Turb equation
#ifdef ISENDIRECVMPI //! MPI_Isend and MPI_Irecv
#ifdef CATOVERLAP
            CATOverlapInterpointIsendIrecvTurb(iinfo);
#else
            CATSeparateInterpointIsendIrecvTurb(iinfo);
#endif
#else //! MPI_Sendrecv
            CATInterpointSendrecvTurb(iinfo);
#endif
        }
    }

    void CATInterpointSendrecvNS(InterpointInformation *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
    }

    void CATInterpointSendrecvTurb(InterpointInformation *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
    }

    void CompressAndTransferInterpolateDataV2(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        if (solverID == 0)
        {            //! NS equation
#ifdef ISENDIRECVMPI //! MPI_Isend and MPI_Irecv
#ifdef CATOVERLAP
            CATOverlapInterfaceIsendIrecvNS(iinfo);
#else
            CATSeparateInterfaceIsendIrecvNS(iinfo);
#endif
#else //! MPI_Sendrecv
            CATInterfaceSendrecvNS(iinfo);
#endif
        }
        else if (solverID == 1)
        {            //! Turb equation
#ifdef ISENDIRECVMPI //! MPI_Isend and MPI_Irecv
#ifdef CATOVERLAP
            CATOverlapInterfaceIsendIrecvTurb(iinfo);
#else
            CATSeparateInterfaceIsendIrecvTurb(iinfo);
#endif
#else //! MPI_Sendrecv
            CATInterfaceSendrecvTurb(iinfo);
#endif
        }
    }

    void CATOverlapInterfaceIsendIrecvTurb(InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! Get total interface size of iZone
        int nIFace = iinfo->GetNIFace();
        //! for string total number of MPI operations
        int countRequest = 0;

#ifndef ORGISENDIRECVORDER //! NOT original MPI operations' order
        //! loop all of neighbor zones for MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthTurb[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameTurbBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferTurb
                int offsetVar = offsetBufferTurb[nameID];
                //! get dim of variable
                int dimData = dimBufferTurb[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                            GPUSendBufferTurb, GPUMPIDataSendTurb);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
            //! Get initial address of send data in GPUMPIDataSendTurb
            double *sendData = GPUMPIDataSendTurb + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifndef CUDA_AWARE_MPI //! just by traditional MPI
            size_t sizeAWARE = sizeof(double) * lengthBufferTurb;
            HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendTurb, GPUMPIDataSendTurb, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(GPUMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#else
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(GPUMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#else
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestTurb)
        {
            printf("countRequest is not equal to nTotalRequestTurb, countRequest = %d, "
                   "nTotalRequestTurb = %d\n",
                   countRequest, nTotalRequestTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerTurb, MPI_STATUSES_IGNORE);
#else
        int nTotalMPIOperations = nTotalRequestTurb;
        for (int mpiID = 0; mpiID < nTotalRequestTurb; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbID[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecv[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcess[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbBufferLengthTurb[ngbID];
            //! int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendTurb
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                {
                    const string &name      = (*nameTurbBuffer)[nameID];
                    int           offsetVar = offsetBufferTurb[nameID];
                    int           dimData   = dimBufferTurb[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                                offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                                GPUSendBufferTurb, GPUMPIDataSendTurb);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }

                double *sendData = GPUMPIDataSendTurb + offsetNgb;
                //! get send tag from iZone to jZone
                int tagIZone = iZone + offsetTagIJ;
#ifndef CUDA_AWARE_MPI
                //! transfer compuressed data from device to host
                size_t sizeAWARE = sizeof(double) * lengthBufferTurb;
                HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendTurb, GPUMPIDataSendTurb, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
                requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(GPUMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
                countRequest++;
#else
                requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(HostMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
                countRequest++;
#endif
            }
            else if (sendOrRecv == 0)
            { //! receive data
                double *recvData = GPUMPIDataRecvTurb + offsetNgb;
                //! set recv tag from jZone to iZone
                int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(GPUMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
                countRequest++;
#else
                requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(HostMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
                countRequest++;
#endif
            } //! end if sendOrRecv
        } //! end loop of mpiID

        //! test total number of MPI operations
        if (countRequest != nTotalRequestTurb)
        {
            printf("countRequest is not equal to nTotalRequestTurb, countRequest = %d, "
                   "nTotalRequestTurb = %d\n",
                   countRequest, nTotalRequestTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerTurb, MPI_STATUSES_IGNORE);

#endif

//! test GPUMPIDataSendTurb and HostSendLargeBufferTurb
#ifdef TESTNEWORGMPI
        double *HostSendBufferTurbAWARE = new double[lengthBufferTurb];
        size_t  sizeAWARETest           = sizeof(double) * lengthBufferTurb;
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferTurbAWARE, GPUMPIDataSendTurb, sizeAWARETest, cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForSend      = iinfo->GetFaceIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                int dimData = dimBufferTurb[nameID];
                //! Get start position for variable in HostSendBufferTurbAWARE.
                TestCompressAndTransferInterpolateDataTurb(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                           offsetNgbLengthBuffer, faceIndexForSend,
                                                           HostSendBufferTurbAWARE, HostSendBufferTurb);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

        delete[] HostSendBufferTurbAWARE;
#endif
    }

    void CATSeparateInterfaceIsendIrecvTurb(InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! Get total interface size of iZone
        int nIFace = iinfo->GetNIFace();
        //! for string total number of MPI operations
        int countRequest = 0;
#ifdef MPIOVERLAP
        int countRequestFirst   = 0;
        int countRequestRemains = 0;
#endif

#ifndef ORGISENDIRECVORDER //! NOT original MPI operations' order
        //! loop all of neighbor zones for Parcel send data
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthTurb[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameTurbBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferTurb
                int offsetVar = offsetBufferTurb[nameID];
                //! get dim of variable
                int dimData = dimBufferTurb[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                            GPUSendBufferTurb, GPUMPIDataSendTurb);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

#ifndef CUDA_AWARE_MPI //! just by traditional MPI
        size_t sizeAWARE = sizeof(double) * lengthBufferTurb;
        HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendTurb, GPUMPIDataSendTurb, sizeAWARE, cudaMemcpyDeviceToHost));
#endif
#ifdef CUDA_AWARE_MPI
        //! cudaDeviceSynchronize();
        cudaStreamSynchronize(0);
#endif

#ifndef MPIOVERLAP
        //! loop all of neighbor zones only for MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthTurb[ngbID];
            //! Get initial address of send data in GPUMPIDataSendTurb
            double *sendData = GPUMPIDataSendTurb + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(GPUMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#else
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(GPUMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#else
            requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestTurb)
        {
            printf("countRequest is not equal to nTotalRequestTurb, countRequest = %d, "
                   "nTotalRequestTurb = %d\n",
                   countRequest, nTotalRequestTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerTurb, MPI_STATUSES_IGNORE);
#else

        //! loop all of neighbor zones only for MPI First operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthTurbFirst[ngbID];
            //! Get initial address of send data in GPUMPIDataSendTurb
            double *sendData = GPUMPIDataSendTurb + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurbFirst[countRequestFirst] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerTurbFirst[countRequestFirst]);
            countRequestFirst++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurbFirst[countRequestFirst] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerTurbFirst[countRequestFirst]);
            countRequestFirst++;
#endif
        }
        //! test total number of MPI operations
        if (countRequestFirst != nTotalRequestTurb)
        {
            printf("countRequestFirst is not equal to nTotalRequestTurb, countRequest = "
                   "%d, nTotalRequestTurb = %d\n",
                   countRequestFirst, nTotalRequestTurb);
            exit(1);
        }
        MPI_Waitall(countRequestFirst, requestContainerTurbFirst, MPI_STATUSES_IGNORE);

        //! loop all of neighbor zones only for MPI Remain operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthTurbRemains[ngbID];
            //! Get initial address of send data in GPUMPIDataSendTurb
            double *sendData = GPUMPIDataSendTurb + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurbRemains[countRequestRemains] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendTurb + ngbBufferLengthTurbFirst[ngbID] + offsetNgb, ngbLength * sizeof(double),
                      MPI_CHAR, ngbProcess, tagIZone, MPI_COMM_WORLD,
                      &requestContainerTurbRemains[countRequestRemains]);
            countRequestRemains++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerTurbRemains[countRequestRemains] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvTurb + ngbBufferLengthTurbFirst[ngbID] + offsetNgb, ngbLength * sizeof(double),
                      MPI_CHAR, ngbProcess, tagJZone, MPI_COMM_WORLD,
                      &requestContainerTurbRemains[countRequestRemains]);
            countRequestRemains++;
#endif
        }
        //! test total number of MPI operations
        if (countRequestRemains != nTotalRequestTurb)
        {
            printf("countRequestRemains is not equal to nTotalRequestTurb, countRequest = "
                   "%d, nTotalRequestTurb = %d\n",
                   countRequestRemains, nTotalRequestTurb);
            exit(1);
        }
//! It should be placed into Decompress NS Interface part
//! MPI_Waitall(countRequestRemains, requestContainerTurbRemains,
//! MPI_STATUSES_IGNORE);
#endif //! end ifndef MPIOVERLAP
#else
        int nTotalMPIOperations = nTotalRequestTurb;
        //! Loop for Parcel Send data
        for (int mpiID = 0; mpiID < nTotalRequestTurb; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbID[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecv[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcess[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbBufferLengthTurb[ngbID];
            //! int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            //! Parcel send data for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendTurb
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                {
                    const string &name      = (*nameTurbBuffer)[nameID];
                    int           offsetVar = offsetBufferTurb[nameID];
                    int           dimData   = dimBufferTurb[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                                offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                                GPUSendBufferTurb, GPUMPIDataSendTurb);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
            }
        } //! end loop of mpiID

#ifndef CUDA_AWARE_MPI
        //! transfer compuressed data from device to host
        size_t sizeAWARE = sizeof(double) * lengthBufferTurb;
        HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendTurb, GPUMPIDataSendTurb, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

        //! Loop for Send/Recv MPI operations
        for (int mpiID = 0; mpiID < nTotalRequestTurb; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbID[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecv[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvTurb
            int offsetNgb = offsetMPIDataTurb[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcess[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbBufferLengthTurb[ngbID];
            //! int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterfaceTurb;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                double *sendData = GPUMPIDataSendTurb + offsetNgb;
                //! get send tag from iZone to jZone
                int tagIZone = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI
                requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(GPUMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
                countRequest++;
#else
                requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(HostMPIDataSendTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
                countRequest++;
#endif
            }
            else if (sendOrRecv == 0)
            { //! receive data
                double *recvData = GPUMPIDataRecvTurb + offsetNgb;
                //! set recv tag from jZone to iZone
                int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(GPUMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
                countRequest++;
#else
                requestContainerTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(HostMPIDataRecvTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerTurb[countRequest]);
                countRequest++;
#endif
            } //! end if sendOrRecv
        } //! end loop of mpiID

        //! test total number of MPI operations
        if (countRequest != nTotalRequestTurb)
        {
            printf("countRequest is not equal to nTotalRequestTurb, countRequest = %d, "
                   "nTotalRequestTurb = %d\n",
                   countRequest, nTotalRequestTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerTurb, MPI_STATUSES_IGNORE);

#endif

//! test GPUMPIDataSendTurb and HostSendLargeBufferTurb
#ifdef TESTNEWORGMPI
        double *HostSendBufferTurbAWARE = new double[lengthBufferTurb];
        size_t  sizeAWARETest           = sizeof(double) * lengthBufferTurb;
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferTurbAWARE, GPUMPIDataSendTurb, sizeAWARETest, cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForSend      = iinfo->GetFaceIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                int dimData = dimBufferTurb[nameID];
                //! Get start position for variable in HostSendBufferTurbAWARE.
                TestCompressAndTransferInterpolateDataTurb(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                           offsetNgbLengthBuffer, faceIndexForSend,
                                                           HostSendBufferTurbAWARE, HostSendBufferTurb);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

        delete[] HostSendBufferTurbAWARE;
#endif
    }

    void CATInterfaceSendrecvTurb(InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
    }

    void CATInterfaceSendrecvNS(InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! rearange data in GPUSendBufferNS
        int iZone          = sendGlobalZone[0];
        int nTotalNgbZones = numNgbZones[0];
        int nIFace         = iinfo->GetNIFace();
        //! vector <MPI_Request> requestContainer;
        int countRequest = 0;

        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int ngbZoneID = recvGlobalZone[ngbID];
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! it may be a bug?? here or creation part
            int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int offsetNgbLengthBuffer = 0;
            int offsetFaces           = offsetFaceIndex[ngbID];
            int offsetNgb             = offsetMPIDataNS[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                const string &name      = (*nameNSBuffer)[nameID];
                int           offsetVar = offsetBufferNS[nameID];
                int           dimData   = dimBufferNS[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                            GPUSendBufferNS, GPUMPIDataSendNS);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
            //! send data
            double *sendData   = GPUMPIDataSendNS + offsetNgb;
            int     ngbProcess = recvProcess[ngbID];
            int     ngbLength  = offsetNgbLengthBuffer;
            //! int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            int tagIZone    = iZone + offsetTagIJ;
#ifndef CUDA_AWARE_MPI
            size_t sizeAWARE = sizeof(double) * lengthBufferNS;
            HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendNS, GPUMPIDataSendNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvNS + offsetNgb;
            //! double ** recvData = &GPUMPIDataRecvNS + offsetNgb;
            int tagJZone = ngbZoneID + offsetTagIJ;

            MPI_Status status;
#ifdef CUDA_AWARE_MPI
            MPI_Sendrecv(GPUMPIDataSendNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess, tagIZone,
                         GPUMPIDataRecvNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess, tagJZone, MPI_COMM_WORLD,
                         &status);
#else
            MPI_Sendrecv(HostMPIDataSendNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess, tagIZone,
                         HostMPIDataRecvNS + offsetNgb, ngbLength, MPI_DOUBLE, ngbProcess, tagJZone, MPI_COMM_WORLD,
                         &status);
#endif
            //! check MPI_Sendrecv results
            int recvElems;
            MPI_Get_count(&status, MPI_DOUBLE, &recvElems);
            if (recvElems != ngbLength)
            {
                printf("MPI_Sendrecv Error in iZone=%d, ngbZoneID = %d, ngbProcess = %d, "
                       "recvElems = %d, ngbLength=%d\n",
                       iZone, ngbZoneID, ngbProcess, recvElems, ngbLength);
                exit(1);
            }
        }

//! test GPUMPIDataSendNS and HostSendLargeBufferNS
#ifdef TESTNEWORGMPI
        double *HostSendBufferNSAWARE = new double[lengthBufferNS];
        size_t  sizeAWARETest         = sizeof(double) * lengthBufferNS;
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNSAWARE, GPUMPIDataSendNS, sizeAWARETest, cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForSend      = iinfo->GetFaceIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                int dimData = dimBufferNS[nameID];
                //! Get start position for variable in HostSendBufferNSAWARE.
                TestCompressAndTransferInterpolateDataNS(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                         offsetNgbLengthBuffer, faceIndexForSend, HostSendBufferNSAWARE,
                                                         HostSendBufferNS);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

        delete[] HostSendBufferNSAWARE;
#endif
    }

    void CATSeparateInterpointIsendIrecvNS(InterpointInformation *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        //! Get global zone index on current process
        int iZone = sendGlobalZoneForPoint[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZonesForPoint[0];
        //! Get total interface point size of iZone
        int nIPoint = iinfo->GetNumberOfInterpoints(); //! for interpoint
        //! for string total number of MPI operations
        int countRequest = 0;

#ifndef ORGISENDIRECVORDER //! NOT original MPI operations' order
        //! loop all of neighbor zones for compressing send data
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total points of jZone's interface
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetPoints = offsetPointIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendInterpointNS
            int offsetNgb = offsetMPIDataInterpointNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbInterpointBufferLengthNS[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIPointOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameNSInterpointBuffer)[nameID];
                //! get initial address of variable in GPUSendInterpointBufferNS
                int offsetVar = offsetInterpointBufferNS[nameID];
                //! get dim of variable
                int dimData = dimInterpointBufferNS[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetPoints, d_pointIndexForSend,
                                                            GPUSendInterpointBufferNS, GPUMPIDataSendInterpointNS);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
            }
        }
//! After Data Compressing, just once from Device to Host
#ifndef CUDA_AWARE_MPI //! just by traditional MPI
        size_t sizeAWARE = sizeof(double) * lengthInterpointBufferNS;
        HANDLE_API_ERR(
            cudaMemcpy(HostMPIDataSendInterpointNS, GPUMPIDataSendInterpointNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif
#ifdef CUDA_AWARE_MPI
        //! cudaDeviceSynchronize();
        cudaStreamSynchronize(0);
#endif
        //! loop all of neighbor zones for send/recv MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total points of jZone's interface
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetPoints = offsetPointIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendInterpointNS
            int offsetNgb = offsetMPIDataInterpointNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbInterpointBufferLengthNS[ngbID];
            //! Get initial address of send data in GPUMPIDataSendInterpointNS
            double *sendData = GPUMPIDataSendInterpointNS + offsetNgb;
            //! get jZone's prcess
            //! int ngbProcess = recvProcess[ngbID];
            int ngbProcess = recvProcessForPoint[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterpointNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI
            requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(GPUMPIDataSendInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagIZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
            countRequest++;
#else
            //! Add for test
            //! printf("Interpoint, iZone = %d, ngbID = %d, ngbProcess=%d\n", iZone,
            //! ngbID, ngbProcess); Add end
            requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagIZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvInterpointNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
            requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(GPUMPIDataRecvInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagJZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
            countRequest++;
#else
            requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagJZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestInterpointNS)
        {
            printf("countRequest is not equal to nTotalRequestInterpointNS, countRequest "
                   "= %d, nTotalRequestInterpointNS = %d\n",
                   countRequest, nTotalRequestInterpointNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerInterpointNS, MPI_STATUSES_IGNORE);
#else
        int nTotalMPIOperations = nTotalRequestInterpointNS;
        for (int mpiID = 0; mpiID < nTotalRequestInterpointNS; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbIDForPoint[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecvForPoint[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetPoints = offsetPointIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvInterpointNS
            int offsetNgb = offsetMPIDataInterpointNS[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcessForPoint[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbInterpointBufferLengthNS[ngbID];
            //! int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterpointNS;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendInterpointNS
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIPointOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSInterpointBuffer)[nameID];
                    int           offsetVar = offsetInterpointBufferNS[nameID];
                    int           dimData   = dimInterpointBufferNS[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor, offsetNgbLengthBuffer, offsetPoints,
                        d_pointIndexForSend, GPUSendInterpointBufferNS, GPUMPIDataSendInterpointNS);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
                }
            }
        } //! end loop of mpiID

#ifndef CUDA_AWARE_MPI
        //! transfer compuressed data from device to host
        size_t sizeAWARE = sizeof(double) * lengthInterpointBufferNS;
        HANDLE_API_ERR(
            cudaMemcpy(HostMPIDataSendInterpointNS, GPUMPIDataSendInterpointNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

        for (int mpiID = 0; mpiID < nTotalRequestInterpointNS; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbIDForPoint[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecvForPoint[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetPoints = offsetPointIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvInterpointNS
            int offsetNgb = offsetMPIDataInterpointNS[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcessForPoint[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbInterpointBufferLengthNS[ngbID];
            //! int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterpointNS;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendInterpointNS
                int offsetNgbLengthBuffer = 0;

                double *sendData = GPUMPIDataSendInterpointNS + offsetNgb;
                //! get send tag from iZone to jZone
                int tagIZone = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI
                requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(GPUMPIDataSendInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagIZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
                countRequest++;
#else
                requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(HostMPIDataSendInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagIZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
                countRequest++;
#endif
            }
            else if (sendOrRecv == 0)
            { //! receive data
                double *recvData = GPUMPIDataRecvInterpointNS + offsetNgb;
                //! set recv tag from jZone to iZone
                int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(GPUMPIDataRecvInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagJZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
                countRequest++;
#else
                requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(HostMPIDataRecvInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagJZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
                countRequest++;
#endif
            } //! end if sendOrRecv
        } //! end loop of mpiID

        //! test total number of MPI operations
        if (countRequest != nTotalRequestInterpointNS)
        {
            printf("countRequest is not equal to nTotalRequestInterpointNS, countRequest "
                   "= %d, nTotalRequestInterpointNS = %d\n",
                   countRequest, nTotalRequestInterpointNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerInterpointNS, MPI_STATUSES_IGNORE);

#endif

//! test GPUMPIDataSendInterpointNS and HostSendLargeBufferInterpointNS
#ifdef TESTNEWORGMPI
        double *HostSendInterpointBufferNSAWARE = new double[lengthInterpointBufferNS];
        size_t  sizeAWARETest                   = sizeof(double) * lengthInterpointBufferNS;
        HANDLE_API_ERR(cudaMemcpy(HostSendInterpointBufferNSAWARE, GPUMPIDataSendInterpointNS, sizeAWARETest,
                                  cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIPointOfNeighbor     = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            int *pointIndexForSend     = iinfo->GetPointIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
            {
                int dimData = dimInterpointBufferNS[nameID];
                //! Get start position for variable in HostSendBufferNSAWARE.
                TestCompressAndTransferInterpointDataNS(iZone, ngbZoneID, nameID, ngbID, nIPoint, nIPointOfNeighbor,
                                                        offsetNgbLengthBuffer, pointIndexForSend,
                                                        HostSendInterpointBufferNSAWARE, HostSendInterpointBufferNS);

                offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
            }
        }

        delete[] HostSendInterpointBufferNSAWARE;
#endif
    }

    void CATSeparateInterpointIsendIrecvTurb(InterpointInformation *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! Get global zone index on current process
        int iZone = sendGlobalZoneForPoint[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZonesForPoint[0];
        //! Get total interface point size of iZone
        int nIPoint = iinfo->GetNumberOfInterpoints(); //! for interpoint
        //! for string total number of MPI operations
        int countRequest = 0;

#ifndef ORGISENDIRECVORDER //! NOT original MPI operations' order
        //! loop all of neighbor zones only for compress send data
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total points of jZone's interface
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetPoints = offsetPointIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendInterpointTurb
            int offsetNgb = offsetMPIDataInterpointTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbInterpointBufferLengthTurb[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIPointOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameTurbInterpointBuffer)[nameID];
                //! get initial address of variable in GPUSendInterpointBufferTurb
                int offsetVar = offsetInterpointBufferTurb[nameID];
                //! get dim of variable
                int dimData = dimInterpointBufferTurb[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetPoints, d_pointIndexForSend,
                                                            GPUSendInterpointBufferTurb, GPUMPIDataSendInterpointTurb);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
            }
        }

#ifndef CUDA_AWARE_MPI //! just by traditional MPI
        size_t sizeAWARE = sizeof(double) * lengthInterpointBufferTurb;
        HANDLE_API_ERR(
            cudaMemcpy(HostMPIDataSendInterpointTurb, GPUMPIDataSendInterpointTurb, sizeAWARE, cudaMemcpyDeviceToHost));
#endif
#ifdef CUDA_AWARE_MPI
        //! cudaDeviceSynchronize();
        cudaStreamSynchronize(0);
#endif

        //! loop all of neighbor zones only for send/recv MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total points of jZone's interface
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetPoints = offsetPointIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendInterpointTurb
            int offsetNgb = offsetMPIDataInterpointTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbInterpointBufferLengthTurb[ngbID];
            //! Get initial address of send data in GPUMPIDataSendInterpointTurb
            double *sendData = GPUMPIDataSendInterpointTurb + offsetNgb;
            //! get jZone's prcess
            //! int ngbProcess = recvProcess[ngbID];
            int ngbProcess = recvProcessForPoint[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterpointTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI
            requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(GPUMPIDataSendInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagIZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
            countRequest++;
#else
            requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagIZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvInterpointTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
            requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(GPUMPIDataRecvInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagJZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
            countRequest++;
#else
            requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagJZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestInterpointTurb)
        {
            printf("countRequest is not equal to nTotalRequestInterpointTurb, "
                   "countRequest = %d, nTotalRequestInterpointTurb = %d\n",
                   countRequest, nTotalRequestInterpointTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerInterpointTurb, MPI_STATUSES_IGNORE);
#else
        int nTotalMPIOperations = nTotalRequestInterpointTurb;
        for (int mpiID = 0; mpiID < nTotalRequestInterpointTurb; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbIDForPoint[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecvForPoint[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetPoints = offsetPointIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvInterpointTurb
            int offsetNgb = offsetMPIDataInterpointTurb[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcessForPoint[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbInterpointBufferLengthTurb[ngbID];
            //! int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterpointTurb;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendInterpointTurb
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIPointOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
                {
                    const string &name      = (*nameTurbInterpointBuffer)[nameID];
                    int           offsetVar = offsetInterpointBufferTurb[nameID];
                    int           dimData   = dimInterpointBufferTurb[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor, offsetNgbLengthBuffer, offsetPoints,
                        d_pointIndexForSend, GPUSendInterpointBufferTurb, GPUMPIDataSendInterpointTurb);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
                }
            } //! end if sendOrRecv
        } //! end loop of mpiID

#ifndef CUDA_AWARE_MPI
        //! transfer compuressed data from device to host
        size_t sizeAWARE = sizeof(double) * lengthInterpointBufferTurb;
        HANDLE_API_ERR(
            cudaMemcpy(HostMPIDataSendInterpointTurb, GPUMPIDataSendInterpointTurb, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

        for (int mpiID = 0; mpiID < nTotalRequestInterpointTurb; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbIDForPoint[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecvForPoint[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetPoints = offsetPointIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvInterpointTurb
            int offsetNgb = offsetMPIDataInterpointTurb[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcessForPoint[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbInterpointBufferLengthTurb[ngbID];
            //! int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterpointTurb;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendInterpointTurb
                double *sendData = GPUMPIDataSendInterpointTurb + offsetNgb;
                //! get send tag from iZone to jZone
                int tagIZone = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI
                requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(GPUMPIDataSendInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagIZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
                countRequest++;
#else
                requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(HostMPIDataSendInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagIZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
                countRequest++;
#endif
            }
            else if (sendOrRecv == 0)
            { //! receive data
                double *recvData = GPUMPIDataRecvInterpointTurb + offsetNgb;
                //! set recv tag from jZone to iZone
                int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(GPUMPIDataRecvInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagJZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
                countRequest++;
#else
                requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(HostMPIDataRecvInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagJZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
                countRequest++;
#endif
            } //! end if sendOrRecv
        } //! end loop of mpiID
        //! test total number of MPI operations
        if (countRequest != nTotalRequestInterpointTurb)
        {
            printf("countRequest is not equal to nTotalRequestInterpointTurb, "
                   "countRequest = %d, nTotalRequestInterpointTurb = %d\n",
                   countRequest, nTotalRequestInterpointTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerInterpointTurb, MPI_STATUSES_IGNORE);

#endif

//! test GPUMPIDataSendInterpointTurb and HostSendLargeBufferInterpointTurb
#ifdef TESTNEWORGMPI
        double *HostSendInterpointBufferTurbAWARE = new double[lengthInterpointBufferTurb];
        size_t  sizeAWARETest                     = sizeof(double) * lengthInterpointBufferTurb;
        HANDLE_API_ERR(cudaMemcpy(HostSendInterpointBufferTurbAWARE, GPUMPIDataSendInterpointTurb, sizeAWARETest,
                                  cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIPointOfNeighbor     = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            int *pointIndexForSend     = iinfo->GetPointIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
            {
                int dimData = dimInterpointBufferTurb[nameID];
                //! Get start position for variable in HostSendBufferTurbAWARE.
                TestCompressAndTransferInterpointDataTurb(
                    iZone, ngbZoneID, nameID, ngbID, nIPoint, nIPointOfNeighbor, offsetNgbLengthBuffer,
                    pointIndexForSend, HostSendInterpointBufferTurbAWARE, HostSendInterpointBufferTurb);

                offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
            }
        }

        delete[] HostSendInterpointBufferTurbAWARE;
#endif
    }

    void CATOverlapInterpointIsendIrecvTurb(InterpointInformation *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! Get global zone index on current process
        int iZone = sendGlobalZoneForPoint[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZonesForPoint[0];
        //! Get total interface point size of iZone
        int nIPoint = iinfo->GetNumberOfInterpoints(); //! for interpoint
        //! for string total number of MPI operations
        int countRequest = 0;

#ifndef ORGISENDIRECVORDER //! NOT original MPI operations' order
        //! loop all of neighbor zones for MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total points of jZone's interface
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetPoints = offsetPointIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendInterpointTurb
            int offsetNgb = offsetMPIDataInterpointTurb[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbInterpointBufferLengthTurb[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIPointOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameTurbInterpointBuffer)[nameID];
                //! get initial address of variable in GPUSendInterpointBufferTurb
                int offsetVar = offsetInterpointBufferTurb[nameID];
                //! get dim of variable
                int dimData = dimInterpointBufferTurb[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetPoints, d_pointIndexForSend,
                                                            GPUSendInterpointBufferTurb, GPUMPIDataSendInterpointTurb);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
            }
            //! Get initial address of send data in GPUMPIDataSendInterpointTurb
            double *sendData = GPUMPIDataSendInterpointTurb + offsetNgb;
            //! get jZone's prcess
            //! int ngbProcess = recvProcess[ngbID];
            int ngbProcess = recvProcessForPoint[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterpointTurb;
            int tagIZone    = iZone + offsetTagIJ;

#ifndef CUDA_AWARE_MPI //! just by traditional MPI
            size_t sizeAWARE = sizeof(double) * lengthInterpointBufferTurb;
            HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendInterpointTurb, GPUMPIDataSendInterpointTurb, sizeAWARE,
                                      cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
            requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(GPUMPIDataSendInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagIZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
            countRequest++;
#else
            requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagIZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvInterpointTurb + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
            requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(GPUMPIDataRecvInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagJZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
            countRequest++;
#else
            requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagJZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestInterpointTurb)
        {
            printf("countRequest is not equal to nTotalRequestInterpointTurb, "
                   "countRequest = %d, nTotalRequestInterpointTurb = %d\n",
                   countRequest, nTotalRequestInterpointTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerInterpointTurb, MPI_STATUSES_IGNORE);
#else
        int nTotalMPIOperations = nTotalRequestInterpointTurb;
        for (int mpiID = 0; mpiID < nTotalRequestInterpointTurb; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbIDForPoint[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecvForPoint[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetPoints = offsetPointIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvInterpointTurb
            int offsetNgb = offsetMPIDataInterpointTurb[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcessForPoint[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbInterpointBufferLengthTurb[ngbID];
            //! int offsetTagIJ = 10000;
            int offsetTagIJ = offsetTagInterpointTurb;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendInterpointTurb
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIPointOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
                {
                    const string &name      = (*nameTurbInterpointBuffer)[nameID];
                    int           offsetVar = offsetInterpointBufferTurb[nameID];
                    int           dimData   = dimInterpointBufferTurb[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor, offsetNgbLengthBuffer, offsetPoints,
                        d_pointIndexForSend, GPUSendInterpointBufferTurb, GPUMPIDataSendInterpointTurb);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
                }

                double *sendData = GPUMPIDataSendInterpointTurb + offsetNgb;
                //! get send tag from iZone to jZone
                int tagIZone = iZone + offsetTagIJ;
#ifndef CUDA_AWARE_MPI
                //! transfer compuressed data from device to host
                size_t sizeAWARE = sizeof(double) * lengthInterpointBufferTurb;
                HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendInterpointTurb, GPUMPIDataSendInterpointTurb, sizeAWARE,
                                          cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
                requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(GPUMPIDataSendInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagIZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
                countRequest++;
#else
                requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(HostMPIDataSendInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagIZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
                countRequest++;
#endif
            }
            else if (sendOrRecv == 0)
            { //! receive data
                double *recvData = GPUMPIDataRecvInterpointTurb + offsetNgb;
                //! set recv tag from jZone to iZone
                int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(GPUMPIDataRecvInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagJZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
                countRequest++;
#else
                requestContainerInterpointTurb[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(HostMPIDataRecvInterpointTurb + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagJZone, MPI_COMM_WORLD, &requestContainerInterpointTurb[countRequest]);
                countRequest++;
#endif
            } //! end if sendOrRecv
        } //! end loop of mpiID

        //! test total number of MPI operations
        if (countRequest != nTotalRequestInterpointTurb)
        {
            printf("countRequest is not equal to nTotalRequestInterpointTurb, "
                   "countRequest = %d, nTotalRequestInterpointTurb = %d\n",
                   countRequest, nTotalRequestInterpointTurb);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerInterpointTurb, MPI_STATUSES_IGNORE);

#endif

//! test GPUMPIDataSendInterpointTurb and HostSendLargeBufferInterpointTurb
#ifdef TESTNEWORGMPI
        double *HostSendInterpointBufferTurbAWARE = new double[lengthInterpointBufferTurb];
        size_t  sizeAWARETest                     = sizeof(double) * lengthInterpointBufferTurb;
        HANDLE_API_ERR(cudaMemcpy(HostSendInterpointBufferTurbAWARE, GPUMPIDataSendInterpointTurb, sizeAWARETest,
                                  cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZoneForPoint[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIPointOfNeighbor     = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            int *pointIndexForSend     = iinfo->GetPointIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
            {
                int dimData = dimInterpointBufferTurb[nameID];
                //! Get start position for variable in HostSendBufferTurbAWARE.
                TestCompressAndTransferInterpointDataTurb(
                    iZone, ngbZoneID, nameID, ngbID, nIPoint, nIPointOfNeighbor, offsetNgbLengthBuffer,
                    pointIndexForSend, HostSendInterpointBufferTurbAWARE, HostSendInterpointBufferTurb);

                offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
            }
        }

        delete[] HostSendInterpointBufferTurbAWARE;
#endif
    }
    void CATOverlapInterpointIsendIrecvNS(InterpointInformation *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! Get global zone index on current process
        int iZone = sendGlobalZoneForPoint[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZonesForPoint[0];
        //! Get total interface point size of iZone
        int nIPoint = iinfo->GetNumberOfInterpoints(); //! for interpoint
        //! for string total number of MPI operations
        int countRequest = 0;

#ifndef ORGISENDIRECVORDER //! NOT original MPI operations' order
        //! loop all of neighbor zones for MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total points of jZone's interface
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetPoints = offsetPointIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendInterpointNS
            int offsetNgb = offsetMPIDataInterpointNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbInterpointBufferLengthNS[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIPointOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameNSInterpointBuffer)[nameID];
                //! get initial address of variable in GPUSendInterpointBufferNS
                int offsetVar = offsetInterpointBufferNS[nameID];
                //! get dim of variable
                int dimData = dimInterpointBufferNS[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetPoints, d_pointIndexForSend,
                                                            GPUSendInterpointBufferNS, GPUMPIDataSendInterpointNS);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
            }
            //! Get initial address of send data in GPUMPIDataSendInterpointNS
            double *sendData = GPUMPIDataSendInterpointNS + offsetNgb;
            //! get jZone's prcess
            //! int ngbProcess = recvProcess[ngbID];
            int ngbProcess = recvProcessForPoint[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 500;
            int offsetTagIJ = offsetTagInterpointNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifndef CUDA_AWARE_MPI //! just by traditional MPI
            size_t sizeAWARE = sizeof(double) * lengthInterpointBufferNS;
            HANDLE_API_ERR(
                cudaMemcpy(HostMPIDataSendInterpointNS, GPUMPIDataSendInterpointNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
            requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(GPUMPIDataSendInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagIZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
            countRequest++;
#else
            requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagIZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvInterpointNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
            requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(GPUMPIDataRecvInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagJZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
            countRequest++;
#else
            requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                      tagJZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestInterpointNS)
        {
            printf("countRequest is not equal to nTotalRequestInterpointNS, countRequest "
                   "= %d, nTotalRequestInterpointNS = %d\n",
                   countRequest, nTotalRequestInterpointNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerInterpointNS, MPI_STATUSES_IGNORE);
#else
        int nTotalMPIOperations = nTotalRequestInterpointNS;
        for (int mpiID = 0; mpiID < nTotalRequestInterpointNS; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbIDForPoint[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecvForPoint[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZoneForPoint[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIPointOfNeighbor = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetPoints = offsetPointIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvInterpointNS
            int offsetNgb = offsetMPIDataInterpointNS[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcessForPoint[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbInterpointBufferLengthNS[ngbID];
            //! int offsetTagIJ = 500;
            int offsetTagIJ = offsetTagInterpointNS;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendInterpointNS
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIPointOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSInterpointBuffer)[nameID];
                    int           offsetVar = offsetInterpointBufferNS[nameID];
                    int           dimData   = dimInterpointBufferNS[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(
                        offsetVar, dimData, nIPoint, offsetNgb, nIPointOfNeighbor, offsetNgbLengthBuffer, offsetPoints,
                        d_pointIndexForSend, GPUSendInterpointBufferNS, GPUMPIDataSendInterpointNS);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
                }

                double *sendData = GPUMPIDataSendInterpointNS + offsetNgb;
                //! get send tag from iZone to jZone
                int tagIZone = iZone + offsetTagIJ;
#ifndef CUDA_AWARE_MPI
                //! transfer compuressed data from device to host
                size_t sizeAWARE = sizeof(double) * lengthInterpointBufferNS;
                HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendInterpointNS, GPUMPIDataSendInterpointNS, sizeAWARE,
                                          cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
                requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(GPUMPIDataSendInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagIZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
                countRequest++;
#else
                requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(HostMPIDataSendInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagIZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
                countRequest++;
#endif
            }
            else if (sendOrRecv == 0)
            { //! receive data
                double *recvData = GPUMPIDataRecvInterpointNS + offsetNgb;
                //! set recv tag from jZone to iZone
                int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(GPUMPIDataRecvInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagJZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
                countRequest++;
#else
                requestContainerInterpointNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(HostMPIDataRecvInterpointNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess,
                          tagJZone, MPI_COMM_WORLD, &requestContainerInterpointNS[countRequest]);
                countRequest++;
#endif
            } //! end if sendOrRecv
        } //! end loop of mpiID

        //! test total number of MPI operations
        if (countRequest != nTotalRequestInterpointNS)
        {
            printf("countRequest is not equal to nTotalRequestInterpointNS, countRequest "
                   "= %d, nTotalRequestInterpointNS = %d\n",
                   countRequest, nTotalRequestInterpointNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerInterpointNS, MPI_STATUSES_IGNORE);

#endif

//! test GPUMPIDataSendInterpointNS and HostSendLargeBufferInterpointNS
#ifdef TESTNEWORGMPI
        double *HostSendInterpointBufferNSAWARE = new double[lengthInterpointBufferNS];
        size_t  sizeAWARETest                   = sizeof(double) * lengthInterpointBufferNS;
        HANDLE_API_ERR(cudaMemcpy(HostSendInterpointBufferNSAWARE, GPUMPIDataSendInterpointNS, sizeAWARETest,
                                  cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZoneForPoint[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIPointOfNeighbor     = iinfo->GetNumberOfInterpointsForNeighbor(iNeighbor);
            int *pointIndexForSend     = iinfo->GetPointIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
            {
                int dimData = dimInterpointBufferNS[nameID];
                //! Get start position for variable in HostSendBufferNSAWARE.
                TestCompressAndTransferInterpointDataNS(iZone, ngbZoneID, nameID, ngbID, nIPoint, nIPointOfNeighbor,
                                                        offsetNgbLengthBuffer, pointIndexForSend,
                                                        HostSendInterpointBufferNSAWARE, HostSendInterpointBufferNS);

                offsetNgbLengthBuffer += dimData * nIPointOfNeighbor;
            }
        }

        delete[] HostSendInterpointBufferNSAWARE;
#endif
    }

    void CATOverlapInterfaceIsendIrecvNS(InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! Get total interface size of iZone
        int nIFace = iinfo->GetNIFace();
        //! for string total number of MPI operations
        int countRequest = 0;

#ifndef ORGISENDIRECVORDER //! NOT original MPI operations' order
        //! loop all of neighbor zones for MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthNS[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameNSBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferNS
                int offsetVar = offsetBufferNS[nameID];
                //! get dim of variable
                int dimData = dimBufferNS[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                            GPUSendBufferNS, GPUMPIDataSendNS);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
            //! Get initial address of send data in GPUMPIDataSendNS
            double *sendData = GPUMPIDataSendNS + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifndef CUDA_AWARE_MPI //! just by traditional MPI
            size_t sizeAWARE = sizeof(double) * lengthBufferNS;
            HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendNS, GPUMPIDataSendNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(GPUMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#else
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(GPUMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#else
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestNS)
        {
            printf("countRequest is not equal to nTotalRequestNS, countRequest = %d, "
                   "nTotalRequestNS = %d\n",
                   countRequest, nTotalRequestNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerNS, MPI_STATUSES_IGNORE);
#else
        int nTotalMPIOperations = nTotalRequestNS;
        for (int mpiID = 0; mpiID < nTotalRequestNS; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbID[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecv[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcess[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbBufferLengthNS[ngbID];
            //! int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendNS
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSBuffer)[nameID];
                    int           offsetVar = offsetBufferNS[nameID];
                    int           dimData   = dimBufferNS[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                                offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                                GPUSendBufferNS, GPUMPIDataSendNS);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }

                double *sendData = GPUMPIDataSendNS + offsetNgb;
                //! get send tag from iZone to jZone
                int tagIZone = iZone + offsetTagIJ;
#ifndef CUDA_AWARE_MPI
                //! transfer compuressed data from device to host
                size_t sizeAWARE = sizeof(double) * lengthBufferNS;
                HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendNS, GPUMPIDataSendNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

#ifdef CUDA_AWARE_MPI
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(GPUMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#else
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#endif
            }
            else if (sendOrRecv == 0)
            { //! receive data
                double *recvData = GPUMPIDataRecvNS + offsetNgb;
                //! set recv tag from jZone to iZone
                int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(GPUMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#else
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#endif
            } //! end if sendOrRecv
        } //! end loop of mpiID

        //! test total number of MPI operations
        if (countRequest != nTotalRequestNS)
        {
            printf("countRequest is not equal to nTotalRequestNS, countRequest = %d, "
                   "nTotalRequestNS = %d\n",
                   countRequest, nTotalRequestNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerNS, MPI_STATUSES_IGNORE);

#endif

//! test GPUMPIDataSendNS and HostSendLargeBufferNS
#ifdef TESTNEWORGMPI
        double *HostSendBufferNSAWARE = new double[lengthBufferNS];
        size_t  sizeAWARETest         = sizeof(double) * lengthBufferNS;
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNSAWARE, GPUMPIDataSendNS, sizeAWARETest, cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForSend      = iinfo->GetFaceIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                int dimData = dimBufferNS[nameID];
                //! Get start position for variable in HostSendBufferNSAWARE.
                TestCompressAndTransferInterpolateDataNS(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                         offsetNgbLengthBuffer, faceIndexForSend, HostSendBufferNSAWARE,
                                                         HostSendBufferNS);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

        delete[] HostSendBufferNSAWARE;
#endif
    }

    void CATSeparateInterfaceIsendIrecvNS(InterfaceInfo *iinfo)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! Get global zone index on current process
        int iZone = sendGlobalZone[0];
        //! Get total number of neighbor zone of iZone
        int nTotalNgbZones = numNgbZones[0];
        //! Get total interface size of iZone
        int nIFace = iinfo->GetNIFace();
        //! for string total number of MPI operations
        int countRequest = 0;
#ifdef MPIOVERLAP
        int countRequestFirst   = 0;
        int countRequestRemains = 0;
#endif

#ifndef ORGISENDIRECVORDER //! NOT original MPI operations' order
        //! loop all of neighbor zones for compressing data
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
            //! ngbID=iNeighbor
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of jZone's interface
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get offset of each variables in buffer of jZone
            int offsetNgbLengthBuffer = 0;
            //! get initial address of jZone in d_send
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthNS[ngbID];

            int gridSize       = 1;
            int blockSize      = 1;
            int loopLen        = nIFaceOfNeighbor;
            int blockSizeLimit = 0;
            KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                //! get name of varialbe
                const string &name = (*nameNSBuffer)[nameID];
                //! get initial address of variable in GPUSendBufferNS
                int offsetVar = offsetBufferNS[nameID];
                //! get dim of variable
                int dimData = dimBufferNS[nameID];
                GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                            offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                            GPUSendBufferNS, GPUMPIDataSendNS);
#ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
#endif
                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

#ifndef CUDA_AWARE_MPI //! just by traditional MPI
        size_t sizeAWARE = sizeof(double) * lengthBufferNS;
        HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendNS, GPUMPIDataSendNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif
#ifdef CUDA_AWARE_MPI //! guarantee GPU data compress finished
        //! cudaDeviceSynchronize();
        cudaStreamSynchronize(0);
#endif
#ifndef MPIOVERLAP
        //! loop all of neighbor zones for MPI operations
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get initial address of jZone in d_send
            //! int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthNS[ngbID];
            //! Get initial address of send data in GPUMPIDataSendNS
            double *sendData = GPUMPIDataSendNS + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(GPUMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#else
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(GPUMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#else
            requestContainerNS[countRequest] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerNS[countRequest]);
            countRequest++;
#endif
        }
        //! test total number of MPI operations
        if (countRequest != nTotalRequestNS)
        {
            printf("countRequest is not equal to nTotalRequestNS, countRequest = %d, "
                   "nTotalRequestNS = %d\n",
                   countRequest, nTotalRequestNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerNS, MPI_STATUSES_IGNORE);
#else
        //! loop all of neighbor zones for MPI operations on First NS buffer
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get initial address of jZone in d_send
            //! int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthNSFirst[ngbID];
            //! Get initial address of send data in GPUMPIDataSendNS
            double *sendData = GPUMPIDataSendNS + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerNSFirst[countRequestFirst] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                      MPI_COMM_WORLD, &requestContainerNSFirst[countRequestFirst]);
            countRequestFirst++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerNSFirst[countRequestFirst] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                      MPI_COMM_WORLD, &requestContainerNSFirst[countRequestFirst]);
            countRequestFirst++;
#endif
        }
        //! test total number of MPI operations
        if (countRequestFirst != nTotalRequestNS)
        {
            printf("countRequestFirst is not equal to nTotalRequestNS, countRequest = %d, "
                   "nTotalRequestNS = %d\n",
                   countRequestFirst, nTotalRequestNS);
            exit(1);
        }
        MPI_Waitall(countRequestFirst, requestContainerNSFirst, MPI_STATUSES_IGNORE);

        //! loop all of neighbor zones for MPI operations on Remain NS buffer
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            //! get neighbor zone's global index (jZone)
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get initial address of jZone in d_send
            //! int offsetFaces = offsetFaceIndex[ngbID];
            //! get inital address of jZone in GPUMPIDataSendNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get the total data size for send and recv MPI ops
            int ngbLength = ngbBufferLengthNSRemains[ngbID];
            //! Get initial address of send data in GPUMPIDataSendNS
            double *sendData = GPUMPIDataSendNS + offsetNgb;
            //! get jZone's prcess
            int ngbProcess = recvProcess[ngbID];
            //! set send tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag. int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            int tagIZone    = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI

#else
            requestContainerNSRemains[countRequestRemains] = MPI_REQUEST_NULL;
            MPI_Isend(HostMPIDataSendNS + ngbBufferLengthNSFirst[ngbID] + offsetNgb, ngbLength * sizeof(double),
                      MPI_CHAR, ngbProcess, tagIZone, MPI_COMM_WORLD, &requestContainerNSRemains[countRequestRemains]);
            countRequestRemains++;
#endif
            //! recv data
            double *recvData = GPUMPIDataRecvNS + offsetNgb;
            //! set recv tag. offsetTagIJ is used for difference from original PHengLEI
            //! tag.
            int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI

#else
            requestContainerNSRemains[countRequestRemains] = MPI_REQUEST_NULL;
            MPI_Irecv(HostMPIDataRecvNS + ngbBufferLengthNSFirst[ngbID] + offsetNgb, ngbLength * sizeof(double),
                      MPI_CHAR, ngbProcess, tagJZone, MPI_COMM_WORLD, &requestContainerNSRemains[countRequestRemains]);
            countRequestRemains++;
#endif
        }
        //! test total number of MPI operations
        if (countRequestRemains != nTotalRequestNS)
        {
            printf("countRequestRemains is not equal to nTotalRequestNS, countRequest = "
                   "%d, nTotalRequestNS = %d\n",
                   countRequestRemains, nTotalRequestNS);
            exit(1);
        }
//! It should be placed on Decompress of Turb interface data.
//! MPI_Waitall(countRequestRemains, requestContainerNSRemains,
//! MPI_STATUSES_IGNORE);
#endif //! end ifndef MPIOVERLAP
#else
        int nTotalMPIOperations = nTotalRequestNS;
        //! just compressing data into GPUMPIDataSendNS
        for (int mpiID = 0; mpiID < nTotalRequestNS; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbID[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecv[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get local zone ID of neighbor zone from 0 to nTotalNgbZones-1, in fact
            //! iNeighbor=ngbID
            int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
            //! get total faces of iZone and jZone
            int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            //! get initial address of jZone in faceSendIndex
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcess[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbBufferLengthNS[ngbID];
            //! int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! compress data into GPUMPIDataSendNS
                int offsetNgbLengthBuffer = 0;

                int gridSize       = 1;
                int blockSize      = 1;
                int loopLen        = nIFaceOfNeighbor;
                int blockSizeLimit = 0;
                KernelLaunchPara((void *)GPUMPIDataCompress, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
                ReportKernelPara((void *)GPUMPIDataCompress, 0, gridSize, blockSize);
#endif
                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                {
                    const string &name      = (*nameNSBuffer)[nameID];
                    int           offsetVar = offsetBufferNS[nameID];
                    int           dimData   = dimBufferNS[nameID];
                    GPUMPIDataCompress<<<gridSize, blockSize>>>(offsetVar, dimData, nIFace, offsetNgb, nIFaceOfNeighbor,
                                                                offsetNgbLengthBuffer, offsetFaces, d_faceIndexForSend,
                                                                GPUSendBufferNS, GPUMPIDataSendNS);
#ifdef KERNELLAUNCHTEST
                    HANDLE_KERNEL_ERR();
#endif
                    offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                }
            }
        } //! end loop of mpiID

#ifndef CUDA_AWARE_MPI
        //! transfer compuressed data from device to host
        size_t sizeAWARE = sizeof(double) * lengthBufferNS;
        HANDLE_API_ERR(cudaMemcpy(HostMPIDataSendNS, GPUMPIDataSendNS, sizeAWARE, cudaMemcpyDeviceToHost));
#endif

        //! just for sending and recving MPI operations
        for (int mpiID = 0; mpiID < nTotalRequestNS; mpiID++)
        {
            //! get neighbor local label zone ID
            int ngbID = sendRecvNgbID[mpiID]; //! from 0 to nTotalNgbZones-1
            //! judge send or recv operation
            int sendOrRecv = isSendOrRecv[mpiID];
            //! get global zone ID (jZone) of neighbor zone
            int ngbZoneID = recvGlobalZone[ngbID];
            //! get initial address of jZone in faceSendIndex
            int offsetFaces = offsetFaceIndex[ngbID];
            //! get initial address of jZone in GPUMPIDataSend/RecvNS
            int offsetNgb = offsetMPIDataNS[ngbID];
            //! get process ID of jZone
            int ngbProcess = recvProcess[ngbID];
            //! get length of send/recv data for jZone
            int ngbLength = ngbBufferLengthNS[ngbID];
            //! int offsetTagIJ = 1000;
            int offsetTagIJ = offsetTagInterfaceNS;
            //! for MPI_Isend operation
            if (sendOrRecv == 1)
            {
                //! get send tag from iZone to jZone
                int tagIZone = iZone + offsetTagIJ;

#ifdef CUDA_AWARE_MPI
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(GPUMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#else
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Isend(HostMPIDataSendNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagIZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#endif
            }
            else if (sendOrRecv == 0)
            { //! receive data
                double *recvData = GPUMPIDataRecvNS + offsetNgb;
                //! set recv tag from jZone to iZone
                int tagJZone = ngbZoneID + offsetTagIJ;
#ifdef CUDA_AWARE_MPI
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(GPUMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#else
                requestContainerNS[countRequest] = MPI_REQUEST_NULL;
                MPI_Irecv(HostMPIDataRecvNS + offsetNgb, ngbLength * sizeof(double), MPI_CHAR, ngbProcess, tagJZone,
                          MPI_COMM_WORLD, &requestContainerNS[countRequest]);
                countRequest++;
#endif
            } //! end if sendOrRecv
        } //! end loop of mpiID

        //! test total number of MPI operations
        if (countRequest != nTotalRequestNS)
        {
            printf("countRequest is not equal to nTotalRequestNS, countRequest = %d, "
                   "nTotalRequestNS = %d\n",
                   countRequest, nTotalRequestNS);
            exit(1);
        }
        MPI_Waitall(countRequest, requestContainerNS, MPI_STATUSES_IGNORE);

#endif

//! test GPUMPIDataSendNS and HostSendLargeBufferNS
#ifdef TESTNEWORGMPI
        double *HostSendBufferNSAWARE = new double[lengthBufferNS];
        size_t  sizeAWARETest         = sizeof(double) * lengthBufferNS;
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNSAWARE, GPUMPIDataSendNS, sizeAWARETest, cudaMemcpyDeviceToHost));
        for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
        {
            int  ngbZoneID             = recvGlobalZone[ngbID];
            int  iNeighbor             = iinfo->FindIthNeighbor(ngbZoneID);
            int  nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
            int *faceIndexForSend      = iinfo->GetFaceIndexForSend(iNeighbor);
            int  offsetNgbLengthBuffer = 0;

            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                int dimData = dimBufferNS[nameID];
                //! Get start position for variable in HostSendBufferNSAWARE.
                TestCompressAndTransferInterpolateDataNS(iZone, ngbZoneID, nameID, ngbID, nIFace, nIFaceOfNeighbor,
                                                         offsetNgbLengthBuffer, faceIndexForSend, HostSendBufferNSAWARE,
                                                         HostSendBufferNS);

                offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
            }
        }

        delete[] HostSendBufferNSAWARE;
#endif
    }

    void CallCUDAStreamSynchronize(cudaStream_t streamProcess) { HANDLE_API_ERR(cudaStreamSynchronize(streamProcess)); }

    void CallDeviceUploadInterfaceValueHostLargeBuffer(const int nTotal, const int neqn, const int nIFace,
                                                       const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        int     solverID;
        int     nameNSID, nameTurbID;
        RFloat *deviceFieldVar;
        double *deviceSendBuffer;
        RFloat *HostSendLargeBuffer;
        int     dimBuffer;
        int     offsetBuffer;

        if ("q" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_q_ns;
            deviceSendBuffer = d_fg_send_q;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("t" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_t_proxy;
            deviceSendBuffer = d_fg_send_t;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("turb::q" == name)
        {
            solverID         = 1;
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            dimBuffer        = dimBufferTurb[nameTurbID];
            deviceFieldVar   = d_q_turb_proxy;
            deviceSendBuffer = d_fg_send_qTurb;
            //! HostSendLargeBuffer = HostSendBufferTurb;
            HostSendLargeBuffer = HostSendBufferTurb + offsetBuffer;
        }
        else if ("dqdx" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_dqdx_proxy;
            deviceSendBuffer = d_fg_send_dqdx;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("dqdy" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_dqdy_proxy;
            deviceSendBuffer = d_fg_send_dqdy;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("dqdz" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_dqdz_proxy;
            deviceSendBuffer = d_fg_send_dqdz;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("limit" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_limit;
            deviceSendBuffer = d_fg_send_limit;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("dtdx" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_dtdx_proxy;
            deviceSendBuffer = d_fg_send_dtdx;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("dtdy" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_dtdy_proxy;
            deviceSendBuffer = d_fg_send_dtdy;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("dtdz" == name)
        {
            solverID         = 0;
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            dimBuffer        = dimBufferNS[nameNSID];
            deviceFieldVar   = d_dtdz_proxy;
            deviceSendBuffer = d_fg_send_dtdz;
            //! HostSendLargeBuffer = HostSendBufferNS;
            HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
        }
        else if ("turb::dqdx" == name)
        {
            solverID         = 1;
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            dimBuffer        = dimBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdx_proxy;
            deviceSendBuffer = d_fg_send_dqTurbdx;
            //! HostSendLargeBuffer = HostSendBufferTurb;
            HostSendLargeBuffer = HostSendBufferTurb + offsetBuffer;
        }
        else if ("turb::dqdy" == name)
        {
            solverID         = 1;
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            dimBuffer        = dimBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdy_proxy;
            deviceSendBuffer = d_fg_send_dqTurbdy;
            //! HostSendLargeBuffer = HostSendBufferTurb;
            HostSendLargeBuffer = HostSendBufferTurb + offsetBuffer;
        }
        else if ("turb::dqdz" == name)
        {
            solverID         = 1;
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            dimBuffer        = dimBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdz_proxy;
            deviceSendBuffer = d_fg_send_dqTurbdz;
            //! HostSendLargeBuffer = HostSendBufferTurb;
            HostSendLargeBuffer = HostSendBufferTurb + offsetBuffer;
        }
        else
        {
            cout << "Error:" << name << " is not included in CallGPUUploadInterfaceValue" << endl;
            exit(1);
        }
        //! From device varaible to device send buffer
        CallGPUFieldVarToIFVar(deviceSendBuffer, deviceFieldVar, nTotal, neqn, nIFace);
        //! record and wait for D2D data transfer
        if (solverID == 0)
        { //! NS equation
            //! record finish of DToD data transfer
            HANDLE_API_ERR(cudaEventRecord(downloadDToDNS_done, 0));
            //! wait for event downlaodDTODNS_done
            HANDLE_API_ERR(cudaStreamWaitEvent(dataTransferDToH, downloadDToDNS_done, 0));
        }
        else if (solverID == 1)
        { //! Turbulence equation
            HANDLE_API_ERR(cudaEventRecord(downloadDToDTurb_done, 0));
            HANDLE_API_ERR(cudaStreamWaitEvent(dataTransferDToH, downloadDToDTurb_done, 0));
        }
        //! From device send buffer to host send large buffer
        size_t sizeArray = dimBuffer * nIFace * sizeof(double);
        /*
  //!Add for test
  printf("nIFace = %d, dimBuffer = %d\n", nIFace, dimBuffer);
  cout<<"name ="<<name<<endl;
  //!Add end
  */
        HANDLE_API_ERR(cudaMemcpyAsync(HostSendLargeBuffer, deviceSendBuffer, sizeArray, cudaMemcpyDeviceToHost,
                                       dataTransferDToH));
    }

    void GPUSendLargeBufferToHostSendLargeBufferInterpoint(const int solverID, InterpointInformation *iinfo)
    {
        //! just for 1 process case
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        if (solverID == 0)
        { //! NS euqation
            size_t bufferSize = lengthInterpointBufferNS * sizeof(double);
            HANDLE_API_ERR(
                cudaMemcpy(HostSendInterpointBufferNS, GPUSendInterpointBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        }
        else if (solverID == 1)
        { //! Turb equation
            size_t bufferSize = lengthInterpointBufferTurb * sizeof(double);
            HANDLE_API_ERR(cudaMemcpy(HostSendInterpointBufferTurb, GPUSendInterpointBufferTurb, bufferSize,
                                      cudaMemcpyDeviceToHost));
        }
    }

    void HostSendLargeBufferToOrgSendBufferInterpoint(const int solverID, InterpointInformation *iinfo)
    {
        //! just for 1 process case
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        Data_ParamFieldSuite *sendData = iinfo->GetSendDataStorage(0);
        int                   nIPoint  = iinfo->GetNumberOfInterpoints();
        //! printf("solverID = %d\n", solverID);
        if (solverID == 0)
        { //! NS euqation
            for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
            {
                const string &name    = (*nameNSInterpointBuffer)[nameID];
                int           offset  = offsetInterpointBufferNS[nameID];
                int           dimData = dimInterpointBufferNS[nameID];
                RFloat      **fg      = reinterpret_cast<RFloat **>(sendData->GetDataPtr(name));
#ifdef CUDAUNITTEST
                for (int pointID = 0; pointID < nIPoint; pointID++)
                {
                    for (int m = 0; m < dimData; ++m)
                    {
                        //! add for test
                        double error =
                            fabs(fg[m][pointID] - HostSendInterpointBufferNS[offset + m * nIPoint + pointID]);
                        if (error > 1.0e-16)
                        {
                            printf("name = %s, fg = %.30e, HostSendBuffer = %.30e, err = %.30e\n", name.c_str(),
                                   fg[m][pointID], HostSendInterpointBufferNS[offset + m * nIPoint + pointID], error);
                            exit(0);
                        }
                        //! add end
                        //! fg[m][faceID] = HostSendBufferNS[offset + m*nIFace + faceID];
                    }
                }
#endif
                memcpy(fg[0], HostSendInterpointBufferNS + offset, dimData * nIPoint * sizeof(double));
            }
        }
        else if (solverID == 1)
        { //! Turbulence equation
            for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
            {
                const string &name    = (*nameTurbInterpointBuffer)[nameID];
                int           offset  = offsetInterpointBufferTurb[nameID];
                int           dimData = dimInterpointBufferTurb[nameID];
                RFloat      **fg      = reinterpret_cast<RFloat **>(sendData->GetDataPtr(name));
#ifdef CUDAUNITTEST
                for (int pointID = 0; pointID < nIPoint; pointID++)
                {
                    for (int m = 0; m < dimData; ++m)
                    {
                        //! add for test
                        double error =
                            fabs(fg[m][pointID] - HostSendInterpointBufferTurb[offset + m * nIPoint + pointID]);
                        if (error > 1.0e-16)
                        {
                            printf("name = %s, fg = %.30e, HostSendBuffer = %.30e, err = %.30e\n", name.c_str(),
                                   fg[m][pointID], HostSendInterpointBufferTurb[offset + m * nIPoint + pointID], error);
                            exit(0);
                        }
                        //! add end
                        //! fg[m][faceID] = HostSendBufferNS[offset + m*nIFace + faceID];
                    }
                }
#endif
                memcpy(fg[0], HostSendInterpointBufferTurb + offset, dimData * nIPoint * sizeof(double));
            }
        }
    }

    void HostRecvLargeBufferToGPURecvSmallBufferInterpoint(const int solverID, InterpointInformation *iinfo)
    {
        //! just for 1 process case
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        double *HostRecvLargeBuffer;
        RFloat *deviceIFVar;
        size_t  arraySize;
        int     nIPoint = iinfo->GetNumberOfInterpoints();

        if (solverID == 0)
        { //! NS euqation
            for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
            {
                const string &name    = (*nameNSInterpointBuffer)[nameID];
                int           offset  = offsetInterpointBufferNS[nameID];
                int           dimData = dimInterpointBufferNS[nameID];
                HostRecvLargeBuffer   = HostRecvInterpointBufferNS + offset;
                arraySize             = dimData * nIPoint * sizeof(double);
                if ("qnode" == name)
                {
                    deviceIFVar = d_fg_recv_qNode;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToDInterpoint));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSqNode_done, dataTransferHToDInterpoint));
                }
                else if ("tnode" == name)
                {
                    deviceIFVar = d_fg_recv_tNode;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToDInterpoint));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNStNode_done, dataTransferHToDInterpoint));
                }
            }
        }
        else if (solverID == 1)
        { //! Turb equation
            for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
            {
                const string &name    = (*nameTurbInterpointBuffer)[nameID];
                int           offset  = offsetInterpointBufferTurb[nameID];
                int           dimData = dimInterpointBufferTurb[nameID];
                HostRecvLargeBuffer   = HostRecvInterpointBufferTurb + offset;
                arraySize             = dimData * nIPoint * sizeof(double);

                if ("qTurbNode" == name)
                {
                    deviceIFVar = d_fg_recv_qTurbNode;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToDInterpoint));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDTurbqNode_done, dataTransferHToDInterpoint));
                }
            }
        }
    }

    void HostRecvLargeBufferToGPURecvLargeBufferInterpoint(const int solverID, InterpointInformation *iinfo)
    {
        //! just for 1 process case
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        if (solverID == 0)
        { //! NS euqation
            size_t bufferSize = lengthInterpointBufferNS * sizeof(double);
            HANDLE_API_ERR(
                cudaMemcpy(GPURecvInterpointBufferNS, HostRecvInterpointBufferNS, bufferSize, cudaMemcpyHostToDevice));
        }
        else if (solverID == 1)
        { //! Turb equation
            size_t bufferSize = lengthInterpointBufferTurb * sizeof(double);
            HANDLE_API_ERR(cudaMemcpy(GPURecvInterpointBufferTurb, HostRecvInterpointBufferTurb, bufferSize,
                                      cudaMemcpyHostToDevice));
        }
    }

    void OrgRecvBufferToHostRecvLargeBufferInterpoint(const int solverID, InterpointInformation *iinfo)
    {
        //! just for 1 process case
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        Data_ParamFieldSuite *recvData = iinfo->GetReceiveDataStorage(0);
        int                   nIPoint  = iinfo->GetNumberOfInterpoints();
        if (solverID == 0)
        { //! NS euqation
            for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
            {
                const string &name    = (*nameNSInterpointBuffer)[nameID];
                int           offset  = offsetInterpointBufferNS[nameID];
                int           dimData = dimInterpointBufferNS[nameID];

                RFloat **fg = reinterpret_cast<RFloat **>(recvData->GetDataPtr(name));
                /*
      for (int faceID = 0; faceID < nIFace; faceID++){
              for ( int m = 0; m < dimData; ++ m ){
                     //!add for test
                      //!printf("name = %s, fg[%d][%d] = %.30e\n", name.c_str(),
      m, faceID, fg[m][faceID]);
                      //!add end
                      //!HostRecvBufferNS[offset + m*nIFace + faceID] =
      fg[m][faceID];
              }
      }
      */
                memcpy(HostRecvInterpointBufferNS + offset, fg[0], dimData * nIPoint * sizeof(double));
            }
        }
        else if (solverID == 1)
        {
            for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
            {
                const string &name    = (*nameTurbInterpointBuffer)[nameID];
                int           offset  = offsetInterpointBufferTurb[nameID];
                int           dimData = dimInterpointBufferTurb[nameID];

                RFloat **fg = reinterpret_cast<RFloat **>(recvData->GetDataPtr(name));
                /*
      for (int faceID = 0; faceID < nIFace; faceID++){
              for ( int m = 0; m < dimData; ++ m ){
                      //!HostRecvBufferTurb[offset + m*nIFace + faceID] =
      fg[m][faceID];
              }
      }
      */

                memcpy(HostRecvInterpointBufferTurb + offset, fg[0], dimData * nIPoint * sizeof(double));
            }
        }
    }
    void GPUSendLargeBufferToHostSendLargeBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        //! just for 1 process case
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        if (solverID == 0)
        { //! NS euqation
            size_t bufferSize = lengthBufferNS * sizeof(double);
            HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        }
        else if (solverID == 1)
        { //! Turb equation
            size_t bufferSize = lengthBufferTurb * sizeof(double);
            HANDLE_API_ERR(cudaMemcpy(HostSendBufferTurb, GPUSendBufferTurb, bufferSize, cudaMemcpyDeviceToHost));
        }
    }

    void HostSendLargeBufferToOrgSendBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        //! just for 1 process case
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        Data_ParamFieldSuite *sendData = iinfo->GetSendDataStorage(0);
        int                   nIFace   = iinfo->GetNIFace();
        if (solverID == 0)
        { //! NS euqation
            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                const string &name    = (*nameNSBuffer)[nameID];
                int           offset  = offsetBufferNS[nameID];
                int           dimData = dimBufferNS[nameID];
                //! add for test
                //! printf("solverID = %d, name = %s\n", solverID, name.c_str());
                //! add end
                //! Data_ParamFieldSuite * sendData = iinfo->GetSendDataStorage(0);
                RFloat **fg = reinterpret_cast<RFloat **>(sendData->GetDataPtr(name));
#ifdef CUDAUNITTEST
                for (int faceID = 0; faceID < nIFace; faceID++)
                {
                    for (int m = 0; m < dimData; ++m)
                    {
                        //! add for test
                        double error = fabs(fg[m][faceID] - HostSendBufferNS[offset + m * nIFace + faceID]);
                        if (error > 1.0e-16)
                        {
                            printf("name = %s, fg = %.30e, HostSendBuffer = %.30e, err = %.30e\n", name.c_str(),
                                   fg[m][faceID], HostSendBufferNS[offset + m * nIFace + faceID], error);
                            exit(0);
                        }
                        //! add end
                        //! fg[m][faceID] = HostSendBufferNS[offset + m*nIFace + faceID];
                    }
                }
#endif
                memcpy(fg[0], HostSendBufferNS + offset, dimData * nIFace * sizeof(double));
                //! add for test
                //! if (fg) printf("solverID = %d, name = %s\n", solverID, name.c_str());
                //! add end
            }
        }
        else if (solverID == 1)
        { //! Turbulence equation
            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                const string &name    = (*nameTurbBuffer)[nameID];
                int           offset  = offsetBufferTurb[nameID];
                int           dimData = dimBufferTurb[nameID];
                //! add for test
                //! printf("solverID = %d, name = %s\n", solverID, name.c_str());
                //! add end
                //! Data_ParamFieldSuite * sendData = iinfo->GetSendDataStorage(0);
                RFloat **fg = reinterpret_cast<RFloat **>(sendData->GetDataPtr(name));
#ifdef CUDAUNITTEST
                for (int faceID = 0; faceID < nIFace; faceID++)
                {
                    for (int m = 0; m < dimData; ++m)
                    {
                        //! add for test
                        double error = fabs(fg[m][faceID] - HostSendBufferTurb[offset + m * nIFace + faceID]);
                        if (error > 1.0e-16)
                        {
                            printf("name = %s, fg = %.30e, HostSendBuffer = %.30e, err = %.30e\n", name.c_str(),
                                   fg[m][faceID], HostSendBufferTurb[offset + m * nIFace + faceID], error);
                            exit(0);
                        }
                        //! add end
                        //! fg[m][faceID] = HostSendBufferTurb[offset + m*nIFace + faceID];
                    }
                }
#endif
                memcpy(fg[0], HostSendBufferTurb + offset, dimData * nIFace * sizeof(double));

                //! add for test
                //! if (fg) printf("solverID = %d, name = %s\n", solverID, name.c_str());
                //! add end
            }
        }
        else
        {
            printf("Error: solverID = %d is not considere\n", solverID);
            exit(0);
        }
    }

    void HostRecvLargeBufferToGPURecvSmallBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        size_t  arraySize;
        double *HostRecvLargeBuffer;
        double *deviceIFVar;
        int     nIFace = iinfo->GetNIFace();
        if (solverID == 0)
        { //! NS euqation
            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                const string &name    = (*nameNSBuffer)[nameID];
                int           offset  = offsetBufferNS[nameID];
                int           dimData = dimBufferNS[nameID];
                HostRecvLargeBuffer   = HostRecvBufferNS + offset;
                arraySize             = dimData * nIFace * sizeof(double);
                if ("q" == name)
                {
                    deviceIFVar = d_fg_recv_q;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    //! record the finish of data transfered from host large recv buffer to
                    //! device small recv buffer
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSq_done, dataTransferHToD));
                    //! Next D2D (data transfer from device small recv buffer to device
                    //! variable) is also put in stream dataTransferH2D. So, at the moment,
                    //! cudaEventRecord is not required. However, it is not good, all of H2D
                    //! async operations will be blocked.
                }
                else if ("t" == name)
                {
                    deviceIFVar = d_fg_recv_t;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    //! record the finish of data transfered from host large recv buffer to
                    //! device small recv buffer
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSt_done, dataTransferHToD));
                    //! Next D2D (data transfer from device small recv buffer to device
                    //! variable) is also put in stream dataTransferH2D. So, at the moment,
                    //! cudaEventRecord is not required. However, it is not efficient.
                }
                else if ("dqdx" == name)
                {
                    deviceIFVar = d_fg_recv_dqdx;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    //! record the finish of data transfered from host large recv buffer to
                    //! device small recv buffer
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSdqdx_done, dataTransferHToD));
                }
                else if ("dqdy" == name)
                {
                    deviceIFVar = d_fg_recv_dqdy;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSdqdy_done, dataTransferHToD));
                }
                else if ("dqdz" == name)
                {
                    deviceIFVar = d_fg_recv_dqdz;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSdqdz_done, dataTransferHToD));
                }
                else if ("limit" == name)
                {
                    deviceIFVar = d_fg_recv_limit;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSlimit_done, dataTransferHToD));
                }
                else if ("dtdx" == name)
                {
                    deviceIFVar = d_fg_recv_dtdx;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSdtdx_done, dataTransferHToD));
                }
                else if ("dtdy" == name)
                {
                    deviceIFVar = d_fg_recv_dtdy;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSdtdy_done, dataTransferHToD));
                }
                else if ("dtdz" == name)
                {
                    deviceIFVar = d_fg_recv_dtdz;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDNSdtdz_done, dataTransferHToD));
                }
                else
                {
                    //! no operations, because of visl, vist
                }
            }
        }
        else if (solverID == 1)
        { //! Turb equation
            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                const string &name    = (*nameTurbBuffer)[nameID];
                int           offset  = offsetBufferTurb[nameID];
                int           dimData = dimBufferTurb[nameID];
                HostRecvLargeBuffer   = HostRecvBufferTurb + offset;
                arraySize             = dimData * nIFace * sizeof(double);

                if ("turb::q" == name)
                {
                    deviceIFVar = d_fg_recv_qTurb;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDTurbq_done, dataTransferHToD));
                    //! Next D2D (data transfer from device small recv buffer to device
                    //! variable) is also put in stream dataTransferH2D. So, at the moment,
                    //! cudaEventRecord is not required. However, it is not efficient.
                }
                else if ("turb::dqdx" == name)
                {
                    deviceIFVar = d_fg_recv_dqTurbdx;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDTurbdqdx_done, dataTransferHToD));
                }
                else if ("turb::dqdy" == name)
                {
                    deviceIFVar = d_fg_recv_dqTurbdy;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDTurbdqdy_done, dataTransferHToD));
                }
                else if ("turb::dqdz" == name)
                {
                    deviceIFVar = d_fg_recv_dqTurbdz;
                    HANDLE_API_ERR(cudaMemcpyAsync(deviceIFVar, HostRecvLargeBuffer, arraySize, cudaMemcpyHostToDevice,
                                                   dataTransferHToD));
                    HANDLE_API_ERR(cudaEventRecord(downloadHToDTurbdqdz_done, dataTransferHToD));
                }
                else
                {
                    //! no operations, because of visl, vist
                }
            }
        }
    }

    void HostRecvLargeBufferToGPURecvLargeBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        if (solverID == 0)
        { //! NS euqation
            size_t bufferSize = lengthBufferNS * sizeof(double);
            HANDLE_API_ERR(cudaMemcpy(GPURecvBufferNS, HostRecvBufferNS, bufferSize, cudaMemcpyHostToDevice));
        }
        else if (solverID == 1)
        { //! Turb equation
            size_t bufferSize = lengthBufferTurb * sizeof(double);
            HANDLE_API_ERR(cudaMemcpy(GPURecvBufferTurb, HostRecvBufferTurb, bufferSize, cudaMemcpyHostToDevice));
        }
    }
    //! void HostRecvLargeBufferToOrgBuffer(const int solverID, InterfaceInfo *
    //! iinfo){
    void OrgRecvBufferToHostRecvLargeBuffer(const int solverID, InterfaceInfo *iinfo)
    {
        //! just for 1 process case, without it, bugs will be induced.
        if (!iinfo) return;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        Data_ParamFieldSuite *recvData = iinfo->GetRecvDataStorage(0);
        int                   nIFace   = iinfo->GetNIFace();
        if (solverID == 0)
        { //! NS euqation
            for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                const string &name    = (*nameNSBuffer)[nameID];
                int           offset  = offsetBufferNS[nameID];
                int           dimData = dimBufferNS[nameID];

                RFloat **fg = reinterpret_cast<RFloat **>(recvData->GetDataPtr(name));
                /*
      for (int faceID = 0; faceID < nIFace; faceID++){
              for ( int m = 0; m < dimData; ++ m ){
                     //!add for test
                      //!printf("name = %s, fg[%d][%d] = %.30e\n", name.c_str(),
      m, faceID, fg[m][faceID]);
                      //!add end
                      //!HostRecvBufferNS[offset + m*nIFace + faceID] =
      fg[m][faceID];
              }
      }
      */
                memcpy(HostRecvBufferNS + offset, fg[0], dimData * nIFace * sizeof(double));
            }
        }
        else if (solverID == 1)
        {
            for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                const string &name    = (*nameTurbBuffer)[nameID];
                int           offset  = offsetBufferTurb[nameID];
                int           dimData = dimBufferTurb[nameID];

                RFloat **fg = reinterpret_cast<RFloat **>(recvData->GetDataPtr(name));
                /*
      for (int faceID = 0; faceID < nIFace; faceID++){
              for ( int m = 0; m < dimData; ++ m ){
                      //!HostRecvBufferTurb[offset + m*nIFace + faceID] =
      fg[m][faceID];
              }
      }
      */

                memcpy(HostRecvBufferTurb + offset, fg[0], dimData * nIFace * sizeof(double));
            }
        }
    }

    void CallDeviceUploadInterpointValueHostLargeBuffer(const int nTotalNode, const int neqn, const int nIPoint,
                                                        const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        int     solverID;
        int     nameNSID, nameTurbID;
        int     offsetBuffer;
        RFloat *deviceSendBuffer;
        RFloat *deviceFieldVar;
        RFloat *HostSendLargeBuffer;
        int     dimBuffer;

        if ("qnode" == name)
        {
            solverID            = 0;
            nameNSID            = GetIndexInterpointBufferNS(name);
            offsetBuffer        = offsetInterpointBufferNS[nameNSID];
            dimBuffer           = dimInterpointBufferNS[nameNSID];
            deviceFieldVar      = d_qNode;
            deviceSendBuffer    = d_fg_send_qNode;
            HostSendLargeBuffer = HostSendInterpointBufferNS + offsetBuffer;
        }
        else if ("tnode" == name)
        {
            solverID            = 0;
            nameNSID            = GetIndexInterpointBufferNS(name);
            offsetBuffer        = offsetInterpointBufferNS[nameNSID];
            dimBuffer           = dimInterpointBufferNS[nameNSID];
            deviceFieldVar      = d_tNode;
            deviceSendBuffer    = d_fg_send_tNode;
            HostSendLargeBuffer = HostSendInterpointBufferNS + offsetBuffer;
        }
        else if ("qTurbNode" == name)
        {
            solverID            = 1;
            nameTurbID          = GetIndexInterpointBufferTurb(name);
            offsetBuffer        = offsetInterpointBufferTurb[nameTurbID];
            dimBuffer           = dimInterpointBufferTurb[nameTurbID];
            deviceFieldVar      = d_qTurbNode;
            deviceSendBuffer    = d_fg_send_qTurbNode;
            HostSendLargeBuffer = HostSendInterpointBufferTurb + offsetBuffer;
        }
        else
        {
            cout << "Error:" << name << "is not included in CallGPUUploadInterpointValue" << endl;
            exit(1);
        }
        //! D2D:From device varaible to device send buffer
        CallGPUFieldVarToIFVarInterpoint(deviceSendBuffer, deviceFieldVar, nTotalNode, neqn, nIPoint);
        //! record and wait for D2D data transfer
        if (solverID == 0)
        { //! NS equation
            //! record finish of DToD data transfer
            HANDLE_API_ERR(cudaEventRecord(downloadDToDNS_done, 0));
            //! wait for event downlaodDTODNS_done
            HANDLE_API_ERR(cudaStreamWaitEvent(dataTransferDToHInterpoint, downloadDToDNS_done, 0));
        }
        else if (solverID == 1)
        { //! Turbulence equation
            HANDLE_API_ERR(cudaEventRecord(downloadDToDTurb_done, 0));
            HANDLE_API_ERR(cudaStreamWaitEvent(dataTransferDToHInterpoint, downloadDToDTurb_done, 0));
        }
        //! D2H:From device send buffer to host send large buffer
        size_t sizeArray = dimBuffer * nIPoint * sizeof(double);
        HANDLE_API_ERR(cudaMemcpyAsync(HostSendLargeBuffer, deviceSendBuffer, sizeArray, cudaMemcpyDeviceToHost,
                                       dataTransferDToHInterpoint));
    }

    void CallGPUUploadInterpointValueLargeBuffer(const int nTotalNode, const int neqn, const int nIPoint,
                                                 const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        int     nameNSID, nameTurbID;
        int     offsetBuffer;
        RFloat *deviceSendBuffer;
        RFloat *deviceFieldVar;
        if ("qnode" == name)
        {
            nameNSID         = GetIndexInterpointBufferNS(name);
            offsetBuffer     = offsetInterpointBufferNS[nameNSID];
            deviceFieldVar   = d_qNode;
            deviceSendBuffer = GPUSendInterpointBufferNS;
        }
        else if ("tnode" == name)
        {
            nameNSID         = GetIndexInterpointBufferNS(name);
            offsetBuffer     = offsetInterpointBufferNS[nameNSID];
            deviceFieldVar   = d_tNode;
            deviceSendBuffer = GPUSendInterpointBufferNS;
        }
        else if ("qTurbNode" == name)
        {
            nameTurbID       = GetIndexInterpointBufferTurb(name);
            offsetBuffer     = offsetInterpointBufferTurb[nameTurbID];
            deviceFieldVar   = d_qTurbNode;
            deviceSendBuffer = GPUSendInterpointBufferTurb;
        }
        else
        {
            cout << "Error:" << name << "is not included in CallGPUUploadInterpointValue" << endl;
            exit(1);
        }
        CallGPUFieldToInterpointBuffer(deviceFieldVar, deviceSendBuffer, nTotalNode, neqn, nIPoint, offsetBuffer);
    }

    void CallGPUFieldToInterpointBuffer(const RFloat *deviceFieldVar, double *deviceSendBuffer, const int nTotalNode,
                                        const int neqn, const int nIPoint, const int offsetBuffer)
    {
        using namespace GPUMemory;
        using namespace GPUGeomVariables;
        int blockSize, gridSize, loopLength;
        loopLength         = nIPoint;
        int blockSizeLimit = 0;
        blockSize          = 1;
        gridSize           = 1;
        KernelLaunchPara((void *)GPUFieldToInterpointBuffer, loopLength, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUFieldToInterpointBuffer, 0, gridSize, blockSize);
#endif
        GPUFieldToInterpointBuffer<<<gridSize, blockSize>>>(deviceFieldVar, deviceSendBuffer, d_interPoint2GlobalPoint,
                                                            nTotalNode, neqn, nIPoint, offsetBuffer);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUFieldToInterpointBuffer(const RFloat *deviceFieldVar, double *deviceSendBuffer,
                                               const int *interPoint2GlobalPoint, const int nTotalNode, const int neqn,
                                               const int nIPoint, const int offset)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int iPoint;
        int sourcePoint;
        for (iPoint = threadID; iPoint < nIPoint; iPoint += blockDim.x * gridDim.x)
        {
            sourcePoint = interPoint2GlobalPoint[iPoint];
            for (int m = 0; m < neqn; ++m)
            {
                deviceSendBuffer[offset + m * nIPoint + iPoint] = deviceFieldVar[m * nTotalNode + sourcePoint];
            }
        }
    }

    void CallGPUDownloadInterpointValueSmallBuffer(const int nTotalNode, const int neqn, const int nIPoint,
                                                   const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        RFloat *deviceIFVar;
        RFloat *deviceFieldVar;
        size_t  sizeIF = neqn * nIPoint * sizeof(RFloat);

        if ("qnode" == name)
        {
            deviceIFVar    = d_fg_recv_qNode;
            deviceFieldVar = d_qInterPoint;
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSqNode_done, 0));
            CallGPUIFVarToFieldVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else if ("tnode" == name)
        {
            deviceIFVar    = d_fg_recv_tNode;
            deviceFieldVar = d_tInterPoint;
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNStNode_done, 0));
            CallGPUIFVarToFieldVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else if ("qTurbNode" == name)
        {
            deviceIFVar    = d_fg_recv_qTurbNode;
            deviceFieldVar = d_qTurbInterPoint;
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDTurbqNode_done, 0));
            CallGPUIFVarToFieldVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else
        {
            cout << "Error:" << name << "is not included in CallGPUDownloadInterpointValue" << endl;
            exit(1);
        }
    }

    void CallGPUDownloadInterpointValueLargeBuffer(const int nTotalNode, const int neqn, const int nIPoint,
                                                   const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        double *deviceBuffer;
        RFloat *deviceFieldVar;
        int     nameNSID, nameTurbID;
        int     offsetBuffer;

        if ("qnode" == name)
        {
            nameNSID       = GetIndexInterpointBufferNS(name);
            offsetBuffer   = offsetInterpointBufferNS[nameNSID];
            deviceBuffer   = GPURecvInterpointBufferNS;
            deviceFieldVar = d_qInterPoint;
        }
        else if ("tnode" == name)
        {
            nameNSID       = GetIndexInterpointBufferNS(name);
            offsetBuffer   = offsetInterpointBufferNS[nameNSID];
            deviceBuffer   = GPURecvInterpointBufferNS;
            deviceFieldVar = d_tInterPoint;
        }
        else if ("qTurbNode" == name)
        {
            nameTurbID     = GetIndexInterpointBufferTurb(name);
            offsetBuffer   = offsetInterpointBufferTurb[nameTurbID];
            deviceBuffer   = GPURecvInterpointBufferTurb;
            deviceFieldVar = d_qTurbInterPoint;
        }
        else
        {
            cout << "Error:" << name << "is not included in CallGPUDownloadInterpointValue" << endl;
            exit(1);
        }

        CallGPUInterpointBufferToField(deviceFieldVar, deviceBuffer, nTotalNode, neqn, nIPoint, offsetBuffer);
    }

    void CallGPUInterpointBufferToField(RFloat *deviceFieldVar, const double *deviceRecvBuffer, const int nTotalNode,
                                        const int neqn, const int nIPoint, const int offsetBuffer)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int blockSize, gridSize, loopLength;
        loopLength         = nIPoint;
        int blockSizeLimit = 0;
        blockSize          = 1;
        gridSize           = 1;
        KernelLaunchPara((void *)GPUInterpointBufferToField, loopLength, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUInterpointBufferToField, 0, gridSize, blockSize);
#endif
        GPUInterpointBufferToField<<<gridSize, blockSize>>>(deviceFieldVar, deviceRecvBuffer, nTotalNode, neqn, nIPoint,
                                                            offsetBuffer);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUInterpointBufferToField(RFloat *deviceFieldVar, const double *deviceBuffer, const int nTotalNode,
                                               const int neqn, const int nIPoint, const int offset)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int iPoint;
        int sourcePoint;
        for (iPoint = threadID; iPoint < nIPoint; iPoint += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < neqn; ++m)
            {
                deviceFieldVar[m * nIPoint + iPoint] = deviceBuffer[offset + m * nIPoint + iPoint];
            }
        }
    }

    void CallGPUUploadInterpointValue(const int nTotalNode, const int neqn, const int nIPoint, const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        RFloat *deviceIFVar;
        RFloat *deviceFieldVar;
        if ("qnode" == name)
        {
            deviceIFVar    = d_fg_send_qNode;
            deviceFieldVar = d_qNode;
            CallGPUFieldVarToIFVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else if ("tnode" == name)
        {
            deviceIFVar    = d_fg_send_tNode;
            deviceFieldVar = d_tNode;
            CallGPUFieldVarToIFVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else if ("qTurbNode" == name)
        {
            deviceIFVar    = d_fg_send_qTurbNode;
            deviceFieldVar = d_qTurbNode;
            CallGPUFieldVarToIFVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else
        {
            cout << "Error:" << name << "is not included in CallGPUUploadInterpointValue" << endl;
            exit(1);
        }
    }
    void CallGPUDownloadInterpointValue(RFloat **fg, const int nTotalNode, const int neqn, const int nIPoint,
                                        const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        RFloat       *deviceIFVar;
        RFloat       *deviceFieldVar;
        size_t        sizeIF = neqn * nIPoint * sizeof(RFloat);
        const RFloat *fg0    = fg[0];

        if ("qnode" == name)
        {
            deviceIFVar = d_fg_recv_qNode;
            //! deviceFieldVar = d_qNode;
            deviceFieldVar = d_qInterPoint;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else if ("tnode" == name)
        {
            deviceIFVar = d_fg_recv_tNode;
            //! deviceFieldVar = d_tNode;
            deviceFieldVar = d_tInterPoint;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else if ("qTurbNode" == name)
        {
            deviceIFVar = d_fg_recv_qTurbNode;
            //! deviceFieldVar = d_qTurbNode;
            deviceFieldVar = d_qTurbInterPoint;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVarInterpoint(deviceIFVar, deviceFieldVar, nTotalNode, neqn, nIPoint);
        }
        else
        {
            cout << "Error:" << name << "is not included in CallGPUDownloadInterpointValue" << endl;
            exit(1);
        }
    }
    void CallGPUDownloadInterfaceValue(RFloat **fg, const int nTotal, const int neqn, const int nIFace,
                                       const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        RFloat       *deviceIFVar;
        RFloat       *deviceFieldVar;
        size_t        sizeIF = neqn * nIFace * sizeof(RFloat);
        const RFloat *fg0    = fg[0];
        if ("q" == name)
        {
            deviceIFVar    = d_fg_recv_q;
            deviceFieldVar = d_q_ns;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("t" == name)
        {
            deviceIFVar    = d_fg_recv_t;
            deviceFieldVar = d_t_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::q" == name)
        {
            deviceIFVar    = d_fg_recv_qTurb;
            deviceFieldVar = d_q_turb_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dqdx" == name)
        {
            deviceIFVar    = d_fg_recv_dqdx;
            deviceFieldVar = d_dqdx_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dqdy" == name)
        {
            deviceIFVar    = d_fg_recv_dqdy;
            deviceFieldVar = d_dqdy_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dqdz" == name)
        {
            deviceIFVar    = d_fg_recv_dqdz;
            deviceFieldVar = d_dqdz_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("limit" == name)
        {
            deviceIFVar    = d_fg_recv_limit;
            deviceFieldVar = d_limit;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdx" == name)
        {
            deviceIFVar    = d_fg_recv_dtdx;
            deviceFieldVar = d_dtdx_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdy" == name)
        {
            deviceIFVar    = d_fg_recv_dtdy;
            deviceFieldVar = d_dtdy_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdz" == name)
        {
            deviceIFVar    = d_fg_recv_dtdz;
            deviceFieldVar = d_dtdz_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdx" == name)
        {
            deviceIFVar    = d_fg_recv_dqTurbdx;
            deviceFieldVar = d_dq_turbdx_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdy" == name)
        {
            deviceIFVar    = d_fg_recv_dqTurbdy;
            deviceFieldVar = d_dq_turbdy_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdz" == name)
        {
            deviceIFVar    = d_fg_recv_dqTurbdz;
            deviceFieldVar = d_dq_turbdz_proxy;
            HANDLE_API_ERR(cudaMemcpy(deviceIFVar, fg0, sizeIF, cudaMemcpyHostToDevice));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else
        {
            cout << "Error:" << name << " is not included in CallGPUDownloadInterfaceValue" << endl;
            exit(1);
        }
    }

    void CallGPUDownloadInterfaceValueSmallBuffer(const int nTotal, const int neqn, const int nIFace,
                                                  const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        RFloat *deviceFieldVar;
        double *deviceIFVar;

        if ("q" == name)
        {
            deviceIFVar    = d_fg_recv_q;
            deviceFieldVar = d_q_ns;
            //! stream 0 should wait for the finish of downloadHToDNSq_done
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSq_done, 0));
            //! date transfered from device small buffer to device variable
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
            //! On stream dataTransferHToD, data from small device recv buffer to device
            //! variable CallGPUIFVarToFieldVarHToDAsync(deviceIFVar, deviceFieldVar,
            //! nTotal, neqn, nIFace);
            //! HANDLE_API_ERR(cudaEventRecord(downloadHToDNSq_done, dataTransferHToD));
            //! Stream 0 should wait for finish of HToD. Due to UpdateInterfaceData, it
            //! is put in the t transfer. It can be referenced by
            //! HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSq_done, 0));
        }
        else if ("t" == name)
        {
            deviceIFVar    = d_fg_recv_t;
            deviceFieldVar = d_t_proxy;
            //! stream 0 should wait for the finish of downloadHToDNSt_done
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSt_done, 0));
            //! date transfered from device small buffer to device variable
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
            //! On stream dataTransferHToD, data from small device recv buffer to device
            //! variable CallGPUIFVarToFieldVarHToDAsync(deviceIFVar, deviceFieldVar,
            //! nTotal, neqn, nIFace);
            //! HANDLE_API_ERR(cudaEventRecord(downloadHToDNSt_done, dataTransferHToD));
            //! Stream 0 should wait for finish of HToD. The position can be changed. At
            //! the current stage, the HToD is not hiden.
            //! HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSt_done, 0));
        }
        else if ("turb::q" == name)
        {
            deviceIFVar    = d_fg_recv_qTurb;
            deviceFieldVar = d_q_turb_proxy;
            //! stream 0 should wait for the finish of downloadHToDTurbq_done
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDTurbq_done, 0));
            //! date transfered from device small buffer to device variable
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
            //! On stream dataTransferHToD, data from small device recv buffer to device
            //! variable CallGPUIFVarToFieldVarHToDAsync(deviceIFVar, deviceFieldVar,
            //! nTotal, neqn, nIFace);
            //! HANDLE_API_ERR(cudaEventRecord(downloadHToDTurbq_done,
            //! dataTransferHToD)); Stream 0 should wait for finish of HToD. The position
            //! can be changed. At the current stage, the HToD is not hiden.
            //! HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDTurbq_done, 0));
        }
        else if ("dqdx" == name)
        {
            deviceIFVar    = d_fg_recv_dqdx;
            deviceFieldVar = d_dqdx_proxy;
            //! downloadHToDNSdqdx_done record the finish of data from Host large recv
            //! buffer to device small recv buffer
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSdqdx_done, 0));
            //! On default stream 0, data transfer from small recv device buffer to
            //! device variable
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dqdy" == name)
        {
            deviceIFVar    = d_fg_recv_dqdy;
            deviceFieldVar = d_dqdy_proxy;
            //! downloadHToDNSdqdy_done record the finish of data from Host large recv
            //! buffer to device small recv buffer
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSdqdy_done, 0));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dqdz" == name)
        {
            deviceIFVar    = d_fg_recv_dqdz;
            deviceFieldVar = d_dqdz_proxy;
            //! downloadHToDNSdqdz_done record the finish of data from Host large recv
            //! buffer to device small recv buffer
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSdqdz_done, 0));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("limit" == name)
        {
            deviceIFVar    = d_fg_recv_limit;
            deviceFieldVar = d_limit;
            //! downloadHToDNSlimit_done record the finish of data from Host large recv
            //! buffer to device small recv buffer
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSlimit_done, 0));
            //! On default stream 0, data transfer from small recv device buffer to
            //! device variable
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdx" == name)
        {
            deviceIFVar    = d_fg_recv_dtdx;
            deviceFieldVar = d_dtdx_proxy;
            //! downloadHToDNSdtdx_done record the finish of data from Host large recv
            //! buffer to device small recv buffer
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSdtdx_done, 0));
            //! On default stream 0, data transfer from small recv device buffer to
            //! device variable
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdy" == name)
        {
            deviceIFVar    = d_fg_recv_dtdy;
            deviceFieldVar = d_dtdy_proxy;
            //! downloadHToDNSdtdy_done record the finish of data from Host large recv
            //! buffer to device small recv buffer
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSdtdy_done, 0));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdz" == name)
        {
            deviceIFVar    = d_fg_recv_dtdz;
            deviceFieldVar = d_dtdz_proxy;
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDNSdtdz_done, 0));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdx" == name)
        {
            deviceIFVar    = d_fg_recv_dqTurbdx;
            deviceFieldVar = d_dq_turbdx_proxy;
            //! downloadHToDTurbdqdx_done record the finish of data from Host large recv
            //! buffer to device small recv buffer
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDTurbdqdx_done, 0));
            //! On default stream 0, data transfer from small recv device buffer to
            //! device variable
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdy" == name)
        {
            deviceIFVar    = d_fg_recv_dqTurbdy;
            deviceFieldVar = d_dq_turbdy_proxy;
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDTurbdqdy_done, 0));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdz" == name)
        {
            deviceIFVar    = d_fg_recv_dqTurbdz;
            deviceFieldVar = d_dq_turbdz_proxy;
            HANDLE_API_ERR(cudaStreamWaitEvent(0, downloadHToDTurbdqdz_done, 0));
            CallGPUIFVarToFieldVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else
        {
            cout << "Error:" << name << " is not included in CallGPUDownloadInterfaceValue" << endl;
            exit(1);
        }
    }

    void CallGPUDownloadInterfaceValueLargeBuffer(RFloat **fg, const int nTotal, const int neqn, const int nIFace,
                                                  const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        int     nameNSID, nameTurbID;
        int     offsetBuffer;
        RFloat *deviceFieldVar;
        double *deviceRecvBuffer;

        if ("q" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_q_ns;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("t" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_t_proxy;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("turb::q" == name)
        {
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            deviceFieldVar   = d_q_turb_proxy;
            deviceRecvBuffer = GPURecvBufferTurb;
        }
        else if ("dqdx" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dqdx_proxy;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("dqdy" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dqdy_proxy;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("dqdz" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dqdz_proxy;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("limit" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_limit;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("dtdx" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dtdx_proxy;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("dtdy" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dtdy_proxy;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("dtdz" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dtdz_proxy;
            deviceRecvBuffer = GPURecvBufferNS;
        }
        else if ("turb::dqdx" == name)
        {
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdx_proxy;
            deviceRecvBuffer = GPURecvBufferTurb;
        }
        else if ("turb::dqdy" == name)
        {
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdy_proxy;
            deviceRecvBuffer = GPURecvBufferTurb;
        }
        else if ("turb::dqdz" == name)
        {
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdz_proxy;
            deviceRecvBuffer = GPURecvBufferTurb;
        }
        else
        {
            cout << "Error:" << name << " is not included in CallGPUDownloadInterfaceValue" << endl;
            exit(1);
        }

        CallGPUBufferToField(deviceFieldVar, deviceRecvBuffer, nTotal, neqn, nIFace, offsetBuffer);
    }

    void CallGPUBufferToField(RFloat *deviceFieldVar, const double *deviceRecvBuffer, const int nTotal, const int neqn,
                              const int nIFace, const int offsetBuffer)
    {
        using namespace GPUMemory;
        using namespace GPUGeomVariables;
        int gridSize, blockSize, loopLength;
        loopLength         = nIFace;
        int blockSizeLimit = 0;
        blockSize          = 1;
        gridSize           = 1;
        KernelLaunchPara((void *)GPUBufferToField, loopLength, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUBufferToField, 0, gridSize, blockSize);
#endif
        GPUBufferToField<<<gridSize, blockSize>>>(deviceFieldVar, deviceRecvBuffer, d_interFace2BoundaryFace,
                                                  d_right_cell_of_face, nTotal, neqn, nIFace, offsetBuffer);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUBufferToField(RFloat *deviceFieldVar, const double *deviceRecvBuffer,
                                     const int *interFace2BoundaryFace, const int *rightCellOfFace, const int nTotal,
                                     const int neqn, const int nIFace, const int offset)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int iFace, targetCell, m;
        for (iFace = threadID; iFace < nIFace; iFace += blockDim.x * gridDim.x)
        {
            targetCell = rightCellOfFace[interFace2BoundaryFace[iFace]];
            for (m = 0; m < neqn; ++m)
            {
                //! f[m][targetCell] = fg[m][iFace];
                deviceFieldVar[m * nTotal + targetCell] = deviceRecvBuffer[offset + m * nIFace + iFace];
                /*
      //!add for test
      if ((targetCell == 0)&&(m == 0)) printf("iFace = %d, deviceRecvBuffer =
      %.30e\n", iFace, deviceRecvBuffer[offset + m*nIFace+iFace]); if ((m ==
      1)&&(iFace == 0)) { printf("sourceCell = %d, deviceFieldVar = %.30e\n",
      sourceCell, deviceFieldVar[m * nTotal + sourceCell]);
      }
      //!add end
      */
            }
        }
    }

    void CallGPUUploadInterfaceValueLargeBuffer(const int nTotal, const int neqn, const int nIFace, const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        int     nameNSID, nameTurbID;
        RFloat *deviceFieldVar;
        double *deviceSendBuffer;
        int     offsetBuffer;
        /*
  //!add for test
  nameNSID = GetIndexBufferNS(name);
  if (nameNSID == 666) {
          nameTurbID = GetIndexBufferTurb(name);
          if (nameTurbID == 666) {
                  printf("Error: %s is not in NS buffer or Turb buffer\n",
  name.c_str()); exit(0);
          }
          printf("%s in TurbSolver by %d\n", name.c_str(), nameTurbID);
  }
  else {
          printf("%s in NSSolver by %d\n", name.c_str(), nameNSID);
  }
  //!add end
  */
        if ("q" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_q_ns;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("t" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_t_proxy;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("turb::q" == name)
        {
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            deviceFieldVar   = d_q_turb_proxy;
            deviceSendBuffer = GPUSendBufferTurb;
        }
        else if ("dqdx" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dqdx_proxy;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("dqdy" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dqdy_proxy;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("dqdz" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dqdz_proxy;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("limit" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_limit;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("dtdx" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dtdx_proxy;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("dtdy" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dtdy_proxy;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("dtdz" == name)
        {
            nameNSID         = GetIndexBufferNS(name);
            offsetBuffer     = offsetBufferNS[nameNSID];
            deviceFieldVar   = d_dtdz_proxy;
            deviceSendBuffer = GPUSendBufferNS;
        }
        else if ("turb::dqdx" == name)
        {
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdx_proxy;
            deviceSendBuffer = GPUSendBufferTurb;
        }
        else if ("turb::dqdy" == name)
        {
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdy_proxy;
            deviceSendBuffer = GPUSendBufferTurb;
        }
        else if ("turb::dqdz" == name)
        {
            nameTurbID       = GetIndexBufferTurb(name);
            offsetBuffer     = offsetBufferTurb[nameTurbID];
            deviceFieldVar   = d_dq_turbdz_proxy;
            deviceSendBuffer = GPUSendBufferTurb;
        }
        else
        {
            cout << "Error:" << name << " is not included in CallGPUUploadInterfaceValue" << endl;
            exit(1);
        }

        CallGPUFieldToBuffer(deviceSendBuffer, deviceFieldVar, nTotal, neqn, nIFace, offsetBuffer);
    }

    void CallGPUFieldToBuffer(double *DeviceSendBuffer, const RFloat *deviceFieldVar, const int nTotal, const int neqn,
                              const int nIFace, const int offset)
    {
        using namespace GPUMemory;
        using namespace GPUGeomVariables;
        int gridSize, blockSize, loopLength;
        loopLength         = nIFace;
        int blockSizeLimit = 0;
        blockSize          = 1;
        gridSize           = 1;
        KernelLaunchPara((void *)GPUFieldToBuffer, loopLength, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUFieldToBuffer, 0, gridSize, blockSize);
#endif
        GPUFieldToBuffer<<<gridSize, blockSize>>>(DeviceSendBuffer, deviceFieldVar, d_interFace2BoundaryFace,
                                                  d_left_cell_of_face, nTotal, neqn, nIFace, offset);
//! GPUFieldToBuffer<<<1, 1>>>(DeviceSendBuffer,
//! deviceFieldVar,d_interFace2BoundaryFace, d_left_cell_of_face, nTotal, neqn,
//! nIFace, offset);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUFieldToBuffer(double *DeviceSendBuffer, const RFloat *deviceFieldVar,
                                     const int *interFace2BoundaryFace, const int *leftCellOfFace, const int nTotal,
                                     const int neqn, const int nIFace, const int offset)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int iFace, sourceCell, m;
        for (iFace = threadID; iFace < nIFace; iFace += blockDim.x * gridDim.x)
        {
            sourceCell = leftCellOfFace[interFace2BoundaryFace[iFace]];
            for (m = 0; m < neqn; ++m)
            {
                //! fg[m * nIFace + iFace] = f[m * nTotal + sourceCell];
                DeviceSendBuffer[offset + m * nIFace + iFace] = deviceFieldVar[m * nTotal + sourceCell];
                /*
      //!add for test
      if ((m == 1)&&(iFace == 0)) {
              printf("sourceCell = %d, deviceFieldVar = %.30e\n", sourceCell,
      deviceFieldVar[m * nTotal + sourceCell]);
      }
      //!add end
      */
            }
        }
    }

    void CallGPUUploadInterfaceValue(const int nTotal, const int neqn, const int nIFace, const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        RFloat *deviceIFVar;
        RFloat *deviceFieldVar;
        if ("q" == name)
        {
            deviceIFVar    = d_fg_send_q;
            deviceFieldVar = d_q_ns;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("t" == name)
        {
            deviceIFVar    = d_fg_send_t;
            deviceFieldVar = d_t_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::q" == name)
        {
            deviceIFVar    = d_fg_send_qTurb;
            deviceFieldVar = d_q_turb_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dqdx" == name)
        {
            deviceIFVar    = d_fg_send_dqdx;
            deviceFieldVar = d_dqdx_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dqdy" == name)
        {
            deviceIFVar    = d_fg_send_dqdy;
            deviceFieldVar = d_dqdy_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dqdz" == name)
        {
            deviceIFVar    = d_fg_send_dqdz;
            deviceFieldVar = d_dqdz_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("limit" == name)
        {
            deviceIFVar    = d_fg_send_limit;
            deviceFieldVar = d_limit;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdx" == name)
        {
            deviceIFVar    = d_fg_send_dtdx;
            deviceFieldVar = d_dtdx_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdy" == name)
        {
            deviceIFVar    = d_fg_send_dtdy;
            deviceFieldVar = d_dtdy_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("dtdz" == name)
        {
            deviceIFVar    = d_fg_send_dtdz;
            deviceFieldVar = d_dtdz_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdx" == name)
        {
            deviceIFVar    = d_fg_send_dqTurbdx;
            deviceFieldVar = d_dq_turbdx_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdy" == name)
        {
            deviceIFVar    = d_fg_send_dqTurbdy;
            deviceFieldVar = d_dq_turbdy_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else if ("turb::dqdz" == name)
        {
            deviceIFVar    = d_fg_send_dqTurbdz;
            deviceFieldVar = d_dq_turbdz_proxy;
            CallGPUFieldVarToIFVar(deviceIFVar, deviceFieldVar, nTotal, neqn, nIFace);
        }
        else
        {
            cout << "Error:" << name << " is not included in CallGPUUploadInterfaceValue" << endl;
            exit(1);
        }

        /*
  //!printf("program is running in CallGPUUploadInterfaceValue, name = %s\n",
  name); cout<<"program is running in CallGPUUploadInterfaceValue, name =
  "<<name<<" iGhostLayer = "<<iGhostLayer<<"in globalRank "<<globalRank<<endl;
  if ("turb::q" == name) printf("name is equal to turb::q in globalRank %d\n",
  globalRank); else printf("name is not equal to turb::q in globalRank %d\n",
  globalRank);
  */
        /*
  //!error: expression must be an integral constant expression
  switch (name) {
          case "q":
                  cout<<"program is running in CallGPUUploadInterfaceValue, name
  = "<<name<<" iGhostLayer = "<<iGhostLayer<<"in globalRank "<<globalRank<<endl;
                  break;
          case "t":
                  cout<<"program is running in CallGPUUploadInterfaceValue, name
  = "<<name<<" iGhostLayer = "<<iGhostLayer<<"in globalRank "<<globalRank<<endl;
                  break;
          case "turb::q":
                  cout<<"program is running in CallGPUUploadInterfaceValue, name
  = "<<name<<" iGhostLayer = "<<iGhostLayer<<"in globalRank "<<globalRank<<endl;
                  break;
          default:
                  cout<<"name is NOT q, t or turb::q"<<endl;
                  break;
  }
  */
    }
    void CallGPUIFVarToFieldVarInterpoint(RFloat *deviceIFVar, RFloat *deviceFieldVar, const int nTotal, const int neqn,
                                          const int nIPoint)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nIPoint;
        KernelLaunchPara((void *)GPUIFVarToFieldVarInterpoint, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUIFVarToFieldVarInterpoint, 0, gridSize, blockSize);
#endif
        GPUIFVarToFieldVarInterpoint<<<gridSize, blockSize>>>(deviceIFVar, deviceFieldVar, d_interPoint2GlobalPoint,
                                                              nTotal, neqn, nIPoint);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    void CallGPUIFVarToFieldVarHToDAsync(RFloat *deviceIFVar, RFloat *deviceFieldVar, const int nTotal, const int neqn,
                                         const int nIFace)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nIFace;
        KernelLaunchPara((void *)GPUIFVarToFieldVar, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUIFVarToFieldVar, 0, gridSize, blockSize);
#endif
        GPUIFVarToFieldVar<<<gridSize, blockSize, 0, dataTransferHToD>>>(
            deviceIFVar, deviceFieldVar, d_interFace2BoundaryFace, d_right_cell_of_face, nTotal, neqn, nIFace);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    void CallGPUIFVarToFieldVar(RFloat *deviceIFVar, RFloat *deviceFieldVar, const int nTotal, const int neqn,
                                const int nIFace)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nIFace;
        KernelLaunchPara((void *)GPUIFVarToFieldVar, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUIFVarToFieldVar, 0, gridSize, blockSize);
#endif
        GPUIFVarToFieldVar<<<gridSize, blockSize>>>(deviceIFVar, deviceFieldVar, d_interFace2BoundaryFace,
                                                    d_right_cell_of_face, nTotal, neqn, nIFace);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUFieldVarToIFVarInterpoint(RFloat *deviceIFVar, const RFloat *deviceFieldVar, const int nTotalNode,
                                          const int neqn, const int nIPoint)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nIPoint;
        KernelLaunchPara((void *)GPUFieldVarToIFVarInterpoint, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUFieldVarToIFVarInterpoint, 0, gridSize, blockSize);
#endif
        GPUFieldVarToIFVarInterpoint<<<gridSize, blockSize>>>(deviceIFVar, deviceFieldVar, d_interPoint2GlobalPoint,
                                                              nTotalNode, neqn, nIPoint);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUFieldVarToIFVar(RFloat *deviceIFVar, const RFloat *deviceFieldVar, const int nTotal, const int neqn,
                                const int nIFace)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nIFace;
        KernelLaunchPara((void *)GPUFieldVarToIFVar, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUFieldVarToIFVar, 0, gridSize, blockSize);
#endif
        GPUFieldVarToIFVar<<<gridSize, blockSize>>>(deviceIFVar, deviceFieldVar, d_interFace2BoundaryFace,
                                                    d_left_cell_of_face, nTotal, neqn, nIFace);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallUploadInterpointValueDeviceToHost(RFloat **fg, const int neqn, const int nIPoint, const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        RFloat *deviceIFVar;
        size_t  sizeIF = neqn * nIPoint * sizeof(RFloat);
        RFloat *h_fg0  = fg[0];
        if ("qnode" == name)
        {
            deviceIFVar = d_fg_send_qNode;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("tnode" == name)
        {
            deviceIFVar = d_fg_send_tNode;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("qTurbNode" == name)
        {
            deviceIFVar = d_fg_send_qTurbNode;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else
        {
            cout << "Error:" << name << "is not included in CallUploadInterpointValueDeviceToHost" << endl;
            exit(1);
        }
    }

    void CallUploadInterfaceValueDeviceToHost(RFloat **fg, const int neqn, const int nIFace, const string &name)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        RFloat *deviceIFVar;
        size_t  sizeIF = neqn * nIFace * sizeof(RFloat);
        RFloat *h_fg0  = fg[0];
        if ("q" == name)
        {
            deviceIFVar = d_fg_send_q;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("t" == name)
        {
            deviceIFVar = d_fg_send_t;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("turb::q" == name)
        {
            deviceIFVar = d_fg_send_qTurb;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("dqdx" == name)
        {
            deviceIFVar = d_fg_send_dqdx;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("dqdy" == name)
        {
            deviceIFVar = d_fg_send_dqdy;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("dqdz" == name)
        {
            deviceIFVar = d_fg_send_dqdz;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("limit" == name)
        {
            deviceIFVar = d_fg_send_limit;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("dtdx" == name)
        {
            deviceIFVar = d_fg_send_dtdx;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("dtdy" == name)
        {
            deviceIFVar = d_fg_send_dtdy;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("dtdz" == name)
        {
            deviceIFVar = d_fg_send_dtdz;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("turb::dqdx" == name)
        {
            deviceIFVar = d_fg_send_dqTurbdx;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("turb::dqdy" == name)
        {
            deviceIFVar = d_fg_send_dqTurbdy;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else if ("turb::dqdz" == name)
        {
            deviceIFVar = d_fg_send_dqTurbdz;
            HANDLE_API_ERR(cudaMemcpy(h_fg0, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
        }
        else
        {
            cout << "Error:" << name << " is not included in CallUploadInterfaceValueDeviceToHost" << endl;
            exit(1);
        }
    }
    __global__ void GPUFieldVarToIFVar(RFloat *fg, const RFloat *f, const int *interFace2BoundaryFace,
                                       const int *leftCellOfFace, int nTotal, int neqn, int nIFace)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int iFace    = 0;
        for (iFace = threadID; iFace < nIFace; iFace += blockDim.x * gridDim.x)
        {
            int sourceCell = leftCellOfFace[interFace2BoundaryFace[iFace]];
            int m;
            for (m = 0; m < neqn; ++m)
            {
                fg[m * nIFace + iFace] = f[m * nTotal + sourceCell];
            }
        }
    }
    __global__ void GPUIFVarToFieldVar(const RFloat *fg, RFloat *f, const int *interFace2BoundaryFace,
                                       const int *rightCellOfFace, int nTotal, int neqn, int nIFace)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int iFace    = 0;
        for (iFace = threadID; iFace < nIFace; iFace += blockDim.x * gridDim.x)
        {
            int sourceCell = rightCellOfFace[interFace2BoundaryFace[iFace]];
            int m;
            for (m = 0; m < neqn; ++m)
            {
                f[m * nTotal + sourceCell] = fg[m * nIFace + iFace];
            }
        }
    }
    __global__ void GPUFieldVarToIFVarInterpoint(RFloat *fg, const RFloat *f, const int *interPoint2GlobalPoint,
                                                 const int nTotalNode, const int neqn, const int nIPoint)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int iPoint;
        int sourcePoint;
        for (iPoint = threadID; iPoint < nIPoint; iPoint += blockDim.x * gridDim.x)
        {
            sourcePoint = interPoint2GlobalPoint[iPoint];
            //! just for test
            //! if (sourcePoint != iPoint) printf("sourcePoint = %d, iPoint = %d\n",
            //! sourcePoint, iPoint);
            //! test end
            for (int m = 0; m < neqn; ++m)
            {
                fg[m * nIPoint + iPoint] = f[m * nTotalNode + sourcePoint];
            }
        }
    }
    __global__ void GPUIFVarToFieldVarInterpoint(const RFloat *fg, RFloat *f, const int *interPoint2GlobalPoint,
                                                 const int nTotalNode, const int neqn, const int nIPoint)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int iPoint;
        int sourcePoint;
        for (iPoint = threadID; iPoint < nIPoint; iPoint += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < neqn; ++m)
            {
                //! f[m * nTotalNode + iPoint] = fg[m * nIPoint + iPoint];
                f[m * nIPoint + iPoint] = fg[m * nIPoint + iPoint];
            }
        }
    }

    __global__ void GPUMPIDataCompress(const int offsetVar, const int neqn, const int nIFaceTotal,
                                       const int offsetMPIDataNS, const int nIFaceOfNeighbor, const int offsetNgbVar,
                                       const int offsetFaceIndex, const int *faceIndexForSend,
                                       const double *GPUSendBufferNS, double *GPUMPIDataSendNS)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int faceID;
        int equationID;
        int ngbZoneFaceID;
        int largeStorageID;
        int faceIndexLargeStore;
        for (faceID = threadID; faceID < nIFaceOfNeighbor; faceID += blockDim.x * gridDim.x)
        {
            /*
                            ngbZoneFaceID = offsetMPIDataNS + offsetNgbVar +
       faceID; faceIndexLargeStore = faceIndexForSend[offsetFaceIndex+faceID];
                            for (equationID = 0; equationID < neqn;
       equationID++) { largeStorageID = offsetVar + equationID * nIFaceTotal +
       faceIndexLargeStore; GPUMPIDataSendNS[ngbZoneFaceID] =
       GPUSendBufferNS[largeStorageID];
                            }
    */

            faceIndexLargeStore = faceIndexForSend[offsetFaceIndex + faceID];
            for (equationID = 0; equationID < neqn; equationID++)
            {
                ngbZoneFaceID  = offsetMPIDataNS + offsetNgbVar + equationID * nIFaceOfNeighbor + faceID;
                largeStorageID = offsetVar + equationID * nIFaceTotal + faceIndexLargeStore;
                GPUMPIDataSendNS[ngbZoneFaceID] = GPUSendBufferNS[largeStorageID];
            }
        }
    }

    //!__global__ void GPUMPIDataDecompress(const int offsetVar, const int neqn,
    //!const int nIFaceTotal, const int offsetMPIDataNS, const int nIFaceOfNeighbor,
    //!const int offsetNgbVar, const int offsetFaceIndex, const int *
    //!faceIndexForRecv, const double * GPUMPIDataRecvNS, double * GPURecvBufferNS){
    __global__ void GPUMPIDataDecompress(const int offsetVar, const int neqn, const int nIFaceTotal,
                                         const int offsetMPIData, const int nIFaceOfNeighbor, const int offsetNgbVar,
                                         const int offsetFaceIndex, const int *faceIndexForRecv,
                                         const double *GPUMPIDataRecv, double *GPURecvBuffer)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int faceID;
        int equationID;
        int ngbZoneFaceID;
        int largeStorageID;
        int faceIndexLargeStore;
        for (faceID = threadID; faceID < nIFaceOfNeighbor; faceID += blockDim.x * gridDim.x)
        {
            //! ngbZoneFaceID = offsetMPIDataNS + offsetNgbVar + faceID;
            faceIndexLargeStore = faceIndexForRecv[offsetFaceIndex + faceID];
            for (equationID = 0; equationID < neqn; equationID++)
            {
                //! ngbZoneFaceID = offsetMPIDataNS + offsetNgbVar +
                //! equationID*nIFaceOfNeighbor + faceID;
                ngbZoneFaceID  = offsetMPIData + offsetNgbVar + equationID * nIFaceOfNeighbor + faceID;
                largeStorageID = offsetVar + equationID * nIFaceTotal + faceIndexLargeStore;
                //! GPURecvBufferNS[largeStorageID] = GPUMPIDataRecvNS[ngbZoneFaceID];
                GPURecvBuffer[largeStorageID] = GPUMPIDataRecv[ngbZoneFaceID];
            }
        }
    }
} //! namespace GPUKernels
