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
//! @file      GPUMPICommunication.h
//! @brief     Communication functions.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
#include "Geo_Interface.h"
using namespace PHSPACE;
namespace GPUKernels
{
    //!GPUSMALLBUFFER
    void CompressInterfaceTurbDataSmallGPUBufferOverlap(const int nameID, const int nTotal, const int neqn,
                                                        const int nIFace, const int nTotalNgbZones,
                                                        InterfaceInfo *iinfo, RFloat **GPUDataSendSmallBuffer);
    void CompressInterfaceNSDataSmallGPUBufferOverlap(const int nameID, const int nTotal, const int neqn,
                                                      const int nIFace, const int nTotalNgbZones, InterfaceInfo *iinfo,
                                                      RFloat **GPUDataSendSmallBuffer);
    void CompressInterfaceNSDataSmallGPUBufferOverlap(const int nameID, const int nTotal, const int neqn,
                                                      const int nIFace, const int nTotalNgbZones, InterfaceInfo *iinfo,
                                                      RFloat **GPUDataSendSmallBuffer);
    void CallCompressInterfaceDataSmallGPUBufferOverlap(const int nTotal, const int neqn, const int nIFace,
                                                        const string &name, InterfaceInfo *iinfo);
    //!MPI_Allreduce for globalMinDt
    void ComputeMinimumDtAllreduce(double &globalMinDt);
    void SetGlobalMinDt(double &globalMinDt);
    //!Decompuresss MPI data into interface
    void DecompressInterfaceNSDataGPUSmallBuffer(const int solverID, InterfaceInfo *iinfo);
    void DecompressInterfaceDataSmallGPUBuffer(const int solverID, InterfaceInfo *iinfo);

    __global__ void GPUMPIDataDecompressSmallBuf(const int offsetVar, const int neqn, const int nIFaceTotal,
                                                 const int nIFaceOfNeighbor, const int offsetFaceIndex,
                                                 const int *faceIndexForRecv, const double *GPUMPIDataRecvSmall,
                                                 double *GPURecvBufferLarge);

    void DecompressInterpolateData(const int solverID, InterfaceInfo *iinfo);
    void DecompressInterpolateDataMPIOverlap(const int solverID, InterfaceInfo *iinfo);
    void DecompressInterpolateDataV2(const int solverID, InterfaceInfo *iinfo);
    void DecompressInterpointData(const int solverID, InterpointInformation *iinfo);
    //!Transfer interface data by CUDA_AWARE_MPI or traditional MPI method
    void TransferInterfaceDataNSMPI(const int solverID, InterfaceInfo *iinfo);
    void TransferInterfaceDataTurbMPI(const int solverID, InterfaceInfo *iinfo);
    void TransferInterfaceDataMPI(const int solverID, InterfaceInfo *iinfo);
    //!Compress interface data by zone connectivity informaion
    //!solverID determines NS or Turb
    void CompressInterfaceDataSmallGPUBuffer(const int solverID, InterfaceInfo *iinfo);
    void CompressInterfaceDataNSGPUSmallBuffer(const int solverID, InterfaceInfo *iinfo);
    void CompressInterfaceDataTurbGPUSmallBuffer(const int solverID, InterfaceInfo *iinfo);

    __global__ void GPUMPIDataCompressSmallBuf(const int offsetVar, const int neqn, const int nIFaceTotal,
                                               const int nIFaceOfNeighbor, const int offsetFaceIndex,
                                               const int *faceIndexForSend, const double *GPUSendBufferLarge,
                                               double *GPUMPIDataSendSmall);

    void CompressAndTransferInterpolateData(const int solverID, InterfaceInfo *iinfo);
    void CompressAndTransferInterpolateDataV2(const int solverID, InterfaceInfo *iinfo);
    void CompressAndTransferInterpointData(const int solverID, InterpointInformation *iinfo);
    void CATOverlapInterfaceIsendIrecvNS(InterfaceInfo *iinfo);
    void CATSeparateInterfaceIsendIrecvNS(InterfaceInfo *iinfo);
    void CATInterfaceSendrecvNS(InterfaceInfo *iinfo);
    void CATOverlapInterfaceIsendIrecvTurb(InterfaceInfo *iinfo);
    void CATSeparateInterfaceIsendIrecvTurb(InterfaceInfo *iinfo);
    void CATInterfaceSendrecvTurb(InterfaceInfo *iinfo);

    void CATOverlapInterpointIsendIrecvNS(InterpointInformation *iinfo);
    void CATSeparateInterpointIsendIrecvNS(InterpointInformation *iinfo);
    void CATOverlapInterpointIsendIrecvTurb(InterpointInformation *iinfo);
    void CATSeparateInterpointIsendIrecvTurb(InterpointInformation *iinfo);
    void CATInterpointSendrecvNS(InterpointInformation *iinfo);
    void CATInterpointSendrecvTurb(InterpointInformation *iinfo);
    //!make cuda stream synchronize for a cuda stream
    void CallCUDAStreamSynchronize(cudaStream_t streamProcess);
    void CallDeviceUploadInterfaceValueHostLargeBuffer(const int nTotal, const int neqn, const int nIFace,
                                                       const string &name);
    //!translate data from GPUSendInterpointBufferNS into HostSendInterpointBufferNS. Similarly, translate data from GPUSendInterpointBufferTurb into HostSendInterpointBufferTurb. GPUSendLargeBufferToHostSendLargeBufferInterpoint is used in CFDSolver::UpdateInterpointData.
    void GPUSendLargeBufferToHostSendLargeBufferInterpoint(const int solverID, InterpointInformation *iinfo);
    //!translate data from HostSendInterpointBufferNS and HostSendInterpointBufferTurb into sendData buffer. nameNSBuffer, offsetBufferNS, and dimBufferNS (also Turb) are required. HostSendLargeBufferToOrgSendBufferInterpoint is used in CFDSolver::UpdateInterpointData.
    void HostSendLargeBufferToOrgSendBufferInterpoint(const int solverID, InterpointInformation *iinfo);
    void OrgRecvBufferToHostRecvLargeBufferInterpoint(const int solverID, InterpointInformation *iinfo);
    void HostRecvLargeBufferToGPURecvSmallBufferInterpoint(const int solverID, InterpointInformation *iinfo);
    void HostRecvLargeBufferToGPURecvLargeBufferInterpoint(const int solverID, InterpointInformation *iinfo);
    //!translate data from GPUSendInterfaceBufferNS into HostSendInterfaceBufferNS. Similarly, translate data from GPUSendInterfaceBufferTurb into HostSendInterfaceBufferTurb. The function is used in CFDSolver::UpdateInterfaceData
    void GPUSendLargeBufferToHostSendLargeBuffer(const int solverID, InterfaceInfo *iinfo);
    void HostSendLargeBufferToOrgSendBuffer(const int solverID, InterfaceInfo *iinfo);
    void OrgRecvBufferToHostRecvLargeBuffer(const int solverID, InterfaceInfo *iinfo);
    void HostRecvLargeBufferToGPURecvSmallBuffer(const int solverID, InterfaceInfo *iinfo);
    void HostRecvLargeBufferToGPURecvLargeBuffer(const int solverID, InterfaceInfo *iinfo);
    void CallGPUUploadInterfaceValue(const int nTotal, const int neqn, const int nIFace, const string &name);
    void CallGPUUploadInterfaceValueLargeBuffer(const int nTotal, const int neqn, const int nIFace, const string &name);
    void CallGPUFieldToBuffer(double *DeviceSendBuffer, const RFloat *deviceFieldVar, const int nTotal, const int neqn,
                              const int nIFace, const int offset);
    void CallGPUFieldVarToIFVar(RFloat *deviceIFVar, const RFloat *deviceFieldVar, const int nTotal, const int neqn,
                                const int nIFace);

    void CallUploadInterpointValueDeviceToHost(RFloat **fg, const int neqn, const int nIPoint, const string &name);
    void CallUploadInterfaceValueDeviceToHost(RFloat **fg, const int neqn, const int nIFace, const string &name);
    //!void CallGPUDownloadInterfaceValue(RFloat ** fg, const int neqn, const int nIFace, const string & name, const int iGhostLayer);
    void CallGPUDownloadInterfaceValueSmallBuffer(const int nTotal, const int neqn, const int nIFace,
                                                  const string &name);
    void CallGPUDownloadInterfaceValueLargeBuffer(RFloat **fg, const int nTotal, const int neqn, const int nIFace,
                                                  const string &name);
    void CallGPUBufferToField(RFloat *deviceFieldVar, const double *deviceRecvBuffer, const int nTotal, const int neqn,
                              const int nIFace, const int offsetBuffer);
    __global__ void GPUBufferToField(RFloat *deviceFieldVar, const double *deviceRecvBuffer,
                                     const int *interFace2BoundaryFace, const int *rightCellOfFace, const int nTotal,
                                     const int neqn, const int nIFace, const int offset);

    void CallGPUDownloadInterfaceValue(RFloat **fg, const int nTotal, const int neqn, const int nIFace,
                                       const string &name);
    void CallGPUDownloadInterpointValue(RFloat **fg, const int nTotalNode, const int neqn, const int nIPoint,
                                        const string &name);
    void CallGPUDownloadInterpointValueSmallBuffer(const int nTotalNode, const int neqn, const int nIPoint,
                                                   const string &name);
    void CallGPUDownloadInterpointValueLargeBuffer(const int nTotalNode, const int neqn, const int nIPoint,
                                                   const string &name);
    void CallGPUInterpointBufferToField(RFloat *deviceFieldVar, const double *deviceRecvBuffer, const int nTotalNode,
                                        const int neqn, const int nIPoint, const int offsetBuffer);
    __global__ void GPUInterpointBufferToField(RFloat *deviceFieldVar, const double *deviceBuffer, const int nTotalNode,
                                               const int neqn, const int nIPoint, const int offset);
    void CallGPUIFVarToFieldVarHToDAsync(RFloat *deviceIFVar, RFloat *deviceFieldVar, const int nTotal, const int neqn,
                                         const int nIFace);
    void CallGPUIFVarToFieldVar(RFloat *deviceIFVar, RFloat *deviceFieldVar, const int nTotal, const int neqn,
                                const int nIFace);
    void CallGPUUploadInterpointValue(const int nTotalNode, const int neqn, const int nIPoint, const string &name);
    void CallDeviceUploadInterpointValueHostLargeBuffer(const int nTotalNode, const int neqn, const int nIPoint,
                                                        const string &name);
    void CallGPUUploadInterpointValueLargeBuffer(const int nTotalNode, const int neqn, const int nIPoint,
                                                 const string &name);
    void CallGPUFieldToInterpointBuffer(const RFloat *deviceFieldVar, double *deviceSendBuffer, const int nTotalNode,
                                        const int neqn, const int nIPoint, const int offsetBuffer);
    __global__ void GPUFieldToInterpointBuffer(const RFloat *deviceFieldVar, double *deviceSendBuffer,
                                               const int *interPoint2GlobalPoint, const int nTotalNode, const int neqn,
                                               const int nIPoint, const int offset);
    void CallGPUFieldVarToIFVarInterpoint(RFloat *deviceIFVar, const RFloat *deviceFieldVar, const int nTotalNode,
                                          const int neqn, const int nIPoint);
    void CallGPUIFVarToFieldVarInterpoint(RFloat *deviceIFVar, RFloat *deviceFieldVar, const int nTotal, const int neqn,
                                          const int nIPoint);

    __global__ void GPUFieldVarToIFVar(RFloat *fg, const RFloat *fm, const int *interFace2BoundaryFace,
                                       const int *leftCellOfFace, int nTotal, int neqn, int nIFace);
    __global__ void GPUIFVarToFieldVar(const RFloat *fg, RFloat *f, const int *interFace2BoundaryFace,
                                       const int *rightCellOfFace, int nTotal, int neqn, int nIFace);
    __global__ void GPUFieldVarToIFVarInterpoint(RFloat *fg, const RFloat *f, const int *interPoint2GlobalPoint,
                                                 const int nTotalNode, const int neqn, const int nIPoint);
    __global__ void GPUIFVarToFieldVarInterpoint(const RFloat *fg, RFloat *f, const int *interPoint2GlobalPoint,
                                                 const int nTotalNode, const int neqn, const int nIPoint);
    __global__ void GPUFieldToBuffer(double *DeviceSendBuffer, const RFloat *deviceFieldVar,
                                     const int *interFace2BoundaryFace, const int *rightCellOfFace, const int nTotal,
                                     const int neqn, const int nIFace, const int offset);
    //!set GPUMPIDataSendNS by GPUSendBufferNS. In fact, it is regarded as a knid of reorder.offsetVar is offset of variabels in GPUSendBufferNS.
    //!neqn is dimData of variable. nIFaceTotal is nIFace, which is the total faces of boudnary in iZone.
    //!offsetMPIDataNS is offset of the component of jzone in GPUMPIDataSendNS.
    //!nIFaceOfNeighbor is faces between iZone and jzone.
    //!offsetNgbVar is offset of the variables in the component of jzone, which should be used with offsetMPIDataNS.
    //!offsetFaceIndex is offset of faceIndexForSend for jzone. faceIndexForSend is the reflecting between faces in the component of jzone in GPUMPIDataSendNS and faces in neighbour between iZone and jzone.
    //!GPUSendBufferNS is the large storage for variables of NS solver.
    //!GPUMPIDataSendNS is the send buffer for NS Solver.
    __global__ void GPUMPIDataCompress(const int offsetVar, const int neqn, const int nIFaceTotal,
                                       const int offsetMPIDataNS, const int nIFaceOfNeighbor, const int offsetNgbVar,
                                       const int offsetFaceIndex, const int *faceIndexForSend,
                                       const double *GPUSendBufferNS, double *GPUMPIDataSendNS);
    __global__ void GPUMPIDataDecompress(const int offsetVar, const int neqn, const int nIFaceTotal,
                                         const int offsetMPIDataNS, const int nIFaceOfNeighbor, const int offsetNgbVar,
                                         const int offsetFaceIndex, const int *faceIndexForSend,
                                         const double *GPUMPIDataRecvNS, double *GPURecvBufferNS);
} //! namespace GPUKernels
