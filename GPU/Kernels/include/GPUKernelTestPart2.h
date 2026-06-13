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
//! @file      GPUKernelTestPart2.h
//! @brief     Kernel test function for GPU results.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

//!#include <stdio.h>
#include <iostream>
//!#include "FYGrid.h"
//!#include "FYUtil.h"
//!#include <cstdio>
#include "Precision.h"
#include "cudaErrorHandle.h"

using namespace PHSPACE;
//! TestGPUCompNodeVarByGradientLoop1 is just for the test of initialization of
//! d_q_n_double and d_n_count_tmp, which is the same for all of conditions
void TestGPUCompNodeVarByGradientLoop1(const string d_gradient_field_proxy, const int nTotalNode, RFloat *org_q_n,
                                       const int *org_n_count);
void TestGPUNodePosiFaceAllocCopy(const long long int *h_nodePosiFace, const long long int *d_nodePosiFace,
                                  const int nTotalNodePosiFace);
void TestGPUBoundaryTypeAllocCopy(const int *h_boundaryType, const int *d_boundaryType, const int nBoundFace);
//! TestGPUCompNodeVarByGradientLoop2 is more general, which is used for the
//! other loops in CompNodeVar. Thus, many conditions should be considered in the
//! test
void TestGPUCompNodeVarByGradientLoop2(const string d_gradient_field_proxy, const int nTotalNode, RFloat *org_q_n,
                                       const int *org_n_count);
void TestGPUCompNodeVarByGradientFinal(const string d_gradient_field_proxy, const int nTotalNode, RFloat *org_q_n);
void TestGPUCompGradientGGNode(const string q_name, const int index, const int nTotalCell, const int nBoundFace,
                               RFloat *org_dqdx, RFloat *org_dqdy, RFloat *org_dqdz);
void TestGPUCompGradientGGCell(const string q_name, const int index, const int nTotalCell, const int nBoundFace,
                               RFloat *org_dqdx, RFloat *org_dqdy, RFloat *org_dqdz);
void TestGPUStoreBoundGrad(const int nBoundFace, const int nTVar, RFloat **bdqx, RFloat **bdqy, RFloat **bdqz);
void TestGPUFaceWeight(const int nst, const int ned, const RDouble *org_deltl, const RDouble *org_deltr);
void TestGPUViscousCoefCopy(const RFloat *org_visl, const RFloat *org_vist, const int nTotal);
void TestGPUGetVisFaceValueMult(const int nst, const int ned, const RFloat *mul, const RFloat *mut);
void TestGPUGetVisFaceValuePrimTmKcp(const int nst, const int ned, RFloat **org_prim, RFloat **org_tmid,
                                     const RFloat *org_kcp);
void TESTGPUCompVisfluxTEST(const int nst, RFloat **org_flux, const int nl);
void TestCompInviscid(const int nst, const int ned, const int n_turb, RFloat **org_flux_turb,
                      RFloat **org_flux_sub_turb, const int flag);
void TestGPUFlux_turb();
void TestGPULoadFlux(const int nst, const int ned, const int nl, const int nTotalCell, const int nBoundFace,
                     RFloat **org_res);
void TestGPUResCopy(RFloat **res, const int nl, const int nTotal);
void TestGPUGetTurbVisFaceValueMulMut(const int nst, const int ned, const RFloat *org_mul, const RFloat *org_mut);
void TestGPUGetTurbVisFaceValuePrim(const int nst, const int ned, const int n_turb, RFloat **org_prim);
void TestGPUGetTurbVisFaceValueMlt(const int nst, const int ned, const int n_turb, RFloat **org_mlt);
void TestGPUParametersMemAllocCopy(const int neqn, RDouble *org_turboo);
void TestCallGPUCompTurbVisfluxTEST(const int nst, const int ned, const int n_turb, RFloat **org_flux,
                                    RFloat **org_flux_sub);
void TestGPUQ_turbCopy(RFloat **org_q_turb, const int nTotal, const int n_turb);
void TestGPULoadTurbFlux(RFloat **res_turb, const int nTotalCell, const int nBoundFace, const int n_turb);
void TestGPUResTurbCopy(RFloat **res_turb, const int nTotalCell, const int nBoundFace, const int n_turb);
void TestGPUSpecTurbCell(const int n_turb, const int nTotalCell, const int nBoundFace, RFloat **spec_turb);
void TestGPUSpecTurbMatTurbFaces(const int n_turb, const int nTotalFace, const int nBoundFace, const int nTotalCell,
                                 RFloat **spec_turb, RFloat **mat_turbl, RFloat **mat_turbr);
void TestGPUSourceFlux_1eq_Original(RFloat **res_turb, RFloat **spec_turb, const int nTotalCell, const int nBoundFace,
                                    const int n_turb);
void TestGPUTurbBoundary(RFloat **org_q_turb, const int nTotal, const int n_turb);
void TestGPUTurbLoadResiduals(RFloat **org_res, const int nTotal, const int n_turb);
void TestGPUTurbFillField(RFloat **org_field, const int nTotalCell, const int nTotal, const int neqn);
void TestGPUTurbLoadQ(RFloat **org_qq_turb, const int nTotalCell, const int nTotal, const int n_turb);
void TestGPUTurbStoreRhsByResidual(RFloat **rhs, const int nTotal, const int n_turb);
void TestGPUTurbZeroResiduals(RFloat **res, const int numberOfGeneralCells, const int numberOfEquations);
void TestGPUSpecTurbInit(RFloat **spec_turb, const int nTotal, const int n_turb);
void TestGPUTurbLhs(const int nTotalCell, const int nTotal, const int n_turb, RFloat **dq_turb);
void TestGPUTurbUpdateFlowField(const int nTotalCell, const int nTotal, const int n_turb, RFloat **q, const int n_neg,
                                const int n_pos);
void TestCallGPURecoverResidual(const int nTotal, const int n_turb, RFloat **res);
void TestGPUMinMax(const int nTotal, const int number, const RFloat *dmin, const RFloat *dmax);
void TestGPULimit(const int nTotal, const RFloat *limit);
void TestGPUVisSpectrum(const int nTotalCell, const RFloat *visSpectrum);
void TestGPUInvSpectrum(const int nTotalCell, const RFloat *org_invSpectrum);
void TestComputeNodeValue(const int nTotalNode, const int nBoundFace, const int nEquation, RFloat **org_qNode,
                          RFloat **org_tNode, const int *org_nodeBC, RFloat **org_qOnBCFace);

void TestNodeValue(const int nTotalNode, const int nEquation, RFloat **org_qNode, RFloat **org_tNode);
void TestCellValue(const int nTotal, const RDouble *h_data, const RDouble *d_data);
void TestNodeWeight(const int nTotalNode, RFloat *org_NodeWeight);

void TestComputeQTurbNodeValue(const int nTotalNode, const int nEquation, RFloat **qNode, const int *nCount);
void TestGPUInterPointsAllocCopy(int nIpoint, const int *org_interPoint2GlobalPoint,
                                 const int *org_cellNumberOfInterPoint, const int *org_labelOfInterPoint);
void TestGPUInterFace2BoundaryFace(const int nIFace, const int *org_interFace2BoundaryFace);
//! noting that TestGPUUploadInterfaceValue should be used with TestFGValue
void TestGPUUploadInterfaceValue(RFloat **org_fg, const int neqn, const int nIFace, const string &name);
void TestFGValue(const RFloat *org_fg0, const RFloat *deviceIFVar, const int nIFace, const int neqn,
                 const string &name);
void TestDeviceUploadInterfaceValueHostLargeBuffer(RFloat **org_fg, const int neqn, const int nIFace,
                                                   const string &name);
void TestGPUUploadInterfaceValueLargeBuffer(RFloat **org_fg, const int neqn, const int nIFace, const string &name);
void TestFGValueLargeBuffer(const RFloat *org_fg0, const double *HostSendBuffer, const int nIFace, const int neqn,
                            const string &name, const int offset);
void TestUploadInterfaceValueDeviceToHost(const int *interFace2BoundaryFace, const int *leftCellOfFace, RFloat **fg,
                                          RFloat **f, const int nIFace, const int neqn, const string &name);
void TestGPUDownloadInterfaceValue(const int *interFace2BoundaryFace, const int *rightCellOfFace, RFloat **f,
                                   const int nIFace, const int neqn, const int nTotal, const string &name);
void TestFieldValue(RFloat **f, const RFloat *deviceFieldVar, const int *interFace2BoundaryFace,
                    const int *rightCellOfFace, const int nIFace, const int nTotal, const int neqn, const string &name);
//! noting that  TestGPUUploadInterpointValue should be used with
//! TestFGValueInterpoint
void TestGPUUploadInterpointValue(RFloat **org_fg, const int neqn, const int nIPoint, const string &name);
void TestDeviceUploadInterpointValueHostLargeBuffer(RFloat **org_fg, const int neqn, const int nIPoint,
                                                    const string &name);
void TestGPUUploadInterpointValueLargeBuffer(RFloat **org_fg, const int neqn, const int nIPoint, const string &name);
void TestFGValueInterpointLargeBuffer(const RFloat *org_fg0, const double *HostSendBuffer, const int nIPoint,
                                      const int neqn, const string &name, const int offset);
void TestFGValueInterpoint(const RFloat *org_fg0, const RFloat *deviceIFVar, const int nIFace, const int neqn,
                           const string &name);
void TestUploadInterpointValueDeviceToHost(const int *interPoint2GlobalPoint, RFloat **fg, RFloat **f, const int neqn,
                                           const int nIPoint, const string &name);
void TestGPUDownloadInterpointValue(const int *interPoint2GlobalPoint, RFloat **f, const int nIPoint, const int neqn,
                                    const int nTotalNode, const string &name);
void TestFieldValueInterpoint(RFloat **f, const RFloat *deviceFieldVar, const int *interPoint2GlobalPoint,
                              const int nIPoint, const int nTotalNode, const int neqn, const string &name);
void TestGPUModifyNodeValue(RFloat **org_qNode, RFloat **org_tNode, const int *interPoint2GlobalPoint,
                            const int nIPoint, const int neqn, const int nTotalNode);
void TestModifyQTurbNodeValue(const int *interPoint2GlobalPoint, RFloat **org_qTurbNode, const int nIPoint,
                              const int neqn, const int nTotalNode);
void TestCompressAndTransferInterpolateDataNS(const int iZone, const int jZone, const int nameID, const int ngbID,
                                              const int nIFace, const int nIFaceOfNeighbor,
                                              const int offsetNgbLengthBuffer, const int *faceIndexForSend,
                                              const double *sendBufferNSAWARE, const double *sendBufferNS);
void TestCompressAndTransferInterpointDataNS(const int iZone, const int jZone, const int nameID, const int ngbID,
                                             const int nIPoint, const int nIPointOfNeighbor,
                                             const int offsetNgbLengthBuffer, const int *pointIndexForSend,
                                             const double *sendBufferNSAWARE, const double *sendBufferNS);
void TestCompressAndTransferInterpointDataTurb(const int iZone, const int jZone, const int nameID, const int ngbID,
                                               const int nIPoint, const int nIPointOfNeighbor,
                                               const int offsetNgbLengthBuffer, const int *pointIndexForSend,
                                               const double *sendBufferNSAWARE, const double *sendBufferNS);
void TestCompressAndTransferInterpolateDataTurb(const int iZone, const int jZone, const int nameID, const int ngbID,
                                                const int nIFace, const int nIFaceOfNeighbor,
                                                const int offsetNgbLengthBuffer, const int *faceIndexForSend,
                                                const double *sendBufferNSAWARE, const double *sendBufferNS);
void TestGPUMPIDataDecompressInterpointNS(const int iZone, const int jZone, const int nameID, const int ngbID,
                                          const int nIPoint, const int nIPointOfNeighbor,
                                          const int offsetNgbLengthBuffer, const int *pointIndexForRecv,
                                          const double *recvBufferNSAWARE, const double *recvBufferNS);
void TestGPUMPIDataDecompressInterpointTurb(const int iZone, const int jZone, const int nameID, const int ngbID,
                                            const int nIPoint, const int nIPointOfNeighbor,
                                            const int offsetNgbLengthBuffer, const int *pointIndexForRecv,
                                            const double *recvBufferNSAWARE, const double *recvBufferNS);
void TestGPUMPIDataDecompressNS(const int iZone, const int jZone, const int nameID, const int ngbID, const int nIFace,
                                const int nIFaceOfNeighbor, const int offsetNgbLengthBuffer,
                                const int *faceIndexForRecv, const double *HostRecvBufferNSAWARE,
                                const double *HostRecvLargeBuffer);
void TestGPUMPIDataDecompressTurb(const int iZone, const int jZone, const int nameID, const int ngbID, const int nIFace,
                                  const int nIFaceOfNeighbor, const int offsetNgbLengthBuffer,
                                  const int *faceIndexForRecv, const double *HostRecvBufferNSAWARE,
                                  const double *HostRecvLargeBuffer);
void TestAllreduceMinDt(const double h_localMinDt, const double globalMinDt);
void TestGPUInitSpec(const RFloat *org_spec, const int nTotal);
void TestGPUDq(const int nTotal, const int nl, RFloat **org_dq);
void TestGPUTurbDq(const int nTotal, const int nl, RFloat **org_dq);
void TestSetNSResTmpbyRes(const int nTotal, const int nTotalCell, const int nEquation, RFloat **resTmp);
void TestGPUNSDualTimeSourceRes(const int nTotal, const int nTotalCell, const int nEquation, RFloat **res);
void TestGPUTurbSetResTmpByRes(const int nTotal, const int nTotalCell, const int nEquation, RFloat **resTmp);
void TestGPUTurbDualTimeSourceRes(const int nTotal, const int nTotalCell, const int nEquation, RFloat **res);
void TestGPUNSUnsteadyConvergence(const RFloat sum1, const RFloat sum2);
void TestGPUTurbUnsteadyConvergence(const RFloat sum1, const RFloat sum2);
void TestCallGPUNSUpdateUnsteadyFlow(const int nTotal, const int nEquation, RFloat **qn1, RFloat **qn2, RFloat **resn1,
                                     RFloat **resn2);
void TestCallGPUTurbUpdateUnsteadyFlow(const int nTotal, const int nEquation, RFloat **qn1, RFloat **qn2,
                                       RFloat **resn1, RFloat **resn2);
