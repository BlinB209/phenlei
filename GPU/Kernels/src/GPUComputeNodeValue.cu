#include "GPUComputeNodeValue.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#include "Geo_SimpleBC.h"
#include "GPUBasicFunctionsPart2.h"
#include "GPUKernelTestPart2.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif
using namespace std;
using namespace PHENGLEI;
namespace GPUKernels
{
    void CallGPUComputeNodeValue(Grid *gridIn, int nEquation)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUControlVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalNode = grid->GetNTotalNode();

        InterpointInformation *interPointInfor = grid->GetInterpointInfo();

        int gridSize      = 1;
        int blockSize     = 1;
        int loopLen       = nTotalNode;
        int sizeqOnBCFace = (nEquation + 1) * nBoundFace * sizeof(RFloat);
        HANDLE_API_ERR(cudaMalloc((void **)&d_qOnBCFace, sizeqOnBCFace));
        HANDLE_API_ERR(cudaMalloc((void **)&d_nodeBC, nTotalNode * sizeof(int)));

        cudaMemset(d_qOnBCFace, 0, sizeqOnBCFace);
        cudaMemset(d_nodeBC, 0, nTotalNode * sizeof(int));

        loopLen = nTotalNode;
        KernelLaunchPara((void *)GPUNcountQNodeTNodeInit, loopLen, 0, gridSize, blockSize);
        GPUNodeWeightQNodeTNodeInit<<<gridSize, blockSize>>>(nTotalNode, nEquation, d_qNode, d_tNode, d_nodeWeight);

        //! Node BC type.
        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUnodeBCForOtherBC, loopLen, 0, gridSize, blockSize);

        GPUnodeBCForSoildSurface<<<gridSize, blockSize>>>(nBoundFace, d_node_number_of_each_face, d_face2node,
                                                          d_nodePosiFace, d_nodeBC, d_boundaryType,
                                                          PHENGLEI::SOLID_SURFACE);

        GPUnodeBCForFarField<<<gridSize, blockSize>>>(nBoundFace, d_node_number_of_each_face, d_face2node,
                                                      d_nodePosiFace, d_nodeBC, d_boundaryType, PHENGLEI::SOLID_SURFACE,
                                                      PHENGLEI::FARFIELD);

        GPUnodeBCForOtherBC<<<gridSize, blockSize>>>(nBoundFace, d_node_number_of_each_face, d_face2node,
                                                     d_nodePosiFace, d_nodeBC, d_boundaryType, PHENGLEI::SOLID_SURFACE,
                                                     PHENGLEI::FARFIELD, PHENGLEI::SYMMETRY, PHENGLEI::INTERFACE);

        //! Init q on BC faces.
        GPUqOnBCFaceIni<<<gridSize, blockSize>>>(nEquation, nBoundFace, nTotalCell, d_left_cell_of_face,
                                                 d_right_cell_of_face, d_qOnBCFace, d_q_ns, d_t_proxy, d_boundaryType,
                                                 PHENGLEI::SYMMETRY, PHENGLEI::INTERFACE);

        //! Step1: solid surface.
        GPUqtWeightNodeForSolidSurface<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, d_node_number_of_each_face, d_qOnBCFace, d_face2node, d_nodePosiFace,
            d_qNode, d_tNode, d_nodeWeight, d_boundaryType, PHENGLEI::SOLID_SURFACE, d_x, d_y, d_z, d_xfc, d_yfc,
            d_zfc);

        //! Step2: Far-field.
        GPUqtWeightNodeForFarField<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, d_node_number_of_each_face, d_qOnBCFace, d_face2node, d_nodePosiFace,
            d_qNode, d_tNode, d_nodeWeight, d_nodeBC, d_boundaryType, PHENGLEI::SOLID_SURFACE, PHENGLEI::FARFIELD, d_x,
            d_y, d_z, d_xfc, d_yfc, d_zfc);

        //! Step3: other BC except solid surface/far field/symmetry/Interface.
        GPUqtWeightNodeForOtherBC<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, d_node_number_of_each_face, d_qOnBCFace, d_face2node, d_nodePosiFace,
            d_qNode, d_tNode, d_nodeWeight, d_nodeBC, d_boundaryType, PHENGLEI::SOLID_SURFACE, PHENGLEI::FARFIELD,
            PHENGLEI::SYMMETRY, PHENGLEI::INTERFACE, d_x, d_y, d_z, d_xfc, d_yfc, d_zfc);

        //! Step4: Now the interior points.
        loopLen = nTotalCell;
        KernelLaunchPara((void *)GPUqtWeightNodeForInterioFace, loopLen, 0, gridSize, blockSize);
        GPUqtWeightNodeForInterioFace<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, d_qNode, d_tNode, d_nodeWeight, d_nodeBC, d_q_ns, d_t_proxy,
            d_node_number_of_each_cell, d_cell2Node, d_posiCell2Node, d_x, d_y, d_z, d_xcc, d_ycc, d_zcc);

        //! Change the number of points divided by the interpoint to the number of points before partition,because it is a face loop,
        //! each point is calculated twice,so the number of points needs to be multiplied by 2.
        if (interPointInfor)
        {
            int nIPoint = interPointInfor->GetNumberOfInterpoints();
            loopLen     = nIPoint;
            KernelLaunchPara((void *)GPUqNodetNodeReset, loopLen, 0, gridSize, blockSize);

            GPUqNodetNodeReset<<<gridSize, blockSize>>>(nIPoint, nEquation, nTotalNode, d_qNode, d_tNode,
                                                        d_interPoint2GlobalPoint, d_cellNumberOfInterPoint,
                                                        d_labelOfInterPoint);
            GPUqIPtIPSetZero<<<gridSize, blockSize>>>(nIPoint, nEquation, nTotalNode, d_qInterPoint, d_tInterPoint,
                                                      d_interPoint2GlobalPoint, d_cellNumberOfInterPoint,
                                                      d_labelOfInterPoint);
        }

        HANDLE_API_ERR(cudaFree(d_nodeBC));
        HANDLE_API_ERR(cudaFree(d_qOnBCFace));
    }

    __global__ void GPUnodeBCForSoildSurface(const int nBoundFace, const int *node_number_of_each_face,
                                             const int *face2node, const long long int *nodepos, int *d_nodeBC,
                                             const int *boundaryType, const int SOLID_SURFACE)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            int bcType = boundaryType[iface];
            if (bcType == SOLID_SURFACE)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; j++)
                {
                    int pt = face2node[nodepos[iface] + j];
                    atomicExch(&d_nodeBC[pt], bcType);
                }
            }
        }
    }
    __global__ void GPUnodeBCForFarField(const int nBoundFace, const int *node_number_of_each_face,
                                         const int *face2node, const long long int *nodepos, int *d_nodeBC,
                                         const int *boundaryType, const int SOLID_SURFACE, const int FARFIELD)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            int bcType = boundaryType[iface];
            if (bcType == FARFIELD)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; j++)
                {
                    int pt = face2node[nodepos[iface] + j];
                    if (d_nodeBC[pt] == SOLID_SURFACE)
                    {
                        continue;
                    }
                    atomicExch(&d_nodeBC[pt], bcType);
                }
            }
        }
    }
    __global__ void GPUnodeBCForOtherBC(const int nBoundFace, const int *node_number_of_each_face, const int *face2node,
                                        const long long int *nodepos, int *d_nodeBC, const int *boundaryType,
                                        const int SOLID_SURFACE, const int FARFIELD, const int SYMMETRY,
                                        const int INTERFACE)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            int bcType = boundaryType[iface];
            if (bcType != SOLID_SURFACE && bcType != FARFIELD && bcType != SYMMETRY && bcType != INTERFACE)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; j++)
                {
                    int pt = face2node[nodepos[iface] + j];
                    if (d_nodeBC[pt] == SOLID_SURFACE || d_nodeBC[pt] == FARFIELD)
                    {
                        continue;
                    }
                    atomicExch(&d_nodeBC[pt], bcType);
                }
            }
        }
    }

    __global__ void GPUqOnBCFaceIni(const int nEquation, const int nBoundFace, const int nTotalCell,
                                    const int *left_cell_of_face, const int *right_cell_of_face, RFloat *d_qOnBCFace,
                                    const RFloat *q, const RFloat *t, const int *boundaryType, const int SYMMETRY,
                                    const int INTERFACE)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            int le     = left_cell_of_face[iface];
            int re     = right_cell_of_face[iface];
            int bcType = boundaryType[iface];
            if (bcType == INTERFACE || bcType == SYMMETRY)
            {
                continue;
            }

            for (int m = 0; m < nEquation; ++m)
            {
                d_qOnBCFace[m * nBoundFace + iface] =
                    0.5 * (q[re + m * (nTotalCell + nBoundFace)] + q[le + m * (nTotalCell + nBoundFace)]);
            }
            d_qOnBCFace[nEquation * nBoundFace + iface] =
                0.5 * (t[re + 0 * (nTotalCell + nBoundFace)] + t[le + 0 * (nTotalCell + nBoundFace)]);
        }
    }

    __global__ void GPUqtWeightNodeForSolidSurface(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                   const int *node_number_of_each_face, RFloat *d_qOnBCFace,
                                                   const int *face2node, const long long int *nodepos, RFloat *qNode,
                                                   RFloat *tNode, RDouble *nodeWeight, const int *boundaryType,
                                                   const int SOLID_SURFACE, const RDouble *d_x, const RDouble *d_y,
                                                   const RDouble *d_z, const RDouble *d_xfc, const RDouble *d_yfc,
                                                   const RDouble *d_zfc)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            int bcType = boundaryType[iface];
            if (bcType == SOLID_SURFACE)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; j++)
                {
                    int     pt = face2node[nodepos[iface] + j];
                    RDouble dx = d_x[pt] - d_xfc[iface];
                    RDouble dy = d_y[pt] - d_yfc[iface];
                    RDouble dz = d_z[pt] - d_zfc[iface];

                    RDouble dist       = sqrt(GPUSQR(dx, dy, dz));
                    RDouble weightTemp = 1.0 / dist;

                    for (int m = 0; m < nEquation; ++m)
                    {
                        RDouble faceQtemp = d_qOnBCFace[m * nBoundFace + iface];
#if __CUDA_ARCH__ < 600
                        atomicAddTest(qNode + m * nTotalNode + pt, faceQtemp * weightTemp);
#else
                        atomicAdd(qNode + m * nTotalNode + pt, faceQtemp * weightTemp);
#endif
                    }
#if __CUDA_ARCH__ < 600
                    atomicAddTest(tNode + 0 * nTotalNode + pt,
                                  d_qOnBCFace[nEquation * nBoundFace + iface] * weightTemp);
#else
                    atomicAdd(tNode + 0 * nTotalNode + pt, d_qOnBCFace[nEquation * nBoundFace + iface] * weightTemp);
#endif
                    atomicAdd(nodeWeight + pt, weightTemp);
                }
            }
        }
    }

    __global__ void GPUqtWeightNodeForFarField(const int nTotalNode, const int nEquation, const int nBoundFace,
                                               const int *node_number_of_each_face, RFloat *d_qOnBCFace,
                                               const int *face2node, const long long int *nodepos, RFloat *qNode,
                                               RFloat *tNode, RDouble *nodeWeight, const int *d_nodeBC,
                                               const int *boundaryType, const int SOLID_SURFACE, const int FARFIELD,
                                               const RDouble *d_x, const RDouble *d_y, const RDouble *d_z,
                                               const RDouble *d_xfc, const RDouble *d_yfc, const RDouble *d_zfc)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            int bcType = boundaryType[iface];
            if (bcType == FARFIELD)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; j++)
                {
                    int pt = face2node[nodepos[iface] + j];
                    if (d_nodeBC[pt] == SOLID_SURFACE)
                    {
                        continue;
                    }
                    RDouble dx = d_x[pt] - d_xfc[iface];
                    RDouble dy = d_y[pt] - d_yfc[iface];
                    RDouble dz = d_z[pt] - d_zfc[iface];

                    RDouble dist       = sqrt(GPUSQR(dx, dy, dz));
                    RDouble weightTemp = 1.0 / dist;

                    for (int m = 0; m < nEquation; ++m)
                    {
                        RDouble faceQtemp = d_qOnBCFace[m * nBoundFace + iface];
#if __CUDA_ARCH__ < 600
                        atomicAddTest(qNode + m * nTotalNode + pt, faceQtemp * weightTemp);
#else
                        atomicAdd(qNode + m * nTotalNode + pt, faceQtemp * weightTemp);
#endif
                    }
#if __CUDA_ARCH__ < 600
                    atomicAddTest(tNode + 0 * nTotalNode + pt,
                                  d_qOnBCFace[nEquation * nBoundFace + iface] * weightTemp);
#else
                    atomicAdd(tNode + 0 * nTotalNode + pt, d_qOnBCFace[nEquation * nBoundFace + iface] * weightTemp);
#endif
                    atomicAdd(nodeWeight + pt, weightTemp);
                }
            }
        }
    }
    __global__ void GPUqtWeightNodeForOtherBC(const int nTotalNode, const int nEquation, const int nBoundFace,
                                              const int *node_number_of_each_face, RFloat *d_qOnBCFace,
                                              const int *face2node, const long long int *nodepos, RFloat *qNode,
                                              RFloat *tNode, RDouble *nodeWeight, const int *d_nodeBC,
                                              const int *boundaryType, const int SOLID_SURFACE, const int FARFIELD,
                                              const int SYMMETRY, const int INTERFACE, const RDouble *d_x,
                                              const RDouble *d_y, const RDouble *d_z, const RDouble *d_xfc,
                                              const RDouble *d_yfc, const RDouble *d_zfc)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            int bcType = boundaryType[iface];
            if (bcType != SOLID_SURFACE && bcType != FARFIELD && bcType != SYMMETRY && bcType != INTERFACE)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; j++)
                {
                    int pt = face2node[nodepos[iface] + j];
                    if (d_nodeBC[pt] == SOLID_SURFACE || d_nodeBC[pt] == FARFIELD)
                    {
                        continue;
                    }
                    RDouble dx = d_x[pt] - d_xfc[iface];
                    RDouble dy = d_y[pt] - d_yfc[iface];
                    RDouble dz = d_z[pt] - d_zfc[iface];

                    RDouble dist       = sqrt(GPUSQR(dx, dy, dz));
                    RDouble weightTemp = 1.0 / dist;

                    for (int m = 0; m < nEquation; ++m)
                    {
                        RDouble faceQtemp = d_qOnBCFace[m * nBoundFace + iface];
#if __CUDA_ARCH__ < 600
                        atomicAddTest(qNode + m * nTotalNode + pt, faceQtemp * weightTemp);
#else
                        atomicAdd(qNode + m * nTotalNode + pt, faceQtemp * weightTemp);
#endif
                    }
#if __CUDA_ARCH__ < 600
                    atomicAddTest(tNode + 0 * nTotalNode + pt,
                                  d_qOnBCFace[nEquation * nBoundFace + iface] * weightTemp);
#else
                    atomicAdd(tNode + 0 * nTotalNode + pt, d_qOnBCFace[nEquation * nBoundFace + iface] * weightTemp);
#endif
                    atomicAdd(nodeWeight + pt, weightTemp);
                }
            }
        }
    }

    __global__ void GPUqtWeightNodeForInterioFace(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                  const int nTotalCell, RFloat *qNode, RFloat *tNode,
                                                  RDouble *nodeWeight, const int *d_nodeBC, RFloat *q, RFloat *t,
                                                  const int *nodeNumberOfEachCell, const int *cell2Node,
                                                  const int *cell2NodePosition, const RDouble *d_x, const RDouble *d_y,
                                                  const RDouble *d_z, const RDouble *d_xcc, const RDouble *d_ycc,
                                                  const RDouble *d_zcc)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iCell = i; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            int nNode        = nodeNumberOfEachCell[iCell];
            int nodePosition = cell2NodePosition[iCell];
            for (int jNode = 0; jNode < nNode; ++jNode)
            {
                int point = cell2Node[nodePosition + jNode];
                if (d_nodeBC[point] != 0)
                {
                    continue;
                }

                //! NO_BC && Symmetry && Interface?
                RDouble dx         = d_x[point] - d_xcc[iCell];
                RDouble dy         = d_y[point] - d_ycc[iCell];
                RDouble dz         = d_z[point] - d_zcc[iCell];
                RDouble dist       = sqrt(GPUSQR(dx, dy, dz));
                RDouble weightTemp = 1.0 / dist;

                for (int m = 0; m < nEquation; ++m)
                {
                    RDouble cellQtemp = q[m * (nTotalCell + nBoundFace) + iCell];
#if __CUDA_ARCH__ < 600
                    atomicAddTest(qNode + m * nTotalNode + point, cellQtemp * weightTemp);
#else
                    atomicAdd(qNode + m * nTotalNode + point, cellQtemp * weightTemp);
#endif
                }
#if __CUDA_ARCH__ < 600
                atomicAddTest(tNode + 0 * nTotalNode + point, t[0 * (nTotalCell + nBoundFace) + iCell] * weightTemp);
#else
                atomicAdd(tNode + 0 * nTotalNode + point, t[0 * (nTotalCell + nBoundFace) + iCell] * weightTemp);
#endif
                atomicAdd(nodeWeight + point, weightTemp);
            }
        }
    }

    void CallGPUModifyNodeValue(Grid *gridIn, int nEquation)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUControlVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalNode = grid->GetNTotalNode();

        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nTotalNode;

        loopLen = nTotalNode;
        KernelLaunchPara((void *)GPUModifyqNodetNode, loopLen, 0, gridSize, blockSize);
        GPUModifyqNodetNode<<<gridSize, blockSize>>>(nTotalNode, nEquation, d_qNode, d_tNode, d_nodeWeight);
    }
    __global__ void GPUModifyqNodetNode(const int nTotalNode, const int nEquation, RFloat *qNode, RFloat *tNode,
                                        RDouble *nodeWeight)
    {
        int           i     = blockDim.x * blockIdx.x + threadIdx.x;
        const RDouble small = 1.0e-36f;
        for (int iNode = i; iNode < nTotalNode; iNode += gridDim.x * blockDim.x)
        {
            if (nodeWeight[iNode] > small)
            {
                for (int m = 0; m < nEquation; m++)
                {
                    qNode[m * nTotalNode + iNode] /= nodeWeight[iNode];
                }
                tNode[0 * nTotalNode + iNode] /= nodeWeight[iNode];
            }
            else
            {
                for (int m = 0; m < nEquation; m++)
                {
                    qNode[m * nTotalNode + iNode] = 0.0;
                }
                tNode[0 * nTotalNode + iNode] = 0.0;
            }
        }
    }

    void CallGPUComputeNodeValueInterpolation(const int nTotalNode, const int nBoundFace, const int nTotalFace,
                                              const int nTotalCell, const int nEquation, const int nIPoint)
    {
#ifdef KERNELLAUNCHTEST
        cout << "Test in CallGPUCompNodeVarByGradient with gradient_field_proxy=" << d_gradient_field_proxy
             << "and d_gradient_var_index=" << d_gradient_var_index << endl;
#endif

        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUControlVariables;
        using namespace GPUGeomVariables;
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nTotalNode;

        KernelLaunchPara((void *)GPUNcountQNodeTNodeInit, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNcountQNodeTNodeInit, 0, gridSize, blockSize);
#endif

        GPUNcountQNodeTNodeInit<<<gridSize, blockSize>>>(nTotalNode, nEquation, d_qNode, d_tNode, d_nCount);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUBoundaryFaceNCountQNodeTNodeCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUBoundaryFaceNCountQNodeTNodeCal, 0, gridSize, blockSize);
#endif

#ifndef CUDADEBUGATOMIC
        GPUBoundaryFaceNCountQNodeTNodeCal<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, d_left_cell_of_face, d_right_cell_of_face,
            d_node_number_of_each_face, d_face2node, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount,
            d_boundaryType);
#else
        GPUBoundaryFaceNCountQNodeTNodeCal<<<1, 1>>>(nTotalNode, nEquation, nBoundFace, nTotalCell, d_left_cell_of_face,
                                                     d_right_cell_of_face, d_node_number_of_each_face, d_face2node,
                                                     d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount,
                                                     d_boundaryType);
#endif

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
#ifndef LOOPMODEOPT
        loopLen = nTotalFace - nBoundFace;
        KernelLaunchPara((void *)GPUInteriorFaceNCountQNodeTNodeCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUInteriorFaceNCountQNodeTNodeCal, 0, gridSize, blockSize);
#endif
#else
        loopLen = nTotalCell;
        /*
                    KernelLaunchPara((void*)GPUCellLoopNCountQNodeTNodeCalFinal, loopLen, 0, gridSize, blockSize);
                    #ifdef KERNELLAUNCHTEST
                            ReportKernelPara((void*)GPUCellLoopNCountQNodeTNodeCalFinal, 0, gridSize, blockSize);
                    #endif
            */
        KernelLaunchPara((void *)GPUCellLoopQNodeCalFinalSep, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCellLoopQNodeCalFinalSep, 0, gridSize, blockSize);
#endif
#endif

#ifndef CUDADEBUGATOMIC
#ifndef LOOPMODEOPT
        GPUInteriorFaceNCountQNodeTNodeCal<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, nTotalFace, d_left_cell_of_face, d_right_cell_of_face,
            d_node_number_of_each_face, d_face2node, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount);
#else
        //!GPUCellLoopNCountQNodeTNodeCalFinal<<<gridSize, blockSize>>>(nTotalNode, nEquation, nTotalCell+nBoundFace, nTotalCell, nTotalFace, d_node_number_of_each_cell, d_cell2Node, d_posiCell2Node, d_cell2NodeCount, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount);
        //!seperation loop
        for (int equationID = 0; equationID < nEquation; equationID++)
        {
            GPUCellLoopQNodeCalFinalSep<<<gridSize, blockSize>>>(
                nTotalNode, equationID, nTotalCell + nBoundFace, nTotalCell, nTotalFace, d_node_number_of_each_cell,
                d_cell2Node, d_posiCell2Node, d_cell2NodeCount, d_q_ns, d_qNode);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
        GPUCellLoopTNodeCalFinalSep<<<gridSize, blockSize>>>(nTotalNode, nTotalCell + nBoundFace, nTotalCell,
                                                             nTotalFace, d_node_number_of_each_cell, d_cell2Node,
                                                             d_posiCell2Node, d_cell2NodeCount, d_t_proxy, d_tNode);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        GPUCellLoopNCountCalFinalSep<<<gridSize, blockSize>>>(nTotalNode, nTotalCell + nBoundFace, nTotalCell,
                                                              nTotalFace, d_node_number_of_each_cell, d_cell2Node,
                                                              d_posiCell2Node, d_cell2NodeCount, d_nCount);
#endif //!for ifndef LOOPMODEOPT
#else
#ifndef LOOPMODEOPT
        GPUInteriorFaceNCountQNodeTNodeCal<<<1, 1>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, nTotalFace, d_left_cell_of_face, d_right_cell_of_face,
            d_node_number_of_each_face, d_face2node, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount);
#else
        //!for whole loop
        //!GPUCellLoopNCountQNodeTNodeCalFinal<<<1, 1>>>(nTotalNode, nEquation, nTotalCell+nBoundFace, nTotalCell, nTotalFace, d_node_number_of_each_cell, d_cell2Node, d_posiCell2Node, d_cell2NodeCount, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount);
        //!seperation loop
        for (int equationID = 0; equationID < nEquation; equationID++)
        {
            GPUCellLoopQNodeCalFinalSep<<<1, 1>>>(nTotalNode, equationID, nTotalCell + nBoundFace, nTotalCell,
                                                  nTotalFace, d_node_number_of_each_cell, d_cell2Node, d_posiCell2Node,
                                                  d_cell2NodeCount, d_q_ns, d_qNode);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
        GPUCellLoopTNodeCalFinalSep<<<1, 1>>>(nTotalNode, nTotalCell + nBoundFace, nTotalCell, nTotalFace,
                                              d_node_number_of_each_cell, d_cell2Node, d_posiCell2Node,
                                              d_cell2NodeCount, d_t_proxy, d_tNode);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        GPUCellLoopNCountCalFinalSep<<<1, 1>>>(nTotalNode, nTotalCell + nBoundFace, nTotalCell, nTotalFace,
                                               d_node_number_of_each_cell, d_cell2Node, d_posiCell2Node,
                                               d_cell2NodeCount, d_nCount);
#endif //!for ifndef LOOPMODEOPT
#endif //!for ifndef CUDADEBUGATOMIC

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        CallGPUComputeNodeValueBoundaryTreatMent(nEquation, nTotalNode, nBoundFace, nTotalCell);
        //!for multi-process case, nIPoint > 0, more treatment for interpolation points is required.
        if (nIPoint) CallGPUMultiProcessReset(nIPoint, nTotalNode, nEquation);

        loopLen = nTotalNode;
        KernelLaunchPara((void *)GPUComputeNodeValueAvr, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeNodeValueAvr, 0, gridSize, blockSize);
#endif
        GPUComputeNodeValueAvr<<<gridSize, blockSize>>>(nEquation, nTotalNode, d_qNode, d_tNode, d_nCount);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    void CallGPUMultiProcessReset(const int nIPoint, const int nTotalNode, const int nEquation)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nIPoint;

        KernelLaunchPara((void *)GPUNCountReset, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNCountReset, 0, gridSize, blockSize);
#endif
        GPUNCountReset<<<gridSize, blockSize>>>(nIPoint, d_nCount, d_interPoint2GlobalPoint, d_cellNumberOfInterPoint);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        KernelLaunchPara((void *)GPUqNodetNodeReset, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)qNodetNodeReset, 0, gridSize, blockSize);
#endif
        GPUqNodetNodeReset<<<gridSize, blockSize>>>(nIPoint, nEquation, nTotalNode, d_qNode, d_tNode,
                                                    d_interPoint2GlobalPoint, d_cellNumberOfInterPoint,
                                                    d_labelOfInterPoint);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        KernelLaunchPara((void *)GPUqIPtIPSetZero, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUqIPtIPSetZero, 0, gridSize, blockSize);
#endif
        GPUqIPtIPSetZero<<<gridSize, blockSize>>>(nIPoint, nEquation, nTotalNode, d_qInterPoint, d_tInterPoint,
                                                  d_interPoint2GlobalPoint, d_cellNumberOfInterPoint,
                                                  d_labelOfInterPoint);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUNCountReset(const int nIPoint, int *nCount, const int *interPoint2GlobalPoint,
                                   const int *cellNumberOfInterPoint)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int iPoint   = 0;
        int globalPoint;
        for (iPoint = threadID; iPoint < nIPoint; iPoint += blockDim.x * gridDim.x)
        {
            globalPoint         = interPoint2GlobalPoint[iPoint];
            nCount[globalPoint] = cellNumberOfInterPoint[iPoint];
        }
    }

    __global__ void GPUqNodetNodeReset(const int nIPoint, const int nEquation, const int nTotalNode, RFloat *qNode,
                                       RFloat *tNode, const int *interPoint2GlobalPoint,
                                       const int *cellNumberOfInterPoint, const int *labelOfInterPoint)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int iPoint   = 0;
        int globalPoint;
        int m = 0;
        for (iPoint = threadID; iPoint < nIPoint; iPoint += blockDim.x * gridDim.x)
        {
            if (labelOfInterPoint[iPoint] != 1)
            {
                globalPoint = interPoint2GlobalPoint[iPoint];
                for (m = 0; m < nEquation; ++m)
                {
                    qNode[m * nTotalNode + globalPoint] = 0;
                }
                tNode[0 * nTotalNode + globalPoint] = 0;
            }
        }
    }

    __global__ void GPUqIPtIPSetZero(const int nIPoint, const int nEquation, const int nTotalNode, RFloat *qInterPoint,
                                     RFloat *tInterPoint, const int *interPoint2GlobalPoint,
                                     const int *cellNumberOfInterPoint, const int *labelOfInterPoint)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int iPoint   = 0;
        int m        = 0;
        for (iPoint = threadID; iPoint < nIPoint; iPoint += blockDim.x * gridDim.x)
        {
            for (m = 0; m < nEquation; ++m)
            {
                qInterPoint[m * nIPoint + iPoint] = 0;
            }
            tInterPoint[0 * nIPoint + iPoint] = 0;
        }
    }
    void CallGPUComputeNodeValueBoundaryTreatMent(const int nEquation, const int nTotalNode, const int nBoundFace,
                                                  const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        using namespace GPUControlVariables;
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nBoundFace;
        int blockSizeLimit = 0;

        //!On Symmetry boundary
        KernelLaunchPara((void *)GPUComputeNodeValueBoundaryInit, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeNodeValueBoundaryInit, 0, gridSize, blockSize);
#endif

        GPUComputeNodeValueBoundaryInit<<<gridSize, blockSize>>>(
            nEquation, nBoundFace, nTotalNode, d_face2node, d_node_number_of_each_face, d_boundaryType, d_nodePosiFace,
            d_qNode, d_tNode, d_nCount, PHENGLEI::SYMMETRY);

        KernelLaunchPara((void *)GPUComputeNodeValueBoundaryCal, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeNodeValueBoundaryCal, 0, gridSize, blockSize);
#endif

#ifndef CUDADEBUGATOMIC
        GPUComputeNodeValueBoundaryCal<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, d_left_cell_of_face, d_face2node, d_node_number_of_each_face,
            d_boundaryType, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount, PHENGLEI::SYMMETRY);
#else
        GPUComputeNodeValueBoundaryCal<<<1, 1>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, d_left_cell_of_face, d_face2node, d_node_number_of_each_face,
            d_boundaryType, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount, PHENGLEI::SYMMETRY);
#endif

        //!On FARFIELD boundary
        KernelLaunchPara((void *)GPUComputeNodeValueBoundaryInit, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeNodeValueBoundaryInit, 0, gridSize, blockSize);
#endif

        GPUComputeNodeValueBoundaryInit<<<gridSize, blockSize>>>(
            nEquation, nBoundFace, nTotalNode, d_face2node, d_node_number_of_each_face, d_boundaryType, d_nodePosiFace,
            d_qNode, d_tNode, d_nCount, PHENGLEI::FARFIELD);

        KernelLaunchPara((void *)GPUComputeNodeValueBoundaryCal, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeNodeValueBoundaryCal, 0, gridSize, blockSize);
#endif

#ifndef CUDADEBUGATOMIC
        GPUComputeNodeValueBoundaryCal<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, d_left_cell_of_face, d_face2node, d_node_number_of_each_face,
            d_boundaryType, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount, PHENGLEI::FARFIELD);
        //!Add by Zhang Xi for checking effects of grid&block size on results: it is found that enlarging the grid or block size, the precision goes down.
#else
        GPUComputeNodeValueBoundaryCal<<<1, 1>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, d_left_cell_of_face, d_face2node, d_node_number_of_each_face,
            d_boundaryType, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount, PHENGLEI::FARFIELD);
#endif

        KernelLaunchPara((void *)GPUComputeNodeValueBoundaryInit, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeNodeValueBoundaryInit, 0, gridSize, blockSize);
#endif

        GPUComputeNodeValueBoundaryInit<<<gridSize, blockSize>>>(
            nEquation, nBoundFace, nTotalNode, d_face2node, d_node_number_of_each_face, d_boundaryType, d_nodePosiFace,
            d_qNode, d_tNode, d_nCount, PHENGLEI::SOLID_SURFACE);

        KernelLaunchPara((void *)GPUComputeNodeValueBoundaryCal, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeNodeValueBoundaryCal, 0, gridSize, blockSize);
#endif

#ifndef CUDADEBUGATOMIC
        GPUComputeNodeValueBoundaryCal<<<gridSize, blockSize>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, d_left_cell_of_face, d_face2node, d_node_number_of_each_face,
            d_boundaryType, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount, PHENGLEI::SOLID_SURFACE);
#else
        GPUComputeNodeValueBoundaryCal<<<1, 1>>>(
            nTotalNode, nEquation, nBoundFace, nTotalCell, d_left_cell_of_face, d_face2node, d_node_number_of_each_face,
            d_boundaryType, d_nodePosiFace, d_q_ns, d_t_proxy, d_qNode, d_tNode, d_nCount, PHENGLEI::SOLID_SURFACE);
#endif
    }

    __global__ void GPUComputeNodeValueBoundaryInit(const int nEquation, const int nBoundFace, const int nTotalNode,
                                                    const int *face2node, const int *node_number_of_each_face,
                                                    const int *boundaryType, const long long int *nodepos,
                                                    RFloat *qNode, RFloat *tNode, int *nCount, const int bctype_in)
    {
        int pt, bctype, m, j;
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            bctype = boundaryType[iface];
            if (bctype == bctype_in)
            {
                for (j = 0; j < node_number_of_each_face[iface]; ++j)
                {
                    pt = face2node[nodepos[iface] + j];
                    //!The judge is deleted because INTERFACE does not exist in the case
                    for (m = 0; m < nEquation; ++m)
                    {
                        qNode[m * nTotalNode + pt] = 0;
                    }
                    tNode[0 * nTotalNode + pt] = 0;
                    nCount[pt]                 = 0;
                }
            }
        }
    }

    __global__ void GPUComputeNodeValueBoundaryCal(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                   const int nTotalCell, const int *left_cell_of_face,
                                                   const int *face2node, const int *node_number_of_each_face,
                                                   const int *boundaryType, const long long int *nodepos,
                                                   const RFloat *q, const RFloat *t, RFloat *qNode, RFloat *tNode,
                                                   int *nCount, const int bctype_in)
    {
        //!Try different type on q_n
        int le, re, pt, bctype;
        int iface = 0;
        int j     = 0;
        int m     = 0;
        //!double qle, qre;
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            le     = left_cell_of_face[iface];
            re     = iface + nTotalCell;
            bctype = boundaryType[iface];
            if (bctype == bctype_in)
            {
                for (j = 0; j < node_number_of_each_face[iface]; ++j)
                {
                    pt = face2node[nodepos[iface] + j];
                    for (m = 0; m < nEquation; ++m)
                    {
//!qNode[m][point] += q[m][le];
#if __CUDA_ARCH__ < 600
                        atomicAddTest(qNode + m * nTotalNode + pt, q[le + m * (nTotalCell + nBoundFace)]);
#else
                        atomicAdd(qNode + m * nTotalNode + pt, q[le + m * (nTotalCell + nBoundFace)]);
#endif
                    }

//!tNode[0][point] += t[0][le];
#if __CUDA_ARCH__ < 600
                    atomicAddTest(tNode + 0 * nTotalNode + pt, t[le + 0 * (nTotalCell + nBoundFace)]);
#else
                    atomicAdd(tNode + 0 * nTotalNode + pt, t[le + 0 * (nTotalCell + nBoundFace)]);
#endif

                    //!nCount[point] += 1;
                    atomicAdd(nCount + pt, 1);

                    for (int m = 0; m < nEquation; ++m)
                    {
//!qNode[m][point] += q[m][re];
#if __CUDA_ARCH__ < 600
                        atomicAddTest(qNode + m * nTotalNode + pt, q[re + m * (nTotalCell + nBoundFace)]);
#else
                        atomicAdd(qNode + m * nTotalNode + pt, q[re + m * (nTotalCell + nBoundFace)]);
#endif
                    }

//!tNode[0][point] += t[0][re];
#if __CUDA_ARCH__ < 600
                    atomicAddTest(tNode + 0 * nTotalNode + pt, t[re + 0 * (nTotalCell + nBoundFace)]);
#else
                    atomicAdd(tNode + 0 * nTotalNode + pt, t[re + 0 * (nTotalCell + nBoundFace)]);
#endif
                    atomicAdd(nCount + pt, 1);
                }
            }
        }
    }

    __global__ void GPUNcountQNodeTNodeInit(const int nTotalNode, const int nEquation, RFloat *qNode, RFloat *tNode,
                                            int *nCount)
    {
        int i     = blockIdx.x * blockDim.x + threadIdx.x;
        int iNode = 0;

        for (iNode = i; iNode < nTotalNode; iNode += blockDim.x * gridDim.x)
        {
            nCount[iNode] = 0;
            for (int m = 0; m < nEquation; ++m)
            {
                //!qNode[m][iNode] = 0;
                qNode[m * nTotalNode + iNode] = 0;
            }
            //!tNode[0][iNode] = 0;
            tNode[0 * nTotalNode + iNode] = 0;
        }
    }
    __global__ void GPUNodeWeightQNodeTNodeInit(const int nTotalNode, const int nEquation, RFloat *qNode, RFloat *tNode,
                                                RFloat *nodeWeight)
    {
        int i     = blockIdx.x * blockDim.x + threadIdx.x;
        int iNode = 0;

        for (iNode = i; iNode < nTotalNode; iNode += blockDim.x * gridDim.x)
        {
            nodeWeight[iNode] = 0;
            for (int m = 0; m < nEquation; ++m)
            {
                //!qNode[m][iNode] = 0;
                qNode[m * nTotalNode + iNode] = 0;
            }
            //!tNode[0][iNode] = 0;
            tNode[0 * nTotalNode + iNode] = 0;
        }
    }

    __global__ void GPUBoundaryFaceNCountQNodeTNodeCal(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                       const int nTotalCell, const int *left_cell_of_face,
                                                       const int *right_cell_of_face,
                                                       const int *node_number_of_each_face, const int *face2node,
                                                       const long long int *nodepos, const RFloat *q, const RFloat *t,
                                                       RFloat *qNode, RFloat *tNode, int *nCount,
                                                       const int *boundaryType)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int le, re, pt, j, bcType, m;
        for (iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            le     = left_cell_of_face[iface];
            re     = right_cell_of_face[iface];
            bcType = boundaryType[iface];
            if (bcType != PHENGLEI::INTERFACE)
            {
                for (j = 0; j < node_number_of_each_face[iface]; j++)
                {
                    pt = face2node[nodepos[iface] + j];
                    for (m = 0; m < nEquation; ++m)
                    {
#if __CUDA_ARCH__ < 600
                        atomicAddTest(qNode + m * nTotalNode + pt, q[le + m * (nTotalCell + nBoundFace)]);
#else
                        atomicAdd(qNode + m * nTotalNode + pt, q[le + m * (nTotalCell + nBoundFace)]);
#endif
                    }

#if __CUDA_ARCH__ < 600
                    atomicAddTest(tNode + 0 * nTotalNode + pt, t[le + 0 * (nTotalCell + nBoundFace)]);
#else
                    atomicAdd(tNode + 0 * nTotalNode + pt, t[le + 0 * (nTotalCell + nBoundFace)]);
#endif
                    atomicAdd(nCount + pt, 1);

                    //! From right
                    for (int m = 0; m < nEquation; ++m)
                    {
#if __CUDA_ARCH__ < 600
                        atomicAddTest(qNode + m * nTotalNode + pt, q[re + m * (nTotalCell + nBoundFace)]);
#else
                        atomicAdd(qNode + m * nTotalNode + pt, q[re + m * (nTotalCell + nBoundFace)]);
#endif
                    }

#if __CUDA_ARCH__ < 600
                    atomicAddTest(tNode + 0 * nTotalNode + pt, t[re + 0 * (nTotalCell + nBoundFace)]);
#else
                    atomicAdd(tNode + 0 * nTotalNode + pt, t[re + 0 * (nTotalCell + nBoundFace)]);
#endif
                    atomicAdd(nCount + pt, 1);
                }
            }
            else
            {
                for (j = 0; j < node_number_of_each_face[iface]; j++)
                {
                    pt = face2node[nodepos[iface] + j];

                    for (m = 0; m < nEquation; ++m)
                    {
#if __CUDA_ARCH__ < 600
                        atomicAddTest(qNode + m * nTotalNode + pt, q[le + m * (nTotalCell + nBoundFace)]);
#else
                        atomicAdd(qNode + m * nTotalNode + pt, q[le + m * (nTotalCell + nBoundFace)]);
#endif
                    }

#if __CUDA_ARCH__ < 600
                    atomicAddTest(tNode + 0 * nTotalNode + pt, t[le + 0 * (nTotalCell + nBoundFace)]);
#else
                    atomicAdd(tNode + 0 * nTotalNode + pt, t[le + 0 * (nTotalCell + nBoundFace)]);
#endif
                    atomicAdd(nCount + pt, 1);
                }
            }
        }
    }

    __global__ void GPUCellLoopQNodeCalFinalSep(const int nTotalNode, const int equationID, const int nTotal,
                                                const int nTotalCell, const int nTotalFace,
                                                const int *nodeNumberOfEachCell, const int *cell2Node,
                                                const int *cell2NodePosition, const int *cell2NodeCount,
                                                const RFloat *qNS, RFloat *qNode)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int nodePosition, numNodesInCell, nodeOffset, cellID;
        int nodeID;
        int accessFrequency;
        for (cellID = threadID; cellID < nTotalCell; cellID += blockDim.x * gridDim.x)
        {
            nodePosition   = cell2NodePosition[cellID];
            numNodesInCell = nodeNumberOfEachCell[cellID];

            for (nodeOffset = 0; nodeOffset < numNodesInCell; nodeOffset++)
            {
                nodeID          = cell2Node[nodePosition + nodeOffset];
                accessFrequency = cell2NodeCount[nodePosition + nodeOffset];
#if __CUDA_ARCH__ < 600
                atomicAddTest(qNode + equationID * nTotalNode + nodeID,
                              qNS[cellID + equationID * nTotal] * accessFrequency);
#else
                atomicAdd(qNode + equationID * nTotalNode + nodeID,
                          qNS[cellID + equationID * nTotal] * accessFrequency);
#endif
            }
        }
    }

    __global__ void GPUCellLoopTNodeCalFinalSep(const int nTotalNode, const int nTotal, const int nTotalCell,
                                                const int nTotalFace, const int *nodeNumberOfEachCell,
                                                const int *cell2Node, const int *cell2NodePosition,
                                                const int *cell2NodeCount, const RFloat *tCell, RFloat *tNode)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int nodePosition, numNodesInCell, nodeOffset, cellID;
        int nodeID;
        int accessFrequency;
        for (cellID = threadID; cellID < nTotalCell; cellID += blockDim.x * gridDim.x)
        {
            nodePosition   = cell2NodePosition[cellID];
            numNodesInCell = nodeNumberOfEachCell[cellID];

            for (nodeOffset = 0; nodeOffset < numNodesInCell; nodeOffset++)
            {
                nodeID          = cell2Node[nodePosition + nodeOffset];
                accessFrequency = cell2NodeCount[nodePosition + nodeOffset];

#if __CUDA_ARCH__ < 600
                atomicAddTest(tNode + 0 * nTotalNode + nodeID, tCell[cellID + 0 * nTotal] * accessFrequency);
#else
                atomicAdd(tNode + 0 * nTotalNode + nodeID, tCell[cellID + 0 * nTotal] * accessFrequency);
#endif
            }
        }
    }

    __global__ void GPUCellLoopNCountCalFinalSep(const int nTotalNode, const int nTotal, const int nTotalCell,
                                                 const int nTotalFace, const int *nodeNumberOfEachCell,
                                                 const int *cell2Node, const int *cell2NodePosition,
                                                 const int *cell2NodeCount, int *nCount)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int nodePosition, numNodesInCell, nodeOffset, cellID;
        int nodeID;
        int accessFrequency;

        for (cellID = threadID; cellID < nTotalCell; cellID += blockDim.x * gridDim.x)
        {
            nodePosition   = cell2NodePosition[cellID];
            numNodesInCell = nodeNumberOfEachCell[cellID];

            for (nodeOffset = 0; nodeOffset < numNodesInCell; nodeOffset++)
            {
                nodeID          = cell2Node[nodePosition + nodeOffset];
                accessFrequency = cell2NodeCount[nodePosition + nodeOffset];
                atomicAdd(nCount + nodeID, 1 * accessFrequency);
            }
        }
    }

    __global__ void GPUCellLoopNCountQNodeTNodeCalFinal(const int nTotalNode, const int nEquation, const int nTotal,
                                                        const int nTotalCell, const int nTotalFace,
                                                        const int *nodeNumberOfEachCell, const int *cell2Node,
                                                        const int *cell2NodePosition, const int *cell2NodeCount,
                                                        const RFloat *qNS, const RFloat *tCell, RFloat *qNode,
                                                        RFloat *tNode, int *nCount)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int nodePosition, numNodesInCell, nodeOffset, cellID;
        int nodeID;
        int accessFrequency;

        for (cellID = threadID; cellID < nTotalCell; cellID += blockDim.x * gridDim.x)
        {
            nodePosition   = cell2NodePosition[cellID];
            numNodesInCell = nodeNumberOfEachCell[cellID];
            for (nodeOffset = 0; nodeOffset < numNodesInCell; nodeOffset++)
            {
                nodeID          = cell2Node[nodePosition + nodeOffset];
                accessFrequency = cell2NodeCount[nodePosition + nodeOffset];
                for (int m = 0; m < nEquation; m++)
                {
//!qNode[m][nodeID] += qNS[m][cellID];
#if __CUDA_ARCH__ < 600
                    atomicAddTest(qNode + m * nTotalNode + nodeID, qNS[cellID + m * nTotal] * accessFrequency);
#else
                    atomicAdd(qNode + m * nTotalNode + nodeID, qNS[cellID + m * nTotal] * accessFrequency);
#endif
                }

//!tNode[0][nodeID] += tCell[0][cellID];
#if __CUDA_ARCH__ < 600
                atomicAddTest(tNode + 0 * nTotalNode + nodeID, tCell[cellID + 0 * nTotal] * accessFrequency);
#else
                atomicAdd(tNode + 0 * nTotalNode + nodeID, tCell[cellID + 0 * nTotal] * accessFrequency);
#endif
                atomicAdd(nCount + nodeID, 1 * accessFrequency);
            }
        }
    }

    __global__ void GPUInteriorFaceNCountQNodeTNodeCal(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                       const int nTotalCell, const int nTotalFace,
                                                       const int *left_cell_of_face, const int *right_cell_of_face,
                                                       const int *node_number_of_each_face, const int *face2node,
                                                       const long long int *nodepos, const RFloat *q, const RFloat *t,
                                                       RFloat *qNode, RFloat *tNode, int *nCount)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int le, re, pt, j, m;

        for (iface = i + nBoundFace; iface < nTotalFace; iface += gridDim.x * blockDim.x)
        {
            le = left_cell_of_face[iface];
            re = right_cell_of_face[iface];
            for (j = 0; j < node_number_of_each_face[iface]; j++)
            {
                pt = face2node[nodepos[iface] + j];

                for (m = 0; m < nEquation; m++)
                {
//!qNode[m][point] += q[m][le];
#if __CUDA_ARCH__ < 600
                    atomicAddTest(qNode + m * nTotalNode + pt, q[le + m * (nTotalCell + nBoundFace)]);
#else
                    atomicAdd(qNode + m * nTotalNode + pt, q[le + m * (nTotalCell + nBoundFace)]);
#endif
                }

//!tNode[0][point] += t[0][le];
#if __CUDA_ARCH__ < 600
                atomicAddTest(tNode + 0 * nTotalNode + pt, t[le + 0 * (nTotalCell + nBoundFace)]);
#else
                atomicAdd(tNode + 0 * nTotalNode + pt, t[le + 0 * (nTotalCell + nBoundFace)]);
#endif
                //!nCount[point] += 1;
                atomicAdd(nCount + pt, 1);

                for (m = 0; m < nEquation; ++m)
                {
//!qNode[m][point] += q[m][re];
#if __CUDA_ARCH__ < 600
                    atomicAddTest(qNode + m * nTotalNode + pt, q[re + m * (nTotalCell + nBoundFace)]);
#else
                    atomicAdd(qNode + m * nTotalNode + pt, q[re + m * (nTotalCell + nBoundFace)]);
#endif
                }

                //!tNode[0][point] += t[0][re];
#if __CUDA_ARCH__ < 600
                atomicAddTest(tNode + 0 * nTotalNode + pt, t[re + 0 * (nTotalCell + nBoundFace)]);
#else
                atomicAdd(tNode + 0 * nTotalNode + pt, t[re + 0 * (nTotalCell + nBoundFace)]);
#endif
                //!nCount[point] += 1;
                atomicAdd(nCount + pt, 1);
            }
        }
    }

    __global__ void GPUComputeNodeValueAvr(const int nEquation, const int nTotalNode, RFloat *qNode, RFloat *tNode,
                                           int *nCount)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int inode = 0;

        for (inode = i; inode < nTotalNode; inode += blockDim.x * gridDim.x)
        {
            if (nCount[inode] != 0)
            {
                for (int m = 0; m < nEquation; m++)
                {
                    //!qNode[m][iNode] /= nCount[iNode];
                    qNode[m * nTotalNode + inode] /= nCount[inode];
                }

                //!tNode[0][iNode] /= nCount[iNode];
                tNode[0 * nTotalNode + inode] /= nCount[inode];
            }
            else
            {
                for (int m = 0; m < nEquation; m++)
                {
                    //!qNode[m][iNode] = 0.0;
                    qNode[m * nTotalNode + inode] = 0.0;
                }

                //!tNode[0][iNode] = 0.0;
                tNode[0 * nTotalNode + inode] = 0.0;
            }
        }
    }

    void CallGPUModifyNodeValue(const int nIPoint, const int nEquation, const int nTotalNode)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        int gridSize  = 1;
        int blockSize = 1;
        //!int loopLen= nTotalNode;

        int loopLen = nIPoint;
        KernelLaunchPara((void *)GPUModifyNodeValue, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUModifyNodeValue, 0, gridSize, blockSize);
#endif

#ifdef CUDADEBUGATOMIC
        GPUModifyNodeValue<<<1, 1>>>(nIPoint, nEquation, nTotalNode, d_qNode, d_tNode, d_qInterPoint, d_tInterPoint,
                                     d_interPoint2GlobalPoint);
#else
        GPUModifyNodeValue<<<gridSize, blockSize>>>(nIPoint, nEquation, nTotalNode, d_qNode, d_tNode, d_qInterPoint,
                                                    d_tInterPoint, d_interPoint2GlobalPoint);
#endif

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUModifyNodeValue(const int nIPoint, const int nEquation, const int nTotalNode, RFloat *qNode,
                                       RFloat *tNode, RFloat *qInterPoint, RFloat *tInterPoint,
                                       const int *interPoint2GlobalPoint)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        int iPoint   = 0;
        int GlobalPoint, m;

        for (iPoint = threadID; iPoint < nIPoint; iPoint += blockDim.x * gridDim.x)
        {
            GlobalPoint = interPoint2GlobalPoint[iPoint];
            for (m = 0; m < nEquation; m++)
            {
//!qNode[m * nTotalNode + GlobalPoint] += qInterPoint[m * nIPoint + iPoint];
#if __CUDA_ARCH__ < 600
                atomicAddTest(qNode + m * nTotalNode + GlobalPoint, qInterPoint[m * nIPoint + iPoint]);
#else
                atomicAdd(qNode + m * nTotalNode + GlobalPoint, qInterPoint[m * nIPoint + iPoint]);
#endif
            }

//!tNode[0 * nTotalNode + GlobalPoint] += tInterPoint[0 * nIPoint + iPoint];
#if __CUDA_ARCH__ < 600
            atomicAddTest(tNode + 0 * nTotalNode + GlobalPoint, tInterPoint[0 * nIPoint + iPoint]);
#else
            atomicAdd(tNode + 0 * nTotalNode + GlobalPoint, tInterPoint[0 * nIPoint + iPoint]);
#endif
        }
    }
} //! namespace GPUKernels
