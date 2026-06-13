#include "GPUCompGradientGGNodeNew.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif

using namespace std;
namespace GPUKernels
{
    void CallGPUCompGradientGGNode_NEW(const string q_name, const int d_index, const int nTotalNode,
                                       const int nTotalCell, const int nBoundFace, const int nTotalFace)
    {
        RFloat *dqdx;
        RFloat *dqdy;
        RFloat *dqdz;
        RFloat *q_n;
        RFloat *q;

        int index;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        if (q_name == "Empty")
        {
            return;
            //！Just temporary code block
        }

        if (q_name == "d_q_proxy")
        {
            dqdx  = d_dqdx_proxy;
            dqdy  = d_dqdy_proxy;
            dqdz  = d_dqdz_proxy;
            q_n   = d_qNode;
            index = d_index;
            q     = d_q_ns;
        }
        else if (q_name == "d_t_proxy")
        {
            dqdx  = d_dtdx_proxy;
            dqdy  = d_dtdy_proxy;
            dqdz  = d_dtdz_proxy;
            q_n   = d_tNode;
            index = d_index;
            q     = d_t_proxy;
        }
        else if (q_name == "d_q_turb_proxy")
        {
            dqdx  = d_dq_turbdx_proxy;
            dqdy  = d_dq_turbdy_proxy;
            dqdz  = d_dq_turbdz_proxy;
            q_n   = d_qTurbNode;
            index = d_index;
        }
        else if (q_name == "d_velocity_proxy")
        {
            dqdx = d_dveldx_proxy;
            dqdy = d_dveldy_proxy;
            dqdz = d_dveldz_proxy;
            //！set d_qVelNode by d_qNode firstly
            CallGPUSetVelocityNode(nTotalNode);
            q_n   = d_qVelNode;
            index = d_index;
        }
        else
        {
            return;
        }
        int nTotal = nTotalCell + nBoundFace;
        /* original kernel size 
        gridSize=2092;
        blockSize=640;
        */
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUCompGradientGGNodeNewInit, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeNewInit, 0, gridSize, blockSize);
#endif

        //！set zero into dqdx, dqdy, dqdz
        GPUCompGradientGGNodeNewInit<<<gridSize, blockSize>>>(index, nTotal, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGNodeNewBoundaryFaceCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeNewBoundaryFaceCal, 0, gridSize, blockSize);
#endif
#ifndef CUDADEBUGATOMIC
        GPUCompGradientGGNodeNewBoundaryFaceCal<<<gridSize, blockSize>>>(
            nTotalNode, index, nTotalCell, nBoundFace, d_left_cell_of_face, d_node_number_of_each_face, d_face2node,
            d_xfn, d_yfn, d_zfn, d_area, d_nodePosiFace, q_n, dqdx, dqdy, dqdz, q, d_boundaryType, PHENGLEI::SYMMETRY,
            PHENGLEI::INTERFACE);
#else
        GPUCompGradientGGNodeNewBoundaryFaceCal<<<1, 1>>>(
            nTotalNode, index, nTotalCell, nBoundFace, d_left_cell_of_face, d_node_number_of_each_face, d_face2node,
            d_xfn, d_yfn, d_zfn, d_area, d_nodePosiFace, q_n, dqdx, dqdy, dqdz, q, d_boundaryType, PHENGLEI::SYMMETRY,
            PHENGLEI::INTERFACE);
#endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        loopLen = nTotalFace - nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGNodeNewInteriorFaceCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeNewInteriorFaceCal, 0, gridSize, blockSize);
#endif

#ifndef CUDADEBUGATOMIC
        GPUCompGradientGGNodeNewInteriorFaceCal<<<gridSize, blockSize>>>(
            nTotalNode, index, nTotalCell, nTotalFace, nBoundFace, d_left_cell_of_face, d_right_cell_of_face,
            d_node_number_of_each_face, d_face2node, d_xfn, d_yfn, d_zfn, d_area, d_nodePosiFace, q_n, dqdx, dqdy,
            dqdz);
#else
        GPUCompGradientGGNodeNewInteriorFaceCal<<<1, 1>>>(nTotalNode, index, nTotalCell, nTotalFace, nBoundFace,
                                                          d_left_cell_of_face, d_right_cell_of_face,
                                                          d_node_number_of_each_face, d_face2node, d_xfn, d_yfn, d_zfn,
                                                          d_area, d_nodePosiFace, q_n, dqdx, dqdy, dqdz);
#endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = nTotalCell;
        KernelLaunchPara((void *)GPUCompGradientGGNodeNewDiVol, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeNewDiVol, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGNodeNewDiVol<<<gridSize, blockSize>>>(index, nTotalCell, nTotal, dqdx, dqdy, dqdz, d_vol);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGNodeNewGhostSet, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeNewGhostSet, 0, gridSize, blockSize);
#endif
        GPUCompGradientGGNodeNewGhostSet<<<gridSize, blockSize>>>(index, nTotalCell, nTotal, nBoundFace,
                                                                  d_left_cell_of_face, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUSetVelocityNode(const int nTotalNode)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalNode;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUSetVelocityNode, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSetVelocityNode, 0, gridSize, blockSize);
#endif
        GPUSetVelocityNode<<<gridSize, blockSize>>>(nTotalNode, d_qNode, d_qVelNode);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPUSetVelocityNode(const int nTotalNode, const RFloat *qNode, RFloat *qVelNode)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iNode = 0;
        for (iNode = i; iNode < nTotalNode; iNode += blockDim.x * gridDim.x)
        {
            qVelNode[0 * nTotalNode + iNode] = qNode[1 * nTotalNode + iNode];
            qVelNode[1 * nTotalNode + iNode] = qNode[2 * nTotalNode + iNode];
            qVelNode[2 * nTotalNode + iNode] = qNode[3 * nTotalNode + iNode];
        }
    }
    __global__ void GPUCompGradientGGNodeNewInit(const int index, const int nTotal, RFloat *dqdx, RFloat *dqdy,
                                                 RFloat *dqdz)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            dqdx[icell + index * nTotal] = 0.0;
            dqdy[icell + index * nTotal] = 0.0;
            dqdz[icell + index * nTotal] = 0.0;
        }
    }

    __global__ void GPUCompGradientGGNodeNewBoundaryFaceCal(const int nTotalNode, const int index, const int nTotalCell,
                                                            const int nBoundFace, const int *left_cell_of_face,
                                                            const int *node_number_of_each_face, const int *face2node,
                                                            const RDouble *nxs, const RDouble *nys, const RDouble *nzs,
                                                            const RDouble *ns, const long long int *nodepos,
                                                            const double *q_n, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz,
                                                            const RFloat *q, const int *boundaryType,
                                                            const int d_SYMMETRY, const int d_INTERFACE)
    {
        int le, re, pt;
        //！double qfc, qfcnx, qfcny, qfcnz;
        RFloat  qfc, qfcnx, qfcny, qfcnz; //！Try different varialbe type
        RDouble nx, ny, nz;
        int     i     = blockDim.x * blockIdx.x + threadIdx.x;
        int     iface = 0;
        for (iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            le         = left_cell_of_face[iface];
            re         = iface + nTotalCell;
            int bcType = boundaryType[iface];
            int nTotal = nTotalCell + nBoundFace;
            qfc        = 0.0;
            if (bcType == d_SYMMETRY || bcType == d_INTERFACE)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; ++j)
                {
                    pt = face2node[nodepos[iface] + j];
                    //！qfc += q_n[pt];
                    qfc += q_n[index * nTotalNode + pt];
                }
                qfc /= node_number_of_each_face[iface];
            }
            else
            {
                qfc = 0.5 * (q[index * nTotal + le] + q[index * nTotal + re]);
            }

            nx = nxs[iface] * ns[iface];
            ny = nys[iface] * ns[iface];
            nz = nzs[iface] * ns[iface];

            qfcnx = qfc * nx;
            qfcny = qfc * ny;
            qfcnz = qfc * nz;

#if __CUDA_ARCH__ < 600
            atomicAddTest(dqdx + (nTotalCell + nBoundFace) * index + le, qfcnx);
            atomicAddTest(dqdy + (nTotalCell + nBoundFace) * index + le, qfcny);
            atomicAddTest(dqdz + (nTotalCell + nBoundFace) * index + le, qfcnz);
#else
            atomicAdd(dqdx + (nTotalCell + nBoundFace) * index + le, qfcnx);
            atomicAdd(dqdy + (nTotalCell + nBoundFace) * index + le, qfcny);
            atomicAdd(dqdz + (nTotalCell + nBoundFace) * index + le, qfcnz);
#endif
        }
    }

    __global__ void GPUCompGradientGGNodeNewInteriorFaceCal(const int nTotalNode, const int index, const int nTotalCell,
                                                            const int nTotalFace, const int nBoundFace,
                                                            const int *left_cell_of_face, const int *right_cell_of_face,
                                                            const int *node_number_of_each_face, const int *face2node,
                                                            const RDouble *nxs, const RDouble *nys, const RDouble *nzs,
                                                            const RDouble *ns, const long long int *nodepos,
                                                            const double *q_n, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz)
    {
        int le, re, pt;
        //！double qfc, qfcnx, qfcny, qfcnz;
        RFloat  qfc, qfcnx, qfcny, qfcnz; //！Try different type
        RDouble nx, ny, nz;
        int     i     = blockDim.x * blockIdx.x + threadIdx.x;
        int     iface = 0;
        for (iface = i + nBoundFace; iface < nTotalFace; iface += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iface];
            re = right_cell_of_face[iface];

            qfc = 0.0;
            for (int j = 0; j < node_number_of_each_face[iface]; ++j)
            {
                pt = face2node[nodepos[iface] + j];
                //！qfc += q_n[pt];
                qfc += q_n[index * nTotalNode + pt];
            }

            qfc /= node_number_of_each_face[iface];

            nx = nxs[iface] * ns[iface];
            ny = nys[iface] * ns[iface];
            nz = nzs[iface] * ns[iface];

            qfcnx = qfc * nx;
            qfcny = qfc * ny;
            qfcnz = qfc * nz;

#if __CUDA_ARCH__ < 600
            atomicAddTest(dqdx + (nTotalCell + nBoundFace) * index + le, qfcnx);
            atomicAddTest(dqdy + (nTotalCell + nBoundFace) * index + le, qfcny);
            atomicAddTest(dqdz + (nTotalCell + nBoundFace) * index + le, qfcnz);

            atomicAddTest(dqdx + (nTotalCell + nBoundFace) * index + re, -qfcnx);
            atomicAddTest(dqdy + (nTotalCell + nBoundFace) * index + re, -qfcny);
            atomicAddTest(dqdz + (nTotalCell + nBoundFace) * index + re, -qfcnz);
#else
            atomicAdd(dqdx + (nTotalCell + nBoundFace) * index + le, qfcnx);
            atomicAdd(dqdy + (nTotalCell + nBoundFace) * index + le, qfcny);
            atomicAdd(dqdz + (nTotalCell + nBoundFace) * index + le, qfcnz);

            atomicAdd(dqdx + (nTotalCell + nBoundFace) * index + re, -qfcnx);
            atomicAdd(dqdy + (nTotalCell + nBoundFace) * index + re, -qfcny);
            atomicAdd(dqdz + (nTotalCell + nBoundFace) * index + re, -qfcnz);
#endif
        }
    }

    __global__ void GPUCompGradientGGNodeNewDiVol(const int index, const int nTotalCell, const int nTotal, RFloat *dqdx,
                                                  RFloat *dqdy, RFloat *dqdz, RDouble *vol)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotalCell; icell += blockDim.x * gridDim.x)
        {
            //！just for test
            //！if ((icell == 0)||(icell == 1)) printf("vol[%d]=%e\n", icell, vol[icell]);
            dqdx[icell + index * nTotal] /= vol[icell];
            dqdy[icell + index * nTotal] /= vol[icell];
            dqdz[icell + index * nTotal] /= vol[icell];
        }
    }

    __global__ void GPUCompGradientGGNodeNewGhostSet(const int index, const int nTotalCell, const int nTotal,
                                                     const int nBoundFace, const int *left_cell_of_face, RFloat *dqdx,
                                                     RFloat *dqdy, RFloat *dqdz)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        for (iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            int le = left_cell_of_face[iface];
            int re = iface + nTotalCell;

            dqdx[re + index * nTotal] = dqdx[le + index * nTotal];
            dqdy[re + index * nTotal] = dqdy[le + index * nTotal];
            dqdz[re + index * nTotal] = dqdz[le + index * nTotal];
        }
    }
} //！ namespace GPUKernels
