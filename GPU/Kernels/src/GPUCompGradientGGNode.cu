#include "GPUCompGradientGGNode.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif

using namespace std;
namespace GPUKernels
{
    void CallGPUCompGradientGGNode(const string q_name, const int d_index, const int nTotalCell, const int nBoundFace,
                                   const int nTotalFace)
    {
        RFloat *dqdx;
        RFloat *dqdy;
        RFloat *dqdz;
        RFloat *q_n;
        int     index;
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
            dqdx = d_dqdx_proxy;
            dqdy = d_dqdy_proxy;
            dqdz = d_dqdz_proxy;
            //！q_n = d_q_n_double;
            //！q_n = d_q_n_tmp;
            q_n   = d_qNode;
            index = d_index;
        }
        else if (q_name == "d_t_proxy")
        {
            dqdx = d_dtdx_proxy;
            dqdy = d_dtdy_proxy;
            dqdz = d_dtdz_proxy;
            //！q_n = d_q_n_double;
            //！q_n = d_q_n_tmp;
            q_n   = d_tNode;
            index = d_index;
        }
        else if (q_name == "d_q_turb_proxy")
        {
            dqdx = d_dq_turbdx_proxy;
            dqdy = d_dq_turbdy_proxy;
            dqdz = d_dq_turbdz_proxy;
            //！q_n = d_q_n_double;
            q_n   = d_q_n_tmp;
            index = d_index;
        }
        else if (q_name == "d_velocity_proxy")
        {
            dqdx = d_dveldx_proxy;
            dqdy = d_dveldy_proxy;
            dqdz = d_dveldz_proxy;
            //！q_n = d_q_n_double;
            q_n   = d_q_n_tmp;
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
        KernelLaunchPara((void *)GPUCompGradientGGNodeInit, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeInit, 0, gridSize, blockSize);
#endif

        //！set zero into dqdx, dqdy, dqdz
        GPUCompGradientGGNodeInit<<<gridSize, blockSize>>>(index, nTotal, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGNodeBoundaryFaceCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeBoundaryFaceCal, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGNodeBoundaryFaceCal<<<gridSize, blockSize>>>(
            index, nTotalCell, nBoundFace, d_left_cell_of_face, d_node_number_of_each_face, d_face2node, d_xfn, d_yfn,
            d_zfn, d_area, d_nodePosiFace, q_n, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = nTotalFace - nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGNodeInteriorFaceCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeInteriorFaceCal, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGNodeInteriorFaceCal<<<gridSize, blockSize>>>(
            index, nTotalCell, nTotalFace, nBoundFace, d_left_cell_of_face, d_right_cell_of_face,
            d_node_number_of_each_face, d_face2node, d_xfn, d_yfn, d_zfn, d_area, d_nodePosiFace, q_n, dqdx, dqdy,
            dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = nTotalCell;
        KernelLaunchPara((void *)GPUCompGradientGGNodeDiVol, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeDiVol, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGNodeDiVol<<<gridSize, blockSize>>>(index, nTotalCell, nTotal, dqdx, dqdy, dqdz, d_vol);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGNodeGhostSet, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGNodeGhostSet, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGNodeGhostSet<<<gridSize, blockSize>>>(index, nTotalCell, nTotal, nBoundFace,
                                                               d_left_cell_of_face, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUCompGradientGGNodeInit(const int index, const int nTotal, RFloat *dqdx, RFloat *dqdy,
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

    __global__ void GPUCompGradientGGNodeBoundaryFaceCal(const int index, const int nTotalCell, const int nBoundFace,
                                                         const int *left_cell_of_face,
                                                         const int *node_number_of_each_face, const int *face2node,
                                                         const RDouble *nxs, const RDouble *nys, const RDouble *nzs,
                                                         const RDouble *ns, const long long int *nodepos,
                                                         const double *q_n, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz)
    {
        int le, re, pt;
        //！double qfc, qfcnx, qfcny, qfcnz;
        RFloat  qfc, qfcnx, qfcny, qfcnz; //！Try different varialbe type
        RDouble nx, ny, nz;
        int     i     = blockDim.x * blockIdx.x + threadIdx.x;
        int     iface = 0;
        for (iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iface];
            re = iface + nTotalCell;

            qfc = 0.0;

            for (int j = 0; j < node_number_of_each_face[iface]; ++j)
            {
                pt = face2node[nodepos[iface] + j];
                qfc += q_n[pt];
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
#else
            atomicAdd(dqdx + (nTotalCell + nBoundFace) * index + le, qfcnx);
            atomicAdd(dqdy + (nTotalCell + nBoundFace) * index + le, qfcny);
            atomicAdd(dqdz + (nTotalCell + nBoundFace) * index + le, qfcnz);
#endif
        }
    }

    __global__ void GPUCompGradientGGNodeInteriorFaceCal(const int index, const int nTotalCell, const int nTotalFace,
                                                         const int nBoundFace, const int *left_cell_of_face,
                                                         const int *right_cell_of_face,
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
                qfc += q_n[pt];
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

    __global__ void GPUCompGradientGGNodeDiVol(const int index, const int nTotalCell, const int nTotal, RFloat *dqdx,
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

    __global__ void GPUCompGradientGGNodeGhostSet(const int index, const int nTotalCell, const int nTotal,
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
